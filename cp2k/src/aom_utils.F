!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

!!!!!!!!!! CHANGE_AC !!!!!!!!!!!!!

! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2016.01.15: creation, put only ffsh_aom subroutines
! *****************************************************************************

! *****************************************************************************
!>      jacob: comment
! construct_ham_X_SH_full_degen
!
! This program is designed to propagate a single excess electronic charge
! along a molecular dynamics trajectory, either during the trajectory
! (to do surface hopping simulations) or afterwards (classical path
! approximation).
!
! This program takes molecular coordinates as inputs.
! From these coordinate inputs, it calculates the electronic Hamiltonian
! in a molecular orbital basis, using the overlapFINAL module to calculate
! the off-diagonals and either accepting the diagonals as inputs or calculating
! them in some other way.
!
! With the electronic Hamiltonian, the program then propagates the electron
! using the RK4 method.
! At the final electronic state, the probabilities for hopping between adiabatic
! surfaces are calculated.
! The program then decides what adiabatic surface the system is on and
! calculates
! the appropriate atomic forces.
!
! In the classical path approximation, the above steps are repeated until the
! end of the input molecular trajectory.  For surface-hopping, the data here is
! fed back into
! the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE aom_utils
   USE kinds,                          ONLY: default_string_length,&
                                               dp
   USE machine,                        ONLY: m_getpid
   USE aom_types,                      ONLY: sh_print_env_type, &
                                             aom_env_type
   USE aom_print,                      ONLY: adiab_print_sh_info

  USE cp_parser_methods,               ONLY: read_float_object

  USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                             cp_print_key_finished_output
  USE cp_log_handling,                 ONLY: cp_logger_type,&
                                             cp_logger_release,&
                                             cp_get_default_logger,&
                                             cp_to_string
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_list_get,&
                                             section_vals_val_get
!  USE kinds,                           ONLY: default_string_length,&
!                                            dp
  USE mixed_energy_types,              ONLY: mixed_force_type
  USE particle_list_types,             ONLY: particle_list_type
  USE virial_types,                    ONLY: virial_p_type,&
                                             virial_type,&
                                             zero_virial

  USE mixed_environment_utils,         ONLY: get_subsys_map_index
  USE mixed_environment_types,         ONLY: mixed_environment_type
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE cp_linked_list_val,              ONLY: &
       cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
       cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
       cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
       cp_sll_val_set_el_at, cp_sll_val_type


  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

  USE distribution_1d_types,           ONLY: distribution_1d_type

  USE force_env_types,                 ONLY: force_env_get, &
                                             force_env_type

  USE non_adiabatic_types,             ONLY: sh_env_type

  USE parallel_rng_types,              ONLY: next_random_number 

  USE aom_types,                       ONLY: aom_env_type

!  USE non_adiabatic_utils,             ONLY: simple_vel_rescaling, &
!                                             translate_com,        &
!                                             cal_delta_energy_nacv,&
!                                             diab_to_adiab

#include "./base/base_uses.f90"




IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_utils'

   PUBLIC :: makeX,            &
             calc_coms_and_nlist_aom, &
             make_adiabX,      &
             calcforces,           &
             read_listcoeff,       & 
             calc_adiab_nace, &
             calc_diabatic_nace, &
             read_listdiabcoeff, & 
             calc_adiab_nace_old, &
             matinvsqrt
 

CONTAINS


!**********************************************************************************
!> \brief 
!**********************************************************************************
SUBROUTINE calc_coms_and_nlist_aom(aom_env, atoms, mol, cutoff_sites, orb)

     TYPE(aom_env_type), INTENT(INOUT)     :: aom_env
     INTEGER, INTENT(IN)                   :: mol, atoms, orb
     REAL(KIND=dp), intent (IN)            :: cutoff_sites

     CHARACTER(len=*), PARAMETER                     :: routineN='calc_coms_and_nlist_aom', &
     routineP = moduleN//':'//routineN

     integer                               :: i, j, m, d
  

     !CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     aom_env%coms = 0
     do i=1,mol
          do m=1,atoms
               do d=1,3
                    aom_env%coms(d,i) = aom_env%coms(d,i) + (aom_env%crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do

     !CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               aom_env%dist(i,j) = ((aom_env%coms(1,i) - aom_env%coms(1,j))*(aom_env%coms(1,i) &
                - aom_env%coms(1,j)) + (aom_env%coms(2,i) - aom_env%coms(2,j))*&
               &(aom_env%coms(2,i) - aom_env%coms(2,j)) + (aom_env%coms(3,i) & 
                - aom_env%coms(3,j))*(aom_env%coms(3,i) - aom_env%coms(3,j)))
               if (aom_env%dist(i,j) < cutoff_sites*cutoff_sites) then
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 1
                         end do
                    end do
               else
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 0
                         end do
                    end do
               end if                    
          end do
     end do
END SUBROUTINE calc_coms_and_nlist_aom

  
!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*aom_env%ctin
        !makeX constructs an X matrix from inputs of S and H
        ! the X matrix is = H - d^{diab}, note that elements of H are pure
        ! imaginary = -iH_{kl} due to the construction of S matrix

        implicit none

        integer, intent(in) :: mol
        COMPLEX(KIND=dp), dimension(:,:), intent(in) :: H
        REAL(KIND=dp), dimension(:,:), intent(in) :: Stt, Sttp, Sttm, Stptp, Stmtm
        COMPLEX(KIND=dp), dimension(:,:), intent(out) :: X
        REAL(KIND=dp), intent(in) :: deltat
        REAL(KIND=dp), dimension(:,:), intent(out), optional    :: NACE   

       CHARACTER(len=*), PARAMETER :: routineN ='makeX', &
                               routineP = moduleN//':'//routineN

        REAL(KIND=dp), dimension(:,:), allocatable :: Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt
        REAL(KIND=dp), dimension(:,:), allocatable :: Sterm1, Sterm2, Splace1, Splace2
        REAL(KIND=dp), dimension(:,:), allocatable     :: my_NACE
!TEST
!        integer                                        :: i, j 

        !allocate matrices
        allocate(Sttinvsqrt(mol,mol))
        allocate(Sttsqrt(mol,mol))
        allocate(Stptpinvsqrt(mol,mol))
        allocate(Stmtminvsqrt(mol,mol))
        allocate(Sterm1(mol,mol))
        allocate(Sterm2(mol,mol))
        allocate(Splace1(mol,mol))
        allocate(Splace2(mol,mol))
        allocate(my_NACE(mol,mol))

        Sttsqrt = Stt
        call matinvsqrt(mol, .FALSE., Sttsqrt)
        Sttinvsqrt = Stt
        call matinvsqrt(mol, .TRUE., Sttinvsqrt)
        Stptpinvsqrt = Stptp
        call matinvsqrt(mol, .TRUE., Stptpinvsqrt)
        Stmtminvsqrt = Stmtm
        call matinvsqrt(mol, .TRUE., Stmtminvsqrt)

        !subroutine DGEMM multiplies two matrices, from the BLAS
        !use to build the S terms in the X matrix
       !centering difference are used through this subroutine

       !first term: Sttsqrt*(Stptpinvsqrt - Stmtmtinvsqrt) this the second term
       !in Eq 13 Spencer16 calculated constistently at time t-dt as
       !S^{1/2}_{TT}*(S^{-1/2}_{TpTp} - S^{-1/2}_{TmTm}) 
        Splace1 = Stptpinvsqrt - Stmtminvsqrt
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttsqrt, mol, Splace1, mol, 0.0D0, Sterm1, mol)

       !second term: Sttinvsqrt*(Sttp - Sttm)*Sttinvsqrt
       ! this is actually the first term in Eq 13 Spencer16 calculated
       ! constistently at time t-dt (indeicated by T) as S^{-1/2}_{TT}*D'_{T}*S^{-1/2}_{TT}
        Splace1 = Sttp - Sttm
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Splace1, mol, Sttinvsqrt, mol, 0.0D0, Splace2, mol)
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttinvsqrt, mol, Splace2, mol, 0.0D0, Sterm2, mol)
        !ADD: SYMMETRIZE NACE
        my_NACE = (0.5D0/deltat)*(Sterm1 + Sterm2)

        my_NACE = (my_NACE - TRANSPOSE(my_NACE)) / 2
        !X = H - (0.5D0/deltat)*(Sterm1 + Sterm2)

!!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "DIABATIC NACE:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", my_NACE(i,j), "  "
!   END DO
!   PRINT *, " "
!END DO

!!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "DIABATIC H:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", -AIMAG(H(i,j)), "  "
!   END DO
!   PRINT *, " "
!END DO

        X = H - my_NACE
!TEST 
!        X = H 
!END TEST
        IF (PRESENT(NACE)) THEN
           NACE = my_NACE
        ENDIF

        deallocate(Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt, Sterm1, Sterm2, Splace1, Splace2)
        deallocate(my_NACE)
end subroutine makeX


!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine calc_diabatic_nace(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, deltat, NACE, X)
       !X is the name I give to the matrix in the solution to the Schroedinger equation:
       !dc/dt = X*aom_env%ctin
       !makeX constructs an X matrix from inputs of S and H

       integer, intent(in) :: mol
       COMPLEX(KIND=dp), dimension(:,:), intent(in) :: H
       REAL(KIND=dp), dimension(:,:), intent(in) :: Stt, Sttp, Sttm, Stptp, Stmtm
       REAL(KIND=dp), intent(in) :: deltat
       REAL(KIND=dp), dimension(:,:), intent(out), optional    :: NACE   
       COMPLEX(KIND=dp), dimension(:,:), intent(out), optional :: X

       CHARACTER(len=*), PARAMETER :: routineN ='calc_diabatic_nace', &
                               routineP = moduleN//':'//routineN

       REAL(KIND=dp), dimension(:,:), allocatable :: Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt
       REAL(KIND=dp), dimension(:,:), allocatable :: Sterm1, Sterm2, Splace1, Splace2
       REAL(KIND=dp), dimension(:,:), allocatable     :: diab_NACE   

       !allocate matrices
       allocate(Sttinvsqrt(mol,mol))
       allocate(Sttsqrt(mol,mol))
       allocate(Stptpinvsqrt(mol,mol))
       allocate(Stmtminvsqrt(mol,mol))
       allocate(Sterm1(mol,mol))
       allocate(Sterm2(mol,mol))
       allocate(Splace1(mol,mol))
       allocate(Splace2(mol,mol))
       allocate(diab_NACE(mol,mol))

       Sttsqrt = Stt
       call matinvsqrt(mol, .FALSE., Sttsqrt)
       Sttinvsqrt = Stt
       call matinvsqrt(mol, .TRUE., Sttinvsqrt)
       Stptpinvsqrt = Stptp
       call matinvsqrt(mol, .TRUE., Stptpinvsqrt)
       Stmtminvsqrt = Stmtm
       call matinvsqrt(mol, .TRUE., Stmtminvsqrt)

       !subroutine DGEMM multiplies two matrices, from the BLAS
       !use to build the S terms in the X matrix, 
       !centering difference are used in all case for this subroutine

       !first term: Sttsqrt*(Stptpinvsqrt - Stmtmtinvsqrt) this the second term
       !in Eq 13 Spencer16 calculated constistently at time t-dt as
       !S^{1/2}_{TT}*(S^{-1/2}_{TpTp} - S^{-1/2}_{TmTm}) 
       Splace1 = Stptpinvsqrt - Stmtminvsqrt
       call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttsqrt, mol, Splace1, mol, 0.0D0, Sterm1, mol)

       !second term: Sttinvsqrt*(Sttp - Sttm)*Sttinvsqrt
       ! this is actually the first term in Eq 13 Spencer16 calculated
       ! constistently at time t-dt (indeicated by T) as S^{-1/2}_{TT}*D'_{T}*S^{-1/2}_{TT}
       Splace1 = Sttp - Sttm
       call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Splace1, mol, Sttinvsqrt, mol, 0.0D0, Splace2, mol)
       call DGEMM('N', 'N', mol, mol, mol, 1.0D0, Sttinvsqrt, mol, Splace2, mol, 0.0D0, Sterm2, mol)
       !ADD: SYMMETRIZE NACE
       diab_NACE = (0.5D0/deltat)*(Sterm1 + Sterm2)

       diab_NACE = (diab_NACE - TRANSPOSE(diab_NACE)) / 2
!       X = H - (0.5D0/deltat)*(Sterm1 + Sterm2)

       if (PRESENT(X)) then
          X = H - diab_NACE
       end if

       IF (PRESENT(NACE)) THEN
          NACE = diab_NACE
       ENDIF

       deallocate(Sttinvsqrt, Sttsqrt, Stptpinvsqrt, Stmtminvsqrt, Sterm1, Sterm2, Splace1, Splace2)
       deallocate(diab_NACE)

end subroutine calc_diabatic_nace

!!! WARNING::::Adiabatic propgation is DEPRECATED so this may not work!!!
subroutine make_adiabX(mol, E_t, NACE_ad_t, X) !deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*aom_env%ctin
        !makeX constructs an X matrix from inputs of S and H

        integer, intent(in) :: mol
        REAL(KIND=dp), DIMENSION(:), intent(in)     :: E_t
        REAL(KIND=dp), dimension(:,:), intent(in)     :: NACE_ad_t
        COMPLEX(KIND=dp), dimension(:,:), allocatable, intent(out) :: X

       CHARACTER(len=*), PARAMETER :: routineN ='make_adiabX', &
                               routineP = moduleN//':'//routineN

        COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H
        integer :: i, j
        REAL(KIND=dp), dimension(:,:), allocatable     :: my_NACE   

        allocate(H(mol,mol))
        allocate(my_NACE(mol,mol))
        allocate(X(mol,mol))

        do i=1,mol
            do j=1,mol
               if (j .eq. i) then
                    H(i,j) = cmplx(0.0D0,-1.0D0*E_t(i),dp)
               else
                    H(i,j) =  cmplx(0.0D0, 0.0D0, dp)
               end if
            end do
        end do

        my_NACE = (NACE_ad_t - TRANSPOSE(NACE_ad_t)) / 2
                
        X = H - my_NACE
end subroutine make_adiabX


!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine old_diab_to_adiab( aom_env, C2, old_state, mol, cstate, ct, C1, cttm_ad, cstate_tm)
        REAL(KIND=dp), dimension(:,:), intent(in) :: C2
        integer, intent(in) :: old_state
        TYPE(aom_env_type) :: aom_env
        integer, intent(in) :: mol
        COMPLEX(KIND=dp), intent(out) :: cstate
        COMPLEX(KIND=dp), dimension(:), intent(out) :: ct
        REAL(KIND=dp), dimension(:,:), optional, intent(in) :: C1
        COMPLEX(KIND=dp), dimension(:), optional, intent(inout) ::  cttm_ad
        COMPLEX(KIND=dp), optional, intent(out) :: cstate_tm

        CHARACTER(len=*), PARAMETER :: routineN ='old_diab_to_adiab', &
                               routineP = moduleN//':'//routineN

        integer :: i, j!, mol!, k, l

        cstate = cmplx(0.0D0, 0.0D0, dp)
        IF (PRESENT(cstate_tm)) cstate_tm = cmplx(0.0D0, 0.0D0, dp)
        do i=1,mol
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*aom_env%ctin(i)
                IF (PRESENT(C1)) THEN
                     cstate_tm = cstate_tm + cmplx(C1(i,old_state),0.0D0,dp)*aom_env%ctin(i)
                END IF
        end do
        do i=1,mol
                !calculate a and b
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                IF (PRESENT(cttm_ad)) cttm_ad(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
                        ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*aom_env%ctin(j)
                        IF (PRESENT(C1)) THEN
                              cttm_ad(i) = cttm_ad(i) + cmplx(C1(j,i),0.0D0,dp)*aom_env%ctin(j)
                        END IF
                end do
        end do
end subroutine old_diab_to_adiab



!**********************************************************************************
!> \brief Calculate the diabatic NACE first
!**********************************************************************************
subroutine nace1(aom_env, mol, deltat)

    TYPE(aom_env_type)                         :: aom_env
    INTEGER                                    :: mol, i, j, sub, add
    REAL(KIND=dp)                              :: deltat

    CHARACTER(len=*), PARAMETER :: routineN ='nace1', &
                routineP = moduleN//':'//routineN

    !> Based on Eq. 19 of paper: Spencer et al., J. Chem. Phys. 145, 064102 (2016)

    !> 1. Calculate the term  D  (first term in Eq. 19)
    !many allocations
    !use the matinvsqrt subroutine to find T = S^(-1/2)
    aom_env%adiab_nace_env%T1 = aom_env%Stmtm
    aom_env%adiab_nace_env%T2 = aom_env%Stptp
    call matinvsqrt(mol, .TRUE., aom_env%adiab_nace_env%T1)
    call matinvsqrt(mol, .TRUE., aom_env%adiab_nace_env%T2)
    !and do some straightforward first-order numerical derivatives
    !Time deriv of overlap
    aom_env%adiab_nace_env%DT = (0.5D0/deltat)*(aom_env%adiab_nace_env%T2 - aom_env%adiab_nace_env%T1)
    !Time deriv of Lowdin Mat
    aom_env%adiab_nace_env%DS = (0.5D0/deltat)*(aom_env%Sttp - aom_env%Sttm)

!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "TST Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(6, "(A2, e12.4, A2)", advance='no') "| ", aom_env%adiab_nace_env%TST(i,j), " "
!   END DO
!   PRINT *, " "
!END DO

    !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
    !by transforming from the diabatic overlaps
    !TST = T1*Stt*DT. This is done in 2 steps.
    !Note that T2 is overwritten in the following DGEMM operation 
    !1. T2 = Stt*DT
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
                aom_env%Stt,               mol, &
                aom_env%adiab_nace_env%DT, mol, &
                0.0_DP, aom_env%adiab_nace_env%T2, mol)

    !2. TST = T1^{dagg}*T2
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                aom_env%adiab_nace_env%T1,  mol, &
                aom_env%adiab_nace_env%T2, mol, &       ! at this point diab_D = TST
                0.0_DP, aom_env%adiab_nace_env%diab_D, mol)

    ! This part takes advantage of the sparse nature of the DS matrix
    !IS ONLY VALID FOR CHAINS WITH ONLY OFF DIAGONALS ELEMENTS
    ! Again the mat-mat mult is done in 2 steps:
    !1. T2 = DS*T1
   DO j=1, mol
     aom_env%adiab_nace_env%T2(1,j) = &
      aom_env%adiab_nace_env%DS(1,2)*aom_env%adiab_nace_env%T1(2,j)
   END DO

   DO i=2, mol-1
      DO j=1, mol
        sub = i-1
        add = i+1
        aom_env%adiab_nace_env%T2(i,j) = &
          aom_env%adiab_nace_env%DS(i,sub)*aom_env%adiab_nace_env%T1(sub,j) &
        + aom_env%adiab_nace_env%DS(i,add)*aom_env%adiab_nace_env%T1(add,j)
      END DO
   END DO
   sub = mol-1
   DO j=1, mol
     aom_env%adiab_nace_env%T2(mol,j) = &
      aom_env%adiab_nace_env%DS(mol,sub)*aom_env%adiab_nace_env%T1(sub,j)
   END DO


!   ! If the loop above doesn't work then try using this.
!    !1. T2 = DS*T1
!    CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
!               aom_env%adiab_nace_env%DS, mol, &
!               aom_env%adiab_nace_env%T1, mol, &
!               0.0_DP, aom_env%adiab_nace_env%T2, mol)

    !2. TDT = T1^{dagg}*T2 + TST
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                aom_env%adiab_nace_env%T1,  mol, &
                aom_env%adiab_nace_env%T2, mol, &       ! Now diab_D becomes = TDT + TST
                1.0_DP, aom_env%adiab_nace_env%diab_D, mol)
!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
PRINT *, "D Bef:"
DO i=1, mol
   DO j=1, mol
      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", aom_env%adiab_nace_env%diab_D(i,j), "  "
   END DO
   PRINT *, " "
END DO
end subroutine nace1


!**********************************************************************************
!> \brief Calculate the diabatic NACE first
!> NB: this subroutine calculates diabatic NACE (d_kl) in the othogonal basis and in
!> practice it does exactly the same job as makeX subroutine with the exception
!> that does not calculate the term H_kl - d_kl
!**********************************************************************************
subroutine nace1_modified(aom_env, mol, deltat)

    TYPE(aom_env_type)                         :: aom_env
    INTEGER                                    :: mol!, i, j!, sub, add
    REAL(KIND=dp)                              :: deltat

    CHARACTER(len=*), PARAMETER :: routineN ='nace1', &
                routineP = moduleN//':'//routineN
!TEST
!    REAL(KIND=dp), dimension(:,:), allocatable     :: my_NACE

    !> Based on Eq. 19 of paper: Spencer et al., J. Chem. Phys. 145, 064102 (2016)

    !> 1. Calculate the term  D  (first term in Eq. 19)
    !many allocations
    !use the matinvsqrt subroutine to find T = S^(-1/2)
    aom_env%adiab_nace_env%T1 = aom_env%Stmtm
    aom_env%adiab_nace_env%T2 = aom_env%Stptp
    call matinvsqrt(mol, .TRUE., aom_env%adiab_nace_env%T1)
    call matinvsqrt(mol, .TRUE., aom_env%adiab_nace_env%T2)
    !and do some straightforward first-order numerical derivatives

    aom_env%adiab_nace_env%Ttt = aom_env%Stt
    call matinvsqrt(mol, .TRUE., aom_env%adiab_nace_env%Ttt)
    aom_env%adiab_nace_env%Sttsqrt_var = aom_env%Stt
    call matinvsqrt(mol, .FALSE., aom_env%adiab_nace_env%Sttsqrt_var)

    !Time deriv of overlap Tdot
    aom_env%adiab_nace_env%DT = (0.5D0/deltat)*(aom_env%adiab_nace_env%T2 - aom_env%adiab_nace_env%T1)
    !Time deriv of Lowdin Mat Dprime
    aom_env%adiab_nace_env%DS = (0.5D0/deltat)*(aom_env%Sttp - aom_env%Sttm)

!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "TST Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(6, "(A2, e12.4, A2)", advance='no') "| ", aom_env%adiab_nace_env%TST(i,j), " "
!   END DO
!   PRINT *, " "
!END DO

!   !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
!   !by transforming from the diabatic overlaps
!   !TST = T1*Stt*DT. This is done in 2 steps.
!   !Note that T2 is now overwritten in the following DGEMM operation 
!   !1. T2 = Stt*DT
!   CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
!               aom_env%Stt,               mol, &
!               aom_env%adiab_nace_env%DT, mol, &
!               0.0_DP, aom_env%adiab_nace_env%T2, mol)
!
!   !2. TST = T1^{dagg}*T2
!   CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
!               aom_env%adiab_nace_env%Ttt,  mol, &
!               aom_env%adiab_nace_env%T2, mol, &       ! at this point diab_D = TST
!               0.0_DP, aom_env%adiab_nace_env%diab_D, mol)

    !alternative to the double step above
    !2. TST = S^(1/2)*DT
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
                aom_env%adiab_nace_env%Sttsqrt_var, mol, &
                aom_env%adiab_nace_env%DT, mol, &
                0.0_DP, aom_env%adiab_nace_env%diab_D, mol)

    ! This part takes advantage of the sparse nature of the DS matrix
    ! THIS PART WORKS JUST FOR CHAINS AND IT'S NOT GENERAL
    ! Again the mat-mat mult is done in 2 steps:
    !1. T2 = DS*T1
!   DO j=1, mol
!     aom_env%adiab_nace_env%T2(1,j) = &
!      aom_env%adiab_nace_env%DS(1,2)*aom_env%adiab_nace_env%T1(2,j)
!   END DO
!   DO i=2, mol-1
!      DO j=1, mol
!        sub = i-1
!        add = i+1
!        aom_env%adiab_nace_env%T2(i,j) = &
!          aom_env%adiab_nace_env%DS(i,sub)*aom_env%adiab_nace_env%T1(sub,j) &
!        + aom_env%adiab_nace_env%DS(i,add)*aom_env%adiab_nace_env%T1(add,j)
!      END DO
!   END DO
!   sub = mol-1
!   DO j=1, mol
!     aom_env%adiab_nace_env%T2(mol,j) = &
!      aom_env%adiab_nace_env%DS(mol,sub)*aom_env%adiab_nace_env%T1(sub,j)
!   END DO


   ! If the loop above doesn't work then try using this.
    !1. T2 = DS*T1
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
               aom_env%adiab_nace_env%DS, mol, &
               aom_env%adiab_nace_env%Ttt, mol, &
               0.0_DP, aom_env%adiab_nace_env%T2, mol)


    !2. TDT = T1^{dagg}*T2 + TST
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                aom_env%adiab_nace_env%Ttt,  mol, &
                aom_env%adiab_nace_env%T2, mol, &       ! Now diab_D becomes = TDT + TST
                1.0_DP, aom_env%adiab_nace_env%diab_D, mol)
!!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "D Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", aom_env%adiab_nace_env%diab_D(i,j), "  "
!   END DO
!   PRINT *, " "
!END DO

!NB: becuase of the way finite differences are defined actually here
!aom_env%diabatic_nace_t refers at time T (namely t-dt) not Tp
   aom_env%diabatic_nace_t = (aom_env%adiab_nace_env%diab_D - & 
                              TRANSPOSE(aom_env%adiab_nace_env%diab_D)) / 2

!aom_env%adiab_nace_env%diab_D = aom_env%diabatic_nace !this would
!antisimmetrize NACE even in the calculation of adiab nace, 
!We have shown that it makes no difference at all!

end subroutine nace1_modified


!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 2.
!**********************************************************************************
subroutine nace2(aom_env, mol, print_here)

        TYPE(aom_env_type)                         :: aom_env
        INTEGER                                    :: mol
        LOGICAL                                    :: print_here

        CHARACTER(len=*), PARAMETER :: routineN ='nace2', &
                    routineP = moduleN//':'//routineN

        INTEGER                                    :: i, j, k, new_index
        REAL(KIND=dp)                              :: max_val, threshold
        LOGICAL                                    :: go_on

        !> 2. Check reordering of eigenvectors between the different timestep.
        aom_env%adiab_nace_env%map_new2old = 0
        aom_env%adiab_nace_env%map_old2new = 0
        aom_env%adiab_nace_env%not_associated = 0
        aom_env%adiab_nace_env%zeros = 0
        aom_env%adiab_nace_env%basic = (/( i, i=1, mol )/)
        threshold = 0.1
        !> First, we associate the new state i with the old state k if overlap is above
        !> (1 - threshold)
        do i=1, mol !loop over new index
           max_val = 0.0
           do j= 1, mol !loop over old index
              aom_env%adiab_overlap(i,j) = 0.0
              do k = 1, mol
                 aom_env%adiab_overlap(i,j) = aom_env%adiab_overlap(i,j) + aom_env%C_tp(k,i)*aom_Env%C_t(k,j)
              enddo
              if (abs(aom_env%adiab_overlap(i,j)).GE.abs(max_val)) THEN
                 max_val = aom_env%adiab_overlap(i,j)
                 new_index = j
              endif
           enddo
           if ( abs( max_val**2 - 1 ).LT.threshold) then
              aom_env%adiab_nace_env%map_new2old(i) = new_index
              aom_env%adiab_nace_env%map_old2new( aom_env%adiab_nace_env%map_new2old(i) ) = i
           endif
           !print_here = .TRUE.
           if (print_here) then
              if (i.NE.new_index) then
                 if ( abs( max_val**2 - 1 ).LT.threshold) then
                    print*, "FOBSH| Reordering: for state", i, "old index was", new_index, &
                            "with overlap", max_val, "and former", aom_env%adiab_overlap(i,i)
                 else
                    print*, "FOBSH! Reorgering: for state", i, " we stay in with overlap", &
                            aom_env%adiab_overlap(i,i), " instead of goind to ", new_index, "with overlap", max_val
                 endif
              endif
           endif
        enddo
        IF (ANY(aom_env%adiab_nace_env%map_new2old.NE.aom_env%adiab_nace_env%basic) .and. print_here) THEN
            print*, "FOBSH| New map: ", aom_env%adiab_nace_env%map_new2old
        END IF

        !> If we've associated state i with state j, we try to associate j with if if the
        !> overlap is above 0.5
        IF (ANY(aom_env%adiab_nace_env%map_new2old.EQ.aom_env%adiab_nace_env%zeros)) THEN
         do i=1,mol
            IF ((aom_env%adiab_nace_env%map_new2old(i).EQ.0).AND.(aom_env%adiab_nace_env%map_old2new(i).NE.0)) THEN
               IF ( abs( (aom_env%adiab_overlap(i, aom_env%adiab_nace_env%map_old2new(i)))**2 - 1 ).LT.0.5 ) THEN
                  aom_env%adiab_nace_env%map_new2old(i) = aom_env%adiab_nace_env%map_old2new(i)
                  aom_env%adiab_nace_env%map_old2new( aom_env%adiab_nace_env%map_new2old(i) ) = i
               ENDIF
            ENDIF
         enddo
        IF (ANY(aom_env%adiab_nace_env%map_new2old.NE.aom_env%adiab_nace_env%basic) .and. print_here) THEN
            print*, "FOBSH! New map again:", aom_env%adiab_nace_env%map_new2old
        END IF
       ENDIF

        !> Finally, we associate the not associate j with the remaining j state by order
        !> of energy
        IF (ANY(aom_env%adiab_nace_env%map_new2old.EQ.aom_env%adiab_nace_env%zeros)) THEN
         k = 1
         do i=1, mol
               j = 1
               go_on = .TRUE.
               do while ((go_on).AND.(j.LE.mol))
                  if (i.eq.aom_env%adiab_nace_env%map_new2old(j)) then
                     go_on = .FALSE.
                  else
                     j = j + 1
                  endif
               enddo
               if (go_on) then
                  aom_env%adiab_nace_env%not_associated(k) = i
                  k = k +1
               endif
         enddo
        j = 1
        do i=1, mol
           if (aom_env%adiab_nace_env%map_new2old(i).EQ.0) then
              aom_env%adiab_nace_env%map_new2old(i) = aom_env%adiab_nace_env%not_associated(j)
              aom_env%adiab_nace_env%map_old2new( aom_env%adiab_nace_env%map_new2old(i) ) = i
              j = j+1
           endif
        enddo
        IF (ANY(aom_env%adiab_nace_env%map_new2old.NE.aom_env%adiab_nace_env%basic) .and. print_here) THEN
           print*, "FOBSH| New map final:", aom_env%adiab_nace_env%map_new2old
           print*, "FOBSH| New map final, from old to new:", aom_env%adiab_nace_env%map_old2new
        ENDIF
       ENDIF
       IF (SUM(aom_env%adiab_nace_env%map_new2old).NE.( mol*(mol+1)/2 )) THEN
           CPABORT("Problem in the mapping between old and new states.")
       ENDIF
!> End of 2.

end subroutine nace2

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 3.
!**********************************************************************************
subroutine nace3(aom_env)

        TYPE(aom_env_type)                         :: aom_env

        CHARACTER(len=*), PARAMETER :: routineN ='nace3', &
                    routineP = moduleN//':'//routineN

        ! 3. Decide map
        IF (.NOT.(aom_env%do_reorder)) THEN
           aom_env%adiab_nace_env%map_new2old = aom_env%adiab_nace_env%basic
           aom_env%adiab_nace_env%map_old2new = aom_env%adiab_nace_env%basic
        ENDIF
        ! End of 3        

end subroutine nace3

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 4. This
!         does correct the phase as it would be arbitrary otherwise.
!**********************************************************************************
subroutine nace4(aom_env, mol)

        TYPE(aom_env_type)                         :: aom_env
        INTEGER                                    :: mol
        REAL(KIND=dp)                              :: sum_

        CHARACTER(len=*), PARAMETER :: routineN ='nace4', &
                    routineP = moduleN//':'//routineN

        LOGICAL                                    :: do_sign_change
        INTEGER                                    :: i, j

        ! 4. Change sign to keep constant convention
        do_sign_change = .True.
        IF (do_sign_change) THEN
           DO i=1, mol
              !PRINT*, aom_env%adiab_overlap(i, aom_env%adiab_nace_env%map_new2old(i))
              if (aom_env%adiab_overlap(i, aom_env%adiab_nace_env%map_new2old(i)).lt.0) then
                 aom_env%C_tp(:,i) = - aom_env%C_tp(:, i)
                 aom_env%adiab_overlap(i,:) = - aom_env%adiab_overlap(i,:)
              endif
              sum_ = 0  
              do j=1, mol
                 sum_ = sum_ + aom_env%C_tp(j, i) * aom_env%C_t(j, i)
              enddo
              !print*, "AC_PHASE", i, sum_
              !print*, "aom_env%adiab_overlap", aom_env%adiab_overlap(i, :)
              !if (i.eq.aom_env%active_state) print*, "AC_STATE_WARNING"          
           ENDDO
        ENDIF
        ! End of 4
end subroutine nace4


!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 5.
!**********************************************************************************
subroutine nace5(aom_env)
   
        TYPE(aom_env_type)                         :: aom_env

        CHARACTER(len=*), PARAMETER :: routineN ='nace5', &
                    routineP = moduleN//':'//routineN

        !> 5. Reassociate state active_state
        aom_env%active_state = aom_env%adiab_nace_env%map_old2new(aom_env%old_state)
        print*, "FOBSH| Put former active state", aom_env%old_state, "in new active state", aom_env%active_state 
        ! End of 5

end subroutine nace5

!**********************************************************************************
!> \brief Transform the diabatic NACE and apply the map.
!> Note that the diab_D matrix refers to time T whereas the U matrices refer to
!time Tp instead (this incosistency needs to be investigated further even though
!it will probably not afftect the results as UDU is a nagligible term)
!**********************************************************************************
subroutine nace6(aom_env, mol, deltat)

    TYPE(aom_env_type)                         :: aom_env
    INTEGER                                    :: mol
    REAL(KIND=dp)                              :: deltat, inv_deltat
 
    CHARACTER(len=*), PARAMETER :: routineN ='nace6', &
                routineP = moduleN//':'//routineN
 
    INTEGER                                    :: i,j
 
    inv_deltat = 1.0/deltat
    ! 6. Calculate the full NACE, Eq. (19) in paper 1.
    DO i=1, mol
       DO j=1, mol
          aom_env%adiab_nace_tp(i,j) = - aom_env%adiab_overlap(i, aom_env%adiab_nace_env%map_new2old(j))
          IF (i ==j) THEN 
             aom_env%adiab_nace_tp(i,j) = aom_env%adiab_nace_tp(i,j) + 1
          ENDIF
       ENDDO
    ENDDO
    ! End of 6

   if (aom_env%no_diab_nace) then 
     ! the adiabatic nace are just d = U^{dagger}Udot
     aom_env%adiab_nace_tp(:,:) = inv_deltat*aom_env%adiab_nace_tp(:,:)
   else
     !NB read description about time inconsistency
     !1. T1 <= D*U
     CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
          aom_env%adiab_nace_env%diab_D, mol, &    
          aom_env%C_tp,                  mol, &     ! C_tp = U
          0.0_DP, aom_env%adiab_nace_env%T1, mol) 

      !U^{dagger}U^dot = (1/dt)*adiab_nace_tp
      !2. adiab_nace_tp <= U^{dagg}*T1 + (1/dt)*adiab_nace_tp
      CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
          aom_env%C_tp,                  mol, &     ! C_tp = U
          aom_env%adiab_nace_env%T1,   mol, &    
          inv_deltat, aom_env%adiab_nace_tp, mol) 
   end if 
     
end subroutine nace6

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 7.
!**********************************************************************************
subroutine nace7(aom_env, mol) 
   
        TYPE(aom_env_type)                         :: aom_env
        INTEGER                                    :: mol

        CHARACTER(len=*), PARAMETER :: routineN ='nace7', &
                    routineP = moduleN//':'//routineN

        INTEGER                                    :: i,j

        ! 7. If reordering, we use the anti-symmetry of the NACE
        IF (aom_env%do_reorder) THEN 
           Do i=1, mol
              Do j=1, mol  
                 !print*, i,j, aom_env%adiab_nace_tp(i,j), aom_env%adiab_nace_tp(j,i)
                 !print*, i,j, ( aom_env%adiab_nace_tp(i,j) + aom_env%adiab_nace_tp(j,i) )
                 IF (j.lt.i) THEN 
                     aom_env%adiab_nace_tp(i,j) = - aom_env%adiab_nace_tp(j,i)
                 ELSE
                     aom_env%adiab_nace_tp(i,j) = ( aom_env%adiab_nace_tp(i,j) - aom_env%adiab_nace_tp(j,i) ) / 2
                 ENDIF
              Enddo
           Enddo
           !Do i=1, mol
           !   Do j=1, mol 
           !      print*, i,j, aom_env%adiab_nace_tp(i,j), aom_env%adiab_nace_tp(j,i)
           !   Enddo
           !Enddo
        ENDIF
        ! End of 7

end subroutine nace7

!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine calc_adiab_nace(aom_env, adiab_section)
        !this subroutine calculates the adiabatic nace used in the calculation
        !of the hopping probability and also in the RK adiabatic propagation 
        TYPE(aom_env_type)                         :: aom_env
        TYPE(section_vals_type), POINTER           :: adiab_section

        CHARACTER(len=*), PARAMETER :: routineN ='calc_adiab_nace', &
                                       routineP = moduleN//':'//routineN

        INTEGER                                    :: mol
        LOGICAL                                    :: print_here

        mol = size(aom_env%Stt,1)

        CALL section_vals_val_get(adiab_section,"PRINT_MORE", &
                                                   l_val=print_here)

        if (.NOT.aom_env%no_diab_nace) then
            CALL nace1_modified(aom_env, mol, aom_env%deltat)
        end if
        !CALL nace1(aom_env, mol, aom_env%deltat)
        CALL nace2(aom_env, mol, print_here)
        CALL nace3(aom_env)
        CALL nace4(aom_env, mol)
        CALL nace5(aom_env)
        CALL nace6(aom_env, mol, aom_env%deltat)
        CALL nace7(aom_env, mol)
 
end subroutine calc_adiab_nace

!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine matinvsqrt(mol, invert, S)

        !calculates either the square root or the inverse of the square root of a square matrix S
        !depending on whether invert is false or true respectively
        !to do this: diagonalise matrix, root eigenvalues, transform back

        integer, intent(in) :: mol
        logical, intent(in) :: invert
        REAL(KIND=dp), dimension(:,:), intent(inout) :: S
       
        CHARACTER(len=*), PARAMETER :: routineN ='matinvsqrt', &
                               routineP = moduleN//':'//routineN

        integer :: i, j, stat
        REAL(KIND=dp), dimension(:,:), allocatable :: P, Pinv
        REAL(KIND=dp), dimension(:), allocatable :: D, DGETRI_WORK
        character(len=1) :: job, up

        allocate(D(mol))
        allocate(DGETRI_WORK(3*mol))
        allocate(P(mol,mol))
        allocate(Pinv(mol,mol))

        !to calculate inverse square root of matrix, first diagonalise it
        job = 'V'
        up = 'U'
        !DSYEV is from the BLAS. Returns vector of eigenvalues D and overwrites S with matrix of eigenvectors
        call DSYEV(job, up, mol, S, mol, D, DGETRI_WORK, 3*mol, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if

        P = S
        !multiply P by diagonal matrix D^(-1/2)
        do i=1,mol
                if (invert) then
                        !D -> D^(-1/2)
                        D(i) = 1.0D0/sqrt(D(i))
                else
                        !D -> D^(1/2)
                        D(i) = sqrt(D(i))
                end if
                do j=1,mol
                        !Pinv = D*(P     )
                        Pinv(i,j) = D(i)*P(j,i)
                end do
        end do

        !call DGEMM from the BLAS to return S^(-1/2)
        !real symmetric matrix -> S = P*Pinv = P*D*(P     )
        call DGEMM('N', 'N', mol, mol, mol, 1.0D0, P, mol, Pinv, mol, 0.0D0, S, mol)

        deallocate(D, P, Pinv, DGETRI_WORK)
end subroutine matinvsqrt


! COMMENT_SG: this subroutine is still a bit tight up to AOM and should be
! easily generalized for other Hamiltonians 
!**********************************************************************************
!> \brief Calculate adiabatic forces according to Spencer 2016 and Carof 2019
!**********************************************************************************
subroutine calcforces(state, Hin, C, nacv, atoms, scaling, cbar, Fin, forces_adiab, & 
                  ext_pres)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom
        integer, intent(in)                         :: state
        REAL(KIND=dp), intent(in), dimension(:,:)   :: Hin, C
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: nacv
        integer, intent(in)                         :: atoms
        REAL(KIND=dp), intent(in)                   :: scaling, cbar
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), dimension(:,:), POINTER      :: forces_adiab
        LOGICAL, dimension(:), INTENT(IN)           :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN ='calcforces', &
                               routineP = moduleN//':'//routineN

        integer                                     :: i, j, mat, d
        INTEGER                                     :: natom
        INTEGER                                     :: num_mol, nadiab, num_j
        REAL(KIND=dp)                               :: convert_nacv 

        natom = SIZE(Fin,2)
        nadiab= SIZE(Fin,3)
        mat = size(Hin,1)

        ! Based on paper [Spencer16]: Spencer et al., J. Chem. Phys. 145, 064102 (2016)
        ! Use Eq. (21) of paper [Spencer16]
        convert_nacv = - scaling / cbar ! to convert the FOB NACV
        num_mol = 0                     ! molecular index
        do d=1,3                        ! sum of cartesian 
                num_j = 0               ! active atom index
                do j=1,natom            ! sum over all atoms
                   forces_adiab(d,j) = 0.0D0
                   IF (ext_pres(j)) THEN   
                      num_j   = num_j + 1                 ! update active atom index
                      num_mol = int((num_j-1)/atoms) + 1  ! update molecular index
                      ! WARNING: we suppose that all molecules have the same
                      ! number of atoms
                   ENDIF
                   do i=1,mat
                      ! Diagonal term, it depends on the forces
                      forces_adiab(d,j) = forces_adiab(d,j) + (abs(C(i,state))**2)*Fin(d,j,i)
                      ! Of diagonal term, only for active atom j present in
                      ! molecule num_mol
                      IF (ext_pres(j)) THEN
                         IF (i.ne.num_mol) THEN
                            forces_adiab(d,j) = forces_adiab(d,j) + &
                                        2.0D0*((C(num_mol,state))*nacv(d,i, num_j)*C(i,state))*convert_nacv
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
        ENDDO
end subroutine calcforces

!**********************************************************************************
!> \brief This reads the wavefunction coefficient section from the restart file, 
!> the whole restart procedure most probably not completely general and it would
!> be still to improve if a more general Hamiltonian is use. That's why this
!> particular subroutine is still in aom_utils.F
!**********************************************************************************
  SUBROUTINE read_listdiabcoeff(listdiabcoeff, norbitals, nadiab, ctin)
   TYPE(cp_sll_val_type), POINTER                         :: listdiabcoeff
   INTEGER, INTENT(IN)                                    :: norbitals, nadiab
   COMPLEX(KIND=dp), dimension(:), allocatable, INTENT(OUT):: ctin

   CHARACTER(len=*), PARAMETER :: routineN = 'read_listdiabcoeff', &
     routineP = moduleN//':'//routineN

   CHARACTER(len=2*default_string_length)   :: line_att
   CHARACTER(len=default_string_length)     :: error_message
   INTEGER                                  :: default_id, end_c, &
                                               iatom, j,  &
                                               start_c, wrd, mat
   LOGICAL                                  :: is_ok
   TYPE(val_type), POINTER                  :: val
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: realcoeff, imagcoeff
!   INTEGER                                  :: my_pres
!   REAL(KIND=dp)                            :: my_ext_elements, my_pres2

   NULLIFY(val)

   mat = norbitals*nadiab
!   ALLOCATE(ext_elements(natom))
!   ALLOCATE(ext_pres(natom))
   ALLOCATE(realcoeff(mat))
   ALLOCATE(imagcoeff(mat))
   ALLOCATE(ctin(mat))


      DO iatom=1,mat
         ! we use only the first default_string_length characters of each line
         is_ok=cp_sll_val_next(listdiabcoeff,val) !logical 
         CALL val_get(val,c_val=line_att)
         default_id = str2id(s2s(""))
         ! Read name and atomic coordinates
         start_c=1
         DO wrd=1,4  !loop over entries in DIAB_COEFF.inc
            DO j=start_c,LEN(line_att)
               IF (line_att(j:j)/=' ') THEN
                  start_c=j
                  EXIT
               END IF
            END DO
            end_c=LEN(line_att)+1
            DO j=start_c,LEN(line_att)
               IF (line_att(j:j)==' ') THEN
                  end_c=j
                  EXIT
               END IF
            END DO
            IF(LEN_TRIM(line_att(start_c:end_c-1)) == 0)& !length of the number
                                                           ! without space
               CALL cp_abort(__LOCATION__,&
                    "Incorrectly formatted input line for atom "//&
                    TRIM(ADJUSTL(cp_to_string(iatom)))//&
                    " found in COORD section. Input line: <"//&
                    TRIM(line_att)//"> ")
            SELECT CASE (wrd)
! We ignore case 1 and 2 since they report state and orbital per state
! this may be used in case of multiple orbitals per state
!!             CASE (1)
!!                 = str2id(s2s(line_att(start_c:end_c-1)))
!            CASE (2)
!               CALL read_float_object(line_att(start_c:end_c-1),&
!                                      my_ext_elements,error_message)
!              ext_elements(iatom) = INT(my_ext_elements)
!              IF(LEN_TRIM(error_message) /= 0)&
!                 CALL cp_abort(__LOCATION__,&
!                      "Incorrectly formatted input line for atom "//&
!                      TRIM(ADJUSTL(cp_to_string(iatom)))//&
!                      " found in COORD section. "//TRIM(error_message)//&
!                      " Input line: <"//TRIM(line_att)//"> ")
           CASE (3)
              CALL read_float_object(line_att(start_c:end_c-1),&
                                     realcoeff(iatom),error_message)
              IF(LEN_TRIM(error_message) /= 0)&
                 CALL cp_abort(__LOCATION__,&
                      "Incorrectly formatted input line for atom "//&
                      TRIM(ADJUSTL(cp_to_string(iatom)))//&
                      " found in COORD section. "//TRIM(error_message)//&
                      " Input line: <"//TRIM(line_att)//"> ")
           CASE (4)
              CALL read_float_object(line_att(start_c:end_c-1),&
                                     imagcoeff(iatom),error_message)
              IF(LEN_TRIM(error_message) /= 0)&
                 CALL cp_abort(__LOCATION__,&
                      "Incorrectly formatted input line for atom "//&
                      TRIM(ADJUSTL(cp_to_string(iatom)))//&
                      " found in COORD section. "//TRIM(error_message)//&
                      " Input line: <"//TRIM(line_att)//"> ")
           END SELECT
           start_c = end_c
            IF (start_c > LEN_TRIM(line_att)) EXIT
         END DO
print *, "real", realcoeff(iatom)
print *, "imag", imagcoeff(iatom)
!store diabatic coeff as a complex number
         ctin(iatom) = CMPLX(realcoeff(iatom), imagcoeff(iatom),dp)
     END DO
print *, "CTIN", ctin

     !DEALLOCATE(val)
  END SUBROUTINE read_listdiabcoeff

!**********************************************************************************
!> \brief Read the aom coefficients given as input file in AOM.include
!**********************************************************************************
   SUBROUTINE read_listcoeff(listcoeff, ext_elements, ext_pres, &
                            ext_scoeff, ext_pcoeff, natom)
    TYPE(cp_sll_val_type), POINTER                         :: listcoeff
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: ext_elements
    LOGICAL, DIMENSION(:), POINTER,  INTENT(OUT)           :: ext_pres
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: ext_scoeff, ext_pcoeff
    INTEGER, INTENT(IN)                                    :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'read_listcoeff', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)                 :: line_att
    CHARACTER(len=default_string_length)                   :: error_message
    INTEGER                                                :: default_id, end_c, &
                                                              iatom, j,  &
                                                              start_c, wrd
    LOGICAL                                                :: is_ok
    TYPE(val_type), POINTER                                :: val
    INTEGER                                                :: my_pres
    REAL(KIND=dp)                                          :: my_ext_elements, my_pres2

    NULLIFY(val)

    ALLOCATE(ext_elements(natom))
    ALLOCATE(ext_pres(natom))
    ALLOCATE(ext_scoeff(natom))
    ALLOCATE(ext_pcoeff(natom))

       DO iatom=1,natom
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(listcoeff,val) !logical 
          CALL val_get(val,c_val=line_att)
          default_id = str2id(s2s(""))
          ! Read name and atomic coordinates
          start_c=1
          DO wrd=1,5
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)/=' ') THEN
                   start_c=j
                   EXIT
                END IF
             END DO
             end_c=LEN(line_att)+1
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)==' ') THEN
                   end_c=j
                   EXIT
                END IF
             END DO
             IF(LEN_TRIM(line_att(start_c:end_c-1)) == 0)& !length of the number
                                                            ! without space
                CALL cp_abort(__LOCATION__,&
                     "Incorrectly formatted input line for atom "//&
                     TRIM(ADJUSTL(cp_to_string(iatom)))//&
                     " found in COORD section. Input line: <"//&
                     TRIM(line_att)//"> ")
             SELECT CASE (wrd)
             CASE (2)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_ext_elements,error_message)
                ext_elements(iatom) = INT(my_ext_elements)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (3)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_pres2,error_message)
                my_pres = INT(my_pres2)
                IF (my_pres.eq.0) THEN
                   ext_pres(iatom) = .FALSE.
                ELSE
                   ext_pres(iatom) = .TRUE.
                ENDIF
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (4)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       ext_scoeff(iatom),error_message)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (5)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       ext_pcoeff(iatom),error_message)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")

             END SELECT
             start_c = end_c
             IF (start_c > LEN_TRIM(line_att)) EXIT
          END DO
      END DO

   END SUBROUTINE read_listcoeff



!******************************************************************************
!!!!!!!!!!!!!!   DEPRECATED SUB. BELOW  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

subroutine calc_adiab_nace_old(aom_env, adiab_section)
        !this subroutine calculates the adiabatic nace used in the calculation
        !of the hopping probability and also in the RK adiabatic propagation 
!WARNING: THIS SUBROUTINE HAS TIME INCONSISTENCIES IN THE TIME DERIVATIVES THAT
!HAVE BEEN CORRECTED IN THE CHOPPED VERSION ABOVE
        TYPE(aom_env_type)                         :: aom_env
        TYPE(section_vals_type), POINTER           :: adiab_section

       CHARACTER(len=*), PARAMETER :: routineN ='calc_adiab_nace_old', &
       routineP = moduleN//':'//routineN

        integer :: i, j, mol, k, l, new_index
        REAL(KIND=dp), dimension(:,:), allocatable :: T1, T2, DT, TST, TDT, DS!, DC
        INTEGER, dimension(:),  allocatable        :: map_new2old, map_old2new, basic
        INTEGER, dimension(:),  allocatable        :: not_associated, zeros
        REAL(KIND=dp)                              :: deltat, max_val, threshold
        LOGICAL                                    :: go_on, print_here, do_sign_change
        REAL(kind=dp)  :: sum_

        mol = size(aom_env%Stt,1)
        deltat = aom_env%deltat

        CALL section_vals_val_get(adiab_section,"PRINT_MORE", &
                                                   l_val=print_here)

        !> Based on Eq. 19 of paper: Spencer et al., J. Chem. Phys. 145, 064102 (2016)

        !> 1. Calculate the term  D  (first term in Eq. 19)
        !many allocations
        allocate(T1(mol,mol))
        allocate(T2(mol,mol))
        allocate(DT(mol,mol))
        allocate(TST(mol,mol))
        allocate(TDT(mol,mol))
        allocate(DS(mol,mol))
        !use the matinvsqrt subroutine to find T = S^(-1/2)
        T1 = aom_env%Stmtm
        T2 = aom_env%Stptp
        call matinvsqrt(mol, .TRUE., T1)
        call matinvsqrt(mol, .TRUE., T2)
        !and do some straightforward first-order numerical derivatives
        DT = (0.5D0/deltat)*(T2 - T1)
        DS = (0.5D0/deltat)*(aom_env%Sttp - aom_env%Sttm)
        !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
        !by transforming from the diabatic overlaps
        do i=1,mol
                do j=1,mol
                        TST(i,j) = 0.0D0
                        TDT(i,j) = 0.0D0
                        do k=1,mol
                                do l=1,mol
                                        TST(i,j) = TST(i,j) + T1(k,i)*aom_env%Stt(k,l)*DT(l,j)
                                        TDT(i,j) = TDT(i,j) + T1(k,i)*DS(k,l)*T1(l,j)
                                end do
                        end do
                end do
        end do
        !> End of 1.

 
        !> 2. Check reordering of eigenvectors between the different timestep.
        ALLOCATE(map_new2old(mol))
        ALLOCATE(map_old2new(mol))
        ALLOCATE(basic(mol))
        ALLOCATE(zeros(mol))
        ALLOCATE(not_associated(mol))
        map_new2old = 0
        map_old2new = 0
        not_associated = 0
        zeros = 0
        basic = (/( i, i=1, mol )/)
        threshold = 0.1
        !> First, we associate the new state i with the old state k if overlap is above
        !> (1 - threshold)
        do i=1, mol !loop over new index
           max_val = 0.0
           do j= 1, mol !loop over old index
              aom_env%adiab_overlap(i,j) = 0.0
              do k = 1, mol
                 aom_env%adiab_overlap(i,j) = aom_env%adiab_overlap(i,j) + aom_env%C_tp(k,i)*aom_Env%C_t(k,j)
              enddo
              if (abs(aom_env%adiab_overlap(i,j)).GE.abs(max_val)) THEN
                 max_val = aom_env%adiab_overlap(i,j)
                 new_index = j
              endif
           enddo
           if ( abs( max_val**2 - 1 ).LT.threshold) then
              map_new2old(i) = new_index
              map_old2new( map_new2old(i) ) = i
           endif
           !print_here = .TRUE.
           if (print_here) then
              if (i.NE.new_index) then
                 if ( abs( max_val**2 - 1 ).LT.threshold) then
                    print*, "FOBSH| Reordering: for state", i, "old index was", new_index, &
                            "with overlap", max_val, "and former", aom_env%adiab_overlap(i,i)
                 else 
                    print*, "FOBSH! Reorgering: for state", i, " we stay in with overlap", &
                            aom_env%adiab_overlap(i,i), " instead of goind to ", new_index, "with overlap", max_val
                 endif
              endif
           endif
        enddo
        IF (ANY(map_new2old.NE.basic) .and. print_here) print*, "FOBSH| New map: ", map_new2old

        !> If we've associated state i with state j, we try to associate j with if if the
        !> overlap is above 0.5
        IF (ANY(map_new2old.EQ.zeros)) THEN
         do i=1,mol
            IF ((map_new2old(i).EQ.0).AND.(map_old2new(i).NE.0)) THEN
               IF ( abs( (aom_env%adiab_overlap(i, map_old2new(i)))**2 - 1 ).LT.0.5 ) THEN
                  map_new2old(i) = map_old2new(i)
                  map_old2new( map_new2old(i) ) = i
               ENDIF
            ENDIF
         enddo
        IF (ANY(map_new2old.NE.basic) .and. print_here) print*, "FOBSH! New map again:", map_new2old
       ENDIF

        !> Finally, we associate the not associate j with the remaining j state by order
        !> of energy
        IF (ANY(map_new2old.EQ.zeros)) THEN
         k = 1
         do i=1, mol
               j = 1
               go_on = .TRUE.
               do while ((go_on).AND.(j.LE.mol))
                  if (i.eq.map_new2old(j)) then
                     go_on = .FALSE.
                  else
                     j = j + 1
                  endif
               enddo
               if (go_on) then
                  not_associated(k) = i
                  k = k +1
               endif
         enddo
        j = 1
        do i=1, mol
           if (map_new2old(i).EQ.0) then
              map_new2old(i) = not_associated(j)
              map_old2new( map_new2old(i) ) = i
              j = j+1
           endif
        enddo 
        IF (ANY(map_new2old.NE.basic) .and. print_here) THEN 
           print*, "FOBSH| New map final:", map_new2old
           print*, "FOBSH| New map final, from old to new:", map_old2new
        ENDIF
       ENDIF
       IF (SUM(map_new2old).NE.( mol*(mol+1)/2 )) THEN
           CPABORT("Problem in the mapping between old and new states.")
       ENDIF
!> End of 2.


! 3. Decide map
       IF (.NOT.(aom_env%do_reorder)) THEN
          map_new2old = basic
          map_old2new = basic
       ENDIF
! End of 3        

! 4. Change sign to keep constant convention
       do_sign_change = .True.
       IF (do_sign_change) THEN
          DO i=1, mol
             !PRINT*, aom_env%adiab_overlap(i, map_new2old(i))
             if (aom_env%adiab_overlap(i, map_new2old(i)).lt.0) then
                aom_env%C_tp(:,i) = - aom_env%C_tp(:, i)
                aom_env%adiab_overlap(i,:) = - aom_env%adiab_overlap(i,:)
             endif
             sum_ = 0 
             do j=1, mol
                sum_ = sum_ + aom_env%C_tp(j, i) * aom_env%C_t(j, i)
             enddo
             !print*, "AC_PHASE", i, sum_
             !print*, "aom_env%adiab_overlap", aom_env%adiab_overlap(i, :)
             !if (i.eq.aom_env%active_state) print*, "AC_STATE_WARNING"          
          ENDDO
       ENDIF
! End of 4


!> 5. Reassociate state active_state
        aom_env%active_state = map_old2new(aom_env%old_state)
        print*, "FOBSH| Put former active state", aom_env%old_state, "in new active state", aom_env%active_state 
! End of 5

! 6. Calculate the full NACE, Eq. (19) in paper 1.
        DO i=1, mol
           DO j=1, mol
              aom_env%adiab_nace_tp(i,j) = - aom_env%adiab_overlap(i, map_new2old(j))
              IF (i ==j) THEN
                 aom_env%adiab_nace_tp(i,j) = aom_env%adiab_nace_tp(i,j) + 1
              ENDIF
              aom_env%adiab_nace_tp(i,j) = aom_env%adiab_nace_tp(i,j) / deltat
              DO k=1, mol
                 DO l=1, mol
                    aom_env%adiab_nace_tp(i,j) = & 
                                    aom_env%adiab_nace_tp(i,j) + aom_env%C_tp(k,i)*aom_env%C_tp(l,j)*(TST(k,l) + TDT(k,l))
                 ENDDO
              ENDDO  
           ENDDO
        ENDDO
! End of 6

! 7. If reordering, we use the anti-symmetry of the NACE
        IF (aom_env%do_reorder) THEN
           Do i=1, mol
              Do j=1, mol 
                 !print*, i,j, aom_env%adiab_nace_tp(i,j), aom_env%adiab_nace_tp(j,i)
                 !print*, i,j, ( aom_env%adiab_nace_tp(i,j) + aom_env%adiab_nace_tp(j,i) )
                 IF (j.lt.i) THEN
                     aom_env%adiab_nace_tp(i,j) = - aom_env%adiab_nace_tp(j,i)
                 ELSE
                     aom_env%adiab_nace_tp(i,j) = ( aom_env%adiab_nace_tp(i,j) - aom_env%adiab_nace_tp(j,i) ) / 2
                 ENDIF
              Enddo
           Enddo
           !Do i=1, mol
           !   Do j=1, mol 
           !      print*, i,j, aom_env%adiab_nace_tp(i,j), aom_env%adiab_nace_tp(j,i)
           !   Enddo
           !Enddo
        ENDIF
! End of 7

        deallocate(T1)
        deallocate(T2)
        deallocate(DT)
        deallocate(TST)
        deallocate(TDT)
        deallocate(DS)
        !deallocate(DC)
end subroutine calc_adiab_nace_old

subroutine surfacehop(deltat, old_state, C2, &
              new_state, attempt_hop, decohere_criterion, sh_print_env, &
                                                     aom_env, adiab_section)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface
        REAL(KIND=dp), intent(in) :: deltat
        integer, intent(in) :: old_state
        REAL(KIND=dp), dimension(:,:), intent(in) :: C2
        INTEGER, INTENT(OUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: attempt_hop
        REAL(KIND=dp), INTENT(IN) :: decohere_criterion
        TYPE(sh_print_env_type), INTENT(OUT) :: sh_print_env
        TYPE(aom_env_type)                   :: aom_env
        TYPE(section_vals_type), POINTER               :: adiab_section

        CHARACTER(len=*), PARAMETER :: routineN ='surfacehop', &
                   routineP = moduleN//':'//routineN

        CHARACTER(LEN=default_string_length)         :: wfpropag_keyword
        integer :: i, j, mol !, k, l
        REAL(KIND=dp) :: hop, sumg
        REAL(KIND=dp), dimension(:), allocatable :: b, g
        COMPLEX(KIND=dp), dimension(:), allocatable :: ct
        COMPLEX(KIND=dp), dimension(:), allocatable :: ak
        COMPLEX(KIND=dp) :: cstate
        
        mol = size(aom_env%Stt,1)

        !many allocations
        allocate(ak(mol))
        allocate(b(mol))
        allocate(g(mol))
        ALLOCATE(sh_print_env%proba(mol))
        allocate(ct(mol))

        CALL section_vals_val_get(adiab_section,"REPRESENTATION", &
                                                     c_val=wfpropag_keyword)

        SELECT CASE(wfpropag_keyword)
        CASE("DIABATIC_BASIS")
            CALL old_diab_to_adiab(aom_env, C2, old_state, mol, cstate, ct)
        CASE("ADIABATIC_BASIS")
            ct(:) = aom_env%ctin(:)
            cstate = aom_env%ctin_t(old_state)
        CASE DEFAULT
           CPABORT("There is no method to propagate the wavefunction.")
        END SELECT

!NOT YET READY FOR ADIAB PRPAGATION

        do i=1,mol
                ak(i) = ct(i)*conjg(cstate)
                if (aom_env%HS_propagation .and. &
                              (i .eq. aom_env%class_forbidden_hop_t(i))) then
                   b(i) = 0.0D0 
                else 
                b(i) = -2.0D0*REAL(conjg(ak(i))*aom_env%adiab_nace_tp(i,old_state))
                end if
        end do

        sumg = 0.0D0
        !calculate g and sum(g)
        do i=1,mol
                g(i) = (deltat*b(i))/abs(ak(old_state))
                !as in Tully, set negative probabilities to 0
                if (g(i) < 0.0D0) then
                        g(i) = 0.0D0
                end if
                if (i .ne. old_state) then
                        sumg = sumg + g(i)
                end if
        end do

        !due to finite timestep, sometimes the probabilities may sum to greater than 1
        !if this happens, rescale them down to sum = 1
        if (sumg > 1.0D0) then
                do i=1,mol
                        g(i) = g(i)/sumg
                        g(old_state) = 0.0D0
                end do
!> CHANGE_AC
                sh_print_env%renorm = .TRUE.
                sh_print_env%sumprob   = sumg        
        else
!> CHANGE_AC
                sh_print_env%renorm = .FALSE.
                sh_print_env%sumprob   = sumg        
                g(old_state) = 1.0D0 - sumg
        end if

!print *, "g_AFTER", g
!> CHANGE_AC
        sh_print_env%proba     = g

        hop = next_random_number(aom_env%rng_stream)

        sh_print_env%hop       = hop
        sh_print_env%old_state = old_state
        !and see if this random number causes a surface hop
        sumg = 0.0D0
        do i=1,mol
                sumg = sumg + g(i)
                if (hop <= sumg) then
                        if (i .ne. old_state) then
                                new_state = i
                               ! decohere = .FALSE.
                                aom_env%decohere = .FALSE.
                                attempt_hop = .TRUE.
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
                                EXIT
                        elseif (i == old_state) then
                                new_state = old_state
                                attempt_hop = .FALSE.
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
!                                write(24,*) state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                aom_env%decohere = .TRUE.
                                do j=1,mol
                                        if (abs(aom_env%adiab_nace_tp(old_state,j)) & 
                                                 > decohere_criterion .and. old_state .ne. j) then
                                               ! decohere = .FALSE.
                                                aom_env%decohere = .FALSE.
                                                EXIT
                                        end if
                                end do
                               ! decohere = .TRUE.
                                EXIT
                        end if
                end if
        end do


        CALL adiab_print_sh_info(sh_print_env, keyword = 'TULLY_TEST') 

        DEALLOCATE(sh_print_env%proba)
        deallocate(ak)
        deallocate(b)
        deallocate(G)
end subroutine surfacehop

!this subroutine might be redundant with the next one
SUBROUTINE rotation_diab_to_adiab(C_tp, ctin, mat, adiab_ctin, pop)

        REAL(KIND=dp), dimension(:,:), intent(in) :: C_tp
        COMPLEX(KIND=dp), dimension(:), intent(in) :: ctin
        integer, intent(in) :: mat
        COMPLEX(KIND=dp), dimension(:), allocatable, intent(out) :: adiab_ctin
        REAL(KIND=dp), dimension(:), allocatable, intent(out)    :: pop
         
        integer :: i, j!, mol!, k, l

        ALLOCATE(adiab_ctin(mat))
        ALLOCATE(pop(mat))

        pop(:) = 0.0D0
        !adiab_ctin =  cmplx(0.0D0, 0.0D0, dp)
        do i=1,mat
            adiab_ctin(i) = cmplx(0.0D0, 0.0D0,dp)
            do j=1,mat
               adiab_ctin(i) = adiab_ctin(i) + cmplx(C_tp(j,i), 0.0D0,dp)*ctin(j)
            end do
            pop(i) = real(adiab_ctin(i)*conjg(adiab_ctin(i)))
        enddo 
        

END SUBROUTINE rotation_diab_to_adiab

SUBROUTINE old_adiab_overlap_check(aom_env, changed_state, itimes)
    TYPE(aom_env_type)                         :: aom_env
    INTEGER, OPTIONAL, intent(in)               :: changed_state, itimes  

    CHARACTER(len=*), PARAMETER :: routineN = 'old_adiab_overlap_check', &
              routineP = moduleN//':'//routineN

    REAL(KIND=dp), dimension(:), allocatable    :: overlaps
    integer                                     :: mat, i, j, res_mod, state, &
                                                     past_state
    REAL(KIND=dp)                               :: sum_ov, max_val, resid 
    LOGICAL :: overlap_over_thresh
!!ADD OPTIONAL

    mat    = SIZE(aom_env%Stptp, 1)

    allocate(overlaps(2*mat))

    overlaps = 0.0D0
    
    
    if (PRESENT(changed_state)) then
        past_state = changed_state
!print *, "OLD_STATE", changed_state
    else
        past_state = aom_env%active_state
!print *, "STATE_CHANGED", past_state
!print *, "NEW_STATE_here", aom_env%active_state    
    end if

!print *, "NEW_STATE", aom_env%active_state    

     do i=1, mat
         overlaps(1) = overlaps(1) + aom_env%C_t(i, past_state)*aom_env%C_tp(i, aom_env%active_state)
     end do
!print *, "OVERLAP_NOTSQUARED", overlaps(1)

     overlaps(1) = abs(overlaps(1))*abs(overlaps(1))
print *, "OVERALP_SQUARED", past_state, aom_env%active_state, overlaps(1)

     overlap_over_thresh = .False.
  
     if ( overlaps(1) .gt. 0.5) then
!        print *, "CARRY_ON_THE_SAME_STATE"
         overlap_over_thresh = .True.
     else
!        print *, "POSSIBLE_TRIVIAL_CROSSING"
        sum_ov = overlaps(1)
     !   max_val = maxval(overlaps(:))
        do j =2, 2*mat  !I want to be sure to count all the states
           res_mod = ((-1)**(modulo(j ,2)))*(j/2)   !this is +1,-1,+2,-2
           state = aom_env%active_state+res_mod
           !two conditions in case state is negative or is larger than mat
           if (state .le. 0 .or. state .gt. mat) then
!              print *, "skipped", j, state
              cycle
           else  
              do i = 1, mat
                 overlaps(j) = overlaps(j) + aom_env%C_t(i, past_state)*&
                                             aom_env%C_tp(i, state)
              end do
              overlaps(j) = abs(overlaps(j))*abs(overlaps(j))
!              print *, j, state
           end if
!           print *, "OV_CALC_NOW", overlaps(j)
           max_val = maxval(overlaps(:))
!           print *, "MAX_VAL_UPTONOW", max_val 
           sum_ov = sum_ov +  overlaps(j)
!           print *, "SUM_OV", sum_ov
           resid = 1.0D0 - sum_ov
           if (resid .gt. max_val) then 
!               print *, "REPEAT"
               cycle
           else
!               print *, "MAX FOUND", max_val
               if (max_val .le. 0.5) then 
                  print *, "WARNING: MORE THAN TWO STATES CROSS"
               end if 
               exit
           end if
        end do
     end if

!print *, "OVERALP_CALCC", overlaps
if (present(itimes)) then 
print *, "MAX_OVERLAP", itimes, overlaps(1), maxval(overlaps(:))
else
print *, "MAX_OVERLAP_HOPPING", overlaps(1), maxval(overlaps(:))
end if

if ( overlaps(1) .gt. 0.5) then
print *, "LOCALIZATION", maxloc(abs(aom_env%C_t(:, past_state)*&
                            aom_env%C_t(:, past_state))), &
      maxloc(abs(aom_env%C_tp(:,aom_env%active_state)*aom_env%C_tp(:,aom_env%active_state)))
print *, "STATE_TOGOTO", past_state, aom_env%active_state
else
print *, "LOCALIZATION", maxloc(abs(aom_env%C_t(:, past_state)*&
                            aom_env%C_t(:, past_state)))
!            maxloc(abs(aom_env%C_tp(:,aom_env%enforced_active_state)*&
!                            aom_env%C_tp(:,aom_env%enforced_active_state)))
!print *, "STATE_TOGOTO", past_state, aom_env%enforced_active_state
end if

    deallocate(overlaps)
     
END SUBROUTINE old_adiab_overlap_check

END MODULE aom_utils
