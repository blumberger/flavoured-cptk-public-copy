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
! propagate_namd_all_steps_full_degen
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
   USE sh_types,                      ONLY: sh_print_env_type, &
                                             sh_env_type
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

 ! USE non_adiabatic_types,             ONLY: sh_env_type

  USE parallel_rng_types,              ONLY: next_random_number 

!  USE sh_types,                       ONLY: sh_env_type

  USE aom_overlapfinal,               ONLY:  connect_list2, &
                                              calc_sab, &
                                              sab_screen_new, &
                                              calc_dRSab, &
                                              calc_rvecs, &
                                              calc_sab_pete
!  USE non_adiabatic_utils,             ONLY: simple_vel_rescaling, &
!                                             cal_delta_energy_nacv,&
!                                             diab_to_adiab

#include "./base/base_uses.f90"




IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_utils'

   PUBLIC :: makeX,                                &
             calc_hop_nacv_aom,                    &
             adiabaticity,                         &
             calc_multiphase_Xmatrix_elements,     &
             calc_aom_Xmatrix_elements,            & 
             calc_aom_Xmatrix_elements_first_step, & 
             make_adiabX,                          &
             read_listcoeff,                       & 
             calc_adiab_nace,                      &
             calc_diabatic_nace,                   &
             read_listdiabcoeff,                   & 
             read_list_tresp,                      &
             read_list_frz_couplings,              &
             calc_adiab_nace_old,                  &
             matinvsqrt,                           &
             nace1_time_t                   
 

CONTAINS

!**********************************************************************************
!> \brief Diagonalize the Hamiltonian and calculate eigevlaues and eigenvectors.
!>        Although this subroutine is not strictly related to AOM (but is is
!>        general, it has to stay here to avoid compilation issue 
!>        non_adiabatic_utils.F -> aom_utils.F -> non_adiabatic_utils.F
!> \author AC [Edited SG]
!**********************************************************************************
subroutine adiabaticity(H, C, E)
        !this subroutine will allow calculation of surface hopping probabilities
        !between N adiabatic surfaces

        REAL(KIND=dp), dimension(:,:), intent(in) :: H
        REAL(KIND=dp), dimension(:,:), intent(out), allocatable :: C
        REAL(KIND=dp), dimension(:), intent(out), allocatable :: E
        !logical, intent(in) :: store_C
        !TYPE(sh_env_type), INTENT(INOUT), OPTIONAL               ::    sh_env

        CHARACTER(len=*), PARAMETER :: routineN ='adiabaticity', &
                               routineP = moduleN//':'//routineN

        integer :: mol
        integer :: stat, LWORK
        character(len=1) :: job, up
        REAL(KIND=dp), dimension(:), allocatable :: DGETRI_WORK

        mol = SIZE(H, 1) !this is actually mat because H could be degenerate
 
        allocate(C(mol,mol))
        allocate(E(mol))
        allocate(DGETRI_WORK(3*size(C)-1))
 
        C = H

        job = 'V'
        up = 'U'
        !DSYEV is from the BLAS. Returns vector of eigenvalues D and overwrites S with matrix of eigenvectors
        !Original call establishes optimum size for speed purposes
        LWORK = -1
        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if
        LWORK = max(3*size(C)-1, int(DGETRI_WORK(1)))

        deallocate(DGETRI_WORK)
        allocate(DGETRI_WORK(LWORK))

        call DSYEV(job, up, mol, C, mol, E, DGETRI_WORK, LWORK, stat)
        if (stat /= 0) then
                CPABORT("Error with subroutine DSYEV: error status")
        end if
        
        !if (AND.PRESENT(sh_env)) then 
        !allocate(sh_env%C_tp(mol,mol))
        !        sh_env%C_tp=C
        !end if 
        
         
        deallocate(DGETRI_WORK)
end subroutine adiabaticity


!**********************************************************************************
!> \brief This calculates the adiabatic nacv according to Eq 11 in Carof 17
!paper. This subroutine is stricly bound to AOM. Other Hamiltonians would
!require a different implementation.
!> \author AC [edited SG] 
!**********************************************************************************
SUBROUTINE calc_hop_nacv_aom(sh_env, diagonal_F, hop_nacv, ext_old_state, ext_new_state, symmetry, fast, &
                         exF_11, exF_1, exF_2, exF_3)
        TYPE(sh_env_type)                                       :: sh_env
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE  :: hop_nacv
        INTEGER, OPTIONAL                                        :: ext_old_state, ext_new_state
        CHARACTER(len=*), OPTIONAL                               :: symmetry
        LOGICAL, OPTIONAL                                        :: fast
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE, OPTIONAL  :: exF_11, exF_1, exF_2, exF_3
        REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F


        CHARACTER(len=*), PARAMETER :: routineN = 'calc_hop_nacv_aom', &
        routineP = moduleN//':'//routineN

        INTEGER                                      :: num_mol, nadiab, num_j, natom
        INTEGER                                      :: old_state, new_state
        REAL(KIND=dp)                                :: scaling, cbar, diff_energy
        REAL(KIND=dp)                                :: E_old_state, E_new_state
        !REAL(KIND=dp), dimension(:,:,:), allocatable :: Fin
        !REAL(KIND=dp), dimension(:,:,:), allocatable   :: nacv
        REAL(KIND=dp), dimension(:,:,:,:), allocatable :: Wprime
        integer :: i, j, k, mat, d, l, atoms
        !REAL(KIND=dp), dimension(:,:,:,:), allocatable :: DH, my_DP
        REAL(KIND=dp), dimension(:,:), allocatable :: F_1, F_11, F_2, F_22, F_3
        REAL(KIND=dp), dimension(:), allocatable ::   Sdiag
        REAL(KIND=dp), dimension(:,:), allocatable :: S, C, T, X, LL, Y
        LOGICAL, dimension(:), allocatable         :: ext_pres
        LOGICAL                                    :: my_fast
        REAL(KIND=dp) :: total_norm, error
        REAL(KIND=dp)                           :: convert_nacv 

        ! Based on paper Carof, Giannini, Blumberger, J. Chem. Phys. 147,
        ! 214113 (2017), [Carof17]


        mat    = SIZE(sh_env%Stptp, 1)
        natom  = SIZE(diagonal_F, 2)
        nadiab = SIZE(diagonal_F, 3)
        atoms  = SIZE(sh_env%connlist1, 2)

        my_fast = .FALSE.
        IF (PRESENT(fast)) my_fast = fast

        allocate(T(mat,mat))
        allocate(S(mat,mat))
        allocate(C(mat,mat))
        allocate(X(mat,mat))
        allocate(Y(mat,mat))
        allocate(LL(mat,mat))

        allocate(Sdiag(mat))
        allocate(ext_pres(natom))

        !allocate(DH(mat,mat,natom,3))
        !allocate(my_DP(mat,mat,natom,3)) 
        !allocate(nacv(3,mat, mat*atoms)) 
        allocate(F_1(3,natom))
        allocate(F_11(3,natom))
        allocate(F_2(3,natom))
        allocate(F_22(3,natom))
        allocate(F_3(3,natom))
        !allocate(Fin(3,natom, mat))
        allocate(hop_nacv(3, natom)) 

        scaling = sh_env%scaling
        cbar = sh_env%cbar
        S = sh_env%Stptp
        C = sh_env%C_tp
        
        !nacv = sh_env%nacv
        !Fin = sh_env%F
        ext_pres = sh_env%ext_pres       
        IF (PRESENT(ext_old_state)) THEN
           old_state = ext_old_state
        ELSE
           old_state = sh_env%active_state
        ENDIF
        IF (PRESENT(ext_new_state)) THEN
           new_state = ext_new_state 
        ELSE
           new_state = sh_env%new_state
        ENDIF
        E_old_state = sh_env%E_tp(old_state)
        E_new_state = sh_env%E_tp(new_state)
        IF (E_old_state.EQ.E_new_state) THEN
           diff_energy = 0.0D0
        ELSE
           diff_energy = 1 / (E_new_state - E_old_state)
        ENDIF  

     IF (.NOT.my_fast) THEN
           ! Use the full Eq. (14) of [Carof17]

        T = S
        CALL matinvsqrt(mat, .FALSE., T)
        X = matmul(T, C)                   ! X is T*U
        !CALL adiabaticity(T, LL, Sdiag)     
        CALL adiabaticity(S, LL, Sdiag)    ! LL is L, the matrix that
        ! diagonalizes S and T. Sdiag is the eigenvalues of S.
        Y = matmul( transpose(LL), C)     ! Y is L^{dagger}C

        !CALCULATE W' with Eq. (21) if [Carof17]
        allocate(Wprime(mat,mat,natom,3)) 
        do d=1,3
           num_j = -1
           do j=1,natom
              IF (ext_pres(j)) THEN
                  num_j   = num_j + 1
                  num_mol = int(num_j/atoms) + 1
              ENDIF
              do k=1, mat
                 do l=1, mat
                    Wprime(k,l,j,d) = 0.0d0
                    do i=1,mat 
                       IF (ext_pres(j)) THEN
                          IF (i.ne.num_mol) THEN
                             !Wprime(k,l,j,d) = Wprime(k,l,j,d) + my_DP(i,num_mol,j,d) &
                             Wprime(k,l,j,d) = Wprime(k,l,j,d) + sh_env%nacv(d, i, num_j) &
                                               *( LL(num_mol,k)*LL(i,l)  + LL(num_mol,l) * LL(i, k))
                          ENDIF
                       ENDIF
                    ENDDO
                    Wprime(k, l, j, d) = - Wprime( k, l, j, d) / ( SQRT(Sdiag(k))*SQRT(Sdiag(l)) + Sdiag(l) )
                 ENDDO
              ENDDO
           ENDDO
        ENDDO          
     ENDIF



       convert_nacv = - scaling / cbar
       IF (.NOT.my_fast) THEN
           num_mol = 0
           do d=1,3
              num_j = -1
              do j=1,natom
                 F_1(d,j) = 0.0D0
                 F_11(d,j) = 0.0D0
                 F_2(d,j) = 0.0D0
                 F_3(d,j) = 0.0D0
                 IF (ext_pres(j)) THEN
                    num_j   = num_j + 1
                    num_mol = int(num_j/atoms) + 1 
                 ENDIF
                 do i=1,mat
                    !F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                    F_1(d,j) = F_1(d,j) + C(i,old_state)*diagonal_F(d,j,i)*C(i,new_state)*diff_energy
                    IF (ext_pres(j)) THEN
                       !F_2(d, j) = F_2(d,j) + X(i, old_state) * nacv(d, i, num_mol, j) * X( num_mol, new_state) 
                       F_2(d, j) = F_2(d,j) + X(i, old_state) * sh_env%nacv(d, i, num_j) * X( num_mol, new_state) 
                       IF (i.ne.num_mol) THEN
                          !F_11(d,j) = F_11(d,j) + diff_energy * DH(i,num_mol,j,d) * (              &
                          F_11(d,j) = F_11(d,j) + diff_energy * convert_nacv*sh_env%nacv(d, i, num_j) * (              &
                                                 C(num_mol, new_state) * C(i, old_state) +   &
                                                 C(num_mol, old_state) * C(i,new_state) )
                       ENDIF
                       DO k = 1, mat
                          F_3(d,j) = F_3(d,j) + Y(i, old_state) * Wprime(i, k, j, d) * Y(k, new_state)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO 

           ! Carefull with the sign! 
           hop_nacv = -(F_1 + F_11) + F_2 + F_3
           IF (PRESENT(symmetry)) THEN
              DO j=1, natom
                 DO d=1, 3
                    print*, "NACV", symmetry, hop_nacv(d,j) 
                 ENDDO
              ENDDO
           ENDIF
           IF (PRESENT(exF_1))  exF_1  = F_1
           IF (PRESENT(exF_11)) exF_11 = F_11
           IF (PRESENT(exF_2))  exF_2  = F_2
           IF (PRESENT(exF_3))  exF_3  = F_3
           IF (PRESENT(exF_1)) THEN
              F_22 = F_2 + F_3
              total_norm = 0.0D0
              error = 0.0D0
              DO j=1, 12
                 !DO d=1,3
                 DO d=1, 3
                       error = error + ( hop_nacv(d,j) - (F_1(d,j) + F_11(d,j)))**2 
                       total_norm = total_norm +  (hop_nacv(d,j))**2
                       print*, "hop_nacv", hop_nacv(d,j)
                       print*, "d_1_term", F_1(d,j) + F_11(d,j)
                       print*, "total_norm = ", total_norm
                       print*, "error = ", error
                 ENDDO
                 IF (( mod(j,3).eq.1 ).AND. (j < 12)) THEN
                       print*, "NACVTOTCARBON = ", j, NORM2(hop_nacv(:,j))
                       print*, "NACVFF1CARBON   = ", j, NORM2(F_1(:,j) + F_11(:,j))
                 ELSE 
                       print*, "NACVTOTHYDROGEN = ", j, NORM2(hop_nacv(:,j))
                       print*, "NACVFF1HYDROGEN   = ", j, NORM2(F_1(:,j) + F_11(:,j))
                 ENDIF 
              ENDDO
              print*, "TOTALNORM = ", sqrt(total_norm)
              print*, "ERROR =  ", sqrt(error/total_norm)
           ENDIF
        ELSE
           ! Only calculation of term a_{I, ij} in Eq. (14) of [Carof17], i.e.
           ! Eq. (15)
           num_mol = 0   ! Molecular index 
           do d=1,3      ! Cartesian coordinates
              num_j = 0  ! Active atom index
              do j=1,natom ! Sum over all the atoms
                 F_1(d,j) = 0.0D0
                 F_11(d,j) = 0.0D0
                 IF (ext_pres(j)) THEN
                    num_j   = num_j + 1                  ! Update the active atom index
                    num_mol = int((num_j-1)/atoms) + 1   ! Update the molecular index
                    ! WARNING: we suppose that all molecules have the same
                    ! number of atoms
                 ENDIF
                 do i=1,mat
                    !F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                    ! Diagonal term, it depends on the forces
                    F_1(d,j) = F_1(d,j) + C(i,old_state)*diagonal_F(d,j,i)*C(i,new_state)*diff_energy
                    IF (ext_pres(j)) THEN
                       IF (i.ne.num_mol) THEN
                         ! Of diagonal term, only for active atom j present in
                         ! molecule num_mol
                         F_11(d,j) = F_11(d,j) + diff_energy * convert_nacv*sh_env%nacv(d, i, num_j) * (&
                                                 C(num_mol, new_state) * C(i, old_state) +   &
                                                 C(num_mol, old_state) * C(i,new_state) )
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO 

           ! minus sign to take the same convention of the paper.
           hop_nacv = -(F_1 + F_11) 
        ENDIF
 
        deallocate(T)
        !deallocate(DH)
        !deallocate(my_DP)
        deallocate(F_1)
        deallocate(F_11)
        deallocate(F_2)
        deallocate(F_22)
        deallocate(F_3)
        deallocate(C)
END SUBROUTINE calc_hop_nacv_aom



FUNCTION NORM2(array) RESULT(norm)
         REAL(kind=dp), dimension(:), INTENT(IN) :: array
        
         INTEGER       :: i
         REAL(KIND=dp) :: norm
 
         norm = 0
         DO i=1, SIZE(array)
            norm = norm + array(i)*array(i)
         ENDDO
         norm = SQRT(norm)
END FUNCTION NORM2


!**********************************************************************************
!> \brief 
!**********************************************************************************
SUBROUTINE calc_multiphase_Xmatrix_elements(sh_env, orb, energies_tp, adiab_section, positions_com_tp)
     TYPE(sh_env_type)                         :: sh_env
     REAL(kind=dp), DIMENSION(:),  POINTER     :: energies_tp
     TYPE(section_vals_type), POINTER          :: adiab_section
     INTEGER                                   :: orb ! state space dimension including degeneracy
     REAL(kind=dp), DIMENSION(:,:),POINTER, &
                                   INTENT(IN)        :: positions_com_tp

     CHARACTER(len=*), PARAMETER               :: routineN ='calc_multiphase_Xmatrix_elements', &
     routineP = moduleN//':'//routineN

     !INTEGER, DIMENSION(:,:), ALLOCATABLE      :: active_system_info
     INTEGER                                   :: state, total_states, state1, state2
     !INTEGER :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
     !INTEGER :: init_atom_mol_j, final_atom_mol_j, atoms_per_mol_j
     INTEGER                                   ::ind1, ind2
     INTEGER                                   :: nadiab, i, j!, temp
     !INTEGER                                   :: start_i,end_i,start_j,end_j
     INTEGER, DIMENSION(:), POINTER            :: donor_index, acceptor_index
     INTEGER, DIMENSION(:), ALLOCATABLE        :: first_ind, second_ind
     !INTEGER, DIMENSION(:), ALLOCATABLE        :: my_mol_index !this is the decomp index
     TYPE(section_vals_type), POINTER          :: site_energy_section, ct_block
     INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CT_states_info
     !real(kind=dp), allocatable, dimension(:,:) :: Atp, Btp
     !REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: coeffap, coeffbp
     !REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: coeffas, coeffbs
     !REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE  :: nacva, nacvb
     REAL(KIND=dp)                            :: mol1, mol2
     !INTEGER, DIMENSION(:,:), ALLOCATABLE      ::  bneighbour, aneighbour
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE  :: Sab_list_dphase
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE  :: Sab_list_aphase
     integer                                   :: atoms_per_donor, atoms_per_acceptor, temp, count_


     NULLIFY(site_energy_section,ct_block)

     ! TOTAL number states
     CALL section_vals_val_get(adiab_section,"NUMBER_DIABATIC_STATES",&
                          i_val=nadiab)
     
     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
     ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
     
     CALL section_vals_val_get(ct_block,"INDEX_DONOR_MOLS",&
                          i_vals=donor_index)
     CALL section_vals_val_get(ct_block,"INDEX_ACCEPTOR_MOLS",&
                          i_vals=acceptor_index)

     ! not stricly necessary as this info is already elsewhere
     CALL section_vals_val_get(ct_block,"N_DONOR_ATOMS",&
              i_val=atoms_per_donor)
     CALL section_vals_val_get(ct_block,"N_ACCEPTOR_ATOMS",&
              i_val=atoms_per_acceptor)
  
     ! some checks for the indeces are done in multi_block_energies_forces
     ! subroutine and could be repeated here   
     total_states = SIZE(donor_index)*SIZE(acceptor_index)
     IF (total_states .ne. nadiab) THEN
         CPABORT("NUMBER_DIABATIC_STATES DIFFERENT FROM N*N STATES")
     END IF
     IF (total_states .ne. nadiab*orb) THEN
         CPABORT("DEGENERACY NOT IMPLEMENTED")
     END IF
     IF (total_states .ne. SIZE(energies_tp)) THEN
         CPABORT("PROBLEM WITH SITE ENERGIES NUMBER OF STATES")
     END IF

     !active_system_info = sh_env%active_system_info
     !get the CT index in the 1,2,3... base (active region basis)
     !independently on the actual molecule
     !numbering
     !my_mol_index  = sh_env%index_active_mols
     
     ALLOCATE(second_ind(SIZE(acceptor_index)))
     do i=1,SIZE(acceptor_index)
       do j=1,SIZE(sh_env%index_active_mols)
           if (acceptor_index(i) .eq. sh_env%index_active_mols(j)) then
              second_ind(i) = j
              EXIT
           end if
       end do
     end do
   
     ALLOCATE(first_ind(SIZE(donor_index)))
     do i=1,SIZE(donor_index)
       do j=1,SIZE(sh_env%index_active_mols)
           if (donor_index(i) .eq. sh_env%index_active_mols(j)) then
              first_ind(i) = j
              EXIT
           end if
       end do
     end do     
    
print *, "DECOMP", sh_env%index_active_mols
print *, "DONOR IND", donor_index
print *, "ACCEPTOR IND", acceptor_index
print *, "FIRST IND", first_ind
print *, "SEC IND", second_ind

 
print *, "CALC DONOR PHASE STUFF"
  CALL calc_aom_coup_nacv_single_phase(sh_env, atoms_per_donor, sh_env%ext_elements_d, sh_env%pcoeffs_d, sh_env%scoeffs_d, &
                                          first_ind, sh_env%connlist1_d, sh_env%connlist2_d, sh_env%atomlist_d, &
                                          positions_com_tp, sh_env%active_system_info, sh_env%nacv_dphase, Sab_list_dphase)
print *, "CALC ACCEPTOR PHASE STUFF"
  CALL calc_aom_coup_nacv_single_phase(sh_env, atoms_per_acceptor, sh_env%ext_elements_a, sh_env%pcoeffs_a, sh_env%scoeffs_a, &
                                          second_ind, sh_env%connlist1_a, sh_env%connlist2_a, sh_env%atomlist_a, &
                                          positions_com_tp, sh_env%active_system_info, sh_env%nacv_aphase, Sab_list_aphase)

!TEMP XXX
do temp=1, SIZE(Sab_list_dphase,2)
print *, "Sab_list", Sab_list_dphase(:, temp)
enddo

do temp=1, SIZE(Sab_list_aphase,2)
print *, "Sab_list", Sab_list_aphase(:, temp)
enddo

!PRINT OUT TEST
OPEN (UNIT = 27, FILE = "nacv_formatted-donor.dat",  position="append",action="write")
do j=1, SIZE(first_ind)
   write(27,*) "state", j
   do i=1,atoms_per_donor*SIZE(first_ind) !! HARDCODEC NUMBER OF ACTIVE ATOMS XXX
     write(27,*)  sh_env%nacv_dphase(1,j,i), sh_env%nacv_dphase(2,j,i), sh_env%nacv_dphase(3,j,i)
  end do
end do
close (27)

!PRINT OUT TEST
OPEN (UNIT = 28, FILE = "nacv_formatted-acceptor.dat",  position="append",action="write")
do j=1, SIZE(second_ind)
   write(28,*) "state", j
   do i=1,atoms_per_acceptor*SIZE(second_ind) !! HARDCODEC NUMBER OF ACTIVE ATOMS XXX
     write(28,*)  sh_env%nacv_aphase(1,j,i), sh_env%nacv_aphase(2,j,i), sh_env%nacv_aphase(3,j,i)
  end do
end do
close (28)

   !> loop over states: note that it does not matter if donor molecules comes
   !after acceptor in the xyz or viceversa, or the molecules are DDDAAA or ADADA
   !as the index of each molecule is referred to DECOMP file
   ALLOCATE(CT_states_info(2,total_states))
   state = 0
   do i=1,SIZE(first_ind)
      ind1 = first_ind(i)
      do j=1,SIZE(second_ind)
         ind2 = second_ind(j)
         state = state + 1
         CT_states_info(1,state) = ind1
         CT_states_info(2,state) = ind2
      end do
   end do

   sh_env%H_tp(:,:) = cmplx(0.0D0, 0.0D0,dp)

  !loop over states and fill-in the CT Hamiltonian with book-kept couplings
  do state1=1, SIZE(CT_states_info(1,:))
    !loop over upper diagonal
       do state2=1, SIZE(CT_states_info(2,:)) 
           !if connectivity allows XXXX
           if (state2 >  state1) then
               ! look at the acceptor couplings (the donor are the same)
               if (CT_states_info(1,state1) .eq. CT_states_info(1,state2)) then 
                   ! molecules interacting of the acceptor
                   mol1 = CT_states_info(2,state1)
                   mol2 = CT_states_info(2,state2)
!print *, "states", state1, state2
!print *, "STATE1",  CT_states_info(1,state1), CT_states_info(2,state1)
!print *, "STATE2",  CT_states_info(1,state2), CT_states_info(2,state2)
!print *, "intercating molecules", mol1, mol2
                   do count_=1, SIZE(Sab_list_aphase,2) 
                       if ((mol1 .eq. INT(Sab_list_aphase(1, count_))) .AND. (mol2 .eq. INT(Sab_list_aphase(2, count_)))) then 
!print *, "Sab_list_aphase(3, count_)", Sab_list_aphase(3, count_)
                          !Minus to be consistent with diago
                          sh_env%H_tp(state1, state2) = cmplx(0.0D0, -sh_env%scaling_acceptor*Sab_list_aphase(3, count_),dp) 
                          EXIT
                       end if
                   end do 
               ! look at the donor couplings (the acceptors are the same)
               else if (CT_states_info(2,state1) .eq. CT_states_info(2,state2)) then 
                   mol1 = CT_states_info(1,state1)
                   mol2 = CT_states_info(1,state2)
!print *, "states", state1, state2
!print *, "STATE1",  CT_states_info(1,state1), CT_states_info(2,state1)
!print *, "STATE2",  CT_states_info(1,state2), CT_states_info(2,state2)
!print *, "intercating molecules", mol1, mol2
                   do count_=1, SIZE(Sab_list_dphase,2) 
                       if ((mol1 .eq. INT(Sab_list_dphase(1, count_))) .AND. (mol2 .eq. INT(Sab_list_dphase(2, count_)))) then 
!print *, "Sab_list_aphase(3, count_)", Sab_list_dphase(3, count_)
                          !Minus to be consistent with diago
                          sh_env%H_tp(state1, state2) = cmplx(0.0D0, -sh_env%scaling_donor*Sab_list_dphase(3, count_),dp) 
                          EXIT
                       end if
                   end do 
                
               end if !state matching
  
           else !state
                if (state2 .eq. state1) then
                   !minus because we propagate -X
                   sh_env%H_tp(state1, state2) = cmplx(0.0D0,(-energies_tp(state1)),dp) 
                else 
                   !simmetrization 
                   sh_env%H_tp(state1, state2) = sh_env%H_tp(state2, state1)
                end if 
           end if
        end do
    end do
         
! TEST PRINT H 
do state1=1, SIZE(sh_env%H_tp,1)
   write (*,*) "H", -AIMAG(sh_env%H_tp(state1,:))
end do
!200 format(f8.6)

!print *, "START OTHER MORE CONSUMING H CALCULATION"
!
!   ! SG: It would be probably faster just to calculate couplings within the
!   ! phases, store them (book-keeping) and the fill the H afterwards in a second
!   ! loop. Even the normalization of the AOM coeff should be done once and then
!   ! stored like in Jacob Hamiltonian.
!
!   sh_env%H_tp(:,:) = cmplx(0.0D0, 0.0D0,dp)
!
!!do state=1, SIZE(CT_states_info(1,:))
!!print *,  CT_states_info(1,state), CT_states_info(2,state)
!!enddo
!            
!  do state1=1, SIZE(CT_states_info(1,:))
!    !loop over upper diagonal
!       do state2=1, SIZE(CT_states_info(2,:)) 
!!print *, "states", state1, state2
!           !if connectivity allows XXXX
!           if (state2 >  state1) then 
!!print *, "STATE1",  CT_states_info(1,state1), CT_states_info(2,state1)
!!print *, "STATE2",  CT_states_info(1,state2), CT_states_info(2,state2)
!               if (CT_states_info(1,state1) .eq. CT_states_info(1,state2)) then 
!
!                   !> find index mol i (which is part of acceptor phase) of the pair
!                   init_atom_mol_i = active_system_info(1,CT_states_info(2,state1))
!                   final_atom_mol_i = active_system_info(2, CT_states_info(2,state1))
!                   atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
!           
!                   !> find index mol j (which is part of acceptor phase) of the pair
!                   init_atom_mol_j = active_system_info(1,CT_states_info(2,state2))
!                   final_atom_mol_j = active_system_info(2,CT_states_info(2,state2))
!                   atoms_per_mol_j = (final_atom_mol_j-init_atom_mol_j)+1
!!print *, "index don", init_atom_mol_i, final_atom_mol_i
!!print*, "index acc", init_atom_mol_j, final_atom_mol_j
!
!                    !! THIS CAN BE PUT IN ANOTHER SUBROUTINE      !!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   !> this allocation deallocation is maybe inefficient. Since
!                   !for a given phase molecules have the same number of atoms
!                   !the allocation could be done at the beginning
!                   ALLOCATE(Atp(4,atoms_per_mol_i))
!                   ALLOCATE(Btp(4,atoms_per_mol_j))                    
!                   ALLOCATE(coeffap(atoms_per_mol_i))
!                   ALLOCATE(coeffbp(atoms_per_mol_j))
!             
!                   Atp(1,1:atoms_per_mol_i) = dble(sh_env%ext_elements_a)
!                   Btp(1,1:atoms_per_mol_j) = dble(sh_env%ext_elements_a) 
!                   Atp(2:4,1:atoms_per_mol_i) = positions_com_tp(1:3, init_atom_mol_i:final_atom_mol_i)
!                   Btp(2:4,1:atoms_per_mol_j) = positions_com_tp(1:3, init_atom_mol_j:final_atom_mol_j)
!
!                   !Give donor phase coeffiecients to be normalized (note all
!                   !the coefficients could be normalizized only once in a single
!                   !loop)
!                   coeffap = sh_env%pcoeffs_a ! they are in-out of calc_sab
!                   coeffbp = sh_env%pcoeffs_a ! they are in-out of calc_sab
!                   !sh_env%connlist2_d is the same because i and j are the same
!                   !(donor phase)
!                   !connlistsa = sh_env%connlist2_a
!                   !connlistsd =  sh_env%connlist2_a
!                   !Normalisation of the coefficients without the s orbitals
!                   call calc_sab(Atp, sh_env%connlist2_a, coeffap)
!                   call calc_sab(Btp, sh_env%connlist2_a, coeffbp)
!
!                   !Calculate Sab
!                   call calc_sab(Atp, sh_env%connlist2_a, coeffap, & 
!                             S, Btp, sh_env%connlist2_a, coeffbp)
!
!                   sh_env%H_tp(state1, state2) = cmplx(0.0D0, -sh_env%scaling*S,dp) !Minus to be consistent with diago
!!print *, "IT'S ELECTR COUPL", S
!
!                   ALLOCATE(aneighbour(0:size(sh_env%atomlist_a),0:atoms_per_mol_i)) !atom list is the first dim
!                   ALLOCATE(bneighbour(0:size(sh_env%atomlist_a),0:atoms_per_mol_j)) !atom list is the first dim
!                   ALLOCATE(coeffas(atoms_per_mol_i))
!                   ALLOCATE(coeffbs(atoms_per_mol_j))
!                   ALLOCATE(nacva(3, atoms_per_mol_i))
!                   ALLOCATE(nacvb(3, atoms_per_mol_j))
!                   aneighbour = 0
!                   bneighbour = 0
!                   nacva = 0.0d0
!                   nacvb = 0.0d0
!
!
!                   call sab_screen_new(Atp, sh_env%atomlist_a, sh_env%connlist1_a, coeffap,&
!                         sh_env%atomlist_a, Btp, sh_env%atomlist_a, sh_env%connlist1_a, coeffbp,&
!                         aneighbour(:,1:atoms_per_mol_i), &
!                         bneighbour(:,1:atoms_per_mol_j), &
!                         sh_env%minao)
!
!!print *, "aneighbour"
!!do i=1,size(aneighbour,2)
!!    print*, aneighbour(:,i)
!!end do
!!print *, "bneighbour"
!!do i=1,size(bneighbour,2)
!!    print*, bneighbour(:,i)
!!end do
!                  
!                   coeffas = sh_env%scoeffs_a ! they are in-out of calc_dRSab
!                   coeffbs = sh_env%scoeffs_a ! they are in-out of calc_dRSab
!
!                   call calc_dRSab(Atp, sh_env%atomlist_a, sh_env%connlist1_a, coeffap, coeffas, &
!                         nacva, nacvb, Btp, sh_env%atomlist_a, sh_env%connlist1_a, &
!                         coeffbp, coeffbs, aneighbour, bneighbour, sh_env%dR, sh_env%do_random, sh_env)
!
!!!open(20, file = 'nacv_multip_ap.dat', status = 'old', position="append",action="write")
!!write(20,*) "states", state1, state2
!!write(20,*) "index acc", init_atom_mol_j, final_atom_mol_j
!!write(20,*) 'NAC vector components on the nuclei of the acceptor'
!!do j=1,size(nacvb,2)
!!     !write(*,*) nacva(1,j),  nacva(2,j),  nacva(3,j)
!!     write(20,*) nacvb(1,j),  nacvb(2,j),  nacvb(3,j)
!!end do
!!
!!write(20,*) "index don", init_atom_mol_i, final_atom_mol_i
!!write(20,*) 'NAC vector components on the nuclei of the donor'
!!do j=1,size(nacva,2)
!!     write(20,*) nacva(1,j),  nacva(2,j),  nacva(3,j)
!!     !write(*,*) nacvb(1,j),  nacvb(2,j),  nacvb(3,j)
!!end do
!!!close(20)
!
!
!                   DEALLOCATE(aneighbour, bneighbour)
!                   DEALLOCATE(coeffas, coeffbs)
!                   DEALLOCATE(nacva, nacvb)
!                   DEALLOCATE(Atp,Btp)        
!                   DEALLOCATE(coeffap, coeffbp)        
!                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!               else if (CT_states_info(2,state1) .eq. CT_states_info(2,state2)) then 
!
!                   !> find index mol i (which is state i) of the pair
!                   init_atom_mol_i = active_system_info(1,CT_states_info(1,state1))
!                   final_atom_mol_i = active_system_info(2, CT_states_info(1,state1))
!                   atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
!           
!                   !> find index mol j (which is state j) of the pair
!                   init_atom_mol_j = active_system_info(1,CT_states_info(1,state2))
!                   final_atom_mol_j = active_system_info(2,CT_states_info(1,state2))
!                   atoms_per_mol_j = (final_atom_mol_j-init_atom_mol_j)+1
!!print *, "index don", init_atom_mol_i, final_atom_mol_i
!!print*, "index acc", init_atom_mol_j, final_atom_mol_j
!
!                    !! THIS CAN BE PUT IN ANOTHER SUBROUTINE      !!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   !> this allocation deallocation is maybe inefficient. Since
!                   !for a given phase molecules have the same number of atoms
!                   !the allocation could be done at the beginning
!                   ALLOCATE(Atp(4,atoms_per_mol_i))
!                   ALLOCATE(Btp(4,atoms_per_mol_j))                    
!                   ALLOCATE(coeffap(atoms_per_mol_i))
!                   ALLOCATE(coeffbp(atoms_per_mol_j))
!             
!                   Atp(1,1:atoms_per_mol_i) = dble(sh_env%ext_elements_d) !donor
!                   Btp(1,1:atoms_per_mol_j) = dble(sh_env%ext_elements_d) !donor
!                   Atp(2:4,1:atoms_per_mol_i) = positions_com_tp(1:3, init_atom_mol_i:final_atom_mol_i)
!                   Btp(2:4,1:atoms_per_mol_j) = positions_com_tp(1:3, init_atom_mol_j:final_atom_mol_j)
!
!                   !Give donor phase coeffiecients to be normalized (note all
!                   !the coefficients could be normalizized only once in a single
!                   !loop)
!                   coeffap = sh_env%pcoeffs_d ! they are in-out of calc_sab
!                   coeffbp = sh_env%pcoeffs_d ! they are in-out of calc_sab
!                   !sh_env%connlist2_d is the same because i and j are the same
!                   !(donor phase)
!                   !connlistsa = sh_env%connlist2_a
!                   !connlistsd =  sh_env%connlist2_a
!
!!do temp=1,SIZE(coeffap)
!!print *, "coeffap, coeffbp", coeffap(temp), coeffbp(temp)
!!enddo
!!do temp=1, SIZE(Atp,2)
!!print *, "xyza", Atp(1,temp), Atp(2,temp)
!!enddo
!!do temp=1, SIZE(Btp,2)
!!print *, "xyzb", Btp(1,temp), Btp(2,temp)
!!enddo
!
!                   !Normalisation of the coefficients without the s orbitals
!                   call calc_sab(Atp, sh_env%connlist2_d, coeffap) !donor
!                   call calc_sab(Btp, sh_env%connlist2_d, coeffbp) !donor
!
!!print *, "NORMALIZED"
!!do temp=1,SIZE(coeffap)
!!print *, "coeffap, coeffbp", coeffap(temp), coeffbp(temp)
!!enddo
!
!                   !Calculate Sab
!                   call calc_sab(Atp, sh_env%connlist2_d, coeffap, & 
!                             S, Btp, sh_env%connlist2_d, coeffbp)
!
!                   sh_env%H_tp(state1, state2) = cmplx(0.0D0, -sh_env%scaling*S,dp) !Minus consistent with diagonals
!!print *, "IT'S HOLE COUPL", S
!  
!               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   ALLOCATE(aneighbour(0:size(sh_env%atomlist_d),0:atoms_per_mol_i)) !atom list is the first dim
!                   ALLOCATE(bneighbour(0:size(sh_env%atomlist_d),0:atoms_per_mol_j)) !atom list is the first dim
!                   ALLOCATE(coeffas(atoms_per_mol_i))
!                   ALLOCATE(coeffbs(atoms_per_mol_j))
!                   ALLOCATE(nacva(3, atoms_per_mol_i))
!                   ALLOCATE(nacvb(3, atoms_per_mol_j))
!                   aneighbour = 0
!                   bneighbour = 0
!                   nacva = 0.0d0
!                   nacvb = 0.0d0
!
!
!                   call sab_screen_new(Atp, sh_env%atomlist_d, sh_env%connlist1_d, coeffap,&
!                         sh_env%atomlist_d, Btp, sh_env%atomlist_d, sh_env%connlist1_d, coeffbp,&
!                         aneighbour(:,1:atoms_per_mol_i), &
!                         bneighbour(:,1:atoms_per_mol_j), &
!                         sh_env%minao)
!
!!print *, "aneighbour"
!!do i=1,size(aneighbour,2)
!!    print*, aneighbour(:,i)
!!end do
!!print *, "bneighbour"
!!do i=1,size(bneighbour,2)
!!    print*, bneighbour(:,i)
!!end do
!                  
!                   coeffas = sh_env%scoeffs_d ! they are in-out of calc_dRSab
!                   coeffbs = sh_env%scoeffs_d ! they are in-out of calc_dRSab
!
!                   call calc_dRSab(Atp, sh_env%atomlist_d, sh_env%connlist1_d, coeffap, coeffas, &
!                         nacva, nacvb, Btp, sh_env%atomlist_d, sh_env%connlist1_d, &
!                         coeffbp, coeffbs, aneighbour, bneighbour, sh_env%dR, sh_env%do_random, sh_env)
!
!!!open(20, file = 'nacv_multip.dat', status = 'old', position="append",action="write")
!!write(*,*) "states", state1, state2
!!write(*,*) "index acc", init_atom_mol_j, final_atom_mol_j
!!write(*,*) 'NAC vector components on the nuclei of the acceptor'
!!do j=1,size(nacvb,2)
!!     !write(*,*) nacva(1,j),  nacva(2,j),  nacva(3,j)
!!     write(*,*) nacvb(1,j),  nacvb(2,j),  nacvb(3,j)
!!end do
!!
!!write(*,*) "index don", init_atom_mol_i, final_atom_mol_i
!!write(*,*) 'NAC vector components on the nuclei of the donor'
!!do j=1,size(nacva,2)
!!     write(*,*) nacva(1,j),  nacva(2,j),  nacva(3,j)
!!     !write(*,*) nacvb(1,j),  nacvb(2,j),  nacvb(3,j)
!!end do
!!!close(20)
!
!
!                   DEALLOCATE(aneighbour, bneighbour)
!                   DEALLOCATE(coeffas, coeffbs)
!                   DEALLOCATE(nacva, nacvb)
!
!                   DEALLOCATE(Atp,Btp)        
!                   DEALLOCATE(coeffap, coeffbp)        
!
!               else
!!                   print *, "NOT INTERACTING"
!               end if
!!print *, "index don", init_atom_mol_i, final_atom_mol_i
!!print*, "index acc", init_atom_mol_j, final_atom_mol_j
!           else 
!                if (state2 .eq. state1) then
!                   !minus because we propagate -X
!                   sh_env%H_tp(state1, state2) = cmplx(0.0D0,(-energies_tp(state1)),dp) 
!                else 
!!print *, "LOWER PART"
!                   sh_env%H_tp(state1, state2) = sh_env%H_tp(state2, state1)
!                   
!                end if 
!           end if
!     end do
!  end do
!
!do state1=1, SIZE(sh_env%H_tp,1)
!   write (*,*) "H", -AIMAG(sh_env%H_tp(state1,:))
!end do
!!300 format(f8.6)

   DEALLOCATE(CT_states_info) 
   DEALLOCATE(first_ind)
   DEALLOCATE(second_ind)
   DEALLOCATE(Sab_list_dphase)
   DEALLOCATE(Sab_list_aphase)

END SUBROUTINE


!**********************************************************************************
!> \brief 
!**********************************************************************************
SUBROUTINE calc_aom_coup_nacv_single_phase(sh_env, atoms, elements, phase_pcoeffs, phase_scoeffs, phase_ind, & 
                                           phase_connlist1, phase_connlist2, phase_atomlist, &
                                           positions_com_tp, active_system_info, nacv, Sab_list)

     TYPE(sh_env_type)                                   :: sh_env
     integer, intent(IN)                                 :: atoms
     integer, dimension(:), ALLOCATABLE, INTENT(IN)      :: elements
     !integer, intent(in)                     :: mol!, orb
     INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN)      :: phase_ind !is active region-based
     INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)    :: active_system_info
     REAL(kind=dp), DIMENSION(:,:),POINTER, &
                                   INTENT(IN)            :: positions_com_tp
     REAL(KIND=dp), DIMENSION(:), POINTER, &
                                   INTENT(IN)            :: phase_pcoeffs, phase_scoeffs
     integer, dimension(:,:), POINTER, &
                                   INTENT(IN)            :: phase_connlist1, phase_connlist2
     INTEGER, DIMENSION(:), POINTER, INTENT(IN)          :: phase_atomlist
     !INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)   :: Interac_list 
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)    :: Sab_list
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT)  :: nacv    
  
   
     CHARACTER(len=*), PARAMETER             :: routineN ='calc_aom_coup_nacv_single_phase', &
     routineP = moduleN//':'//routineN

     integer :: i, ind_i, init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i 
     integer :: k, ind_k, init_atom_mol_k, final_atom_mol_k, atoms_per_mol_k!, temp
     real(kind=dp), allocatable, dimension(:,:)          :: Atp, Btp
     !real(kind=dp), allocatable                          :: PC_rvecs(:,:)
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE          :: pcoeffsrun 
     INTEGER, DIMENSION(:,:), ALLOCATABLE                :: bneighbour, aneighbour
     REAL(KIND=dp)                                       :: S 
     INTEGER                                             :: count_
     !REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE        :: nacv
      

     !ALLOCATE(PC_allrvecs(5,mol*atoms))
     !ALLOCATE(PC_rvecs(3,atoms))

     ALLOCATE(Atp(4,atoms))
     ALLOCATE(Btp(4,atoms))
     !similar array to AOM single phase H but without the third dimension (here
     !is just tp). Also pcoeffsrun exists because AOM coeffiient need to be
     !renormalized for each structure
     ALLOCATE(pcoeffsrun(atoms,SIZE(phase_ind))) 

     ! SG: both neighbour and nacv could be allocated alredy in
     ! sh_initialization

     !XXX IF MULTISTATE THEN 
     !atom list is the first and second molecule of the dimer 
     ALLOCATE(aneighbour(0:size(phase_atomlist),0:atoms))
     ALLOCATE(bneighbour(0:size(phase_atomlist),0:atoms))
     aneighbour = 0
     bneighbour = 0
     ! END IF XXXX

print *, "ENTER NEW SINGLE PHASE SUB"

     !variable storage of Sab (two molecules may get e.g. further apart and
     !would not be stored if NN !=1) 
     count_ = 0
     do i=1,SIZE(phase_ind)     
        do k=1, SIZE(phase_ind)
              ! XXX TO ADD IF CONNECTIVITY THEN 
!                   if (sh_env%NN(i,k) == 1) then
           ! count only upper part
           if (k>i) then
              count_ = count_ + 1
           end if
        end do
     end do 

     ! Must be deallocated in the mother subroutine
     ALLOCATE(Sab_list(3,count_))
     !ALLOCATE(Interac_list(2,count_))

     count_ = 0
     ! loop over the index of a single phase
     do i=1,SIZE(phase_ind)
               ind_i = phase_ind(i)
            
               !> find index mol i (which is the donor mol a given phase)
               init_atom_mol_i = active_system_info(1, ind_i)
               final_atom_mol_i = active_system_info(2, ind_i)
               atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

               IF (atoms_per_mol_i .ne. atoms) THEN
                   CPABORT("There is a problem with the num. of atoms per mol")
               END IF
        
               !create the array of all its atomic coordinates at tp
               !numerical factor converts from Angstroms into Bohr radii
               Atp(1,1:atoms_per_mol_i) = dble(elements) ! this is the donor of the pair
               Atp(2:4,1:atoms_per_mol_i) = positions_com_tp(1:3,init_atom_mol_i:final_atom_mol_i)

               !then for each other molecule
               do k=1, SIZE(phase_ind)
                    
                    ind_k = phase_ind(k)
                    !> find index mol j (which is the acceptor mol a given phase)
                    init_atom_mol_k = active_system_info(1, ind_k)
                    final_atom_mol_k = active_system_info(2, ind_k)
                    atoms_per_mol_k = (final_atom_mol_k-init_atom_mol_k)+1

!print *, "index don", init_atom_mol_i, final_atom_mol_i
!print*, "index acc", init_atom_mol_k, final_atom_mol_k

                    IF (atoms_per_mol_i .ne. atoms_per_mol_k) THEN
                       CPABORT("This works only for molecules within a given phase")
                    END IF

                    !create the array of its atomic coordinates at tp (the
                    !number of atoms is still the same as before)
                    Btp(1,1:atoms_per_mol_k) = dble(elements)
                    Btp(2:4,1:atoms_per_mol_k) = positions_com_tp(1:3, init_atom_mol_k:final_atom_mol_k)
!do temp=1,SIZE(coeffap)
!print *, "coeffap, coeffbp", coeffap(temp), coeffbp(temp)
!enddo
!do temp=1, SIZE(Atp,2)
!print *, "xyza", Atp(1,temp), Atp(2,temp)
!enddo
!do temp=1, SIZE(Btp,2)
!print *, "xyzb", Btp(1,temp), Btp(2,temp)
!enddo

                    ! the AOM coefficient must be normalized according to the
                    ! given geometry. Since all teh molecules are the same the
                    ! geometry is taken for all the molecules from the second
                    ! loop k=1,mol when the first loop is fixed
                    if (i == 1) then
                         pcoeffsrun(1:atoms_per_mol_k,k) = phase_pcoeffs(1:atoms_per_mol_k)
                         call calc_sab(Btp, phase_connlist2, pcoeffsrun(:,k))
                    end if

                    !XXX MULTISTATE NEEDS ADJUSTMENT
                    !if (mod(itimes-1, sh_env%multiple_time_step) .eq. 0 & 
                    !      .and. sh_env%multiple_time_step .ne. 1) then 
                           call sab_screen_new(Atp, phase_atomlist, &
                            phase_connlist1, pcoeffsrun(:,i),&
                           & phase_atomlist, Btp, phase_atomlist, &
                            phase_connlist1, pcoeffsrun(:,k),&
                           & aneighbour(:,1:atoms_per_mol_k), bneighbour(:,1:atoms_per_mol_k), & 
                                                                                sh_env%minao)
                     !end if

                   ! TO DO: optimize with neighbouring list 
                   !only call overlap if neighbours list says so
!                   if (sh_env%NN(i,k) == 1) then
                        ! we only do the upper diagonal to save time
                        if (k>i) then 
                           count_  = count_+1
                           call calc_sab(Atp, phase_connlist2, &
                               pcoeffsrun(1:atoms_per_mol_i, i), &
                               S, Btp, phase_connlist2, &
                               pcoeffsrun(1:atoms_per_mol_k, k),&
                               print_rvecs='yes', molA_index=i, molB_index=k)                           
!print *, "IT'S COUPL", phase_ind(i), phase_ind(i), S
                               Sab_list(1, count_) = dble(phase_ind(i))
                               Sab_list(2, count_) = dble(phase_ind(k))
                               Sab_list(3, count_) = dble(S)
!print *, "SEED", sh_env%do_random
                          !IF MULTISTATE THEN 
                           call calc_dRSab(Atp, phase_atomlist, &
                           phase_connlist1, pcoeffsrun(1:atoms_per_mol_i,i), phase_scoeffs, &
                           nacv(1:3, k,(i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i), &
                           nacv(1:3, i,(k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k), &
                           Btp, phase_atomlist, phase_connlist1, pcoeffsrun(1:atoms_per_mol_k,k), &
                            !sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), &
                            !         sh_env%neighbourlist(i,k,:,:), dR, &
                            phase_scoeffs, aneighbour(:,:), bneighbour(:,:), sh_env%dR, &
                            sh_env%do_random, sh_env)                         
                         
                          !END IF
                         
!!open(20, file = 'nacv_multip_inter.dat', status = 'old',
!!position="append",action="write")
!write (*,*) "states", i, k
!write (*,*) "index acc", i, "from to", (k-1)*atoms+1, k*atoms
!write (*,*) "NAC vector components on the nuclei of the acceptor"
!do temp=(k-1)*atoms+1, k*atoms
!write (*,*) nacv(1, i,temp), nacv(2, i,temp), nacv(3, i,temp)
!end do
!
!write (*,*) "index donor", k, "from to", (i-1)*atoms+1, i*atoms
!write (*,*) "NAC vector components on the nuclei of the donor"
!do temp=(i-1)*atoms+1, i*atoms
!write (*,*) nacv(1, k,temp), nacv(2, k,temp), nacv(3, k ,temp)
!enddo
!!close(20)

                       end if !over the upper diag

               end do !loop over k
 
              ! call calc_sab_pete(Atp, phase_connlist2, PC_rvecs)
              ! do o=1, size(PC_rvecs,2)
              !    PC_allrvecs(1,((i-1)*atoms)+o) = i
              !    PC_allrvecs(2,((i-1)*atoms)+o) = o
              !    PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
              ! end do

     end do !loop over i

!!PRINT OUT TEST
!OPEN (UNIT = 27, FILE = "nacv_formatted-newsub.dat", STATUS="NEW")
!
!do k=1, SIZE(phase_ind)
!   write(27,*) "state", k
!   do i=1,atoms*SIZE(phase_ind) !! HARDCODEC NUMBER OF ACTIVE ATOMS XXX
!     write(27,*)  nacv(1,k,i), nacv(2,k,i), nacv(3,k,i)
!  end do
!end do
!close (27)


     !DEALLOCATE(PC_rvecs) 
     DEALLOCATE(Atp,Btp)
     DEALLOCATE(pcoeffsrun)
     DEALLOCATE(aneighbour,bneighbour)


END SUBROUTINE calc_aom_coup_nacv_single_phase


!**********************************************************************************
!> \brief 
!**********************************************************************************
SUBROUTINE calc_aom_Xmatrix_elements(sh_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb, minao)

     TYPE(sh_env_type)                      :: sh_env
     INTEGER, POINTER                        :: itimes
     integer, intent(IN)                     :: atoms
     integer, intent(in)                     :: mol, orb
     REAL(KIND=dp), intent(in)                   :: dR  
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in)                   :: minao
     REAL(KIND=dp), intent(in) :: scaling
     real(kind=dp), DIMENSION(:,:), allocatable, intent(inout) :: PC_allrvecs

     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_aom_Xmatrix_elements', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp)                   :: S
     integer :: i, j, k, l, m, d, o
     real(kind=dp), allocatable :: PC_rvecs(:,:)
      
     !ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))
     !COMMENT: I guess Atm and Btm are only need for the first timestep (itimes = 1)
     !one could duplicate this subroutine and create one only for the first timestep
     !for each molecule i
     do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               do j=1,atoms
                    sh_env%At(1,j)  = dble(sh_env%elements(j))
                    sh_env%Atp(1,j) = dble(sh_env%elements(j))
                    sh_env%Atm(1,j) = dble(sh_env%elements(j))
                    do d=2,4
                         !> CHANGE_AC : already coordonates in atomic units
                         sh_env%Atm(d,j) = sh_env%crds(d-1,(i-1)*atoms+j,1)
                         sh_env%At(d,j)  = sh_env%crds(d-1,(i-1)*atoms+j,2)
                         sh_env%Atp(d,j) = sh_env%crds(d-1,(i-1)*atoms+j,3)
                    end do
               end do

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         sh_env%Bt(1,l)  = dble(sh_env%elements(l))
                         sh_env%Btp(1,l) = dble(sh_env%elements(l))
                         sh_env%Btm(1,l) = dble(sh_env%elements(l))
                         do d=2,4
                              !> CHANGE_AC : already coordonates in atomic units
                              sh_env%Btm(d,l) = sh_env%crds(d-1,(k-1)*atoms+l,1)
                              sh_env%Bt(d,l)  = sh_env%crds(d-1,(k-1)*atoms+l,2)
                              sh_env%Btp(d,l) = sh_env%crds(d-1,(k-1)*atoms+l,3)
                         end do
                    end do
                    ! the AOM coefficient must be normalized according to the
                    ! given geometry. Since all teh molecules are the same the
                    ! geometry is taken for all the molecules from the second
                    ! loop k=1,mol
                    if (i == 1) then
                         do d=1,orb
                              sh_env%pcoeffsrun(:,orb*(k-1)+d,1) = sh_env%pcoeffs(:,d)
                              call calc_sab(sh_env%Btm, sh_env%connlist2, &
                                                 sh_env%pcoeffsrun(:,orb*(k-1)+d,1))

                              sh_env%pcoeffsrun(:,orb*(k-1)+d,2) = sh_env%pcoeffs(:,d)
                              call calc_sab(sh_env%Bt, sh_env%connlist2, &
                                                 sh_env%pcoeffsrun(:,orb*(k-1)+d,2))

                              sh_env%pcoeffsrun(:,orb*(k-1)+d,3) = sh_env%pcoeffs(:,d)
                              call calc_sab(sh_env%Btp, sh_env%connlist2, & 
                                               sh_env%pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if
                    !neighbourlist(:,:,:,:) = sh_env%neighbourlist(:,:,:,:)

                    !!recalculate neighbouring lists from scratch
                    !! you could even reduce the times we calculate neighbouring
                    !list even further
                    if (mod(itimes-1, sh_env%multiple_time_step) .eq. 0 & 
                          .and. sh_env%multiple_time_step .ne. 1) then 
                       do d=1,orb
                            do m=1,orb
                                 call sab_screen_new(sh_env%Atp, sh_env%atomlist, &
                                  sh_env%connlist1, sh_env%pcoeffsrun(:,orb*(i-1)+d,3),&
                                 & sh_env%atomlist, sh_env%Btp, sh_env%atomlist, &
                                  sh_env%connlist1, sh_env%pcoeffsrun(:,orb*(k-1)+m,3),&
                                 & sh_env%aneighbour(:,1:atoms), sh_env%bneighbour(:,1:atoms), & 
                                                                                 minao)
                            end do
                       end do
                     end if

                    !only call overlap if neighbours list says so
                    if (sh_env%NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              CALL diagonal_aom_elem(i,k,orb,itimes,sh_env)
                         !now the off-diagonals                         
                         else ! ( i != k)
                              CALL off_diagonal_aom_elem(sh_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)
                         end if !(i!=k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   sh_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 0.0D0
                                   sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
     !> in case no_diab_nace is true nace in the diabatic basis are not used
     !> neither to
     !> propagate_wf nor in the adiabatic nace calculation for the hopping
     !> probability
     !> so that all the Sttp and co matrices are also not used and not calculated
     !> NB probably the same could be done for sh_env%Stt
                                   if (.NOT. sh_env%no_diab_nace) then
                                       sh_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 0.0D0
                                       sh_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       sh_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       sh_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       sh_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   end if
                              end do
                         end do
                    end if !(NN(i,k) == 1

               end do !loop over k
               call calc_sab_pete(sh_env%Atp, sh_env%connlist2, PC_rvecs)
               do o=1, size(PC_rvecs,2)
                  PC_allrvecs(1,((i-1)*atoms)+o) = i
                  PC_allrvecs(2,((i-1)*atoms)+o) = o
                  PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
               end do
     end do !loop over i
    DEALLOCATE(PC_rvecs) 
END SUBROUTINE calc_aom_Xmatrix_elements

SUBROUTINE diagonal_aom_elem(i,k,orb,itimes,sh_env)

     integer, intent(IN) :: i, k,orb
     INTEGER, POINTER                               :: itimes
     TYPE(sh_env_type)              :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='diagonal_aom_elem', &
     routineP = moduleN//':'//routineN

     integer ::  m, d

     do d=1,orb
          do m=1,orb
               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
               sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
     !> in case no_diab_nace is true nace in the diabatic basis are not used
     !> neither to
     !> propagate_wf nor in the adiabatic nace calculation for the hopping
     !> probability
     !> so that all the Sttp and co matrices are also not used and not calculated
     !> NB probably the same could be done for sh_env%Stt
               if (.NOT. sh_env%no_diab_nace) then
                  sh_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                  sh_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
               end if
               sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, &
                                         (-sh_env%siteenergies(3,k)),dp)
          end do
     end do
     if (itimes.eq.1 .or. sh_env%restart_key)  then
          do d=1,orb
               do m=1,orb
               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
                    sh_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 1.0D0
                    if (.NOT. sh_env%no_diab_nace) then
                      sh_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                      sh_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                      sh_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                    end if
               end do
          end do     
     end if ! (.not.inrun)
END SUBROUTINE diagonal_aom_elem

SUBROUTINE off_diagonal_aom_elem(sh_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)

     TYPE(sh_env_type)                           :: sh_env
     integer, intent(IN) :: i, k,orb, atoms
     INTEGER, POINTER                               :: itimes
     REAL(KIND=dp), intent(in)                   :: dR  
     REAL(KIND=dp), intent(out)                  :: S
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in) :: scaling

     CHARACTER(len=*), PARAMETER :: routineN ='off_diagonal_aom_elem', &
     routineP = moduleN//':'//routineN

     integer ::  m, d
     !symmetric matrices
     !no sense calling calc_sab twice as often as you need to
     if (k < i) then
          do m=1,orb
               do d=1,orb
                    sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = sh_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                    sh_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = &
                                             sh_env%Stt(orb*(k-1)+m,orb*(i-1)+d)
                    sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = &
                                            sh_env%Stptp(orb*(k-1)+m,orb*(i-1)+d)
     !>in case no_diab_nace is true nace in the diabatic basis are not used
     !>neither to
     !>propagate_wf nor in the adiabatic nace calculation for the hopping
     !>probability
     !>so that all the Sttp and co matrices are also not used and not calculated
     !>NB probably the same could be done for sh_env%Stt
                    if (.NOT. sh_env%no_diab_nace) then
                       sh_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = &
                                                sh_env%Stmtm(orb*(k-1)+m,orb*(i-1)+d)
                       !A cross symmetry for Sttp, Stpt,
                       !Sttm, Stmt
                       sh_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                sh_env%Stpt(orb*(k-1)+m,orb*(i-1)+d)
                       sh_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = & 
                                                 sh_env%Sttp(orb*(k-1)+m,orb*(i-1)+d)
                       if (itimes.eq.1 .or. sh_env%restart_key)  then
                          sh_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                 sh_env%Stmt(orb*(k-1)+m,orb*(i-1)+d)
                          sh_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                 sh_env%Sttm(orb*(k-1)+m,orb*(i-1)+d)
                       ENDIF
                    end if
               end do
          end do
     else
          do d=1,orb
               do m=1,orb
               !PC added print_rvecs =yes so p vectors and coords would be printed 
                    call calc_sab(sh_env%Atp, sh_env%connlist2, &
                        sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                     S, sh_env%Btp, sh_env%connlist2, &
                               sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3),& 
                     print_rvecs='yes', molA_index=i, molB_index=k)
                    sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                    sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

                    !> In the following diabatic non-orthog. nacv are
                    !calculated. Since calculating them is very time consuming
                    !and expensive we opted for a multiple time step scheme.
                    !nacv within the counter(chosen in the input 
                    !file) are considered a
                    !constant and taken from the prevous time step.
                    !Calcul of the NACV
                    if (mod(itimes-1, sh_env%multiple_time_step) .eq. 0 &
                          .and. sh_env%multiple_time_step .ne. 1) then 
                        call calc_dRSab(sh_env%Atp, sh_env%atomlist, &
                                      sh_env%connlist1, &
                        sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), sh_env%scoeffs, &
                        sh_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                        sh_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                        sh_env%Btp, sh_env%atomlist, sh_env%connlist1, &
                        sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3), &
                         !sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), &
                         !         sh_env%neighbourlist(i,k,:,:), dR, &
                         sh_env%scoeffs, sh_env%aneighbour(:,:), &
                                          sh_env%bneighbour(:,:), dR, &
                         do_random, sh_env)
                        !> note that calc_dRSab uses a radom number to randomize
                        !forwards derivative in the NACV (this improves energy
                        !conservation by one/two order of magnitude. To be
                        !consistent with random sequence, one should:
                       !else
                             !call random generator function

                   else if (sh_env%multiple_time_step .eq. 1) then 
                        call calc_dRSab(sh_env%Atp, sh_env%atomlist, &
                                      sh_env%connlist1, &
                        sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), sh_env%scoeffs, &
                        sh_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                        sh_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                        sh_env%Btp, sh_env%atomlist, sh_env%connlist1, &
                        sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3), &
                         sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), &
                                  sh_env%neighbourlist(i,k,:,:), dR, &
                         do_random, sh_env)
                        !> note that calc_dRSab uses a radom number to randomize
                        !forwards derivative in the NACV (this improves energy
                        !conservation by one/two order of magnitude. To be
                        !consistent with random sequence, one should:
                       !else
                             !call random generator function
          
                    end if
               end do
          end do
     ! in case no_diab_nace is true nace in the diabatic basis are not used
     ! neither to
     ! propagate_wf nor in the adiabatic nace calculation for the hopping
     ! probability
     ! so that all the Sttp and co matrices are also not used and not calculated
     ! NB probably the same could be done for sh_env%Stt and sh_env%Stptp
         if (.NOT. sh_env%no_diab_nace) then
             do d=1,orb
                  do m=1,orb
                       call calc_sab(sh_env%At, sh_env%connlist2, &
                        sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                        sh_env%Btp,&
                       & sh_env%connlist2, sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                       sh_env%Sttp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
                       call calc_sab(sh_env%Atp, sh_env%connlist2, &
                            sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), S,& 
                        sh_env%Bt,&
                       & sh_env%connlist2, sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                       sh_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                  end do
             end do
             if (itimes.eq.1 .or. sh_env%restart_key)  then
                  do d=1,orb
                       do m=1,orb
                            call calc_sab(sh_env%Atm, sh_env%connlist2, &
                            &sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                                sh_env%Btm, sh_env%connlist2,&
                            &sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                            sh_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
    
                            call calc_sab(sh_env%At, sh_env%connlist2, &
                            &sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                                       sh_env%Btm, sh_env%connlist2,&
                            &sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                            sh_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
                            call calc_sab(sh_env%Atm, sh_env%connlist2, &
                            &sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                              sh_env%Bt, sh_env%connlist2, &
                            &sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                            sh_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                       end do
                  end do
    
             end if ! (.not.inrun)
         end if 
     end if ! ( k < i)
END SUBROUTINE off_diagonal_aom_elem

!**********************************************************************************
!> \brief 
!**********************************************************************************
SUBROUTINE calc_aom_Xmatrix_elements_first_step(sh_env, atoms, mol, dR, do_random,cbar, &
                      scaling, PC_allrvecs, orb, minao)

     TYPE(sh_env_type)                          :: sh_env
     integer, intent(IN)                         :: atoms
     integer, intent(in)                         :: mol, orb
     REAL(KIND=dp), intent(in)                   :: dR
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in)                   :: minao
     REAL(KIND=dp), intent(in)                   :: scaling
     real(kind=dp), DIMENSION(:,:), allocatable, intent(inout) :: PC_allrvecs

     CHARACTER(len=*), PARAMETER                 :: routineN='calc_aom_Xmatrix_elements_first_step', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp)                               :: S
     integer                                     :: i, j, k, l, m, d, o, temp
     real(kind=dp), allocatable                  :: PC_rvecs(:,:)


     !ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))

     !CALCULATE THE HAMILTONIAN AND THE NACV
     do i=1,mol
print *, "state i", i
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               sh_env%Atp(1,1:atoms) = dble(sh_env%elements((i-1)*atoms+1:i*atoms))
               !sh_env%Atp(1,1:atoms) = dble(sh_env%elements(:))
               do j=1,atoms
                    sh_env%Atp(2:4,j) = sh_env%crds(1:3,(i-1)*atoms+j,3)
               end do

               !then for each other molecule
               do k=1,mol
print *, "state k", k
                    !create the array of its atomic coordinates at both t and t+dt
                    sh_env%Btp(1,1:atoms) = dble(sh_env%elements((k-1)*atoms+1:k*atoms))
                    !sh_env%Btp(1,1:atoms) = dble(sh_env%elements(:))
                    do l=1,atoms
                        sh_env%Btp(2:4,l) = sh_env%crds(1:3,(k-1)*atoms+l,3)
                    end do


                    if (i == 1) then
!do temp=1,atoms
!print *, "sh_env%pcoeffs", sh_env%pcoeffs(temp,1)
!enddo
                         ! SG: renormalize coefficient (atoms, state, tstep) of all k (meaning all molecules)
                         ! this is a trick to do it once only for all molecules
                         ! and store them in a more general pcoeffsrun that
                         ! contains them all. This could be extended by
                         ! considering that each molecule may have different n
                         ! atoms.
                         do d=1,orb
                              sh_env%pcoeffsrun(1:atoms,orb*(k-1)+d,3) = sh_env%pcoeffs(1:atoms,d)
                              call calc_sab(sh_env%Btp, sh_env%connlist2, sh_env%pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
!do temp=1,atoms
!print  *, "sh_env%pcoeffsrun", sh_env%pcoeffsrun(temp, k, 3)
!enddo
                    end if  ! (i == 1)

!print *, "STATES", i, k
!do temp=1,atoms
!print  *, "sh_env%pcoeffsrun entering", sh_env%pcoeffsrun(temp, i, 3), sh_env%pcoeffsrun(temp, k, 3)
!enddo

                    do d=1,orb
                         do m=1,orb
                              call sab_screen_new(sh_env%Atp, sh_env%atomlist, sh_env%connlist1, &
                                                                      sh_env%pcoeffsrun(:,orb*(i-1)+d,3),&
                              & sh_env%atomlist, sh_env%Btp, sh_env%atomlist, sh_env%connlist1, &
                                                                     sh_env%pcoeffsrun(:,orb*(k-1)+m,3),&
                              & sh_env%aneighbour(:,1:atoms), sh_env%bneighbour(:,1:atoms), minao)
                         end do
                    end do
                    !store neighbouring list when we do not use MTS feature
                    if (sh_env%multiple_time_step .eq. 1) then
                        sh_env%neighbourlist(i,k,:,1:atoms) = sh_env%bneighbour(:,1:atoms)
                        sh_env%neighbourlist(k,i,:,1:atoms) = sh_env%aneighbour(:,1:atoms)
                    end if 
!print *, "sh_env%aneighbour(:,1:atoms)", sh_env%aneighbour(:,1:atoms)
!print *, "sh_env%aneighbour(:,1:atoms)", SIZE(sh_env%aneighbour(1,1:atoms))
!print *, "sh_env%aneighbour(:,1:atoms)", SIZE(sh_env%aneighbour(:,1))
!print *, "sh_env%bneighbour(:,1:atoms)", sh_env%bneighbour(:,1:atoms)
                    if (sh_env%NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                                        sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-sh_env%siteenergies(3,k)),dp)
                                   end do
                              end do

                         !now the off-diagonals                         
                         else ! ( i == k)
                              !symmetric matrices
                              !no sense calling calc_sab twice as often as you need to
                              if (k < i) then
                                   do m=1,orb
                                        do d=1,orb
                                             sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = sh_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                                             sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = sh_env%Stptp(orb*(k-1)+m,orb*(i-1)+d)
                     
                             !call calc_sab(sh_env%Atp, sh_env%connlist2, &
                             !                             sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                             !                             S, sh_env%Btp, sh_env%connlist2, & 
                             !                             sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                             !                sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                             !                sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, -scaling*S,dp)

                                           ! call calc_dRSab(sh_env%Atp, sh_env%atomlist, sh_env%connlist1, &
                                           ! sh_env%pcoeffsrun(1:atoms,orb*(k-1)+d,3), sh_env%scoeffs, &
                                           ! sh_env%nacv_tot(1:3,orb*(i-1)+d,orb*(k-1)+m,(i-1)*atoms+1:i*atoms), &
                                           ! sh_env%nacv_tot(1:3,orb*(k-1)+m,orb*(i-1)+d,(k-1)*atoms+1:k*atoms), &
                                           ! !sh_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                                           ! !sh_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                                           ! sh_env%Btp, sh_env%atomlist, sh_env%connlist1, & 
                                           ! sh_env%pcoeffsrun(1:atoms,orb*(i-1)+m,3),&
                                           ! !sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), sh_env%neighbourlist(i,k,:,:), dR, &
                                           ! sh_env%scoeffs, sh_env%aneighbour(:,:), sh_env%bneighbour(:,:), dR, &
                                           ! do_random, sh_env)
                                           !

                                        end do
                                   end do
                              else
print *, "OFF-DIAGONAL"
                                   do d=1,orb
                                        do m=1,orb
!do temp=1,atoms
!print *, "sh_env%Atp just x", sh_env%Atp(1,temp), sh_env%Atp(2,temp)
!print *, "sh_env%Btp just x", sh_env%Btp(1,temp), sh_env%Btp(2,temp)
!enddo
!print *, "sh_env%connlist2", sh_env%connlist2
!do temp=1,atoms
!print  *, "sh_env%pcoeffsrun(temp, i, 3)", sh_env%pcoeffsrun(temp, i, 3)
!enddo

                                             call calc_sab(sh_env%Atp, sh_env%connlist2, &
                                                          sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                                          S, sh_env%Btp, sh_env%connlist2, & 
                                                          sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, -scaling*S,dp)

print *, "SEED", do_random
                                             call calc_dRSab(sh_env%Atp, sh_env%atomlist, sh_env%connlist1, &
                                             sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), sh_env%scoeffs, &
                                             !nacv(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                             !nacv(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                             sh_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                                             sh_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                                             sh_env%Btp, sh_env%atomlist, sh_env%connlist1, & 
                                             sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3),&
                                             !sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), sh_env%neighbourlist(i,k,:,:), dR, &
                                             sh_env%scoeffs, sh_env%aneighbour(:,:), sh_env%bneighbour(:,:), dR, &
                                             do_random, sh_env)

!open(20, file = 'nacv_multip.dat', status = 'old', position="append",action="write")
write (*,*) "states", i, k
write (*,*) "index acc", orb*(i-1)+m, "from to", (k-1)*atoms+1, k*atoms
write (*,*) "NAC vector components on the nuclei of the acceptor"
do temp=(k-1)*atoms+1, k*atoms
write (*,*) sh_env%nacv(1,orb*(i-1)+m,temp), sh_env%nacv(2,orb*(i-1)+m,temp), sh_env%nacv(3,orb*(i-1)+m,temp)
end do

write (*,*) "index donor", orb*(k-1)+d, "from to", (i-1)*atoms+1, i*atoms
write (*,*) "NAC vector components on the nuclei of the donor"
do temp=(i-1)*atoms+1, i*atoms
write (*,*) sh_env%nacv(1,orb*(k-1)+d,temp), sh_env%nacv(2,orb*(k-1)+d,temp), sh_env%nacv(3,orb*(k-1)+d,temp)
enddo 
!close(20)

                                           !!TEST XXXX NACV TOT 
                                           ! call calc_dRSab(sh_env%Atp, sh_env%atomlist, sh_env%connlist1, &
                                           ! sh_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), sh_env%scoeffs, &
                                           ! sh_env%nacv_tot(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                           ! sh_env%nacv_tot(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                           ! !sh_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                                           ! !sh_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                                           ! sh_env%Btp, sh_env%atomlist, sh_env%connlist1, & 
                                           ! sh_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3),&
                                           ! !sh_env%scoeffs, sh_env%neighbourlist(k,i,:,:), sh_env%neighbourlist(i,k,:,:), dR, &
                                           ! sh_env%scoeffs, sh_env%aneighbour(:,:), sh_env%bneighbour(:,:), dR, &
                                           ! do_random, sh_env)


                                        end do
                                   end do
                              endif 
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   sh_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1
               end do !loop over k 

               call calc_sab_pete(sh_env%Atp, sh_env%connlist2, PC_rvecs)
               do o=1, size(PC_rvecs,2)
                  PC_allrvecs(1,((i-1)*atoms)+o) = i
                  PC_allrvecs(2,((i-1)*atoms)+o) = o
                  PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
               end do

     enddo !loop-over i
     DEALLOCATE(PC_rvecs)

!PRINT OUT TEST 
!OPEN (UNIT = 24, FILE = "nacv_formatted.dat", STATUS="NEW")

!do k=1, mol
!   write(24,*) "state", k
!   do i=1,SIZE(sh_env%nacv,3) !! HARDCODEC NUMBER OF ACTIVE ATOMS XXX
!     write(24,*)  sh_env%nacv(1,k,i), sh_env%nacv(2,k,i), sh_env%nacv(3,k,i)
!  end do  
!end do
!close (24)

!!PRINT OUT TEST 
!OPEN (UNIT = 25, FILE = "nacv_formatted_total.dat", STATUS="NEW")
!
!do k=1, mol
!do i=1, mol
!   write(25,*) "state", k,i
!   do ci=1,54
!     write(25,*)  sh_env%nacv_tot(1,k,i, ci), sh_env%nacv_tot(2,k,i,ci), sh_env%nacv_tot(3,k,i, ci)
!  end do  
!end do
!end do
!close (25)
END SUBROUTINE calc_aom_Xmatrix_elements_first_step


  
!**********************************************************************************
!> \brief 
!> \Note this subroutine could be optimized by putting the local
!allocation/deallocation step only at the first step and then overwriting the
!matrices like it is done in nace1 
!**********************************************************************************
subroutine makeX(mol, H, Stt, Sttp, Sttm, Stptp, Stmtm, X, deltat, NACE)
        !X is the name I give to the matrix in the solution to the Schroedinger equation:
        !dc/dt = X*sh_env%ctin
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
       !dc/dt = X*sh_env%ctin
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
        !dc/dt = X*sh_env%ctin
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
subroutine old_diab_to_adiab( sh_env, C2, old_state, mol, cstate, ct, C1, cttm_ad, cstate_tm)
        REAL(KIND=dp), dimension(:,:), intent(in) :: C2
        integer, intent(in) :: old_state
        TYPE(sh_env_type) :: sh_env
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
                cstate = cstate + cmplx(C2(i,old_state), 0.0D0,dp)*sh_env%ctin(i)
                IF (PRESENT(C1)) THEN
                     cstate_tm = cstate_tm + cmplx(C1(i,old_state),0.0D0,dp)*sh_env%ctin(i)
                END IF
        end do
        do i=1,mol
                !calculate a and b
                ct(i) = cmplx(0.0D0, 0.0D0,dp)
                IF (PRESENT(cttm_ad)) cttm_ad(i) = cmplx(0.0D0, 0.0D0,dp)
                do j=1,mol
                        ct(i) = ct(i) + cmplx(C2(j,i), 0.0D0,dp)*sh_env%ctin(j)
                        IF (PRESENT(C1)) THEN
                              cttm_ad(i) = cttm_ad(i) + cmplx(C1(j,i),0.0D0,dp)*sh_env%ctin(j)
                        END IF
                end do
        end do
end subroutine old_diab_to_adiab


!**********************************************************************************
!> \brief Calculate the diabatic NACE first at time Tp
!>    NB: this subroutine calculates diabatic NACE (d_kl) in the othogonal basis and in
!>    practice it does exactly the same job as makeX subroutine with the exception
!>    that does not calculate the term H_kl - d_kl
!> \param diabatic_nace_t and diab_D are the outputs of these subroutine they
!>    are the same array but the former is further simmetrized in order to be
!>    transferred to the propagation subroutine with linear interpolation of
!>    diabatic_nace_t
!> \note in all the following nace subroutines: "mol" should be replaced by "mat", 
!>    the latter can deal with degeneracy present in the H   
!**********************************************************************************
subroutine nace1_time_tp(sh_env, mol, deltat)

    TYPE(sh_env_type)                         :: sh_env
    INTEGER                                    :: mol!, i, j!, sub, add
    REAL(KIND=dp)                              :: deltat

    CHARACTER(len=*), PARAMETER :: routineN ='nace1', &
                routineP = moduleN//':'//routineN

    !> Refer to Overleaf Note by SG  
    ! note T1, T2, Sqrt are just supporting matrices they get overwrited
     
    !*********  FIRST PART *****************
    !use the matinvsqrt subroutine to find T = S^(-1/2)
    sh_env%adiab_nace_env%T1 = sh_env%Stt
    sh_env%adiab_nace_env%T2 = sh_env%Stptp
    !find S^{-1/2} at given time
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T1)
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T2)

    !find S^{1/2} at a given time
    sh_env%adiab_nace_env%Ssqrt = sh_env%Stptp
    call matinvsqrt(mol, .FALSE., sh_env%adiab_nace_env%Ssqrt)

    !Time deriv of Tdot(at time Tp) = DT = [STpTp^{-1/2} - STT^{-1/2}]*(1/dt) 
    sh_env%adiab_nace_env%DT = (1.0D0/deltat)*(sh_env%adiab_nace_env%T2 - sh_env%adiab_nace_env%T1)

    !Second term (at time Tp) = T*S*DT = S^{1/2}*DT = S_TpTp^{1/2}*DT  
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
                sh_env%adiab_nace_env%Ssqrt, mol, &
                sh_env%adiab_nace_env%DT, mol, &
                0.0_DP, sh_env%adiab_nace_env%diab_D, mol)


    !*********  SECOND  PART *****************
    !Overwrites T1 to save space so now: T1 = STpTp^{-1/2}
    sh_env%adiab_nace_env%T1 = sh_env%Stptp
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T1)

    !Time deriv Dprime (at time Tp) = [STpTp^{-1/2} - STpT^{-1/2}]*(1/dt) 
    sh_env%adiab_nace_env%DS = (1.0D0/deltat)*(sh_env%Stptp - sh_env%Stpt)

    ! Overwriting T2 to save space 
    !T2 = DS*T1 = DS*STpTp^{-1/2}
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
               sh_env%adiab_nace_env%DS, mol, &
               sh_env%adiab_nace_env%T1, mol, &
               0.0_DP, sh_env%adiab_nace_env%T2, mol)

    !The total orthogonal diabatic NACE (at time Tp)= T*D*T + T*S*DT
    !Diab NACE = T1^{dagg}*T2 + T*S*DT
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                sh_env%adiab_nace_env%T1,  mol, &
                sh_env%adiab_nace_env%T2, mol, &       ! Now diab_D becomes = T*D*T + T*S*DT
                1.0_DP, sh_env%adiab_nace_env%diab_D, mol)

!!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "D Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", sh_env%adiab_nace_env%diab_D(i,j), "  "
!   END DO
!   PRINT *, " "
!END DO

!NB: becuase of the way finite differences are defined actually here
!sh_env%diabatic_nace_t refers at time T (namely t-dt) not Tp 
   sh_env%diabatic_nace_tp = (sh_env%adiab_nace_env%diab_D - & 
                              TRANSPOSE(sh_env%adiab_nace_env%diab_D)) / 2
!> this previous line would
!> antisimmetrize NACE even in the calculation of adiab nace, 
!> We have shown that it makes no difference at all!

end subroutine nace1_time_tp


!**********************************************************************************
!> \brief Calculate the diabatic NACE first
!>    NB: this subroutine calculates diabatic NACE (d_kl) in the othogonal basis and in
!>    practice it does exactly the same job as makeX subroutine with the exception
!>    that does not calculate the term H_kl - d_kl
!> \param diabatic_nace_t and diab_D are the outputs of these subroutine they
!>    are the same array but the former is further simmetrized in order to be
!>    transferred to the propagation subroutine with linear interpolation of
!>    diabatic_nace_t
!> \note in all the following nace subroutines: "mol" should be replaced by "mat", 
!>    the latter can deal with degeneracy present in the H   
!> \note the orthogonal NACE calculated here refer to time T!, like H_t (not Tp) 
!**********************************************************************************
subroutine nace1_time_t(sh_env, mol, deltat)

    TYPE(sh_env_type)                         :: sh_env
    INTEGER                                    :: mol!, i, j!, sub, add
    REAL(KIND=dp)                              :: deltat

    CHARACTER(len=*), PARAMETER :: routineN ='nace1', &
                routineP = moduleN//':'//routineN

    !> Based on Eq. 19 of paper: Spencer et al., J. Chem. Phys. 145, 064102 (2016)
    !> Refer to Overleaf Note by SG  
    ! note T1, T2, Sqrt are just supporting matrices they get overwrited
     
    !*********  FIRST PART *****************
    !use the matinvsqrt subroutine to find T = S^(-1/2)
    sh_env%adiab_nace_env%T1 = sh_env%Stmtm
    sh_env%adiab_nace_env%T2 = sh_env%Stptp
    !find S^{-1/2} at given time
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T1)
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T2)

    !find S^{1/2} at a given time
    sh_env%adiab_nace_env%Ssqrt = sh_env%Stt
    call matinvsqrt(mol, .FALSE., sh_env%adiab_nace_env%Ssqrt)

    !Time deriv of Tdot(at time T) = DT = [STpTp^{-1/2} - STmTm^{-1/2}]*(1/2dt) 
    sh_env%adiab_nace_env%DT = (0.5D0/deltat)*(sh_env%adiab_nace_env%T2 - sh_env%adiab_nace_env%T1)

    !Second term (at time T) = T*S*DT = S^{1/2}*DT = S_TT^{1/2}*DT  
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
                sh_env%adiab_nace_env%Ssqrt, mol, &
                sh_env%adiab_nace_env%DT, mol, &
                0.0_DP, sh_env%adiab_nace_env%diab_D, mol)


    !*********  SECOND  PART *****************
    !Overwrites T1 to save space so now: T1 = STT^{-1/2}
    sh_env%adiab_nace_env%T1 = sh_env%Stt
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T1)

    !Time deriv Dprime (at time T) = [STTp^{-1/2} - STTm^{-1/2}]*(1/2dt) 
    sh_env%adiab_nace_env%DS = (0.5D0/deltat)*(sh_env%Sttp - sh_env%Sttm)

    ! Overwriting T2 to save space 
    !T2 = DS*T1 = DS*STT^{-1/2}
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
               sh_env%adiab_nace_env%DS, mol, &
               sh_env%adiab_nace_env%T1, mol, &
               0.0_DP, sh_env%adiab_nace_env%T2, mol)

    !The total orthogonal diabatic NACE (at time T)= T*D*T + T*S*DT
    !Diab NACE = T1^{dagg}*T2 + T*S*DT
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                sh_env%adiab_nace_env%T1,  mol, &
                sh_env%adiab_nace_env%T2, mol, &       ! Now diab_D becomes = T*D*T + T*S*DT
                1.0_DP, sh_env%adiab_nace_env%diab_D, mol)

!!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "D Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", sh_env%adiab_nace_env%diab_D(i,j), "  "
!   END DO
!   PRINT *, " "
!END DO

!NB: becuase of the way finite differences are defined actually here
!sh_env%diabatic_nace_t refers at time T (namely t-dt) not Tp 
   sh_env%diabatic_nace_t = (sh_env%adiab_nace_env%diab_D - & 
                              TRANSPOSE(sh_env%adiab_nace_env%diab_D)) / 2
!> this previous line would
!> antisimmetrize NACE even in the calculation of adiab nace, 
!> We have shown that it makes no difference at all!

end subroutine nace1_time_t


!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 2.
!**********************************************************************************
subroutine nace2(sh_env, mol, print_here)

        TYPE(sh_env_type)                         :: sh_env
        INTEGER                                    :: mol
        LOGICAL                                    :: print_here

        CHARACTER(len=*), PARAMETER :: routineN ='nace2', &
                    routineP = moduleN//':'//routineN

        INTEGER                                    :: i, j, k, new_index
        REAL(KIND=dp)                              :: max_val, threshold
        LOGICAL                                    :: go_on

        !> 2. Check reordering of eigenvectors between the different timestep.
        sh_env%adiab_nace_env%map_new2old = 0
        sh_env%adiab_nace_env%map_old2new = 0
        sh_env%adiab_nace_env%not_associated = 0
        sh_env%adiab_nace_env%zeros = 0
        sh_env%adiab_nace_env%basic = (/( i, i=1, mol )/)
        threshold = 0.1
        !> First, we associate the new state i with the old state k if overlap is above
        !> (1 - threshold)
        do i=1, mol !loop over new index
           max_val = 0.0
           do j= 1, mol !loop over old index
              sh_env%adiab_overlap(i,j) = 0.0
              do k = 1, mol
                 sh_env%adiab_overlap(i,j) = sh_env%adiab_overlap(i,j) + sh_env%C_tp(k,i)*sh_env%C_t(k,j)
              enddo
              if (abs(sh_env%adiab_overlap(i,j)).GE.abs(max_val)) THEN
                 max_val = sh_env%adiab_overlap(i,j)
                 new_index = j
              endif
           enddo
           if ( abs( max_val**2 - 1 ).LT.threshold) then
              sh_env%adiab_nace_env%map_new2old(i) = new_index
              sh_env%adiab_nace_env%map_old2new( sh_env%adiab_nace_env%map_new2old(i) ) = i
           endif
           !print_here = .TRUE.
           if (print_here) then
              if (i.NE.new_index) then
                 if ( abs( max_val**2 - 1 ).LT.threshold) then
                    print*, "FOBSH| Reordering: for state", i, "old index was", new_index, &
                            "with overlap", max_val, "and former", sh_env%adiab_overlap(i,i)
                 else
                    print*, "FOBSH! Reorgering: for state", i, " we stay in with overlap", &
                            sh_env%adiab_overlap(i,i), " instead of goind to ", new_index, "with overlap", max_val
                 endif
              endif
           endif
        enddo
        IF (ANY(sh_env%adiab_nace_env%map_new2old.NE.sh_env%adiab_nace_env%basic) .and. print_here) THEN
            print*, "FOBSH| New map: ", sh_env%adiab_nace_env%map_new2old
        END IF

        !> If we've associated state i with state j, we try to associate j with if if the
        !> overlap is above 0.5
        IF (ANY(sh_env%adiab_nace_env%map_new2old.EQ.sh_env%adiab_nace_env%zeros)) THEN
         do i=1,mol
            IF ((sh_env%adiab_nace_env%map_new2old(i).EQ.0).AND.(sh_env%adiab_nace_env%map_old2new(i).NE.0)) THEN
               IF ( abs( (sh_env%adiab_overlap(i, sh_env%adiab_nace_env%map_old2new(i)))**2 - 1 ).LT.0.5 ) THEN
                  sh_env%adiab_nace_env%map_new2old(i) = sh_env%adiab_nace_env%map_old2new(i)
                  sh_env%adiab_nace_env%map_old2new( sh_env%adiab_nace_env%map_new2old(i) ) = i
               ENDIF
            ENDIF
         enddo
        IF (ANY(sh_env%adiab_nace_env%map_new2old.NE.sh_env%adiab_nace_env%basic) .and. print_here) THEN
            print*, "FOBSH! New map again:", sh_env%adiab_nace_env%map_new2old
        END IF
       ENDIF

        !> Finally, we associate the not associate j with the remaining j state by order
        !> of energy
        IF (ANY(sh_env%adiab_nace_env%map_new2old.EQ.sh_env%adiab_nace_env%zeros)) THEN
         k = 1
         do i=1, mol
               j = 1
               go_on = .TRUE.
               do while ((go_on).AND.(j.LE.mol))
                  if (i.eq.sh_env%adiab_nace_env%map_new2old(j)) then
                     go_on = .FALSE.
                  else
                     j = j + 1
                  endif
               enddo
               if (go_on) then
                  sh_env%adiab_nace_env%not_associated(k) = i
                  k = k +1
               endif
         enddo
        j = 1
        do i=1, mol
           if (sh_env%adiab_nace_env%map_new2old(i).EQ.0) then
              sh_env%adiab_nace_env%map_new2old(i) = sh_env%adiab_nace_env%not_associated(j)
              sh_env%adiab_nace_env%map_old2new( sh_env%adiab_nace_env%map_new2old(i) ) = i
              j = j+1
           endif
        enddo
        IF (ANY(sh_env%adiab_nace_env%map_new2old.NE.sh_env%adiab_nace_env%basic) .and. print_here) THEN
           print*, "FOBSH| New map final:", sh_env%adiab_nace_env%map_new2old
           print*, "FOBSH| New map final, from old to new:", sh_env%adiab_nace_env%map_old2new
        ENDIF
       ENDIF
       IF (SUM(sh_env%adiab_nace_env%map_new2old).NE.( mol*(mol+1)/2 )) THEN
           CPABORT("Problem in the mapping between old and new states.")
       ENDIF
!> End of 2.

end subroutine nace2

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 3.
!**********************************************************************************
subroutine nace3(sh_env)

        TYPE(sh_env_type)                         :: sh_env

        CHARACTER(len=*), PARAMETER :: routineN ='nace3', &
                    routineP = moduleN//':'//routineN

        ! 3. Decide map
        IF (.NOT.(sh_env%do_reorder)) THEN
           sh_env%adiab_nace_env%map_new2old = sh_env%adiab_nace_env%basic
           sh_env%adiab_nace_env%map_old2new = sh_env%adiab_nace_env%basic
        ENDIF
        ! End of 3        

end subroutine nace3

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 4. This
!         does correct the phase as it would be arbitrary otherwise.
!**********************************************************************************
subroutine nace4(sh_env, mol)

        TYPE(sh_env_type)                         :: sh_env
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
              !PRINT*, sh_env%adiab_overlap(i, sh_env%adiab_nace_env%map_new2old(i))
              if (sh_env%adiab_overlap(i, sh_env%adiab_nace_env%map_new2old(i)).lt.0) then
                 sh_env%C_tp(:,i) = - sh_env%C_tp(:, i)
                 sh_env%adiab_overlap(i,:) = - sh_env%adiab_overlap(i,:)
              endif
              sum_ = 0  
              do j=1, mol
                 sum_ = sum_ + sh_env%C_tp(j, i) * sh_env%C_t(j, i)
              enddo
              !print*, "AC_PHASE", i, sum_
              !print*, "sh_env%adiab_overlap", sh_env%adiab_overlap(i, :)
              !if (i.eq.sh_env%active_state) print*, "AC_STATE_WARNING"          
           ENDDO
        ENDIF
        ! End of 4
end subroutine nace4


!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 5.
!**********************************************************************************
subroutine nace5(sh_env)
   
        TYPE(sh_env_type)                         :: sh_env

        CHARACTER(len=*), PARAMETER :: routineN ='nace5', &
                    routineP = moduleN//':'//routineN

        !> 5. Reassociate state active_state
        sh_env%active_state = sh_env%adiab_nace_env%map_old2new(sh_env%old_state)
        print*, "FOBSH| Put former active state", sh_env%old_state, "in new active state", sh_env%active_state 
        ! End of 5

end subroutine nace5

!**********************************************************************************
!> \brief Transform the diabatic NACE and apply the map.
!> Note that the diab_D matrix refers to time T whereas the U matrices refer to
!time Tp instead (this incosistency needs to be investigated further even though
!it will probably not afftect the results as UDU is a nagligible term)
!**********************************************************************************
subroutine nace6(sh_env, mol, deltat)

    TYPE(sh_env_type)                         :: sh_env
    INTEGER                                    :: mol
    REAL(KIND=dp)                              :: deltat, inv_deltat
 
    CHARACTER(len=*), PARAMETER :: routineN ='nace6', &
                routineP = moduleN//':'//routineN
 
    INTEGER                                    :: i,j
 
    inv_deltat = 1.0/deltat
    ! 6. Calculate the full NACE, Eq. (19) in paper 1.
    DO i=1, mol
       DO j=1, mol
          sh_env%adiab_nace_tp(i,j) = - sh_env%adiab_overlap(i, sh_env%adiab_nace_env%map_new2old(j))
          IF (i ==j) THEN 
             sh_env%adiab_nace_tp(i,j) = sh_env%adiab_nace_tp(i,j) + 1
          ENDIF
       ENDDO
    ENDDO
    ! End of 6

   if (sh_env%no_diab_nace) then 
     ! the adiabatic nace are just d = U^{dagger}Udot
     sh_env%adiab_nace_tp(:,:) = inv_deltat*sh_env%adiab_nace_tp(:,:)
   else
     !NB read description about time inconsistency
     !1. T1 <= D*U
     CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
          sh_env%adiab_nace_env%diab_D, mol, &    
          sh_env%C_tp,                  mol, &     ! C_tp = U
          0.0_DP, sh_env%adiab_nace_env%T1, mol) 

      !U^{dagger}U^dot = (1/dt)*adiab_nace_tp
      !2. adiab_nace_tp <= U^{dagg}*T1 + (1/dt)*adiab_nace_tp
      CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
          sh_env%C_tp,                  mol, &     ! C_tp = U
          sh_env%adiab_nace_env%T1,   mol, &    
          inv_deltat, sh_env%adiab_nace_tp, mol) 
   end if 
     
end subroutine nace6

!**********************************************************************************
!> \brief I've divided the calc_adiab_nace subroutine into discrete chunks in
!         order to profiling it with a finer resolution. This is chunk 7.
!**********************************************************************************
subroutine nace7(sh_env, mol) 
   
        TYPE(sh_env_type)                         :: sh_env
        INTEGER                                    :: mol

        CHARACTER(len=*), PARAMETER :: routineN ='nace7', &
                    routineP = moduleN//':'//routineN

        INTEGER                                    :: i,j

        ! 7. If reordering, we use the anti-symmetry of the NACE
        IF (sh_env%do_reorder) THEN 
           Do i=1, mol
              Do j=1, mol  
                 !print*, i,j, sh_env%adiab_nace_tp(i,j), sh_env%adiab_nace_tp(j,i)
                 !print*, i,j, ( sh_env%adiab_nace_tp(i,j) + sh_env%adiab_nace_tp(j,i) )
                 IF (j.lt.i) THEN 
                     sh_env%adiab_nace_tp(i,j) = - sh_env%adiab_nace_tp(j,i)
                 ELSE
                     sh_env%adiab_nace_tp(i,j) = ( sh_env%adiab_nace_tp(i,j) - sh_env%adiab_nace_tp(j,i) ) / 2
                 ENDIF
              Enddo
           Enddo
           !Do i=1, mol
           !   Do j=1, mol 
           !      print*, i,j, sh_env%adiab_nace_tp(i,j), sh_env%adiab_nace_tp(j,i)
           !   Enddo
           !Enddo
        ENDIF
        ! End of 7

end subroutine nace7

!**********************************************************************************
!> \brief 
!**********************************************************************************
subroutine calc_adiab_nace(sh_env, adiab_section)
        !this subroutine calculates the adiabatic nace used in the calculation
        !of the hopping probability and also in the RK adiabatic propagation 
        TYPE(sh_env_type)                         :: sh_env
        TYPE(section_vals_type), POINTER           :: adiab_section

        CHARACTER(len=*), PARAMETER :: routineN ='calc_adiab_nace', &
                                       routineP = moduleN//':'//routineN

        INTEGER                                    :: mat
        LOGICAL                                    :: print_here

        mat = size(sh_env%H_tp,1)

        CALL section_vals_val_get(adiab_section,"PRINT_MORE", &
                                                   l_val=print_here)

        if (.NOT.sh_env%no_diab_nace) then
            CALL nace1_time_tp(sh_env, mat, sh_env%deltat)
        end if
        !CALL nace1(sh_env, mat, sh_env%deltat)
        CALL nace2(sh_env, mat, print_here)
        CALL nace3(sh_env)
        CALL nace4(sh_env, mat)
        CALL nace5(sh_env)
        CALL nace6(sh_env, mat, sh_env%deltat)
        CALL nace7(sh_env, mat)
 
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

!**********************************************************************************
!> \brief 
!**********************************************************************************
   SUBROUTINE read_list_frz_couplings(ninteractions, listcoeff, h1, frz_couplings)

    TYPE(cp_sll_val_type), POINTER                         :: listcoeff
    INTEGER, DIMENSION(:,:), POINTER, INTENT(OUT)        :: h1
    !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: h2
    REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(OUT)  :: frz_couplings
    INTEGER, INTENT(IN)                                    :: ninteractions

    CHARACTER(len=*), PARAMETER :: routineN = 'read_listcoeff', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)                 :: line_att
    CHARACTER(len=default_string_length)                   :: error_message
    INTEGER                                                :: default_id, end_c, &
                                                              iatom, j,  &
                                                              start_c, wrd
    LOGICAL                                                :: is_ok
    TYPE(val_type), POINTER                                :: val
    REAL(KIND=dp)                                          :: my_pres

    NULLIFY(val)

    ALLOCATE(h1(ninteractions,2))
    ALLOCATE(frz_couplings(ninteractions))

    !this reads line by line the section listcoeff, and for each line it has
    !the different CASES 
       DO iatom=1,ninteractions
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(listcoeff,val) !logical 
          CALL val_get(val,c_val=line_att)
          default_id = str2id(s2s(""))
          ! Read name and atomic coordinates
          start_c=1
          DO wrd=1,3
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
             !IF(LEN_TRIM(line_att(start_c:end_c-1)) == 0)& !length of the number
             !                                               ! without space
             !   CALL cp_abort(__LOCATION__,&
             !        "Incorrectly formatted input line for atom "//&
             !        TRIM(ADJUSTL(cp_to_string(iatom)))//&
             !        " found in COORD section. Input line: <"//&
             !        TRIM(line_att)//"> ")
             SELECT CASE (wrd)
             CASE (1)  ! read active or not atom (carrying a charge)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_pres,error_message)
                h1(iatom,1) = INT(my_pres)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (2)  ! read active or not atom (carrying a charge)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_pres,error_message)
                h1(iatom,2) = INT(my_pres)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (3) ! transition charge (this must be divided by sqrt(2) if it
                      !comes from gaussian fitting of the electrostatic potential 
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       frz_couplings(iatom),error_message)
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

   END SUBROUTINE read_list_frz_couplings

!**********************************************************************************
!> \brief 
!**********************************************************************************
   SUBROUTINE read_list_tresp(natom_mol, listcoeff, ext_elements, ext_pres, tresp_charges)

    TYPE(cp_sll_val_type), POINTER                         :: listcoeff
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: ext_elements
    LOGICAL, DIMENSION(:), POINTER,  INTENT(OUT)           :: ext_pres
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: tresp_charges
    INTEGER, INTENT(IN)                                    :: natom_mol

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

    ALLOCATE(ext_pres(natom_mol))
    ALLOCATE(ext_elements(natom_mol))
    ALLOCATE(tresp_charges(natom_mol))

    !this reads line by line the section listcoeff, and for each line it has
    !the different CASES 
       DO iatom=1,natom_mol
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
             CASE (2)  ! read atomic number 
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       my_ext_elements,error_message)
                ext_elements(iatom) = INT(my_ext_elements)
                IF(LEN_TRIM(error_message) /= 0)&
                   CALL cp_abort(__LOCATION__,&
                        "Incorrectly formatted input line for atom "//&
                        TRIM(ADJUSTL(cp_to_string(iatom)))//&
                        " found in COORD section. "//TRIM(error_message)//&
                        " Input line: <"//TRIM(line_att)//"> ")
             CASE (3)  ! read active or not atom (carrying a charge)
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
             CASE (4) ! transition charge (this must be divided by sqrt(2) if it
                      !comes from gaussian fitting of the electrostatic potential 
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       tresp_charges(iatom),error_message)
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

   END SUBROUTINE read_list_tresp

!******************************************************************************
!!!!!!!!!!!!!!   DEPRECATED SUB. BELOW  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

!**********************************************************************************
!> \brief Calculate adiabatic forces according to Spencer 2016 and Carof 2019
!for each atom in the system on the state passed as input.
! COMMENT_SG: I generalized this subroutine and kept this just for testing
! purposes
!**********************************************************************************
subroutine calcforces_aom(state, Hin, C, nacv, atoms, scaling, cbar, Fin, forces_adiab, & 
                  ext_pres)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom
        integer, intent(in)                         :: state
        REAL(KIND=dp), intent(in), dimension(:,:)   :: Hin, C
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: nacv
        integer, intent(in)                         :: atoms !this is atoms per mol
        REAL(KIND=dp), intent(in)                   :: scaling, cbar
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), dimension(:,:), POINTER      :: forces_adiab
        LOGICAL, dimension(:), INTENT(IN)           :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN ='calcforces_aom', &
                               routineP = moduleN//':'//routineN

        integer                                     :: i, j, mat, d
        INTEGER                                     :: natom
        INTEGER                                     :: num_mol, num_j
        REAL(KIND=dp)                               :: convert_nacv 

        natom = SIZE(Fin,2)
        !nadiab= SIZE(Fin,3)
        mat = size(Hin,1) ! COMMENT_SG: In general this array could have
                            !degenerate states (becuase the construction of AOM
                            !Hamiltonian can be done with multiple orbitals per
                            !site), however the following subroutine cannot
                            !handle the degenerate case becuase Fin(d,j,i) is
                            !the array of diagonal forces and its dimension
                            !stops at the number of total force_evals
                            !considered. To sum up if you want this subroutine
                            !to be able to handle degenerate cases you should
                            !add a further loop over "orb" like is done for site
                            !energies of the Hamiltoninan. This is the best
                            !strategy becuase it would not even require to
                            !change the dimension of Fin(d,j,i) that is one of
                            !the biggest arrays in the code.

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
end subroutine calcforces_aom

!**********************************************************************************
!> \brief Calculate the diabatic NACE first
! DEPRECATED BECAUSE TIME INCONSISTENT
!**********************************************************************************
subroutine nace1(sh_env, mol, deltat)

    TYPE(sh_env_type)                         :: sh_env
    INTEGER                                    :: mol, i, j, sub, add
    REAL(KIND=dp)                              :: deltat

    CHARACTER(len=*), PARAMETER :: routineN ='nace1', &
                routineP = moduleN//':'//routineN

    !> Based on Eq. 19 of paper: Spencer et al., J. Chem. Phys. 145, 064102 (2016)

    !> 1. Calculate the term  D  (first term in Eq. 19)
    !many allocations
    !use the matinvsqrt subroutine to find T = S^(-1/2)
    sh_env%adiab_nace_env%T1 = sh_env%Stmtm
    sh_env%adiab_nace_env%T2 = sh_env%Stptp
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T1)
    call matinvsqrt(mol, .TRUE., sh_env%adiab_nace_env%T2)
    !and do some straightforward first-order numerical derivatives
    !Time deriv of overlap
    sh_env%adiab_nace_env%DT = (0.5D0/deltat)*(sh_env%adiab_nace_env%T2 - sh_env%adiab_nace_env%T1)
    !Time deriv of Lowdin Mat
    sh_env%adiab_nace_env%DS = (0.5D0/deltat)*(sh_env%Sttp - sh_env%Sttm)

!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
!PRINT *, "TST Bef:"
!DO i=1, mol
!   DO j=1, mol
!      WRITE(6, "(A2, e12.4, A2)", advance='no') "| ", sh_env%adiab_nace_env%TST(i,j), " "
!   END DO
!   PRINT *, " "
!END DO

    !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
    !by transforming from the diabatic overlaps
    !TST = T1*Stt*DT. This is done in 2 steps.
    !Note that T2 is overwritten in the following DGEMM operation 
    !1. T2 = Stt*DT
    CALL DGEMM('N', 'N', mol, mol, mol, 1.0D0, &
                sh_env%Stt,               mol, &
                sh_env%adiab_nace_env%DT, mol, &
                0.0_DP, sh_env%adiab_nace_env%T2, mol)

    !2. TST = T1^{dagg}*T2
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                sh_env%adiab_nace_env%T1,  mol, &
                sh_env%adiab_nace_env%T2, mol, &       ! at this point diab_D = TST
                0.0_DP, sh_env%adiab_nace_env%diab_D, mol)

    ! This part takes advantage of the sparse nature of the DS matrix
    !IS ONLY VALID FOR CHAINS WITH ONLY OFF DIAGONALS ELEMENTS
    ! Again the mat-mat mult is done in 2 steps:
    !1. T2 = DS*T1
   DO j=1, mol
     sh_env%adiab_nace_env%T2(1,j) = &
      sh_env%adiab_nace_env%DS(1,2)*sh_env%adiab_nace_env%T1(2,j)
   END DO

   DO i=2, mol-1
      DO j=1, mol
        sub = i-1
        add = i+1
        sh_env%adiab_nace_env%T2(i,j) = &
          sh_env%adiab_nace_env%DS(i,sub)*sh_env%adiab_nace_env%T1(sub,j) &
        + sh_env%adiab_nace_env%DS(i,add)*sh_env%adiab_nace_env%T1(add,j)
      END DO
   END DO
   sub = mol-1
   DO j=1, mol
     sh_env%adiab_nace_env%T2(mol,j) = &
      sh_env%adiab_nace_env%DS(mol,sub)*sh_env%adiab_nace_env%T1(sub,j)
   END DO


!   ! If the loop above doesn't work then try using this.
!    !1. T2 = DS*T1
!    CALL DGEMM('N', 'N', mol, mol, mol, 1.0_DP, &
!               sh_env%adiab_nace_env%DS, mol, &
!               sh_env%adiab_nace_env%T1, mol, &
!               0.0_DP, sh_env%adiab_nace_env%T2, mol)

    !2. TDT = T1^{dagg}*T2 + TST
    CALL DGEMM('T', 'N', mol, mol, mol, 1.0_DP, &
                sh_env%adiab_nace_env%T1,  mol, &
                sh_env%adiab_nace_env%T2, mol, &       ! Now diab_D becomes = TDT + TST
                1.0_DP, sh_env%adiab_nace_env%diab_D, mol)
!   THE FOLLOWING CODE IS USEFUL FOR PRINTING NxN MATRICIES       
PRINT *, "D Bef:"
DO i=1, mol
   DO j=1, mol
      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", sh_env%adiab_nace_env%diab_D(i,j), "  "
   END DO
   PRINT *, " "
END DO
end subroutine nace1

subroutine calc_adiab_nace_old(sh_env, adiab_section)
        !this subroutine calculates the adiabatic nace used in the calculation
        !of the hopping probability and also in the RK adiabatic propagation 
!WARNING: THIS SUBROUTINE HAS TIME INCONSISTENCIES IN THE TIME DERIVATIVES THAT
!HAVE BEEN CORRECTED IN THE CHOPPED VERSION ABOVE
        TYPE(sh_env_type)                         :: sh_env
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

        mol = size(sh_env%Stt,1)
        deltat = sh_env%deltat

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
        T1 = sh_env%Stmtm
        T2 = sh_env%Stptp
        call matinvsqrt(mol, .TRUE., T1)
        call matinvsqrt(mol, .TRUE., T2)
        !and do some straightforward first-order numerical derivatives
        DT = (0.5D0/deltat)*(T2 - T1)
        DS = (0.5D0/deltat)*(sh_env%Sttp - sh_env%Sttm)
        !with these two loops, calculate the nacv in the time domain between adiabatic wavefunctions
        !by transforming from the diabatic overlaps
        do i=1,mol
                do j=1,mol
                        TST(i,j) = 0.0D0
                        TDT(i,j) = 0.0D0
                        do k=1,mol
                                do l=1,mol
                                        TST(i,j) = TST(i,j) + T1(k,i)*sh_env%Stt(k,l)*DT(l,j)
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
              sh_env%adiab_overlap(i,j) = 0.0
              do k = 1, mol
                 sh_env%adiab_overlap(i,j) = sh_env%adiab_overlap(i,j) + sh_env%C_tp(k,i)*sh_env%C_t(k,j)
              enddo
              if (abs(sh_env%adiab_overlap(i,j)).GE.abs(max_val)) THEN
                 max_val = sh_env%adiab_overlap(i,j)
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
                            "with overlap", max_val, "and former", sh_env%adiab_overlap(i,i)
                 else 
                    print*, "FOBSH! Reorgering: for state", i, " we stay in with overlap", &
                            sh_env%adiab_overlap(i,i), " instead of goind to ", new_index, "with overlap", max_val
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
               IF ( abs( (sh_env%adiab_overlap(i, map_old2new(i)))**2 - 1 ).LT.0.5 ) THEN
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
       IF (.NOT.(sh_env%do_reorder)) THEN
          map_new2old = basic
          map_old2new = basic
       ENDIF
! End of 3        

! 4. Change sign to keep constant convention
       do_sign_change = .True.
       IF (do_sign_change) THEN
          DO i=1, mol
             !PRINT*, sh_env%adiab_overlap(i, map_new2old(i))
             if (sh_env%adiab_overlap(i, map_new2old(i)).lt.0) then
                sh_env%C_tp(:,i) = - sh_env%C_tp(:, i)
                sh_env%adiab_overlap(i,:) = - sh_env%adiab_overlap(i,:)
             endif
             sum_ = 0 
             do j=1, mol
                sum_ = sum_ + sh_env%C_tp(j, i) * sh_env%C_t(j, i)
             enddo
             !print*, "AC_PHASE", i, sum_
             !print*, "sh_env%adiab_overlap", sh_env%adiab_overlap(i, :)
             !if (i.eq.sh_env%active_state) print*, "AC_STATE_WARNING"          
          ENDDO
       ENDIF
! End of 4


!> 5. Reassociate state active_state
        sh_env%active_state = map_old2new(sh_env%old_state)
        print*, "FOBSH| Put former active state", sh_env%old_state, "in new active state", sh_env%active_state 
! End of 5

! 6. Calculate the full NACE, Eq. (19) in paper 1.
        DO i=1, mol
           DO j=1, mol
              sh_env%adiab_nace_tp(i,j) = - sh_env%adiab_overlap(i, map_new2old(j))
              IF (i ==j) THEN
                 sh_env%adiab_nace_tp(i,j) = sh_env%adiab_nace_tp(i,j) + 1
              ENDIF
              sh_env%adiab_nace_tp(i,j) = sh_env%adiab_nace_tp(i,j) / deltat
              DO k=1, mol
                 DO l=1, mol
                    sh_env%adiab_nace_tp(i,j) = & 
                                    sh_env%adiab_nace_tp(i,j) + sh_env%C_tp(k,i)*sh_env%C_tp(l,j)*(TST(k,l) + TDT(k,l))
                 ENDDO
              ENDDO  
           ENDDO
        ENDDO
! End of 6

! 7. If reordering, we use the anti-symmetry of the NACE
        IF (sh_env%do_reorder) THEN
           Do i=1, mol
              Do j=1, mol 
                 !print*, i,j, sh_env%adiab_nace_tp(i,j), sh_env%adiab_nace_tp(j,i)
                 !print*, i,j, ( sh_env%adiab_nace_tp(i,j) + sh_env%adiab_nace_tp(j,i) )
                 IF (j.lt.i) THEN
                     sh_env%adiab_nace_tp(i,j) = - sh_env%adiab_nace_tp(j,i)
                 ELSE
                     sh_env%adiab_nace_tp(i,j) = ( sh_env%adiab_nace_tp(i,j) - sh_env%adiab_nace_tp(j,i) ) / 2
                 ENDIF
              Enddo
           Enddo
           !Do i=1, mol
           !   Do j=1, mol 
           !      print*, i,j, sh_env%adiab_nace_tp(i,j), sh_env%adiab_nace_tp(j,i)
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
                                                     sh_env, adiab_section)

        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface
        REAL(KIND=dp), intent(in) :: deltat
        integer, intent(in) :: old_state
        REAL(KIND=dp), dimension(:,:), intent(in) :: C2
        INTEGER, INTENT(OUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: attempt_hop
        REAL(KIND=dp), INTENT(IN) :: decohere_criterion
        TYPE(sh_print_env_type), INTENT(OUT) :: sh_print_env
        TYPE(sh_env_type)                   :: sh_env
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
        
        mol = size(sh_env%Stt,1)

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
            CALL old_diab_to_adiab(sh_env, C2, old_state, mol, cstate, ct)
        CASE("ADIABATIC_BASIS")
            ct(:) = sh_env%ctin(:)
            cstate = sh_env%ctin_t(old_state)
        CASE DEFAULT
           CPABORT("There is no method to propagate the wavefunction.")
        END SELECT

!NOT YET READY FOR ADIAB PRPAGATION

        do i=1,mol
                ak(i) = ct(i)*conjg(cstate)
                if (sh_env%HS_propagation .and. &
                              (i .eq. sh_env%class_forbidden_hop_t(i))) then
                   b(i) = 0.0D0 
                else 
                b(i) = -2.0D0*REAL(conjg(ak(i))*sh_env%adiab_nace_tp(i,old_state))
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

        hop = next_random_number(sh_env%rng_stream)

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
                                sh_env%decohere = .FALSE.
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
                                sh_env%decohere = .TRUE.
                                do j=1,mol
                                        if (abs(sh_env%adiab_nace_tp(old_state,j)) & 
                                                 > decohere_criterion .and. old_state .ne. j) then
                                               ! decohere = .FALSE.
                                                sh_env%decohere = .FALSE.
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

SUBROUTINE old_adiab_overlap_check(sh_env, changed_state, itimes)
    TYPE(sh_env_type)                         :: sh_env
    INTEGER, OPTIONAL, intent(in)               :: changed_state, itimes  

    CHARACTER(len=*), PARAMETER :: routineN = 'old_adiab_overlap_check', &
              routineP = moduleN//':'//routineN

    REAL(KIND=dp), dimension(:), allocatable    :: overlaps
    integer                                     :: mat, i, j, res_mod, state, &
                                                     past_state
    REAL(KIND=dp)                               :: sum_ov, max_val, resid 
    LOGICAL :: overlap_over_thresh
!!ADD OPTIONAL

    mat    = SIZE(sh_env%Stptp, 1)

    allocate(overlaps(2*mat))

    overlaps = 0.0D0
    
    
    if (PRESENT(changed_state)) then
        past_state = changed_state
!print *, "OLD_STATE", changed_state
    else
        past_state = sh_env%active_state
!print *, "STATE_CHANGED", past_state
!print *, "NEW_STATE_here", sh_env%active_state    
    end if

!print *, "NEW_STATE", sh_env%active_state    

     do i=1, mat
         overlaps(1) = overlaps(1) + sh_env%C_t(i, past_state)*sh_env%C_tp(i, sh_env%active_state)
     end do
!print *, "OVERLAP_NOTSQUARED", overlaps(1)

     overlaps(1) = abs(overlaps(1))*abs(overlaps(1))
print *, "OVERALP_SQUARED", past_state, sh_env%active_state, overlaps(1)

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
           state = sh_env%active_state+res_mod
           !two conditions in case state is negative or is larger than mat
           if (state .le. 0 .or. state .gt. mat) then
!              print *, "skipped", j, state
              cycle
           else  
              do i = 1, mat
                 overlaps(j) = overlaps(j) + sh_env%C_t(i, past_state)*&
                                             sh_env%C_tp(i, state)
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
print *, "LOCALIZATION", maxloc(abs(sh_env%C_t(:, past_state)*&
                            sh_env%C_t(:, past_state))), &
      maxloc(abs(sh_env%C_tp(:,sh_env%active_state)*sh_env%C_tp(:,sh_env%active_state)))
print *, "STATE_TOGOTO", past_state, sh_env%active_state
else
print *, "LOCALIZATION", maxloc(abs(sh_env%C_t(:, past_state)*&
                            sh_env%C_t(:, past_state)))
!            maxloc(abs(sh_env%C_tp(:,sh_env%enforced_active_state)*&
!                            sh_env%C_tp(:,sh_env%enforced_active_state)))
!print *, "STATE_TOGOTO", past_state, sh_env%enforced_active_state
end if

    deallocate(overlaps)
     
END SUBROUTINE old_adiab_overlap_check

END MODULE aom_utils
