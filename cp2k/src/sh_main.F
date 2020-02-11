!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2016.01.15: Main subroutine for FSSH Hamiltonian creations and
!>                          propagation of electronic equation of motions
!> Edited SG
! *****************************************************************************

! *****************************************************************************
! \brief below a very OLD comment of Jacob (who was the PhD first to write the code) 
!  the code has hugely changed since then:
!
!>   jacob: OLD comment 
!    propagate_namd_all_steps_full_degen
!
!    This program is designed to propagate a single excess electronic charge
!    along a molecular dynamics trajectory, either during the trajectory
!    (to do surface hopping simulations) or afterwards (classical path approximation).
!
!    This program takes molecular coordinates as inputs.
!    From these coordinate inputs, it calculates the electronic Hamiltonian
!    in a molecular orbital basis, using the overlapFINAL module to calculate
!    the off-diagonals and either accepting the diagonals as inputs or calculating
!    them in some other way.
!
!    With the electronic Hamiltonian, the program then propagates the electron using the RK4 method.
!    At the final electronic state, the probabilities for hopping between adiabatic surfaces are calculated.
!    The program then decides what adiabatic surface the system is on and calculates
!    the appropriate atomic forces.
!
!    In the classical path approximation, the above steps are repeated until the
!    end of the input molecular trajectory.  For surface-hopping, the data here is fed back into
!    the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE sh_main 

     USE kinds,                        ONLY: default_string_length,&
                                             dp


     USE aom_utils,                    ONLY: makeX,                     &
                                             calc_hop_nacv_aom,         &
                                             adiabaticity,              &
                                             calc_aom_Xmatrix_elements_first_step, & 
                                             calc_aom_Xmatrix_elements, &
                                             make_adiabX,               &
                                             calcforces_aom,                &
                                             read_listcoeff,            &
                                             calc_adiab_nace,           &
                                             calc_adiab_nace_old

     USE sh_types,                    ONLY: sh_env_type, &
                                                sh_print_env_type

     USE mixed_energy_types,           ONLY: mixed_force_type
     USE input_section_types,          ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_list_get,&
                                              section_vals_val_get

     USE cp_linked_list_val,           ONLY: &
         cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
         cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
         cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
         cp_sll_val_set_el_at, cp_sll_val_type

     USE aom_print,                    ONLY: adiab_print_coeff,              &
                                             adiab_print_adiab_pop,                &
                                             adiab_print_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_site_energies,   &
                                             adiab_print_off_diagonals,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian, &
                                             adiab_print_pseudo_hamiltonian, &
                                             adiab_print_nacv, &
                                             adiab_print_pvecs



     USE atomic_kind_list_types,       ONLY: atomic_kind_list_type

     USE atomic_kind_types,            ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

     USE distribution_1d_types,        ONLY: distribution_1d_type

     USE non_adiabatic_utils,          ONLY: decohere_wf,                     &
                                             calc_coms_and_nlist,             &
                                             calc_frozen_H,                   &
                                             calcforces_frozenH,              &
                                             diab_to_adiab,                   & 
                                             surf_hop_type,                   &
                                             choose_firstadiabat,             &
                                             print_matrix_standard_out,       & 
                                             translate_com,                   &
                                             propagate_wf,                    &
                                             propagate_wf_fast,               &
                                             propagate_adiab_wf,              &
                                             propagate_wf_diab_nace_interpol, &
                                             prop_wf_HS_corr,                 & 
                                             classically_forbidden_hops,      &
                                             define_active_region,            &
                                             calc_psi_after_correction

     USE parallel_rng_types,           ONLY: UNIFORM,&
                                          create_rng_stream,&
                                          delete_rng_stream

     USE cp_log_handling,              ONLY: cp_get_default_logger, &
                                             cp_logger_type

     USE cp_output_handling,           ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_should_output


#include "./base/base_uses.f90"


     IMPLICIT NONE

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='sh_main'   

        PUBLIC :: propagate_namd_all_steps, & 
                  propagate_namd_first_step

CONTAINS

!**********************************************************************************
!> \brief  This construct the Hamiltonian and Initialize electronic coefficients
!at the first time step. Note that it could be merged with the subroutine for
!the remaining step in future (there is no particular reason why should be
!separated from that, except for the fact that it was easier to implement like
!this) 
!> authors AC and SG
!> edited 02/02/20
!**********************************************************************************
  SUBROUTINE propagate_namd_first_step(nadiab, norbitals,  &
                                  positions_tp,  &
                                  velocities_tp,  &
                                  energies_tp,  &
                                  sh_env, itimes, &
                                  adiab_section, my_bo_keyword,  &
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword, method_keyword, diagonal_F)

     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: velocities_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, POINTER                               :: itimes
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section,optimization_section
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     CHARACTER(LEN=default_string_length)           :: method_keyword

     CHARACTER(len=*), PARAMETER :: routineN ='propagate_namd_first_step', &
     routineP = moduleN//':'//routineN

     integer                                        :: coord, mol,  & 
                                                       first_diabat, active_state, mat, orb 
     REAL(KIND=dp)                                  :: deltatMD, deltatE
     COMPLEX(KIND=dp), dimension(:), allocatable    :: ctin
     INTEGER                                        :: new_state , ELsteps   
     REAL(KIND=dp)                                  :: times 
     INTEGER                                        :: imol,  natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE     :: exact_forces
     LOGICAL                                        :: keyword_translate
     TYPE(sh_print_env_type)                        :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)           :: keyword_init, keyword_basis_repr, &
                                                        prop_wf_keyword
     real(kind=dp), DIMENSION(:,:), allocatable     :: PC_allrvecs
     REAL(kind=dp), DIMENSION(:), ALLOCATABLE       :: site_energies, couplings
     REAL(KIND=dp), dimension(:), allocatable       :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable    :: adiab_c
     TYPE(cp_logger_type), POINTER                  :: logger
     INTEGER                                        :: unit_nr
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F
     TYPE(section_vals_type), POINTER               :: coupling_section
     CHARACTER(LEN=default_string_length)           :: coupling_keyword
     !INTEGER                                        :: i, k

     NULLIFY(positions_com_tp, velocities_com_tp)
     NULLIFY(logger)
     NULLIFY(optimization_section)
     NULLIFY(coupling_section)

     logger => cp_get_default_logger()

     CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)

     !> Calculate the orbital overlap in the COM reference sistem of the
     !molecules
     CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS", l_val=keyword_translate)

     CALL section_vals_val_get(adiab_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ELsteps)

     !atoms = sh_env%atoms_per_site !Number of atoms per site
     natom = sh_env%natom !Total number of atoms
     deltatMD = sh_env%deltat !Nuclear timestep

     !> Number of sites
     mol = nadiab ! CHANGE_AC
     !> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
     !> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
     !> electronic timestep is a fraction of the molecular timestep
     deltatE = deltatMD/ELsteps !CHANGE_AC

     sh_print_env%adiab_section => adiab_section
     sh_print_env%nadiab = mat
     sh_print_env%norbitals = orb
     sh_print_env%deltatMD = deltatMD

     ALLOCATE(sh_print_env%adiab_energies_t(mat))
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))
     allocate(ctin(mat))      !ctin goes into the subroutine ZGEMV
                              !these are carrier wf coefficient in the diab
                              !basis

     !READ THE COORDINATES, ENERGIES AND STATE
     IF (method_keyword.eq."FROZEN_HAMILTONIAN") THEN
        velocities_com_tp = 0.0D0
     ELSE
        CALL translate_com(velocities_tp, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     ENDIF
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  

print *, "positions_tp", positions_tp(:,1)
print *, "positions_tp", positions_tp(:,2)
print *, "positions_tp", positions_tp(:,3)

print *, "positions_com_tp", positions_com_tp(:,1)
print *, "positions_com_tp", positions_com_tp(:,2)
print *, "positions_com_tp", positions_com_tp(:,3)

     IF (sh_env%aom_couplings_keyword) THEN
         !------------------  AOM section -------------------------------------
         ALLOCATE(PC_allrvecs(5,mol*sh_env%atoms_per_site))
         do coord=1,3
            sh_env%crds(coord,:,1)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
            !At time 0, we must guess the positions at time -1 with the inital velocities
            sh_env%crds(coord,:,2)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres) -  &
                              PACK(velocities_com_tp(coord,:), sh_env%ext_pres)*deltatMD
            sh_env%crds(coord,:,3)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
            sh_env%siteenergies(2,:) = 0.0D0 ! CHANGE_AC
            sh_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
         enddo
print *, "sh_env%ext_pres", sh_env%ext_pres
print *, "SIZE sh_env%ext_pres", SIZE(sh_env%ext_pres)         
print *, "positions_com_tp(coord,:)", SIZE(positions_com_tp(1,:))
         !> CALCULATE THE CENTER OF MASS AND CONSTRUCT NEIGHBOURS LIST
         !> COMMENT_SG: it assumes same number of atoms per mol
         CALL calc_coms_and_nlist(mol, sh_env%cutoff_sites, orb, crds=sh_env%crds, &
                                    active_system_info=sh_env%active_system_info, NN=sh_env%NN)

         ! Allocate neighbouring list just at the first step as Jacob suggested 
         !only if sh_env%multiple_time_step .eq. 1 otherwise we use Multiple time
         !step feature and to save memory we recalculate bneighbour and aneighbour
         !every N steps. This is done also for the NACV giving a huge speed-up. NOTE
         !that in general is not convenient to recalculate neighboring list at each
         !step because is to expensive to do and not useful.
         if (sh_env%multiple_time_step .eq. 1) then 
           allocate(sh_env%neighbourlist(mol,mol,0:sh_env%atoms_per_site,0:sh_env%atoms_per_site))
           sh_env%neighbourlist = 0
         end if 

         !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
         CALL calc_aom_Xmatrix_elements_first_step(sh_env, sh_env%atoms_per_site, mol, sh_env%dR, sh_env%do_random, sh_env%cbar, &
                          sh_env%scaling, PC_allrvecs, orb, sh_env%minao)

print *, "sh_env%cutoff_sites",sh_env%cutoff_sites
print *, "NN", sh_env%NN 

         !> CONVERT NACV AND CALCULATE EIGENPROPERTIES 
         !> THESE ASSOCIATIONS ARE NEEDED FOR THE NACV CALC
         sh_env%nacv = sh_env%cbar*sh_env%nacv
         !-----------------------------------------------------------------------------
     ELSE IF (sh_env%frozen_couplings_keyword) THEN 
         !> FROZEN COUPLINGS 
         CALL calc_frozen_H(sh_env, energies_tp, orb)

     ELSE IF (sh_env%tresp_couplings_keyword) THEN 
         !> TRESP COUPLINGS
         do coord=1,3
            sh_env%crds(coord,:,1)   = 0.0D0 !time tm not used
            sh_env%crds(coord,:,2)   = 0.0D0 !time t not used
            sh_env%crds(coord,:,3)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
            !sh_env%crds(coord,:,3)   = 0.0D0
         enddo
         
print *, "sh_env%ext_pres", sh_env%ext_pres
print *, "SIZE sh_env%ext_pres", SIZE(sh_env%ext_pres)         
print *, "positions_com_tp(coord,:)", SIZE(positions_com_tp(1,:))
         !> CALCULATE THE CENTER OF MASS AND CONSTRUCT NEIGHBOURS LIST
         !this is needed to calculate the NN list
         CALL calc_coms_and_nlist(mol, sh_env%cutoff_sites, orb, crds=sh_env%crds, &
                                    active_system_info=sh_env%active_system_info, NN=sh_env%NN)

print *, "sh_env%cutoff_sites",sh_env%cutoff_sites
print *, "NN", sh_env%NN 

         CALL calc_frozen_H(sh_env, energies_tp, orb)

     ELSE
         !> In a possible extension one should add here 
         CPABORT("Another coupling method to constrcut H should be chosen")
     END IF

!!TEST FROZEN 
!    do i =1, size(sh_env%H_tp,1)
!      do k = 1, size(sh_env%H_tp,1)
!         if (i .ne. k) then
!            if ((k+1 == i) .or.(k-1 == i)) then
!               sh_env%H_tp(i,k) = cmplx(0.0D0,0.001102479,dp) 
!            else
!               sh_env%H_tp(i,k) = cmplx(0.0D0,0.0D0,dp)
!            end if
!         end if
!     end do
!    end do
!sh_env%nacv = 0.0D0
!!print *, "sh_env%H_tp(:,:)", sh_env%H_tp

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian H_tp
     !It gives us back the eigenvalues E_tp and eigenvectors C_tp
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(sh_env%H_tp), sh_env%C_tp, sh_env%E_tp)
     
     !> PRINT
     !> print site energies 
     allocate(site_energies(mat))
     allocate(couplings(mat-1))
     do imol=1,mat
         site_energies(imol) = -AIMAG(sh_env%H_tp(imol,imol))
         if (imol .ne. mat) then
           couplings(imol) = -AIMAG(sh_env%H_tp(imol,imol+1))
         end if
     end do
   
     times = DBLE(itimes)

     !> print info in the output only if we do not restart the file
     sh_print_env%times = times
     sh_print_env%times = DBLE(itimes)
     if (.NOT. sh_env%restart_key) then
        CALL adiab_print_site_energies(adiab_section, site_energies, nadiab,deltatMD, times) 
        CALL adiab_print_off_diagonals(adiab_section, couplings, nadiab-1, deltatMD, times)
        CALL adiab_print_adiabat_energies(adiab_section, sh_env%E_tp, nadiab, deltatMD, times)
        CALL adiab_print_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        CALL adiab_print_pseudo_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        IF (sh_env%aom_couplings_keyword) THEN
            CALL adiab_print_nacv(adiab_section, sh_env%nacv, nadiab, norbitals, sh_env%atoms_per_site, itimes, deltatMD)
            call adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, sh_env%atoms_per_site, nadiab)
        END IF
     end if

     !> INITIALIZE THE WAVEFUNCTION (WE DO NOT PROPAGATE IT AT THE FIRST STEP) 
     IF  (.NOT. my_test_keyword) THEN
         sh_env%decohere = .FALSE.
         sh_env%attempt_hop = .FALSE.
         !if we want to restart wavefunction then
         if (sh_env%restart_key) then
             keyword_init = "DIABATIC_RESTART"
         else
             CALL section_vals_val_get(adiab_section, "INITIALIZATION", c_val=keyword_init)
         end if

         CALL section_vals_val_get(adiab_section, "REPRESENTATION", c_val=keyword_basis_repr)
         SELECT CASE(keyword_init)
         CASE("DIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_DIABAT", &
                                       i_val=first_diabat)
             CALL choose_firstadiabat(sh_env%C_tp, first_diabat, sh_print_env, active_state, sh_env) 
             !an active state in needed for the nuclear propagation but the wafunction has to
             !start form an eigenvector of the hamiltonian in the diab_basis
             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(first_diabat) = CMPLX(1.0D0, 0.0D0,dp)  
             CASE("ADIABATIC_BASIS")
                 !!!!! THIS SECTION IS DEPRECATED !!!!!!
                 !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
                 CALL nacv_first_ts(sh_env,diagonal_F, natom, mat, velocities_com_tp) 
                 ctin(:) = sh_env%C_tp(first_diabat,:)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("ADIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_ADIABAT",&
                                      i_val=active_state)
             !> print init
             sh_print_env%first_adiabat    = active_state
             ALLOCATE(sh_print_env%pop_init(SIZE(sh_env%C_tp(:,active_state))))
             !all pop is initialized in the active state, this is a hack to
             !correct the printing
             sh_print_env%initrandom      = 0.0D0 ! COMMENT_SG: hack for regtests
             sh_print_env%pop_init(:)      = 0.0D0
             sh_print_env%pop_init(active_state)  = 1.0D0
             CALL adiab_print_sh_info(sh_print_env, keyword = 'START')
             DEALLOCATE(sh_print_env%pop_init)

             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = sh_env%C_tp(:, active_state)
             CASE("ADIABATIC_BASIS")
                 !!!!! THIS SECTION IS DEPRECATED !!!!!!
                 !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
                 CALL nacv_first_ts(sh_env, diagonal_F, natom, mat, velocities_com_tp) 
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(active_state) = CMPLX(1.0D0, 0.0D0,dp)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("DIABATIC_RESTART")
                 ctin(:) = sh_env%readcoeff(:)
                 active_state = sh_env%active_state_res
                 ! do not print this duplicated step in run-sh.log
                 !sh_print_env%first_adiabat = active_state 
                 !CALL adiab_print_sh_info(sh_print_env, keyword = 'START') 
         CASE DEFAULT
             CPABORT("INITIALIZATION implemented: DIABATIC, ADIABATIC")
         END SELECT

         allocate(adiab_c(mat)) 
         allocate(adiab_pop(mat))

         !> print info in the output only if we do not restart the file
         if (.NOT. sh_env%restart_key) then
             !NB: this subroutine is wasted computational work if adiab_pop is not
             !printed out, same probably for AOM.main

             unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIAB_POPULATIONS",&
               extension=".xyz", middle_name="adiab_pop",log_filename=.FALSE.)
             if (unit_nr > 0) then 
                CALL diab_to_adiab( sh_env%C_tp, ctin, adiab_c, adiab_pop)     
                CALL adiab_print_adiab_pop(adiab_section, adiab_pop, & 
                                  nadiab,norbitals,itimes,deltatMD, times)
             end if 
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, times)
             CALL adiab_print_pop(adiab_section, abs(ctin),nadiab,norbitals,itimes,deltatMD, times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes, &
                                        deltatMD, times)
         end if
     ELSE 
         sh_env%attempt_hop = .TRUE.
         active_state = 1
         new_state = mod(active_state + 1, mat) 
         IF ((new_state).EQ.0) new_state = mat
     ENDIF 


     !> ######## FORCES CALCULATION   ###########################################################
     print*, "FOBSH| ACTIVE STATE = ", active_state
     IF (sh_env%aom_couplings_keyword) THEN
        !------------------  AOM section --------------------------------------------------------
        !> second part CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE
        !RESULT OF THE SURFACE HOP ROUTINE IS.
        call calcforces_aom(state=active_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, &
                                                            nacv=sh_env%nacv, atoms=sh_env%atoms_per_site, & 
                        scaling=sh_env%scaling, cbar=sh_env%cbar, Fin=diagonal_F, forces_adiab=sh_env%forces_adiab, &
                        ext_pres=sh_env%ext_pres)

        IF (sh_env%analytics) THEN
           times = DBLE(itimes)
           CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
        ENDIF
        !-----------------------------------------------------------------------------------------------
     ELSE IF (sh_env%frozen_couplings_keyword) THEN 
        !frozen Hamiltonian case
        CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, state=active_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab)

     ELSE IF (sh_env%tresp_couplings_keyword) THEN 
        !frozen Hamiltonian case
        CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, state=active_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab)
     ELSE
         !> In a possible extension one should add here 
         CPABORT("Another way to calculate adiabatic forces should be developed")
     END IF

     !> This selects the most general DIABATIC_RK_FAST that does not include
     !diabatic nace calculation
     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
     CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
     if (prop_wf_keyword .eq. "DIABATIC_RK_FAST") then 
         ! in this case diabatic nace are never calculated, neither in the
         ! propagate_wf nor in the adiabatic nace calculation for the hopping
         ! probability (thus d_ij = Udot\dagger U)
         ! so that all the Sttp matrices are also not used and so the memory
         ! allocate for them can be deallocated 
         sh_env%no_diab_nace = .True.
     end if

     
     !> STORE THE COORDINATES IF NOT SUCCESSFUL RUN
     IF (sh_env%aom_couplings_keyword) THEN
        !------------------  AOM section --------------------------------------------------------
        sh_env%rm_tm      = positions_com_tp - velocities_com_tp*deltatMD
        sh_env%rm_t       = positions_com_tp
        sh_env%Stt = sh_env%Stptp
        !------------------------------------------------------------------------------------- 
        DEALLOCATE(PC_allrvecs)    
     END IF

     sh_env%sh_print_env = sh_print_env
     sh_env%ctin        = ctin
     sh_env%ctin_t        = ctin
     sh_env%H_t = sh_env%H_tp 
     sh_env%diabatic_nace_t = sh_env%diabatic_nace_tp
     sh_env%E_t = sh_env%E_tp
     sh_env%C_t = sh_env%C_tp
     sh_env%new_state = new_state
     sh_env%active_state = active_state
     !> IMPORTANT FOR THE FIRST TIME STEP
     sh_env%old_state = sh_env%active_state
     sh_env%new_state = -1
     sh_env%successful_hop = .FALSE.


    !> DEALLOCATE
    deallocate(ctin)
    deallocate(exact_forces)
    deallocate(site_energies, couplings, adiab_c, adiab_pop)
    DEALLOCATE(positions_com_tp)
    DEALLOCATE(velocities_com_tp)

  END SUBROUTINE propagate_namd_first_step

!**********************************************************************************
!> \brief 
!**********************************************************************************
  !Seems to change the nace not the nacv 
  SUBROUTINE nacv_first_ts(sh_env, diagonal_F, natom, mat, velocities_com_tp)

     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, intent(in)                            :: natom
     integer, intent(in)                            :: mat
     REAL(kind=dp), DIMENSION(:,:), INTENT(IN)      :: velocities_com_tp
     REAL(kind=dp), DIMENSION(:,:,:), POINTER, & 
                                    INTENT(IN)      :: diagonal_F

     CHARACTER(len=*), PARAMETER                    :: routineN ='nacv_first_ts', &
     routineP = moduleN//':'//routineN
     
     integer                                        :: i,j, iparticle, coord
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv

     sh_env%adiab_nace_tp = 0.0D0
     do i =1, mat
         do j = 1, mat
              if (j .ge. i) then
                  CALL calc_hop_nacv_aom(sh_env, diagonal_F, hop_nacv, i, j, fast=.TRUE.)
                  do iparticle=1,natom
                       do coord = 1,3
                            sh_env%adiab_nace_tp(i,j) = sh_env%adiab_nace_tp(i,j) + &
                                                    velocities_com_tp(coord, iparticle)* & 
                                                            hop_nacv(coord, iparticle)
                       end do
                  end do
              else
                  sh_env%adiab_nace_tp(i,j)= -sh_env%adiab_nace_tp(j,i)
              end if 
         end do
     end do
     sh_env%adiab_nace_tp = -sh_env%adiab_nace_tp 
  END SUBROUTINE nacv_first_ts

!**********************************************************************************
!> \brief  This construct the Hamiltonian, propagate the electronic equation of
!motion and defines if a hop is attempted or not.
!> authors AC and SG
!> edited 02/02/20
!**********************************************************************************
    SUBROUTINE propagate_namd_all_steps(nadiab, norbitals,                            &
                                  positions_tp,                                 &
                                  energies_tp,                                  &
                                  sh_env, itimes,                               &
                                  adiab_section, my_bo_keyword,                 & 
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword, diagonal_F)


     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, POINTER                               :: itimes
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     INTEGER                                        :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), POINTER  :: atomic_kind_set
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F

     CHARACTER(len=*), PARAMETER :: routineN ='propagate_namd_all_steps', &
     routineP = moduleN//':'//routineN

     TYPE(section_vals_type), POINTER               :: decoherence_section, optimization_section
     TYPE(section_vals_type), POINTER               :: e_conservation_section
     CHARACTER(LEN=default_string_length)           :: coupling_keyword
     integer                                        :: coord
     integer                                        :: mol,  mat, orb
     REAL(KIND=dp)                                  :: deltatMD, deltatE
     INTEGER, dimension(:), allocatable             :: class_forbidden_hop_tp
     LOGICAL                                        :: do_nace_speedup
     INTEGER                                        :: ELsteps   
     REAL(KIND=dp)                                  :: times 
     REAL(KIND=dp)                                  :: decohere_criterion
     INTEGER                                        :: imol,  natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE     :: exact_forces
     REAL(kind=dp), DIMENSION(:), ALLOCATABLE       :: site_energies, couplings
     LOGICAL                                        :: keyword_translate
     TYPE(sh_print_env_type)                        :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)           :: prop_wf_keyword 
     CHARACTER(LEN=default_string_length)           :: reversal_keyword
     real(kind=dp), DIMENSION(:,:), allocatable     :: PC_allrvecs
     REAL(KIND=dp), dimension(:), allocatable       :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable    :: adiab_c
     TYPE(cp_logger_type), POINTER                  :: logger
     INTEGER                                        :: unit_nr
!     INTEGER                                        :: i, k
     
     NULLIFY(logger)
     NULLIFY(optimization_section)
     NULLIFY(decoherence_section)
     NULLIFY(e_conservation_section)

     logger => cp_get_default_logger()

     CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
     !> Calculate the orbital overlap in the COM reference sistem of the
     !molecules
     CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS", l_val=keyword_translate)

     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
     decoherence_section => section_vals_get_subs_vals(adiab_section,"DECOHERENCE")
     e_conservation_section => section_vals_get_subs_vals(adiab_section,"ENERGY_CONSERVATION")
     
     !COMMENT SG: the following keyword could be removed since new nace are
     !better
     !> uses chunked nace (improved) subroutines as default 
     CALL section_vals_val_get(optimization_section, "DO_SPEEDUP_NACE", l_val=do_nace_speedup)

     CALL section_vals_val_get(adiab_section,"ELECTRONIC_PARTIAL_STEP",&
                               i_val=ELsteps)
     !> Instantaneous decoherence related threshold criteria 
     CALL section_vals_val_get(decoherence_section,"DECOHERENCE_CRITERION",&
             r_val=decohere_criterion)

     !atoms = sh_env%atoms_per_site !Number of atoms per site for AOM 
     natom = sh_env%natom !Total number of atoms
     deltatMD = sh_env%deltat !Nuclear timestep

     !> Number of sites
     mol = nadiab ! CHANGE_AC
     !> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
     !> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
     !> electronic timestep is a fraction of the molecular timestep
     deltatE = deltatMD/ELsteps !CHANGE_AC

     sh_print_env = sh_env%sh_print_env
     times = DBLE(itimes)
     sh_print_env%times = times
     
     ALLOCATE(positions_com_tp(3, natom))
     allocate(exact_forces(3,natom))
     !> OBTAIN COORDINATES, ENERGIES, STATES, etc. FROM PREVIOUS STEP
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     IF (sh_env%aom_couplings_keyword) THEN
         !------------------  AOM section -------------------------------------
         ALLOCATE(PC_allrvecs(5,mol*sh_env%atoms_per_site))
         sh_env%siteenergies = 0.0_dp
         do coord=1,3
            sh_env%crds(coord,:,1) = PACK(sh_env%rm_tm(coord,:), sh_env%ext_pres)
            sh_env%crds(coord,:,2) = PACK(sh_env%rm_t(coord,:), sh_env%ext_pres)
            sh_env%crds(coord,:,3) = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
            sh_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
         enddo

         !> CALCULATE THE CENTER OF MAS AND CONSTRUCT NEIGHBOURS LIST
         CALL calc_coms_and_nlist(mol, sh_env%cutoff_sites, orb, crds=sh_env%crds, &
                                    active_system_info=sh_env%active_system_info, NN=sh_env%NN)

         !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
         CALL calc_aom_Xmatrix_elements(sh_env, itimes, sh_env%atoms_per_site, mol, sh_env%dR, sh_env%do_random, sh_env%cbar, &
                          sh_env%scaling, PC_allrvecs, orb, sh_env%minao)

         !> CONVERT THE diabatic NACV WITH THE GOOD UNITS
         if (mod(itimes-1, sh_env%multiple_time_step) .eq. 0) then
             print *, "FOBSH| *** OFF-DIAGONAL GRADIENT CALCULATION"
             sh_env%nacv = sh_env%cbar*sh_env%nacv
         end if
         !TEST NACV ABSENCE:
         !sh_env%nacv = 0.0D0
         !print *, "NACV tp", sh_env%nacv(1,:,:)
         !-----------------------------------------------------------------------------
     ELSE IF (sh_env%frozen_couplings_keyword) THEN
         !> FROZEN COUPLINGS 
         CALL calc_frozen_H(sh_env, energies_tp, orb)

     ELSE IF (sh_env%tresp_couplings_keyword) THEN
         !> TRESP COUPLINGS 
         do coord=1,3
            sh_env%crds(coord,:,1) = 0.0D0 !time tm is not used 
            sh_env%crds(coord,:,2) = 0.0D0 !time t is not used
            sh_env%crds(coord,:,3) = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
         enddo

         !> CALCULATE THE CENTER OF MAS AND CONSTRUCT NEIGHBOURS LIST
         CALL calc_coms_and_nlist(mol, sh_env%cutoff_sites, orb, crds=sh_env%crds, &
                                    active_system_info=sh_env%active_system_info, NN=sh_env%NN)

         CALL calc_frozen_H(sh_env, energies_tp, orb)

     ELSE
         !> In a possible extension one should add here 
         CPABORT("Another coupling method to constrcut H should be chosen")
     END IF



!!TEST FROZEN 
!    do i =1, size(sh_env%H_tp,1)
!      do k = 1, size(sh_env%H_tp,1)
!         if (i .ne. k) then
!            if ((k+1 == i) .or.(k-1 == i)) then
!               sh_env%H_tp(i,k) = cmplx(0.0D0,0.001102479,dp) 
!            else
!               sh_env%H_tp(i,k) = cmplx(0.0D0,0.0D0,dp)
!            end if
!         end if
!     end do
!    end do
!sh_env%nacv = 0.0D0
!!print *, "sh_env%H_tp(:,:)", sh_env%H_tp

     !> THIS KEYWORD helps treating initial itimes correctly in case of restart file
     !usage 
     sh_env%restart_key = .False.

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian H_tp
     !It gives us back the eigenvalues E_tp and eigenvectors C_tp
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(sh_env%H_tp), sh_env%C_tp, sh_env%E_tp)

     !> store adiabatic nace and diabatic from the previous step
     sh_env%adiab_nace_t = sh_env%adiab_nace_tp

     !> Calculate the overlap between C_tp and C_t, reorder state if required,
     !> calculate the nace
     IF (do_nace_speedup) THEN
        CALL calc_adiab_nace(sh_env, adiab_section)
     ELSE
        CALL calc_adiab_nace_old(sh_env, adiab_section)
     END IF

     !>this if statement is usefull at t=1 becuase sh_env%diabatic_nace_tp
     !> is not present and we would need it in the interpolation scheme during the
     !> wavefunction propagation, so no we approximate that sh_env%diabatic_nace_t =0
     !> and sh_env%diabatic_nace_tp = sh_env%diabatic_nace_t
     if (itimes.eq.1 .or. sh_env%restart_key)  then
       sh_env%diabatic_nace_tp = sh_env%diabatic_nace_t
     end if
     
     IF (sh_env%aom_couplings_keyword) THEN
         !------------------  AOM section    ---------------------------------------------
         !> PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
         !>the CHARGE CARRIER WF is called sh_env%ctin
         sh_env%HS_propagation = .FALSE.
         CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
         SELECT CASE(prop_wf_keyword)                                 
         CASE("DIABATIC_RK")
               CALL propagate_wf(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                                 itimes, adiab_section, deltatMD, deltatE , times, mat)
         CASE("DIABATIC_RK_FAST")
               CALL propagate_wf_fast(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                                 itimes, adiab_section, deltatMD, deltatE , times, mat)
         CASE("DIABATIC_RK_DLI")
               CALL propagate_wf_diab_nace_interpol(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                                 itimes, adiab_section, deltatMD, deltatE , times, mat)
         CASE("ADIABATIC_RK")
               CALL propagate_adiab_wf(sh_env%E_t, sh_env%E_tp, sh_env%ctin, sh_env, nadiab, norbitals, ELsteps, &
                                       itimes, adiab_section, deltatMD, deltatE , times,  mat)
         CASE("HS_CORRECTED_RK")  
               !> ADIABATIC NACV EVALUATION FOR HS CORRECTION
               !> class_forbidden_hop_tp is calculated refferring to tp to be
               !consistent with the hopping probability, which is evaluated at tp.
               !These transition will be set to zero in the hamiltonian evolution
               !between t and tp. However the velocities used are the ones at t+dt/2,
               !this is a bit inconsistent but we could not find a better way. 

               ALLOCATE(velocities_com_tp(3, natom))
               CALL translate_com(sh_env%velocities_tp, keyword_translate, velocities_com_tp, &
                                              local_particles, atomic_kinds = atomic_kinds)  
               CALL classically_forbidden_hops(sh_env%active_state,diagonal_F,  sh_env, mat,  &
                      velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                              atomic_kind_set, class_forbidden_hop_tp)
               CALL prop_wf_HS_corr(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                               itimes, adiab_section, deltatMD, deltatE , times, &
                                mat, class_forbidden_hop_tp, sh_env%active_state)
               CALL section_vals_val_get(e_conservation_section,"METHOD_REVERSAL",c_val=reversal_keyword)
               if (reversal_keyword .eq. "ALWAYS") then
                    sh_env%HS_propagation = .FALSE.
               else 
                    sh_env%HS_propagation = .TRUE.
               end if
               DEALLOCATE(velocities_com_tp)
         CASE DEFAULT                                           
               CPABORT("There is no method to propagate the wavefunction.") 
         END SELECT                                             
         ! -------------------------------------------------------------------------------------------
     ELSE
         !> For now the only way to propagate the WF with an Hamiltoninian
         !different from aom is the FAST propagator, where diabatic NACEa are
         !neglected. If you wanted to still use the linear interpolation of
         !diabatic nace you should code up a different way to calculate diabatic_nace_t
         CALL propagate_wf_fast(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                                 itimes, adiab_section, deltatMD, deltatE , times, mat)
         
     END IF

     !> SURFACE HOPPING TEST
     !The eigenvectors are needed to calculate the surface hopping probabilities
     !The surfacehop subroutine calculates these probabilites and decides on the
     !new (possibly unchanged) surface for simulation
     sh_env%decohere = .FALSE.
     sh_env%attempt_hop = .FALSE.
     IF (my_test_keyword) THEN   
        sh_env%attempt_hop = .TRUE.
        sh_env%active_state = 1
        sh_env%new_state = mod(sh_env%active_state + 1, mat) 
        IF ((sh_env%new_state).EQ.0) sh_env%new_state = mat
     ELSE IF (.NOT.my_bo_keyword) THEN ! THEN IT IS A FSSH RUN
            call surf_hop_type(adiab_section, sh_env%C_t, sh_env%C_tp,   &
                                 deltatMD, &
                                 sh_env%new_state, &
                                 sh_env%attempt_hop, &
                                 decohere_criterion, &
                                 sh_print_env, sh_env)
     ENDIF

     print*, "FOBSH| ACTIVE STATE = ", sh_env%active_state

     IF (my_bo_keyword) THEN
          sh_print_env%first_adiabat = sh_env%active_state
          CALL adiab_print_sh_info(sh_print_env, keyword = 'BOMD')
     END IF 


     IF (sh_env%aom_couplings_keyword) THEN
         !------------------  AOM section    --------------------------------------------------------
         !CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
         !HOP ROUTINE IS. 
         call calcforces_aom(state=sh_env%active_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, nacv=sh_env%nacv, &
                               atoms=sh_env%atoms_per_site, & 
                               scaling=sh_env%scaling, cbar=sh_env%cbar, Fin=diagonal_F, forces_adiab=sh_env%forces_adiab, &
                               ext_pres=sh_env%ext_pres)
         IF (sh_env%attempt_hop) THEN
             call calcforces_aom(state=sh_env%new_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, nacv=sh_env%nacv, &
                                  atoms=sh_env%atoms_per_site, & 
                                  scaling=sh_env%scaling, cbar=sh_env%cbar, Fin=diagonal_F, &
                                                          forces_adiab=sh_env%forces_adiab_new, &
                                                          ext_pres=sh_env%ext_pres)
         ENDIF
         ! --------------------------------------------------------------------------------------------------
     ELSE IF (sh_env%frozen_couplings_keyword) THEN
        CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, state=sh_env%active_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab)
        IF (sh_env%attempt_hop) THEN
                     CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp),C=sh_env%C_tp, state=sh_env%new_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab_new)
        END IF

     ELSE IF (sh_env%tresp_couplings_keyword) THEN
        CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, state=sh_env%active_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab)
        IF (sh_env%attempt_hop) THEN
                     CALL calcforces_frozenH(Hin=-AIMAG(sh_env%H_tp),C=sh_env%C_tp, state=sh_env%new_state, Fin=diagonal_F, &
                                  forces_adiab=sh_env%forces_adiab_new)
        END IF

     ELSE
         !> In a possible extension one should add here 
         CPABORT("Another way to calculate adiabatic forces should be developed")
     END IF
     

     !> PRINT
     CALL adiab_print_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)
     CALL adiab_print_pseudo_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)

     !> print site energies 
     allocate(site_energies(mat))
     allocate(couplings(mat-1))
     allocate(adiab_c(mat))
     allocate(adiab_pop(mat))
 
     do imol=1,mat
         site_energies(imol) = -AIMAG(sh_env%H_tp(imol,imol))
         if (imol .ne. mat) then
           couplings(imol) = -AIMAG(sh_env%H_tp(imol,imol+1))
         end if 
     end do


     unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIAB_POPULATIONS",&
               extension=".xyz", middle_name="adiab_pop",log_filename=.FALSE.)
     if (unit_nr > 0) then
        !CALL rotation_diab_to_adiab(sh_env%C_tp, sh_env%ctin, mat, adiab_ctin, adiab_pop)
        CALL diab_to_adiab( sh_env%C_tp, sh_env%ctin, adiab_c, adiab_pop)
        CALL adiab_print_adiab_pop(adiab_section, adiab_pop, nadiab,norbitals,itimes, deltatMD, times) 
     end if

     CALL adiab_print_site_energies(adiab_section, site_energies, nadiab, deltatMD, times)
     CALL adiab_print_off_diagonals(adiab_section, couplings, nadiab-1, deltatMD, times)
     CALL adiab_print_adiabat_energies(adiab_section, sh_env%E_tp, nadiab, deltatMD, times)
     IF (sh_env%aom_couplings_keyword) THEN
        CALL adiab_print_nacv(adiab_section, sh_env%nacv, nadiab, norbitals, sh_env%atoms_per_site ,itimes, deltatMD)
        CALL adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, sh_env%atoms_per_site, nadiab)
     END IF
     IF (sh_env%analytics) THEN
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF

     !> TRANSLATE THE COORDINATE FOR ONE TIMESTEP
     IF (sh_env%aom_couplings_keyword) THEN
         !----------------------- AOM ----------------------
         sh_env%rm_tm      = sh_env%rm_t
         sh_env%rm_t       = positions_com_tp
         ! --------------------------------------------------
         DEALLOCATE(PC_allrvecs)
     END IF

     sh_env%sh_print_env = sh_print_env

     !> DEALLOCATE
     DEALLOCATE(positions_com_tp) 
     deallocate(exact_forces)
     deallocate(site_energies, couplings, adiab_c, adiab_pop)
    END SUBROUTINE propagate_namd_all_steps

END MODULE sh_main
