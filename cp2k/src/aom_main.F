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
! (to do surface hopping simulations) or afterwards (classical path approximation).
!
! This program takes molecular coordinates as inputs.
! From these coordinate inputs, it calculates the electronic Hamiltonian
! in a molecular orbital basis, using the overlapFINAL module to calculate
! the off-diagonals and either accepting the diagonals as inputs or calculating
! them in some other way.
!
! With the electronic Hamiltonian, the program then propagates the electron using the RK4 method.
! At the final electronic state, the probabilities for hopping between adiabatic surfaces are calculated.
! The program then decides what adiabatic surface the system is on and calculates
! the appropriate atomic forces.
!
! In the classical path approximation, the above steps are repeated until the
! end of the input molecular trajectory.  For surface-hopping, the data here is fed back into
! the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE aom_main 

     USE kinds,                        ONLY: default_string_length,&
                                             dp


!    USE aom_overlapfinal,             ONLY:  connect_list2, &
!                                             calc_sab, &
!                                             sab_screen_new, &
!                                             calc_dRSab, &
!                                             calc_rvecs, &
!                                             calc_sab_pete

     USE aom_utils,                    ONLY: makeX,    &
                                             calc_aom_Xmatrix_elements, &
                                             make_adiabX,      &
                                             calcforces,     &
                                             read_listcoeff, &
                                             calc_adiab_nace, &
                                             calc_adiab_nace_old, &
                                             calc_coms_and_nlist_aom

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

     USE non_adiabatic_utils,          ONLY: decohere_wf,       &
                                             diab_to_adiab,     & 
                                             surf_hop_type,     &
                                             choose_firstadiabat, &
                                             print_matrix_standard_out,& 
                                             translate_com, &
                                             adiabaticity,   &
                                             propagate_wf, &
                                             propagate_wf_fast, &
                                             propagate_adiab_wf, &
                                             propagate_wf_diab_nace_interpol, &
                                             prop_wf_HS_corr, & 
                                             classically_forbidden_hops, &
                                             define_active_region, &
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
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_main'   

        PUBLIC :: construct_ham_X_SH

CONTAINS

    SUBROUTINE construct_ham_X_SH(nadiab, aom_atoms, norbitals,  &
                                  positions_tp,  &
                                  energies_tp,  &
                                  sh_env, itimes, &
                                   coupling_section, &
                                  adiab_section, my_bo_keyword,                 & 
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword, diagonal_F)


     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     INTEGER, INTENT(IN)                            :: aom_atoms
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, POINTER                               :: itimes
     TYPE(section_vals_type), POINTER               :: coupling_section
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     INTEGER                                        :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), POINTER  :: atomic_kind_set
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F

     CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH', &
     routineP = moduleN//':'//routineN

     integer                                        :: coord
     integer                                        :: mol, atoms, c,  mat, orb
     REAL(KIND=dp)                                  :: cutoff_sites, cutoff_connect, &
                                                       scaling, deltatMD, deltatE
     REAL(KIND=dp)                                  :: cbar, minao
     INTEGER, dimension(:), allocatable             :: class_forbidden_hop_tp
     LOGICAL                                        :: analytics, do_nace_speedup
     INTEGER                                        :: ELsteps   
     REAL(KIND=dp)                                  :: dR, times 
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
     LOGICAL                                        :: do_random 
     REAL(KIND=dp), dimension(:), allocatable       :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable    :: adiab_c
     TYPE(cp_logger_type), POINTER                  :: logger
     INTEGER                                        :: unit_nr
!TEST
!     INTEGER                                    :: d, i,j


     NULLIFY(logger)
     logger => cp_get_default_logger()

     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=cutoff_sites)
     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     !> minao is an atomic overlap cutoff
     CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=minao)
     !> dR is the increment for the gradient finite difference calculation 
     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
             r_val=dR)
     !> analytics would probably calculates analytic forces for a 2 state sytems (Spencer 2016)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=analytics)
     !> This remove bias coming from the finite difference calculation of the
     !> gradients of the Hamiltonian from AOM
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)
     !> Multiple time step calculation of aom gradient
     CALL section_vals_val_get(coupling_section,"MULTIPLE_TIME_STEP",&
                               i_val=sh_env%multiple_time_step)

     
     CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS", l_val=keyword_translate)
     !COMMENT SG: the following keyword could be removed since new nace are
     !better
     !> uses chunked nace (improved) subroutines as default 
     CALL section_vals_val_get(adiab_section, "DO_SPEEDUP_NACE", l_val=do_nace_speedup)
     CALL section_vals_val_get(adiab_section,"ELECTRONIC_PARTIAL_STEP",&
                               i_val=ELsteps)
     !> Instantaneous decoherence related threshold criteria 
     CALL section_vals_val_get(adiab_section,"DECOHERENCE_CRITERION",&
             r_val=decohere_criterion)

     !> Scaling factor Hab = scaling*Sab_bar
     scaling = sh_env%scaling
     cbar = sh_env%cbar
     atoms = sh_env%atoms !Number of atoms per site
     natom = sh_env%natom !Total number of atoms
     deltatMD = sh_env%deltat !Nuclear timestep

     !> Number of sites
     mol = nadiab ! CHANGE_AC
     !> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
     !> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
     !> Numbers of active atoms 
     c = mol*atoms
     !> electronic timestep is a fraction of the molecular timestep
     deltatE = deltatMD/ELsteps !CHANGE_AC

     sh_print_env = sh_env%sh_print_env
     times = DBLE(itimes)
     sh_print_env%times = times
     
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))

     !> OBTAIN COORDINATES, ENERGIES, STATES, etc. FROM PREVIOUS STEP
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     CALL translate_com(sh_env%velocities_tp, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     sh_env%siteenergies = 0.0_dp
     do coord=1,3
        sh_env%crds(coord,:,1) = PACK(sh_env%rm_tm(coord,:), sh_env%ext_pres)
        sh_env%crds(coord,:,2) = PACK(sh_env%rm_t(coord,:), sh_env%ext_pres)
        sh_env%crds(coord,:,3) = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
        sh_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo

     !> store adiabatic nace and diabatic from the previous step
     sh_env%adiab_nace_t = sh_env%adiab_nace_tp

     !> CALCULATE THE CENTER OF MAS AND CONSTRUCT NEIGHBOURS LIST
     CALL calc_coms_and_nlist_aom(sh_env, atoms, mol, cutoff_sites, orb)

     !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
     CALL calc_aom_Xmatrix_elements(sh_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb, minao)

     !> CONVERT THE diabatic NACV WITH THE GOOD UNITS
    if (mod(itimes-1, sh_env%multiple_time_step) .eq. 0) then
        print *, "FOBSH| *** OFF-DIAGONAL GRADIENT CALCULATION"
        sh_env%nacv = cbar*sh_env%nacv
    end if
    !TEST NACV ABSENCE:
    !sh_env%nacv = 0.0D0
    !print *, "NACV tp", sh_env%nacv(1,:,:)


     !> THIS KEYWORD helps treating initial itimes correctly in case of restart file
     !usage 
     sh_env%restart_key = .False.

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian H_tp
     !It gives us back the eigenvalues E_tp and eigenvectors C_tp
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(sh_env%H_tp), sh_env%C_tp, sh_env%E_tp)
     !Calculate the overlap between C_tp and C_t, reorder state if required,
     !calculate the nace
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
     
     !> PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
     !>the CHARGE CARRIER WF is called sh_env%ctin
     sh_env%HS_propagation = .FALSE.
     CALL section_vals_val_get(adiab_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
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
           !between t and tp.
           CALL classically_forbidden_hops(sh_env%active_state,diagonal_F,  sh_env, mat,  &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, class_forbidden_hop_tp)
           CALL prop_wf_HS_corr(sh_env%ctin, sh_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times, &
                            mat, class_forbidden_hop_tp, sh_env%active_state)
           CALL section_vals_val_get(adiab_section,"METHOD_REVERSAL",c_val=reversal_keyword)
           if (reversal_keyword .eq. "ALWAYS") then
                sh_env%HS_propagation = .FALSE.
           else 
                sh_env%HS_propagation = .TRUE.
           end if
     CASE DEFAULT                                           
           CPABORT("There is no method to propagate the wavefunction.") 
     END SELECT                                             


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


     !CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
     !HOP ROUTINE IS. 
     call calcforces(state=sh_env%active_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, nacv=sh_env%nacv, &
                           atoms=atoms, & 
                           scaling=scaling, cbar=cbar, Fin=diagonal_F, forces_adiab=sh_env%forces_adiab, &
                           ext_pres=sh_env%ext_pres)
     IF (sh_env%attempt_hop) THEN
         call calcforces(state=sh_env%new_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, nacv=sh_env%nacv, &
                              atoms=atoms, & 
                              scaling=scaling, cbar=cbar, Fin=diagonal_F, forces_adiab=sh_env%forces_adiab_new, &
                              ext_pres=sh_env%ext_pres)
     ENDIF

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
     CALL adiab_print_nacv(adiab_section, sh_env%nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
     CALL adiab_print_adiabat_energies(adiab_section, sh_env%E_tp, nadiab, deltatMD, times)
     CALL adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, atoms, nadiab)
     IF (analytics) THEN
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF

     !> TRANSLATE THE COORDINATE FOR ONE TIMESTEP
     sh_env%rm_tm      = sh_env%rm_t
     sh_env%rm_t       = positions_com_tp
     sh_env%sh_print_env = sh_print_env

     !> DEALLOCATE
     DEALLOCATE(positions_com_tp) 
     DEALLOCATE(velocities_com_tp)
     deallocate(exact_forces)
     deallocate(site_energies, couplings, adiab_c, adiab_pop)
     DEALLOCATE(PC_allrvecs)
    END SUBROUTINE construct_ham_X_SH



END MODULE aom_main
