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

MODULE sh_types
  USE kinds,                           ONLY: default_string_length, &
                                             dp
  USE input_section_types,             ONLY: section_vals_type

  USE parallel_rng_types,              ONLY: rng_stream_type


IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sh_types'


   TYPE real_ES_NN
      INTEGER                                    :: atoma, atomb
   END  TYPE real_ES_NN 

   TYPE sh_print_env_type
!> General variable
     INTEGER                                     :: nadiab, norbitals
     REAL(KIND = dp)                             :: times, deltatMD
     TYPE(section_vals_type), POINTER            :: adiab_section
     REAL(KIND = dp), DIMENSION(:), POINTER      :: adiab_energies_t   
!> Variable for the first timestep
     REAL(KIND = dp), DIMENSION(:), POINTER      :: pop_init
     INTEGER                                     :: first_adiabat
     REAL(KIND = dp)                             :: initrandom
!> Variable for Tully's test
     LOGICAL                                     :: attempt_hop
     REAL(KIND = dp)                             :: sumprob
     REAL(KIND = dp)                             :: enermin, exprob, sumflux 
     LOGICAL                                     :: successful_hop
     LOGICAL                                     :: additional_prints
     INTEGER                                     :: old_state
     INTEGER                                     :: new_state
     LOGICAL                                     :: renorm
     REAL(KIND = dp), DIMENSION(:), POINTER      :: proba
     REAL(KIND = dp), DIMENSION(:), POINTER      :: tau_deco
     REAL(KIND = dp)                             :: hop
!Variable for the rescaling
     REAL(KIND = dp)                             :: delta_energy
     REAL(KIND = dp)                             :: kinetic_energy
     REAL(KIND = dp)                             :: A, B
     REAL(KIND = dp)                             :: rescale_factor
     REAL(KIND = dp), DIMENSION(3)               :: momentum_shift
!> Variable for the reversing
     LOGICAL                                     :: reverse
     CHARACTER(LEN=default_string_length)        :: reversal_keyword
     REAL(KIND = dp)                             :: product_pd, product_dfn, product_dfo
!> Variable for the decoherence
     LOGICAL                                     :: decohere
     LOGICAL                                     :: collapse
!     LOGICAL, DIMENSION(:), ALLOCATABLE          :: pure_deph
!> To store energies in the Jacob's Method
!> The probability array for the hop
!> Random number for test hop
!> If renormalization occurs
!> If renormalization occurs: the sum of probabilities
!> Old and new state number
!> Delta energy between old_state and putative new_state
!> Random number for initial choice
!> Initial population in the adiabatic states
!> First adiabat after the random choice
   END TYPE sh_print_env_type

!***************************************************************************************
!> \brief A storage container for the arrays used in calculating the adiab NACE.
!> \authour Matt
!> \note I have kept this seperate to make it easy to see what changes I've
!        made. Feel free to hard-code it into the sh_env_type.
!***************************************************************************************
TYPE adiab_nace_env_type
   REAL(KIND=dp), dimension(:,:), allocatable :: T1
   REAL(KIND=dp), dimension(:,:), allocatable :: T2
   REAL(KIND=dp), dimension(:,:), allocatable :: DT
   REAL(KIND=dp), dimension(:,:), allocatable :: diab_D
   REAL(KIND=dp), dimension(:,:), allocatable :: DS

   INTEGER, dimension(:),  allocatable        :: map_new2old
   INTEGER, dimension(:),  allocatable        :: map_old2new
   INTEGER, dimension(:),  allocatable        :: basic
   INTEGER, dimension(:),  allocatable        :: zeros
   INTEGER, dimension(:),  allocatable        :: not_associated

   !SG NEW VARIABLES
   REAL(KIND=dp), dimension(:,:), allocatable :: Ssqrt

END TYPE adiab_nace_env_type


!***************************************************************************************
!> \brief A storage container for the arrays used in fssh
!> \authors AC and SG
!> \edited 02/02/20
!***************************************************************************************
   TYPE sh_env_type
!> AOM arrays specific for donor acceptor MULTIPHASE
     REAL(KIND=dp), DIMENSION(:), POINTER     :: pcoeffs_d, pcoeffs_a
     REAL(KIND=dp), DIMENSION(:), POINTER     :: scoeffs_d, scoeffs_a
     INTEGER, DIMENSION(:), POINTER           :: atomlist_d, atomlist_a !List atoms with orbitals
     integer, dimension(:,:), POINTER         :: connlist1_d, connlist2_d, & ! AOM connectivities of individual mols
                                                 connlist1_a, connlist2_a
     integer, dimension(:), ALLOCATABLE       :: ext_elements_d, ext_elements_a ! atomic number
     !integer, dimension(:,:), POINTER         :: bneighbour_d, aneighbour_d
     !integer, dimension(:,:), POINTER         :: bneighbour_a, aneighbour_a
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: nacv_dphase, nacv_aphase   ! off-diagonal gradient of aom Hamilton.
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: nacv_fe, nacv_ct
     REAL(KIND=dp), DIMENSION(:), POINTER     :: pcoeffs_fe, pcoeffs_ct
     REAL(KIND=dp), DIMENSION(:), POINTER     :: scoeffs_fe, scoeffs_ct
     INTEGER, DIMENSION(:), POINTER           :: atomlist_fe, atomlist_ct !should be able to get rid of
     integer, dimension(:,:), POINTER         :: connlist1_fe, connlist2_fe, & !AOM connectivities of individual mols
                                                 connlist1_ct, connlist2_ct
     integer, dimension(:), ALLOCATABLE       :: ext_elements_fe, ext_elements_ct ! atomic number
     REAL(KIND = dp)                          :: scaling_donor, scaling_acceptor, scaling_fect
     LOGICAL                                  :: fe_donor  !True: donor is FE, False: acceptor is FE
     REAL(KIND = dp)                          :: e_offset  !WTP_4.2: E_FE(R_gs) - E_CT(R_gs)

     !SINGLE PHASE XSH AOM PARAMETERS DEFINED HERE
     REAL(KIND = dp)                          :: scaling_hh, scaling_ee, scaling_fect_e, scaling_fect_h
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: nacv_hh, nacv_ee, nacv_fe_h, nacv_ct_h, nacv_fe_e, nacv_ct_e
     REAL(KIND=dp), DIMENSION(:), POINTER       :: pcoeffs_ee, scoeffs_ee, pcoeffs_hh, scoeffs_hh
     REAL(KIND=dp), DIMENSION(:), POINTER       :: pcoeffs_fe_h, scoeffs_fe_h, pcoeffs_ct_h, scoeffs_ct_h
     REAL(KIND=dp), DIMENSION(:), POINTER       :: pcoeffs_fe_e, scoeffs_fe_e, pcoeffs_ct_e, scoeffs_ct_e
     integer, dimension(:), ALLOCATABLE       :: elements_sp

     LOGICAL                                  :: freeze_multip_forces

!> AOM arrays
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: rm_t,&        ! position step min 1
                                                 rm_tm       ! position step min 2
     !> Information for the overlap method: calculate at the first timestep
     INTEGER, DIMENSION(:), POINTER           :: elements, atomlist !List atoms with orbitals
     INTEGER, DIMENSION(:,:,:,:), POINTER     :: neighbourlist !very big list allocated 
                                                               !only if MTS=1
     integer, dimension(:,:), POINTER         :: connlist1, connlist2, & ! AOM neighbouring
                                                 bneighbour, aneighbour
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: pcoeffs ! aom coeffs for single phase H
     REAL(KIND=dp), DIMENSION(:), POINTER     :: scoeffs ! aom coeffs for single phase H
     REAL(KIND = dp)                          :: scaling, cbar, minao !For AOM calculations
     !> Calculate at each timestep
     REAL(KIND=dp), dimension(:,:,:), POINTER :: pcoeffsrun ! normalized AOM coefficients
     REAL(KIND=dp), dimension(:,:), POINTER   :: At, Atm, Atp, Bt, Btm, &
                                                 !A and B are donor and acceptor coordinates at different times  
                                                 Btp, siteenergies
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: nacv! off-diagonal gradient of aom Hamilton.
     !REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER :: nacv_tot! off-diagonal gradient of aom Hamilton.
     REAL(KIND=dp), dimension(:,:), POINTER   :: Stt, Stmt, Sttm, Stmtm, &
                                                 Stptp, Sttp, Stpt
     !> multiple time step
     INTEGER                                  :: multiple_time_step
     REAL(KIND=dp)                            :: cutoff_sites, cutoff_connect
     LOGICAL                                  :: analytics, do_random
     REAL(KIND=dp)                            :: dR

!> TRESP arrays
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tresp_charges
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: grad_tresp! off-diagonal gradient of tresp Hamilton.
     !> this is a test keyword for full storage of tresp gradient
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER :: test_grad_tresp! off-diagonal gradient of tresp Hamilonian
     LOGICAL                                  ::  Gaussian_charges

!> BORDERLINE used for AOM and some other couplings mainly TRESP
     LOGICAL, DIMENSION(:), POINTER           :: ext_pres !T or F array 
                                                          !indicating active atoms
     REAL(KIND=dp), dimension(:,:,:), POINTER :: crds !crds active atoms in COM reference frames
     INTEGER, DIMENSION(:,:), POINTER         :: NN   !intermolacular neighbrouring list
     integer                                  :: atoms_per_site !Number of atoms per site
     REAL(KIND=dp)                            :: cutoff_sites_fe

!> FROZEN_COUPLINGS section 
     REAL(KIND = dp)                          :: first_off_d_Hab
     REAL(KIND=dp), DIMENSION(:), POINTER     :: frz_couplings
     INTEGER, DIMENSION(:,:), POINTER         :: frz_coup_connectivity
     LOGICAL                                  :: connectivity_frz_c
     

!> Couplings section possibilities
     LOGICAL                                  :: aom_couplings_keyword, frozen_couplings_keyword 
     LOGICAL                                  :: tresp_couplings_keyword 
     LOGICAL                                  :: multiphase_couplings_keyword 
     LOGICAL                                  :: multifect_couplings_keyword 
     LOGICAL                                  :: sp_fect_keyword !single-phase fect keyword
 
     INTEGER, DIMENSION(:,:), ALLOCATABLE     :: full_system_info  
     INTEGER, DIMENSION(:,:), ALLOCATABLE     :: active_system_info  

!> Site eneries section
!> Multi state approach in case of a different number of states and active
!molecules
     LOGICAL                                  :: multistate
     LOGICAL                                  :: add_coulomb_pot
     INTEGER                                  :: number_ct_states
     REAL(KIND=DP), DIMENSION(:), POINTER     :: coulomb_barrier
     REAL(KIND=DP)                            :: energy_offset

!> Constant over the simulations
     integer                                  :: natom !Total number of atoms
     REAL(KIND=dp)                            :: deltat !Timestep

     TYPE(rng_stream_type), POINTER           :: rng_stream !the random seed generator

     REAL(KIND = dp), DIMENSION(:,:), POINTER :: forces_adiab, &     
                                                 forces_adiab_new
!> index active molecules
     INTEGER, DIMENSION(:), ALLOCATABLE       :: index_active_mols

!> Use to store the present info in the future timestep
     REAL(KIND = dp), DIMENSION(:,:), POINTER     :: velocities_tp ! velocities step
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE       :: C_t, C_tp !Hamiltoninan eigenvectors
     REAL(KIND=dp), DIMENSION(:,:), POINTER       :: diabatic_nace_tp, & 
                                                     diabatic_nace_t  !othogonal diab NACE
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: E_t, E_tp
     COMPLEX(KIND=dp), DIMENSION(:,:), POINTER    :: H_t, H_tp
     COMPLEX(KIND=dp), DIMENSION(:),POINTER       :: ctin, ctin_t !wavefunction at tp and t
     COMPLEX(KIND=dp), dimension(:), ALLOCATABLE  :: readcoeff ! coeff from the restart file 
     INTEGER                                      :: active_state_res
     LOGICAL                                      :: restart_key

!> Decoherence part
     LOGICAL                                 :: do_deco !should we do decoherence in that run?
     !> arrays used in the decoherence spurious transfer correction
     COMPLEX(KIND=dp), dimension(:), POINTER :: ctin_before_deco 
     COMPLEX(KIND=dp), dimension(:), POINTER :: ctin_after_deco
     COMPLEX(KIND=dp), dimension(:), POINTER :: ctin_after_tranf_correction
     INTEGER, dimension(:), POINTER          :: active_region
     !> Logicals for instant collapse or IDA
     LOGICAL                                 :: decohere !activate the decoherence according to instant collapse criteria
     LOGICAL                                 :: futur_decohere !futur decoherence possible in instant collapse

!> States index and reordering
     INTEGER                                 :: new_state, old_state, active_state
     REAL(KIND=dp), dimension(:,:), POINTER  :: adiab_overlap 
     LOGICAL                                 :: attempt_hop
     LOGICAL                                 :: successful_hop
     LOGICAL                                 :: do_reorder
     LOGICAL                                 :: no_diab_nace

!> Required for HS propagation
     LOGICAL                                 :: HS_propagation
     REAL(KIND=dp), dimension(:,:), POINTER  :: adiab_nace_tp, adiab_nace_t
     INTEGER, dimension(:), POINTER          :: class_forbidden_hop_t

!> To print
     TYPE(sh_print_env_type)                 :: sh_print_env

! Full Coulomb Interactions
     LOGICAL                                 :: finished_coulomb_step

!> CHANGE ME, added nace_allocatables
    TYPE(adiab_nace_env_type)                :: adiab_nace_env

   END TYPE sh_env_type

   PUBLIC :: sh_env_type,       &
             sh_print_env_type, &
             sh_release,        &
             sh_init


CONTAINS

   !**********************************************************
   !> \brief Will initialise the attributes of the sh_env_type
   !>        . This basically just means NULLIFYing pointers.
   !> \author Matt
   !**********************************************************
   SUBROUTINE sh_init(sh_env)
      TYPE(sh_env_type)                      :: sh_env

      CHARACTER(len=*), PARAMETER :: routineN ='sh_init',  &
                                     routineP = moduleN//':'//routineN
      NULLIFY(sh_env%rm_t, sh_env%elements, sh_env%siteenergies,  &
              sh_env%atomlist, sh_env%neighbourlist,              &
              sh_env%connlist1, sh_env%connlist2, sh_env%Btp,     &
              sh_env%pcoeffs, sh_env%scoeffs, sh_env%pcoeffsrun,  &
              sh_env%At, sh_env%Atm, sh_env%Atp, sh_env%Bt,       &
              sh_env%Btm, sh_env%nacv, sh_env%Stt, sh_env%Stmt,   &
              sh_env%Sttm, sh_env%Stmtm, sh_env%grad_tresp,       &
              sh_env%test_grad_tresp,                             &
              sh_env%crds, sh_env%NN, sh_env%frz_couplings,       &
              sh_env%frz_coup_connectivity, sh_env%rng_stream,    &
              sh_env%forces_adiab, sh_env%velocities_tp,          &
              sh_env%diabatic_nace_tp, sh_env%H_t, sh_env%H_tp,   &
              sh_env%ctin, sh_env%ctin_t, sh_env%bneighbour,      &
              sh_env%ctin_before_deco, sh_env%ctin_after_deco,    &
              sh_env%ctin_after_tranf_correction, sh_env%rm_tm,   &
              sh_env%active_region, sh_env%adiab_overlap,         &
              sh_env%adiab_nace_tp, sh_env%adiab_nace_t,          &
              sh_env%class_forbidden_hop_t, sh_env%aneighbour,    &
              sh_env%forces_adiab_new, sh_env%diabatic_nace_t,    &
              sh_env%Stptp, sh_env%Stpt, sh_env%Sttp              )

      
   END SUBROUTINE

   !**********************************************************
   !> \brief Deallocate the sh environment
   !> \params sh_env => The surface hopping environment
   !**********************************************************
   SUBROUTINE sh_release(sh_env)
      TYPE(sh_env_type)                      :: sh_env

      CHARACTER(len=*), PARAMETER :: routineN ='sh_release',  &
                                     routineP = moduleN//':'//routineN
      
      ! First handle the allocatables
      IF (ALLOCATED(sh_env%tresp_charges)) THEN
          DEALLOCATE(sh_env%tresp_charges)
      END IF

      IF (ALLOCATED(sh_env%full_system_info)) THEN
          DEALLOCATE(sh_env%full_system_info)
      END IF

      IF (ALLOCATED(sh_env%active_system_info)) THEN
          DEALLOCATE(sh_env%active_system_info)
      END IF

      IF (ALLOCATED(sh_env%C_t)) THEN
          DEALLOCATE(sh_env%C_t)
      END IF

      IF (ALLOCATED(sh_env%C_tp)) THEN
          DEALLOCATE(sh_env%C_tp)
      END IF

      IF (ALLOCATED(sh_env%E_t)) THEN
          DEALLOCATE(sh_env%E_t)
      END IF

      IF (ALLOCATED(sh_env%E_tp)) THEN
          DEALLOCATE(sh_env%E_tp)
      END IF

      IF (ALLOCATED(sh_env%readcoeff)) THEN
          DEALLOCATE(sh_env%readcoeff)
      END IF

      

      ! Now handle the pointers
      ! N.B. These pointers should mostly be allocatables.
      
      IF (ASSOCIATED(sh_env%forces_adiab_new)) THEN
          DEALLOCATE(sh_env%forces_adiab_new)
      END IF

      IF (ASSOCIATED(sh_env%siteenergies)) THEN
          DEALLOCATE(sh_env%siteenergies)
      END IF

      IF (ASSOCIATED(sh_env%aneighbour)) THEN
          DEALLOCATE(sh_env%aneighbour)
      END IF

      IF (ASSOCIATED(sh_env%bneighbour)) THEN
          DEALLOCATE(sh_env%bneighbour)
      END IF

      IF (ASSOCIATED(sh_env%rm_t)) THEN
          DEALLOCATE(sh_env%rm_t)
      END IF
      
      IF (ASSOCIATED(sh_env%rm_tm)) THEN
          DEALLOCATE(sh_env%rm_tm)
      END IF
      
      IF (ASSOCIATED(sh_env%elements)) THEN
          DEALLOCATE(sh_env%elements)
      END IF
      
      IF (ASSOCIATED(sh_env%atomlist)) THEN
          DEALLOCATE(sh_env%atomlist)
      END IF
      
      IF (ASSOCIATED(sh_env%neighbourlist)) THEN
          DEALLOCATE(sh_env%neighbourlist)
      END IF
      
      IF (ASSOCIATED(sh_env%connlist1)) THEN
          DEALLOCATE(sh_env%connlist1)
      END IF
      
      IF (ASSOCIATED(sh_env%connlist2)) THEN
          DEALLOCATE(sh_env%connlist2)
      END IF
      
      IF (ASSOCIATED(sh_env%pcoeffs)) THEN
          DEALLOCATE(sh_env%pcoeffs)
      END IF
      
      IF (ASSOCIATED(sh_env%scoeffs)) THEN
          DEALLOCATE(sh_env%scoeffs)
      END IF
      
      IF (ASSOCIATED(sh_env%pcoeffsrun)) THEN
          DEALLOCATE(sh_env%pcoeffsrun)
      END IF
      
      IF (ASSOCIATED(sh_env%At)) THEN
          DEALLOCATE(sh_env%At)
      END IF
      
      IF (ASSOCIATED(sh_env%Atm)) THEN
          DEALLOCATE(sh_env%Atm)
      END IF
      
      IF (ASSOCIATED(sh_env%Atp)) THEN
          DEALLOCATE(sh_env%Atp)
      END IF
      
      IF (ASSOCIATED(sh_env%Bt)) THEN
          DEALLOCATE(sh_env%Bt)
      END IF
      
      IF (ASSOCIATED(sh_env%Btp)) THEN
          DEALLOCATE(sh_env%Btp)
      END IF
      
      IF (ASSOCIATED(sh_env%Btm)) THEN
          DEALLOCATE(sh_env%Btm)
      END IF
      
      IF (ASSOCIATED(sh_env%nacv)) THEN
          DEALLOCATE(sh_env%nacv)
      END IF
      
      IF (ASSOCIATED(sh_env%Stt)) THEN
          DEALLOCATE(sh_env%Stt)
      END IF
      
      IF (ASSOCIATED(sh_env%Stpt)) THEN
          DEALLOCATE(sh_env%Stpt)
      END IF
      
      IF (ASSOCIATED(sh_env%Sttp)) THEN
          DEALLOCATE(sh_env%Sttp)
      END IF
      
      IF (ASSOCIATED(sh_env%Stptp)) THEN
          DEALLOCATE(sh_env%Stptp)
      END IF
      
      IF (ASSOCIATED(sh_env%Stmt)) THEN
          DEALLOCATE(sh_env%Stmt)
      END IF
      
      IF (ASSOCIATED(sh_env%Sttm)) THEN
          DEALLOCATE(sh_env%Sttm)
      END IF
      
      IF (ASSOCIATED(sh_env%Stmtm)) THEN
          DEALLOCATE(sh_env%Stmtm)
      END IF
      
      IF (ASSOCIATED(sh_env%grad_tresp)) THEN
          DEALLOCATE(sh_env%grad_tresp)
      END IF
      
      IF (ASSOCIATED(sh_env%test_grad_tresp)) THEN
          DEALLOCATE(sh_env%test_grad_tresp)
      END IF
      
      IF (ASSOCIATED(sh_env%crds)) THEN
          DEALLOCATE(sh_env%crds)
      END IF
      
      IF (ASSOCIATED(sh_env%NN)) THEN
          DEALLOCATE(sh_env%NN)
      END IF
      
      IF (ASSOCIATED(sh_env%frz_couplings)) THEN
          DEALLOCATE(sh_env%frz_couplings)
      END IF
      
      IF (ASSOCIATED(sh_env%frz_coup_connectivity)) THEN
          DEALLOCATE(sh_env%frz_coup_connectivity)
      END IF
      
      IF (ASSOCIATED(sh_env%rng_stream)) THEN
          DEALLOCATE(sh_env%rng_stream)
      END IF
      
      IF (ASSOCIATED(sh_env%forces_adiab)) THEN
          DEALLOCATE(sh_env%forces_adiab)
      END IF
      
      IF (ASSOCIATED(sh_env%velocities_tp)) THEN
          DEALLOCATE(sh_env%velocities_tp)
      END IF
      
      IF (ASSOCIATED(sh_env%diabatic_nace_t)) THEN
          DEALLOCATE(sh_env%diabatic_nace_t)
      END IF

      IF (ASSOCIATED(sh_env%diabatic_nace_tp)) THEN
          DEALLOCATE(sh_env%diabatic_nace_tp)
      END IF
      
      IF (ASSOCIATED(sh_env%H_t)) THEN
          DEALLOCATE(sh_env%H_t)
      END IF
      
      IF (ASSOCIATED(sh_env%H_tp)) THEN
          DEALLOCATE(sh_env%H_tp)
      END IF
      
      IF (ASSOCIATED(sh_env%ctin)) THEN
          DEALLOCATE(sh_env%ctin)
      END IF
      
      IF (ASSOCIATED(sh_env%ctin_t)) THEN
          DEALLOCATE(sh_env%ctin_t)
      END IF
      
      IF (ASSOCIATED(sh_env%ctin_before_deco)) THEN
          DEALLOCATE(sh_env%ctin_before_deco)
      END IF
      
      IF (ASSOCIATED(sh_env%ctin_after_deco)) THEN
          DEALLOCATE(sh_env%ctin_after_deco)
      END IF
      
      IF (ASSOCIATED(sh_env%ctin_after_tranf_correction)) THEN
          DEALLOCATE(sh_env%ctin_after_tranf_correction)
      END IF
      
      IF (ASSOCIATED(sh_env%active_region)) THEN
          DEALLOCATE(sh_env%active_region)
      END IF
      
      IF (ASSOCIATED(sh_env%adiab_overlap)) THEN
          DEALLOCATE(sh_env%adiab_overlap)
      END IF
      
      IF (ASSOCIATED(sh_env%adiab_nace_tp)) THEN
          DEALLOCATE(sh_env%adiab_nace_tp)
      END IF
      
      IF (ASSOCIATED(sh_env%adiab_nace_t)) THEN
          DEALLOCATE(sh_env%adiab_nace_t)
      END IF
      
      IF (ASSOCIATED(sh_env%class_forbidden_hop_t)) THEN
          DEALLOCATE(sh_env%class_forbidden_hop_t)
      END IF


      IF (ALLOCATED(sh_env%readcoeff)) THEN
         DEALLOCATE(sh_env%readcoeff)
      END IF

      CALL sh_print_release(sh_env%sh_print_env)

   END SUBROUTINE

   !**********************************************************
   !> \brief Deallocate the sh print environment
   !> \params sh_env => The surface hopping environment
   !**********************************************************
   SUBROUTINE sh_print_release(sh_print_env)
      TYPE(sh_print_env_type)                      :: sh_print_env

      CHARACTER(len=*), PARAMETER :: routineN ='sh_print_release',  &
                                     routineP = moduleN//':'//routineN

      DEALLOCATE(sh_print_env%adiab_energies_t)   

   END SUBROUTINE

END MODULE sh_types
