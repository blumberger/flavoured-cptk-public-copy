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

MODULE aom_types
  USE kinds,                           ONLY: default_string_length, &
                                             dp
  USE input_section_types,             ONLY: section_vals_type

  USE parallel_rng_types,              ONLY: rng_stream_type


IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aom_types'



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
!        made. Feel free to hard-code it into the aom_env_type.
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
   REAL(KIND=dp), dimension(:,:), allocatable :: Ttt
   REAL(KIND=dp), dimension(:,:), allocatable :: Sttsqrt_var

END TYPE adiab_nace_env_type

!COMMENT: one should check that all those variables are actually used somewhere
!in the code
   TYPE aom_env_type
!> Constant over the simulations
     integer                                :: atoms !Number of active atoms
     integer                                :: natom !Total number of atoms
     REAL(KIND=dp)                          :: deltat !Timestep
     REAL(KIND = dp)                        :: scaling, cbar !For AOM calculations

!> Information for the overlap method: calculate at the first timestep
     INTEGER, DIMENSION(:), POINTER           :: elements, atomlist
     INTEGER, DIMENSION(:,:,:,:), POINTER     :: neighbourlist
     integer, dimension(:,:), POINTER         :: connlist1, connlist2, bneighbour, aneighbour
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: pcoeffs
     REAL(KIND=dp), DIMENSION(:), POINTER     :: scoeffs
     LOGICAL, DIMENSION(:), POINTER           :: ext_pres
     TYPE(rng_stream_type), POINTER           :: rng_stream !the random seed generator

!> Calculate at each timestep
     REAL(KIND=dp), dimension(:,:), POINTER         :: coms
     REAL(KIND=dp), dimension(:,:,:), POINTER       :: crds
     REAL(KIND=dp), dimension(:,:), POINTER         :: dist
     INTEGER, DIMENSION(:,:), POINTER               :: NN
     REAL(KIND=dp), dimension(:,:,:), allocatable   :: pcoeffsrun
     REAL(KIND=dp), dimension(:,:), POINTER         :: At, Atm, Atp, Bt, Btm, & 
                                                       Btp, siteenergies
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: nacv!, nacv_t
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: F
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F
     REAL(KIND = dp), DIMENSION(:,:), POINTER       :: forces_adiab, &      ! forces at step 0
                                                       forces_adiab_new!, &

!> Use to store the present info in the future timestep
     REAL(KIND = dp), DIMENSION(:,:), POINTER     :: rm_t,&        ! position step min 1
                                                     rm_tm,&       ! position step min 2
                                                     velocities_tp ! velocities step
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: C_t, C_tp
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: diabatic_nace_tp, & 
                                                    diabatic_nace_t  !othogonal diab NACE
     REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE   :: E_t, E_tp
     COMPLEX(KIND=dp), DIMENSION(:,:), POINTER    :: H_t, H_tp
     REAL(KIND=dp), dimension(:,:), POINTER       :: Stt, Stmt, Sttm, Stmtm, &
                                                     Stptp, Sttp, Stpt
     COMPLEX(KIND=dp), DIMENSION(:),POINTER  :: ctin, ctin_t !wavefunction at tp and t
     COMPLEX(KIND=dp), dimension(:), allocatable  :: readcoeff
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
     INTEGER                                :: new_state, old_state, active_state
     REAL(KIND=dp), dimension(:,:), ALLOCATABLE  :: adiab_overlap 
     LOGICAL                                :: attempt_hop
     LOGICAL                                :: successful_hop
     LOGICAL                                :: do_reorder
     LOGICAL                                :: no_diab_nace
!> multiple time step
     INTEGER                                      :: multiple_time_step

!> Required for HS propagation
     LOGICAL                                :: HS_propagation
     REAL(KIND=dp), dimension(:,:), POINTER :: adiab_nace_tp, adiab_nace_t
     INTEGER, dimension(:), POINTER         :: class_forbidden_hop_t

!> To print
     TYPE(sh_print_env_type)                :: sh_print_env

!> CHANGE ME, added nace_allocatables
    TYPE(adiab_nace_env_type)                   :: adiab_nace_env

   END TYPE aom_env_type

   PUBLIC :: aom_env_type, &
             sh_print_env_type


END MODULE aom_types
