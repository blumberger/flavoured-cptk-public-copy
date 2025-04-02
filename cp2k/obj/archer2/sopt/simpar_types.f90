# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/simpar_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/simpar_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Type for storing MD parameters
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! *****************************************************************************
MODULE simpar_types
  
  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/./base/base_uses.f90" 1
! Basic use statements and preprocessor macros
! should be included in the use statements

  USE base_hooks,                      ONLY: cp__a,&
                                             cp__b,&
                                             cp__w,&
                                             cp__l,&
                                             cp_abort,&
                                             cp_warn,&
                                             timeset,&
                                             timestop


! Dangerous: Full path can be arbitrarily long and might overflow Fortran line.









! The MARK_USED macro can be used to mark an argument/variable as used.
! It is intended to make it possible to switch on -Werror=unused-dummy-argument,
! but deal elegantly with e.g. library wrapper routines that take arguments only used if the library is linked in. 
! This code should be valid for any Fortran variable, is always standard conforming,
! and will be optimized away completely by the compiler
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/simpar_types.F" 2

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
!> \brief Simulation parameter type for molecular dynamics
!> \par History
!>         created [CJM]
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! *****************************************************************************
  TYPE simpar_type
     INTEGER        :: nsteps
     REAL (KIND=dp) :: dt
     REAL (KIND=dp) :: dt_fact
     REAL (KIND=dp) :: dr_tol
     REAL (KIND=dp) :: dsc_tol
     REAL (KIND=dp) :: temp_ext
     REAL (KIND=dp) :: temp_baro_ext
     REAL (KIND=dp) :: temp_baro
     REAL (KIND=dp) :: temp_tol
     REAL (KIND=dp) :: temp_baro_tol
     REAL (KIND=dp) :: p_ext
     REAL (KIND=dp) :: cmass
     REAL (KIND=dp) :: cmass_nph
     REAL (KIND=dp) :: v0
     REAL (KIND=dp) :: e0
     REAL (KIND=dp) :: v_shock
     REAL (KIND=dp) :: p0
     REAL (KIND=dp) :: f_annealing
     REAL (KIND=dp) :: f_annealing_cell
     REAL (KIND=dp) :: gamma_nph
     INTEGER        :: ensemble
     LOGICAL        :: constraint
     LOGICAL        :: annealing
     LOGICAL        :: annealing_cell
     LOGICAL        :: dump_lm
     LOGICAL        :: angvel_zero
     LOGICAL        :: variable_dt
     INTEGER        :: nfree, nfree_rot_transl
     INTEGER        :: info_constraint
     INTEGER        :: lagrange_multipliers
     REAL (KIND=dp) :: tau_cell
     ! Constraints Parameters
     REAL (KIND=dp) :: shake_tol, roll_tol
     ! Langevin Parameters
     REAL (KIND=dp) :: gamma
     REAL (KIND=dp) :: noisy_gamma
     REAL (KIND=dp) :: shadow_gamma
     REAL (KIND=dp) :: var_w
     ! RESPA Parameters
     LOGICAL        :: multi_time_switch, do_respa
     INTEGER        :: n_time_steps
     ! SHELL parameters
     REAL (KIND=dp) :: temp_sh_ext
     REAL (KIND=dp) :: temp_sh_tol
     LOGICAL        :: temperature_per_kind
     LOGICAL        :: scale_temperature_per_kind
     LOGICAL        :: do_thermal_region
     ! ADIABATIC parameters
     REAL (KIND=dp) :: temp_slow
     REAL (KIND=dp) :: temp_fast
     REAL (KIND=dp) :: temp_tol_fast, temp_tol_slow
     INTEGER :: n_resp_fast
     ! Velocity softening Parameters
     INTEGER        :: soften_nsteps
     REAL (KIND=dp) :: soften_alpha
     REAL (KIND=dp) :: soften_delta
  END TYPE simpar_type


  PUBLIC :: simpar_type,&
            create_simpar_type,&
            release_simpar_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'simpar_types'

CONTAINS
! *****************************************************************************
!> \brief Creates the simulation parameters type
!> \param simpar ...
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_simpar_type(simpar)
    TYPE(simpar_type), POINTER               :: simpar

    CHARACTER(len=*), PARAMETER :: routineN = 'create_simpar_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(simpar)))CALL cp__a("simpar_types.F",105)
    ALLOCATE(simpar)
  END SUBROUTINE create_simpar_type

! *****************************************************************************
!> \brief Releases the simulation parameters type
!> \param simpar ...
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE release_simpar_type(simpar)
    TYPE(simpar_type), POINTER               :: simpar

    CHARACTER(len=*), PARAMETER :: routineN = 'release_simpar_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(simpar)))CALL cp__a("simpar_types.F",120)
    DEALLOCATE(simpar)
  END SUBROUTINE release_simpar_type

END MODULE simpar_types
