# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_energy_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_energy_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author fschiff
!> \date   11.06
! *****************************************************************************
MODULE mixed_energy_types

  
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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_energy_types.F" 2

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_energy_types'
  PRIVATE

! *****************************************************************************
  TYPE mixed_energy_type
     REAL ( kind = dp ) :: pot
     REAL ( kind = dp ) :: kin
  END TYPE mixed_energy_type

! *****************************************************************************
  TYPE mixed_force_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER         :: forces
  END TYPE mixed_force_type

! Public data types
  PUBLIC :: mixed_energy_type,&
            mixed_force_type

! Public subroutines
  PUBLIC :: allocate_mixed_energy,&
            deallocate_mixed_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a mixed energy data structure.
!> \param mixed_energy ...
!> \date    11.06
!> \author  fschiff
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_mixed_energy', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(mixed_energy)) THEN
      ALLOCATE (mixed_energy)
    END IF
    CALL init_mixed_energy(mixed_energy)
  END SUBROUTINE allocate_mixed_energy

! *****************************************************************************
!> \brief   Deallocate a mixed energy data structure.
!> \param mixed_energy ...
!> \date    11.06
!> \author  fschiff
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_mixed_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(mixed_energy)) THEN
      DEALLOCATE (mixed_energy)
    END IF
  END SUBROUTINE deallocate_mixed_energy

! *****************************************************************************
!> \brief ...
!> \param mixed_energy ...
! *****************************************************************************
  SUBROUTINE init_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mixed_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(mixed_energy)) THEN
      mixed_energy%pot = 0.0_dp
    ELSE
      CALL cp_abort(cp__l("mixed_energy_types.F",91),&
           "The mixed_energy pointer is not associated "//&
           "and cannot be initialised")
    END IF
  END SUBROUTINE init_mixed_energy

END MODULE mixed_energy_types
