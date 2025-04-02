# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermal_region_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermal_region_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Thermal regions type: to initialize and control the temperature of
!>        different regions
!> \par History
!>   - Added support for langevin regions (2014/01/08, LT)
!> \author MI
! *****************************************************************************
MODULE thermal_region_types

  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/../base/base_uses.f90" 1
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
# 18 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermal_region_types.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: thermal_regions_type,&
            thermal_region_type,&
            allocate_thermal_regions,&
            release_thermal_regions,&
            retain_thermal_regions

  TYPE thermal_regions_type
     INTEGER :: id_nr, ref_count, nregions
     LOGICAL :: force_rescaling
     REAL(KIND=dp) :: temp_reg0
     LOGICAL, DIMENSION(:), POINTER                   :: do_langevin
     TYPE(section_vals_type), POINTER                 :: section
     TYPE(thermal_region_type), DIMENSION(:), POINTER :: thermal_region
  END TYPE thermal_regions_type

  TYPE thermal_region_type
     INTEGER :: region_index, npart
     INTEGER, DIMENSION(:), POINTER :: part_index
     REAL(KIND=dp) :: ekin, temperature, temp_expected, temp_tol
  END TYPE thermal_region_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermal_region_types'
CONTAINS

! *****************************************************************************
!> \brief allocate thermal_regions
!> \param thermal_regions ...
!> \author
! *****************************************************************************
  SUBROUTINE allocate_thermal_regions(thermal_regions)
    TYPE(thermal_regions_type), POINTER      :: thermal_regions

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_thermal_regions', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check

    check = .NOT.ASSOCIATED(thermal_regions)
    IF(.NOT.(check))CALL cp__a("motion/thermal_region_types.F",60)

    ALLOCATE(thermal_regions)
    thermal_regions%ref_count =  1
    thermal_regions%nregions = 0
    NULLIFY(thermal_regions%thermal_region)
    NULLIFY(thermal_regions%do_langevin)

  END SUBROUTINE allocate_thermal_regions
! *****************************************************************************
!> \brief retains  thermal_regions
!> \param thermal_regions ...
!> \author
! *****************************************************************************
  SUBROUTINE retain_thermal_regions(thermal_regions)

    TYPE(thermal_regions_type), POINTER      :: thermal_regions

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_thermal_regions', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(thermal_regions)) THEN
       IF(.NOT.(thermal_regions%ref_count>0))CALL cp__a("motion/thermal_region_types.F",82)
       thermal_regions%ref_count=thermal_regions%ref_count+1
    END IF

  END SUBROUTINE retain_thermal_regions

! *****************************************************************************
!> \brief release thermal_regions
!> \param thermal_regions ...
!> \author
! *****************************************************************************
  SUBROUTINE release_thermal_regions(thermal_regions)

    TYPE(thermal_regions_type), POINTER      :: thermal_regions

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermal_regions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ireg
    LOGICAL                                  :: check

    check = ASSOCIATED(thermal_regions)
    IF (check) THEN
       check = thermal_regions%ref_count>0
       IF(.NOT.(check))CALL cp__a("motion/thermal_region_types.F",106)
       thermal_regions%ref_count=thermal_regions%ref_count-1
       IF (thermal_regions%ref_count<1) THEN
          IF (ASSOCIATED(thermal_regions%thermal_region)) THEN
            DO ireg = 1,SIZE(thermal_regions%thermal_region)
              DEALLOCATE(thermal_regions%thermal_region(ireg)%part_index)
            END DO
            DEALLOCATE(thermal_regions%thermal_region)
          END IF
          IF (ASSOCIATED(thermal_regions%do_langevin)) THEN
             DEALLOCATE(thermal_regions%do_langevin)
          END IF
          DEALLOCATE(thermal_regions)
       END IF
    END IF

  END SUBROUTINE release_thermal_regions

END MODULE thermal_region_types
