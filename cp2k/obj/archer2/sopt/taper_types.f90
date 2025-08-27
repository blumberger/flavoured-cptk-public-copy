# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/taper_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/taper_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the semi empirical parameter types.
!> \author Teodoro Laino [tlaino] - 10.2008 University of Zurich
! *****************************************************************************
MODULE taper_types
  
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
# 14 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/taper_types.F" 2

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'taper_types'

! *****************************************************************************
!> \brief Taper type
! *****************************************************************************
  TYPE taper_type
     LOGICAL                               :: apply_taper
     REAL(KIND=dp)                         :: r0, rscale
  END TYPE taper_type

  PUBLIC :: taper_type, taper_create, taper_release, taper_eval, dtaper_eval

CONTAINS

! *****************************************************************************
!> \brief Creates taper type
!> \param taper ...
!> \param rc ...
!> \param range ...
! *****************************************************************************
  SUBROUTINE taper_create(taper, rc, range)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rc, range

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(taper)))CALL cp__a("taper_types.F",48)
    ALLOCATE (taper)
    IF (range > EPSILON(0.0_dp)) THEN
       taper%apply_taper = .TRUE.
       IF(.NOT.(range>0.0_dp))CALL cp__a("taper_types.F",52)
       taper%r0     = 2.0_dp*rc - 20.0_dp * range
       taper%rscale = 1.0_dp/range
    ELSE
       taper%apply_taper = .FALSE.
    END IF

  END SUBROUTINE taper_create

! *****************************************************************************
!> \brief Releases taper type
!> \param taper ...
! *****************************************************************************
  SUBROUTINE taper_release(taper)
    TYPE(taper_type), POINTER                :: taper

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(taper)) THEN
       DEALLOCATE (taper)
    END IF
  END SUBROUTINE taper_release

! *****************************************************************************
!> \brief Taper functions
!> \param taper ...
!> \param rij ...
!> \retval ft ...
! *****************************************************************************
  FUNCTION taper_eval (taper, rij) RESULT(ft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    REAL(KIND=dp)                            :: ft

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr

    ft = 1._dp
    IF (taper%apply_taper) THEN
       dr = taper%rscale*(rij-taper%r0)
       ft = 0.5_dp*(1.0_dp-TANH(dr))
    END IF
  END FUNCTION taper_eval

! *****************************************************************************
!> \brief Analytical derivatives for taper function
!> \param taper ...
!> \param rij ...
!> \retval dft ...
! *****************************************************************************
  FUNCTION dtaper_eval (taper, rij) RESULT(dft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    REAL(KIND=dp)                            :: dft

    CHARACTER(len=*), PARAMETER :: routineN = 'dtaper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr

    dft = 0.0_dp
    IF (taper%apply_taper) THEN
       dr  = taper%rscale*(rij-taper%r0)
       dft = -0.5_dp*(1.0_dp-TANH(dr)**2)*taper%rscale
    END IF
  END FUNCTION dtaper_eval

END MODULE taper_types
