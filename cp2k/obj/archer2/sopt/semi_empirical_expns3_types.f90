# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_expns3_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_expns3_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the type to handle the 1/R^3 residual integral part
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
MODULE semi_empirical_expns3_types
  
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
# 14 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_expns3_types.F" 2

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \brief 1/R^3 expansion type
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  TYPE semi_empirical_expns3_type
     REAL(KIND=dp)                                  :: core_core
     REAL(KIND=dp), DIMENSION(9)                    :: e1b, e2a
     REAL(KIND=dp), DIMENSION(81)                   :: w
  END TYPE semi_empirical_expns3_type

! *****************************************************************************
!> \brief 1/R^3 expansion type: array of pointers
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  TYPE semi_empirical_expns3_p_type
     TYPE(semi_empirical_expns3_type), POINTER      :: expns3
  END TYPE semi_empirical_expns3_p_type

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_expns3_types'

  PUBLIC :: semi_empirical_expns3_p_type,&
            semi_empirical_expns3_create,&
            semi_empirical_expns3_release

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical 1/R^3 expansion type
!> \param expns3 ...
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE semi_empirical_expns3_create(expns3)
    TYPE(semi_empirical_expns3_type), &
      POINTER                                :: expns3

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_expns3_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(expns3)))CALL cp__a("semi_empirical_expns3_types.F",57)
    ALLOCATE (expns3)
    expns3%core_core = 0.0_dp
    expns3%e1b       = 0.0_dp
    expns3%e2a       = 0.0_dp
    expns3%w         = 0.0_dp
  END SUBROUTINE semi_empirical_expns3_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical type
!> \param expns3 ...
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE semi_empirical_expns3_release(expns3)
    TYPE(semi_empirical_expns3_type), &
      POINTER                                :: expns3

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_expns3_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(expns3)) THEN
       DEALLOCATE (expns3)
    END IF
  END SUBROUTINE semi_empirical_expns3_release

END MODULE semi_empirical_expns3_types
