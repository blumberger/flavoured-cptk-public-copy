# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/et_coupling_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/et_coupling_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the et_coupling data type.
!> \author Florian Schiffmann (01.2007,fschiff)
! *****************************************************************************
MODULE et_coupling_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/et_coupling_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_types'

! *** Public data types ***

  PUBLIC :: et_coupling_type

! *** Public subroutines ***

  PUBLIC :: et_coupling_create,&
            et_coupling_release,&
            set_et_coupling_type

! *****************************************************************************
!> \par History
!>      01.2007 created [Florian Schiffmann]
!> \author fschiff
! *****************************************************************************
  TYPE et_coupling_type
    TYPE(cp_fm_p_type),DIMENSION(:), POINTER           :: et_mo_coeff
    TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER     :: rest_mat
    LOGICAL                                            :: first_run
    LOGICAL                                            :: keep_matrix
    REAL(KIND = dp)                                    :: energy,e1,order_p
  END TYPE

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
! *****************************************************************************
    SUBROUTINE et_coupling_create(et_coupling)
    TYPE(et_coupling_type), POINTER          :: et_coupling

    CHARACTER(len=*), PARAMETER :: routineN = 'et_coupling_create', &
      routineP = moduleN//':'//routineN

      ALLOCATE(et_coupling)

      NULLIFY(et_coupling%et_mo_coeff)
      NULLIFY(et_coupling%rest_mat)
      et_coupling%first_run=.TRUE.
      et_coupling%keep_matrix=.FALSE.
      ALLOCATE(et_coupling%rest_mat(2))

    END SUBROUTINE et_coupling_create

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
! *****************************************************************************
    SUBROUTINE get_et_coupling_type(et_coupling,et_mo_coeff,rest_mat)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: et_mo_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rest_mat

    CHARACTER(len=*), PARAMETER :: routineN = 'get_et_coupling_type', &
      routineP = moduleN//':'//routineN

      IF(PRESENT(et_mo_coeff))et_mo_coeff => et_coupling%et_mo_coeff
      IF(PRESENT(rest_mat))rest_mat => et_coupling%rest_mat

    END SUBROUTINE get_et_coupling_type

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
! *****************************************************************************
    SUBROUTINE set_et_coupling_type(et_coupling,et_mo_coeff,rest_mat)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: et_mo_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rest_mat

    CHARACTER(len=*), PARAMETER :: routineN = 'set_et_coupling_type', &
      routineP = moduleN//':'//routineN

      IF(PRESENT(et_mo_coeff))  et_coupling%et_mo_coeff  = et_mo_coeff
      IF(PRESENT(rest_mat)) et_coupling%rest_mat => rest_mat

    END SUBROUTINE set_et_coupling_type

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
! *****************************************************************************
    SUBROUTINE et_coupling_release(et_coupling)
    TYPE(et_coupling_type), POINTER          :: et_coupling

    CHARACTER(LEN=*), PARAMETER :: routineN = 'et_coupling_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

      IF(ASSOCIATED(et_coupling%et_mo_coeff))THEN
         DO i=1,SIZE(et_coupling%et_mo_coeff)
            CALL cp_fm_release(et_coupling%et_mo_coeff(i)%matrix)
         END DO
         DEALLOCATE(et_coupling%et_mo_coeff)
      END IF
      IF(ASSOCIATED(et_coupling%rest_mat))THEN
!         CALL deallocate_matrix_set(et_coupling%rest_mat)
         DEALLOCATE(et_coupling%rest_mat)
      END IF

      DEALLOCATE(et_coupling)
    END SUBROUTINE et_coupling_release

END MODULE et_coupling_types

