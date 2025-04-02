# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_period_efield_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_period_efield_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief type for berry phase efield matrices. At the moment only used for
!>        cosmat and sinmat
!> \par History
!>      none
!> \author fschiff (06.2010)
! *****************************************************************************

MODULE qs_period_efield_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
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
# 20 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_period_efield_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_period_efield_types'

  PUBLIC :: efield_berry_type, efield_berry_release, init_efield_matrices,&
            set_efield_matrices

  TYPE efield_berry_type
     REAL(KIND=dp)                                         :: field_energy
     REAL(KIND=dp), DIMENSION(3)                           :: polarisation
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: cosmat
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: sinmat
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: dipmat
  END TYPE efield_berry_type

CONTAINS


! *****************************************************************************
!> \brief ...
!> \param efield ...
! *****************************************************************************
  SUBROUTINE init_efield_matrices(efield)
    TYPE(efield_berry_type), POINTER         :: efield

    CHARACTER(len=*), PARAMETER :: routineN = 'init_efield_matrices', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: field_energy
    REAL(KIND=dp), DIMENSION(3)              :: polarisation

! retain possible values for energy and polarisation

    IF(ASSOCIATED(efield)) THEN
       field_energy = efield%field_energy
       polarisation = efield%polarisation
       CALL efield_berry_release(efield)
    ELSE
       field_energy = 0.0_dp
       polarisation = 0.0_dp
    END IF

    ALLOCATE(efield)
    NULLIFY(efield%cosmat)
    NULLIFY(efield%sinmat)
    NULLIFY(efield%dipmat)

    efield%field_energy = field_energy
    efield%polarisation = polarisation

  END SUBROUTINE init_efield_matrices

! *****************************************************************************
!> \brief ...
!> \param efield ...
!> \param sinmat ...
!> \param cosmat ...
!> \param dipmat ...
! *****************************************************************************
  SUBROUTINE set_efield_matrices(efield,sinmat,cosmat,dipmat)

    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sinmat, cosmat, dipmat

    CHARACTER(len=*), PARAMETER :: routineN = 'set_efield_matrices', &
      routineP = moduleN//':'//routineN

     IF(PRESENT(cosmat))efield%cosmat=>cosmat
     IF(PRESENT(sinmat))efield%sinmat=>sinmat
     IF(PRESENT(dipmat))efield%dipmat=>dipmat

  END SUBROUTINE set_efield_matrices

! *****************************************************************************
!> \brief ...
!> \param efield ...
! *****************************************************************************
  SUBROUTINE efield_berry_release(efield)
    TYPE(efield_berry_type), POINTER         :: efield

    CHARACTER(len=*), PARAMETER :: routineN = 'efield_berry_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(efield))THEN
       IF(ASSOCIATED(efield%sinmat).AND.ASSOCIATED(efield%cosmat))THEN
          CALL cp_dbcsr_deallocate_matrix_set ( efield%cosmat)
          CALL cp_dbcsr_deallocate_matrix_set ( efield%sinmat)
       END IF
       IF(ASSOCIATED(efield%dipmat))THEN
          CALL cp_dbcsr_deallocate_matrix_set ( efield%dipmat)
       END IF
       DEALLOCATE(efield)
    END IF
  END SUBROUTINE efield_berry_release

END MODULE qs_period_efield_types
