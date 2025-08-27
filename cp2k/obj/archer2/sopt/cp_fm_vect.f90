# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_vect.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_vect.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routine to handle vectors of full matrixes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/../base/base_uses.f90" 1
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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_vect.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_vect'

  PUBLIC :: cp_fm_vect_dealloc, cp_fm_vect_copy
!***
CONTAINS

! *****************************************************************************
!> \brief deallocate an array of pointers to blacs matrixes
!> \param matrixes the array of matrixes to deallocate
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_fm_vect_dealloc(matrixes)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(matrixes)) THEN
       DO i=1,SIZE(matrixes)
          CALL cp_fm_release(matrixes(i)%matrix)
       END DO
       DEALLOCATE(matrixes)
    END IF
  END SUBROUTINE cp_fm_vect_dealloc

! *****************************************************************************
!> \brief Does a shallow copy of an array of full matrices (i.e. just retains
!>      the matrices)
!> \param matrixes the matrixes to copy
!> \param copy ...
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_vect_copy(matrixes, copy)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: matrixes
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: copy

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

  ALLOCATE(copy(SIZE(matrixes)))
  DO i=1,SIZE(matrixes)
     copy(i)%matrix => matrixes(i)%matrix
     CALL cp_fm_retain(matrixes(i)%matrix)
  END DO
END SUBROUTINE cp_fm_vect_copy

END MODULE cp_fm_vect
