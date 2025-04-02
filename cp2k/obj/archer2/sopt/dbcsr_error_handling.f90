# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/error/dbcsr_error_handling.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/error/dbcsr_error_handling.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
!> \author VW
!>
!> <b>Modification history:</b>
!> - Created Feb 2010
! *****************************************************************************
MODULE dbcsr_error_handling


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/error/../../base/base_uses.f90" 1
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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/error/dbcsr_error_handling.F" 2
  IMPLICIT NONE
  PRIVATE

  ! procedures
  PUBLIC :: dbcsr_assert, dbcsr_abort
  ! parameters
  PUBLIC :: dbcsr_failure_level, dbcsr_fatal_level
  PUBLIC :: dbcsr_caller_error, dbcsr_wrong_args_error,&
            dbcsr_internal_error,&
            dbcsr_unimplemented_error_nr

  ! interfaces
  INTERFACE dbcsr_assert
     MODULE PROCEDURE dbcsr_int_assert
     MODULE PROCEDURE dbcsr_logical_assert
     MODULE PROCEDURE dbcsr_not_assert
     MODULE PROCEDURE dbcsr_true_assert
  END INTERFACE

  INTEGER, PARAMETER :: dbcsr_error_stack_size = 30
  !! level of an error
  INTEGER, PARAMETER :: dbcsr_fatal_level=3
  !! level of a failure
  INTEGER, PARAMETER :: dbcsr_failure_level=2
  !! level of a note
  INTEGER, PARAMETER :: dbcsr_note_level=0
  !! error number: no error
  INTEGER, PARAMETER :: dbcsr_no_error = 0
  !! error number: generic error on the side of the caller
  INTEGER, PARAMETER :: dbcsr_caller_error = 1
  !! error number: one or more arguments have and invalid value
  INTEGER, PARAMETER :: dbcsr_wrong_args_error = 100
  !! error number: precondition failed
  INTEGER, PARAMETER :: dbcsr_precondition_failed = 200
  !! error number: generic error inside the routine
  INTEGER, PARAMETER :: dbcsr_internal_error = -1
  !! error number: postcondition failed
  INTEGER, PARAMETER :: dbcsr_postcondition_failed = -200
  !! error number: invariant failed
  INTEGER, PARAMETER :: dbcsr_invariant_failed = -100
  !! error number: assertion failure
  INTEGER, PARAMETER :: dbcsr_assertion_failed = -300
  !! error number: not implemented
  INTEGER, PARAMETER :: dbcsr_unimplemented_error_nr = -1000


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_error_handling'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param routine ...
!> \param line ...
!> \param msg ...
!> \param[inout]
! *****************************************************************************
  SUBROUTINE dbcsr_abort(routine, line,  msg)
    CHARACTER(*), INTENT(in)                 :: routine
    INTEGER, INTENT(in)                      :: line
    CHARACTER(*), INTENT(in)                 :: msg

    CALL dbcsr_abort_low(dbcsr_fatal_level, dbcsr_internal_error, routine, msg, line)

  END SUBROUTINE dbcsr_abort

! *****************************************************************************
!> \brief ...
!> \param level ...
!> \param etype ...
!> \param routine ...
!> \param msg ...
!> \param line ...
!> \param[inout]
! *****************************************************************************
  SUBROUTINE dbcsr_abort_low (level, etype, routine, msg, line)
    INTEGER, INTENT(in)                      :: level, etype
    CHARACTER(*), INTENT(in)                 :: routine, msg
    INTEGER, INTENT(in)                      :: line

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(level))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(etype))==-1) EXIT ;  END DO ; ENDIF
    CALL cp_abort(cp__l(routine, line), msg)

  END SUBROUTINE dbcsr_abort_low

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
! *****************************************************************************
  SUBROUTINE dbcsr_int_assert(left, rel, right, level, etype, routine, msg, line)
    INTEGER, INTENT(IN)                      :: left
    CHARACTER(len=2), INTENT(IN)             :: rel
    INTEGER, INTENT(IN)                      :: right, level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQ")
       l = left .EQ. right
    CASE ("LT")
       l = left .LT. right
    CASE ("LE")
       l = left .LE. right
    CASE ("GT")
       l = left .GT. right
    CASE ("GE")
       l = left .GE. right
    CASE ("NE")
       l = left .NE. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, 139)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,I9,A4,I9)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line)
    ENDIF
  END SUBROUTINE dbcsr_int_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
! *****************************************************************************
  SUBROUTINE dbcsr_logical_assert(left, rel, right, level, etype, routine, msg, line)
    LOGICAL, INTENT(IN)                      :: left
    CHARACTER(len=*), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQV")
       l = left .EQV. right
    CASE ("NEQV")
       l = left .NEQV. right
    CASE ("OR")
       l = left .OR. right
    CASE ("AND")
       l = left .AND. right
    CASE ("IMP")
       l = .NOT. left .OR. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, 185)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,L1,A,L1)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line)
    ENDIF
  END SUBROUTINE dbcsr_logical_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
! *****************************************************************************
  SUBROUTINE dbcsr_true_assert(right, level, etype, routine, msg, line)
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    l = right
    IF (.NOT.l) THEN
       WRITE(*,'(1X,A,1X,L1)')"ASSERTION FAILED:",&
            right
       CALL dbcsr_abort_low (level, etype, routine, msg, line)
    ENDIF
  END SUBROUTINE dbcsr_true_assert


! *****************************************************************************
!> \brief Assertion
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
! *****************************************************************************
  SUBROUTINE dbcsr_not_assert(rel, right, level, etype, routine, msg, line)
    CHARACTER(len=3), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("NOT")
       l = .NOT. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, 249)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,A,L1)')"ASSERTION FAILED:",&
            "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line)
    ENDIF
  END SUBROUTINE dbcsr_not_assert

END MODULE dbcsr_error_handling
