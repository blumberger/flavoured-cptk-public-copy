# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_event



  USE acc_stream,                      ONLY: acc_stream_cptr,&
                                             acc_stream_type

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/../base/base_uses.f90" 1
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
# 18 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_event'

  PUBLIC :: acc_event_type
  PUBLIC :: acc_event_create, acc_event_destroy
  PUBLIC :: acc_event_record, acc_event_query
  PUBLIC :: acc_stream_wait_event, acc_event_synchronize


  TYPE acc_event_type
    PRIVATE



    INTEGER :: dummy = 1

  END TYPE acc_event_type

# 96 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"

CONTAINS

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamWaitEvent.
!>        Because of fortran circular dependency restriction this can not go into acc_stream.F
!> \param[in] stream stream
!> \param[in] event event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_wait_event(stream, event)
    TYPE(acc_stream_type), INTENT(IN) :: stream
    TYPE(acc_event_type), INTENT(IN)  :: event


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(event))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",113,"__ACC not compiled in.")
# 127 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END SUBROUTINE acc_stream_wait_event


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventRecord.
!> \param[in] this event
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_record(this, stream)
    TYPE(acc_event_type), INTENT(IN)  :: this
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",143,"__ACC not compiled in.")
# 157 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END SUBROUTINE acc_event_record


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventCreate.
!> \param[in,out] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_create(this)
    TYPE(acc_event_type), &
      INTENT(INOUT)                          :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",171,"__ACC not compiled in.")
# 181 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END SUBROUTINE acc_event_create


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventDestroy.
!> \param[in,out] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_destroy(this)
    TYPE(acc_event_type), &
      INTENT(INOUT)                          :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",195,"__ACC not compiled in.")
# 205 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END SUBROUTINE acc_event_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventQuery.
!> \param[in] this event
!> \retval res true if event has occured, false otherwise
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_event_query(this) RESULT(res)
    TYPE(acc_event_type), INTENT(IN)         :: this
    LOGICAL                                  :: res


    res = .FALSE.
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",221,"__ACC not compiled in.")
# 231 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END FUNCTION acc_event_query


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventSynchronize.
!> \param[in] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_synchronize(this)
    TYPE(acc_event_type), INTENT(IN)  :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_event.F",244,"__ACC not compiled in.")
# 253 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_event.F"
END SUBROUTINE acc_event_synchronize

END MODULE acc_event
