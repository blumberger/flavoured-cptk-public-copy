# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_stream




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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_stream'

  PUBLIC :: acc_stream_type
  PUBLIC :: acc_stream_create, acc_stream_destroy
  PUBLIC :: acc_stream_synchronize
  PUBLIC :: acc_stream_priority_range
  PUBLIC :: acc_stream_equal, acc_stream_associated
  PUBLIC :: acc_stream_cptr

  TYPE acc_stream_type
    PRIVATE



    INTEGER :: dummy = 1

  END TYPE acc_stream_type

# 79 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
CONTAINS




! *****************************************************************************
!> \brief Returns C-pointer of given stream.
!> \param[in] this stream ID
!> \retval res false (accelerator support is not enabled)
!> \author  Ole Schuett
! *****************************************************************************

FUNCTION acc_stream_cptr(this) RESULT(res)
    INTEGER, INTENT(in)                      :: this
    LOGICAL                                  :: res

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    res = .FALSE.
END FUNCTION acc_stream_cptr

# 113 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamCreate.
!> \param[out] this stream
!> \param[in] name stream name
!> \param[in] priority (optional) stream priority
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_create(this, name, priority)
    TYPE(acc_stream_type),INTENT(OUT) :: this
    CHARACTER(LEN=*), INTENT(IN)             :: name
    INTEGER, INTENT(IN), OPTIONAL            :: priority


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(name))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(priority))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_stream.F",131,"__ACC not compiled in.")
# 147 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
END SUBROUTINE acc_stream_create

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamDestroy.
!> \param[in,out] this stream
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_destroy(this)
    TYPE(acc_stream_type), &
      INTENT(INOUT)                          :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_stream.F",160,"__ACC not compiled in.")
# 170 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
END SUBROUTINE acc_stream_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamSynchronize.
!> \param[in] this stream
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_synchronize(this)
    TYPE(acc_stream_type), &
      INTENT(IN)                             :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_stream.F",184,"__ACC not compiled in.")
# 193 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_stream.F"
END SUBROUTINE acc_stream_synchronize


! *****************************************************************************
!> \brief Fortran-wrapper for cudaDeviceGetStreamPriorityRange.
!> \param[out] least least priority
!> \param[out] greatest greatest priroity
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_priority_range(least, greatest)
    INTEGER, INTENT(OUT)                     :: least, greatest


    least=-1; greatest=-1 ! assign intent-out arguments to silence compiler warnings
    CALL cp__b("acc/acc_stream.F",207,"__ACC not compiled in.")






END SUBROUTINE acc_stream_priority_range


! *****************************************************************************
!> \brief Checks if two streams are equal
!> \param[in] this first stream
!> \param[in] other second stream
!> \retval res true if equal, false otherwise
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_stream_equal(this, other) RESULT(res)
    TYPE(acc_stream_type), INTENT(IN) :: this, other
    LOGICAL                                  :: res

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(other))==-1) EXIT ;  END DO ; ENDIF
    res = .TRUE.



END FUNCTION acc_stream_equal


! *****************************************************************************
!> \brief Checks if a streams is associated
!> \param[in] this stream
!> \retval res true if associated, false otherwise
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_stream_associated(this) RESULT(res)
    TYPE(acc_stream_type), INTENT(IN) :: this
    LOGICAL                                  :: res

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    res = .FALSE.



END FUNCTION acc_stream_associated

END MODULE acc_stream
