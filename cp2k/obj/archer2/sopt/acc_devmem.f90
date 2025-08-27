# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
! *****************************************************************************
MODULE acc_devmem



  USE acc_kinds,                       ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_cptr,&
                                             acc_stream_synchronize,&
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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_devmem'


  PUBLIC :: acc_devmem_type
  PUBLIC :: acc_devmem_allocate_bytes, acc_devmem_deallocate
  PUBLIC :: acc_devmem_setzero_bytes
  PUBLIC :: acc_devmem_allocated
  PUBLIC :: acc_devmem_dev2host, acc_devmem_host2dev
  PUBLIC :: acc_devmem_size_in_bytes
  PUBLIC :: acc_devmem_ensure_size_bytes
  PUBLIC :: acc_devmem_cptr

 INTERFACE acc_devmem_dev2host
  MODULE PROCEDURE dev2host_i4_1D
  MODULE PROCEDURE dev2host_i8_1D
  MODULE PROCEDURE dev2host_r4_1D
  MODULE PROCEDURE dev2host_r8_1D
  MODULE PROCEDURE dev2host_c4_1D
  MODULE PROCEDURE dev2host_c8_1D
 END INTERFACE acc_devmem_dev2host

 INTERFACE acc_devmem_host2dev
  MODULE PROCEDURE host2dev_i4_1D
  MODULE PROCEDURE host2dev_i8_1D
  MODULE PROCEDURE host2dev_r4_1D
  MODULE PROCEDURE host2dev_r8_1D
  MODULE PROCEDURE host2dev_c4_1D
  MODULE PROCEDURE host2dev_c8_1D
  MODULE PROCEDURE host2dev_i4_2D
  MODULE PROCEDURE host2dev_i8_2D
  MODULE PROCEDURE host2dev_r4_2D
  MODULE PROCEDURE host2dev_r8_2D
  MODULE PROCEDURE host2dev_c4_2D
  MODULE PROCEDURE host2dev_c8_2D
 END INTERFACE acc_devmem_host2dev

 TYPE acc_devmem_type
   PRIVATE
   INTEGER                      :: size_in_bytes = -1



 END TYPE acc_devmem_type

# 163 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"

CONTAINS


! *****************************************************************************
!> \brief Ensures that given devmem has at least the requested size.
!> \param[in,out] this device memory
!> \param[in] stream on which zeroing and memcopying is performed
!> \param[in] requested_size_in_bytes requested size in bytes
!> \param[in] nocopy (optional) if after growin old content should NOT be copied over. Default: false. 
!> \param[in] zero_pad (optional) if after growing the new memory should be zeroed. Default: false. 
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_ensure_size_bytes(this, stream, requested_size_in_bytes, nocopy, zero_pad)
    TYPE(acc_devmem_type), &
      INTENT(INOUT)                          :: this
    TYPE(acc_stream_type), INTENT(IN) :: stream
    INTEGER, INTENT(IN)                      :: requested_size_in_bytes
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(requested_size_in_bytes))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nocopy))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(zero_pad))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",189,"__ACC not compiled in.")
# 237 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
  END SUBROUTINE acc_devmem_ensure_size_bytes

! *****************************************************************************
!> \brief Returns a logical, which indicates if the given devmem is allocated.
!> \param[in] this device memory
!> \retval res true if device memory is allocated, false otherwise
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION acc_devmem_allocated(this) RESULT(res)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    LOGICAL                                  :: res

     res = this%size_in_bytes >= 0
  END FUNCTION acc_devmem_allocated


! *****************************************************************************
!> \brief Returns size of given devmem in terms of item count (not bytes!)
!> \param[in] this device memory
!> \retval res size of device memory (item count)
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION acc_devmem_size_in_bytes(this) RESULT(res)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    INTEGER                                  :: res

     IF(this%size_in_bytes < 0)&
        CALL cp__b("acc/acc_devmem.F",264,"acc_devmem_len: not allocated")
     res = this%size_in_bytes
  END FUNCTION acc_devmem_size_in_bytes



! *****************************************************************************
!> \brief Returns C-pointer to data of given devmem.
!> \param[in] this device memory
!> \retval res false (accelerator support is not enabled)
!> \author  Ole Schuett
! *****************************************************************************

  FUNCTION acc_devmem_cptr(this) RESULT(res)
    INTEGER, INTENT(IN)                      :: this
    LOGICAL                                  :: res

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    res = .FALSE.
  END FUNCTION acc_devmem_cptr
# 300 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"


! *****************************************************************************
!> \brief Fortran-wrapper for cudaMemGetInfo.
!> \param[out] free free device memory
!> \param[out] avail available device memory
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE acc_devmem_info(free, avail)
    INTEGER, INTENT(OUT)                     :: free, avail


    avail=-1; free=-1 ! assign intent-out arguments to silence compiler warnings
    CALL cp__b("acc/acc_devmem.F",313,"__ACC not compiled in.")






 END SUBROUTINE acc_devmem_info


! *****************************************************************************
!> \brief Allocates a given devmem.
!> \param[in,out] this device memory
!> \param[in] size_in_bytes size in bytes
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_allocate_bytes(this, size_in_bytes)
    TYPE(acc_devmem_type), INTENT(INOUT)     :: this
    INTEGER, INTENT(IN)                      :: size_in_bytes


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(size_in_bytes))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",336,"__ACC not compiled in.")
# 352 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
  END SUBROUTINE acc_devmem_allocate_bytes


! *****************************************************************************
!> \brief Deallocates a given devmem.
!> \param[in,out] this device memory
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_deallocate(this)
    TYPE(acc_devmem_type), INTENT(INOUT) :: this


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",365,"__ACC not compiled in.")
# 383 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
  END SUBROUTINE acc_devmem_deallocate


! *****************************************************************************
!> \brief Sets entries in given devmem to zero, asynchronously.
!> \param[in,out] this device memory
!> \param[in] first_byte (optional) begin of region to zero, defaults to 1 if not given. 
!> \param[in] last_byte (optional) end of region to zero, defaults to size if not given. 
!> \param[in] stream stream on which zeroing is performed.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_setzero_bytes(this, first_byte, last_byte, stream)
    TYPE(acc_devmem_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN), OPTIONAL        :: first_byte, last_byte
    TYPE(acc_stream_type), INTENT(IN)    :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(first_byte))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(last_byte))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",404,"__ACC not compiled in.")
# 426 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"
  END SUBROUTINE acc_devmem_setzero_bytes



! *****************************************************************************
!> \brief Helper-routine performing actuall host2dev transfers.
!> \param[in] this device memory
!> \param hostmem_cptr host memory
!> \param[in] n_bytes number of bytes
!> \param[in] stream stream used for memory transfer
!> \author  Ole Schuett
! *****************************************************************************
# 463 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"


! *****************************************************************************
!> \brief Helper-routine performing actual dev2host transfers.
!> \param[in] this device memory
!> \param hostmem_cptr host memory
!> \param[in] n_bytes number of bytes
!> \param[in] stream stream used for memory transfer
!> \author  Ole Schuett
! *****************************************************************************
# 499 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F"



# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_i4.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_i4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_4), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_i4_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_i4_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_4), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_i4_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_i4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_4), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_i4_1D

# 502 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_i8.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_i8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_8), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_i8_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_i8_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_8), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_i8_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_i8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    INTEGER(kind=int_8), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_i8_1D

# 503 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_r4.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_r4_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r4_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_r4_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_r4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_r4_1D

# 504 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_r8.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_8), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_r8_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r8_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_8), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_r8_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_r8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_8), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_r8_1D

# 505 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_c4.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_c4_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c4_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_c4_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_c4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_c4_1D

# 506 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem_c8.f90" 1
! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",17,"__ACC not compiled in.")



 END SUBROUTINE host2dev_c8_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c8_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",40,"__ACC not compiled in.")



 END SUBROUTINE host2dev_c8_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_c8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(this))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hostmem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_devmem.F",63,"__ACC not compiled in.")



 END SUBROUTINE dev2host_c8_1D

# 507 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_devmem.F" 2

END MODULE acc_devmem
