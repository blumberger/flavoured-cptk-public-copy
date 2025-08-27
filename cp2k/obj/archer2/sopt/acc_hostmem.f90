# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_hostmem



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
# 27 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_hostmem'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.


  PUBLIC :: acc_hostmem_allocate, acc_hostmem_deallocate


  INTERFACE acc_hostmem_allocate
     MODULE PROCEDURE acc_hostmem_alloc_i,    acc_hostmem_alloc_l
     MODULE PROCEDURE acc_hostmem_alloc_r,    acc_hostmem_alloc_d
     MODULE PROCEDURE acc_hostmem_alloc_c,    acc_hostmem_alloc_z
     MODULE PROCEDURE acc_hostmem_alloc_i_2D, acc_hostmem_alloc_l_2D
     MODULE PROCEDURE acc_hostmem_alloc_r_2D, acc_hostmem_alloc_d_2D
     MODULE PROCEDURE acc_hostmem_alloc_c_2D, acc_hostmem_alloc_z_2D
  END INTERFACE

  INTERFACE acc_hostmem_deallocate
     MODULE PROCEDURE acc_hostmem_dealloc_i,    acc_hostmem_dealloc_l
     MODULE PROCEDURE acc_hostmem_dealloc_r,    acc_hostmem_dealloc_d
     MODULE PROCEDURE acc_hostmem_dealloc_c,    acc_hostmem_dealloc_z
     MODULE PROCEDURE acc_hostmem_dealloc_i_2D, acc_hostmem_dealloc_l_2D
     MODULE PROCEDURE acc_hostmem_dealloc_r_2D, acc_hostmem_dealloc_d_2D
     MODULE PROCEDURE acc_hostmem_dealloc_c_2D, acc_hostmem_dealloc_z_2D
  END INTERFACE



# 84 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F"


CONTAINS


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_i.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_i (host_mem, n, stream)
    INTEGER(KIND=int_4), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_i: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_i



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_i_2D (host_mem, n1, n2, stream)
    INTEGER(KIND=int_4), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_i.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_i_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_i_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i (host_mem, stream)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_i: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_i


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i_2D (host_mem, stream)
    INTEGER(KIND=int_4), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_i: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_i_2D
# 89 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_l.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_l (host_mem, n, stream)
    INTEGER(KIND=int_8), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_l: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_l



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_l_2D (host_mem, n1, n2, stream)
    INTEGER(KIND=int_8), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_l.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_l_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_l_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_l (host_mem, stream)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_l', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_l: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_l


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_l_2D (host_mem, stream)
    INTEGER(KIND=int_8), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_l_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_l: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_l_2D
# 90 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_r.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_r (host_mem, n, stream)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_r: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_r



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_r_2D (host_mem, n1, n2, stream)
    REAL(kind=real_4), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_r.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_r_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_r_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r (host_mem, stream)
    REAL(kind=real_4), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_r: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_r


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r_2D (host_mem, stream)
    REAL(kind=real_4), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_r: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_r_2D
# 91 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_d.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_d (host_mem, n, stream)
    REAL(kind=real_8), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_d: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_d



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_d_2D (host_mem, n1, n2, stream)
    REAL(kind=real_8), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_d.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_d_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_d_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_d (host_mem, stream)
    REAL(kind=real_8), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_d', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_d: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_d


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_d_2D (host_mem, stream)
    REAL(kind=real_8), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_d_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_d: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_d_2D
# 92 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_c.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c (host_mem, n, stream)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_c: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_c



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c_2D (host_mem, n1, n2, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_c.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_c_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_c_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_c: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_c


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c_2D (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_c: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_c_2D
# 93 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_z.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_z (host_mem, n, stream)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream






    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",27,"acc_hostmem_alloc_z: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_z



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_z_2D (host_mem, n1, n2, stream)
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
# 53 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem_z.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",57,"acc_hostmem_alloc_z_2D: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_alloc_z_2D


! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_z (host_mem, stream)
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_z', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",81,"acc_hostmem_dealloc_z: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_z


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_z_2D (host_mem, stream)
    COMPLEX(kind=real_8), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_z_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN



    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(host_mem))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_hostmem.F",105,"acc_hostmem_dealloc_z: ACC not compiled in.")

  END SUBROUTINE acc_hostmem_dealloc_z_2D
# 94 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F" 2


! *****************************************************************************
!> \brief Helper-routine performing allocation of host-pinned cuda memory.
!> \param host_mem_c_ptr pointer to allocated memory
!> \param n_bytes number of bytes to allocate
!> \param stream ...
! *****************************************************************************
# 125 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F"

# 152 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_hostmem.F"

END MODULE acc_hostmem
