# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-06
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-06
!> - 2014-04, Ole Schuett: generalized into acc-framework
! *****************************************************************************
MODULE dbcsr_acc_operations
  USE ISO_C_BINDING,                   ONLY: C_INT,&
                                             C_PTR
  USE acc_devmem,                      ONLY: acc_devmem_cptr,&
                                             acc_devmem_type
  USE acc_stream,                      ONLY: acc_stream_cptr,&
                                             acc_stream_type
  USE dbcsr_error_handling,            ONLY: dbcsr_abort
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/../../base/base_uses.f90" 1
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
# 26 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_acc_operations'

  LOGICAL, PARAMETER :: careful_mod = .FALSE.


  PUBLIC :: dbcsr_acc_do_mm_stack, dbcsr_acc_transpose


# 73 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F"

CONTAINS

! *****************************************************************************
!> \brief Launch an accelerated kernel for processing a stack.
!> \param param_stack ...
!> \param stack_size ...
!> \param datatype ...
!> \param a_data ...
!> \param b_data ...
!> \param c_data ...
!> \param m_max ...
!> \param n_max ...
!> \param k_max ...
!> \param def_mnk ...
!> \param stream ...
!> \param success ...
! *****************************************************************************
  SUBROUTINE dbcsr_acc_do_mm_stack(param_stack, stack_size, datatype, &
       a_data, b_data, c_data, m_max, n_max, k_max, def_mnk, stream, success)
    TYPE(acc_devmem_type), INTENT(IN)        :: param_stack
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, INTENT(IN)                      :: datatype
    TYPE(acc_devmem_type), INTENT(IN)        :: a_data, b_data
    TYPE(acc_devmem_type), &
      INTENT(INOUT)                          :: c_data
    INTEGER, INTENT(IN)                      :: m_max, n_max, k_max
    LOGICAL, INTENT(IN)                      :: def_mnk
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    LOGICAL, INTENT(INOUT)                   :: success


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(param_stack))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_size))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(datatype))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(a_data))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(b_data))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(c_data))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(m_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(k_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(def_mnk))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(success))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("dbcsr/mm/dbcsr_acc_operations.F",117,"__DBCSR_ACC not compiled in.")
# 149 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F"
  END SUBROUTINE dbcsr_acc_do_mm_stack


! *****************************************************************************
!> \brief Launch an accelerated transpose kernel
!> \param trs_stack ...
!> \param offset ...
!> \param nblks ...
!> \param datatype ...
!> \param buffer ...
!> \param m ...
!> \param n ...
!> \param stream ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_acc_transpose(trs_stack, offset, nblks, datatype, buffer, m, n, stream)
    TYPE(acc_devmem_type), INTENT(IN)        :: trs_stack
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, INTENT(IN)                      :: datatype
    TYPE(acc_devmem_type), INTENT(IN)        :: buffer
    INTEGER, INTENT(IN)                      :: m, n
    TYPE(acc_stream_type), INTENT(IN)        :: stream


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(trs_stack))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(offset))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nblks))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(datatype))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(buffer))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(m))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stream))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("dbcsr/mm/dbcsr_acc_operations.F",182,"__DBCSR_ACC not compiled in.")
# 204 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_acc_operations.F"
  END SUBROUTINE dbcsr_acc_transpose

END MODULE dbcsr_acc_operations
