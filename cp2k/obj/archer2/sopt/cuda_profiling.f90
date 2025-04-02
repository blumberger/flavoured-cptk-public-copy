# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  routines for profiling cuda
!> \par History
!>      05.2013 created
!> \author Ole Schuett
! *****************************************************************************
MODULE cuda_profiling
  USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                             C_INT,&
                                             C_NULL_CHAR,&
                                             C_SIZE_T
  USE kinds,                           ONLY: default_string_length,&
                                             int_8

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/../base/base_uses.f90" 1
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
# 20 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F" 2

 !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

 IMPLICIT NONE

 PRIVATE

 PUBLIC  :: cuda_nvtx_init, cuda_nvtx_range_push, cuda_nvtx_range_pop, cuda_mem_info

# 71 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F"


 CONTAINS


! *****************************************************************************
!> \brief ...
! *****************************************************************************
 SUBROUTINE cuda_nvtx_init()
# 91 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F"
  END SUBROUTINE cuda_nvtx_init

! *****************************************************************************
!> \brief ...
!> \param routineN ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_push(routineN)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN




    CALL cp_abort(cp__l("common/cuda_profiling.F",103),"cuda_nvtx_range_push: "//&
          "__CUDA_PROFILING not compiled in, but called with:"//TRIM(routineN))

  END SUBROUTINE cuda_nvtx_range_push

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_pop()




     CALL cp__b("common/cuda_profiling.F",116,"cuda_nvtx_range_push: __CUDA_PROFILING not compiled in.")

  END SUBROUTINE cuda_nvtx_range_pop

! *****************************************************************************
!> \brief ...
!> \param free ...
!> \param total ...
! *****************************************************************************
  SUBROUTINE cuda_mem_info(free, total)
    INTEGER(KIND=int_8), INTENT(OUT)         :: free, total
# 136 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cuda_profiling.F"
    free  = 0
    total = 0

  END SUBROUTINE cuda_mem_info

END MODULE cuda_profiling
