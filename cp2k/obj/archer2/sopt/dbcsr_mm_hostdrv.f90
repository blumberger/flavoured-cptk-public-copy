# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Stacks of small matrix multiplications
!> \author  Urban Borstnik
!> \date    2011-09-26
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2011-09-26 Split dbcsr_internal_operations
! *****************************************************************************
MODULE dbcsr_mm_hostdrv
  USE dbcsr_config,                    ONLY: has_acc,&
                                             mm_driver,&
                                             mm_driver_blas,&
                                             mm_driver_matmul,&
                                             mm_driver_smm,&
                                             mm_driver_xsmm
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width,&
                                             p_a_first,&
                                             p_b_first,&
                                             p_c_first,&
                                             p_k,&
                                             p_m,&
                                             p_n,&
                                             stack_descriptor_type
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_work_type
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp

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
# 47 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads


  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_hostdrv'

  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

  PUBLIC :: dbcsr_mm_hostdrv_lib_init, dbcsr_mm_hostdrv_lib_finalize
  PUBLIC :: dbcsr_mm_hostdrv_process
  PUBLIC :: dbcsr_mm_hostdrv_type
  PUBLIC :: dbcsr_mm_hostdrv_init

  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  TYPE dbcsr_mm_hostdrv_type
      TYPE(dbcsr_data_obj)          :: data_area
  END TYPE dbcsr_mm_hostdrv_type

CONTAINS

! *****************************************************************************
!> \brief Initialize the library
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_lib_init()







  END SUBROUTINE dbcsr_mm_hostdrv_lib_init


! *****************************************************************************
!> \brief Finalize the library
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_lib_finalize()







  END SUBROUTINE dbcsr_mm_hostdrv_lib_finalize


! *****************************************************************************
!> \brief Initialize the library
!> \param this ...
!> \param product_wm ...
! *****************************************************************************
 SUBROUTINE dbcsr_mm_hostdrv_init(this, product_wm)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_work_type), POINTER           :: product_wm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    this%data_area = product_wm%data_area

    CALL timestop(handle)

  END SUBROUTINE dbcsr_mm_hostdrv_init


! *****************************************************************************
!> \brief Calls the various drivers that process the stack.
!>
!> \param this ...
!> \param[in] left Left-matrix data 
!> \param[in] right Right-matrix data
!> \param[in] params           Stack of GEMM parameters
!> \param stack_size ...
!> \param stack_descr ...
!> \param success ...
!> \param used_smm ...
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_process(this, left, right, params, stack_size, &
       stack_descr, success, used_smm)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(1:dbcsr_ps_width, &
      stack_size), INTENT(INOUT)             :: params
    TYPE(stack_descriptor_type), INTENT(IN)  :: stack_descr
    LOGICAL, INTENT(OUT)                     :: success, used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_process', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.

    INTEGER                                  :: error_handle, sp
    REAL(KIND=dp)                            :: rnd

    IF(has_acc) & !for cpu-only runs this is called too often
       CALL timeset (routineN, error_handle)

    success = .TRUE. !host driver never failes...hopefully
    used_smm = .FALSE.

    IF (dbg) THEN
       CALL RANDOM_NUMBER (rnd)
       IF (rnd < 0.01_dp) THEN
          WRITE(*,*)routineN//" Stack size", stack_size, dbcsr_ps_width
          CALL print_gemm_parameters(params(:,1:stack_size))
       ENDIF
    ENDIF

    ! Verify stack consistency.  Only the upper bound is verified.
    IF (careful) THEN
       DO sp = 1, stack_size
          CALL dbcsr_assert (params(p_a_first,sp)&
               + params(p_m,sp) * params(p_k,sp) - 1,&
               "LE", dbcsr_data_get_size (left%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "A data out of bounds.", 180)
          CALL dbcsr_assert (params(p_b_first,sp)&
               + params(p_k,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (right%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "B data out of bounds.", 185)
          CALL dbcsr_assert (params(p_c_first,sp)&
               + params(p_m,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (this%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "C data out of bounds.", 190)
       ENDDO
    ENDIF

    SELECT CASE (mm_driver)
    CASE (mm_driver_matmul)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL internal_process_mm_stack_s (params, &
               stack_size, &
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp)
       CASE (dbcsr_type_real_8)
          CALL internal_process_mm_stack_d (params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp)
       CASE (dbcsr_type_complex_4)
          CALL internal_process_mm_stack_c (params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp)
       CASE (dbcsr_type_complex_8)
          CALL internal_process_mm_stack_z (params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp)
       CASE default
          CALL cp__b("dbcsr/mm/dbcsr_mm_hostdrv.F",214,"Invalid data type")
       END SELECT
    CASE (mm_driver_smm)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL smm_process_mm_stack_s (stack_descr, params, &
               stack_size, &
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp, used_smm)
       CASE (dbcsr_type_real_8)
          CALL smm_process_mm_stack_d (stack_descr, params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp, used_smm)
       CASE (dbcsr_type_complex_4)
          CALL smm_process_mm_stack_c (stack_descr, params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp, used_smm)
       CASE (dbcsr_type_complex_8)
          CALL smm_process_mm_stack_z (stack_descr, params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp, used_smm)
       CASE default
          CALL cp__b("dbcsr/mm/dbcsr_mm_hostdrv.F",235,"Invalid data type")
       END SELECT

    CASE (mm_driver_xsmm)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL xsmm_process_mm_stack_s (stack_descr, params, stack_size, &
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp, used_smm)
       CASE (dbcsr_type_real_8)
          CALL xsmm_process_mm_stack_d (stack_descr, params, stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp, used_smm)
       CASE (dbcsr_type_complex_4)
          CALL xsmm_process_mm_stack_c (stack_descr, params, stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp, used_smm)
       CASE (dbcsr_type_complex_8)
          CALL xsmm_process_mm_stack_z (stack_descr, params, stack_size, &
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp, used_smm)
       CASE default
          CALL cp__b("dbcsr/mm/dbcsr_mm_hostdrv.F",253,"Invalid data type")
       END SELECT

    CASE (mm_driver_blas)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL blas_process_mm_stack_s (params,&
               stack_size,&
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp)
       CASE (dbcsr_type_real_8)
          CALL blas_process_mm_stack_d (params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp)
       CASE (dbcsr_type_complex_4)
          CALL blas_process_mm_stack_c (params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp)
       CASE (dbcsr_type_complex_8)
          CALL blas_process_mm_stack_z (params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp)
       CASE default
          CALL cp__b("dbcsr/mm/dbcsr_mm_hostdrv.F",275,"Invalid data type")
       END SELECT
    CASE default
       CALL cp__b("dbcsr/mm/dbcsr_mm_hostdrv.F",278,"Invalid multiplication driver")
    END SELECT


    IF(has_acc) & !for cpu-only runs this is called too often
       CALL timestop(error_handle)


  END SUBROUTINE dbcsr_mm_hostdrv_process


! *****************************************************************************
!> \brief Helper-routine used by dbcsr_mm_hostdrv_process to print debug info.
!> \param params ...
! *****************************************************************************
  SUBROUTINE print_gemm_parameters(params)
    INTEGER, DIMENSION(:, :), INTENT(in)     :: params

    INTEGER                                  :: sp

    DO sp = 1, SIZE(params,2)
       WRITE(*,'(1X,A,1X,I7,":",3(1X,I4,","),".",3(1X,I12,","))')&
            "GEMM PARAMETERS",&
               sp,&
               params(p_m,sp),&
               params(p_k,sp),&
               params(p_n,sp),&
               params(p_a_first,sp),&
               params(p_b_first,sp),&
               params(p_c_first,sp)
    ENDDO
  END SUBROUTINE print_gemm_parameters



# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_d.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE blas_process_mm_stack_d(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL DGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            1.0_real_8,& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            1.0_real_8,& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_d

! *****************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE internal_process_mm_stack_d(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_d_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_d


! *****************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE smm_process_mm_stack_d(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

# 137 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_d.f90"
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_d(params, stack_size,a_data, b_data, c_data)


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
  END SUBROUTINE smm_process_mm_stack_d


! *****************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE xsmm_process_mm_stack_d(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

# 178 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_d.f90"

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    REAL(kind=real_8), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

# 288 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_d.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_d(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.


  END SUBROUTINE xsmm_process_mm_stack_d

! *****************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! *****************************************************************************
  PURE SUBROUTINE internal_mm_d_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    REAL(kind=real_8), INTENT(IN)                      :: B(K,N)
    REAL(kind=real_8), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_d_nn
# 313 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_z.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE blas_process_mm_stack_z(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL ZGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_8),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_8),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_z

! *****************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE internal_process_mm_stack_z(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_z_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_z


! *****************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE smm_process_mm_stack_z(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

# 137 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_z.f90"
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_z(params, stack_size,a_data, b_data, c_data)


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
  END SUBROUTINE smm_process_mm_stack_z


! *****************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE xsmm_process_mm_stack_z(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

# 178 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_z.f90"

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    COMPLEX(kind=real_8), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

# 288 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_z.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_z(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.


  END SUBROUTINE xsmm_process_mm_stack_z

! *****************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! *****************************************************************************
  PURE SUBROUTINE internal_mm_z_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_z_nn
# 314 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_s.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE blas_process_mm_stack_s(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL SGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            1.0_real_4,& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            1.0_real_4,& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_s

! *****************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE internal_process_mm_stack_s(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_s_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_s


! *****************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE smm_process_mm_stack_s(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

# 137 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_s.f90"
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_s(params, stack_size,a_data, b_data, c_data)


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
  END SUBROUTINE smm_process_mm_stack_s


! *****************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE xsmm_process_mm_stack_s(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

# 178 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_s.f90"

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    REAL(kind=real_4), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

# 288 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_s.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_s(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.


  END SUBROUTINE xsmm_process_mm_stack_s

! *****************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! *****************************************************************************
  PURE SUBROUTINE internal_mm_s_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(kind=real_4), INTENT(INOUT)                   :: C(M,N)
    REAL(kind=real_4), INTENT(IN)                      :: B(K,N)
    REAL(kind=real_4), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_s_nn
# 315 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_c.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE blas_process_mm_stack_c(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL CGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_4),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_4),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_c

! *****************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE internal_process_mm_stack_c(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_c_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_c


! *****************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! *****************************************************************************
  SUBROUTINE smm_process_mm_stack_c(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

# 137 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_c.f90"
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_c(params, stack_size,a_data, b_data, c_data)


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
  END SUBROUTINE smm_process_mm_stack_c


! *****************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE xsmm_process_mm_stack_c(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

# 178 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_c.f90"

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    COMPLEX(kind=real_4), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

# 288 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv_c.f90"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(stack_descr))==-1) EXIT ;  END DO ; ENDIF
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_c(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.


  END SUBROUTINE xsmm_process_mm_stack_c

! *****************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! *****************************************************************************
  PURE SUBROUTINE internal_mm_c_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_4), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_4), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_c_nn
# 316 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

END MODULE dbcsr_mm_hostdrv
