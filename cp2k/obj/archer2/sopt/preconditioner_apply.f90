# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/preconditioner_apply.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/preconditioner_apply.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner_apply
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_copy, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_type
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_s_inverse,&
                                             ot_precond_solver_direct,&
                                             ot_precond_solver_inv_chol,&
                                             ot_precond_solver_update
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type

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
# 36 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/preconditioner_apply.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

  PUBLIC :: apply_preconditioner_fm, apply_preconditioner_dbcsr

CONTAINS

! *****************************************************************************
!> \brief applies a previously created preconditioner to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
  
    SELECT CASE (preconditioner_env%in_use)
    CASE (0)
       CALL cp__b("preconditioner_apply.F",66,"No preconditioner in use")
    CASE (ot_precond_full_single)
       CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
    CASE (ot_precond_full_all)
       CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
    CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
       SELECT CASE (preconditioner_env%solver)
       CASE(ot_precond_solver_inv_chol,ot_precond_solver_update)
          CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
       CASE(ot_precond_solver_direct)
          CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
       CASE DEFAULT
          CALL cp__b("preconditioner_apply.F",78,"Solver not implemented")
       END SELECT
    CASE DEFAULT
       CALL cp__b("preconditioner_apply.F",81,"Unknown preconditioner")
    END SELECT
  
    CALL timestop(handle)
  
  END SUBROUTINE apply_preconditioner_fm

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_preconditioner_dbcsr(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
  
    SELECT CASE (preconditioner_env%in_use)
    CASE (0)
       CALL cp__b("preconditioner_apply.F",108,"No preconditioner in use")
    CASE (ot_precond_full_single)
       CALL apply_single(preconditioner_env, matrix_in, matrix_out)
    CASE (ot_precond_full_all)
       CALL apply_all(preconditioner_env, matrix_in, matrix_out)
    CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
       SELECT CASE (preconditioner_env%solver)
       CASE(ot_precond_solver_inv_chol,ot_precond_solver_update)
          CALL apply_single(preconditioner_env, matrix_in, matrix_out)
       CASE(ot_precond_solver_direct)
          CALL cp__b("preconditioner_apply.F",118,"Apply_full_direct not supported with ot")
          !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
       CASE DEFAULT
          CALL cp__b("preconditioner_apply.F",121,"Wrong solver")
       END SELECT
    CASE DEFAULT
       CALL cp__b("preconditioner_apply.F",124,"Wrong preconditioner")
    END SELECT
  
    CALL timestop(handle)
  
  END SUBROUTINE apply_preconditioner_dbcsr

! *****************************************************************************
!> \brief apply to full matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,handle)
 
    CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
    CALL cp_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                    matrix_in,0.0_dp,matrix_out)
    CALL timestop(handle)

  END SUBROUTINE apply_full_single

! *****************************************************************************
!> \brief apply to dbcsr matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_single(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF (.NOT.ASSOCIATED(preconditioner_env%dbcsr_matrix)) &
       CALL cp__b("preconditioner_apply.F",175,"NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
    CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%dbcsr_matrix,matrix_in,&
         0.0_dp,matrix_out)

    CALL timestop(handle)

  END SUBROUTINE apply_single

! *****************************************************************************
!> \brief preconditioner contains the factorization, application done by
!>        solving the linear system
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_direct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n
    TYPE(cp_fm_type), POINTER                :: work

    CALL timeset(routineN,handle)
  
    CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
    CALL cp_fm_create(work,matrix_in%matrix_struct,name="apply_full_single",&
                      use_sp=matrix_in%use_sp)
    CALL cp_fm_cholesky_restore(matrix_in,k,preconditioner_env%fm,work,&
         &                      "SOLVE",transa="T")
    CALL cp_fm_cholesky_restore(work,k,preconditioner_env%fm,matrix_out,&
         &                      "SOLVE",transa="N")
    CALL cp_fm_release(work)
  
    CALL timestop(handle)
  
  END SUBROUTINE apply_full_direct

! *****************************************************************************
!> \brief full all to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, n, &
                                                ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_type), POINTER                :: matrix_tmp

    CALL timeset(routineN,handle)
  
    CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  
    CALL cp_fm_create(matrix_tmp,matrix_in%matrix_struct,name="apply_full_all")
    CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                               row_indices=row_indices, col_indices=col_indices, local_data=local_data)
  
    !
    CALL cp_gemm('T','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                    matrix_in,0.0_dp,matrix_tmp)
  
    ! do the right scaling
    DO j=1,ncol_local
    DO i=1,nrow_local
       dum=1.0_dp/MAX(preconditioner_env%energy_gap, &
               preconditioner_env%full_evals(row_indices(i))-preconditioner_env%occ_evals(col_indices(j)))
       local_data(i,j)=local_data(i,j)*dum
    ENDDO
    ENDDO
  
    ! mult back
    CALL cp_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                    matrix_tmp,0.0_dp,matrix_out)
  
    CALL cp_fm_release(matrix_tmp)
  
    CALL timestop(handle)
  
  END SUBROUTINE apply_full_all

! *****************************************************************************
!> \brief full all to a dbcsr matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! *****************************************************************************
  SUBROUTINE apply_all(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_offset, col_size, &
                                                handle, i, j, row, &
                                                row_offset, row_size
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: DATA
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)
  
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_copy(matrix_tmp,matrix_in,name="apply_full_all")
    CALL cp_dbcsr_multiply('T','N',1.0_dp,preconditioner_env%dbcsr_matrix, &
                    matrix_in,0.0_dp,matrix_tmp)
    ! do the right scaling
    CALL cp_dbcsr_iterator_start(iter, matrix_tmp)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, DATA, &
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO j=1,col_size
       DO i=1,row_size
          dum=1.0_dp/MAX(preconditioner_env%energy_gap, &
               preconditioner_env%full_evals( row_offset+i-1 )&
               -preconditioner_env%occ_evals( col_offset+j-1 ))
          DATA(i,j)=DATA(i,j)*dum
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    ! mult back
    CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%dbcsr_matrix, &
                    matrix_tmp,0.0_dp,matrix_out)
    CALL cp_dbcsr_release(matrix_tmp)
    CALL timestop(handle)
  
  END SUBROUTINE apply_all

END MODULE preconditioner_apply
