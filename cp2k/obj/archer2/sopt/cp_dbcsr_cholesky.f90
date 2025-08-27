# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Interface to (sca)lapack for the Cholesky based procedures
!> \author  VW
!> \date    2009-09-08
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-09-08
! *****************************************************************************
MODULE cp_dbcsr_cholesky
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_info,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             sp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/base/base_uses.f90" 1
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
# 32 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F" 2

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_cholesky'

  PUBLIC :: cp_dbcsr_cholesky_decompose, cp_dbcsr_cholesky_invert,&
            cp_dbcsr_cholesky_restore

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_decompose(matrix,n,para_env,blacs_env)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_cholesky_decompose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_n, &
                                                nfullcols_total, &
                                                nfullrows_total
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix




    CALL timeset(routineN,handle)


    NULLIFY(fm_matrix, fm_struct)
    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix")
    CALL cp_fm_struct_release(fm_struct)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix)

    my_n = MIN(fm_matrix%matrix_struct%nrow_global,&
         fm_matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       IF(.NOT.(n<=my_n))CALL cp__a("cp_dbcsr_cholesky.F",93)
       my_n=n
    END IF

    a => fm_matrix%local_data
    a_sp => fm_matrix%local_data_sp

# 110 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F"

    IF(fm_matrix%use_sp) THEN
       CALL spotrf('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotrf('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF



    IF(info/=0) &
      CALL cp__b("cp_dbcsr_cholesky.F",120,"Cholesky decomposition failed. Matrix ill conditioned ?")

    CALL copy_fm_to_dbcsr(fm_matrix, matrix)

    CALL cp_fm_release(fm_matrix)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_decompose

! *****************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \param upper_to_full ...
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_invert(matrix,n,para_env,blacs_env,upper_to_full)
    TYPE(cp_dbcsr_type)                           :: matrix
    INTEGER, INTENT(in), OPTIONAL             :: n
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(cp_blacs_env_type), POINTER          :: blacs_env
    LOGICAL, INTENT(IN)                       :: upper_to_full

    CHARACTER(len=*), PARAMETER :: routineN='dbcsr_cholesky_invert',&
         routineP=moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp
    INTEGER                                   :: info,handle
    INTEGER                                   :: my_n, nfullrows_total, nfullcols_total
    TYPE(cp_fm_type), POINTER                 :: fm_matrix, fm_matrix_tmp
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct




    CALL timeset(routineN,handle)


    NULLIFY(fm_matrix, fm_struct)
    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullrows_total,para_env=para_env)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix")
    CALL cp_fm_struct_release(fm_struct)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix)

    my_n = MIN(fm_matrix%matrix_struct%nrow_global,&
         fm_matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       IF(.NOT.(n<=my_n))CALL cp__a("cp_dbcsr_cholesky.F",177)
       my_n=n
    END IF

    a => fm_matrix%local_data
    a_sp => fm_matrix%local_data_sp

# 195 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F"

    IF(fm_matrix%use_sp) THEN
       CALL spotri('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotri('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF



    IF(.NOT.(info==0))CALL cp__a("cp_dbcsr_cholesky.F",204)

    IF(upper_to_full) THEN
       CALL cp_fm_create(fm_matrix_tmp,fm_matrix%matrix_struct,name="fm_matrix_tmp")
       CALL cp_fm_upper_to_full(fm_matrix, fm_matrix_tmp)
       CALL cp_fm_release(fm_matrix_tmp)
    ENDIF

    CALL copy_fm_to_dbcsr(fm_matrix, matrix)

    CALL cp_fm_release(fm_matrix)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_invert

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param neig ...
!> \param matrixb ...
!> \param matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
!> \param para_env ...
!> \param blacs_env ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_restore(matrix,neig,matrixb,matrixout,op,pos,transa,&
       para_env,blacs_env)
    TYPE(cp_dbcsr_type)                                :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                            :: neig
    CHARACTER ( LEN = * ), INTENT ( IN )           :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: transa
    TYPE(cp_para_env_type), POINTER                :: para_env
    TYPE(cp_blacs_env_type), POINTER               :: blacs_env

    CHARACTER(len=*), PARAMETER :: routineN='dbcsr_cholesky_restore',&
         routineP=moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,b,out
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp,b_sp,out_sp
    INTEGER                                   :: itype,handle
    INTEGER                                   :: n
    REAL(KIND = dp)                           :: alpha
    INTEGER                                   :: myprow, mypcol, nfullrows_total, &
         nfullcols_total
    TYPE(cp_blacs_env_type), POINTER          :: context
    CHARACTER                                 :: chol_pos,chol_transa
    TYPE(cp_fm_type), POINTER                 :: fm_matrix,fm_matrixb,fm_matrixout
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct





    CALL timeset(routineN,handle)


    NULLIFY(fm_matrix, fm_matrixb, fm_matrixout, fm_struct)

    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix")
    CALL cp_fm_struct_release(fm_struct)

    CALL cp_dbcsr_get_info(matrixb,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrixb,fm_struct,name="fm_matrixb")
    CALL cp_fm_struct_release(fm_struct)

    CALL cp_dbcsr_get_info(matrixout,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrixout,fm_struct,name="fm_matrixout")
    CALL cp_fm_struct_release(fm_struct)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix)
    CALL copy_dbcsr_to_fm(matrixb, fm_matrixb)
    !CALL copy_dbcsr_to_fm(matrixout, fm_matrixout)

    context => fm_matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = fm_matrix%matrix_struct%nrow_global
    itype = 1
    IF(op /= "SOLVE" .AND. op /= "MULTIPLY")&
       CALL cp__b("cp_dbcsr_cholesky.F",294,"wrong argument op")

    IF (PRESENT(pos)) THEN
       SELECT CASE(pos)
       CASE("LEFT")
         chol_pos='L'
       CASE("RIGHT")
         chol_pos='R'
       CASE DEFAULT
          CALL cp__b("cp_dbcsr_cholesky.F",303,"wrong argument pos")
       END SELECT
    ELSE
       chol_pos='L'
    ENDIF

    chol_transa='N'
    IF (PRESENT(transa)) chol_transa=transa

    IF((fm_matrix%use_sp.NEQV.fm_matrixb%use_sp).OR.(fm_matrix%use_sp.NEQV.fm_matrixout%use_sp))&
       CALL cp__b("cp_dbcsr_cholesky.F",313,"not the same precision")

    ! notice b is the cholesky guy
    a => fm_matrix%local_data
    b => fm_matrixb%local_data
    out => fm_matrixout%local_data
    a_sp => fm_matrix%local_data_sp
    b_sp => fm_matrixb%local_data_sp
    out_sp => fm_matrixout%local_data_sp

# 352 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_dbcsr_cholesky.F"

    alpha=1.0_dp
    IF(fm_matrix%use_sp) THEN
       CALL scopy(neig*n,a_sp(1,1),1,out_sp(1,1),1)
    ELSE
       CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    ENDIF
    IF (op.EQ."SOLVE") THEN
       IF(fm_matrix%use_sp) THEN
          CALL strsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),SIZE(b_sp,1),out_sp(1,1),n)
       ELSE
          CALL dtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),SIZE(b,1),out(1,1),n)
       ENDIF
    ELSE
       IF(fm_matrix%use_sp) THEN
          CALL strmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),n,out_sp(1,1),n)
       ELSE
          CALL dtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),n,out(1,1),n)
       ENDIF
    ENDIF



    CALL copy_fm_to_dbcsr(fm_matrixout, matrixout)

    CALL cp_fm_release(fm_matrix)
    CALL cp_fm_release(fm_matrixb)
    CALL cp_fm_release(fm_matrixout)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_restore

END MODULE cp_dbcsr_cholesky

