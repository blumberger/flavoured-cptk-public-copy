# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief various cholesky decomposition related routines
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_cholesky
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_log_handling,                 ONLY: cp_to_string
  USE kinds,                           ONLY: dp,&
                                             sp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/../base/base_uses.f90" 1
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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_cholesky'

  PUBLIC :: cp_fm_cholesky_decompose, cp_fm_cholesky_invert,&
       cp_fm_cholesky_reduce, cp_fm_cholesky_restore

!***
CONTAINS

! *****************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param info_out ...
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_decompose(matrix,n,info_out)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    INTEGER, INTENT(out), OPTIONAL           :: info_out

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_decompose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_n
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp




    CALL timeset(routineN,handle)


    my_n = MIN(matrix%matrix_struct%nrow_global,&
         matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       IF(.NOT.(n<=my_n))CALL cp__a("fm/cp_fm_cholesky.F",65)
       my_n=n
    END IF

    a => matrix%local_data
    a_sp => matrix%local_data_sp

# 82 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"

    IF(matrix%use_sp) THEN
       CALL spotrf('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotrf('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF



    IF (PRESENT(info_out)) THEN
      info_out = info
    ELSE
      IF (info/=0) &
        CALL cp__b("fm/cp_fm_cholesky.F",95,"Cholesky failed: the matrix is not positive definite or ill-conditioned.")
    END IF

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_decompose

! *****************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_invert(matrix,n)
    TYPE(cp_fm_type), POINTER           :: matrix
    INTEGER, INTENT(in), OPTIONAL                :: n

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_invert',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp
    INTEGER                                   :: info,handle
    INTEGER                                   :: my_n




    CALL timeset(routineN,handle)


    my_n = MIN(matrix%matrix_struct%nrow_global,&
         matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       IF(.NOT.(n<=my_n))CALL cp__a("fm/cp_fm_cholesky.F",130)
       my_n=n
    END IF

    a => matrix%local_data
    a_sp => matrix%local_data_sp

# 148 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"

    IF(matrix%use_sp) THEN
       CALL spotri('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotri('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF



    IF(.NOT.(info==0))CALL cp__a("fm/cp_fm_cholesky.F",157)

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_invert

! *****************************************************************************
!> \brief reduce a matrix pencil A,B to normal form
!>      B has to be cholesky decomposed with  cp_fm_cholesky_decompose
!>      before calling this routine
!>      A,B -> inv(U^T)*A*inv(U),1
!>      (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!> \param matrix the symmetric matrix A
!> \param matrixb the cholesky decomposition of matrix B
!> \param itype ...
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_reduce(matrix,matrixb, itype)
   TYPE(cp_fm_type), POINTER           :: matrix, matrixb
   INTEGER, OPTIONAL                   :: itype

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_reduce',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a, b
    INTEGER                                   :: info, handle
    INTEGER                                   :: n, my_itype





    CALL timeset(routineN,handle)

    n = matrix%matrix_struct%nrow_global

    my_itype =1
    IF( PRESENT(itype) ) my_itype = itype

    a => matrix%local_data
    b => matrixb%local_data

# 212 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"

    CALL dsygst(my_itype,'U',n,a(1,1),n,b(1,1),n,info)



    IF(.NOT.(info==0))CALL cp__a("fm/cp_fm_cholesky.F",217)

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_reduce

!
! op can be "SOLVE" (out = U^-1 * in ) or "MULTIPLY"   (out = U * in )
! pos can be "LEFT" or "RIGHT" (U at the left or at the right)
!
! DEPRECATED, see cp_fm_basic_linalg:cp_fm_triangular_multiply
!
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param neig ...
!> \param matrixb ...
!> \param matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_restore(matrix,neig,matrixb,matrixout,op,pos,transa)
    TYPE(cp_fm_type), POINTER          :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                         :: neig
    CHARACTER ( LEN = * ), INTENT ( IN )        :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: transa

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_restore',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: a,b,out
    REAL(KIND = sp), DIMENSION(:,:), POINTER         :: a_sp,b_sp,out_sp
    INTEGER                                   :: itype, handle
    INTEGER                                   :: n
    REAL(KIND = dp)                           :: alpha
    INTEGER                                   :: myprow, mypcol
    TYPE(cp_blacs_env_type), POINTER          :: context
    CHARACTER                                 :: chol_pos,chol_transa





    CALL timeset(routineN,handle)

    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    n = matrix%matrix_struct%nrow_global
    itype = 1
    IF(op /= "SOLVE" .AND. op /= "MULTIPLY")&
       CALL cp__b("fm/cp_fm_cholesky.F",269,"wrong argument op")

    IF (PRESENT(pos)) THEN
       SELECT CASE(pos)
       CASE("LEFT")
         chol_pos='L'
       CASE("RIGHT")
         chol_pos='R'
       CASE DEFAULT
          CALL cp__b("fm/cp_fm_cholesky.F",278,"wrong argument pos")
       END SELECT
    ELSE
       chol_pos='L'
    ENDIF

    chol_transa='N'
    IF (PRESENT(transa)) chol_transa=transa

    IF((matrix%use_sp.NEQV.matrixb%use_sp).OR.(matrix%use_sp.NEQV.matrixout%use_sp))&
       CALL cp__b("fm/cp_fm_cholesky.F",288,"not the same precision")

    ! notice b is the cholesky guy
    a => matrix%local_data
    b => matrixb%local_data
    out => matrixout%local_data
    a_sp => matrix%local_data_sp
    b_sp => matrixb%local_data_sp
    out_sp => matrixout%local_data_sp

# 327 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_cholesky.F"

    alpha=1.0_dp
    IF(matrix%use_sp) THEN
       CALL scopy(neig*n,a_sp(1,1),1,out_sp(1,1),1)
    ELSE
       CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    ENDIF
    IF (op.EQ."SOLVE") THEN
       IF(matrix%use_sp) THEN
          CALL strsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),SIZE(b_sp,1),out_sp(1,1),n)
       ELSE
          CALL dtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),SIZE(b,1),out(1,1),n)
       ENDIF
    ELSE
       IF(matrix%use_sp) THEN
          CALL strmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),n,out_sp(1,1),n)
       ELSE
          CALL dtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),n,out(1,1),n)
       ENDIF
    ENDIF



    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_restore

END MODULE cp_fm_cholesky
