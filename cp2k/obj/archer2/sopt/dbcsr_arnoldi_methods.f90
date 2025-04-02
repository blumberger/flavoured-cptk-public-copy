# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief interface to the low level methods for arnoldi iteration
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi_methods
  USE dbcsr_arnoldi_types,             ONLY: &
       arnoldi_control, arnoldi_data_c, arnoldi_data_d, arnoldi_data_s, &
       arnoldi_data_z, dbcsr_arnoldi_data, get_control, get_data_c, &
       get_data_d, get_data_s, get_data_z, has_d_cmplx, has_d_real, &
       has_s_cmplx, has_s_real, m_x_v_vectors
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p
  USE dbcsr_geev,                      ONLY: dbcsr_general_local_diag,&
                                             dbcsr_symm_local_diag,&
                                             dbcsr_tridiag_local_diag
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_get_info,&
                                             dbcsr_scale
  USE dbcsr_types,                     ONLY: dbcsr_iterator,&
                                             dbcsr_obj,&
                                             dbcsr_obj_type_p
  USE dbcsr_vector_operations,         ONLY: dbcsr_matrix_colvec_multiply
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/../../base/base_uses.f90" 1
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
# 40 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi_methods'

  PUBLIC :: arnoldi_init, build_subspace, compute_evals, arnoldi_iram,&
            gev_arnoldi_init, gev_build_subspace, gev_update_data

  INTERFACE convert_matrix
     MODULE PROCEDURE convert_matrix_z_to_d, convert_matrix_s_to_c
     MODULE PROCEDURE convert_matrix_d_to_z, convert_matrix_c_to_s
     MODULE PROCEDURE convert_matrix_z_to_z, convert_matrix_c_to_c
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Interface for the routine calcualting the implicit restarts
!>        Currently all based on lapack
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE arnoldi_iram(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL arnoldi_iram_d(arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL arnoldi_iram_s(arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL arnoldi_iram_z(arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL arnoldi_iram_c(arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE arnoldi_iram

! *****************************************************************************
!> \brief Interface to compute the eigenvalues of a nonsymmetric matrix
!>        This is only the serial version
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE compute_evals(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL compute_evals_d(arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL compute_evals_s(arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL compute_evals_z(arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL compute_evals_c(arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE compute_evals

! *****************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        currently it can only setup a random vector but can be improved to 
!>        various types of restarts easily
!> \param matrix pointer to the matrices as described in main interface
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_data all data concerning the subspace
! *****************************************************************************
  SUBROUTINE arnoldi_init(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL arnoldi_init_d (matrix, vectors, arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL arnoldi_init_s (matrix, vectors, arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL arnoldi_init_z (matrix, vectors, arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL arnoldi_init_c (matrix, vectors, arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE arnoldi_init

! *****************************************************************************
!> \brief Interface for the initialization of the arnoldi subspace creation
!>        for the generalized eigenvalue problem
!> \param matrix pointer to the matrices as described in main interface
!> \param matrix_arnoldi ...
!> \param vectors work vectors for the matrix vector multiplications
!> \param arnoldi_data all data concerning the subspace
! *****************************************************************************
  SUBROUTINE gev_arnoldi_init(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix, matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_arnoldi_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL gev_arnoldi_init_d (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL gev_arnoldi_init_s (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL gev_arnoldi_init_z (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL gev_arnoldi_init_c (matrix, matrix_arnoldi, vectors, arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE gev_arnoldi_init

! *****************************************************************************
!> \brief here the iterations are performed and the krylov space is constructed
!> \param matrix see above
!> \param vectors see above
!> \param arnoldi_data see above
! *****************************************************************************
  SUBROUTINE build_subspace(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL build_subspace_d (matrix, vectors, arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL build_subspace_s (matrix, vectors, arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL build_subspace_z (matrix, vectors, arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL build_subspace_c (matrix, vectors, arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE build_subspace

! *****************************************************************************
!> \brief here the iterations are performed and the krylov space for the generalized
!>        eigenvalue probelm is created
!> \param matrix see above
!> \param vectors see above
!> \param arnoldi_data see above
! *****************************************************************************
  SUBROUTINE gev_build_subspace(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_build_subspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL gev_build_subspace_d (matrix, vectors, arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL gev_build_subspace_s (matrix, vectors, arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL gev_build_subspace_z (matrix, vectors, arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL gev_build_subspace_c (matrix, vectors, arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE gev_build_subspace

! *****************************************************************************
!> \brief in the generalized eigenvalue the matrix depende on the projection
!>        therefore the outer loop has to build a new set of matrices for the 
!>        inner loop
!> \param matrix see above
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data see above
! *****************************************************************************
  SUBROUTINE gev_update_data(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix, matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(has_d_real(arnoldi_data))  CALL gev_update_data_d (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL gev_update_data_s (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL gev_update_data_z (matrix, matrix_arnoldi, vectors, arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL gev_update_data_c (matrix, matrix_arnoldi, vectors, arnoldi_data)

    CALL timestop(handle)

  END SUBROUTINE gev_update_data

! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_z_to_d(m_out,m_in)
    REAL(real_8), DIMENSION(:, :)            :: m_out
    COMPLEX(real_8), DIMENSION(:, :)         :: m_in

     m_out(:,:)=REAL(m_in(:,:),KIND=real_8)
  END SUBROUTINE convert_matrix_z_to_d

! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_c_to_s(m_out,m_in)
    REAL(real_4), DIMENSION(:, :)            :: m_out
    COMPLEX(real_4), DIMENSION(:, :)         :: m_in

     m_out(:,:)=REAL(m_in,KIND=real_4)
  END SUBROUTINE convert_matrix_c_to_s

! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_d_to_z(m_out,m_in)
    COMPLEX(real_8), DIMENSION(:, :)         :: m_out
    REAL(real_8), DIMENSION(:, :)            :: m_in

     m_out(:,:)=CMPLX(m_in,0.0,KIND=real_8)
  END SUBROUTINE convert_matrix_d_to_z

! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_s_to_c(m_out,m_in)
    COMPLEX(real_4), DIMENSION(:, :)         :: m_out
    REAL(real_4), DIMENSION(:, :)            :: m_in

     m_out(:,:)=CMPLX(m_in,0.0,KIND=real_4)
  END SUBROUTINE convert_matrix_s_to_c

! I kno that one is stupid but like this it simplifies the template
! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_z_to_z(m_out,m_in)
    COMPLEX(real_8), DIMENSION(:, :)         :: m_out, m_in

     m_out(:,:)=m_in
  END SUBROUTINE convert_matrix_z_to_z

! *****************************************************************************
!> \brief ...
!> \param m_out ...
!> \param m_in ...
! *****************************************************************************
  SUBROUTINE convert_matrix_c_to_c(m_out,m_in)
    COMPLEX(real_4), DIMENSION(:, :)         :: m_out, m_in

     m_out(:,:)=m_in
  END SUBROUTINE convert_matrix_c_to_c


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods_low_d.f90" 1

! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE compute_evals_d(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals_d', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_8), DIMENSION(:, :), ALLOCATABLE   :: levec
    TYPE(arnoldi_data_d), POINTER            :: ar_data
    INTEGER                                  :: ndim
    TYPE(arnoldi_control), POINTER           :: control

    ar_data=>get_data_d(arnoldi_data)
    control=> get_control(arnoldi_data)
    ndim=control%current_step
    ALLOCATE(levec(ndim, ndim))

! Needs antoher interface as the calls to real and complex geev differ (sucks!)
! only perform the diagonalization on processors which hold data
    IF(control%generalized_ev)THEN
       CALL dbcsr_symm_local_diag('V',ar_data%Hessenberg(1:ndim, 1:ndim), ndim,&
                                  ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
    ELSE
       IF(control%symmetric)THEN
          CALL dbcsr_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       ELSE
          CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       END IF
    END IF

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_d

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_init_d (matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_8)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_8)                         :: norm 
    REAL(kind=real_8), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_d), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp
    
    ar_data=>get_data_d(arnoldi_data)

   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_8 ; w_vec = 0.0_real_8
    ar_data%Hessenberg=0.0_real_8

    IF(control%has_initial_vector)THEN
       ! after calling the set routine the initial vector is stored in f_vec
       CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
       ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL dlarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF

    CALL transfer_dbcsr_to_local_array_d(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
    CALL compute_norms_d(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_8)/rnorm)

    ! Everything prepared, initialize the Arnoldi iteration
    CALL transfer_dbcsr_to_local_array_d(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! This permits to compute the subspace of a matrix which is a product of multiple matrices
    DO i=1, SIZE(matrix)
       CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                        0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL dbcsr_copy(vectors%input_vec, vectors%result_vec)
    END DO

    CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, w_vec, nrow_local, control%local_comp)

    ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
    ar_data%Hessenberg(1, 1)=DOT_PRODUCT(v_vec, w_vec)
    CALL mp_sum(ar_data%Hessenberg(1, 1), pcol_group)
    ar_data%f_vec=w_vec-v_vec*ar_data%Hessenberg(1, 1)

    ar_data%local_history(:, 1)=v_vec(:)

    ! We did the first step in here so we should set the current step for the subspace generation accordingly
    control%current_step=1

    DEALLOCATE(v_vec, w_vec)

  END SUBROUTINE arnoldi_init_d

! *****************************************************************************
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementaion which scales subspace size^4
!>        by replacing the lapack calls with direct math the 
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_iram_d(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram_d', &
         routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_d), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER                     :: control
    COMPLEX(real_8), DIMENSION(:,:), ALLOCATABLE  :: tmp_mat, safe_mat, Q, tmp_mat1
    COMPLEX(real_8), DIMENSION(:), ALLOCATABLE    :: work, tau, work_measure
    INTEGER                                   :: msize, lwork, i, j, info, nwant
    REAL(kind=real_8)                          :: beta, sigma
    REAL(kind=real_8),DIMENSION(:,:),ALLOCATABLE :: Qdata


    ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
    ar_data=>get_data_d(arnoldi_data)
    control=>get_control(arnoldi_data)
    msize=control%current_step
    nwant=control%nval_out
    ALLOCATE(tmp_mat(msize,msize)); ALLOCATE(safe_mat(msize,msize))
    ALLOCATE(Q(msize,msize)); ALLOCATE(tmp_mat1(msize,msize))
    ALLOCATE(work_measure(1))
    ALLOCATE(tau(msize)); ALLOCATE(Qdata(msize,msize))
    !make Q identity
    Q=CMPLX(0.0, 0.0, real_8)
    DO i=1,msize
       Q(i,i)=CMPLX(1.0, 0.0, real_8)
    END DO

    ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
    CALL convert_matrix(tmp_mat,ar_data%Hessenberg(1:msize,1:msize))
    safe_mat(:,:)=tmp_mat(:,:)

    DO i=1,msize
       ! A bit a strange check but in the end we only want to shift the unwanted evals
       IF(ANY(control%selected_ind==i))CYCLE
       ! Here we shift the matrix by subtracting unwanted evals from the diagonal
       DO j=1,msize
          tmp_mat(j,j)=tmp_mat(j,j)-ar_data%evals(i)
       END DO
       ! Now we repair the damage by QR factorizing
       lwork=-1
       CALL zgeqrf(msize,msize,tmp_mat,msize,tau,work_measure,lwork,info)
       lwork=INT(work_measure(1))
       IF (ALLOCATED(work)) THEN
          IF (SIZE(work).LT.lwork) THEN
             DEALLOCATE(work)
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(work)) ALLOCATE(work(lwork))
       CALL zgeqrf(msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
       CALL zungqr(msize,msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! update Q=Q*Q_current
       tmp_mat1(:,:)=Q(:,:)
       CALL zgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         Q,msize)       
       ! Update H=(Q*)HQ
       CALL zgemm('C','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat,msize,safe_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         tmp_mat1,msize)
       CALL zgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         safe_mat,msize)

       ! this one is crucial for numerics not to accumulate noise in the subdiagonals
       DO j=1,msize
          safe_mat(j+2:msize,j)=CMPLX(0.0, 0.0, real_8)
       END DO
       tmp_mat(:,:)=safe_mat(:,:)
    END DO

    ! Now we can compute our restart quantities
    ar_data%Hessenberg=0.0_real_8
    CALL convert_matrix(ar_data%Hessenberg(1:msize,1:msize),safe_mat)
    CALL convert_matrix(Qdata,Q)
  
    beta=ar_data%Hessenberg(nwant+1,nwant); sigma=Qdata(msize,nwant)

    !update the residuum and the basis vectors
    IF(control%local_comp)THEN
       ar_data%f_vec=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,nwant+1))*beta+ar_data%f_vec(:)*sigma
       ar_data%local_history(:,1:nwant)=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,1:nwant))
    END IF
    ! Set the current step to nwant so the subspace build knows where to start
    control%current_step=nwant
    
    DEALLOCATE(tmp_mat,safe_mat,Q,Qdata,tmp_mat1,work,tau,work_measure)
    
  END SUBROUTINE arnoldi_iram_d

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE build_subspace_d(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors), TARGET              :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_d', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_8)                        :: rnorm
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_d), POINTER  :: ar_data
    REAL(kind=real_8)                         :: norm
    REAL(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    TYPE(dbcsr_obj), POINTER                 :: input_vec, result_vec, swap_vec

    ar_data=>get_data_d(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_8 ; w_vec = 0.0_real_8
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))

    DO j=control%current_step, control%max_iter-1

       ! compute the vector norm of the residuum, get it real valued as well (rnorm)
       CALL compute_norms_d(ar_data%f_vec, norm, rnorm, control%pcol_group)

       ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
       IF(control%myproc==0)control%converged=rnorm.LT.REAL(control%threshold, real_8)
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT

       ! transfer normalized residdum to history and its norm to the Hessenberg matrix
       IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
       v_vec(:)=ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j+1)=v_vec(:); ar_data%Hessenberg(j+1, j)=norm

       input_vec=>vectors%input_vec
       result_vec=>vectors%result_vec
       CALL transfer_local_array_to_dbcsr_d(input_vec, v_vec, nrow_local, control%local_comp)
 
       ! This permits to compute the subspace of a matrix which is a product of two matrices
       DO i=1, SIZE(matrix)
          CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, input_vec, result_vec, 1.0_real_8, &
                                            0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
          swap_vec=>input_vec
          input_vec=>result_vec
          result_vec=>swap_vec
       END DO

       CALL transfer_dbcsr_to_local_array_d(input_vec, w_vec, nrow_local, control%local_comp)

       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_d(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j+1, &
                               ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
       CALL DGKS_ortho_d(h_vec, ar_data%f_vec, s_vec, nrow_local, j+1, ar_data%local_history, &
                                   ar_data%local_history, control%local_comp, control%pcol_group)
       ! Finally we can put the projections into our Hessenberg matrix
       ar_data%Hessenberg(1:j+1, j+1)= h_vec(1:j+1)
       control%current_step=j+1
    END DO

    ! compute the vector norm of the final residuum and put it in to Hessenberg
    CALL compute_norms_d(ar_data%f_vec, norm, rnorm, control%pcol_group)
    ar_data%Hessenberg(control%current_step+1, control%current_step)=norm
    
    ! broadcast the Hessenberg matrix so we don't need to care later on
    CALL mp_bcast(ar_data%Hessenberg, 0, control%mp_group)

    DEALLOCATE(v_vec, w_vec, h_vec, s_vec)

  END SUBROUTINE  build_subspace_d

! *****************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_local_array_d(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_8), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=0.0_real_8)
    IF(is_local)array(1:n)=data_vec(1:n)

  END SUBROUTINE transfer_dbcsr_to_local_array_d

! *****************************************************************************
!> \brief The inverse routine transfering data back from an array to a dbcsr
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_d(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_8), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=0.0_real_8)
    IF(is_local)data_vec(1:n)=array(1:n)

! *****************************************************************************

  END SUBROUTINE transfer_local_array_to_dbcsr_d

! *****************************************************************************
!> \brief Gram-Schmidt in matrix vector form
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param w_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_d(h_vec, f_vec, s_vec, w_vec, nrow_local,&
                                            j, local_history, reorth_mat, local_comp, pcol_group)
    REAL(kind=real_8), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    REAL(kind=real_8), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho_d', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=0.0_real_8; f_vec=0.0_real_8; s_vec=0.0_real_8
    IF(local_comp)CALL dgemv('T', nrow_local, j, 1.0_real_8, local_history, &
                                      nrow_local, w_vec, 1, 0.0_real_8, h_vec, 1)
    CALL mp_sum(h_vec(1:j), pcol_group)

    IF(local_comp)CALL dgemv('N', nrow_local, j, 1.0_real_8, reorth_mat, &
                                      nrow_local, h_vec, 1, 0.0_real_8, f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

    CALL timestop(handle)

  END SUBROUTINE Gram_Schmidt_ortho_d

! *****************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE DGKS_ortho_d(h_vec, f_vec, s_vec, nrow_local, j, &
                                    local_history, reorth_mat, local_comp, pcol_group)
    REAL(kind=real_8), DIMENSION(:)      :: h_vec, s_vec, f_vec
    REAL(kind=real_8), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DGKS_ortho_d', &
         routineP = moduleN//':'//routineN

    LOGICAL                                          :: local_comp
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(local_comp)CALL dgemv('T', nrow_local, j, 1.0_real_8, local_history, &
                                       nrow_local, f_vec, 1, 0.0_real_8, s_vec, 1)
    CALL mp_sum(s_vec(1:j), pcol_group)
    IF(local_comp)CALL dgemv('N', nrow_local, j, -1.0_real_8, reorth_mat, &
                                       nrow_local, s_vec, 1, 1.0_real_8, f_vec, 1)
    h_vec(1:j)=h_vec(1:j)+s_vec(1:j)

    CALL timestop(handle)

  END SUBROUTINE DGKS_ortho_d

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE compute_norms_d(vec, norm, rnorm, pcol_group)
    REAL(kind=real_8), DIMENSION(:)           :: vec
    REAL(real_8)                        :: rnorm
    REAL(kind=real_8)                         :: norm
    INTEGER                                  :: pcol_group

    ! the norm is computed along the processor column
    norm=DOT_PRODUCT(vec, vec)
    CALL mp_sum(norm, pcol_group)
    rnorm=SQRT(REAL(norm, real_8))
    norm=rnorm

  END SUBROUTINE compute_norms_d

! *****************************************************************************
!> \brief Computes the intial guess for the solution of the generalized eigenvalue 
!>        using the arnoldi method
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_arnoldi_init_d (matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_8)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_8)                         :: norm, denom
    REAL(kind=real_8), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_d), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp

    ar_data=>get_data_d(arnoldi_data)
    
   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_8 ; w_vec = 0.0_real_8
    ar_data%Hessenberg=0.0_real_8

    IF(control%has_initial_vector)THEN
    ! after calling the set routine the initial vector is stored in f_vec
        CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL dlarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF   

    CALL transfer_dbcsr_to_local_array_d(vectors%input_vec, v_vec, nrow_local, control%local_comp)
   
    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_d(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_8)/rnorm)

    CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                      0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, w_vec, nrow_local, control%local_comp)
   
    ar_data%rho_scale=0.0_real_8
    ar_data%rho_scale=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(ar_data%rho_scale, pcol_group)

    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                      0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, w_vec, nrow_local, control%local_comp)
    
    denom=0.0_real_8
    denom=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(denom, pcol_group)
    IF(control%myproc==0) ar_data%rho_scale=ar_data%rho_scale/denom
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

    ! if the maximum ev is requested we need to optimize with -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_8, -ar_data%rho_scale)
   
    ar_data%x_vec=v_vec

  END SUBROUTINE gev_arnoldi_init_d

! *****************************************************************************
!> \brief builds the basis rothogonal wrt. teh metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and 
!>        matrix_vector products are very differently defined. Therefore it is 
!>        cleaner to put it in a seperate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_build_subspace_d(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_d', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ncol_local, nrow_local, pcol_group
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_d), POINTER            :: ar_data
    REAL(kind=real_8)                         :: norm
    REAL(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    REAL(kind=real_8), ALLOCATABLE, DIMENSION(:,:)    :: Zmat, CZmat , BZmat 

    ar_data=>get_data_d(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.
    pcol_group=control%pcol_group

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_8 ; w_vec = 0.0_real_8
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))
    ALLOCATE(Zmat(nrow_local,control%max_iter)); ALLOCATE(CZmat(nrow_local,control%max_iter))
    ALLOCATE(BZmat(nrow_local,control%max_iter))

    CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                      0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, BZmat(:,1), nrow_local, control%local_comp)
    
    norm=0.0_real_8 
    norm=DOT_PRODUCT(ar_data%x_vec,BZmat(:,1)) 
    CALL mp_sum(norm, pcol_group)
    IF(control%local_comp)THEN
       Zmat(:,1)=ar_data%x_vec/SQRT(norm);  BZmat(:,1)= BZmat(:,1)/SQRT(norm)
    END IF

    DO j=1,control%max_iter
       control%current_step=j
       CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, Zmat(:,j), nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                        0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, CZmat(:,j), nrow_local, control%local_comp)
       w_vec(:)=CZmat(:,j)                         
       
       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_d(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               BZmat, Zmat, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can becom a problem later on, there is probably a good check, but we don't perform it
       CALL DGKS_ortho_d(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                                    Zmat, control%local_comp, control%pcol_group)
    
       CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                        0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, v_vec, nrow_local, control%local_comp)
       norm=0.0_real_8
       norm=DOT_PRODUCT(ar_data%f_vec,v_vec)
       CALL mp_sum(norm, pcol_group)
      
       IF(control%myproc==0)control%converged=REAL(norm,real_8).LT.EPSILON(REAL(1.0,real_8))
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT 
       IF(j==control%max_iter-1)EXIT

       IF(control%local_comp)THEN
          Zmat(:,j+1)=ar_data%f_vec/SQRT(norm);  BZmat(:,j+1)= v_vec(:)/SQRT(norm)
       END IF
    END DO

! getting a bit more complicated here as the final matrix is again a product which has to be computed with the 
! ditributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
    ar_data%Hessenberg=0.0_real_8
    IF(control%local_comp)THEN
       ar_data%Hessenberg(1:control%current_step,1:control%current_step)=&
          MATMUL(TRANSPOSE(CZmat(:,1:control%current_step)),Zmat(:,1:control%current_step))
    END IF
    CALL mp_sum(ar_data%Hessenberg,control%mp_group)

    ar_data%local_history=Zmat
    ! broadcast the Hessenberg matrix so we don't need to care later on

    DEALLOCATE(v_vec); DEALLOCATE(w_vec); DEALLOCATE(s_vec); DEALLOCATE(h_vec); DEALLOCATE(CZmat);
    DEALLOCATE(Zmat); DEALLOCATE(BZmat)

  END SUBROUTINE gev_build_subspace_d

! *****************************************************************************
!> \brief Updates all data after an inner loop of the generalized ev arnoldi. 
!>        Updates rho and C=A-rho*B accordingly.
!>        As an update scheme is used for he ev, the output ev has to be replaced
!>        with the updated one.
!>        Furthermore a convergence check is performed. The mv product could
!>        be skiiped by making clever use of the precomputed data, However,
!>        it is most likely not worth the effort.
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_update_data_d(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data_d', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_d), POINTER            :: ar_data
    INTEGER                                  :: ncol_local, nrow_local, ind, i
    REAL(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: v_vec
    REAL(real_8)                        :: rnorm
    REAL(kind=real_8)                         :: norm
    COMPLEX(real_8)                     :: val 

    control=>get_control(arnoldi_data)

    ar_data=>get_data_d(arnoldi_data)

! compute the new shift, hack around the problem templating the conversion
    val=ar_data%evals(control%selected_ind(1))
    ar_data%rho_scale=ar_data%rho_scale+REAL(val,real_8)
! compute the new eigenvector / initial guess for the next arnoldi loop    
    ar_data%x_vec=0.0_real_8
    DO i=1,control%current_step
       val=ar_data%revec(i,control%selected_ind(1))
       ar_data%x_vec(:)=ar_data%x_vec(:)+ar_data%local_history(:,i)*REAL(val,real_8)
    END DO
!      ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
!                        ar_data%revec(1:control%current_step,control%selected_ind(1)))

! update the C-matrix (A-rho*B), if teh maximum value is requested we have to use -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_8, -ar_data%rho_scale)

! compute convergence and set the correct eigenvalue and eigenvector
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    IF (ncol_local>0) THEN
       ALLOCATE(v_vec(nrow_local))
       CALL compute_norms_d(ar_data%x_vec, norm, rnorm, control%pcol_group)
       v_vec(:)=ar_data%x_vec(:)/rnorm
    ENDIF
    CALL transfer_local_array_to_dbcsr_d(vectors%input_vec, v_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_8, &
                                            0.0_real_8, vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_d(vectors%result_vec, v_vec, nrow_local, control%local_comp)
    IF (ncol_local>0) THEN
       CALL compute_norms_d(v_vec, norm, rnorm, control%pcol_group)
       ! check convergence 
       control%converged=rnorm.LT.control%threshold
       DEALLOCATE(v_vec)
    ENDIF
    ! and broadcast the real eigenvalue    
    CALL mp_bcast(control%converged,0,control%mp_group)
    ind=control%selected_ind(1)
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
    ar_data%evals(ind)=ar_data%rho_scale


  END SUBROUTINE gev_update_data_d

# 321 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods_low_s.f90" 1

! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE compute_evals_s(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals_s', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_4), DIMENSION(:, :), ALLOCATABLE   :: levec
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    INTEGER                                  :: ndim
    TYPE(arnoldi_control), POINTER           :: control

    ar_data=>get_data_s(arnoldi_data)
    control=> get_control(arnoldi_data)
    ndim=control%current_step
    ALLOCATE(levec(ndim, ndim))

! Needs antoher interface as the calls to real and complex geev differ (sucks!)
! only perform the diagonalization on processors which hold data
    IF(control%generalized_ev)THEN
       CALL dbcsr_symm_local_diag('V',ar_data%Hessenberg(1:ndim, 1:ndim), ndim,&
                                  ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
    ELSE
       IF(control%symmetric)THEN
          CALL dbcsr_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       ELSE
          CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       END IF
    END IF

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_s

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_init_s (matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_4)                         :: norm 
    REAL(kind=real_4), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp
    
    ar_data=>get_data_s(arnoldi_data)

   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ar_data%Hessenberg=0.0_real_4

    IF(control%has_initial_vector)THEN
       ! after calling the set routine the initial vector is stored in f_vec
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
       ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL slarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF

    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
    CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm)

    ! Everything prepared, initialize the Arnoldi iteration
    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! This permits to compute the subspace of a matrix which is a product of multiple matrices
    DO i=1, SIZE(matrix)
       CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL dbcsr_copy(vectors%input_vec, vectors%result_vec)
    END DO

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)

    ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
    ar_data%Hessenberg(1, 1)=DOT_PRODUCT(v_vec, w_vec)
    CALL mp_sum(ar_data%Hessenberg(1, 1), pcol_group)
    ar_data%f_vec=w_vec-v_vec*ar_data%Hessenberg(1, 1)

    ar_data%local_history(:, 1)=v_vec(:)

    ! We did the first step in here so we should set the current step for the subspace generation accordingly
    control%current_step=1

    DEALLOCATE(v_vec, w_vec)

  END SUBROUTINE arnoldi_init_s

! *****************************************************************************
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementaion which scales subspace size^4
!>        by replacing the lapack calls with direct math the 
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_iram_s(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram_s', &
         routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER                     :: control
    COMPLEX(real_4), DIMENSION(:,:), ALLOCATABLE  :: tmp_mat, safe_mat, Q, tmp_mat1
    COMPLEX(real_4), DIMENSION(:), ALLOCATABLE    :: work, tau, work_measure
    INTEGER                                   :: msize, lwork, i, j, info, nwant
    REAL(kind=real_4)                          :: beta, sigma
    REAL(kind=real_4),DIMENSION(:,:),ALLOCATABLE :: Qdata


    ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    msize=control%current_step
    nwant=control%nval_out
    ALLOCATE(tmp_mat(msize,msize)); ALLOCATE(safe_mat(msize,msize))
    ALLOCATE(Q(msize,msize)); ALLOCATE(tmp_mat1(msize,msize))
    ALLOCATE(work_measure(1))
    ALLOCATE(tau(msize)); ALLOCATE(Qdata(msize,msize))
    !make Q identity
    Q=CMPLX(0.0, 0.0, real_4)
    DO i=1,msize
       Q(i,i)=CMPLX(1.0, 0.0, real_4)
    END DO

    ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
    CALL convert_matrix(tmp_mat,ar_data%Hessenberg(1:msize,1:msize))
    safe_mat(:,:)=tmp_mat(:,:)

    DO i=1,msize
       ! A bit a strange check but in the end we only want to shift the unwanted evals
       IF(ANY(control%selected_ind==i))CYCLE
       ! Here we shift the matrix by subtracting unwanted evals from the diagonal
       DO j=1,msize
          tmp_mat(j,j)=tmp_mat(j,j)-ar_data%evals(i)
       END DO
       ! Now we repair the damage by QR factorizing
       lwork=-1
       CALL cgeqrf(msize,msize,tmp_mat,msize,tau,work_measure,lwork,info)
       lwork=INT(work_measure(1))
       IF (ALLOCATED(work)) THEN
          IF (SIZE(work).LT.lwork) THEN
             DEALLOCATE(work)
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(work)) ALLOCATE(work(lwork))
       CALL cgeqrf(msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
       CALL cungqr(msize,msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! update Q=Q*Q_current
       tmp_mat1(:,:)=Q(:,:)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         Q,msize)       
       ! Update H=(Q*)HQ
       CALL cgemm('C','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat,msize,safe_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         tmp_mat1,msize)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         safe_mat,msize)

       ! this one is crucial for numerics not to accumulate noise in the subdiagonals
       DO j=1,msize
          safe_mat(j+2:msize,j)=CMPLX(0.0, 0.0, real_4)
       END DO
       tmp_mat(:,:)=safe_mat(:,:)
    END DO

    ! Now we can compute our restart quantities
    ar_data%Hessenberg=0.0_real_4
    CALL convert_matrix(ar_data%Hessenberg(1:msize,1:msize),safe_mat)
    CALL convert_matrix(Qdata,Q)
  
    beta=ar_data%Hessenberg(nwant+1,nwant); sigma=Qdata(msize,nwant)

    !update the residuum and the basis vectors
    IF(control%local_comp)THEN
       ar_data%f_vec=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,nwant+1))*beta+ar_data%f_vec(:)*sigma
       ar_data%local_history(:,1:nwant)=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,1:nwant))
    END IF
    ! Set the current step to nwant so the subspace build knows where to start
    control%current_step=nwant
    
    DEALLOCATE(tmp_mat,safe_mat,Q,Qdata,tmp_mat1,work,tau,work_measure)
    
  END SUBROUTINE arnoldi_iram_s

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE build_subspace_s(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors), TARGET              :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_s', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_4)                        :: rnorm
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER  :: ar_data
    REAL(kind=real_4)                         :: norm
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    TYPE(dbcsr_obj), POINTER                 :: input_vec, result_vec, swap_vec

    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))

    DO j=control%current_step, control%max_iter-1

       ! compute the vector norm of the residuum, get it real valued as well (rnorm)
       CALL compute_norms_s(ar_data%f_vec, norm, rnorm, control%pcol_group)

       ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
       IF(control%myproc==0)control%converged=rnorm.LT.REAL(control%threshold, real_4)
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT

       ! transfer normalized residdum to history and its norm to the Hessenberg matrix
       IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
       v_vec(:)=ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j+1)=v_vec(:); ar_data%Hessenberg(j+1, j)=norm

       input_vec=>vectors%input_vec
       result_vec=>vectors%result_vec
       CALL transfer_local_array_to_dbcsr_s(input_vec, v_vec, nrow_local, control%local_comp)
 
       ! This permits to compute the subspace of a matrix which is a product of two matrices
       DO i=1, SIZE(matrix)
          CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, input_vec, result_vec, 1.0_real_4, &
                                            0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
          swap_vec=>input_vec
          input_vec=>result_vec
          result_vec=>swap_vec
       END DO

       CALL transfer_dbcsr_to_local_array_s(input_vec, w_vec, nrow_local, control%local_comp)

       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_s(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j+1, &
                               ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
       CALL DGKS_ortho_s(h_vec, ar_data%f_vec, s_vec, nrow_local, j+1, ar_data%local_history, &
                                   ar_data%local_history, control%local_comp, control%pcol_group)
       ! Finally we can put the projections into our Hessenberg matrix
       ar_data%Hessenberg(1:j+1, j+1)= h_vec(1:j+1)
       control%current_step=j+1
    END DO

    ! compute the vector norm of the final residuum and put it in to Hessenberg
    CALL compute_norms_s(ar_data%f_vec, norm, rnorm, control%pcol_group)
    ar_data%Hessenberg(control%current_step+1, control%current_step)=norm
    
    ! broadcast the Hessenberg matrix so we don't need to care later on
    CALL mp_bcast(ar_data%Hessenberg, 0, control%mp_group)

    DEALLOCATE(v_vec, w_vec, h_vec, s_vec)

  END SUBROUTINE  build_subspace_s

! *****************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_local_array_s(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=0.0_real_4)
    IF(is_local)array(1:n)=data_vec(1:n)

  END SUBROUTINE transfer_dbcsr_to_local_array_s

! *****************************************************************************
!> \brief The inverse routine transfering data back from an array to a dbcsr
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_s(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    REAL(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=0.0_real_4)
    IF(is_local)data_vec(1:n)=array(1:n)

! *****************************************************************************

  END SUBROUTINE transfer_local_array_to_dbcsr_s

! *****************************************************************************
!> \brief Gram-Schmidt in matrix vector form
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param w_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_s(h_vec, f_vec, s_vec, w_vec, nrow_local,&
                                            j, local_history, reorth_mat, local_comp, pcol_group)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho_s', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=0.0_real_4; f_vec=0.0_real_4; s_vec=0.0_real_4
    IF(local_comp)CALL sgemv('T', nrow_local, j, 1.0_real_4, local_history, &
                                      nrow_local, w_vec, 1, 0.0_real_4, h_vec, 1)
    CALL mp_sum(h_vec(1:j), pcol_group)

    IF(local_comp)CALL sgemv('N', nrow_local, j, 1.0_real_4, reorth_mat, &
                                      nrow_local, h_vec, 1, 0.0_real_4, f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

    CALL timestop(handle)

  END SUBROUTINE Gram_Schmidt_ortho_s

! *****************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE DGKS_ortho_s(h_vec, f_vec, s_vec, nrow_local, j, &
                                    local_history, reorth_mat, local_comp, pcol_group)
    REAL(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec
    REAL(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DGKS_ortho_s', &
         routineP = moduleN//':'//routineN

    LOGICAL                                          :: local_comp
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(local_comp)CALL sgemv('T', nrow_local, j, 1.0_real_4, local_history, &
                                       nrow_local, f_vec, 1, 0.0_real_4, s_vec, 1)
    CALL mp_sum(s_vec(1:j), pcol_group)
    IF(local_comp)CALL sgemv('N', nrow_local, j, -1.0_real_4, reorth_mat, &
                                       nrow_local, s_vec, 1, 1.0_real_4, f_vec, 1)
    h_vec(1:j)=h_vec(1:j)+s_vec(1:j)

    CALL timestop(handle)

  END SUBROUTINE DGKS_ortho_s

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE compute_norms_s(vec, norm, rnorm, pcol_group)
    REAL(kind=real_4), DIMENSION(:)           :: vec
    REAL(real_4)                        :: rnorm
    REAL(kind=real_4)                         :: norm
    INTEGER                                  :: pcol_group

    ! the norm is computed along the processor column
    norm=DOT_PRODUCT(vec, vec)
    CALL mp_sum(norm, pcol_group)
    rnorm=SQRT(REAL(norm, real_4))
    norm=rnorm

  END SUBROUTINE compute_norms_s

! *****************************************************************************
!> \brief Computes the intial guess for the solution of the generalized eigenvalue 
!>        using the arnoldi method
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_arnoldi_init_s (matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    REAL(kind=real_4)                         :: norm, denom
    REAL(kind=real_4), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp

    ar_data=>get_data_s(arnoldi_data)
    
   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ar_data%Hessenberg=0.0_real_4

    IF(control%has_initial_vector)THEN
    ! after calling the set routine the initial vector is stored in f_vec
        CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL slarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF   

    CALL transfer_dbcsr_to_local_array_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)
   
    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm)

    CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)
   
    ar_data%rho_scale=0.0_real_4
    ar_data%rho_scale=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(ar_data%rho_scale, pcol_group)

    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, w_vec, nrow_local, control%local_comp)
    
    denom=0.0_real_4
    denom=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(denom, pcol_group)
    IF(control%myproc==0) ar_data%rho_scale=ar_data%rho_scale/denom
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

    ! if the maximum ev is requested we need to optimize with -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_4, -ar_data%rho_scale)
   
    ar_data%x_vec=v_vec

  END SUBROUTINE gev_arnoldi_init_s

! *****************************************************************************
!> \brief builds the basis rothogonal wrt. teh metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and 
!>        matrix_vector products are very differently defined. Therefore it is 
!>        cleaner to put it in a seperate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_build_subspace_s(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_s', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ncol_local, nrow_local, pcol_group
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    REAL(kind=real_4)                         :: norm
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:,:)    :: Zmat, CZmat , BZmat 

    ar_data=>get_data_s(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.
    pcol_group=control%pcol_group

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = 0.0_real_4 ; w_vec = 0.0_real_4
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))
    ALLOCATE(Zmat(nrow_local,control%max_iter)); ALLOCATE(CZmat(nrow_local,control%max_iter))
    ALLOCATE(BZmat(nrow_local,control%max_iter))

    CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                      0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, BZmat(:,1), nrow_local, control%local_comp)
    
    norm=0.0_real_4 
    norm=DOT_PRODUCT(ar_data%x_vec,BZmat(:,1)) 
    CALL mp_sum(norm, pcol_group)
    IF(control%local_comp)THEN
       Zmat(:,1)=ar_data%x_vec/SQRT(norm);  BZmat(:,1)= BZmat(:,1)/SQRT(norm)
    END IF

    DO j=1,control%max_iter
       control%current_step=j
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, Zmat(:,j), nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, CZmat(:,j), nrow_local, control%local_comp)
       w_vec(:)=CZmat(:,j)                         
       
       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_s(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               BZmat, Zmat, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can becom a problem later on, there is probably a good check, but we don't perform it
       CALL DGKS_ortho_s(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                                    Zmat, control%local_comp, control%pcol_group)
    
       CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                        0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, v_vec, nrow_local, control%local_comp)
       norm=0.0_real_4
       norm=DOT_PRODUCT(ar_data%f_vec,v_vec)
       CALL mp_sum(norm, pcol_group)
      
       IF(control%myproc==0)control%converged=REAL(norm,real_4).LT.EPSILON(REAL(1.0,real_4))
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT 
       IF(j==control%max_iter-1)EXIT

       IF(control%local_comp)THEN
          Zmat(:,j+1)=ar_data%f_vec/SQRT(norm);  BZmat(:,j+1)= v_vec(:)/SQRT(norm)
       END IF
    END DO

! getting a bit more complicated here as the final matrix is again a product which has to be computed with the 
! ditributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
    ar_data%Hessenberg=0.0_real_4
    IF(control%local_comp)THEN
       ar_data%Hessenberg(1:control%current_step,1:control%current_step)=&
          MATMUL(TRANSPOSE(CZmat(:,1:control%current_step)),Zmat(:,1:control%current_step))
    END IF
    CALL mp_sum(ar_data%Hessenberg,control%mp_group)

    ar_data%local_history=Zmat
    ! broadcast the Hessenberg matrix so we don't need to care later on

    DEALLOCATE(v_vec); DEALLOCATE(w_vec); DEALLOCATE(s_vec); DEALLOCATE(h_vec); DEALLOCATE(CZmat);
    DEALLOCATE(Zmat); DEALLOCATE(BZmat)

  END SUBROUTINE gev_build_subspace_s

! *****************************************************************************
!> \brief Updates all data after an inner loop of the generalized ev arnoldi. 
!>        Updates rho and C=A-rho*B accordingly.
!>        As an update scheme is used for he ev, the output ev has to be replaced
!>        with the updated one.
!>        Furthermore a convergence check is performed. The mv product could
!>        be skiiped by making clever use of the precomputed data, However,
!>        it is most likely not worth the effort.
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_update_data_s(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_s), POINTER            :: ar_data
    INTEGER                                  :: ncol_local, nrow_local, ind, i
    REAL(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: v_vec
    REAL(real_4)                        :: rnorm
    REAL(kind=real_4)                         :: norm
    COMPLEX(real_4)                     :: val 

    control=>get_control(arnoldi_data)

    ar_data=>get_data_s(arnoldi_data)

! compute the new shift, hack around the problem templating the conversion
    val=ar_data%evals(control%selected_ind(1))
    ar_data%rho_scale=ar_data%rho_scale+REAL(val,real_4)
! compute the new eigenvector / initial guess for the next arnoldi loop    
    ar_data%x_vec=0.0_real_4
    DO i=1,control%current_step
       val=ar_data%revec(i,control%selected_ind(1))
       ar_data%x_vec(:)=ar_data%x_vec(:)+ar_data%local_history(:,i)*REAL(val,real_4)
    END DO
!      ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
!                        ar_data%revec(1:control%current_step,control%selected_ind(1)))

! update the C-matrix (A-rho*B), if teh maximum value is requested we have to use -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, 1.0_real_4, -ar_data%rho_scale)

! compute convergence and set the correct eigenvalue and eigenvector
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    IF (ncol_local>0) THEN
       ALLOCATE(v_vec(nrow_local))
       CALL compute_norms_s(ar_data%x_vec, norm, rnorm, control%pcol_group)
       v_vec(:)=ar_data%x_vec(:)/rnorm
    ENDIF
    CALL transfer_local_array_to_dbcsr_s(vectors%input_vec, v_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, 1.0_real_4, &
                                            0.0_real_4, vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_s(vectors%result_vec, v_vec, nrow_local, control%local_comp)
    IF (ncol_local>0) THEN
       CALL compute_norms_s(v_vec, norm, rnorm, control%pcol_group)
       ! check convergence 
       control%converged=rnorm.LT.control%threshold
       DEALLOCATE(v_vec)
    ENDIF
    ! and broadcast the real eigenvalue    
    CALL mp_bcast(control%converged,0,control%mp_group)
    ind=control%selected_ind(1)
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
    ar_data%evals(ind)=ar_data%rho_scale


  END SUBROUTINE gev_update_data_s

# 322 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods_low_z.f90" 1

! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE compute_evals_z(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals_z', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_8), DIMENSION(:, :), ALLOCATABLE   :: levec
    TYPE(arnoldi_data_z), POINTER            :: ar_data
    INTEGER                                  :: ndim
    TYPE(arnoldi_control), POINTER           :: control

    ar_data=>get_data_z(arnoldi_data)
    control=> get_control(arnoldi_data)
    ndim=control%current_step
    ALLOCATE(levec(ndim, ndim))

! Needs antoher interface as the calls to real and complex geev differ (sucks!)
! only perform the diagonalization on processors which hold data
    IF(control%generalized_ev)THEN
       CALL dbcsr_symm_local_diag('V',ar_data%Hessenberg(1:ndim, 1:ndim), ndim,&
                                  ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
    ELSE
       IF(control%symmetric)THEN
          CALL dbcsr_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       ELSE
          CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       END IF
    END IF

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_z

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_init_z (matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_8)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    COMPLEX(kind=real_8)                         :: norm 
    COMPLEX(kind=real_8), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_z), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp
    
    ar_data=>get_data_z(arnoldi_data)

   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_8) ; w_vec = CMPLX(0.0, 0.0, real_8)
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_8)

    IF(control%has_initial_vector)THEN
       ! after calling the set routine the initial vector is stored in f_vec
       CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
       ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL zlarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF

    CALL transfer_dbcsr_to_local_array_z(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
    CALL compute_norms_z(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_8)/rnorm)

    ! Everything prepared, initialize the Arnoldi iteration
    CALL transfer_dbcsr_to_local_array_z(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! This permits to compute the subspace of a matrix which is a product of multiple matrices
    DO i=1, SIZE(matrix)
       CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                        CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL dbcsr_copy(vectors%input_vec, vectors%result_vec)
    END DO

    CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, w_vec, nrow_local, control%local_comp)

    ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
    ar_data%Hessenberg(1, 1)=DOT_PRODUCT(v_vec, w_vec)
    CALL mp_sum(ar_data%Hessenberg(1, 1), pcol_group)
    ar_data%f_vec=w_vec-v_vec*ar_data%Hessenberg(1, 1)

    ar_data%local_history(:, 1)=v_vec(:)

    ! We did the first step in here so we should set the current step for the subspace generation accordingly
    control%current_step=1

    DEALLOCATE(v_vec, w_vec)

  END SUBROUTINE arnoldi_init_z

! *****************************************************************************
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementaion which scales subspace size^4
!>        by replacing the lapack calls with direct math the 
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_iram_z(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram_z', &
         routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_z), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER                     :: control
    COMPLEX(real_8), DIMENSION(:,:), ALLOCATABLE  :: tmp_mat, safe_mat, Q, tmp_mat1
    COMPLEX(real_8), DIMENSION(:), ALLOCATABLE    :: work, tau, work_measure
    INTEGER                                   :: msize, lwork, i, j, info, nwant
    COMPLEX(kind=real_8)                          :: beta, sigma
    COMPLEX(kind=real_8),DIMENSION(:,:),ALLOCATABLE :: Qdata


    ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
    ar_data=>get_data_z(arnoldi_data)
    control=>get_control(arnoldi_data)
    msize=control%current_step
    nwant=control%nval_out
    ALLOCATE(tmp_mat(msize,msize)); ALLOCATE(safe_mat(msize,msize))
    ALLOCATE(Q(msize,msize)); ALLOCATE(tmp_mat1(msize,msize))
    ALLOCATE(work_measure(1))
    ALLOCATE(tau(msize)); ALLOCATE(Qdata(msize,msize))
    !make Q identity
    Q=CMPLX(0.0, 0.0, real_8)
    DO i=1,msize
       Q(i,i)=CMPLX(1.0, 0.0, real_8)
    END DO

    ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
    CALL convert_matrix(tmp_mat,ar_data%Hessenberg(1:msize,1:msize))
    safe_mat(:,:)=tmp_mat(:,:)

    DO i=1,msize
       ! A bit a strange check but in the end we only want to shift the unwanted evals
       IF(ANY(control%selected_ind==i))CYCLE
       ! Here we shift the matrix by subtracting unwanted evals from the diagonal
       DO j=1,msize
          tmp_mat(j,j)=tmp_mat(j,j)-ar_data%evals(i)
       END DO
       ! Now we repair the damage by QR factorizing
       lwork=-1
       CALL zgeqrf(msize,msize,tmp_mat,msize,tau,work_measure,lwork,info)
       lwork=INT(work_measure(1))
       IF (ALLOCATED(work)) THEN
          IF (SIZE(work).LT.lwork) THEN
             DEALLOCATE(work)
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(work)) ALLOCATE(work(lwork))
       CALL zgeqrf(msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
       CALL zungqr(msize,msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! update Q=Q*Q_current
       tmp_mat1(:,:)=Q(:,:)
       CALL zgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         Q,msize)       
       ! Update H=(Q*)HQ
       CALL zgemm('C','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat,msize,safe_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         tmp_mat1,msize)
       CALL zgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_8),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_8),&
                         safe_mat,msize)

       ! this one is crucial for numerics not to accumulate noise in the subdiagonals
       DO j=1,msize
          safe_mat(j+2:msize,j)=CMPLX(0.0, 0.0, real_8)
       END DO
       tmp_mat(:,:)=safe_mat(:,:)
    END DO

    ! Now we can compute our restart quantities
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_8)
    CALL convert_matrix(ar_data%Hessenberg(1:msize,1:msize),safe_mat)
    CALL convert_matrix(Qdata,Q)
  
    beta=ar_data%Hessenberg(nwant+1,nwant); sigma=Qdata(msize,nwant)

    !update the residuum and the basis vectors
    IF(control%local_comp)THEN
       ar_data%f_vec=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,nwant+1))*beta+ar_data%f_vec(:)*sigma
       ar_data%local_history(:,1:nwant)=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,1:nwant))
    END IF
    ! Set the current step to nwant so the subspace build knows where to start
    control%current_step=nwant
    
    DEALLOCATE(tmp_mat,safe_mat,Q,Qdata,tmp_mat1,work,tau,work_measure)
    
  END SUBROUTINE arnoldi_iram_z

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE build_subspace_z(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors), TARGET              :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_z', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_8)                        :: rnorm
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_z), POINTER  :: ar_data
    COMPLEX(kind=real_8)                         :: norm
    COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    TYPE(dbcsr_obj), POINTER                 :: input_vec, result_vec, swap_vec

    ar_data=>get_data_z(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_8) ; w_vec = CMPLX(0.0, 0.0, real_8)
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))

    DO j=control%current_step, control%max_iter-1

       ! compute the vector norm of the residuum, get it real valued as well (rnorm)
       CALL compute_norms_z(ar_data%f_vec, norm, rnorm, control%pcol_group)

       ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
       IF(control%myproc==0)control%converged=rnorm.LT.REAL(control%threshold, real_8)
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT

       ! transfer normalized residdum to history and its norm to the Hessenberg matrix
       IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
       v_vec(:)=ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j+1)=v_vec(:); ar_data%Hessenberg(j+1, j)=norm

       input_vec=>vectors%input_vec
       result_vec=>vectors%result_vec
       CALL transfer_local_array_to_dbcsr_z(input_vec, v_vec, nrow_local, control%local_comp)
 
       ! This permits to compute the subspace of a matrix which is a product of two matrices
       DO i=1, SIZE(matrix)
          CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, input_vec, result_vec, CMPLX(1.0, 0.0, real_8), &
                                            CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
          swap_vec=>input_vec
          input_vec=>result_vec
          result_vec=>swap_vec
       END DO

       CALL transfer_dbcsr_to_local_array_z(input_vec, w_vec, nrow_local, control%local_comp)

       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_z(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j+1, &
                               ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
       CALL DGKS_ortho_z(h_vec, ar_data%f_vec, s_vec, nrow_local, j+1, ar_data%local_history, &
                                   ar_data%local_history, control%local_comp, control%pcol_group)
       ! Finally we can put the projections into our Hessenberg matrix
       ar_data%Hessenberg(1:j+1, j+1)= h_vec(1:j+1)
       control%current_step=j+1
    END DO

    ! compute the vector norm of the final residuum and put it in to Hessenberg
    CALL compute_norms_z(ar_data%f_vec, norm, rnorm, control%pcol_group)
    ar_data%Hessenberg(control%current_step+1, control%current_step)=norm
    
    ! broadcast the Hessenberg matrix so we don't need to care later on
    CALL mp_bcast(ar_data%Hessenberg, 0, control%mp_group)

    DEALLOCATE(v_vec, w_vec, h_vec, s_vec)

  END SUBROUTINE  build_subspace_z

! *****************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_local_array_z(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    COMPLEX(kind=real_8), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_8))
    IF(is_local)array(1:n)=data_vec(1:n)

  END SUBROUTINE transfer_dbcsr_to_local_array_z

! *****************************************************************************
!> \brief The inverse routine transfering data back from an array to a dbcsr
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_z(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    COMPLEX(kind=real_8), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_8))
    IF(is_local)data_vec(1:n)=array(1:n)

! *****************************************************************************

  END SUBROUTINE transfer_local_array_to_dbcsr_z

! *****************************************************************************
!> \brief Gram-Schmidt in matrix vector form
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param w_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_z(h_vec, f_vec, s_vec, w_vec, nrow_local,&
                                            j, local_history, reorth_mat, local_comp, pcol_group)
    COMPLEX(kind=real_8), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    COMPLEX(kind=real_8), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho_z', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=CMPLX(0.0, 0.0, real_8); f_vec=CMPLX(0.0, 0.0, real_8); s_vec=CMPLX(0.0, 0.0, real_8)
    IF(local_comp)CALL zgemv('T', nrow_local, j, CMPLX(1.0, 0.0, real_8), local_history, &
                                      nrow_local, w_vec, 1, CMPLX(0.0, 0.0, real_8), h_vec, 1)
    CALL mp_sum(h_vec(1:j), pcol_group)

    IF(local_comp)CALL zgemv('N', nrow_local, j, CMPLX(1.0, 0.0, real_8), reorth_mat, &
                                      nrow_local, h_vec, 1, CMPLX(0.0, 0.0, real_8), f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

    CALL timestop(handle)

  END SUBROUTINE Gram_Schmidt_ortho_z

! *****************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE DGKS_ortho_z(h_vec, f_vec, s_vec, nrow_local, j, &
                                    local_history, reorth_mat, local_comp, pcol_group)
    COMPLEX(kind=real_8), DIMENSION(:)      :: h_vec, s_vec, f_vec
    COMPLEX(kind=real_8), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DGKS_ortho_z', &
         routineP = moduleN//':'//routineN

    LOGICAL                                          :: local_comp
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(local_comp)CALL zgemv('T', nrow_local, j, CMPLX(1.0, 0.0, real_8), local_history, &
                                       nrow_local, f_vec, 1, CMPLX(0.0, 0.0, real_8), s_vec, 1)
    CALL mp_sum(s_vec(1:j), pcol_group)
    IF(local_comp)CALL zgemv('N', nrow_local, j, CMPLX(-1.0, 0.0, real_8), reorth_mat, &
                                       nrow_local, s_vec, 1, CMPLX(1.0, 0.0, real_8), f_vec, 1)
    h_vec(1:j)=h_vec(1:j)+s_vec(1:j)

    CALL timestop(handle)

  END SUBROUTINE DGKS_ortho_z

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE compute_norms_z(vec, norm, rnorm, pcol_group)
    COMPLEX(kind=real_8), DIMENSION(:)           :: vec
    REAL(real_8)                        :: rnorm
    COMPLEX(kind=real_8)                         :: norm
    INTEGER                                  :: pcol_group

    ! the norm is computed along the processor column
    norm=DOT_PRODUCT(vec, vec)
    CALL mp_sum(norm, pcol_group)
    rnorm=SQRT(REAL(norm, real_8))
    norm=CMPLX(rnorm, 0.0, real_8)

  END SUBROUTINE compute_norms_z

! *****************************************************************************
!> \brief Computes the intial guess for the solution of the generalized eigenvalue 
!>        using the arnoldi method
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_arnoldi_init_z (matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_8)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    COMPLEX(kind=real_8)                         :: norm, denom
    COMPLEX(kind=real_8), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_z), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp

    ar_data=>get_data_z(arnoldi_data)
    
   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_8) ; w_vec = CMPLX(0.0, 0.0, real_8)
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_8)

    IF(control%has_initial_vector)THEN
    ! after calling the set routine the initial vector is stored in f_vec
        CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL zlarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF   

    CALL transfer_dbcsr_to_local_array_z(vectors%input_vec, v_vec, nrow_local, control%local_comp)
   
    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_z(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_8)/rnorm)

    CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                      CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, w_vec, nrow_local, control%local_comp)
   
    ar_data%rho_scale=CMPLX(0.0, 0.0, real_8)
    ar_data%rho_scale=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(ar_data%rho_scale, pcol_group)

    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                      CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, w_vec, nrow_local, control%local_comp)
    
    denom=CMPLX(0.0, 0.0, real_8)
    denom=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(denom, pcol_group)
    IF(control%myproc==0) ar_data%rho_scale=ar_data%rho_scale/denom
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

    ! if the maximum ev is requested we need to optimize with -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, CMPLX(1.0, 0.0, real_8), -ar_data%rho_scale)
   
    ar_data%x_vec=v_vec

  END SUBROUTINE gev_arnoldi_init_z

! *****************************************************************************
!> \brief builds the basis rothogonal wrt. teh metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and 
!>        matrix_vector products are very differently defined. Therefore it is 
!>        cleaner to put it in a seperate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_build_subspace_z(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_z', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ncol_local, nrow_local, pcol_group
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_z), POINTER            :: ar_data
    COMPLEX(kind=real_8)                         :: norm
    COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:,:)    :: Zmat, CZmat , BZmat 

    ar_data=>get_data_z(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.
    pcol_group=control%pcol_group

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_8) ; w_vec = CMPLX(0.0, 0.0, real_8)
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))
    ALLOCATE(Zmat(nrow_local,control%max_iter)); ALLOCATE(CZmat(nrow_local,control%max_iter))
    ALLOCATE(BZmat(nrow_local,control%max_iter))

    CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                      CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, BZmat(:,1), nrow_local, control%local_comp)
    
    norm=CMPLX(0.0, 0.0, real_8) 
    norm=DOT_PRODUCT(ar_data%x_vec,BZmat(:,1)) 
    CALL mp_sum(norm, pcol_group)
    IF(control%local_comp)THEN
       Zmat(:,1)=ar_data%x_vec/SQRT(norm);  BZmat(:,1)= BZmat(:,1)/SQRT(norm)
    END IF

    DO j=1,control%max_iter
       control%current_step=j
       CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, Zmat(:,j), nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                        CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, CZmat(:,j), nrow_local, control%local_comp)
       w_vec(:)=CZmat(:,j)                         
       
       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_z(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               BZmat, Zmat, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can becom a problem later on, there is probably a good check, but we don't perform it
       CALL DGKS_ortho_z(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                                    Zmat, control%local_comp, control%pcol_group)
    
       CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                        CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, v_vec, nrow_local, control%local_comp)
       norm=CMPLX(0.0, 0.0, real_8)
       norm=DOT_PRODUCT(ar_data%f_vec,v_vec)
       CALL mp_sum(norm, pcol_group)
      
       IF(control%myproc==0)control%converged=REAL(norm,real_8).LT.EPSILON(REAL(1.0,real_8))
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT 
       IF(j==control%max_iter-1)EXIT

       IF(control%local_comp)THEN
          Zmat(:,j+1)=ar_data%f_vec/SQRT(norm);  BZmat(:,j+1)= v_vec(:)/SQRT(norm)
       END IF
    END DO

! getting a bit more complicated here as the final matrix is again a product which has to be computed with the 
! ditributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_8)
    IF(control%local_comp)THEN
       ar_data%Hessenberg(1:control%current_step,1:control%current_step)=&
          MATMUL(TRANSPOSE(CZmat(:,1:control%current_step)),Zmat(:,1:control%current_step))
    END IF
    CALL mp_sum(ar_data%Hessenberg,control%mp_group)

    ar_data%local_history=Zmat
    ! broadcast the Hessenberg matrix so we don't need to care later on

    DEALLOCATE(v_vec); DEALLOCATE(w_vec); DEALLOCATE(s_vec); DEALLOCATE(h_vec); DEALLOCATE(CZmat);
    DEALLOCATE(Zmat); DEALLOCATE(BZmat)

  END SUBROUTINE gev_build_subspace_z

! *****************************************************************************
!> \brief Updates all data after an inner loop of the generalized ev arnoldi. 
!>        Updates rho and C=A-rho*B accordingly.
!>        As an update scheme is used for he ev, the output ev has to be replaced
!>        with the updated one.
!>        Furthermore a convergence check is performed. The mv product could
!>        be skiiped by making clever use of the precomputed data, However,
!>        it is most likely not worth the effort.
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_update_data_z(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data_z', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_z), POINTER            :: ar_data
    INTEGER                                  :: ncol_local, nrow_local, ind, i
    COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:)      :: v_vec
    REAL(real_8)                        :: rnorm
    COMPLEX(kind=real_8)                         :: norm
    COMPLEX(real_8)                     :: val 

    control=>get_control(arnoldi_data)

    ar_data=>get_data_z(arnoldi_data)

! compute the new shift, hack around the problem templating the conversion
    val=ar_data%evals(control%selected_ind(1))
    ar_data%rho_scale=ar_data%rho_scale+val
! compute the new eigenvector / initial guess for the next arnoldi loop    
    ar_data%x_vec=CMPLX(0.0, 0.0, real_8)
    DO i=1,control%current_step
       val=ar_data%revec(i,control%selected_ind(1))
       ar_data%x_vec(:)=ar_data%x_vec(:)+ar_data%local_history(:,i)*val
    END DO
!      ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
!                        ar_data%revec(1:control%current_step,control%selected_ind(1)))

! update the C-matrix (A-rho*B), if teh maximum value is requested we have to use -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, CMPLX(1.0, 0.0, real_8), -ar_data%rho_scale)

! compute convergence and set the correct eigenvalue and eigenvector
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    IF (ncol_local>0) THEN
       ALLOCATE(v_vec(nrow_local))
       CALL compute_norms_z(ar_data%x_vec, norm, rnorm, control%pcol_group)
       v_vec(:)=ar_data%x_vec(:)/rnorm
    ENDIF
    CALL transfer_local_array_to_dbcsr_z(vectors%input_vec, v_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_8), &
                                            CMPLX(0.0, 0.0, real_8), vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_z(vectors%result_vec, v_vec, nrow_local, control%local_comp)
    IF (ncol_local>0) THEN
       CALL compute_norms_z(v_vec, norm, rnorm, control%pcol_group)
       ! check convergence 
       control%converged=rnorm.LT.control%threshold
       DEALLOCATE(v_vec)
    ENDIF
    ! and broadcast the real eigenvalue    
    CALL mp_bcast(control%converged,0,control%mp_group)
    ind=control%selected_ind(1)
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
    ar_data%evals(ind)=ar_data%rho_scale


  END SUBROUTINE gev_update_data_z

# 323 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods_low_c.f90" 1

! *****************************************************************************
!> \brief Call the correct eigensolver, in the arnoldi method only the right
!>        eigenvectors are used. Lefts are created here but dumped immediatly 
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE compute_evals_c(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_evals_c', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_4), DIMENSION(:, :), ALLOCATABLE   :: levec
    TYPE(arnoldi_data_c), POINTER            :: ar_data
    INTEGER                                  :: ndim
    TYPE(arnoldi_control), POINTER           :: control

    ar_data=>get_data_c(arnoldi_data)
    control=> get_control(arnoldi_data)
    ndim=control%current_step
    ALLOCATE(levec(ndim, ndim))

! Needs antoher interface as the calls to real and complex geev differ (sucks!)
! only perform the diagonalization on processors which hold data
    IF(control%generalized_ev)THEN
       CALL dbcsr_symm_local_diag('V',ar_data%Hessenberg(1:ndim, 1:ndim), ndim,&
                                  ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim))
    ELSE
       IF(control%symmetric)THEN
          CALL dbcsr_tridiag_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       ELSE
          CALL dbcsr_general_local_diag('N', 'V', ar_data%Hessenberg(1:ndim, 1:ndim), ndim, &
                                         ar_data%evals(1:ndim), ar_data%revec(1:ndim, 1:ndim), levec)
       END IF
    END IF

    DEALLOCATE(levec)

  END SUBROUTINE compute_evals_c

! *****************************************************************************
!> \brief Initialization of the arnoldi vector. Here a random vector is used,
!>        might be generalized in the future 
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_init_c (matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    COMPLEX(kind=real_4)                         :: norm 
    COMPLEX(kind=real_4), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_c), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp
    
    ar_data=>get_data_c(arnoldi_data)

   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_4) ; w_vec = CMPLX(0.0, 0.0, real_4)
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_4)

    IF(control%has_initial_vector)THEN
       ! after calling the set routine the initial vector is stored in f_vec
       CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
       ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL clarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF

    CALL transfer_dbcsr_to_local_array_c(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! compute the vector norm of the random vectorm, get it real valued as well (rnorm)
    CALL compute_norms_c(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm)

    ! Everything prepared, initialize the Arnoldi iteration
    CALL transfer_dbcsr_to_local_array_c(vectors%input_vec, v_vec, nrow_local, control%local_comp)

    ! This permits to compute the subspace of a matrix which is a product of multiple matrices
    DO i=1, SIZE(matrix)
       CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                        CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL dbcsr_copy(vectors%input_vec, vectors%result_vec)
    END DO

    CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, w_vec, nrow_local, control%local_comp)

    ! Put the projection into the Hessenberg matrix, and make the vectors orthonormal
    ar_data%Hessenberg(1, 1)=DOT_PRODUCT(v_vec, w_vec)
    CALL mp_sum(ar_data%Hessenberg(1, 1), pcol_group)
    ar_data%f_vec=w_vec-v_vec*ar_data%Hessenberg(1, 1)

    ar_data%local_history(:, 1)=v_vec(:)

    ! We did the first step in here so we should set the current step for the subspace generation accordingly
    control%current_step=1

    DEALLOCATE(v_vec, w_vec)

  END SUBROUTINE arnoldi_init_c

! *****************************************************************************
!> \brief Alogorithm for the implicit restarts in the arnoldi method
!>        this is an early implementaion which scales subspace size^4
!>        by replacing the lapack calls with direct math the 
!>        QR and  gemms can be made linear and a N^2 sacling will be acchieved
!>        however this already sets the framework but should be used with care
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE arnoldi_iram_c(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_iram_c', &
         routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_c), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER                     :: control
    COMPLEX(real_4), DIMENSION(:,:), ALLOCATABLE  :: tmp_mat, safe_mat, Q, tmp_mat1
    COMPLEX(real_4), DIMENSION(:), ALLOCATABLE    :: work, tau, work_measure
    INTEGER                                   :: msize, lwork, i, j, info, nwant
    COMPLEX(kind=real_4)                          :: beta, sigma
    COMPLEX(kind=real_4),DIMENSION(:,:),ALLOCATABLE :: Qdata


    ! This is just a terribly inefficient implementation but I hope it is correct and might serve as a reference
    ar_data=>get_data_c(arnoldi_data)
    control=>get_control(arnoldi_data)
    msize=control%current_step
    nwant=control%nval_out
    ALLOCATE(tmp_mat(msize,msize)); ALLOCATE(safe_mat(msize,msize))
    ALLOCATE(Q(msize,msize)); ALLOCATE(tmp_mat1(msize,msize))
    ALLOCATE(work_measure(1))
    ALLOCATE(tau(msize)); ALLOCATE(Qdata(msize,msize))
    !make Q identity
    Q=CMPLX(0.0, 0.0, real_4)
    DO i=1,msize
       Q(i,i)=CMPLX(1.0, 0.0, real_4)
    END DO

    ! Looks a bit odd, but safe_mat will contain the result in the end, while tmpmat gets violated by lapack
    CALL convert_matrix(tmp_mat,ar_data%Hessenberg(1:msize,1:msize))
    safe_mat(:,:)=tmp_mat(:,:)

    DO i=1,msize
       ! A bit a strange check but in the end we only want to shift the unwanted evals
       IF(ANY(control%selected_ind==i))CYCLE
       ! Here we shift the matrix by subtracting unwanted evals from the diagonal
       DO j=1,msize
          tmp_mat(j,j)=tmp_mat(j,j)-ar_data%evals(i)
       END DO
       ! Now we repair the damage by QR factorizing
       lwork=-1
       CALL cgeqrf(msize,msize,tmp_mat,msize,tau,work_measure,lwork,info)
       lwork=INT(work_measure(1))
       IF (ALLOCATED(work)) THEN
          IF (SIZE(work).LT.lwork) THEN
             DEALLOCATE(work)
          ENDIF
       ENDIF
       IF (.NOT.ALLOCATED(work)) ALLOCATE(work(lwork))
       CALL cgeqrf(msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! Ask Lapack to reconstruct Q from its own way of storing data (tmpmat will contain Q)
       CALL cungqr(msize,msize,msize,tmp_mat,msize,tau,work,lwork,info)
       ! update Q=Q*Q_current
       tmp_mat1(:,:)=Q(:,:)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         Q,msize)       
       ! Update H=(Q*)HQ
       CALL cgemm('C','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat,msize,safe_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         tmp_mat1,msize)
       CALL cgemm('N','N',msize,msize,msize,CMPLX(1.0, 0.0, real_4),tmp_mat1,msize,tmp_mat,msize,CMPLX(0.0, 0.0, real_4),&
                         safe_mat,msize)

       ! this one is crucial for numerics not to accumulate noise in the subdiagonals
       DO j=1,msize
          safe_mat(j+2:msize,j)=CMPLX(0.0, 0.0, real_4)
       END DO
       tmp_mat(:,:)=safe_mat(:,:)
    END DO

    ! Now we can compute our restart quantities
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_4)
    CALL convert_matrix(ar_data%Hessenberg(1:msize,1:msize),safe_mat)
    CALL convert_matrix(Qdata,Q)
  
    beta=ar_data%Hessenberg(nwant+1,nwant); sigma=Qdata(msize,nwant)

    !update the residuum and the basis vectors
    IF(control%local_comp)THEN
       ar_data%f_vec=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,nwant+1))*beta+ar_data%f_vec(:)*sigma
       ar_data%local_history(:,1:nwant)=MATMUL(ar_data%local_history(:,1:msize),Qdata(1:msize,1:nwant))
    END IF
    ! Set the current step to nwant so the subspace build knows where to start
    control%current_step=nwant
    
    DEALLOCATE(tmp_mat,safe_mat,Q,Qdata,tmp_mat1,work,tau,work_measure)
    
  END SUBROUTINE arnoldi_iram_c

! *****************************************************************************
!> \brief Here we create the Krylov subspace and fill the Hessenberg matrix
!>        convergence check is only performed on subspace convergence
!>        Gram Schidt is used to orthonogonalize. 
!>        If this is numericall not sufficient a Daniel, Gragg, Kaufman and Steward
!>        correction is performed
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE build_subspace_c(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors), TARGET              :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_c', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncol_local, nrow_local
    REAL(real_4)                        :: rnorm
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_c), POINTER  :: ar_data
    COMPLEX(kind=real_4)                         :: norm
    COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    TYPE(dbcsr_obj), POINTER                 :: input_vec, result_vec, swap_vec

    ar_data=>get_data_c(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_4) ; w_vec = CMPLX(0.0, 0.0, real_4)
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))

    DO j=control%current_step, control%max_iter-1

       ! compute the vector norm of the residuum, get it real valued as well (rnorm)
       CALL compute_norms_c(ar_data%f_vec, norm, rnorm, control%pcol_group)

       ! check convergence and inform everybody about it, a bit annoying to talk to everybody because of that
       IF(control%myproc==0)control%converged=rnorm.LT.REAL(control%threshold, real_4)
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT

       ! transfer normalized residdum to history and its norm to the Hessenberg matrix
       IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
       v_vec(:)=ar_data%f_vec(:)/rnorm; ar_data%local_history(:, j+1)=v_vec(:); ar_data%Hessenberg(j+1, j)=norm

       input_vec=>vectors%input_vec
       result_vec=>vectors%result_vec
       CALL transfer_local_array_to_dbcsr_c(input_vec, v_vec, nrow_local, control%local_comp)
 
       ! This permits to compute the subspace of a matrix which is a product of two matrices
       DO i=1, SIZE(matrix)
          CALL dbcsr_matrix_colvec_multiply(matrix(i)%matrix, input_vec, result_vec, CMPLX(1.0, 0.0, real_4), &
                                            CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
          swap_vec=>input_vec
          input_vec=>result_vec
          result_vec=>swap_vec
       END DO

       CALL transfer_dbcsr_to_local_array_c(input_vec, w_vec, nrow_local, control%local_comp)

       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_c(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j+1, &
                               ar_data%local_history, ar_data%local_history, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can become a problem later on, there is probably a good check whether it's necessary, but we don't perform it
       CALL DGKS_ortho_c(h_vec, ar_data%f_vec, s_vec, nrow_local, j+1, ar_data%local_history, &
                                   ar_data%local_history, control%local_comp, control%pcol_group)
       ! Finally we can put the projections into our Hessenberg matrix
       ar_data%Hessenberg(1:j+1, j+1)= h_vec(1:j+1)
       control%current_step=j+1
    END DO

    ! compute the vector norm of the final residuum and put it in to Hessenberg
    CALL compute_norms_c(ar_data%f_vec, norm, rnorm, control%pcol_group)
    ar_data%Hessenberg(control%current_step+1, control%current_step)=norm
    
    ! broadcast the Hessenberg matrix so we don't need to care later on
    CALL mp_bcast(ar_data%Hessenberg, 0, control%mp_group)

    DEALLOCATE(v_vec, w_vec, h_vec, s_vec)

  END SUBROUTINE  build_subspace_c

! *****************************************************************************
!> \brief Helper routine to transfer the all data of a dbcsr matrix to a local array
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_local_array_c(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    COMPLEX(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_4))
    IF(is_local)array(1:n)=data_vec(1:n)

  END SUBROUTINE transfer_dbcsr_to_local_array_c

! *****************************************************************************
!> \brief The inverse routine transfering data back from an array to a dbcsr
!> \param vec ...
!> \param array ...
!> \param n ...
!> \param is_local ...
! *****************************************************************************
  SUBROUTINE transfer_local_array_to_dbcsr_c(vec, array, n, is_local)
    TYPE(dbcsr_obj)                          :: vec
    COMPLEX(kind=real_4), DIMENSION(:)           :: array
    INTEGER                                  :: n
    LOGICAL                                  :: is_local
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec

    data_vec => dbcsr_get_data_p (vec%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_4))
    IF(is_local)data_vec(1:n)=array(1:n)

! *****************************************************************************

  END SUBROUTINE transfer_local_array_to_dbcsr_c

! *****************************************************************************
!> \brief Gram-Schmidt in matrix vector form
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param w_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE Gram_Schmidt_ortho_c(h_vec, f_vec, s_vec, w_vec, nrow_local,&
                                            j, local_history, reorth_mat, local_comp, pcol_group)
    COMPLEX(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec, w_vec
    COMPLEX(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group
    LOGICAL                                          :: local_comp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'Gram_Schmidt_ortho_c', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    ! Let's do the orthonormalization, first try the Gram Schmidt scheme
    h_vec=CMPLX(0.0, 0.0, real_4); f_vec=CMPLX(0.0, 0.0, real_4); s_vec=CMPLX(0.0, 0.0, real_4)
    IF(local_comp)CALL cgemv('T', nrow_local, j, CMPLX(1.0, 0.0, real_4), local_history, &
                                      nrow_local, w_vec, 1, CMPLX(0.0, 0.0, real_4), h_vec, 1)
    CALL mp_sum(h_vec(1:j), pcol_group)

    IF(local_comp)CALL cgemv('N', nrow_local, j, CMPLX(1.0, 0.0, real_4), reorth_mat, &
                                      nrow_local, h_vec, 1, CMPLX(0.0, 0.0, real_4), f_vec, 1)
    f_vec(:)=w_vec(:)-f_vec(:)

    CALL timestop(handle)

  END SUBROUTINE Gram_Schmidt_ortho_c

! *****************************************************************************
!> \brief Compute the  Daniel, Gragg, Kaufman and Steward correction
!> \param h_vec ...
!> \param f_vec ...
!> \param s_vec ...
!> \param nrow_local ...
!> \param j ...
!> \param local_history ...
!> \param reorth_mat ...
!> \param local_comp ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE DGKS_ortho_c(h_vec, f_vec, s_vec, nrow_local, j, &
                                    local_history, reorth_mat, local_comp, pcol_group)
    COMPLEX(kind=real_4), DIMENSION(:)      :: h_vec, s_vec, f_vec
    COMPLEX(kind=real_4), DIMENSION(:, :)    :: local_history, reorth_mat
    INTEGER                                          :: nrow_local, j, pcol_group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DGKS_ortho_c', &
         routineP = moduleN//':'//routineN

    LOGICAL                                          :: local_comp
    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

    IF(local_comp)CALL cgemv('T', nrow_local, j, CMPLX(1.0, 0.0, real_4), local_history, &
                                       nrow_local, f_vec, 1, CMPLX(0.0, 0.0, real_4), s_vec, 1)
    CALL mp_sum(s_vec(1:j), pcol_group)
    IF(local_comp)CALL cgemv('N', nrow_local, j, CMPLX(-1.0, 0.0, real_4), reorth_mat, &
                                       nrow_local, s_vec, 1, CMPLX(1.0, 0.0, real_4), f_vec, 1)
    h_vec(1:j)=h_vec(1:j)+s_vec(1:j)

    CALL timestop(handle)

  END SUBROUTINE DGKS_ortho_c

! *****************************************************************************
!> \brief Compute the norm of a vector distributed along proc_col
!>        as local arrays. Always return the real part next to the complex rep.
!> \param vec ...
!> \param norm ...
!> \param rnorm ...
!> \param pcol_group ...
! *****************************************************************************
  SUBROUTINE compute_norms_c(vec, norm, rnorm, pcol_group)
    COMPLEX(kind=real_4), DIMENSION(:)           :: vec
    REAL(real_4)                        :: rnorm
    COMPLEX(kind=real_4)                         :: norm
    INTEGER                                  :: pcol_group

    ! the norm is computed along the processor column
    norm=DOT_PRODUCT(vec, vec)
    CALL mp_sum(norm, pcol_group)
    rnorm=SQRT(REAL(norm, real_4))
    norm=CMPLX(rnorm, 0.0, real_4)

  END SUBROUTINE compute_norms_c

! *****************************************************************************
!> \brief Computes the intial guess for the solution of the generalized eigenvalue 
!>        using the arnoldi method
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_arnoldi_init_c (matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'arnoldi_init_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iseed(4), row_size, col_size, &
                                                nrow_local, ncol_local, pcol_group, &
                                                row, col
    REAL(real_4)                        :: rnorm
    TYPE(dbcsr_iterator)                     :: iter
    COMPLEX(kind=real_4)                         :: norm, denom
    COMPLEX(kind=real_4), DIMENSION(:), ALLOCATABLE :: &
                                                v_vec, w_vec
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    LOGICAL                                  :: transposed, local_comp
    TYPE(arnoldi_data_c), POINTER            :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    pcol_group=control%pcol_group
    local_comp=control%local_comp

    ar_data=>get_data_c(arnoldi_data)
    
   ! create a local data copy to store the vectors and make Gram Schmidt a bit simpler
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local))
    ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_4) ; w_vec = CMPLX(0.0, 0.0, real_4)
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_4)

    IF(control%has_initial_vector)THEN
    ! after calling the set routine the initial vector is stored in f_vec
        CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
    ELSE
    ! Setup the initial normalized random vector (sufficient if it only happens on proc_col 0)
       CALL dbcsr_iterator_start(iter, vectors%input_vec)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, data_vec, transposed, row_size=row_size, col_size=col_size)
          iseed(1)=2; iseed(2)=MOD(row, 4095); iseed(3)=MOD(col, 4095); iseed(4)=11
          CALL clarnv(2, iseed, row_size*col_size, data_vec)
       END DO
       CALL dbcsr_iterator_stop(iter)
    END IF   

    CALL transfer_dbcsr_to_local_array_c(vectors%input_vec, v_vec, nrow_local, control%local_comp)
   
    ! compute the vector norm of the reandom vectorm, get it real valued as well (rnorm)
    CALL compute_norms_c(v_vec, norm, rnorm, control%pcol_group)

    IF (rnorm==0) rnorm=1 ! catch case where this rank has no actual data
    CALL dbcsr_scale(vectors%input_vec, REAL(1.0, real_4)/rnorm)

    CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                      CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, w_vec, nrow_local, control%local_comp)
   
    ar_data%rho_scale=CMPLX(0.0, 0.0, real_4)
    ar_data%rho_scale=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(ar_data%rho_scale, pcol_group)

    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                      CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)

    CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, w_vec, nrow_local, control%local_comp)
    
    denom=CMPLX(0.0, 0.0, real_4)
    denom=DOT_PRODUCT(v_vec,w_vec)
    CALL mp_sum(denom, pcol_group)
    IF(control%myproc==0) ar_data%rho_scale=ar_data%rho_scale/denom
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

    ! if the maximum ev is requested we need to optimize with -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, CMPLX(1.0, 0.0, real_4), -ar_data%rho_scale)
   
    ar_data%x_vec=v_vec

  END SUBROUTINE gev_arnoldi_init_c

! *****************************************************************************
!> \brief builds the basis rothogonal wrt. teh metric.
!>        The structure looks similar to normal arnoldi but norms, vectors and 
!>        matrix_vector products are very differently defined. Therefore it is 
!>        cleaner to put it in a seperate subroutine to avoid confusion
!> \param matrix ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_build_subspace_c(matrix, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_c', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: j, ncol_local, nrow_local, pcol_group
    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_c), POINTER            :: ar_data
    COMPLEX(kind=real_4)                         :: norm
    COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: h_vec, s_vec, v_vec, w_vec
    COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:,:)    :: Zmat, CZmat , BZmat 

    ar_data=>get_data_c(arnoldi_data)
    control=>get_control(arnoldi_data)
    control%converged=.FALSE.
    pcol_group=control%pcol_group

    ! create the vectors required during the iterations
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ALLOCATE(v_vec(nrow_local));  ALLOCATE(w_vec(nrow_local))
    v_vec = CMPLX(0.0, 0.0, real_4) ; w_vec = CMPLX(0.0, 0.0, real_4)
    ALLOCATE(s_vec(control%max_iter)); ALLOCATE(h_vec(control%max_iter))
    ALLOCATE(Zmat(nrow_local,control%max_iter)); ALLOCATE(CZmat(nrow_local,control%max_iter))
    ALLOCATE(BZmat(nrow_local,control%max_iter))

    CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, ar_data%x_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                      CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, BZmat(:,1), nrow_local, control%local_comp)
    
    norm=CMPLX(0.0, 0.0, real_4) 
    norm=DOT_PRODUCT(ar_data%x_vec,BZmat(:,1)) 
    CALL mp_sum(norm, pcol_group)
    IF(control%local_comp)THEN
       Zmat(:,1)=ar_data%x_vec/SQRT(norm);  BZmat(:,1)= BZmat(:,1)/SQRT(norm)
    END IF

    DO j=1,control%max_iter
       control%current_step=j
       CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, Zmat(:,j), nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                        CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, CZmat(:,j), nrow_local, control%local_comp)
       w_vec(:)=CZmat(:,j)                         
       
       ! Let's do the orthonormalization, to get the new f_vec. First try the Gram Schmidt scheme
       CALL Gram_Schmidt_ortho_c(h_vec, ar_data%f_vec, s_vec, w_vec, nrow_local, j, &
                               BZmat, Zmat, control%local_comp, control%pcol_group)

       ! A bit more expensive but simpliy always top up with a DGKS correction, otherwise numerics
       ! can becom a problem later on, there is probably a good check, but we don't perform it
       CALL DGKS_ortho_c(h_vec, ar_data%f_vec, s_vec, nrow_local, j, BZmat, &
                                    Zmat, control%local_comp, control%pcol_group)
    
       CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, ar_data%f_vec, nrow_local, control%local_comp)
       CALL dbcsr_matrix_colvec_multiply(matrix(2)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                        CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
       CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, v_vec, nrow_local, control%local_comp)
       norm=CMPLX(0.0, 0.0, real_4)
       norm=DOT_PRODUCT(ar_data%f_vec,v_vec)
       CALL mp_sum(norm, pcol_group)
      
       IF(control%myproc==0)control%converged=REAL(norm,real_4).LT.EPSILON(REAL(1.0,real_4))
       CALL mp_bcast(control%converged, 0, control%mp_group)
       IF(control%converged)EXIT 
       IF(j==control%max_iter-1)EXIT

       IF(control%local_comp)THEN
          Zmat(:,j+1)=ar_data%f_vec/SQRT(norm);  BZmat(:,j+1)= v_vec(:)/SQRT(norm)
       END IF
    END DO

! getting a bit more complicated here as the final matrix is again a product which has to be computed with the 
! ditributed vectors, therefore a sum along the first proc_col is necessary. As we want that matrix everywhere,
! we set it to zero before and compute the distributed product only on the first col and then sum over the full grid
    ar_data%Hessenberg=CMPLX(0.0, 0.0, real_4)
    IF(control%local_comp)THEN
       ar_data%Hessenberg(1:control%current_step,1:control%current_step)=&
          MATMUL(TRANSPOSE(CZmat(:,1:control%current_step)),Zmat(:,1:control%current_step))
    END IF
    CALL mp_sum(ar_data%Hessenberg,control%mp_group)

    ar_data%local_history=Zmat
    ! broadcast the Hessenberg matrix so we don't need to care later on

    DEALLOCATE(v_vec); DEALLOCATE(w_vec); DEALLOCATE(s_vec); DEALLOCATE(h_vec); DEALLOCATE(CZmat);
    DEALLOCATE(Zmat); DEALLOCATE(BZmat)

  END SUBROUTINE gev_build_subspace_c

! *****************************************************************************
!> \brief Updates all data after an inner loop of the generalized ev arnoldi. 
!>        Updates rho and C=A-rho*B accordingly.
!>        As an update scheme is used for he ev, the output ev has to be replaced
!>        with the updated one.
!>        Furthermore a convergence check is performed. The mv product could
!>        be skiiped by making clever use of the precomputed data, However,
!>        it is most likely not worth the effort.
!> \param matrix ...
!> \param matrix_arnoldi ...
!> \param vectors ...
!> \param arnoldi_data ...
! *****************************************************************************

  SUBROUTINE gev_update_data_c(matrix, matrix_arnoldi, vectors, arnoldi_data)
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix_arnoldi
    TYPE(m_x_v_vectors)                      :: vectors
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gev_update_data_c', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(arnoldi_data_c), POINTER            :: ar_data
    INTEGER                                  :: ncol_local, nrow_local, ind, i
    COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:)      :: v_vec
    REAL(real_4)                        :: rnorm
    COMPLEX(kind=real_4)                         :: norm
    COMPLEX(real_4)                     :: val 

    control=>get_control(arnoldi_data)

    ar_data=>get_data_c(arnoldi_data)

! compute the new shift, hack around the problem templating the conversion
    val=ar_data%evals(control%selected_ind(1))
    ar_data%rho_scale=ar_data%rho_scale+val
! compute the new eigenvector / initial guess for the next arnoldi loop    
    ar_data%x_vec=CMPLX(0.0, 0.0, real_4)
    DO i=1,control%current_step
       val=ar_data%revec(i,control%selected_ind(1))
       ar_data%x_vec(:)=ar_data%x_vec(:)+ar_data%local_history(:,i)*val
    END DO
!      ar_data%x_vec(:)=MATMUL(ar_data%local_history(:,1:control%current_step),&
!                        ar_data%revec(1:control%current_step,control%selected_ind(1)))

! update the C-matrix (A-rho*B), if teh maximum value is requested we have to use -A-rho*B
    CALL dbcsr_copy(matrix_arnoldi(1)%matrix,matrix(1)%matrix)
    CALL dbcsr_add(matrix_arnoldi(1)%matrix, matrix(2)%matrix, CMPLX(1.0, 0.0, real_4), -ar_data%rho_scale)

! compute convergence and set the correct eigenvalue and eigenvector
    CALL dbcsr_get_info(matrix=vectors%input_vec, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    IF (ncol_local>0) THEN
       ALLOCATE(v_vec(nrow_local))
       CALL compute_norms_c(ar_data%x_vec, norm, rnorm, control%pcol_group)
       v_vec(:)=ar_data%x_vec(:)/rnorm
    ENDIF
    CALL transfer_local_array_to_dbcsr_c(vectors%input_vec, v_vec, nrow_local, control%local_comp)
    CALL dbcsr_matrix_colvec_multiply(matrix_arnoldi(1)%matrix, vectors%input_vec, vectors%result_vec, CMPLX(1.0, 0.0, real_4), &
                                            CMPLX(0.0, 0.0, real_4), vectors%rep_row_vec, vectors%rep_col_vec)
    CALL transfer_dbcsr_to_local_array_c(vectors%result_vec, v_vec, nrow_local, control%local_comp)
    IF (ncol_local>0) THEN
       CALL compute_norms_c(v_vec, norm, rnorm, control%pcol_group)
       ! check convergence 
       control%converged=rnorm.LT.control%threshold
       DEALLOCATE(v_vec)
    ENDIF
    ! and broadcast the real eigenvalue    
    CALL mp_bcast(control%converged,0,control%mp_group)
    ind=control%selected_ind(1)
    CALL mp_bcast(ar_data%rho_scale,0,control%mp_group)

! Again the maximum value request is done on -A therefore the eigenvalue needs the opposite sign
    ar_data%evals(ind)=ar_data%rho_scale


  END SUBROUTINE gev_update_data_c

# 324 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_methods.F" 2

END MODULE dbcsr_arnoldi_methods 
