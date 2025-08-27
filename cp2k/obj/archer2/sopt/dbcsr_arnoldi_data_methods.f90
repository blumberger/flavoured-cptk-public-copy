# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The methods which allow to analyze and manipulate the arnoldi procedure
!>        The main routine and this should eb the only public access point for the 
!>        method
!> \par History
!>       2014.09 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_arnoldi_data_methods
  USE dbcsr_arnoldi_types,             ONLY: &
       arnoldi_control, arnoldi_data_c, arnoldi_data_d, arnoldi_data_s, &
       arnoldi_data_z, dbcsr_arnoldi_data, get_control, get_data_c, &
       get_data_d, get_data_s, get_data_z, get_evals_c, get_evals_d, &
       get_evals_s, get_evals_z, get_sel_ind, has_d_cmplx, has_d_real, &
       has_s_cmplx, has_s_real, set_control, set_data_c, set_data_d, &
       set_data_s, set_data_z
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p
  USE dbcsr_error_handling,            ONLY: dbcsr_abort
  USE dbcsr_methods,                   ONLY: dbcsr_get_matrix_type,&
                                             dbcsr_is_initialized,&
                                             dbcsr_release
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_grid_setup
  USE dbcsr_operations,                ONLY: dbcsr_get_info
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_obj,&
                                             dbcsr_obj_type_p,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_type_symmetric
  USE dbcsr_vector_operations,         ONLY: create_col_vec_from_matrix
  USE kinds,                           ONLY: real_4,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_arnoldi_data_methods'

  PUBLIC :: select_evals, get_selected_ritz_val,arnoldi_is_converged, &
            dbcsr_arnoldi_data, get_nrestart, set_arnoldi_initial_vector, &
            setup_arnoldi_data, deallocate_arnoldi_data, get_selected_ritz_vector

CONTAINS 

! *****************************************************************************
!> \brief This routine sets the environment for the arnoldi iteration and 
!>        the krylov subspace creation. All simulation parameters have to be given
!>        at this stage so the rest can run fully automated
!>        In addition, this routine allocates the data necessary for
!> \param arnoldi_data this type which gets filled with information and 
!>        on output contains all information necessary to extract
!>        whatever the user desires
!> \param matrix vector of matrices, only the first gets used to get some dimensions
!>        and parallel information needed later on
!> \param max_iter maximum dimension of the krylov subspace
!> \param threshold convergence threshold, this is used for both subspace and eigenval 
!> \param selection_crit integer defining according to which criterion the 
!>        eigenvalues are selected for the subspace
!> \param nval_request for some sel_crit useful, how many eV to select
!> \param nrestarts ...
!> \param generalized_ev ...
!> \param iram ...
! *****************************************************************************
  SUBROUTINE setup_arnoldi_data(arnoldi_data, matrix, max_iter, threshold, selection_crit,&
                                nval_request, nrestarts, generalized_ev,iram)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter
    REAL(real_8)                             :: threshold
    INTEGER                                  :: selection_crit, nval_request, &
                                                nrestarts
    LOGICAL                                  :: generalized_ev, iram

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_arnoldi_data', &
      routineP = moduleN//':'//routineN

    CALL setup_arnoldi_control(arnoldi_data, matrix, max_iter, threshold, selection_crit,&
                               nval_request, nrestarts, generalized_ev, iram)

    SELECT CASE (matrix(1)%matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL setup_arnoldi_data_d (arnoldi_data, matrix, max_iter)
    CASE (dbcsr_type_real_4)
       CALL setup_arnoldi_data_s (arnoldi_data, matrix, max_iter)
    CASE (dbcsr_type_complex_8)
       CALL setup_arnoldi_data_z (arnoldi_data, matrix, max_iter)
    CASE (dbcsr_type_complex_4)
       CALL setup_arnoldi_data_c (arnoldi_data, matrix, max_iter)
    END SELECT

  END SUBROUTINE setup_arnoldi_data

! *****************************************************************************
!> \brief Creates the control type for arnoldi, see above for details
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
!> \param threshold ...
!> \param selection_crit ...
!> \param nval_request ...
!> \param nrestarts ...
!> \param generalized_ev ...
!> \param iram ...
! *****************************************************************************
  SUBROUTINE setup_arnoldi_control(arnoldi_data, matrix, max_iter, threshold, selection_crit, nval_request,&
                                   nrestarts, generalized_ev, iram)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter
    REAL(real_8)                             :: threshold
    INTEGER                                  :: selection_crit, nval_request, &
                                                nrestarts
    LOGICAL                                  :: generalized_ev, iram

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_arnoldi_control', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control
    TYPE(dbcsr_distribution_obj)             :: distri

    ALLOCATE(control)
! Fill the information which will later on control the behavior of the arnoldi method and allow synchronization
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, distribution=distri)
    control%mp_group=distri%d%mp_env%mp%mp_group
    control%myproc=distri%d%mp_env%mp%mynode
    CALL dbcsr_mp_grid_setup(distri%d%mp_env)
    IF(.NOT. distri%d%mp_env%mp%subgroups_defined)&
       CALL dbcsr_abort(routineP,138,"arnoldi only with subgroups")
    control%pcol_group=distri%d%mp_env%mp%pcol_group
    control%prow_group=distri%d%mp_env%mp%prow_group

    control%symmetric=.FALSE.
! Will need a fix for complex because there it has to be hermitian
    IF(SIZE(matrix)==1)control%symmetric=&
                       dbcsr_get_matrix_type(matrix(1)%matrix)==dbcsr_type_symmetric

! Set the control parameters
    control%max_iter=max_iter
    control%current_step=0
    control%selection_crit=selection_crit
    control%nval_req=nval_request
    control%threshold=threshold
    control%converged=.FALSE.
    control%has_initial_vector=.FALSE.
    control%iram=iram
    control%nrestart=nrestarts
    control%generalized_ev=generalized_ev

    IF(control%nval_req>1.AND.control%nrestart>0.AND..NOT.control%iram)&
      CALL dbcsr_abort(routineP,160,'with more than one eigenvalue requested '//&
      'internal restarting with a previous EVEC is a bad idea, set IRAM or nrestsart=0')


! some checks for the generalized EV mode
    IF(control%generalized_ev.AND.selection_crit==1)&
       CALL dbcsr_abort(routineP,166,'generalized ev can only highest OR lowest EV')
    IF(control%generalized_ev.AND.nval_request.NE.1)&
       CALL dbcsr_abort(routineP,168,'generalized ev can only compute one EV at the time')
    IF(control%generalized_ev.AND.control%nrestart==0)&
       CALL dbcsr_abort(routineP,170,'outer loops are mandatory for generalized EV, set nrestart appropriatly')
    IF(SIZE(matrix).NE.2.AND.control%generalized_ev)&
       CALL dbcsr_abort(routineP,172,'generalized ev needs exactly two matrices as input (2nd is the metric)')

    ALLOCATE(control%selected_ind(max_iter))
    CALL set_control(arnoldi_data,control)

  END SUBROUTINE setup_arnoldi_control

! *****************************************************************************
!> \brief Deallocate the data in dbcsr_arnoldi_data
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector(arnoldi_data,ind,matrix,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix, vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector', &
      routineP = moduleN//':'//routineN

    IF(has_d_real(arnoldi_data))  CALL get_selected_ritz_vector_d(arnoldi_data,ind,matrix,vector)
    IF(has_s_real(arnoldi_data))  CALL get_selected_ritz_vector_s(arnoldi_data,ind,matrix,vector)
    IF(has_d_cmplx(arnoldi_data)) CALL get_selected_ritz_vector_z(arnoldi_data,ind,matrix,vector)
    IF(has_s_cmplx(arnoldi_data)) CALL get_selected_ritz_vector_c(arnoldi_data,ind,matrix,vector)

  END SUBROUTINE get_selected_ritz_vector

! *****************************************************************************
!> \brief Deallocate the data in dbcsr_arnoldi_data
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_control), POINTER           :: control

    IF(has_d_real(arnoldi_data))  CALL deallocate_arnoldi_data_d (arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL deallocate_arnoldi_data_s (arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL deallocate_arnoldi_data_z (arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL deallocate_arnoldi_data_c (arnoldi_data)

    control=>get_control(arnoldi_data)
    DEALLOCATE(control%selected_ind)
    DEALLOCATE(control)

  END SUBROUTINE deallocate_arnoldi_data

! *****************************************************************************
!> \brief perform the selection of eigenvalues, fills the selected_ind array
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE select_evals(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_evals', &
      routineP = moduleN//':'//routineN

    IF(has_d_real(arnoldi_data))  CALL select_evals_d(arnoldi_data)
    IF(has_s_real(arnoldi_data))  CALL select_evals_s(arnoldi_data)
    IF(has_d_cmplx(arnoldi_data)) CALL select_evals_z(arnoldi_data)
    IF(has_s_cmplx(arnoldi_data)) CALL select_evals_c(arnoldi_data)

  END SUBROUTINE select_evals

! *****************************************************************************
!> \brief set a new selection type, if you notice you didn't like the initial one
!> \param ar_data ...
!> \param itype ...
! *****************************************************************************
  SUBROUTINE set_eval_selection(ar_data,itype)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    INTEGER                                  :: itype

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(ar_data)
    control%selection_crit=itype
  END SUBROUTINE set_eval_selection

! *****************************************************************************
!> \brief returns the number of restarts allowed for arnoldi
!> \param arnoldi_data ...
!> \retval nrestart ...
! *****************************************************************************
  FUNCTION get_nrestart(arnoldi_data)         RESULT(nrestart)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: nrestart

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    nrestart=control%nrestart

  END FUNCTION get_nrestart

! *****************************************************************************
!> \brief get the number of eigenvalues matching the search criterion
!> \param ar_data ...
!> \retval nval_out ...
! *****************************************************************************
  FUNCTION get_nval_out(ar_data)     RESULT(nval_out)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    INTEGER                                  :: nval_out

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(ar_data)
    nval_out=control%nval_out

  END FUNCTION get_nval_out

! *****************************************************************************
!> \brief get the dimension of the krylov space. Can be less than max_iter 
!>        if subspace converged early
!> \param ar_data ...
!> \retval current_step ...
! *****************************************************************************
  FUNCTION get_subsp_size(ar_data)     RESULT(current_step)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    INTEGER                                  :: current_step

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(ar_data)
    current_step=control%current_step

  END FUNCTION get_subsp_size

! *****************************************************************************
!> \brief Find out whether the method with the current search criterion is converged
!> \param ar_data ...
!> \retval converged ...
! *****************************************************************************
  FUNCTION arnoldi_is_converged(ar_data)      RESULT(converged)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    LOGICAL                                  :: converged

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(ar_data)
    converged=control%converged

  END FUNCTION

! *****************************************************************************
!> \brief get a single specific Ritz value from the set of selected
!> \param ar_data ...
!> \param ind ...
!> \retval eval_out ...
! *****************************************************************************
  FUNCTION get_selected_ritz_val(ar_data, ind) RESULT(eval_out)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    INTEGER                                  :: ind
    COMPLEX(real_8)                          :: eval_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_val', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_4), DIMENSION(:), POINTER   :: evals_s
    COMPLEX(real_8), DIMENSION(:), POINTER   :: evals_d
    INTEGER                                  :: ev_ind
    INTEGER, DIMENSION(:), POINTER           :: selected_ind

    IF(ind.GT.get_nval_out(ar_data)) CALL dbcsr_abort(routineP,340,'outside range of indexed evals')

    selected_ind=>get_sel_ind(ar_data)
    ev_ind=selected_ind(ind)
    IF(has_d_real(ar_data)) THEN
       evals_d=>get_evals_d(ar_data);   eval_out=evals_d(ev_ind)     
    ELSE IF(has_s_real(ar_data)) THEN
       evals_s=>get_evals_s(ar_data);   eval_out=CMPLX(evals_s(ev_ind),KIND=real_8)     
    ELSE IF(has_d_cmplx(ar_data)) THEN
       evals_d=>get_evals_z(ar_data);   eval_out=evals_d(ev_ind)     
    ELSE IF(has_s_cmplx(ar_data)) THEN
       evals_s=>get_evals_c(ar_data);   eval_out=CMPLX(evals_s(ev_ind),KIND=real_8)     
    END IF

  END FUNCTION get_selected_ritz_val

! *****************************************************************************
!> \brief Get all Ritz values of the selected set. eval_out has to be allocated 
!>        at least the size of get_neval_out()
!> \param ar_data ...
!> \param eval_out ...
! *****************************************************************************
  SUBROUTINE get_all_selected_ritz_val(ar_data, eval_out)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    COMPLEX(real_8), DIMENSION(:)            :: eval_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_all_selected_ritz_val', &
      routineP = moduleN//':'//routineN

    COMPLEX(real_4), DIMENSION(:), POINTER   :: evals_s
    COMPLEX(real_8), DIMENSION(:), POINTER   :: evals_d
    INTEGER                                  :: ev_ind, ind
    INTEGER, DIMENSION(:), POINTER           :: selected_ind

    NULLIFY(evals_d,evals_s)
    IF(SIZE(eval_out).LT.get_nval_out(ar_data))&
       CALL dbcsr_abort(routineP,376,'array for eval output too small')
    selected_ind=>get_sel_ind(ar_data)

    IF(has_d_real(ar_data))  evals_d=>get_evals_d(ar_data)
    IF(has_s_real(ar_data))  evals_s=>get_evals_s(ar_data)
    IF(has_d_cmplx(ar_data)) evals_d=>get_evals_d(ar_data)
    IF(has_s_cmplx(ar_data)) evals_s=>get_evals_c(ar_data)

    DO ind=1, get_nval_out(ar_data)
       ev_ind=selected_ind(ind)
       IF(ASSOCIATED(evals_d))eval_out(ind)=evals_d(ev_ind)
       IF(ASSOCIATED(evals_s))eval_out(ind)=CMPLX(evals_s(ev_ind),KIND=real_8)
    END DO

  END SUBROUTINE get_all_selected_ritz_val

! *****************************************************************************
!> \brief ...
!> \param ar_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_arnoldi_initial_vector(ar_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: ar_data
    TYPE(dbcsr_obj)                          :: vector

    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(ar_data)
    control%has_initial_vector=.TRUE.

    IF(has_d_real(ar_data))  CALL set_initial_vector_d(ar_data,vector) 
    IF(has_s_real(ar_data))  CALL set_initial_vector_s(ar_data,vector) 
    IF(has_d_cmplx(ar_data)) CALL set_initial_vector_z(ar_data,vector) 
    IF(has_s_cmplx(ar_data)) CALL set_initial_vector_c(ar_data,vector) 

  END SUBROUTINE set_arnoldi_initial_vector



# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_selection_d.f90" 1

!!! Here come the methods handling the selection of eigenvalues and eigenvectors !!!
!!! If you want a personal method, simply created a Subroutine returning the index
!!! array selected ind which contains as the first nval_out entries the index of the evals

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE select_evals_d(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_evals_d', &
      routineP = moduleN//':'//routineN 
    
    INTEGER                                  :: my_crit, last_el, my_ind, i
    REAL(real_8)                        :: convergence
    TYPE(arnoldi_data_d),POINTER   :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control => get_control(arnoldi_data)
    ar_data => get_data_d(arnoldi_data)

    last_el=control%current_step
    convergence=REAL(0.0, real_8)
    my_crit=control%selection_crit
    control%nval_out=MIN(control%nval_req, control%current_step)
    SELECT CASE(my_crit)
    ! minimum and maximum real eval
    CASE(1)
       CALL index_min_max_real_eval_d(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n maximum real eval
    CASE(2)
       CALL index_nmax_real_eval_d(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n minimum real eval
    CASE(3)
       CALL index_nmin_real_eval_d(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    CASE DEFAULT
       CALL dbcsr_abort(routineP,39,"unknown selection index")
    END SELECT
    ! test whether we are converged
    DO i=1, control%nval_out
       my_ind=control%selected_ind(i)
       convergence=MAX(convergence, &
                   ABS(ar_data%revec(last_el, my_ind)*ar_data%Hessenberg(last_el+1, last_el)))
    END DO
    control%converged=convergence.LT.control%threshold

  END SUBROUTINE select_evals_d

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_min_max_real_eval_d(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_min_max_real_eval_d', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array
    INTEGER                                   :: i

    neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1,current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(1)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO
    DO i=current_step,1,-1
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(2)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO

  END SUBROUTINE index_min_max_real_eval_d

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmax_real_eval_d(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmax_real_eval_d', &
      routineP = moduleN//':'//routineN
    
    INTEGER                                  :: i, nlimit
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(current_step+1-i))))<EPSILON(0.0_real_8))THEN
          selected_ind(i)=indexing(current_step+1-i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmax_real_eval_d

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmin_real_eval_d(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval,nlimit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmin_real_eval_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(i)=indexing(i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmin_real_eval_d

# 415 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_selection_s.f90" 1

!!! Here come the methods handling the selection of eigenvalues and eigenvectors !!!
!!! If you want a personal method, simply created a Subroutine returning the index
!!! array selected ind which contains as the first nval_out entries the index of the evals

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE select_evals_s(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_evals_s', &
      routineP = moduleN//':'//routineN 
    
    INTEGER                                  :: my_crit, last_el, my_ind, i
    REAL(real_4)                        :: convergence
    TYPE(arnoldi_data_s),POINTER   :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control => get_control(arnoldi_data)
    ar_data => get_data_s(arnoldi_data)

    last_el=control%current_step
    convergence=REAL(0.0, real_4)
    my_crit=control%selection_crit
    control%nval_out=MIN(control%nval_req, control%current_step)
    SELECT CASE(my_crit)
    ! minimum and maximum real eval
    CASE(1)
       CALL index_min_max_real_eval_s(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n maximum real eval
    CASE(2)
       CALL index_nmax_real_eval_s(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n minimum real eval
    CASE(3)
       CALL index_nmin_real_eval_s(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    CASE DEFAULT
       CALL dbcsr_abort(routineP,39,"unknown selection index")
    END SELECT
    ! test whether we are converged
    DO i=1, control%nval_out
       my_ind=control%selected_ind(i)
       convergence=MAX(convergence, &
                   ABS(ar_data%revec(last_el, my_ind)*ar_data%Hessenberg(last_el+1, last_el)))
    END DO
    control%converged=convergence.LT.control%threshold

  END SUBROUTINE select_evals_s

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_min_max_real_eval_s(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_min_max_real_eval_s', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array
    INTEGER                                   :: i

    neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1,current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(1)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO
    DO i=current_step,1,-1
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(2)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO

  END SUBROUTINE index_min_max_real_eval_s

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmax_real_eval_s(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmax_real_eval_s', &
      routineP = moduleN//':'//routineN
    
    INTEGER                                  :: i, nlimit
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(current_step+1-i))))<EPSILON(0.0_real_4))THEN
          selected_ind(i)=indexing(current_step+1-i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmax_real_eval_s

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmin_real_eval_s(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval,nlimit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmin_real_eval_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(i)=indexing(i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmin_real_eval_s

# 416 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_selection_z.f90" 1

!!! Here come the methods handling the selection of eigenvalues and eigenvectors !!!
!!! If you want a personal method, simply created a Subroutine returning the index
!!! array selected ind which contains as the first nval_out entries the index of the evals

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE select_evals_z(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_evals_z', &
      routineP = moduleN//':'//routineN 
    
    INTEGER                                  :: my_crit, last_el, my_ind, i
    REAL(real_8)                        :: convergence
    TYPE(arnoldi_data_z),POINTER   :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control => get_control(arnoldi_data)
    ar_data => get_data_z(arnoldi_data)

    last_el=control%current_step
    convergence=REAL(0.0, real_8)
    my_crit=control%selection_crit
    control%nval_out=MIN(control%nval_req, control%current_step)
    SELECT CASE(my_crit)
    ! minimum and maximum real eval
    CASE(1)
       CALL index_min_max_real_eval_z(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n maximum real eval
    CASE(2)
       CALL index_nmax_real_eval_z(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n minimum real eval
    CASE(3)
       CALL index_nmin_real_eval_z(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    CASE DEFAULT
       CALL dbcsr_abort(routineP,39,"unknown selection index")
    END SELECT
    ! test whether we are converged
    DO i=1, control%nval_out
       my_ind=control%selected_ind(i)
       convergence=MAX(convergence, &
                   ABS(ar_data%revec(last_el, my_ind)*ar_data%Hessenberg(last_el+1, last_el)))
    END DO
    control%converged=convergence.LT.control%threshold

  END SUBROUTINE select_evals_z

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_min_max_real_eval_z(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_min_max_real_eval_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array
    INTEGER                                   :: i

    neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1,current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(1)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO
    DO i=current_step,1,-1
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(2)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO

  END SUBROUTINE index_min_max_real_eval_z

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmax_real_eval_z(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmax_real_eval_z', &
      routineP = moduleN//':'//routineN
    
    INTEGER                                  :: i, nlimit
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(current_step+1-i))))<EPSILON(0.0_real_8))THEN
          selected_ind(i)=indexing(current_step+1-i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmax_real_eval_z

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmin_real_eval_z(evals, current_step, selected_ind, neval)
    COMPLEX(real_8), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval,nlimit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmin_real_eval_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_8), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_8)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_8))THEN
          selected_ind(i)=indexing(i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmin_real_eval_z

# 417 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_selection_c.f90" 1

!!! Here come the methods handling the selection of eigenvalues and eigenvectors !!!
!!! If you want a personal method, simply created a Subroutine returning the index
!!! array selected ind which contains as the first nval_out entries the index of the evals

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE select_evals_c(arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'select_evals_c', &
      routineP = moduleN//':'//routineN 
    
    INTEGER                                  :: my_crit, last_el, my_ind, i
    REAL(real_4)                        :: convergence
    TYPE(arnoldi_data_c),POINTER   :: ar_data
    TYPE(arnoldi_control), POINTER           :: control

    control => get_control(arnoldi_data)
    ar_data => get_data_c(arnoldi_data)

    last_el=control%current_step
    convergence=REAL(0.0, real_4)
    my_crit=control%selection_crit
    control%nval_out=MIN(control%nval_req, control%current_step)
    SELECT CASE(my_crit)
    ! minimum and maximum real eval
    CASE(1)
       CALL index_min_max_real_eval_c(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n maximum real eval
    CASE(2)
       CALL index_nmax_real_eval_c(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    ! n minimum real eval
    CASE(3)
       CALL index_nmin_real_eval_c(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
    CASE DEFAULT
       CALL dbcsr_abort(routineP,39,"unknown selection index")
    END SELECT
    ! test whether we are converged
    DO i=1, control%nval_out
       my_ind=control%selected_ind(i)
       convergence=MAX(convergence, &
                   ABS(ar_data%revec(last_el, my_ind)*ar_data%Hessenberg(last_el+1, last_el)))
    END DO
    control%converged=convergence.LT.control%threshold

  END SUBROUTINE select_evals_c

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_min_max_real_eval_c(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_min_max_real_eval_c', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array
    INTEGER                                   :: i

    neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1,current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(1)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO
    DO i=current_step,1,-1
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(2)=indexing(i)
          neval=neval+1
          EXIT
       END IF
    END DO

  END SUBROUTINE index_min_max_real_eval_c

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmax_real_eval_c(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmax_real_eval_c', &
      routineP = moduleN//':'//routineN
    
    INTEGER                                  :: i, nlimit
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(current_step+1-i))))<EPSILON(0.0_real_4))THEN
          selected_ind(i)=indexing(current_step+1-i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmax_real_eval_c

! *****************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! *****************************************************************************
  SUBROUTINE index_nmin_real_eval_c(evals, current_step, selected_ind, neval)
    COMPLEX(real_4), DIMENSION(:)       :: evals
    INTEGER                                  :: current_step
    INTEGER, DIMENSION(:)                    :: selected_ind
    INTEGER                                  :: neval,nlimit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'index_nmin_real_eval_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    INTEGER, DIMENSION(current_step)         :: indexing
    REAL(real_4), DIMENSION(current_step)        :: tmp_array

    nlimit=neval; neval=0
    selected_ind=0
    tmp_array(1:current_step)=REAL(evals(1:current_step), real_4)
    CALL sort(tmp_array, current_step, indexing)
    DO i=1, current_step
       IF(ABS(AIMAG(evals(indexing(i))))<EPSILON(0.0_real_4))THEN
          selected_ind(i)=indexing(i)
          neval=neval+1
          IF(neval==nlimit)EXIT
       END IF
    END DO

  END SUBROUTINE index_nmin_real_eval_c

# 418 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_manipulation_d.f90" 1
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
! *****************************************************************************
  SUBROUTINE  setup_arnoldi_data_d (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_d), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%x_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_d(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_d

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data_d (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_d', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_d), POINTER            :: ar_data

    ar_data=>get_data_d(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%x_vec))DEALLOCATE(ar_data%x_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_d

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector_d(arnoldi_data,ind,matrix,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector_d', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_d), POINTER           :: ar_data
    INTEGER                                           :: vsize, myind, sspace_size, i
    INTEGER, DIMENSION(:), POINTER           :: selected_ind
    COMPLEX(real_8),DIMENSION(:),ALLOCATABLE       :: ritz_v
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    selected_ind=>get_sel_ind(arnoldi_data)
    ar_data=>get_data_d(arnoldi_data)
    sspace_size=get_subsp_size(arnoldi_data)
    vsize=SIZE(ar_data%f_vec)
    myind=selected_ind(ind)
    ALLOCATE(ritz_v(vsize))
    ritz_v=CMPLX(0.0,0.0,real_8)

    IF(dbcsr_is_initialized (vector))CALL dbcsr_release(vector)
    CALL create_col_vec_from_matrix(vector,matrix,1)
    IF(control%local_comp)THEN
       DO i=1,sspace_size
          ritz_v(:)=ritz_v(:)+ar_data%local_history(:,i)*ar_data%revec(i,myind)
       END DO
       data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=0.0_real_8)
       ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
       ! again I hope the user knows what is required
       data_vec(1:vsize) =REAL(ritz_v(1:vsize),KIND=real_8)
    END IF

    DEALLOCATE(ritz_v)

  END SUBROUTINE get_selected_ritz_vector_d
     
   
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_initial_vector_d(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_initial_vector_d', &
      routineP = moduleN//':'//routineN
    
    TYPE(arnoldi_data_d), POINTER           :: ar_data
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    INTEGER                                           :: nrow_local, ncol_local
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)

    CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ar_data=>get_data_d(arnoldi_data)
    data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=0.0_real_8)
    IF(nrow_local*ncol_local>0)ar_data%f_vec(1:nrow_local)=data_vec(1:nrow_local)

  END SUBROUTINE set_initial_vector_d  

    
# 419 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_manipulation_s.f90" 1
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
! *****************************************************************************
  SUBROUTINE  setup_arnoldi_data_s (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_s), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%x_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_s(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_s

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data_s (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER            :: ar_data

    ar_data=>get_data_s(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%x_vec))DEALLOCATE(ar_data%x_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_s

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector_s(arnoldi_data,ind,matrix,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector_s', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_s), POINTER           :: ar_data
    INTEGER                                           :: vsize, myind, sspace_size, i
    INTEGER, DIMENSION(:), POINTER           :: selected_ind
    COMPLEX(real_4),DIMENSION(:),ALLOCATABLE       :: ritz_v
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    selected_ind=>get_sel_ind(arnoldi_data)
    ar_data=>get_data_s(arnoldi_data)
    sspace_size=get_subsp_size(arnoldi_data)
    vsize=SIZE(ar_data%f_vec)
    myind=selected_ind(ind)
    ALLOCATE(ritz_v(vsize))
    ritz_v=CMPLX(0.0,0.0,real_4)

    IF(dbcsr_is_initialized (vector))CALL dbcsr_release(vector)
    CALL create_col_vec_from_matrix(vector,matrix,1)
    IF(control%local_comp)THEN
       DO i=1,sspace_size
          ritz_v(:)=ritz_v(:)+ar_data%local_history(:,i)*ar_data%revec(i,myind)
       END DO
       data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=0.0_real_4)
       ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
       ! again I hope the user knows what is required
       data_vec(1:vsize) =REAL(ritz_v(1:vsize),KIND=real_4)
    END IF

    DEALLOCATE(ritz_v)

  END SUBROUTINE get_selected_ritz_vector_s
     
   
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_initial_vector_s(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_initial_vector_s', &
      routineP = moduleN//':'//routineN
    
    TYPE(arnoldi_data_s), POINTER           :: ar_data
    REAL(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    INTEGER                                           :: nrow_local, ncol_local
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)

    CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ar_data=>get_data_s(arnoldi_data)
    data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=0.0_real_4)
    IF(nrow_local*ncol_local>0)ar_data%f_vec(1:nrow_local)=data_vec(1:nrow_local)

  END SUBROUTINE set_initial_vector_s  

    
# 420 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_manipulation_z.f90" 1
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
! *****************************************************************************
  SUBROUTINE  setup_arnoldi_data_z (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_z), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%x_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_z(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_z

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data_z (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_z', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_z), POINTER            :: ar_data

    ar_data=>get_data_z(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%x_vec))DEALLOCATE(ar_data%x_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_z

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector_z(arnoldi_data,ind,matrix,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector_z', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_z), POINTER           :: ar_data
    INTEGER                                           :: vsize, myind, sspace_size, i
    INTEGER, DIMENSION(:), POINTER           :: selected_ind
    COMPLEX(real_8),DIMENSION(:),ALLOCATABLE       :: ritz_v
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    selected_ind=>get_sel_ind(arnoldi_data)
    ar_data=>get_data_z(arnoldi_data)
    sspace_size=get_subsp_size(arnoldi_data)
    vsize=SIZE(ar_data%f_vec)
    myind=selected_ind(ind)
    ALLOCATE(ritz_v(vsize))
    ritz_v=CMPLX(0.0,0.0,real_8)

    IF(dbcsr_is_initialized (vector))CALL dbcsr_release(vector)
    CALL create_col_vec_from_matrix(vector,matrix,1)
    IF(control%local_comp)THEN
       DO i=1,sspace_size
          ritz_v(:)=ritz_v(:)+ar_data%local_history(:,i)*ar_data%revec(i,myind)
       END DO
       data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_8))
       ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
       ! again I hope the user knows what is required
       data_vec(1:vsize) =CMPLX(ritz_v(1:vsize),KIND=real_8)
    END IF

    DEALLOCATE(ritz_v)

  END SUBROUTINE get_selected_ritz_vector_z
     
   
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_initial_vector_z(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_initial_vector_z', &
      routineP = moduleN//':'//routineN
    
    TYPE(arnoldi_data_z), POINTER           :: ar_data
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    INTEGER                                           :: nrow_local, ncol_local
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)

    CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ar_data=>get_data_z(arnoldi_data)
    data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_8))
    IF(nrow_local*ncol_local>0)ar_data%f_vec(1:nrow_local)=data_vec(1:nrow_local)

  END SUBROUTINE set_initial_vector_z  

    
# 421 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_manipulation_c.f90" 1
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
! *****************************************************************************
  SUBROUTINE  setup_arnoldi_data_c (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_c), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%x_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_c(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_c

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
! *****************************************************************************
  SUBROUTINE deallocate_arnoldi_data_c (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_c', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_c), POINTER            :: ar_data

    ar_data=>get_data_c(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%x_vec))DEALLOCATE(ar_data%x_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_c

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vector_c(arnoldi_data,ind,matrix,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(dbcsr_obj)                          :: matrix
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_selected_ritz_vector_c', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_c), POINTER           :: ar_data
    INTEGER                                           :: vsize, myind, sspace_size, i
    INTEGER, DIMENSION(:), POINTER           :: selected_ind
    COMPLEX(real_4),DIMENSION(:),ALLOCATABLE       :: ritz_v
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)
    selected_ind=>get_sel_ind(arnoldi_data)
    ar_data=>get_data_c(arnoldi_data)
    sspace_size=get_subsp_size(arnoldi_data)
    vsize=SIZE(ar_data%f_vec)
    myind=selected_ind(ind)
    ALLOCATE(ritz_v(vsize))
    ritz_v=CMPLX(0.0,0.0,real_4)

    IF(dbcsr_is_initialized (vector))CALL dbcsr_release(vector)
    CALL create_col_vec_from_matrix(vector,matrix,1)
    IF(control%local_comp)THEN
       DO i=1,sspace_size
          ritz_v(:)=ritz_v(:)+ar_data%local_history(:,i)*ar_data%revec(i,myind)
       END DO
       data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_4))
       ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
       ! again I hope the user knows what is required
       data_vec(1:vsize) =CMPLX(ritz_v(1:vsize),KIND=real_4)
    END IF

    DEALLOCATE(ritz_v)

  END SUBROUTINE get_selected_ritz_vector_c
     
   
! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE set_initial_vector_c(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj)                          :: vector

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_initial_vector_c', &
      routineP = moduleN//':'//routineN
    
    TYPE(arnoldi_data_c), POINTER           :: ar_data
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: data_vec
    INTEGER                                           :: nrow_local, ncol_local
    TYPE(arnoldi_control), POINTER           :: control

    control=>get_control(arnoldi_data)

    CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
    ar_data=>get_data_c(arnoldi_data)
    data_vec => dbcsr_get_data_p (vector%m%data_area, select_data_type=CMPLX(0.0, 0.0, real_4))
    IF(nrow_local*ncol_local>0)ar_data%f_vec(1:nrow_local)=data_vec(1:nrow_local)

  END SUBROUTINE set_initial_vector_c  

    
# 422 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/arnoldi/dbcsr_arnoldi_data_methods.F" 2

END MODULE dbcsr_arnoldi_data_methods


   


