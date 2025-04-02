# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief used for collecting some of the diagonalization shemes available for
!>      cp_fm_type. cp_fm_power also moved here as it is very related
!> \note
!>      first version : most routines imported
!> \par History
!>      - unused Jacobi routines removed, cosmetics (05.04.06,MK)
!> \author Joost VandeVondele (2003-08)
! *****************************************************************************
MODULE cp_fm_diag
USE cp_log_handling, ONLY: cp_logger_type, cp_get_default_logger, cp_logger_get_default_unit_nr, cp_logger_get_unit_nr
  USE machine,                         ONLY: m_memory
  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinit
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk,&
                                             cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_column_scale,&
                                             cp_fm_triangular_multiply,&
                                             cp_fm_triangular_invert,&
                                             cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_get
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_get_info,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length, dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_free,&
                                             mp_comm_split,&
                                             mp_comm_split_direct,&
                                             mp_sync











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
# 65 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag'

  ! these saved variables are DIAGONALIZATION global
   INTEGER, SAVE :: diag_type = 0
   INTEGER, SAVE :: kernel_type = 0

  ! Public subroutines

  PUBLIC :: choose_eigv_solver,&
            cp_fm_block_jacobi,&
            cp_fm_power,&
            cp_fm_syevd,&
            cp_fm_syevx,&
            cp_fm_geeig,&
            cp_fm_geeig_canon,&
            diag_init

CONTAINS

! *****************************************************************************
!> \brief Setup the diagonalization library to be used
!>         Check of availability not yet fully implemented
!>         It should change library to Scalapack if others are not available
!> \param diag_lib diag_library flag from GLOBAL section in input
!> \param switched ...
!> \param k_elpa ...
!> \author  MI 11.2013       
! *****************************************************************************
  SUBROUTINE diag_init(diag_lib,switched,k_elpa)
    CHARACTER(LEN=*), INTENT(IN)             :: diag_lib
    LOGICAL, INTENT(INOUT)                   :: switched
    INTEGER, INTENT(IN)                      :: k_elpa

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: try_library

! skalapack is the default and always linked

    IF(diag_lib.EQ."SL") THEN
      try_library = .FALSE.
      diag_type = 1
    ELSE
      try_library =.TRUE.
    END IF

    IF(try_library) THEN
        IF(diag_lib.EQ."ELPA") THEN




          ! ELPA library not linked, switch to SL
             diag_type = 1
             switched = .TRUE.

         ELSE IF(diag_lib.EQ."SL2") THEN




          ! SL2 library not linked, switch to SL
             diag_type = 1
             switched = .TRUE.

         END IF
    END IF

    IF(diag_type == 3) THEN
       kernel_type = k_elpa
    END IF

    ! Check that one of the diagonalization type is set
    IF(diag_type<1) THEN
      ! something wrong
      CALL cp__b("fm/cp_fm_diag.F",146,"Unknown DIAG type")
    END IF

  END SUBROUTINE diag_init


! *****************************************************************************
!> \brief   Choose the Eigensolver depending on which library is available
!>           ELPA seems to be unstable for small systems
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for cp_fm_syevd with scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! *****************************************************************************
  SUBROUTINE choose_eigv_solver(matrix,eigenvectors,eigenvalues,info)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(OUT), OPTIONAL           :: info

    CHARACTER(len=*), PARAMETER :: routineN = 'choose_eigv_solver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: myinfo, nmo

    myinfo = 0

    nmo = SIZE(eigenvalues,1)

    !sample peak memory
    CALL m_memory()

    IF(diag_type==3) THEN
      CALL cp_fm_elpa(matrix,eigenvectors,eigenvalues)
    ELSE IF(diag_type==2) THEN
      CALL cp_fm_syevr(matrix,eigenvectors,eigenvalues,1,nmo)
    ELSE IF(diag_type==1) THEN
      CALL cp_fm_syevd(matrix,eigenvectors,eigenvalues,info=myinfo)
    END IF

    IF (PRESENT(info)) info = myinfo 
  END SUBROUTINE choose_eigv_solver

! *****************************************************************************
!> \brief   Computes all eigenvalues and vectors of a real symmetric matrix
!>          significantly faster than syevx, scales also much better.
!>          Needs workspace to allocate all the eigenvectors
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param info ...
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!> \par     info If present returns error code and prevents program stops.
!>               Works currently only for scalapack.
!>               Other solvers will end the program regardless of PRESENT(info).
! *****************************************************************************
  SUBROUTINE cp_fm_syevd(matrix,eigenvectors,eigenvalues,info)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    INTEGER, INTENT(OUT), OPTIONAL           :: info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, myinfo, n, nmo
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig
# 228 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"
    INTEGER                                  :: liwork, lwork
    INTEGER, DIMENSION(:), POINTER           :: iwork
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work


    CALL timeset(routineN,handle)

    myinfo = 0

    n = matrix%matrix_struct%nrow_global
    ALLOCATE(eig(n))


# 335 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    !MK Retrieve the optimal work array sizes first
    myinfo = 0
    lwork = -1
    liwork = -1
    m => matrix%local_data
    eig(:) = 0.0_dp

    ALLOCATE (work(1))
    work(:) = 0.0_dp
    ALLOCATE (iwork(1))
    iwork(:) = 0

    CALL dsyevd('V','U',n,m(1,1),n,eig(1),work(1),lwork,iwork(1),liwork,myinfo)

    IF (myinfo /= 0) THEN
       CALL cp__b("fm/cp_fm_diag.F",351,"ERROR in DSYEVD: Could not retrieve work array sizes")
    END IF

    ! Reallocate work arrays and perform diagonalisation
    lwork = INT(work(1))
    DEALLOCATE (work)
    ALLOCATE(work(lwork))

    liwork = iwork(1)
    DEALLOCATE (iwork)
    ALLOCATE(iwork(liwork))
    iwork(:) = 0

    CALL dsyevd('V','U',n,m(1,1),n,eig(1),work(1),lwork,iwork(1),liwork,myinfo)

    IF (myinfo /= 0) THEN
      CALL cp__b("fm/cp_fm_diag.F",367,"Matrix diagonalization failed")
    END IF

    CALL cp_fm_to_fm(matrix,eigenvectors)

    DEALLOCATE (iwork)
    DEALLOCATE (work)


    IF(PRESENT(info)) myinfo = 0

    nmo = SIZE(eigenvalues,1)
    IF (nmo > n) THEN
      eigenvalues(1:n) = eig(1:n)
    ELSE
      eigenvalues(1:nmo) = eig(1:nmo)
    END IF

    DEALLOCATE (eig)
    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevd

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eig ...
!> \param info ...
! *****************************************************************************
  SUBROUTINE cp_fm_syevd_base(matrix,eigenvectors,eig,info)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eig
    INTEGER, INTENT(OUT)                     :: info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syevd_base', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
# 419 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    CALL timeset(routineN,handle)


# 476 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(matrix))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(eigenvectors))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(eig))==-1) EXIT ;  END DO ; ENDIF
   info = -1
   CALL cp__b("fm/cp_fm_diag.F",480,"cp_fm_syevd_base called with out SCALAPACK compiled in")


    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevd_base

! *****************************************************************************
!> \brief   compute eigenvalues and optionally eigenvectors of a real symmetric matrix using scalapack.
!>          If eigenvectors are required this routine will replicate a full matrix on each CPU...
!>          if more than a handful of vectors are needed, use cp_fm_syevd instead
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param neig ...
!> \param work_syevx ...
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          neig   is the number of vectors needed (default all)
!>          work_syevx evec calculation only, is the fraction of the working buffer allowed (1.0 use full buffer)
!>                     reducing this saves time, but might cause the routine to fail
! *****************************************************************************
  SUBROUTINE cp_fm_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx)

    ! Diagonalise the symmetric n by n matrix using the LAPACK library.

    TYPE(cp_fm_type), POINTER                    :: matrix
    TYPE(cp_fm_type), POINTER, OPTIONAL          :: eigenvectors
    REAL(KIND = dp), OPTIONAL, INTENT(IN)        :: work_syevx
    INTEGER, INTENT(IN), OPTIONAL                :: neig
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevx",&
                                   routineP = moduleN//":"//routineN

    REAL(KIND=dp), PARAMETER                     :: orfac = -1.0_dp,&
                                                    vl = 0.0_dp,&
                                                    vu = 0.0_dp

    TYPE(cp_blacs_env_type), POINTER             :: context
    TYPE(cp_logger_type), POINTER                :: logger
    CHARACTER(LEN=1)                             :: job_type
    REAL(KIND=dp)                                :: abstol, work_syevx_local
    INTEGER                                      :: handle, info, &
                                                    liwork, lwork, m, mypcol, &
                                                    myprow, n, nb, npcol, &
                                                    nprow, output_unit, &
                                                    neig_local
    LOGICAL                                      :: ionode, needs_evecs
    INTEGER, DIMENSION(:), ALLOCATABLE           :: ifail, iwork
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: w, work
    REAL(KIND=dp), DIMENSION(:,:), POINTER       :: a, z

    REAL(KIND=dp), EXTERNAL                      :: dlamch

# 544 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"
    INTEGER, EXTERNAL                            :: ilaenv


    ! by default all
    n = matrix%matrix_struct%nrow_global
    neig_local=n
    IF (PRESENT(neig)) neig_local=neig
    IF (neig_local == 0) RETURN

    CALL timeset(routineN,handle)

    needs_evecs=PRESENT(eigenvectors)

    logger => cp_get_default_logger()
    ionode = logger%para_env%mepos==logger%para_env%source
    n = matrix%matrix_struct%nrow_global

    ! by default allocate all needed space
    work_syevx_local=1.0_dp
    IF (PRESENT(work_syevx)) work_syevx_local=work_syevx

    ! set scalapack job type
    IF (needs_evecs) THEN
       job_type="V"
    ELSE
       job_type="N"
    ENDIF

    ! target the most accurate calculation of the eigenvalues
    abstol = 2.0_dp*dlamch("S")


    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)


    ALLOCATE (w(n))
    eigenvalues(:) = 0.0_dp
# 673 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    a => matrix%local_data
    IF (needs_evecs) THEN
       z => eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => matrix%local_data
    ENDIF

    ! Get the optimal work storage size

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n))
    ifail = 0
    ALLOCATE (iwork(liwork))
    ALLOCATE (work(lwork))
    info = 0

    CALL dsyevx(job_type,"I","U",n,a(1,1),n,vl,vu,1,neig_local,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

    ! Error handling

    IF (info /= 0) THEN
      output_unit = cp_logger_get_unit_nr(logger,local=.FALSE.)
      WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
        "info    = ",info
      IF (info > 0) THEN
        WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
          "ifail   = ",ifail
      END IF
      CALL cp__b("fm/cp_fm_diag.F",709,"Error in dsyevx")
    END IF

    ! Release work storage

    DEALLOCATE (ifail)
    DEALLOCATE (iwork)
    DEALLOCATE (work)


    eigenvalues(1:neig_local) = w(1:neig_local)
    DEALLOCATE (w)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevx

! *****************************************************************************
!> \brief  computes selected eigenvalues and, optionally, eigenvectors of
!>        a real symmetric matrix A distributed in 2D blockcyclic format by
!>       calling the recommended sequence of ScaLAPACK routines.
!>
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param ilow ...
!> \param iup ...
!> \par     matrix is supposed to be in upper triangular form, and overwritten by this routine
!>          subsets of eigenvalues/vectors can be selected by
!>          specifying a range of values or a range of indices for the desired eigenvalues.
! *****************************************************************************
  SUBROUTINE cp_fm_syevr(matrix,eigenvectors,eigenvalues,ilow,iup)

    TYPE(cp_fm_type), POINTER                    :: matrix
    TYPE(cp_fm_type), POINTER, OPTIONAL          :: eigenvectors
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)   :: eigenvalues
    INTEGER, INTENT(IN), OPTIONAL                :: ilow,iup

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevr",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND = dp), PARAMETER  :: vl = 0.0_dp,&
                                   vu = 0.0_dp

    CHARACTER(LEN=1)                           :: job_type
    INTEGER                                    :: handle, ilow_local,&
                                                  info, iup_local,&
                                                  lwork, liwork, mypcol,&
                                                  myprow, n, neig
    INTEGER, DIMENSION(:), ALLOCATABLE         :: iwork
    LOGICAL                                    :: ionode, needs_evecs
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: w, work
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: a, z

    TYPE(cp_blacs_env_type), POINTER           :: context
    TYPE(cp_logger_type), POINTER              :: logger

    REAL(KIND = dp), EXTERNAL :: dlamch







    INTEGER                             :: m, nb
    REAL(dp)                            :: abstol
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ifail
    INTEGER, EXTERNAL                   :: ilaenv


    ! by default all
    n = matrix%matrix_struct%nrow_global
    neig=n
    iup_local = n
    ilow_local = 1
    IF (PRESENT(ilow) .AND. PRESENT(iup)) THEN
        neig=iup-ilow+1
        iup_local = iup
        ilow_local = ilow
    END IF
    IF (neig <= 0) RETURN

    CALL timeset(routineN,handle)

    needs_evecs=PRESENT(eigenvectors)

    logger => cp_get_default_logger()
    ionode = logger%para_env%mepos==logger%para_env%source
    n = matrix%matrix_struct%nrow_global

    ! set scalapack job type
    IF (needs_evecs) THEN
       job_type="V"
    ELSE
       job_type="N"
    ENDIF

    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    ALLOCATE(w(n))

    eigenvalues(:) = 0.0_dp

# 867 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    a => matrix%local_data
    IF (needs_evecs) THEN
       z => eigenvectors%local_data
    ELSE
       ! z will not be referenced
       z => matrix%local_data
    ENDIF

    ! Get the optimal work storage size

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n))
    ifail = 0

    ALLOCATE (iwork(liwork))
    ALLOCATE (work(lwork))

    ! target the most accurate calculation of the eigenvalues
    abstol = 2.0_dp*dlamch("S")

    info = 0
    CALL dsyevx(job_type,"I","U",n,a(1,1),n,vl,vu,ilow_local,iup_local,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

    ! Error handling
    IF(.NOT.(info==0))CALL cp__a("fm/cp_fm_diag.F",898)

    ! Release work storage
    DEALLOCATE (iwork)
    DEALLOCATE (work)



    eigenvalues(ilow_local:iup_local) = w(ilow_local:iup_local)
    DEALLOCATE (w)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_syevr

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! *****************************************************************************
 SUBROUTINE cp_fm_elpa(matrix,eigenvectors,eigenvalues)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eigenvalues

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_elpa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
# 938 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    CALL timeset(routineN,handle)
# 995 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(matrix))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(eigenvectors))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(eigenvalues))==-1) EXIT ;  END DO ; ENDIF

    CALL cp__b("fm/cp_fm_diag.F",1000,"CP2K compiled without the ELPA library.")


    CALL timestop(handle)

  END SUBROUTINE cp_fm_elpa


! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param work ...
!> \param exponent ...
!> \param threshold ...
!> \param n_dependent ...
!> \param verbose ...
! *****************************************************************************
  SUBROUTINE cp_fm_power(matrix,work,exponent,threshold,n_dependent,verbose)

   ! Raise the real symmetric n by n matrix to the power given by
   ! the exponent. All eigenvectors with a corresponding eigenvalue lower
   ! than threshold are quenched. result in matrix

   ! - Creation (29.03.1999, Matthias Krack)
   ! - Parallelised using BLACS and ScaLAPACK (06.06.2001,MK)

   TYPE(cp_fm_type), POINTER                 :: matrix,work
   REAL(KIND = dp), INTENT(IN)               :: exponent,threshold
   INTEGER, INTENT(OUT)                      :: n_dependent
   LOGICAL, INTENT(IN), OPTIONAL             :: verbose

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_power', &
      routineP = moduleN//':'//routineN

   INTEGER                                    :: handle,icol_global,&
                                                 mypcol,myprow,&
                                                 ncol_global,npcol,nprow,&
                                                 nrow_global
   LOGICAL                                    :: my_verbose
   REAL(KIND=dp)                              :: condition_number,f,p
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigenvalues
   REAL(KIND=dp), DIMENSION(:,:), POINTER     :: eigenvectors
   TYPE(cp_blacs_env_type), POINTER           :: context







    CALL timeset(routineN,handle)

    my_verbose = .FALSE.
    IF(PRESENT(verbose)) my_verbose = verbose

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n_dependent = 0
    p = 0.5_dp*exponent

    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    ALLOCATE (eigenvalues(ncol_global))

    eigenvalues(:) = 0.0_dp

    ! Compute the eigenvectors and eigenvalues

    CALL choose_eigv_solver(matrix,work,eigenvalues)

# 1132 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    eigenvectors => work%local_data

    ! Build matrix**exponent with eigenvector quenching

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1
        eigenvectors(1:nrow_global,icol_global) = 0.0_dp

      ELSE

        f = eigenvalues(icol_global)**p
        eigenvectors(1:nrow_global,icol_global) =&
          f*eigenvectors(1:nrow_global,icol_global)

      END IF

    END DO


    CALL cp_fm_syrk("U","N",ncol_global,1.0_dp,work,1,1,0.0_dp,matrix)
    CALL cp_fm_upper_to_full(matrix,work)

    ! Print some warnings/notes

    IF (matrix%matrix_struct%para_env%mepos == 0 .AND. my_verbose) THEN
      condition_number = ABS(eigenvalues(ncol_global)/eigenvalues(1))
      WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,(T2,A,ES15.6))")&
        "CP_FM_POWER: smallest eigenvalue:",eigenvalues(1),&
        "CP_FM_POWER: largest eigenvalue: ",eigenvalues(ncol_global),&
        "CP_FM_POWER: condition number:   ",condition_number
      IF (eigenvalues(1) <= 0.0_dp) THEN
        WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,T2,A)")&
          "WARNING: matrix has a negative eigenvalue, tighten EPS_DEFAULT"
      END IF
      IF (condition_number > 1.0E12_dp) THEN
        WRITE (UNIT=cp_logger_get_default_unit_nr(),FMT="(/,T2,A)")&
          "WARNING: high condition number => possibly ill-conditioned matrix"
      END IF
    END IF

    DEALLOCATE (eigenvalues)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_power

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigval ...
!> \param thresh ...
!> \param start_sec_block ...
! *****************************************************************************
  SUBROUTINE cp_fm_block_jacobi(matrix,eigenvectors,eigval,thresh,&
                                start_sec_block)

    ! Calculates block diagonalizazion from full symmetric matrix
    ! It has its origin in cp_fm_syevx. This routine rotates only elements
    ! which are larger than a threshold thresh.
    ! start_sec_block is the start of the second block.
    ! IT DOES ONLY ONE SWEEP!

    ! - Creation (07.10.2002, Martin Fengler)
    ! - Cosmetics (05.04.06,MK)

    TYPE(cp_fm_type), POINTER                 :: eigenvectors,matrix
    REAL(KIND = dp), DIMENSION(:), INTENT(IN) :: eigval
    INTEGER, INTENT(IN)                       :: start_sec_block
    REAL(KIND = dp), INTENT(IN)               :: thresh

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_block_jacobi', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,ev

    REAL(KIND = dp) :: tan_theta,tau,c,s
    INTEGER  :: q,p,N
    REAL(KIND = dp), DIMENSION (:),ALLOCATABLE :: c_ip

# 1228 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

# 1343 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_diag.F"

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die bearbeitet werden soll

    ALLOCATE(c_ip(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data        ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

    ! Start diagonalizing matrix

    tan_theta = 0.0_dp
    tau       = 0.0_dp

    DO q=start_sec_block,N
      DO p=1,(start_sec_block-1)

          IF(ABS(A(p,q))>thresh) THEN

          tau = (eigval(q)-eigval(p))/(2.0_dp*A(p,q))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV(:,p) - s*EV(:,q)
          !  c_iq = s*EV(:,p) + c*EV(:,q)

          !  EV(:,p)=c_ip
          !  EV(:,q)=c_iq

          CALL dcopy(N,EV(1,p),1,c_ip(1),1)
          CALL dscal(N,c,EV(1,p),1)
          CALL daxpy(N,-s,EV(1,q),1,EV(1,p),1)
          CALL dscal(N,c,EV(1,q),1)
          CALL daxpy(N,s,c_ip(1),1,EV(1,q),1)

         END IF

      END DO
    END DO

    ! Release work storage

    DEALLOCATE (c_ip)



    CALL timestop(handle)

  END SUBROUTINE cp_fm_block_jacobi

! *****************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Single option version: Cholesky decomposition of B
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
! *****************************************************************************
 SUBROUTINE cp_fm_geeig(amatrix,bmatrix,eigenvectors,eigenvalues,work)

    TYPE(cp_fm_type), POINTER                :: amatrix, bmatrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eigenvalues
    TYPE(cp_fm_type), POINTER                :: work

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_geeig', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo

    CALL timeset(routineN,handle)

    CALL cp_fm_get_info(amatrix,nrow_global=nao)
    nmo = SIZE(eigenvalues)
    ! Cholesky decompose S=U(T)U
    CALL cp_fm_cholesky_decompose(bmatrix)
    ! Invert to get U^(-1)
    CALL cp_fm_triangular_invert(bmatrix)
    ! Reduce to get U^(-T) * H * U^(-1)
    CALL cp_fm_triangular_multiply(bmatrix,amatrix,side="R")
    CALL cp_fm_triangular_multiply(bmatrix,amatrix,transpose_tr=.TRUE.)
    ! Diagonalize
    CALL choose_eigv_solver(matrix=amatrix,eigenvectors=work,&
                            eigenvalues=eigenvalues)
    ! Restore vectors C = U^(-1) * C*
    CALL cp_fm_triangular_multiply(bmatrix,work)
    CALL cp_fm_to_fm(work,eigenvectors,nmo)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_geeig

! *****************************************************************************
!> \brief General Eigenvalue Problem  AX = BXE
!>        Use canonical diagonalization : U*s**(-1/2)
!> \param amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param work ...
!> \param epseig ...
! *****************************************************************************
 SUBROUTINE cp_fm_geeig_canon(amatrix,bmatrix,eigenvectors,eigenvalues,work,epseig)

    TYPE(cp_fm_type), POINTER                :: amatrix, bmatrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:)              :: eigenvalues
    TYPE(cp_fm_type), POINTER                :: work
    REAL(KIND=dp), INTENT(IN)                :: epseig

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_geeig_canon', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icol, irow, nao, &
                                                nc, ncol, nmo, nx
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: seigval

    CALL timeset(routineN,handle)

    ! Test sizees
    CALL cp_fm_get_info(amatrix,nrow_global=nao)
    nmo = SIZE(eigenvalues)
    ALLOCATE(seigval(nao))

    ! Diagonalize -S matrix, this way the NULL space is at the end of the spectrum
    CALL cp_fm_scale(-1.0_dp, bmatrix)
    CALL choose_eigv_solver(matrix=bmatrix,eigenvectors=work,eigenvalues=seigval)
    seigval(:) = -seigval(:)
    nc=nao
    DO i=1,nao
       IF (seigval(i) < epseig) THEN
          nc = i-1
          EXIT
       END IF
    END DO
    IF(.NOT.(nc/=0))CALL cp__a("fm/cp_fm_diag.F",1483)

    IF(nc/=nao) THEN
       IF(nc < nmo) THEN
          ! Copy NULL space definition to last vectors of eigenvectors (if needed)
          ncol = nmo - nc
          CALL cp_fm_to_fm(work,eigenvectors,ncol,nc+1,nc+1)
       END IF
       ! Set NULL space in eigenvector matrix of S to zero
       DO icol=nc+1,nao
          DO irow=1,nao
             CALL cp_fm_set_element(work,irow,icol,0.0_dp)
          END DO
       END DO
       ! Set small eigenvalues to a dummy save value
       seigval(nc+1:nao) = 1.0_dp
    END IF
    ! calculate U*s**(-1/2)
    seigval(:) = 1.0_dp/SQRT(seigval(:))
    CALL cp_fm_column_scale(work,seigval)
    ! Reduce to get U^(-T) * H * U^(-1)
    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,work,amatrix,0.0_dp,bmatrix)
    CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,bmatrix,work,0.0_dp,amatrix)
    IF(nc/=nao) THEN
       ! set diagonal values to save large value
       DO icol=nc+1,nao
          CALL cp_fm_set_element(amatrix,icol,icol,10000.0_dp)
       END DO
    END IF
    ! Diagonalize
    CALL choose_eigv_solver(matrix=amatrix,eigenvectors=bmatrix,eigenvalues=eigenvalues)
    nx = MIN(nc,nmo)
    ! Restore vectors C = U^(-1) * C*
    CALL cp_fm_gemm("N","N",nao,nx,nc,1.0_dp,work,bmatrix,0.0_dp,eigenvectors)

    DEALLOCATE(seigval)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_geeig_canon

END MODULE cp_fm_diag
