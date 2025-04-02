# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!






MODULE fftw3_lib

  USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                             C_INT,&
                                             C_INTPTR_T
  USE cp_files,                        ONLY: get_unit_number
  USE fft_kinds,                       ONLY: dp,&
                                             integer8_kind
  USE fft_plan,                        ONLY: fft_plan_type

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/../../base/base_uses.f90" 1
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
# 24 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F" 2

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fftw3_do_init, fftw3_do_cleanup, fftw3_get_lengths, fftw33d, fftw31dm
  PUBLIC :: fftw3_destroy_plan, fftw3_create_plan_1dm, fftw3_create_plan_3d








CONTAINS

! *****************************************************************************
!> \brief ...
!> \param wisdom_file ...
!> \param ionode ...
! *****************************************************************************
SUBROUTINE fftw3_do_cleanup(wisdom_file,ionode)

    CHARACTER(LEN=*), INTENT(IN)             :: wisdom_file
    LOGICAL                                  :: ionode


# 66 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(wisdom_file))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ionode))==-1) EXIT ;  END DO ; ENDIF


END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param wisdom_file ...
! *****************************************************************************
SUBROUTINE fftw3_do_init(wisdom_file)

    CHARACTER(LEN=*), INTENT(IN)             :: wisdom_file

# 146 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(wisdom_file))==-1) EXIT ;  END DO ; ENDIF


END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param DATA ...
!> \param max_length ...
!> \par History
!>      JGH 23-Jan-2006 : initial version
!>      Adapted for new interface
!>      IAB 09-Jan-2009 : Modified to cache plans in fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!>      IAB 09-Oct-2009 : Added OpenMP directives to 1D FFT, and planning routines
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!>      IAB 11-Sep-2012 : OpenMP parallel 3D FFT (Ruyman Reyes, PRACE)
!> \author JGH
! *****************************************************************************
SUBROUTINE fftw3_get_lengths ( DATA, max_length )


    INTEGER, DIMENSION(*)                    :: DATA
    INTEGER, INTENT(INOUT)                   :: max_length

    INTEGER :: h, i, j, k, m, maxn, maxn_elevens, maxn_fives, maxn_sevens, &
      maxn_thirteens, maxn_threes, maxn_twos, ndata, nmax, number
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: dlocal, idx

!------------------------------------------------------------------------------
! compute ndata
!! FFTW can do arbitrary(?) lengths, maybe you want to limit them to some
!!    powers of small prime numbers though...

  maxn_twos = 15
  maxn_threes = 3
  maxn_fives = 2
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = 37748736

  ndata = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
              END DO
           END DO
        END DO
     END DO
  END DO

  ALLOCATE ( dlocal ( ndata ), idx ( ndata ) )

  ndata = 0
  dlocal ( : ) = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
                 dlocal ( ndata ) = number
              END DO
           END DO
        END DO
     END DO
  END DO

  CALL sortint ( dlocal, ndata, idx )
  ndata = MIN ( ndata, max_length )
  DATA(1:ndata) = dlocal(1:ndata)
  max_length = ndata

  DEALLOCATE ( dlocal, idx )

END SUBROUTINE fftw3_get_lengths


! *****************************************************************************
!> \brief ...
!> \param iarr ...
!> \param n ...
!> \param index ...
! *****************************************************************************
SUBROUTINE sortint ( iarr, n, index )

    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(INOUT)                   :: iarr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: a, i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l, temp

!------------------------------------------------------------------------------

  DO i = 1, n
     INDEX(i) = i
  END DO
  jstack = 0
  l = 1
  ir = n
  DO WHILE(.TRUE.)
  IF (ir-l<m) THEN
     DO j = l + 1, ir
        a = iarr(j)
        ib = INDEX(j)
        DO i = j - 1, 0, -1
           IF(i==0) EXIT
           IF (iarr(i)<=a) EXIT
           iarr(i+1) = iarr(i)
           INDEX(i+1) = INDEX(i)
        END DO
        iarr(i+1) = a
        INDEX(i+1) = ib
     END DO
     IF (jstack==0) RETURN
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack - 2
  ELSE
     k = (l+ir)/2
     temp = iarr(k)
     iarr(k) = iarr(l+1)
     iarr(l+1) = temp
     itemp = INDEX(k)
     INDEX(k) = INDEX(l+1)
     INDEX(l+1) = itemp
     IF (iarr(l+1)>iarr(ir)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l)>iarr(ir)) THEN
        temp = iarr(l)
        iarr(l) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l)
        INDEX(l) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l+1)>iarr(l)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(l)
        iarr(l) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(l)
        INDEX(l) = itemp
     END IF
     i = l + 1
     j = ir
     a = iarr(l)
     ib = INDEX(l)
     DO WHILE(.TRUE.)
     i = i + 1
     DO WHILE(iarr(i)<a)
       i = i + 1
     ENDDO
     j = j - 1
     DO WHILE(iarr(j)>a)
       j = j - 1
     ENDDO
     IF (j<i) EXIT
     temp = iarr(i)
     iarr(i) = iarr(j)
     iarr(j) = temp
     itemp = INDEX(i)
     INDEX(i) = INDEX(j)
     INDEX(j) = itemp
     ENDDO
     iarr(l) = iarr(j)
     iarr(j) = a
     INDEX(l) = INDEX(j)
     INDEX(j) = ib
     jstack = jstack + 2
     IF (jstack>nstack) CALL cp__b("pw/fft/fftw3_lib.F",346," Nstack too small in sortr")
     IF (ir-i+1>=j-l) THEN
        istack(jstack) = ir
        istack(jstack-1) = i
        ir = j - 1
     ELSE
        istack(jstack) = j - 1
        istack(jstack-1) = l
        l = i
     END IF
  END IF

  ENDDO

END SUBROUTINE sortint

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param fft_rank ...
!> \param dim_n ...
!> \param dim_istride ...
!> \param dim_ostride ...
!> \param hm_rank ...
!> \param hm_n ...
!> \param hm_istride ...
!> \param hm_ostride ...
!> \param zin ...
!> \param zout ...
!> \param fft_direction ...
!> \param fftw_plan_type ...
!> \param valid ...
! *****************************************************************************
SUBROUTINE fftw3_create_guru_plan(plan, fft_rank, dim_n, &
                                  dim_istride, dim_ostride, hm_rank, &
                                  hm_n, hm_istride, hm_ostride, &
                                  zin, zout, fft_direction, fftw_plan_type, &
                                  valid)


  IMPLICIT NONE

  INTEGER(KIND=integer8_kind), INTENT ( INOUT )      :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zout
  INTEGER, INTENT(IN) :: dim_n(2), dim_istride(2), dim_ostride(2), &
             hm_n(2), hm_istride(2), hm_ostride(2), fft_rank, &
             fft_direction, fftw_plan_type, hm_rank
  LOGICAL, INTENT(OUT)                               :: valid

# 413 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(fft_rank))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(dim_n))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(dim_istride))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(dim_ostride))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hm_rank))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hm_n))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hm_istride))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(hm_ostride))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(fft_direction))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(fftw_plan_type))==-1) EXIT ;  END DO ; ENDIF
  !MARK_USED does not work with assumed size arguments
  IF(.FALSE.)THEN; DO; IF(ABS(zin(1))>ABS(zout(1)))EXIT; ENDDO; ENDIF
  valid = .FALSE.



END SUBROUTINE

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \retval is_mkl ...
! *****************************************************************************
FUNCTION fftw3_is_mkl_wrapper() RESULT(is_mkl)

  IMPLICIT NONE

  LOGICAL :: is_mkl
# 477 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
  is_mkl = .FALSE.


END FUNCTION

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param nrows ...
!> \param nt ...
!> \param rows_per_thread ...
!> \param rows_per_thread_r ...
!> \param th_planA ...
!> \param th_planB ...
! *****************************************************************************
SUBROUTINE fftw3_compute_rows_per_th(nrows,nt,rows_per_thread,rows_per_thread_r,&
                                     th_planA, th_planB)

    INTEGER, INTENT(IN)                      :: nrows, nt
    INTEGER, INTENT(OUT)                     :: rows_per_thread, &
                                                rows_per_thread_r, th_planA, &
                                                th_planB

         IF (MOD(nrows,nt) .EQ. 0) THEN
            rows_per_thread = nrows/nt
            rows_per_thread_r = 0
            th_planA = nt 
            th_planB = 0
         ELSE
            rows_per_thread   = nrows/nt + 1
            rows_per_thread_r = nrows/nt
            th_planA = MOD(nrows,nt)
            th_planB = nt - th_planA
         ENDIF

END SUBROUTINE


! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param plan_r ...
!> \param dim_n ...
!> \param dim_istride ...
!> \param dim_ostride ...
!> \param hm_n ...
!> \param hm_istride ...
!> \param hm_ostride ...
!> \param input ...
!> \param output ...
!> \param fft_direction ...
!> \param fftw_plan_type ...
!> \param rows_per_th ...
!> \param rows_per_th_r ...
! *****************************************************************************
SUBROUTINE fftw3_create_3d_plans(plan, plan_r, dim_n, dim_istride, dim_ostride, &
                                 hm_n, hm_istride, hm_ostride, &
                                 input, output, &
                                 fft_direction, fftw_plan_type, rows_per_th, & 
                                 rows_per_th_r) 


    INTEGER(KIND=integer8_kind), &
      INTENT(INOUT)                          :: plan, plan_r
    INTEGER, INTENT(INOUT)                   :: dim_n(2), dim_istride(2), &
                                                dim_ostride(2), hm_n(2), &
                                                hm_istride(2), hm_ostride(2)
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: input, output
    INTEGER, INTENT(INOUT)                   :: fft_direction, fftw_plan_type
    INTEGER, INTENT(IN)                      :: rows_per_th, rows_per_th_r

    LOGICAL                                  :: valid

! First plans will have an additional row

    hm_n(2) = rows_per_th
    CALL fftw3_create_guru_plan(plan,1, &
                                 dim_n,dim_istride,dim_ostride, &
                                 2,hm_n,hm_istride,hm_ostride, &
                                 input, output, &
                                 fft_direction, fftw_plan_type, valid)

    IF (.NOT. valid) THEN
         CALL cp__b("pw/fft/fftw3_lib.F",564,"fftw3_create_plan")
    ENDIF

    !!!! Remainder
    hm_n(2) = rows_per_th_r
    CALL fftw3_create_guru_plan(plan_r,1, &
                                 dim_n,dim_istride,dim_ostride, &
                                 2,hm_n,hm_istride,hm_ostride, &
                                 input, output, &
                                 fft_direction, fftw_plan_type, valid)
    IF (.NOT. valid) THEN
         CALL cp__b("pw/fft/fftw3_lib.F",575,"fftw3_create_plan (remaining)")
    ENDIF


END SUBROUTINE

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param zin ...
!> \param zout ...
!> \param plan_style ...
! *****************************************************************************
SUBROUTINE fftw3_create_plan_3d(plan, zin, zout, plan_style)

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)      :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)      :: zout
  INTEGER                                            :: plan_style
# 755 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan_style))==-1) EXIT ;  END DO ; ENDIF
   !MARK_USED does not work with assumed size arguments
   IF(.FALSE.)THEN; DO; IF(ABS(zin(1))>ABS(zout(1)))EXIT; ENDDO; ENDIF


END SUBROUTINE fftw3_create_plan_3d


! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param plan_r ...
!> \param split_dim ...
!> \param nt ...
!> \param tid ...
!> \param input ...
!> \param istride ...
!> \param output ...
!> \param ostride ...
! *****************************************************************************
SUBROUTINE fftw3_workshare_execute_dft(plan, plan_r, split_dim, nt, tid, &
                                       input, istride, output, ostride)

  INTEGER, INTENT(IN)                           :: split_dim,nt, tid
  INTEGER, INTENT(IN)                           :: istride, ostride
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT) :: input, output
    INTEGER (KIND=integer8_kind)                :: plan, plan_r
# 822 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan_r))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(split_dim))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nt))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(tid))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(istride))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ostride))==-1) EXIT ;  END DO ; ENDIF
   !MARK_USED does not work with assumed size arguments
   IF(.FALSE.)THEN; DO; IF(ABS(input(1))>ABS(output(1)))EXIT; ENDDO; ENDIF


END SUBROUTINE


! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param scale ...
!> \param zin ...
!> \param zout ...
!> \param stat ...
! *****************************************************************************
SUBROUTINE fftw33d ( plan, scale, zin, zout, stat )

  TYPE(fft_plan_type), INTENT(IN)                      :: plan
  REAL(KIND=dp), INTENT(IN)                            :: scale
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET:: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET:: zout
  INTEGER, INTENT(OUT)                                 :: stat
# 924 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
  !MARK_USED does not work with assumed size arguments
  IF(.FALSE.)THEN; DO; IF(ABS(zin(1))>ABS(zout(1)))EXIT; ENDDO; ENDIF
  stat = 0



END SUBROUTINE fftw33d

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param zin ...
!> \param zout ...
!> \param plan_style ...
! *****************************************************************************
SUBROUTINE fftw3_create_plan_1dm(plan, zin, zout, plan_style)

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zout
  INTEGER, INTENT(IN)                                :: plan_style
# 1028 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan_style))==-1) EXIT ;  END DO ; ENDIF
   !MARK_USED does not work with assumed size arguments
   IF(.FALSE.)THEN; DO; IF(ABS(zin(1))>ABS(zout(1)))EXIT; ENDDO; ENDIF


END SUBROUTINE fftw3_create_plan_1dm

! *****************************************************************************
!> \brief ...
!> \param plan ...
! *****************************************************************************
SUBROUTINE fftw3_destroy_plan ( plan )

  TYPE(fft_plan_type), INTENT (INOUT)   :: plan

# 1063 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF


END SUBROUTINE fftw3_destroy_plan

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param zin ...
!> \param zout ...
!> \param scale ...
!> \param stat ...
! *****************************************************************************
SUBROUTINE fftw31dm ( plan, zin, zout, scale, stat )
    TYPE(fft_plan_type), INTENT(IN)          :: plan
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(IN), TARGET                     :: zin
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT), TARGET                  :: zout
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, INTENT(OUT)                     :: stat

    COMPLEX(KIND=dp), POINTER                :: zin_ptr, zout_ptr, zscal_ptr
    INTEGER                                  :: in_offset, my_id, num_rows, &
                                                out_offset, scal_offset
    INTEGER(KIND=C_INTPTR_T)                 :: fftw_plan

!------------------------------------------------------------------------------

my_id = 0
num_rows = plan%m

!$omp parallel default(none), &
!$omp          private(my_id,num_rows,zin_ptr,zout_ptr,zscal_ptr,in_offset,out_offset,scal_offset,fftw_plan), &
!$omp          shared(zin,zout), &
!$omp          shared(plan,scale,stat)
!$ my_id = omp_get_thread_num()

!$ if (my_id < plan%num_threads_needed) then

fftw_plan = plan%fftw_plan

in_offset = 1
out_offset = 1
scal_offset = 1

!$ in_offset = 1 + plan%num_rows * my_id * plan%n
!$ out_offset = 1 + plan%num_rows * my_id * plan%n
!$ IF ( plan%fsign == +1 .AND. plan%trans ) THEN
!$  in_offset = 1 + plan%num_rows*my_id
!$ ELSEIF ( plan%fsign == -1 .AND. plan%trans ) THEN
!$  out_offset = 1 + plan%num_rows*my_id
!$ ENDIF
!$ scal_offset = 1 + plan%n*plan%num_rows*my_id
!$ IF ( plan%need_alt_plan .AND. my_id .EQ. plan%num_threads_needed - 1 ) THEN
!$   num_rows = plan%alt_num_rows
!$   fftw_plan = plan%alt_fftw_plan
!$ ELSE
!$   num_rows = plan%num_rows
!$ ENDIF

zin_ptr => zin(in_offset)
zout_ptr => zout(out_offset)
zscal_ptr => zout(scal_offset)

# 1141 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/fft/fftw3_lib.F"
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
  !MARK_USED does not work with assumed size arguments
  IF(.FALSE.)THEN; DO; IF(ABS(zin(1))>ABS(zout(1)))EXIT; ENDDO; ENDIF
  stat=0

!$ else
!$ end if



!$omp end parallel

END SUBROUTINE fftw31dm

!     Copyright (c) 2003, 2006 Matteo Frigo
!     Copyright (c) 2003, 2006 Massachusetts Institute of Technology
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     This is an example implementation of Fortran wisdom export/import
!     to/from a Fortran unit (file), exploiting the generic
!     dfftw_export_wisdom/dfftw_import_wisdom functions.
!
!     We cannot compile this file into the FFTW library itself, lest all
!     FFTW-calling programs be required to link to the Fortran I/O
!     libraries.
!
!     adapted to become more standard Fortran 90 [2007-10] Joost VandeVondele
!     and added some namespacing
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! *****************************************************************************
!> \brief ...
!> \param c ...
!> \param iunit ...
! *****************************************************************************
      SUBROUTINE fftw_write_char(c, iunit) BIND(C, name="fftw_write_char")
    CHARACTER(KIND=C_CHAR)                   :: c
    INTEGER(KIND=C_INT)                      :: iunit

         WRITE(iunit,'(a)',ADVANCE="NO") c
      END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param iunit ...
! *****************************************************************************
      SUBROUTINE fftw_export_wisdom_to_file(iunit)
    INTEGER                                  :: iunit




   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(iunit))==-1) EXIT ;  END DO ; ENDIF

      END SUBROUTINE

!     Fortran 77 does not have any portable way to read an arbitrary
!     file one character at a time [needs to wait for stream IO of F2003].
!     The best alternative seems to be to
!     read a whole line into a buffer, since for fftw-exported wisdom we
!     can bound the line length.  (If the file contains longer lines,
!     then the lines will be truncated and the wisdom import should
!     simply fail.)  Ugh (and not thread safe).

! *****************************************************************************
!> \brief ...
!> \param ic ...
!> \param iunit ...
! *****************************************************************************
      SUBROUTINE fftw_read_char(ic, iunit) BIND(C, name="fftw_read_char")
    INTEGER(KIND=C_INT)                      :: ic, iunit

    CHARACTER(LEN=256)                       :: buf

         SAVE buf
         INTEGER ibuf
         DATA ibuf/257/
         SAVE ibuf
         IF (ibuf .LT. 257) THEN
            ic = ICHAR(buf(ibuf:ibuf))
            ibuf = ibuf + 1
            RETURN
         ENDIF
         READ(iunit,123,END=666) buf
         ic = ICHAR(buf(1:1))
         ibuf = 2
         RETURN
 666     ic = -1
         ibuf = 257
 123     FORMAT(a256)
      END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param isuccess ...
!> \param iunit ...
! *****************************************************************************
      SUBROUTINE fftw_import_wisdom_from_file(isuccess, iunit)
    INTEGER                                  :: isuccess, iunit

         isuccess=0



   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(iunit))==-1) EXIT ;  END DO ; ENDIF

      END SUBROUTINE

END MODULE

