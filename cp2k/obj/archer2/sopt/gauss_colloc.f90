# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!>   \brief
!>     Routines to efficently collocate and integrate gaussians on a grid
!>     These use most of Joost's tricks and a couple more...
!>     result is *speed* and genericity
!>   \author Fawzi Mohamed, 2007
!>   \notes original available with BSD style license
! *****************************************************************************
MODULE gauss_colloc
  USE d3_poly,                         ONLY: &
       grad_size3, poly_affine_t3, poly_affine_t3t, poly_p_eval2b, &
       poly_p_eval3b, poly_padd_uneval2b, poly_padd_uneval3b, poly_size1, &
       poly_size2, poly_size3
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE lgrid_types,                     ONLY: lgrid_type

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

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
# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

IMPLICIT NONE
PRIVATE

  PUBLIC :: collocGauss,&
            integrateGaussFull







  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gauss_colloc'

    REAL(dp), PARAMETER :: small=TINY(1.0_dp)

  ! keep prettify happy (does not see the include)
  INTEGER(KIND=int_8), PARAMETER           :: unused_import_of_int_8=1


CONTAINS


! *****************************************************************************
!> \brief collocate a periodically repeated gaussian on a non orthormbic grid
!>
!> this routine has been tested and works well with cells with
!> det(h)/sqrt(tr(dot(h^T,h)))>0.2 (2 angles bigger than 24 deg or one angle
!> bigger than 11 deg).
!> Because of its numerics it might fail badly (infinity or NaN) with
!> with more deformed cells. Avoiding this would be bossible only using
!> IEEE numerics controls, which would also make everything slower and
!> less supported.
!> Still the actual numeric has been carefully tuned, and in normal cases
!> and most anormal it should work.
!> With det(h)/sqrt(tr(dot(h^T,h)))>0.2 I could not find any failure.
!>
!> \param h cell matrix
!> \param h_inv inverse of the cell matrix
!> \param grid the grid
!> \param poly polynomial (d3_poly format)
!> \param alphai exponential coeff
!> \param posi position of the gaussian
!> \param max_r2 maximum radius of collocation squared
!> \param periodic array of 0 or 1 that says which dimensions have pbc (1=pbc)
!> \param gdim dimension of the grid (grid might be a subset)
!> \param local_bounds local bounds of the grid piece that is kept locally
!>   (i.e. of grid) the global grid is assumed to atart at 0,0,0
!> \param local_shift start indexes of the local slice (i.e. of grid)
!> \param poly_shift position of posi in the polynomial reference system.
!>  Set it to posi to use the global reference system.
!> \param scale a global scale factor
!> \param lgrid ...
! *****************************************************************************
SUBROUTINE collocGauss(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,lgrid)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(inout)    :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(lgrid_type), INTENT(inout), &
      OPTIONAL                               :: lgrid

    CHARACTER(len=*), PARAMETER :: routineN = 'collocGauss', &
      routineP = moduleN//':'//routineN

    CALL colloc_int_body(h,h_inv,grid,poly,alphai,posi,max_r2,&
                         periodic,gdim,local_bounds,local_shift,&
                         poly_shift,scale,lgrid,integrate=.FALSE.)

END SUBROUTINE


! *****************************************************************************
!> \brief integrates a gaussian times any polynomial up to a give order.
!>
!> Most things are the same as for collocGauss (see its comments).
!> Returns the integrals of all the monomials in d3 format into poly
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
! *****************************************************************************
SUBROUTINE integrateGaussFull(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(inout)    :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'integrateGaussFull', &
      routineP = moduleN//':'//routineN

    CALL colloc_int_body(h,h_inv,grid,poly,alphai,posi,max_r2,&
                         periodic,gdim,local_bounds,local_shift,&
                         poly_shift,scale,integrate=.TRUE.)


END SUBROUTINE



! *****************************************************************************
!> \brief Common code for collocGauss() and integrateGaussFull()
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param lgrid ...
!> \param integrate ...
! *****************************************************************************
SUBROUTINE colloc_int_body(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,lgrid,integrate)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(inout)    :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(lgrid_type), INTENT(inout), &
      OPTIONAL                               :: lgrid
    LOGICAL                                  :: integrate

    CHARACTER(len=*), PARAMETER :: routineN = 'colloc_int_body', &
        routineP=moduleN//':'//routineN

    INTEGER, DIMENSION(0:2), PARAMETER       :: permut = (/ 2,1,0 /)
    INTEGER :: grad, i, i0, ii, iiShift, iiShift2, iistart, iistart2, ij, &
      ijShift, iJump, ik, ikShift, ikShift2, ikstart, ikstart2, iend,iend2,&
      imax, imax1, imin, imin1, istart, istart2, j, jend, jJump, jmax, jmax1, jmin, &
      jmin1, jstart, k, kend, kend2, kgrad, kJump, kmax, kmax1, kmin, kmin1, &
      kstart, kstart2, max_j, size_jk,size_k,size_ijk,ig,ithread,nthread
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: k_bounds
    INTEGER, DIMENSION(0:2)                  :: cellShift, l_shift, l_ub, &
                                                ndim, period, shiftPos,ldim2
    INTEGER, DIMENSION(2, 0:2)               :: l_bounds
    LOGICAL                                  :: has_overlap
    REAL(dp) :: cci0, cci1, cci2, ccj0, ccj0_i0, ccj0_i1, ccj0_i2, ccj1, &
      ccj1_i0, ccj1_i1, ccj2, cck0, cck0_0, cck0_0_p, cck0_i, cck0_i2, &
      cck0_ij, cck0_j, cck0_j2, cck0_j_p, cck1, cck1_0, cck1_0_p, cck1_i, &
      cck1_j, cck2, delta_i, delta_j, delta_k, g_scale, i_coeffn_i, icoeff0, &
      ii_coeff0, ii_coeff2, ii_coeff2_jump, ii_coeffn, ii_coeffn_jump, &
      ij_coeff0, ij_coeff0_jump, ij_coeff1, ik_coeff0, ik_coeff0_jump, &
      ik_coeff1, j_coeffn_i, j_coeffn_j, jcoeff0, jj_coeff0, jj_coeff2, &
      jj_coeffn, jk_coeff0, jk_coeff1, k_coeffn_i, k_coeffn_j, k_coeffn_k, &
      kcoeff0, kk_coeff0, kk_coeff2, kk_coeffn, m(0:2,0:2), maxr2, p_kk, &
      r_0
    REAL(dp) :: res_0, res_i, res_j, res_k, scaled_h(0:2,0:2), sqDi, sqDj, &
      sqDk
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: poly_ijk, poly_jk,xi
    REAL(dp), DIMENSION(0:2)                 :: l, normD, p_shift, resPos, &
                                                resPosReal, riPos, rpos, wrPos

    REAL(dp) :: det, gval
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: k_vals
    INTEGER, PARAMETER :: npoly=1
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: poly_k
    REAL(dp) :: p_v




    ig = 0
    ithread = 0
!$  ithread = omp_get_thread_num()

    nthread = 1
!$  nthread = omp_get_num_threads()


    IF (integrate) THEN
       poly=0.0_dp
    ELSE
       IF (ALL(poly==0.0_dp)) RETURN
    ENDIF

    IF (PRESENT(poly_shift)) THEN
        p_shift=poly_shift
    ELSE
        p_shift=0.0_dp
    END IF

    ldim2(permut(0))=SIZE(grid,1)
    ldim2(permut(1))=SIZE(grid,2)
    ldim2(permut(2))=SIZE(grid,3)

    IF (PRESENT(gdim)) THEN
        DO i=0,2
            ndim(permut(i))=gdim(i)
        END DO
    ELSE
        ndim=ldim2
    END IF
    g_scale=1.0_dp
    IF (PRESENT(scale)) g_scale=scale

    IF (integrate) THEN
       det=(h(0,0)*(h(1,1)*h(2,2)-h(1,2)*h(2,1))&
           -h(1,0)*(h(0,1)*h(2,2)-h(0,2)*h(2,1))&
           +h(2,0)*(h(0,1)*h(1,2)-h(0,2)*h(1,1)))
       g_scale=g_scale*ABS(det)/REAL(INT(ndim(0),KIND=int_8)*INT(ndim(1),KIND=int_8)*INT(ndim(2),KIND=int_8),dp)
    ENDIF

    IF (PRESENT(local_bounds)) THEN
        DO i=0,2
            l_bounds(:,permut(i))=local_bounds(:,i)
        END DO
    ELSE
        l_bounds(1,:)=0
        l_bounds(2,:)=ndim-1
    END IF
    IF (PRESENT(local_shift)) THEN
        DO i=0,2
            l_shift(permut(i))=local_shift(i)
        END DO
    ELSE
        l_shift=0 ! l_bounds(1,:)
    END IF
    l_ub=l_bounds(2,:)-l_bounds(1,:)+l_shift
    DO i=0,2
        IF(.NOT.(l_ub(i)<ldim2(i)))CALL cp__a("gauss_colloc.F",304)
    END DO

    DO i=0,2
        period(permut(i))=periodic(i)
    END DO
    IF(.NOT.(ALL(l_bounds(2,:)<ndim.or.period(:)==1)))CALL cp__a("gauss_colloc.F",310)
    IF(.NOT.(ALL(l_bounds(1,:)>=0 .or.period(:)==1)))CALL cp__a("gauss_colloc.F",311)
    IF(.NOT.(ALL(l_bounds(2,:)-l_bounds(1,:)<ndim)))CALL cp__a("gauss_colloc.F",312)
    rPos=0.0_dp
    DO j=0,2
        DO i=0,2
            rPos(permut(i))=rPos(permut(i))+h_inv(i,j)*posi(j)
        END DO
    END DO
    cellShift=FLOOR(rPos)*period
    wrPos=rPos-REAL(cellShift,dp)
    riPos=wrPos*ndim
    shiftPos=FLOOR(riPos+0.5_dp)
    resPos=riPos-shiftPos
    normD=1.0_dp/REAL(ndim,dp)
    scaled_h=0.0_dp
    DO j=0,2
        DO i=0,2
            scaled_h(i,permut(j))=h(i,j)*normD(permut(j))
        END DO
    END DO
    resPosReal=0.0_dp
    DO j=0,2
        DO i=0,2
            resPosReal(i)=resPosReal(i)+h(i,j)*(wrPos(permut(j))-normD(permut(j))*REAL(shiftPos(permut(j)),dp))
        END DO
    END DO

    !maxr2=0.0_dp
    !DO j=0,2
    !    DO i=0,2
    !        ! guarantee at least the nearest points (this increases the sphere, increase just the box?)
    !        maxr2=maxr2+(scaled_h(i,j))**2 
    !    END DO
    !END DO
    maxr2=max_r2 !MAX(max_r2,maxr2)

    ! build up quadratic form (ellipsoid)
    m=0.0_dp
    DO j=0,2
        DO i=0,2
            DO k=0,2
                m(i,j)=m(i,j)+scaled_h(k,i)*scaled_h(k,j)
            END DO
        END DO
    END DO

    l=0.0_dp
    DO j=0,2
        DO i=0,2
            l(j)=l(j)-2.0*resPos(i)*m(i,j)
        END DO
    END DO

    r_0=0.0_dp
    DO i=0,2
        r_0=r_0-0.5*resPos(i)*l(i)
    END DO

    ! calc i boundaries
    cci2 = (m(2,2) * m(0,0) * m(1,1) - m(2,2) * m(0,1) ** 2 - m(1,1) * m(0,2) ** 2 &
        + 2.0_dp * m(0,2) * m(0,1) * m(1,2) - m(0,0) * m(1,2) ** 2) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci1 = -(-m(2,2) * l(0) * m(1,1) + m(2,2) * m(0,1) * l(1) + l(2) * m(0,2) * m(1,1) &
        + l(0) * m(1,2) ** 2 - l(2) * m(0,1) * m(1,2) - m(1,2) * l(1) * m(0,2)) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci0 = -((-4.0_dp * m(2,2) * r_0 * m(1,1) + m(2,2) * l(1) ** 2 + l(2) ** 2 * m(1,1) &
        - 2.0_dp * l(1) * m(1,2) * l(2) + 4.0_dp * r_0 * m(1,2) ** 2) &
        / (m(2,2) * m(1,1) - m(1,2) ** 2)) / 4.0_dp-maxr2
    delta_i=cci1*cci1-4.0_dp*cci2*cci0
    IF (delta_i<=0) RETURN
    sqDi=SQRT(delta_i)
    imin=CEILING((-cci1-sqDi)/(2.0_dp*cci2))
    imax=FLOOR((-cci1+sqDi)/(2.0_dp*cci2))
    !! box early return

    IF (period(0)==1) THEN
        has_overlap=imax-imin+1>ndim(0).OR.(l_bounds(1,0)==0.and.l_bounds(2,0)==ndim(0)-1)
        IF (.not.has_overlap) THEN
            imin1=MODULO(imin+shiftPos(0),ndim(0))
            imax1=imin1+imax-imin+1
            IF (imin1<l_bounds(1,0)) THEN
                has_overlap=imax1>=l_bounds(1,0)
            ELSE
                has_overlap=imin1<=l_bounds(2,0).OR.(imax1>=ndim(0).and.l_bounds(1,0)<=imax1+ndim(0))
            END IF
            IF (.not.has_overlap) RETURN
        END IF
    ELSE
        IF (imax+shiftPos(0)<l_bounds(1,0).or.imin+shiftPos(0)>l_bounds(2,0)) RETURN
    END IF

    ! j box bounds
    has_overlap=l_bounds(1,1)==0.and.l_bounds(2,1)==ndim(1)-1
    IF (.not.has_overlap) THEN
        ccj2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
            + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2)
        ccj1 = -(-m(0,0) * l(1) * m(2,2) + m(0,0) * l(2) * m(1,2) + l(0) * m(0,1) * m(2,2) &
            - m(0,2) * l(2) * m(0,1) - l(0) * m(0,2) * m(1,2) + l(1) * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2)
        ccj0 = (4.0_dp * m(0,0) * m(2,2) * r_0 - m(0,0) * l(2) ** 2 - m(2,2) * l(0) ** 2 &
            + 2.0_dp * m(0,2) * l(2) * l(0) - 4.0_dp * r_0 * m(0,2) ** 2) &
            / (m(0,0) * m(2,2) - m(0,2) ** 2) / 4.0_dp-maxr2
        delta_j=ccj1*ccj1-4.0_dp*ccj2*ccj0
        IF (delta_j<=0) RETURN
        sqDj=SQRT(delta_j)
        jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
        jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
        IF (period(1)==1) THEN
            IF (jmax-jmin+1<ndim(1)) THEN
                jmin1=MODULO(jmin+shiftPos(1),ndim(1))
                jmax1=jmin1+jmax-jmin+1
                IF (jmin1<l_bounds(1,1)) THEN
                    has_overlap=jmax1>=l_bounds(1,1)
                ELSE
                    has_overlap=jmin1<=l_bounds(2,1).OR.(jmax1>=ndim(1).AND.(l_bounds(1,1)<=jmax1-ndim(1)))
                END IF
                IF (.not.has_overlap) RETURN
            END IF
        ELSE
            IF (jmax+shiftPos(1)<l_bounds(1,1).or.jmin+shiftPos(1)>l_bounds(2,1)) RETURN
        END IF
    END IF

    ! k box bounds
    has_overlap=l_bounds(1,2)==0.and.l_bounds(2,2)==ndim(2)-1
    IF (.not.has_overlap) THEN
        cck2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
            + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) / (m(0,0) * m(1,1) - m(0,1) ** 2)
        cck1 = (m(0,0) * l(2) * m(1,1) - m(0,0) * m(1,2) * l(1) - m(0,2) * l(0) * m(1,1) &
            + l(0) * m(0,1) * m(1,2) - l(2) * m(0,1) ** 2 + m(0,1) * l(1) * m(0,2)) / (m(0,0) * m(1,1) - m(0,1) ** 2)
        cck0 = (4.0_dp * m(0,0) * m(1,1) * r_0 - m(0,0) * l(1) ** 2 - m(1,1) * l(0) ** 2 &
            + 2.0_dp * m(0,1) * l(1) * l(0) - 4.0_dp * r_0 * m(0,1) ** 2) &
            / (m(0,0) * m(1,1) - m(0,1) ** 2) / 4.0_dp-maxr2
        delta_k=cck1*cck1-4.0_dp*cck2*cck0
        IF (delta_k<=0) RETURN
        sqDk=SQRT(delta_k)
        kmin=CEILING((-cck1-sqDk)/(2.0_dp*cck2))
        kmax=FLOOR((-cck1+sqDk)/(2.0_dp*cck2))

        IF (period(2)==1) THEN
            IF (kmax-kmin+1<ndim(2)) THEN
                kmin1=MODULO(kmin+shiftPos(2),ndim(2))
                kmax1=kmin1+kmax-kmin+1
                IF (kmin1<l_bounds(1,2)) THEN
                    has_overlap=kmax1>=l_bounds(1,2)
                ELSE
                    has_overlap=kmin1<=l_bounds(2,2).OR.&
                        (kmax1>=ndim(2).AND.(l_bounds(1,2)<=MODULO(kmax1,ndim(2))))
                END IF
                IF (.not.has_overlap) RETURN
            END IF
        ELSE
            IF (kmax+shiftPos(2)<l_bounds(1,2).or.kmin+shiftPos(2)>l_bounds(2,2)) RETURN
        END IF
    END IF

    ! k bounds (cache a la cube_info, or inversely integrate in the collocate loop?)
    ccj2   = (m(2,2) * m(1,1) - m(1,2) ** 2) / m(2,2)
    ccj1_i1=(2 * m(2,2) * m(0,1) - 2 * m(0,2) * m(1,2)) / m(2,2)
    ccj1_i0=(-l(2) * m(1,2) + m(2,2) * l(1)) / m(2,2)
    ccj0_i2=(m(2,2) * m(0,0)-m(0,2) ** 2) / m(2,2)
    ccj0_i1=( m(2,2) * l(0) - m(0,2) * l(2) ) / m(2,2)
    ccj0_i0=(m(2,2) * r_0 - 0.25 * l(2) ** 2) / m(2,2) - maxr2
    cck2   = m(2,2)
    cck1_i = 2 * m(0,2)
    cck1_j = 2 * m(1,2)
    cck1_0 = l(2)
    cck0_i2 = m(0,0)
    cck0_ij = 2 * m(0,1)
    cck0_i = l(0)
    cck0_j2 = m(1,1)
    cck0_j = l(1)
    cck0_0 = r_0 - maxr2

    ! find maximum number of j
    max_j=0
    DO i0=0,1
        i=(imin+imax)/2+i0
        ccj1 = ccj1_i1 * i +ccj1_i0
        ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
        delta_j=ccj1*ccj1-4*ccj2*ccj0
        IF (delta_j>=0) THEN
            sqDj=SQRT(delta_j)
            max_j=MAX(max_j,FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2)) &
                        -CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))+1)
        END IF
    END DO
    max_j=max_j+1 ! just to be sure...
    IF (period(1)==0) max_j=MIN(max_j,l_bounds(2,1)-l_bounds(1,1)+1)

    IF (period(0)==0) THEN
        imin=MAX(l_bounds(1,0)-shiftPos(0),imin)
        imax=MIN(l_bounds(2,0)-shiftPos(0),imax)
    END IF

    ! k bounds (cache a la cube_info?)
    has_overlap=.FALSE.
    ALLOCATE(k_bounds(0:1,0:max_j-1,0:MAX(0,imax-imin+1)))
    ! k_bounds=0
    istart=imin
    iiShift=shiftPos(0)-l_bounds(2,0)+istart
    IF (iiShift>0) iiShift=iiShift+ndim(0)-1
    iiShift=(iiShift/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift=CEILING(REAL(shiftPos(0)+istart-l_bounds(2,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0))
    istart=MAX(iiShift+l_bounds(1,0),istart)
    iend=MIN(iiShift+l_bounds(2,0),imax)
    iJump=ndim(0)-l_bounds(2,0)+l_bounds(1,0)-1
    jJump=ndim(1)-l_bounds(2,1)+l_bounds(1,1)-1
    DO
        DO i=istart,iend
            ccj1 = ccj1_i1 * i +ccj1_i0
            ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
            delta_j=ccj1*ccj1-4*ccj2*ccj0
            IF (delta_j<0) CONTINUE
            sqDj=SQRT(delta_j)
            jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
            jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
            cck0_0_p=cck0_0+(cck0_i2*i+cck0_i)*i
            cck0_j_p= cck0_j+cck0_ij*i
            cck1_0_p=cck1_0+cck1_i*i
            IF (period(1)==0) THEN
                jmin=MAX(l_bounds(1,1)-shiftPos(1),jmin)
                jmax=MIN(l_bounds(2,1)-shiftPos(1),jmax)
            END IF
            jstart=jmin
            ijShift=shiftPos(1)+jstart-l_bounds(2,1)
            IF (ijShift>0) ijShift=ijShift+ndim(1)-1
            ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
            ! ijShift=CEILING(REAL(shiftPos(1)+jstart-l_bounds(2,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1)
            jstart=MAX(ijShift+l_bounds(1,1),jstart)
            jend=MIN(ijShift+l_bounds(2,1),jmax)
            DO
                DO j=jstart,jend
                    cck1=cck1_0_p+cck1_j*j
                    cck0=cck0_0_p+(cck0_j_p+cck0_j2*j)*j

                    delta_k=cck1*cck1-4*cck2*cck0
                    IF (delta_k<0) THEN
                        k_bounds(0,j-jmin,i-imin)=0 ! CEILING((-cck1)/(2.0_dp*cck2))
                        k_bounds(1,j-jmin,i-imin)=-1 ! k_bounds(0,j-jmin,i-imin)-1
                    ELSE
                        sqDk=SQRT(delta_k)
                        kmin=CEILING((-cck1-sqDk)/(2.0_dp*cck2))
                        kmax=FLOOR((-cck1+sqDk)/(2.0_dp*cck2))

                        ! ! reduce kmax,kmin
                        ! ! this should be done later if k_bounds are shared by threads with different slices
                        ! ikShift=FLOOR(REAL(shiftPos(2)+kmax-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
                        ! kmax=MIN(kmax,ikShift+l_bounds(2,2))
                        ! ikShift2=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kmin)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
                        ! kmin=MAX(kmin,ikShift2+l_bounds(1,2))

                        k_bounds(0,j-jmin,i-imin)=kmin
                        k_bounds(1,j-jmin,i-imin)=kmax
                        IF (kmax>=kmin) has_overlap=.TRUE.
                    END IF
                END DO
                jstart=jend+jJump+1
                IF (jstart>jmax) EXIT
                jend=MIN(jend+ndim(1),jmax)
            END DO
        END DO
        istart=iend+iJump+1
        IF (istart>imax) EXIT
        iend=MIN(iend+ndim(0),imax)
    END DO
    IF (period(2)==0) THEN
        k_bounds(0,:,:)=MAX(l_bounds(1,2)-shiftPos(2),k_bounds(0,:,:))
        k_bounds(1,:,:)=MIN(l_bounds(2,2)-shiftPos(2),k_bounds(1,:,:))
    END IF

    IF (.not.has_overlap) RETURN

    ! poly x,y,z -> i,j,k
    grad=grad_size3(SIZE(poly)/npoly)
    size_jk=poly_size2(grad)*npoly
    size_k=poly_size1(grad)*npoly
    size_ijk=poly_size3(grad)*npoly
    ALLOCATE(poly_ijk(size_ijk),&
        poly_jk(size_jk),&
        xi(grad+1))

    IF(integrate) THEN
        ALLOCATE(k_vals(0:grad))
    ELSE
        ALLOCATE(poly_k(0:size_k-1))
    ENDIF

    IF (integrate) THEN
       IF(.NOT.(SIZE(poly)==poly_size3(grad)))CALL cp__a("gauss_colloc.F",599)
       poly_ijk=0.0_dp
    ELSE
       CALL poly_affine_t3(poly,scaled_h,-resPosReal+p_shift,poly_ijk,&
           npoly=npoly)
    ENDIF

    ij_coeff0=EXP(-2.0_dp*alphai*m(0,1))
    ik_coeff0=EXP(-2.0_dp*alphai*m(0,2))
    ii_coeff0=EXP(-alphai*m(0,0))
    jk_coeff0=EXP(-2.0_dp*alphai*m(1,2))
    jj_coeff0=EXP(-alphai*m(1,1))
    kk_coeff0=EXP(-alphai*m(2,2))
    jk_coeff1=1.0_dp/jk_coeff0
    ij_coeff1=1.0_dp/ij_coeff0
    ik_coeff1=1.0_dp/ik_coeff0
    ii_coeff2=ii_coeff0*ii_coeff0
    jj_coeff2=jj_coeff0*jj_coeff0
    kk_coeff2=kk_coeff0*kk_coeff0
    icoeff0=EXP(-alphai*l(0))
    jcoeff0=EXP(-alphai*l(1))
    kcoeff0=EXP(-alphai*l(2))
    res_0=EXP(-alphai*r_0)*g_scale

    i_coeffn_i=icoeff0
    j_coeffn_i=jcoeff0
    k_coeffn_i=kcoeff0
    ii_coeffn =i_coeffn_i*ii_coeff0
    res_i=res_0

    iJump=ndim(0)-l_bounds(2,0)+l_bounds(1,0)-1
    istart=MAX(0,imin)
    iiShift=shiftPos(0)-l_bounds(2,0)+istart
    IF (iiShift>0) iiShift=iiShift+ndim(0)-1
    iiShift=(iiShift/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift=CEILING(REAL(shiftPos(0)+istart-l_bounds(2,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0)
    istart=MAX(iiShift+l_bounds(1,0),istart)
    iistart=istart-iiShift-l_bounds(1,0)+l_shift(0)
    istart2=MIN(-1,imax)
    iiShift2=shiftPos(0)+istart2-l_bounds(1,0)
    IF (iiShift2<0) iiShift2=iiShift2-ndim(0)+1
    iiShift2=(iiShift2/ndim(0))*ndim(0)-shiftPos(0)
    !iiShift2=FLOOR(REAL(shiftPos(0)+istart2-l_bounds(1,0))/REAL(ndim(0)))*ndim(0)-shiftPos(0)
    istart2=MIN(iiShift2+l_bounds(2,0),istart2)
    iistart2=istart2-iiShift2-l_bounds(1,0)+l_shift(0)

    IF (iJump/=0.and.(iistart+imax-istart>=ndim(0)+l_shift(0) .OR.&
        iistart2+imin-istart2<=l_ub(0)-ndim(0))) THEN
        ! will wrap
        ij_coeff0_jump =ij_coeff0**(iJump)
        ik_coeff0_jump =ik_coeff0**(iJump)
        ii_coeff2_jump =ii_coeff2**(iJump)
        ii_coeffn_jump =ii_coeff0**((iJump)*(iJump-1))
    ELSE
        ij_coeff0_jump = 1.0_dp
        ik_coeff0_jump = 1.0_dp
        ii_coeff2_jump = 1.0_dp
        ii_coeffn_jump = 1.0_dp
    END IF

    iend=MIN(iiShift+l_bounds(2,0),imax)
    ii=iistart 
    IF (i>0) THEN
        ii_coeffn=i_coeffn_i*ii_coeff0**(2*istart+1)
        j_coeffn_i=jcoeff0*ij_coeff0**istart
        k_coeffn_i=kcoeff0*ik_coeff0**istart
        res_i=res_0*(ii_coeff0**istart*i_coeffn_i)**istart
    END IF
    DO
        DO i=istart,iend
            ! perform j loop
            IF (ABS(res_i)>small) THEN
                CALL j_loop
            END IF
            j_coeffn_i=j_coeffn_i*ij_coeff0
            k_coeffn_i=k_coeffn_i*ik_coeff0
            res_i=res_i*ii_coeffn
            ii_coeffn=ii_coeffn*ii_coeff2
            ii=ii + 1
        END DO
        istart=iend+iJump+1
        IF (istart>imax) EXIT
        iend=MIN(iend+ndim(0),imax)
        ii=l_shift(0) 
        j_coeffn_i=j_coeffn_i*ij_coeff0_jump
        k_coeffn_i=k_coeffn_i*ik_coeff0_jump
        res_i=res_i*ii_coeffn**(iJump)*ii_coeffn_jump
        ii_coeffn=ii_coeffn*ii_coeff2_jump
    END DO

    ! neg i side
    i_coeffn_i=1.0_dp/icoeff0
    j_coeffn_i=jcoeff0
    k_coeffn_i=kcoeff0
    res_i=res_0
    ii_coeffn=i_coeffn_i*ii_coeff0

    iend2=MAX(iiShift2+l_bounds(1,0),imin)
    ii=iistart2 
    IF (istart2<-1) THEN
        ii_coeffn=i_coeffn_i*ii_coeff0**(-(2*istart2+1))
        j_coeffn_i=jcoeff0*ij_coeff0**(istart2+1)
        k_coeffn_i=kcoeff0*ik_coeff0**(istart2+1)
        res_i=res_0*(ii_coeff0**(-istart2-1)*i_coeffn_i)**(-istart2-1)
    END IF
    DO
        DO i=istart2,iend2,-1
            j_coeffn_i=j_coeffn_i*ij_coeff1
            k_coeffn_i=k_coeffn_i*ik_coeff1
            res_i=res_i*ii_coeffn
            ii_coeffn=ii_coeffn*ii_coeff2

            ! perform j loop
            IF (ABS(res_i)>small) THEN
                CALL j_loop
            END IF
            ii=ii-1
        END DO
        istart2=iend2-iJump-1
        IF (istart2<imin) EXIT
        iend2=MAX(iend2-ndim(0),imin)
        ii=l_ub(0) 
        j_coeffn_i=j_coeffn_i/ij_coeff0_jump
        k_coeffn_i=k_coeffn_i/ik_coeff0_jump
        res_i=res_i*ii_coeffn**iJump*ii_coeffn_jump
        ii_coeffn=ii_coeffn*ii_coeff2_jump
    END DO

    ! the final cleanup
    IF(integrate) THEN
       CALL poly_affine_t3t(poly_ijk,scaled_h,-resPosReal+p_shift,poly,&
           npoly=npoly)
    END IF

    ! rely on compiler to deallocate ALLOCATABLEs

    CONTAINS


! *****************************************************************************
!> \brief ...
! *****************************************************************************
SUBROUTINE j_loop()
    ! calculate j bounds
    ccj1 = ccj1_i1 * i +ccj1_i0
    ccj0 = (ccj0_i2*i+ccj0_i1)*i+ccj0_i0
    delta_j=ccj1*ccj1-4*ccj2*ccj0
    IF (delta_j<0) THEN
        RETURN
    END IF
    sqDj=SQRT(delta_j)
    jmin=CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
    jmax=FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))

   IF (integrate) THEN
       poly_jk=0.0_dp
   ELSE
       CALL poly_p_eval3b(poly_ijk(1),size_ijk,REAL(i,dp),&
            poly_jk(1),size_jk,&
               npoly=npoly,grad=grad,xi=xi(1))
   ENDIF

    IF (period(1)==0) THEN
        jmin=MAX(l_bounds(1,1)-shiftPos(1),jmin)
        jmax=MIN(l_bounds(2,1)-shiftPos(1),jmax)
    END IF

    ! pos j side
    j_coeffn_j=j_coeffn_i
    k_coeffn_j=k_coeffn_i
    jj_coeffn=j_coeffn_j*jj_coeff0
    res_j=res_i

    jJump=ndim(1)-l_bounds(2,1)+l_bounds(1,1)
    jstart=MAX(0,jmin)
    ijShift=shiftPos(1)+jstart-l_bounds(2,1)
    IF (ijShift>0) ijShift=ijShift+ndim(1)-1
    ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
    !ijShift=CEILING(REAL(shiftPos(1)+jstart-l_bounds(2,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1)
    jstart=MAX(ijShift+l_bounds(1,1),jstart)
    jend=MIN(ijShift+l_bounds(2,1),jmax)
    ij=(jstart-ijShift-l_bounds(1,1)+l_shift(1))

    IF (jstart>0) THEN
        k_coeffn_j=k_coeffn_i*jk_coeff0**jstart
        jj_coeffn=j_coeffn_j*jj_coeff0**(2*jstart+1)
        res_j=res_i*(jj_coeff0**jstart*j_coeffn_j)**jstart
    END IF
    DO
        DO j=jstart,jend
            kmin=k_bounds(0,j-jmin,i-imin)
            kmax=k_bounds(1,j-jmin,i-imin)
            ! do k loop
            IF (res_j/=0.0_dp.and.k_coeffn_j/=0.0_dp.and.kmin<=kmax &
                .and.ABS(res_j)>small) THEN
                CALL k_loop
            END IF
            k_coeffn_j=k_coeffn_j*jk_coeff0
            res_j=res_j*jj_coeffn
            jj_coeffn=jj_coeffn*jj_coeff2
            ij=ij+1
        END DO
        jstart=jend+jJump
        IF (jstart>jmax) EXIT
        ij=l_shift(1)
        jend=MIN(jend+ndim(1),jmax)
        IF (jJump/=1) THEN ! remove if?
            k_coeffn_j=k_coeffn_i*jk_coeff0**jstart
            jj_coeffn=j_coeffn_j*jj_coeff0**(2*jstart+1)
            res_j=res_i*(jj_coeff0**jstart*j_coeffn_j)**jstart
        END IF
    END DO

    ! neg j side
    j_coeffn_j=1.0_dp/j_coeffn_i
    k_coeffn_j=k_coeffn_i
    jj_coeffn=j_coeffn_j*jj_coeff0
    res_j=res_i

    jstart=MIN(-1,jmax)
    ijShift=shiftPos(1)+jstart-l_bounds(1,1)
    IF (ijShift<0) ijShift=ijShift-ndim(1)+1
    ijShift=(ijShift/ndim(1))*ndim(1)-shiftPos(1)
    !ijShift=FLOOR(REAL(shiftPos(1)+jstart-l_bounds(1,1))/REAL(ndim(1)))*ndim(1)-shiftPos(1))
    jstart=MIN(ijShift+l_bounds(2,1),jstart)
    jend=MAX(ijShift+l_bounds(1,1),jmin)
    ij=(jstart-ijShift-l_bounds(1,1)+l_shift(1))
    IF (jstart<-1) THEN
        k_coeffn_j=k_coeffn_i*jk_coeff0**(jstart+1)
        jj_coeffn=j_coeffn_j*jj_coeff0**(-(2*jstart+1))
        res_j=res_i*(jj_coeff0**(-jstart-1)*j_coeffn_j)**(-jstart-1)
    END IF
    DO
        DO j=jstart,jend,-1
            k_coeffn_j=k_coeffn_j*jk_coeff1
            res_j=res_j*jj_coeffn
            jj_coeffn=jj_coeffn*jj_coeff2

            kmin=k_bounds(0,j-jmin,i-imin)
            kmax=k_bounds(1,j-jmin,i-imin)
            ! perform k loop
            IF (res_j/=0.0_dp.and.k_coeffn_j/=0.0_dp.and.kmin<=kmax &
                .and.ABS(res_j)>small) THEN
                CALL k_loop
            END IF
            ij=ij-1
        END DO
        jstart=jend-jJump
        IF (jstart<jmin) EXIT
        jend=MAX(jend-ndim(1),jmin)
        ij=l_ub(1)
        IF (jJump/=1) THEN ! remove if?
            k_coeffn_j=k_coeffn_i*jk_coeff0**(jstart+1)
            jj_coeffn=j_coeffn_j*jj_coeff0**(-(2*jstart+1))
            res_j=res_i*(jj_coeff0**(-jstart-1)*j_coeffn_j)**(-jstart-1)
        END IF
    END DO

    IF(integrate) THEN
       CALL poly_padd_uneval3b(poly_ijk(1),size_ijk,REAL(i,dp),&
            poly_jk(1),size_jk,&
           npoly=npoly,grad=grad,xi=xi(1))
       !CALL poly_padd_uneval3(poly_ijk,REAL(i,dp),poly_jk,npoly=npoly)
    ENDIF
END SUBROUTINE

! *****************************************************************************
!> \brief ...
! *****************************************************************************
SUBROUTINE k_loop()
    IF(integrate) THEN
       SELECT CASE(grad)
          CASE(1)
              CALL kloop1_int
          CASE(2)
              CALL kloop2_int
          CASE(3)
              CALL kloop3_int
          CASE(4)
              CALL kloop4_int
          CASE(5)
              CALL kloop5_int
          CASE(6)
              CALL kloop6_int
          CASE(7)
              CALL kloop7_int
          CASE(8)
              CALL kloop8_int
          CASE default
              CALL kloopdefault_int
       END SELECT
    ELSE
       SELECT CASE(grad)
          CASE(1)
              CALL kloop1_col
          CASE(2)
              CALL kloop2_col
          CASE(3)
              CALL kloop3_col
          CASE(4)
              CALL kloop4_col
          CASE(5)
              CALL kloop5_col
          CASE(6)
              CALL kloop6_col
          CASE(7)
              CALL kloop7_col
          CASE(8)
              CALL kloop8_col
          CASE default
              CALL kloopdefault_col
       END SELECT
    ENDIF
END SUBROUTINE

! *****************************************************************************
!> \brief ...
! *****************************************************************************
SUBROUTINE kloop1_int()
        INTEGER, PARAMETER :: grad_val=1


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
# 94 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
# 144 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
# 190 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
# 238 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 256 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 921 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop2_int()
        INTEGER, PARAMETER :: grad_val=2


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
# 98 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
# 152 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
# 202 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
# 254 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 273 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 930 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop3_int()
        INTEGER, PARAMETER :: grad_val=3


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
# 102 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
# 160 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
# 214 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
# 270 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 290 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 939 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop4_int()
        INTEGER, PARAMETER :: grad_val=4


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
# 106 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
# 168 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
# 226 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
# 286 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 307 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 948 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop5_int()
        INTEGER, PARAMETER :: grad_val=5


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
# 110 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
# 176 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
# 238 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
# 302 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 324 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 957 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop6_int()
        INTEGER, PARAMETER :: grad_val=6


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
# 114 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
# 184 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
# 250 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
# 318 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 341 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 966 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop7_int()
        INTEGER, PARAMETER :: grad_val=7


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
# 118 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
# 192 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
# 262 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
# 334 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 358 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 975 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop8_int()
        INTEGER, PARAMETER :: grad_val=8


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(8)=k_vals(8)+p_kk
# 122 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(8)=k_vals(8)+p_kk
# 200 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(8)=k_vals(8)+p_kk
# 274 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(1)=k_vals(1)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(2)=k_vals(2)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(3)=k_vals(3)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(4)=k_vals(4)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(5)=k_vals(5)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(6)=k_vals(6)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(7)=k_vals(7)+p_kk
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(8)=k_vals(8)+p_kk
# 350 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 375 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 984 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloopdefault_int()
    INTEGER                                  :: grad_val

        grad_val=grad


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!




    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)


    k_vals=0.0_dp

    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad_val
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
# 98 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad_val
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
# 152 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad_val
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
# 202 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"

                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2

                gval=grid(ik,ij,ii)*res_k
                k_vals(0)=k_vals(0)+gval
                p_kk=gval
                DO kgrad=1,grad_val
                    p_kk=p_kk*REAL(k,dp)
                    k_vals(kgrad)=k_vals(kgrad)+p_kk
                END DO
# 254 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 273 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
    CALL poly_padd_uneval2b(poly_jk(1),size_jk,REAL(j,dp),k_vals(0),&
        size_k,npoly=npoly,grad=grad_val,xi=xi(1))
# 995 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2

    END SUBROUTINE

! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop1_col()
    INTEGER, PARAMETER                       :: grad_val = 1


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop1.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend







                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2







                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend







                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2







                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1005 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop2_col()
    INTEGER, PARAMETER                       :: grad_val = 2


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 80 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 134 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 184 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 236 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop2.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1013 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop3_col()
    INTEGER, PARAMETER                       :: grad_val = 3


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 82 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 140 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 194 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 250 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop3.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1021 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop4_col()
    INTEGER, PARAMETER                       :: grad_val = 4


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 84 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 146 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 204 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 264 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop4.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1029 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop5_col()
    INTEGER, PARAMETER                       :: grad_val = 5


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 86 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 152 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 214 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 278 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop5.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1037 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop6_col()
    INTEGER, PARAMETER                       :: grad_val = 6


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 88 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 158 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 224 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 292 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop6.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1045 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop7_col()
    INTEGER, PARAMETER                       :: grad_val = 7


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 90 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 164 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 234 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 306 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop7.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1053 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloop8_col()
    INTEGER, PARAMETER                       :: grad_val = 8


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 92 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(8)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 170 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(8)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 244 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                    p_v=p_v+poly_k(8)*p_kk
                    p_kk=p_kk*REAL(k,dp)



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 320 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop8.f90"
                p_v=poly_k(0)
                p_kk=k
                    p_v=p_v+poly_k(1)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(2)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(3)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(4)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(5)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(6)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(7)*p_kk
                    p_kk=p_kk*k
                    p_v=p_v+poly_k(8)*p_kk
                    p_kk=p_kk*k



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1061 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE kloopdefault_col()
    INTEGER                                  :: grad_val

        grad_val=grad

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90" 1
! to generate the unrolled loops use the following python scripts
!import os,re
!
!doStart=re.compile(r" *do +kgrad *= *(?P<start>[0-9]+) *, *grad_val *$",re.IGNORECASE)
!doEnd=re.compile(r" *end do",re.IGNORECASE)
!for i in xrange(1,9):
!    fIn=file('colloc_int_kloop.f90')
!    fout=file('colloc_int_kloop%d.f90'%(i),'w')
!    while 1:
!        line=fIn.readline()
!        if not line:
!            break
!        m=doStart.match(line)
!        if m:
!            body=""
!            while 1:
!                line=fIn.readline()
!                if not line:
!                    raise Exception('missing end do')
!                if doEnd.match(line):
!                    print 'unrolling',repr(body)
!                    for j in xrange(int(m.group('start')),i+1):
!                        fout.write(body.replace('kgrad',str(j)))
!                    break
!                body+=line
!        else:
!            fout.write(line)
!    fout.close()
!

    CALL poly_p_eval2b(poly_jk(1),size_jk,REAL(j,dp),&
         poly_k(0),size_k,npoly=npoly,grad=grad_val,xi=xi(1))

    ! starting point
    kJump=ndim(2)-l_bounds(2,2)+l_bounds(1,2)
    kstart=MAX(0,kmin)
    ikShift=shiftPos(2)-l_bounds(2,2)+kstart
    IF (ikShift>0) ikShift=ikShift+ndim(2)-1
    ikShift=(ikShift/ndim(2))*ndim(2)-shiftPos(2)
    ! ikShift=CEILING(REAL(shiftPos(2)-l_bounds(2,2)+kstart)/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart=MAX(ikShift+l_bounds(1,2),kstart)
    kend=MIN(ikShift+l_bounds(2,2),kmax)
    ikstart=kstart-ikShift-l_bounds(1,2)+l_shift(2)
    kstart2=MIN(-1,kmax)
    ikShift2=shiftPos(2)+kstart2-l_bounds(1,2)
    IF (ikShift2<0) ikShift2=ikShift2-ndim(2)+1
    ikShift2=(ikShift2/ndim(2))*ndim(2)-shiftPos(2)
    !ikShift2=FLOOR(REAL(shiftPos(2)+kstart2-l_bounds(1,2))/REAL(ndim(2)))*ndim(2)-shiftPos(2)
    kstart2=MIN(ikShift2+l_bounds(2,2),kstart2)
    kend2=MAX(ikShift2+l_bounds(1,2),kmin)
    ikstart2=kstart2-ikShift2-l_bounds(1,2)+l_shift(2)




    IF (kJump/=1 .AND. (ikstart+kmax-kstart>=ndim(2)+l_shift(2) .OR.&
        ikstart2+kmin-kstart2<=l_ub(2)-ndim(2))) THEN
        ! will wrap
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 80 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                DO kgrad=1,grad_val
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                END DO



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 134 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                p_v=poly_k(0)
                p_kk=k
                DO kgrad=1,grad_val
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*k
                END DO



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END DO
    ELSE
        ! no jump
        ! pos k side
        k_coeffn_k=k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart
        ik=ikstart
        IF (k>0) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(2*kstart+1)
            res_k=res_j*(kk_coeff0**kstart*k_coeffn_k)**kstart
        END IF
        DO
            DO k=kstart,kend
# 184 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                p_v=poly_k(0)
                p_kk=REAL(k,dp)
                DO kgrad=1,grad_val
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*REAL(k,dp)
                END DO



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF



                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
                ik=ik+1
            END DO
            kstart=kend+kJump
            IF (kstart>kmax) EXIT
            kend=MIN(kend+ndim(2),kmax)
            ik=l_shift(2)
        END DO

        ! neg k side
        k_coeffn_k=1.0_dp/k_coeffn_j
        kk_coeffn=k_coeffn_k*kk_coeff0
        res_k=res_j
        k=kstart2
        ik=ikstart2
        IF (k<-1) THEN
            kk_coeffn=k_coeffn_k*kk_coeff0**(-(2*kstart2+1))
            res_k=res_j*(kk_coeff0**(-kstart2-1)*k_coeffn_k)**(-kstart2-1)
        END IF
        DO
            DO k=kstart2,kend2,-1
                res_k=res_k*kk_coeffn
                kk_coeffn=kk_coeffn*kk_coeff2
# 236 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/colloc_int_kloop.f90"
                p_v=poly_k(0)
                p_kk=k
                DO kgrad=1,grad_val
                    p_v=p_v+poly_k(kgrad)*p_kk
                    p_kk=p_kk*k
                END DO



                IF ( PRESENT ( lgrid ) ) THEN
                  ig = lgrid%ldim * ithread + ii * (l_bounds(2,2)-l_bounds(1,2)+1) * (l_bounds(2,1)-l_bounds(1,1)+1) + &
                       ij * (l_bounds(2,2)-l_bounds(1,2)+1) + ik + 1
                  lgrid%r(ig)=lgrid%r(ig) + p_v*res_k
                ELSE
                  grid(ik,ij,ii) = grid(ik,ij,ii) + p_v*res_k
                END IF


                ik=ik-1
            END DO
            kstart2=kend2-kJump
            IF (kstart2<kmin) EXIT
            kend2=MAX(kend2-ndim(2),kmin)
            ik=l_ub(2)
        END DO
    END IF
# 1070 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/gauss_colloc.F" 2
    END SUBROUTINE


END SUBROUTINE colloc_int_body

END MODULE gauss_colloc
