# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factors.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factors.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE structure_factors

  
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: twopi
  USE structure_factor_types,          ONLY: structure_factor_type
  USE termination,                     ONLY: stop_memory

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factors.F" 2

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'structure_factors'

  PRIVATE
  PUBLIC :: structure_factor_evaluate, structure_factor_allocate
  PUBLIC :: structure_factor_deallocate, structure_factor_init

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param exp_igr ...
! *****************************************************************************
  SUBROUTINE structure_factor_init(exp_igr)

    TYPE(structure_factor_type), &
      INTENT(INOUT)                          :: exp_igr

    NULLIFY(exp_igr % ex, exp_igr % ey, exp_igr % ez)
    NULLIFY(exp_igr % shell_ex, exp_igr % shell_ey, exp_igr % shell_ez)
    NULLIFY(exp_igr % core_ex, exp_igr % core_ey, exp_igr % core_ez)
    NULLIFY(exp_igr % centre, exp_igr % shell_centre, exp_igr % core_centre)

  END SUBROUTINE structure_factor_init

! *****************************************************************************
!> \brief ...
!> \param exp_igr ...
! *****************************************************************************
  SUBROUTINE structure_factor_deallocate ( exp_igr )

    TYPE(structure_factor_type), &
      INTENT(INOUT)                          :: exp_igr

    CHARACTER(len=*), PARAMETER :: routineN = 'structure_factor_deallocate', &
      routineP = moduleN//':'//routineN

    DEALLOCATE ( exp_igr % ex)
    DEALLOCATE ( exp_igr % ey)
    DEALLOCATE ( exp_igr % ez)
    IF ( ASSOCIATED ( exp_igr % shell_ex ) ) THEN
       DEALLOCATE ( exp_igr % shell_ex)
       DEALLOCATE ( exp_igr % shell_ey)
       DEALLOCATE ( exp_igr % shell_ez)
    END IF
    IF ( ASSOCIATED ( exp_igr % core_ex ) ) THEN
       DEALLOCATE ( exp_igr % core_ex)
       DEALLOCATE ( exp_igr % core_ey)
       DEALLOCATE ( exp_igr % core_ez)
    END IF
    IF ( ASSOCIATED ( exp_igr % centre ) ) THEN
       DEALLOCATE ( exp_igr % centre)
    END IF
    IF ( ASSOCIATED ( exp_igr % shell_centre ) ) THEN
       DEALLOCATE ( exp_igr % shell_centre)
    END IF
    IF ( ASSOCIATED ( exp_igr % core_centre ) ) THEN
       DEALLOCATE ( exp_igr % core_centre)
    END IF

  END SUBROUTINE structure_factor_deallocate

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param bds ...
!> \param nparts ...
!> \param exp_igr ...
!> \param allocate_centre ...
!> \param allocate_shell_e ...
!> \param allocate_shell_centre ...
!> \param nshell ...
! *****************************************************************************
  SUBROUTINE structure_factor_allocate ( bds, nparts, exp_igr, &
                                      allocate_centre, allocate_shell_e, &
                                      allocate_shell_centre,nshell)

    INTEGER, DIMENSION(:, :), INTENT(IN)     :: bds
    INTEGER, INTENT(IN)                      :: nparts
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    LOGICAL, INTENT(IN), OPTIONAL            :: allocate_centre, &
                                                allocate_shell_e, &
                                                allocate_shell_centre
    INTEGER, INTENT(IN), OPTIONAL            :: nshell

    CHARACTER(len=*), PARAMETER :: routineN = 'structure_factor_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isos

    ALLOCATE ( exp_igr % ex ( bds ( 1, 1 ):bds ( 2, 1 ) + 1, nparts ), &
         STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,114,'ex',0)
    ALLOCATE ( exp_igr % ey ( bds ( 1, 2 ):bds ( 2, 2 ) + 1, nparts ), &
         STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,117,'ey',0)
    ALLOCATE ( exp_igr % ez ( bds ( 1, 3 ):bds ( 2, 3 ) + 1, nparts ), &
        STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,120,'ez',0)
    NULLIFY ( exp_igr % centre )

    exp_igr % lb(1) = LBOUND ( exp_igr % ex, 1 )
    exp_igr % lb(2) = LBOUND ( exp_igr % ey, 1 )
    exp_igr % lb(3) = LBOUND ( exp_igr % ez, 1 )

    IF ( PRESENT ( allocate_centre ) ) THEN
      IF ( allocate_centre ) THEN
        ALLOCATE ( exp_igr % centre ( 3, nparts ), STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,130,&
                                          'centre',dp_size*3*nparts)
      END IF
    END IF

    IF( PRESENT(allocate_shell_e)) THEN
      IF(allocate_shell_e) THEN
      ALLOCATE ( exp_igr % shell_ex ( bds ( 1, 1 ):bds ( 2, 1 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,139,'shell_ex',0)
      ALLOCATE ( exp_igr % shell_ey ( bds ( 1, 2 ):bds ( 2, 2 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,142,'shell_ey',0)
      ALLOCATE ( exp_igr % shell_ez ( bds ( 1, 3 ):bds ( 2, 3 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,145,'shell_ez',0)
      NULLIFY ( exp_igr % shell_centre )

      ALLOCATE ( exp_igr % core_ex ( bds ( 1, 1 ):bds ( 2, 1 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,150,'core_ex',0)
      ALLOCATE ( exp_igr % core_ey ( bds ( 1, 2 ):bds ( 2, 2 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,153,'core_ey',0)
      ALLOCATE ( exp_igr % core_ez ( bds ( 1, 3 ):bds ( 2, 3 ) + 1, nshell ), &
         STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory(routineN,moduleN,156,'core_ez',0)
      NULLIFY ( exp_igr % core_centre )

      IF ( PRESENT ( allocate_shell_centre ) ) THEN
        IF ( allocate_shell_centre ) THEN
          ALLOCATE ( exp_igr % shell_centre ( 3, nshell ), STAT = isos )
          IF ( isos /= 0 ) &
             CALL stop_memory(routineN,moduleN,163,'shell_centre',dp_size*3*nparts)
          ALLOCATE ( exp_igr % core_centre ( 3, nshell ), STAT = isos )
          IF ( isos /= 0 ) &
             CALL stop_memory(routineN,moduleN,166,'core_centre',dp_size*3*nparts)
        END IF
      END IF
      END IF
    ELSE
      NULLIFY(exp_igr % shell_ex , exp_igr % shell_ey , exp_igr % shell_ez)
      NULLIFY(exp_igr % core_ex , exp_igr % core_ey , exp_igr % core_ez)
      NULLIFY(exp_igr % shell_centre, exp_igr % core_centre)
    END IF

  END SUBROUTINE structure_factor_allocate

! *****************************************************************************
!> \brief ...
!> \param delta ...
!> \param lb ...
!> \param ex ...
!> \param ey ...
!> \param ez ...
! *****************************************************************************
  SUBROUTINE structure_factor_evaluate ( delta, lb, ex, ey, ez )

    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: delta
    INTEGER, DIMENSION(3), INTENT(IN)        :: lb
    COMPLEX(KIND=dp), DIMENSION(lb(1):), &
      INTENT(out)                            :: ex
    COMPLEX(KIND=dp), DIMENSION(lb(2):), &
      INTENT(out)                            :: ey
    COMPLEX(KIND=dp), DIMENSION(lb(3):), &
      INTENT(out)                            :: ez

    COMPLEX(KIND=dp)                         :: fm, fp
    INTEGER                                  :: j, l0, l1, m0, m1, n0, n1
    REAL(KIND=dp)                            :: vec( 3 )

    l0 = LBOUND ( ex, 1 )
    l1 = UBOUND ( ex, 1 )
    m0 = LBOUND ( ey, 1 )
    m1 = UBOUND ( ey, 1 )
    n0 = LBOUND ( ez, 1 )
    n1 = UBOUND ( ez, 1 )

    ! delta is in scaled coordinates
    vec ( : ) = twopi * ( delta ( : ) + 0.5_dp  )

    ex ( l0 ) = 1.0_dp
    ey ( m0 ) = 1.0_dp
    ez ( n0 ) = 1.0_dp
    ex ( l1 ) = 1.0_dp
    ey ( m1 ) = 1.0_dp
    ez ( n1 ) = 1.0_dp

    fp = CMPLX ( COS ( vec ( 1 ) ), -SIN ( vec ( 1 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -l0
       ex (  j + l0 ) = ex (  j + l0 - 1 ) * fp
       ex ( -j + l1 ) = ex ( -j + l1 + 1 ) * fm
    END DO

    fp = CMPLX ( COS ( vec ( 2 ) ), -SIN ( vec ( 2 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -m0
       ey (  j + m0 ) = ey (  j + m0 - 1 ) * fp
       ey ( -j + m1 ) = ey ( -j + m1 + 1 ) * fm
    END DO

    fp = CMPLX ( COS ( vec ( 3 ) ), -SIN ( vec ( 3 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -n0
       ez (  j + n0 ) = ez (  j + n0 - 1 ) * fp
       ez ( -j + n1 ) = ez ( -j + n1 + 1 ) * fm
    END DO

  END SUBROUTINE structure_factor_evaluate

END MODULE structure_factors
