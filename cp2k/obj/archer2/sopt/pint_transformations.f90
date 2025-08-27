# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/pint_transformations.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/pint_transformations.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pint_transformations
  USE input_constants,                 ONLY: transformation_stage
  USE kinds,                           ONLY: dp
  USE pint_normalmode,                 ONLY: normalmode_f2uf,&
                                             normalmode_u2x,&
                                             normalmode_x2u
  USE pint_staging,                    ONLY: staging_f2uf,&
                                             staging_u2x,&
                                             staging_x2u
  USE pint_types,                      ONLY: pint_env_type

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/../base/base_uses.f90" 1
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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/pint_transformations.F" 2

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_transformations'

  PUBLIC :: pint_x2u, &
            pint_u2x, &
            pint_f2uf

CONTAINS

! ***************************************************************************
!> \brief Transforms from the x into the u variables
!>      (at the moment a staging transformation for the positions)
!> \param pint_env the path integral environment
!> \param ux will contain the u variable (defaults to pint_env%ux)
!> \param x the positions to transform (defaults to pint_env%x)
!> \par History
!>      Added normal mode transformation [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_x2u(pint_env,ux,x)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_x2u', &
      routineP = moduleN//':'//routineN

    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

    IF(.NOT.(ASSOCIATED(pint_env)))CALL cp__a("motion/pint_transformations.F",52)
    IF(.NOT.(pint_env%ref_count>0))CALL cp__a("motion/pint_transformations.F",53)
    my_x => pint_env%x
    my_ux => pint_env%ux
    IF (PRESENT(x)) my_x => x
    IF (PRESENT(ux)) my_ux => ux
    IF(.NOT.(ASSOCIATED(my_ux)))CALL cp__a("motion/pint_transformations.F",58)
    IF(.NOT.(ASSOCIATED(my_x)))CALL cp__a("motion/pint_transformations.F",59)

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_x2u(pint_env%staging_env,ux=my_ux,x=my_x)
    ELSE
      CALL normalmode_x2u(pint_env%normalmode_env,ux=my_ux,x=my_x)
    END IF
    RETURN
  END SUBROUTINE pint_x2u

! ***************************************************************************
!> \brief transform from the u variable to the x (inverse of x2u)
!> \param pint_env path integral environment
!> \param ux the u variable (positions to be backtransformed)
!> \param x will contain the positions
!> \par History
!>      Added normal mode transformation by hforbert
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_u2x(pint_env,ux,x)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: x

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_u2x', &
      routineP = moduleN//':'//routineN

    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

    IF(.NOT.(ASSOCIATED(pint_env)))CALL cp__a("motion/pint_transformations.F",90)
    IF(.NOT.(pint_env%ref_count>0))CALL cp__a("motion/pint_transformations.F",91)
    my_x => pint_env%x
    my_ux => pint_env%ux
    IF (PRESENT(x)) my_x => x
    IF (PRESENT(ux)) my_ux => ux
    IF(.NOT.(ASSOCIATED(my_ux)))CALL cp__a("motion/pint_transformations.F",96)
    IF(.NOT.(ASSOCIATED(my_x)))CALL cp__a("motion/pint_transformations.F",97)

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_u2x(pint_env%staging_env,ux=my_ux,x=my_x)
    ELSE
      CALL normalmode_u2x(pint_env%normalmode_env,ux=my_ux,x=my_x)
    END IF
    RETURN
  END SUBROUTINE pint_u2x

! ***************************************************************************
!> \brief transformation x to u for the forces
!> \param pint_env the path integral environment
!> \param uf will contain the accelerations for the transformed variables
!>        afterwards
!> \param f the forces to transform
!> \par History
!>      Added normal mode transformation [hforbert]
!>      Divide forces by the number of beads, since the replication
!>        environment (should) give raw forces [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_f2uf(pint_env,uf,f)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: uf
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: f

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_f2uf', &
      routineP = moduleN//':'//routineN

    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_uf

    IF(.NOT.(ASSOCIATED(pint_env)))CALL cp__a("motion/pint_transformations.F",131)
    IF(.NOT.(pint_env%ref_count>0))CALL cp__a("motion/pint_transformations.F",132)
    my_f => pint_env%f
    my_uf => pint_env%uf
    IF (PRESENT(f)) my_f => f
    IF (PRESENT(uf)) my_uf => uf
    IF(.NOT.(ASSOCIATED(my_uf)))CALL cp__a("motion/pint_transformations.F",137)
    IF(.NOT.(ASSOCIATED(my_f)))CALL cp__a("motion/pint_transformations.F",138)

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_f2uf(pint_env%staging_env,uf=my_uf,f=my_f)
    ELSE
      CALL normalmode_f2uf(pint_env%normalmode_env,uf=my_uf,f=my_f)
    END IF

   my_uf=my_uf/pint_env%mass_fict * pint_env%propagator%physpotscale
    RETURN
  END SUBROUTINE pint_f2uf

END MODULE pint_transformations
