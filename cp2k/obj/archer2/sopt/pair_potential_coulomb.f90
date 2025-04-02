# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pair_potential_coulomb.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pair_potential_coulomb.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pair_potential_coulomb

  USE erf_fn,                          ONLY: erf,&
                                             erfc
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: oorootpi
  USE pw_poisson_types,                ONLY: do_ewald_none

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
# 14 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pair_potential_coulomb.F" 2

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_coulomb'

  PUBLIC :: potential_coulomb

CONTAINS

! *****************************************************************************
!> \brief Evaluates the electrostatic energy and force
!> \param r2 ...
!> \param fscalar ...
!> \param qfac ...
!> \param ewald_type ...
!> \param alpha ...
!> \param beta ...
!> \param interaction_cutoff ...
!> \retval potential_coulomb ...
!> \author Toon.Verstraelen@gmail.com
! *****************************************************************************
  FUNCTION potential_coulomb(r2, fscalar, qfac, ewald_type, alpha, beta, &
       interaction_cutoff)

    REAL(KIND=dp), INTENT(IN)                :: r2
    REAL(KIND=dp), INTENT(INOUT)             :: fscalar
    REAL(KIND=dp), INTENT(IN)                :: qfac
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta, &
                                                interaction_cutoff
    REAL(KIND=dp)                            :: potential_coulomb

    REAL(KIND=dp), PARAMETER :: two_over_sqrt_pi = 2.0_dp*oorootpi

    REAL(KIND=dp)                            :: r, x1, x2

    r = SQRT(r2)
    IF (beta > 0.0_dp) THEN
       IF (ewald_type == do_ewald_none) THEN
          x2 = r*beta
          potential_coulomb = erf(x2)/r
          fscalar = fscalar + qfac*(potential_coulomb - &
            two_over_sqrt_pi*EXP(-x2*x2)*beta)/r2
       ELSE
          x1 = alpha*r
          x2 = r*beta
          potential_coulomb = (erf(x2) - erf(x1))/r
          fscalar = fscalar + qfac*(potential_coulomb + &
            two_over_sqrt_pi*(EXP(-x1*x1)*alpha-EXP(-x2*x2)*beta))/r2
       END IF
    ELSE
       IF (ewald_type == do_ewald_none) THEN
          potential_coulomb = 1.0_dp/r
          fscalar = fscalar + qfac*potential_coulomb/r2
       ELSE
          x1 = alpha*r
          potential_coulomb = erfc(x1)/r
          fscalar = fscalar + qfac*(potential_coulomb + &
            two_over_sqrt_pi*EXP(-x1*x1)*alpha)/r2
       END IF
    END IF

    potential_coulomb = qfac*(potential_coulomb - interaction_cutoff)

  END FUNCTION potential_coulomb

END MODULE pair_potential_coulomb
