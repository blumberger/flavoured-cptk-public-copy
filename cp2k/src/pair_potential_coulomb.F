!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pair_potential_coulomb

  USE erf_fn,                          ONLY: erf,            &
                                             erfc
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: oorootpi
  USE pw_poisson_types,                ONLY: do_ewald_none,       &
                                             do_ewald_dsf,        &
                                             do_ewald_ewald,      &
                                             do_ewald_fssh_dsf,   &
                                             do_ewald_fssh_ewald
#include "./base/base_uses.f90"

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
!> \author Toon.Verstraelen@gmail.com (Editted by Matt 95ellismle@gmail.com)
! *****************************************************************************
  FUNCTION potential_coulomb(r2, fscalar, qfac, ewald_type, alpha, beta, &
                             interaction_cutoff, dsf_cutoff, Rc)

    REAL(KIND=dp), INTENT(IN)                :: r2
    REAL(KIND=dp), INTENT(INOUT)             :: fscalar
    REAL(KIND=dp), INTENT(IN)                :: qfac
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta, &
                                                interaction_cutoff
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: dsf_cutoff, Rc
    REAL(KIND=dp)                            :: potential_coulomb

    REAL(KIND=dp), PARAMETER :: two_over_sqrt_pi = 1.1283791670955126_dp

    REAL(KIND=dp)                            :: r, x1, x2, r_rc, dsf_const

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

       potential_coulomb = qfac*(potential_coulomb - interaction_cutoff)
    ELSE
       SELECT CASE(ewald_type)
       ! No Ewald
       CASE(do_ewald_none)
          potential_coulomb = 1.0_dp/r
          fscalar = fscalar + qfac*potential_coulomb/r2
          potential_coulomb = qfac*(potential_coulomb - interaction_cutoff)

       ! Damped Shifted Forces
       CASE (do_ewald_dsf, do_ewald_fssh_dsf)
         IF (.NOT. PRESENT(dsf_cutoff)) &
            CPABORT("Please specify a cutoff for your DSF simulation.")
         IF (r > Rc) THEN
            potential_coulomb = 0.0_dp
            RETURN
         END IF

         x1 = alpha * r
         r_rc = r - rc
         
         potential_coulomb = erfc(x1)/r
         !PRINT *, " "
         !PRINT *, " "
         !PRINT *, "E_dsf contrib = ", qfac*(potential_coulomb - interaction_cutoff &
         !                                    + (dsf_cutoff*r_rc))
         !PRINT *, "erfc(alpha r)/r = ", potential_coulomb
         !PRINT *, "erfc(alpha Rc)/Rc = ", interaction_cutoff
         !PRINT *, "term4 = ", dsf_cutoff
         !PRINT *, "full last term = ", dsf_cutoff * r_rc
         !PRINT *, "(r - Rc) = ", r_rc
         !PRINT *, "fterm1: ", potential_coulomb * qfac / r
         !PRINT *, "fterm2: ", two_over_sqrt_pi*EXP(-x1*x1)*alpha /r2
         !PRINT *, "fterm3: ", dsf_cutoff
         fscalar = fscalar + (qfac*( (potential_coulomb   &
                               + two_over_sqrt_pi*alpha*EXP(-x1*x1))/r2 &
                               - dsf_cutoff/r))

         potential_coulomb = qfac*(potential_coulomb - interaction_cutoff &
                                   + (dsf_cutoff*r_rc))
         
         !PRINT *, "alpha = ", alpha
         !PRINT *, "R = ", r
         !PRINT *, "Rc = ", rc
         !PRINT *, "R -Rc = ", r_rc
         !PRINT *, "dsf_cutoff = ", dsf_cutoff
         !PRINT *, "interaction_cutoff = ", interaction_cutoff
         !PRINT *, "q1q2 = ", qfac
         !PRINT *, " "
         !PRINT *, "E_dsf = ", potential_coulomb
         !PRINT *, " "
         !PRINT *, "fscalar: ", fscalar
         !PRINT *, " "
         !PRINT *, " "


       ! Ewald, PME, SMPE
       CASE DEFAULT
          x1 = alpha*r
          potential_coulomb = erfc(x1)/r
          fscalar = fscalar + qfac*(potential_coulomb + &
            two_over_sqrt_pi*EXP(-x1*x1)*alpha)/r2
          potential_coulomb = qfac*(potential_coulomb - interaction_cutoff)

       END SELECT
    END IF


  END FUNCTION potential_coulomb

END MODULE pair_potential_coulomb
