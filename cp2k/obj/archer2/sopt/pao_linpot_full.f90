# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_linpot_full.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_linpot_full.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Full parametrization of Fock matrix, ie. the identity parametrization.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_linpot_full
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type

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
# 21 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_linpot_full.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_linpot_full'

  PUBLIC :: linpot_full_count_terms, linpot_full_calc_term

CONTAINS

! *****************************************************************************
!> \brief Initialize full linpot parametrization
!> \param qs_env ...
!> \param iatom ...
!> \param nterms ...
! *****************************************************************************
  SUBROUTINE linpot_full_count_terms(qs_env, iatom, nterms)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    INTEGER, INTENT(OUT)                     :: nterms

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_full_count_terms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, n
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

     CALL get_qs_env(qs_env,&
                      particle_set=particle_set,&
                      qs_kind_set=qs_kind_set)

     atomic_kind => particle_set(iatom)%atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind, kind_number=ikind)
     CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
     n = basis_set%nsgf

     nterms = n + n*(n-1)/2
  END SUBROUTINE linpot_full_count_terms


! *****************************************************************************
!> \brief Count number of linear terms in parametrization
!> \param kterm ...
!> \param block_V ...
! *****************************************************************************
  SUBROUTINE linpot_full_calc_term(kterm, block_V)
    INTEGER, INTENT(IN)                      :: kterm
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: block_V

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_full_calc_term', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, i, j, n

    block_V = 0.0
    n = SIZE(block_V,1)
    IF(kterm > n + n*(n-1)/2) CALL cp__b("pao_linpot_full.F",83,"kterm out of bounds")

    a = 0
    outer: &
    DO i=1, SIZE(block_V,1)
    DO j=1, SIZE(block_V,2)
      IF(i>=j) a = a + 1
      IF(a==kterm) EXIT outer
    ENDDO
    ENDDO outer

    block_V(i,j) = 1.0_dp
    block_V(j,i) = 1.0_dp
  END SUBROUTINE linpot_full_calc_term

END MODULE pao_linpot_full
