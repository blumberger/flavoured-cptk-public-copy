# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/linear_systems.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/linear_systems.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides interfaces to LAPACK routines for factorisation and
!>      linear system solving
!> \note
!>      We are using LAPACK interfaces, so please make sure in IBM/AIX you have
!>      the lapack library before essl: "xlf90 ... -llapack -lessl" !!!
!> \par History
!>      none
!> \author JGH (30-5-2001)
! *****************************************************************************
MODULE linear_systems

  
  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_sgesv

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
# 22 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/linear_systems.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'linear_systems'

  PUBLIC :: solve_system

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param mysize ...
!> \param eigenvectors ...
! *****************************************************************************
SUBROUTINE solve_system ( matrix, mysize, eigenvectors )

    REAL(KIND=dp), INTENT(INOUT)             :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    REAL(KIND=dp), INTENT(INOUT)             :: eigenvectors( :, : )

    CHARACTER(len=*), PARAMETER :: routineN = 'solve_system', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, lda, ldb, nrhs, &
                                                ipiv( mysize )

  lda = SIZE ( matrix, 1 )
  ldb = SIZE ( eigenvectors, 1 )
  nrhs = SIZE ( eigenvectors, 2 )

  CALL lapack_sgesv ( mysize, nrhs, matrix, lda, ipiv, &
                      eigenvectors, ldb, info )
  IF ( info /= 0 ) THEN
     CALL cp__b("common/linear_systems.F",58,"Error in inversion")
  END IF

END SUBROUTINE solve_system

END MODULE linear_systems

