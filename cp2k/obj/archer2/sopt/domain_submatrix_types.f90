# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/domain_submatrix_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/domain_submatrix_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types to handle submatrices
!> \par History
!>       2013.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE domain_submatrix_types
  USE kinds,                           ONLY: dp

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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/domain_submatrix_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'domain_submatrix_types'

  INTEGER, PARAMETER, PUBLIC           :: select_row_col = 1
  INTEGER, PARAMETER, PUBLIC           :: select_row = 2

  PUBLIC :: domain_submatrix_type, domain_map_type

  ! submatrix storage with the meta-data necessary to convert 
  ! the submatrix into the DBCSR format
  TYPE domain_submatrix_type
    INTEGER                                       :: domain
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: mdata
    INTEGER                                       :: nbrows
    INTEGER                                       :: nbcols
    INTEGER                                       :: nrows
    INTEGER                                       :: ncols
    INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_row
    INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_col
    INTEGER, DIMENSION(:), ALLOCATABLE            :: size_brow
    INTEGER, DIMENSION(:), ALLOCATABLE            :: size_bcol
    INTEGER                                       :: nnodes
    INTEGER                                       :: groupid
  END TYPE domain_submatrix_type

  TYPE domain_map_type
    INTEGER, DIMENSION(:), ALLOCATABLE     :: index1
    INTEGER, DIMENSION(:,:), ALLOCATABLE   :: pairs
  END TYPE domain_map_type

END MODULE domain_submatrix_types


