# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factor_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factor_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE structure_factor_types

  
  USE kinds,                           ONLY: dp

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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/structure_factor_types.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: structure_factor_type

! *****************************************************************************
  TYPE structure_factor_type
     COMPLEX (KIND=dp), DIMENSION ( :, : ), POINTER :: ex, ey, ez
     COMPLEX (KIND=dp), DIMENSION ( :, : ), POINTER :: shell_ex, shell_ey, shell_ez
     COMPLEX (KIND=dp), DIMENSION ( :, : ), POINTER :: core_ex, core_ey, core_ez
     INTEGER, DIMENSION ( :, : ), POINTER :: centre, core_centre, shell_centre
     INTEGER :: lb ( 3 )
  END TYPE structure_factor_type

END MODULE structure_factor_types

