# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_integrate_potential.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_integrate_potential.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Integrate single or product functions over a potential on a RS grid
!> \par History
!>      Refactored from earlier versions by Joost VandeVondele (2002,2003,2007)
!> \author JGH [04.2014]
! *****************************************************************************
!
! This module acts as a common container for the routines from the low level
! modules
!           qs_integrate_potential_product
!           qs_integrate_potential_single
!           qs_integrate_potential_low
!
! *****************************************************************************
MODULE qs_integrate_potential
  USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
  USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
  USE qs_integrate_potential_single,   ONLY: integrate_ppl_rspace,&
                                             integrate_rho_nlcc,&
                                             integrate_scp_rspace,&
                                             integrate_v_core_rspace,&
                                             integrate_v_rspace_one_center

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
# 29 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_integrate_potential.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential'

! *** Public subroutines ***

  ! included from qs_integrate_potential_product
  PUBLIC :: integrate_v_rspace

  ! included from qs_integrate_potential_single
  PUBLIC :: integrate_v_rspace_one_center,&
            integrate_v_core_rspace,&
            integrate_ppl_rspace,&
            integrate_scp_rspace,&
            integrate_rho_nlcc

  ! included from qs_integrate_potential_low
  PUBLIC :: integrate_pgf_product_rspace

END MODULE qs_integrate_potential
