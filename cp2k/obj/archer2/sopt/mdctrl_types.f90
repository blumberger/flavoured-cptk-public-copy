# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mdctrl_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mdctrl_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A common interface for passing a callback into the md_run loop.
!> \par History
!> \author Ole
! *****************************************************************************
MODULE mdctrl_types

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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mdctrl_types.F" 2

 IMPLICIT NONE
 PRIVATE

 TYPE glbopt_mdctrl_data_type
   INTEGER                                    :: md_bump_counter
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: epot_history
   INTEGER                                    :: output_unit
   INTEGER                                    :: itimes
   INTEGER                                    :: bump_steps_upwards
   INTEGER                                    :: bump_steps_downwards
   INTEGER                                    :: md_bumps_max
 END TYPE glbopt_mdctrl_data_type

 TYPE mdctrl_type
    TYPE(glbopt_mdctrl_data_type), POINTER                 :: glbopt => Null()
    !... and possible more in the future
 END TYPE mdctrl_type


 PUBLIC :: mdctrl_type, glbopt_mdctrl_data_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mdctrl_types'

END MODULE mdctrl_types

