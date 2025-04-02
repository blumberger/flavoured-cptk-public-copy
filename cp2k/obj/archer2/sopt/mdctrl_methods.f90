# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/mdctrl_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/mdctrl_methods.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A common interface (wrapper) for a callback into the md_run loop.
!>        Currently this is only used by the glbopt machinery, but its meant
!>        to be extended if others need to controll the md_run loop, too.
!>
!> \par History
!>      11.2012 created [Ole]
!> \author Ole
! *****************************************************************************
MODULE mdctrl_methods
  USE glbopt_callback,                 ONLY: glbopt_md_callback
  USE md_environment_types,            ONLY: md_environment_type
  USE mdctrl_types,                    ONLY: mdctrl_type

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
# 20 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/mdctrl_methods.F" 2

  IMPLICIT NONE
  PRIVATE

 PUBLIC :: mdctrl_callback

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mdctrl_methods'

CONTAINS

! *****************************************************************************
!> \brief This is called by md_run for each step during during its main-loop.
!> \param mdctrl data which is passed on to the wrapped client-routine
!> \param md_env contains the current state of the md_run
!> \param should_stop can be used to abort the md_run
! *****************************************************************************
  SUBROUTINE mdctrl_callback(mdctrl, md_env, should_stop)
    TYPE(mdctrl_type), POINTER               :: mdctrl
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(inout)                   :: should_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'mdctrl_callback', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(md_env)))CALL cp__a("motion/mdctrl_methods.F",44)
    IF(.NOT.(ASSOCIATED(mdctrl)))CALL cp__a("motion/mdctrl_methods.F",45)

    IF(ASSOCIATED(mdctrl%glbopt)) THEN
      CALL glbopt_md_callback(mdctrl%glbopt, md_env, should_stop)

    !ELSE IF(ASSOCIATED(mdctrl%your_own_hook)) THEN ...

    ELSE
      CALL cp__b("motion/mdctrl_methods.F",53,"mdctrl_callback: No hook found.")
    ENDIF

  END SUBROUTINE mdctrl_callback

END MODULE mdctrl_methods

