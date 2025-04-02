# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/md_util.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/md_util.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Molecular Dynamics
!> \author Teodoro Laino [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
MODULE md_util
  
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE md_energies,                     ONLY: md_write_output
  USE md_environment_types,            ONLY: md_environment_type

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/md_util.F" 2

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_util'

  PUBLIC :: md_output

CONTAINS

! *****************************************************************************
!> \brief collects the part of the MD that, basically, does the output
!> \param md_env ...
!> \param md_section ...
!> \param root_section ...
!> \param forced_io ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE md_output(md_env,md_section,root_section,forced_io)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: md_section, root_section
    LOGICAL, INTENT(IN)                      :: forced_io

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_print
    TYPE(section_vals_type), POINTER         :: print_section

    CALL timeset(routineN,handle)
    do_print = .TRUE.
    IF (forced_io) THEN
       print_section => section_vals_get_subs_vals(md_section,"PRINT")
       CALL section_vals_val_get(print_section,"FORCE_LAST",l_val=do_print)
    END IF
    IF (do_print) THEN
       ! Dumps all files related to the MD run
       CALL md_write_output(md_env)
       CALL write_restart(md_env=md_env,root_section=root_section)
    END IF
    CALL timestop(handle)

  END SUBROUTINE md_output

END MODULE md_util
