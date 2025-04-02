# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_create.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_create.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a QM/MM calculation with Force-Mixing
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type
  USE qmmm_create,                     ONLY: qmmm_env_create
  USE qmmm_types,                      ONLY: qmmm_env_get,&
                                             qmmm_env_release,&
                                             qmmm_env_type
  USE qmmmx_types,                     ONLY: qmmmx_env_type
  USE qmmmx_util,                      ONLY: setup_force_mixing_qmmm_sections,&
                                             update_force_mixing_labels

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
# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_create.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_create'

  PUBLIC :: qmmmx_env_create

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
  SUBROUTINE  qmmmx_env_create(qmmmx_env, root_section, para_env, globenv,&
       force_env_section, subsys_section, use_motion_section)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_create', &
      routineP = moduleN//':'//routineN

    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(qmmm_env_type), POINTER             :: dummy_qmmm_env
    TYPE(section_vals_type), POINTER         :: qmmm_core_section, &
                                                qmmm_extended_section, &
                                                qmmm_section

    NULLIFY(dummy_qmmm_env)

    qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM")

    CALL qmmm_env_create(dummy_qmmm_env, root_section, para_env, globenv,&
                         force_env_section, qmmm_section, subsys_section, use_motion_section, &
                         ignore_outside_box = .TRUE.)
    CALL qmmm_env_get(dummy_qmmm_env, subsys=subsys)

    CALL update_force_mixing_labels(subsys, qmmm_section)

    ! using CUR_INDICES and CUR_LABELS, create appropriate QM_KIND sections for two QM/MM calculations
    CALL setup_force_mixing_qmmm_sections(subsys, qmmm_section, qmmm_core_section, qmmm_extended_section)


    ALLOCATE(qmmmx_env)
    CALL qmmm_env_create(qmmmx_env%core, root_section, para_env, globenv,&
         force_env_section, qmmm_core_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE.)

    CALL qmmm_env_create(qmmmx_env%ext, root_section, para_env, globenv,&
         force_env_section, qmmm_extended_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE.)

    CALL section_vals_release(qmmm_core_section)
    CALL section_vals_release(qmmm_extended_section)
    CALL qmmm_env_release(dummy_qmmm_env)

  END SUBROUTINE qmmmx_env_create

END MODULE qmmmx_create
