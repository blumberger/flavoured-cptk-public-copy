# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/velocity_verlet_control.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/velocity_verlet_control.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides an interface to the velocity-verlet based integrator
!>      routines for all ensembles
!> \author CJM (11-SEPT-2002)
! *****************************************************************************
MODULE velocity_verlet_control

  
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       isokin_ensemble, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble, &
       reftraj_ensemble
  USE integrator,                      ONLY: &
       isokin, langevin, nph_uniaxial, nph_uniaxial_damped, npt_f, npt_i, &
       nve, nve_respa, nvt, nvt_adiabatic, reftraj
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE simpar_types,                    ONLY: simpar_type

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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/velocity_verlet_control.F" 2

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'velocity_verlet_control'
  PUBLIC :: velocity_verlet

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param md_env ...
!> \param globenv ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE velocity_verlet ( md_env, globenv)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'velocity_verlet', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar

    CALL timeset (routineN, handle )

    ! Get force environment
    CALL get_md_env ( md_env, force_env=force_env, simpar=simpar)

    ! RESPA implemented only for NVE
    IF(simpar%do_respa .AND. nve_ensemble.NE.simpar % ensemble) THEN
       CALL cp__b("motion/velocity_verlet_control.F",64,"RESPA integrator not implemented for this ensemble")
    END IF

    ! Choice of the ensemble
    SELECT CASE (simpar%ensemble)
    CASE DEFAULT
       CALL cp__b("motion/velocity_verlet_control.F",70,"Integrator not implemented")
    CASE (nve_ensemble)
       IF(simpar%do_respa)THEN
          CALL nve_respa(md_env)
       ELSE
          CALL nve (md_env, globenv)
       END IF
    CASE (nvt_ensemble)
       CALL nvt (md_env, globenv)
    CASE (nvt_adiabatic_ensemble)
       CALL nvt_adiabatic (md_env, globenv)
    CASE (isokin_ensemble)
       CALL isokin (md_env)
    CASE (npt_i_ensemble)
       CALL npt_i (md_env, globenv)
    CASE (npt_f_ensemble)
       CALL npt_f (md_env, globenv)
    CASE (nph_uniaxial_ensemble)
       CALL nph_uniaxial (md_env)
    CASE (nph_uniaxial_damped_ensemble)
       CALL nph_uniaxial_damped (md_env)
    CASE (reftraj_ensemble)
       CALL reftraj (md_env)
    CASE (langevin_ensemble)
       CALL langevin(md_env)
    CASE (npe_f_ensemble)
       CALL npt_f (md_env, globenv)
    CASE (npe_i_ensemble)
       CALL npt_i (md_env, globenv)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE velocity_verlet

END MODULE velocity_verlet_control
