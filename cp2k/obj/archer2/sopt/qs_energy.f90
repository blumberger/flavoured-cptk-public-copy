# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_energy.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_energy.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Perform a QUICKSTEP wavefunction optimization (single point)
!> \par History
!>      none
!> \author MK (29.10.2002)
! *****************************************************************************
MODULE qs_energy
  USE almo_scf,                        ONLY: almo_entry_scf
  USE cp_control_types,                ONLY: dft_control_type
  USE dm_ls_scf,                       ONLY: ls_scf
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_energy_utils,                 ONLY: qs_energies_compute_matrix_w,&
                                             qs_energies_init,&
                                             qs_energies_mp2,&
                                             qs_energies_properties
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_scf,                          ONLY: scf

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
# 27 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_energy.F" 2

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

  PUBLIC :: qs_energies

CONTAINS

! *****************************************************************************
!> \brief   Driver routine for QUICKSTEP single point wavefunction optimization.
!> \param qs_env ...
!> \param consistent_energies ...
!> \param calc_forces ...
!> \date    29.10.2002
!> \par History
!>          - consistent_energies option added (25.08.2005, TdK)
!>          - introduced driver for energy in order to properly decide between
!>            SCF or RTP (fschiff 02.09)
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE qs_energies (qs_env, consistent_energies, calc_forces)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: consistent_energies, &
                                                calc_forces

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: my_calc_forces, run_rtp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy

    CALL timeset(routineN,handle)

    my_calc_forces = .FALSE.
    IF(PRESENT(calc_forces)) my_calc_forces = calc_forces

    CALL qs_env_rebuild_pw_env(qs_env)

    CALL get_qs_env(qs_env=qs_env,run_rtp=run_rtp)
    IF(.NOT.run_rtp)THEN

      NULLIFY(dft_control, energy)
      CALL qs_energies_init(qs_env, my_calc_forces)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy)

      ! *** Perform a SCF run ***
      IF (dft_control%qs_control%do_ls_scf) THEN
        CALL ls_scf(qs_env=qs_env)
      ELSE IF (dft_control%qs_control%do_almo_scf) THEN
        CALL almo_entry_scf(qs_env=qs_env, calc_forces=my_calc_forces)
      ELSE
        CALL scf(qs_env=qs_env)

        ! Compute MP2 energy
        CALL qs_energies_mp2(qs_env, my_calc_forces)

        ! if calculate forces, time to compute the w matrix
        CALL qs_energies_compute_matrix_w(qs_env,my_calc_forces)

      END IF

      IF (PRESENT(consistent_energies)) THEN
        IF (consistent_energies) THEN
          CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
          ! add MP2 energy if necessary
          IF(ASSOCIATED(qs_env%mp2_env)) THEN
            energy%total = energy%total + energy%mp2
          END IF
        END IF
      END IF

      CALL qs_energies_properties(qs_env)

    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_energies

END MODULE qs_energy
