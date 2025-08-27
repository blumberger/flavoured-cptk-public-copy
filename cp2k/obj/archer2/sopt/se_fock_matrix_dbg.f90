# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_dbg.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_dbg.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE se_fock_matrix_dbg
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_trace
  USE kinds,                           ONLY: dp
  USE qs_energy_types,                 ONLY: init_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE se_fock_matrix_coulomb,          ONLY: build_fock_matrix_coulomb_lr
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_dbg.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_dbg'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: dbg_energy_coulomb_lr

CONTAINS

! *****************************************************************************
!> \brief Debug routine for long-range energy (debug value of EWALD vs VALUE KS)
!> \param energy ...
!> \param ks_matrix ...
!> \param nspins ...
!> \param qs_env ...
!> \param matrix_p ...
!> \param calculate_forces ...
!> \param store_int_env ...
!> \author Teodoro Laino [tlaino] - 04.2009
! *****************************************************************************
  SUBROUTINE dbg_energy_coulomb_lr(energy, ks_matrix, nspins, qs_env, matrix_p,&
       calculate_forces, store_int_env)
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'dbg_energy_coulomb_lr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    REAL(KIND=dp)                            :: ecoul

! Zero structures only for debugging purpose

    CALL init_qs_energy(energy)
    DO ispin = 1, nspins
       CALL cp_dbcsr_set(ks_matrix(ispin)%matrix,0.0_dp)
    END DO

    ! Evaluate Coulomb Long-Range
    CALL build_fock_matrix_coulomb_lr(qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
         store_int_env)

    ! Compute the Hartree energy
    DO ispin=1,nspins
       CALL cp_dbcsr_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul)
       energy%hartree = energy%hartree + ecoul

       WRITE(*,*)ispin,"ECOUL ",ecoul
    END DO
    WRITE(*,*)"ENUC in DBG:",energy%core_overlap

    ! Debug statements
    WRITE(*,*)"TOTAL ENE",0.5_dp*energy%hartree+energy%core_overlap
    CALL cp__b("se_fock_matrix_dbg.F",79,"Debug energy for Coulomb Long-Range")

  END SUBROUTINE dbg_energy_coulomb_lr

END MODULE se_fock_matrix_dbg
