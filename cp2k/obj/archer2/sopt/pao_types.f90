# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types used by the PAO machinery
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE linesearch,                      ONLY: linesearch_type

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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_types'

  PUBLIC :: pao_env_type

  TYPE pao_env_type
     ! input values
     REAL(KIND=dp)                    :: eps_pao = 0.0_dp
     REAL(KIND=dp)                    :: mixing = 0.0_dp
     REAL(KIND=dp)                    :: penalty_dist = 0.0_dp
     REAL(KIND=dp)                    :: penalty_strength = 0.0_dp
     REAL(KIND=dp)                    :: check_unitary_tol = 0.0_dp
     REAL(KIND=dp)                    :: check_grad_tol = 0.0_dp
     REAL(KIND=dp)                    :: num_grad_eps = 0.0_dp
     INTEGER                          :: num_grad_order = -1
     INTEGER                          :: max_pao = -1
     INTEGER                          :: max_cycles = -1
     INTEGER                          :: parameterization = -1
     INTEGER                          :: cg_init_steps = -1
     CHARACTER(LEN=default_string_length) :: preopt_dm_file = ""

     ! output units
     INTEGER                          :: iw = -1
     INTEGER                          :: iw_cg = -1

     ! state variable
     INTEGER                          :: istep = -1
     REAL(KIND=dp)                    :: energy_prev = 0.0_dp
     REAL(KIND=dp)                    :: step_start_time = 0.0_dp
     TYPE(linesearch_type)            :: linesearch
     LOGICAL                          :: need_initial_scf = .FALSE.

     ! matrices
     TYPE(cp_dbcsr_type)              :: matrix_X
     TYPE(cp_dbcsr_type)              :: matrix_U
     TYPE(cp_dbcsr_type)              :: matrix_U0
     TYPE(cp_dbcsr_type)              :: matrix_H0
     TYPE(cp_dbcsr_type)              :: matrix_Y
     TYPE(cp_dbcsr_type)              :: matrix_N
     TYPE(cp_dbcsr_type)              :: matrix_N_inv
     TYPE(cp_dbcsr_type)              :: matrix_X_orig
     TYPE(cp_dbcsr_type)              :: matrix_G
     TYPE(cp_dbcsr_type)              :: matrix_G_prev
     TYPE(cp_dbcsr_type)              :: matrix_D
     TYPE(cp_dbcsr_type)              :: matrix_V_terms

  END TYPE

END MODULE pao_types
