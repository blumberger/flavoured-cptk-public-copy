# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_stati.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_stati.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief tree nodes creation, searching, deallocation, references etc.
!> \par History
!>      11.2012 created [Mandes Schönherr] 
!> \author Mandes
! *****************************************************************************

MODULE tmc_stati

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/../base/base_uses.f90" 1
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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_stati.F" 2
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_stati'

  ! IO
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_default_trajectory_file_name  = "tmc_trajectory.dat"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_default_restart_out_file_name = "tmc_restart.dat"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_default_restart_in_file_name  = "tmc_restart.last"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_energy_worker_out_file_name   = "tmc_E_worker.out"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_NMC_worker_out_file_name      = "tmc_NMC_worker.out"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_master_out_file_name          = "tmc_master.out"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_ana_out_file_name             = "tmc_ana.out"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_default_dot_file_name         = "tmc_tree.dot"
  CHARACTER(LEN=*), PARAMETER, &
    PUBLIC                   :: tmc_default_unspecified_name      = "xxx_unspecified_xxx"

  ! TASK TYPES
  INTEGER, PARAMETER, PUBLIC :: task_type_MC                  = 1
  INTEGER, PARAMETER, PUBLIC :: task_type_ideal_gas           = 2
  INTEGER, PARAMETER, PUBLIC :: task_type_pauling             = 3
  INTEGER, PARAMETER, PUBLIC :: task_type_gaussian_adaptation = 4

  !-- communication status --
  !message tags
  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_OK                    = 0
  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_WAIT_FOR_NEW_TASK     = -42

  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_WORKER_INIT           = 666

  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_CALCULATING           = 1000
  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_FAILED                = 998
  INTEGER, PARAMETER, PUBLIC :: TMC_STATUS_STOP_RECEIPT          = 999

  INTEGER, PARAMETER, PUBLIC :: TMC_MESSAGE_INT                  = 1001
  INTEGER, PARAMETER, PUBLIC :: TMC_MASSAGE_REAL                 = 1002

  INTEGER, PARAMETER, PUBLIC :: TMC_CANCELING_MESSAGE            = 1003
  INTEGER, PARAMETER, PUBLIC :: TMC_CANCELING_RECEIPT            = 1004


!  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_REQUEST_REJECTED        = 1005

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_APPROX_ENERGY_REQUEST   = 1007
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_APPROX_ENERGY_RESULT    = 1008

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_START_CONF_REQUEST      = 1009
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_START_CONF_RESULT       = 1010

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_ENERGY_REQUEST          = 1011
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_ENERGY_RESULT           = 1012

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_NMC_REQUEST             = 1020
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_NMC_RESULT              = 1021
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_NMC_BROADCAST           = 1022

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_MD_REQUEST              = 1030
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_MD_RESULT               = 1031
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_MD_BROADCAST            = 1032

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_SCF_STEP_ENER_RECEIVE   = 2011

  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_INIT_ANALYSIS           = 3000
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_ANALYSIS_REQUEST        = 3001
  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_ANALYSIS_RESULT         = 3002

!  INTEGER, PARAMETER, PUBLIC :: TMC_STAT_SYNC_RND_SEED           = 1040

END MODULE tmc_stati

