# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_move_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_move_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief tree nodes creation, searching, deallocation, references etc.
!> \par History
!>      11.2012 created [Mandes Schönherr] 
!> \author Mandes 11/2012
! *****************************************************************************

MODULE tmc_move_types
  USE kinds,                           ONLY: default_string_length,&
                                             dp

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/tmc/tmc_move_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_move_types'

  !-- list of available move types
  INTEGER, PARAMETER, PUBLIC :: mv_type_none           = 0
  INTEGER, PARAMETER, PUBLIC :: mv_type_swap_conf      = 1 ! swapping of 2 configurations of different temperature
  INTEGER, PARAMETER, PUBLIC :: mv_type_atom_trans     = 2 ! atom translation (done in every posible direction)
  INTEGER, PARAMETER, PUBLIC :: mv_type_mol_trans      = 3 ! molecule translation (done in every posible direction)
  INTEGER, PARAMETER, PUBLIC :: mv_type_mol_rot        = 4 ! molecule rotation
  INTEGER, PARAMETER, PUBLIC :: mv_type_proton_reorder = 5 ! reordering the protons within a chain of molecules
  INTEGER, PARAMETER, PUBLIC :: mv_type_atom_swap      = 6 ! swaps two atoms of different type
  INTEGER, PARAMETER, PUBLIC :: mv_type_MD             = 7 ! certain amount of MD steps
  INTEGER, PARAMETER, PUBLIC :: mv_type_volume_move    = 8 ! volume move for NPT simulations
  INTEGER, PARAMETER, PUBLIC :: mv_type_gausian_adapt  = 9 ! gaussian adaptation
  INTEGER, PARAMETER, PUBLIC :: mv_type_NMC_moves      = 10! indentifies the Nested Monte Carlo move for master
  INTEGER, PARAMETER, PUBLIC :: nr_mv_types            = 10!-- allways update the number of possible types!!

  PUBLIC :: tmc_move_type, move_types_create, move_types_release

  TYPE tmc_move_type
     !-- mv_type, handling indeces to move type (are equal for all several configurations/temperatures)
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mv_weight
     !-- mv_size, moves are normaly done in interval ]-mv_size, mv_size[
     ! 1st dimension are the different types, 2nd dim for configuration/temperature
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: mv_size
     !-- acc_prob, probability of acceptance of a certain move type for a certain temperature
     ! 1st dimension are the different move types, 2nd dim for configuration/temperature
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: acc_prob
     !-- count, remembers the certain amount of moves of certain a move type and temperature
     ! 1st dimension are the different types, 2nd dim for config./Temp
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: mv_count
     !-- count, remembers the certain amount of accepted moves of a certain move type and temperature
     ! 1st dimension are the different types, 2nd dim for config./Temp
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: acc_count
     !-- subbox_prob, probability of acceptance of a certain move type within subbox,
     !   done in Nested Monte Carlo routine
     !   the moves are rejected if atom or center of mass leaves the subbox
     !   1st dimension are the different move types
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: subbox_acc_count
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: subbox_count
     TYPE(list_atoms), DIMENSION(:), POINTER :: atom_lists

     !-- nmc_acc_prob, probability of acceptance of a certain move type,
     !   done in Nested Monte Carlo routine, for different potential
     !   1st dimension are the different move types
!     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: nmc_nr_acc
!     INTEGER, DIMENSION(:), ALLOCATABLE :: nmc_count
  END TYPE tmc_move_type

  TYPE list_atoms
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atoms
  END TYPE list_atoms
CONTAINS

! *****************************************************************************
!> \brief allocating the module variables
!> \param move_types pointer to the structure which should be deallocated
!> \param nr_temp ...
!> \author Mandes 11.2012
!> \note deallocating the module variables
! *****************************************************************************
  SUBROUTINE move_types_create(move_types, nr_temp)
    TYPE(tmc_move_type), POINTER             :: move_types
    INTEGER                                  :: nr_temp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'move_types_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(move_types)))CALL cp__a("tmc/tmc_move_types.F",90)

    ALLOCATE(move_types)
    ALLOCATE(move_types%mv_weight(nr_mv_types))
    move_types%mv_weight(:) = 0.0_dp
    ALLOCATE(move_types%mv_size(nr_mv_types,nr_temp))
    move_types%mv_size(:,:) = 0.0_dp
    ALLOCATE(move_types%acc_prob(0:nr_mv_types,nr_temp))
    move_types%acc_prob(:,:) = 0.0_dp
    ALLOCATE(move_types%mv_count(0:nr_mv_types,nr_temp))
    move_types%mv_count(:,:) = 0
    ALLOCATE(move_types%acc_count(0:nr_mv_types,nr_temp))
    move_types%acc_count(:,:) = 0
    ALLOCATE(move_types%subbox_acc_count(nr_mv_types,nr_temp))
    move_types%subbox_acc_count(:,:) = 0
    ALLOCATE(move_types%subbox_count(nr_mv_types,nr_temp))
    move_types%subbox_count(:,:) = 0
    NULLIFY(move_types%atom_lists)
  END SUBROUTINE move_types_create

! *****************************************************************************
!> \brief deallocating the module variables
!> \param move_types pointer to the structure which should be deallocated
!> \author Mandes 11.2012
!> \note deallocating the module variables
! *****************************************************************************
  SUBROUTINE move_types_release(move_types)
    TYPE(tmc_move_type), POINTER             :: move_types

    CHARACTER(LEN=*), PARAMETER :: routineN = 'move_types_release', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(move_types)))CALL cp__a("tmc/tmc_move_types.F",122)

    IF(ASSOCIATED(move_types%atom_lists)) DEALLOCATE(move_types%atom_lists)
    DEALLOCATE(move_types%mv_weight)
    DEALLOCATE(move_types%mv_size)
    DEALLOCATE(move_types%acc_prob)
    DEALLOCATE(move_types%mv_count)
    DEALLOCATE(move_types%acc_count)
    DEALLOCATE(move_types%subbox_acc_count)
    DEALLOCATE(move_types%subbox_count)
    DEALLOCATE(move_types)    
  END SUBROUTINE move_types_release

END MODULE tmc_move_types
