# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/task_list_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/task_list_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types for task lists
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! *****************************************************************************
MODULE task_list_types
  
  USE kinds,                           ONLY: dp,&
                                             int_8

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/task_list_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

! *****************************************************************************
  TYPE task_list_type
    INTEGER(kind=int_8), DIMENSION(:, :), POINTER :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER       :: dist_ab
    INTEGER(kind=int_8), DIMENSION(:), POINTER    :: atom_pair_send, atom_pair_recv
    INTEGER                                       :: ntasks
    INTEGER, DIMENSION(:,:),POINTER               :: taskstart,taskstop
    INTEGER, DIMENSION(:),POINTER                 :: npairs
  END TYPE task_list_type

  PUBLIC :: task_list_type

  PUBLIC :: allocate_task_list,&
            deallocate_task_list

CONTAINS

! *****************************************************************************
!> \brief allocates and initialised the components of the task_list_type
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE allocate_task_list(task_list)
    TYPE(task_list_type), POINTER            :: task_list

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_task_list', &
      routineP = moduleN//':'//routineN

  ALLOCATE(task_list)

  NULLIFY(task_list%tasks)
  NULLIFY(task_list%dist_ab)
  NULLIFY(task_list%atom_pair_send)
  NULLIFY(task_list%atom_pair_recv)
  NULLIFY(task_list%taskstart)
  NULLIFY(task_list%taskstop)
  NULLIFY(task_list%npairs)
  task_list%ntasks=0
END SUBROUTINE allocate_task_list

! *****************************************************************************
!> \brief deallocates the components and the object itself
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE deallocate_task_list(task_list)
    TYPE(task_list_type), POINTER            :: task_list

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_task_list', &
      routineP = moduleN//':'//routineN

  IF (ASSOCIATED(task_list%tasks)) THEN
     DEALLOCATE(task_list%tasks)
  ENDIF
  IF (ASSOCIATED(task_list%dist_ab)) THEN
     DEALLOCATE(task_list%dist_ab)
  ENDIF
  IF (ASSOCIATED(task_list%atom_pair_send)) THEN
     DEALLOCATE(task_list%atom_pair_send)
  ENDIF
  IF (ASSOCIATED(task_list%atom_pair_recv)) THEN
     DEALLOCATE(task_list%atom_pair_recv)
  ENDIF
  IF (ASSOCIATED(task_list%taskstart)) THEN
     DEALLOCATE(task_list%taskstart)
  ENDIF
  IF (ASSOCIATED(task_list%taskstop)) THEN
     DEALLOCATE(task_list%taskstop)
  ENDIF
  IF (ASSOCIATED(task_list%npairs)) THEN
     DEALLOCATE(task_list%npairs)
  ENDIF

  DEALLOCATE(task_list)
END SUBROUTINE deallocate_task_list
END MODULE task_list_types
