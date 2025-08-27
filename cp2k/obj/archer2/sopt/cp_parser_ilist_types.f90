# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implementation (20.02.2008,AK)
!>      - integration into cp2k (22.02.2008,tlaino)
!>      - variables added (25.02.2008,AK)
!> \author Axel Kohlmeyer [AK] - CMM/UPenn Philadelphia
!> \date 25.02.2008
! *****************************************************************************
MODULE cp_parser_ilist_types
  


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/../base/base_uses.f90" 1
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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_types.F" 2
  IMPLICIT NONE
  PRIVATE

  TYPE ilist_type
     LOGICAL                              :: in_use
     INTEGER                              :: nel_list
     INTEGER                              :: istart, iend
     INTEGER                              :: ipresent
  END TYPE ilist_type

  PUBLIC :: ilist_type, create_ilist_type, release_ilist_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_ilist_types'

CONTAINS

! ****************************************************************************
!> \brief creates the integer listing type
!> \param ilist ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_ilist_type(ilist)
    TYPE(ilist_type), POINTER                :: ilist

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ilist_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(ilist)))CALL cp__a("input/cp_parser_ilist_types.F",46)
    ALLOCATE(ilist)
    ilist%istart   = HUGE(0)
    ilist%iend     = HUGE(0)
    ilist%nel_list = HUGE(0)
    ilist%ipresent = HUGE(0)
    ilist%in_use   = .FALSE.

  END SUBROUTINE create_ilist_type

! ****************************************************************************
!> \brief creates the integer listing type
!> \param ilist ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE release_ilist_type(ilist)
    TYPE(ilist_type), POINTER                :: ilist

    CHARACTER(len=*), PARAMETER :: routineN = 'release_ilist_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(ilist)))CALL cp__a("input/cp_parser_ilist_types.F",68)
    DEALLOCATE(ilist)
  END SUBROUTINE release_ilist_type

END MODULE cp_parser_ilist_types
