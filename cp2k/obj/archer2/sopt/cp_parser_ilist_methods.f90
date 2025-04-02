# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_methods.F"
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
MODULE cp_parser_ilist_methods
  USE cp_log_handling,                 ONLY: cp_to_string
  USE cp_parser_ilist_types,           ONLY: ilist_type

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_ilist_methods.F" 2

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ilist_setup, ilist_update, ilist_reset
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_ilist_methods'

CONTAINS

! ****************************************************************************
!> \brief setup the integer listing type
!> \param ilist ...
!> \param token ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE ilist_setup(ilist, token)
    TYPE(ilist_type), POINTER                :: ilist
    CHARACTER(LEN=*)                         :: token

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind

    IF(.NOT.(ASSOCIATED(ilist)))CALL cp__a("input/cp_parser_ilist_methods.F",44)
    ind = INDEX(token,"..")
    READ (UNIT=token(:ind-1),FMT=*) ilist%istart
    READ (UNIT=token(ind+2:),FMT=*) ilist%iend
    IF(ilist%istart > ilist%iend)&
       CALL cp_abort(cp__l("input/cp_parser_ilist_methods.F",49),&
            "Invalid list range specified: "//&
            TRIM(ADJUSTL(cp_to_string(ilist%istart)))//".."//&
            TRIM(ADJUSTL(cp_to_string(ilist%iend))))
    ilist%nel_list = ilist%iend - ilist%istart + 1
    ilist%ipresent = ilist%istart
    ilist%in_use   = .TRUE.

  END SUBROUTINE ilist_setup

! ****************************************************************************
!> \brief updates the integer listing type
!> \param ilist ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE ilist_update(ilist)
    TYPE(ilist_type), POINTER                :: ilist

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_update', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(ilist)))CALL cp__a("input/cp_parser_ilist_methods.F",71)
    ilist%ipresent = ilist%ipresent + 1
    IF (ilist%ipresent>ilist%iend) THEN
       CALL ilist_reset(ilist)
    END IF
  END SUBROUTINE ilist_update

! ****************************************************************************
!> \brief updates the integer listing type
!> \param ilist ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE ilist_reset(ilist)
    TYPE(ilist_type), POINTER                :: ilist

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_reset', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(ilist)))CALL cp__a("input/cp_parser_ilist_methods.F",90)
    IF (ilist%ipresent==ilist%iend) THEN
       ilist%istart   = HUGE(0)
       ilist%iend     = HUGE(0)
       ilist%nel_list = HUGE(0)
       ilist%ipresent = HUGE(0)
       ilist%in_use   = .FALSE.
    END IF
  END SUBROUTINE ilist_reset

END MODULE cp_parser_ilist_methods
