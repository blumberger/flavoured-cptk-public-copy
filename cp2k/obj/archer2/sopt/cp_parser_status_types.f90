# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_status_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_status_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow the storage of the parser status
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 08.2008
! *****************************************************************************
MODULE cp_parser_status_types
  USE cp_parser_buffer_types,          ONLY: buffer_type,&
                                             create_buffer_type,&
                                             release_buffer_type
  USE kinds,                           ONLY: max_line_length

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_status_types.F" 2

  IMPLICIT NONE
  PRIVATE

  TYPE status_type
     LOGICAL                                        :: in_use
     INTEGER                                        :: old_input_line_number
     INTEGER                                        :: old_icol
     INTEGER                                        :: old_icol1
     INTEGER                                        :: old_icol2
     CHARACTER(LEN=max_line_length)                 :: old_input_line
     ! Store status of the buffer
     TYPE(buffer_type), POINTER                     :: buffer
  END TYPE status_type

  PUBLIC :: status_type, create_status_type, release_status_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_status_types'

CONTAINS

! ****************************************************************************
!> \brief creates the parser status type
!> \param status ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_status_type(status)
    TYPE(status_type), POINTER               :: status

    CHARACTER(len=*), PARAMETER :: routineN = 'create_status_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(status)))CALL cp__a("input/cp_parser_status_types.F",49)
    ALLOCATE(status)
    status%in_use                = .FALSE.
    status%old_input_line        = ""
    status%old_input_line_number = HUGE(0)
    status%old_icol              = HUGE(0)
    status%old_icol1             = HUGE(0)
    status%old_icol2             = HUGE(0)
    NULLIFY(status%buffer)
    CALL create_buffer_type(status%buffer)
  END SUBROUTINE create_status_type

! ****************************************************************************
!> \brief releases the parser status type
!> \param status ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE release_status_type(status)
    TYPE(status_type), POINTER               :: status

    CHARACTER(len=*), PARAMETER :: routineN = 'release_status_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(status)))CALL cp__a("input/cp_parser_status_types.F",73)
    CALL release_buffer_type(status%buffer)
    DEALLOCATE(status)
  END SUBROUTINE release_status_type

END MODULE cp_parser_status_types
