# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_inpp_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_inpp_types.F"
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
MODULE cp_parser_inpp_types
  
  USE kinds,                           ONLY: default_path_length

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_inpp_types.F" 2

  IMPLICIT NONE
  PRIVATE

  TYPE inpp_type
     ! for '@INCLUDE "some_file.inc"'
     ! currently open include file stack pointer
     INTEGER                              :: io_stack_level
     ! include file stack data
     INTEGER, POINTER, DIMENSION(:)       :: io_stack_channel,&
                                             io_stack_lineno
     CHARACTER (len=default_path_length),&
        POINTER, DIMENSION(:)             :: io_stack_filename
     ! for '@SET VAR value' and '${VAR}'
     ! table size
     INTEGER                              :: num_variables
     ! table entries
     CHARACTER (len=default_path_length), &
          POINTER, DIMENSION(:)       :: variable_name
     CHARACTER (len=default_path_length), &
          POINTER, DIMENSION(:)       :: variable_value
  END TYPE inpp_type

  PUBLIC :: inpp_type, create_inpp_type, release_inpp_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_inpp_types'

CONTAINS

! ****************************************************************************
!> \brief creates the internal preprocessing type
!> \param inpp ...
!> \param initial_variables ...
!> \date  22.02.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_inpp_type(inpp, initial_variables)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(len=default_path_length), &
      DIMENSION(:, :), POINTER               :: initial_variables

    CHARACTER(len=*), PARAMETER :: routineN = 'create_inpp_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(inpp)))CALL cp__a("input/cp_parser_inpp_types.F",62)
    ALLOCATE(inpp)

    inpp%io_stack_level = 0
    NULLIFY(inpp%io_stack_channel,&
            inpp%io_stack_lineno,&
            inpp%io_stack_filename)

    inpp%num_variables =0
    NULLIFY(inpp%variable_name,&
            inpp%variable_value)

    IF (ASSOCIATED(initial_variables)) THEN
       inpp%num_variables =SIZE(initial_variables,2)
       ALLOCATE(inpp%variable_name(inpp%num_variables))
       inpp%variable_name=initial_variables(1,:)
       ALLOCATE(inpp%variable_value(inpp%num_variables))
       inpp%variable_value=initial_variables(2,:)
    ENDIF

  END SUBROUTINE create_inpp_type

! ****************************************************************************
!> \brief releases the internal preprocessing type
!> \param inpp ...
!> \date  22.02.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE release_inpp_type(inpp)
    TYPE(inpp_type), POINTER                 :: inpp

    CHARACTER(len=*), PARAMETER :: routineN = 'release_inpp_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(inpp)))CALL cp__a("input/cp_parser_inpp_types.F",96)

    IF (ASSOCIATED(inpp%io_stack_channel)) THEN
       DEALLOCATE(inpp%io_stack_channel)
    END IF
    IF (ASSOCIATED(inpp%io_stack_lineno)) THEN
       DEALLOCATE(inpp%io_stack_lineno)
    END IF
    IF (ASSOCIATED(inpp%io_stack_filename)) THEN
       DEALLOCATE(inpp%io_stack_filename)
    END IF

    IF (ASSOCIATED(inpp%variable_name)) THEN
       DEALLOCATE(inpp%variable_name)
    END IF
    IF (ASSOCIATED(inpp%variable_value)) THEN
       DEALLOCATE(inpp%variable_value)
    END IF

    DEALLOCATE(inpp)
  END SUBROUTINE release_inpp_type

END MODULE cp_parser_inpp_types
