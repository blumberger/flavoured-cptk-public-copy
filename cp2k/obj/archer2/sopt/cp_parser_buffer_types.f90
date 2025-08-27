# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_buffer_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_buffer_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple buffering of read lines of a parser
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 08.2008
! *****************************************************************************
MODULE cp_parser_buffer_types
  
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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/cp_parser_buffer_types.F" 2

  IMPLICIT NONE
  PRIVATE

! ****************************************************************************
!> \brief  Buffer type for speeding-up the parsing in parallel
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  TYPE buffer_type
     INTEGER                              :: size, buffer_id
     INTEGER                              :: present_line_number,&
                                             last_line_number,&
                                             istat
     INTEGER, DIMENSION(:), POINTER       :: input_line_numbers
     CHARACTER(LEN=max_line_length), &
          DIMENSION(:), POINTER           :: input_lines
     TYPE(buffer_type), POINTER           :: sub_buffer
  END TYPE buffer_type

  PUBLIC :: buffer_type, create_buffer_type, release_buffer_type, copy_buffer_type,&
            initialize_sub_buffer, finalize_sub_buffer
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_buffer_types'
  INTEGER, PARAMETER, PRIVATE          :: buffer_size=1000

CONTAINS

! ****************************************************************************
!> \brief  Creates the parser buffer type
!> \param buffer ...
!> \date   08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_buffer_type(buffer)
    TYPE(buffer_type), POINTER               :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'create_buffer_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(buffer)))CALL cp__a("input/cp_parser_buffer_types.F",54)
    ALLOCATE(buffer)
    buffer%size = buffer_size
    ALLOCATE(buffer%input_lines(buffer%size))
    ALLOCATE(buffer%input_line_numbers(buffer%size))
    buffer%buffer_id           = 0
    buffer%input_line_numbers  = 0
    buffer%istat               = 0
    buffer%present_line_number = buffer%size
    buffer%last_line_number    = buffer%size
    NULLIFY(buffer%sub_buffer)
  END SUBROUTINE create_buffer_type

! ****************************************************************************
!> \brief  Releases the parser buffer type
!> \param buffer ...
!> \date   08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  RECURSIVE SUBROUTINE release_buffer_type(buffer)
    TYPE(buffer_type), POINTER               :: buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'release_buffer_type', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(buffer)))CALL cp__a("input/cp_parser_buffer_types.F",79)
    DEALLOCATE(buffer%input_lines)
    DEALLOCATE(buffer%input_line_numbers)
    IF (ASSOCIATED(buffer%sub_buffer)) THEN
       CALL release_buffer_type(buffer%sub_buffer)
    END IF
    DEALLOCATE(buffer)
  END SUBROUTINE release_buffer_type

! ****************************************************************************
!> \brief  Copies  buffer types
!> \param buffer_in ...
!> \param buffer_out ...
!> \param force ...
!> \date   08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  RECURSIVE SUBROUTINE copy_buffer_type(buffer_in, buffer_out, force)
    TYPE(buffer_type), POINTER               :: buffer_in, buffer_out
    LOGICAL, INTENT(IN), OPTIONAL            :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_buffer_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: my_force

    IF(.NOT.(ASSOCIATED(buffer_in)))CALL cp__a("input/cp_parser_buffer_types.F",106)
    IF(.NOT.(ASSOCIATED(buffer_out)))CALL cp__a("input/cp_parser_buffer_types.F",107)
    IF(.NOT.(buffer_in%size==buffer_out%size))CALL cp__a("input/cp_parser_buffer_types.F",108)
    my_force = .FALSE.
    IF (PRESENT(force)) my_force = force
    ! Copy buffer structure
    buffer_out%present_line_number = buffer_in%present_line_number
    buffer_out%last_line_number    = buffer_in%last_line_number
    buffer_out%istat               = buffer_in%istat
    ! This part can be quite expensive.. we do it only when strictly necessary..
    IF ((buffer_out%buffer_id/=buffer_in%buffer_id).OR.(my_force)) THEN
       buffer_out%buffer_id           = buffer_in%buffer_id
       buffer_out%input_line_numbers  = buffer_in%input_line_numbers
       ! Explicit loop: bypass a NAG bug..
       DO i = 1, SIZE(buffer_in%input_lines)
          buffer_out%input_lines(i)   = buffer_in%input_lines(i)
       END DO
    END IF
    IF (ASSOCIATED(buffer_in%sub_buffer).AND.ASSOCIATED(buffer_out%sub_buffer)) THEN
       CALL copy_buffer_type(buffer_in%sub_buffer, buffer_out%sub_buffer, force)
    END IF
  END SUBROUTINE copy_buffer_type

! ****************************************************************************
!> \brief  Initializes sub buffer structure
!> \param sub_buffer ...
!> \param buffer ...
!> \date   08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE initialize_sub_buffer(sub_buffer, buffer)
    TYPE(buffer_type), POINTER               :: sub_buffer, buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_sub_buffer', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(buffer)))CALL cp__a("input/cp_parser_buffer_types.F",142)
    IF(.NOT.(.NOT.ASSOCIATED(sub_buffer)))CALL cp__a("input/cp_parser_buffer_types.F",143)
    CALL create_buffer_type(sub_buffer)
    CALL copy_buffer_type(buffer, sub_buffer)
    sub_buffer%present_line_number = 0
  END SUBROUTINE initialize_sub_buffer


! ****************************************************************************
!> \brief  Finalizes sub buffer structure
!> \param sub_buffer ...
!> \param buffer ...
!> \date   08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE finalize_sub_buffer(sub_buffer, buffer)
    TYPE(buffer_type), POINTER               :: sub_buffer, buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'finalize_sub_buffer', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(buffer)))CALL cp__a("input/cp_parser_buffer_types.F",163)
    IF(.NOT.(ASSOCIATED(sub_buffer)))CALL cp__a("input/cp_parser_buffer_types.F",164)
    CALL copy_buffer_type(sub_buffer,buffer)
    CALL release_buffer_type(sub_buffer)
  END SUBROUTINE finalize_sub_buffer

END MODULE cp_parser_buffer_types
