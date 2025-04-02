# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/input_enumeration_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/input_enumeration_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represents an enumeration, i.e. a mapping between integers and strings
!> \par History
!>      08.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_enumeration_types
  
  USE cp_log_handling,                 ONLY: cp_to_string
  USE kinds,                           ONLY: default_string_length
  USE string_utilities,                ONLY: a2s,&
                                             uppercase

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input/input_enumeration_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_enumeration_types'
  INTEGER, SAVE, PRIVATE :: last_enumeration_id=0

  PUBLIC :: enumeration_type
  PUBLIC :: enum_create, enum_retain, enum_release, enum_i2c, enum_c2i

! *****************************************************************************
!> \brief represents an enumaration, i.e. a mapping between strings and numbers
!> \param id_nr identification number (unique)
!> \param ref_count reference count
!> \param c_vals string values
!> \param i_vals integer values
!> \param strict if integer values not in the list should be accepted
!> \author fawzi
! *****************************************************************************
  TYPE char_array
     CHARACTER, DIMENSION(:), POINTER :: chars => Null()
  END TYPE char_array

  TYPE enumeration_type
     INTEGER :: id_nr, ref_count
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: c_vals
     TYPE(char_array), DIMENSION(:), POINTER :: desc => Null()
     INTEGER, DIMENSION(:), POINTER :: i_vals
     LOGICAL :: strict
  END TYPE enumeration_type

CONTAINS

! *****************************************************************************
!> \brief creates an enumeration
!> \param enum the enumeration to be created
!> \param c_vals string values
!> \param i_vals integer values
!> \param desc ...
!> \param strict if integer values not in the list should be accepted,
!>        defaults defaults to true
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_create(enum,c_vals,i_vals,desc,strict)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: c_vals
    INTEGER, DIMENSION(:), INTENT(in)        :: i_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: desc
    LOGICAL, INTENT(in), OPTIONAL            :: strict

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, n

  IF(.NOT.(.NOT.ASSOCIATED(enum)))CALL cp__a("input/input_enumeration_types.F",77)
  IF(.NOT.(SIZE(c_vals)==SIZE(i_vals)))CALL cp__a("input/input_enumeration_types.F",78)
  ALLOCATE(enum)
  last_enumeration_id=last_enumeration_id+1
  enum%id_nr=last_enumeration_id
  enum%ref_count=1
  ALLOCATE(enum%c_vals(SIZE(c_vals)))
  DO i=1,SIZE(enum%c_vals)
     enum%c_vals(i)=c_vals(i)
     CALL uppercase(enum%c_vals(i))
  END DO
  ALLOCATE(enum%i_vals(SIZE(i_vals)))
  enum%i_vals=i_vals
  enum%strict=.TRUE.
  IF (PRESENT(strict)) enum%strict=strict
  ALLOCATE(enum%desc(SIZE(c_vals)))
  IF (PRESENT(desc)) THEN
     IF(.NOT.(SIZE(enum%desc)==SIZE(desc)))CALL cp__a("input/input_enumeration_types.F",94)
     DO i=1,SIZE(enum%desc)
        n = LEN_TRIM(desc(i))
        ALLOCATE(enum%desc(i)%chars(n))
        DO j=1, n
          enum%desc(i)%chars(j) = desc(i)(j:j)
        END DO
    END DO
  ELSE
     DO i=1,SIZE(enum%desc)
        ALLOCATE(enum%desc(i)%chars(1))
        enum%desc(i)%chars(1:1) = ' '
     END DO
  END IF
END SUBROUTINE enum_create

! *****************************************************************************
!> \brief retains the given enumeration
!> \param enum the obect to retain
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_retain(enum)
    TYPE(enumeration_type), POINTER          :: enum

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_retain', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(enum)))CALL cp__a("input/input_enumeration_types.F",121)
  IF(.NOT.(enum%ref_count>0))CALL cp__a("input/input_enumeration_types.F",122)
  enum%ref_count=enum%ref_count+1
  END SUBROUTINE enum_retain

! *****************************************************************************
!> \brief releases the given enumeration
!> \param enum the obect to release
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_release(enum)
    TYPE(enumeration_type), POINTER          :: enum

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

  IF (ASSOCIATED(enum)) THEN
       IF(.NOT.(enum%ref_count>0))CALL cp__a("input/input_enumeration_types.F",140)
       enum%ref_count=enum%ref_count-1
       IF (enum%ref_count==0) THEN
          DEALLOCATE(enum%c_vals)
          DEALLOCATE(enum%i_vals)
          DO i=1, SIZE(enum%desc)
             DEALLOCATE(enum%desc(i)%chars)
          END DO
          DEALLOCATE(enum%desc)
          DEALLOCATE(enum)
       END IF
    END IF
    NULLIFY(enum)
  END SUBROUTINE enum_release

! *****************************************************************************
!> \brief maps an integer to a string
!> \param enum the enumeration to use for the mapping
!> \param i the value to map
!> \retval res ...
!> \author fawzi
! *****************************************************************************
FUNCTION enum_i2c(enum,i) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    INTEGER, INTENT(in)                      :: i
    CHARACTER(len=default_string_length)     :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_i2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j
    LOGICAL                                  :: found

  IF(.NOT.(ASSOCIATED(enum)))CALL cp__a("input/input_enumeration_types.F",173)
  IF(.NOT.(enum%ref_count>0))CALL cp__a("input/input_enumeration_types.F",174)
  res=" "
  found=.FALSE.
  DO j=1,SIZE(enum%i_vals)
     IF (enum%i_vals(j)==i) THEN
        res=enum%c_vals(j)
        found=.TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.found) THEN
     IF (enum%strict) THEN
      DO j=1,SIZE(enum%desc)
       PRINT *, TRIM(a2s(enum%desc(j)%chars))
       PRINT *, TRIM(enum%c_vals(j))
      ENDDO
       PRINT *, enum%i_vals
     END IF
     IF(enum%strict)&
        CALL cp__b("input/input_enumeration_types.F",193,"invalid value for enumeration:"//cp_to_string(i))
     res=ADJUSTL(cp_to_string(i))
  END IF
END FUNCTION enum_i2c

! *****************************************************************************
!> \brief maps a string to an integer
!> \param enum the enumeration to use for the mapping
!> \param c the value to map
!> \retval res ...
!> \author fawzi
! *****************************************************************************
FUNCTION enum_c2i(enum,c) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), INTENT(in)             :: c
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_c2i', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: upc
    INTEGER                                  :: iostat, j
    LOGICAL                                  :: found

  IF(.NOT.(ASSOCIATED(enum)))CALL cp__a("input/input_enumeration_types.F",217)
  IF(.NOT.(enum%ref_count>0))CALL cp__a("input/input_enumeration_types.F",218)
  upc=c
  CALL uppercase(upc)
  found=.FALSE.
  DO j=1,SIZE(enum%c_vals)
     IF (enum%c_vals(j)==upc) THEN
        res=enum%i_vals(j)
        found=.TRUE.
        EXIT
     END IF
  END DO

  IF (.NOT.found) THEN
     IF(enum%strict)&
        CALL cp__b("input/input_enumeration_types.F",232,"invalid value for enumeration:"//TRIM(c))
     READ(c,"(i10)",iostat=iostat) res
     IF(iostat/=0)&
        CALL cp__b("input/input_enumeration_types.F",235,"invalid value for enumeration2:"//TRIM(c))
  END IF
END FUNCTION enum_c2i

END MODULE input_enumeration_types
