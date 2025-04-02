# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cp_result_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cp_result_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  set of type/routines to handle the storage of results in force_envs
!> \author fschiff (12.2007)
!> \par    History
!>         - 10.2008 Teodoro Laino [tlaino] - University of Zurich
!>                   major rewriting:
!>                   - information stored in a proper type (not in a character!)
!>                   - module more lean
!>                   - splitting types and creating methods for cp_results
! *****************************************************************************
  MODULE cp_result_types

  USE kinds,                           ONLY: default_string_length,&
                                             dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/../base/base_uses.f90" 1
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
# 21 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/cp_result_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_result_types'

  INTEGER, PARAMETER, PUBLIC :: result_type_logical=1,&
                                result_type_integer=2,&
                                result_type_real=3

! *** Public data types ***
  PUBLIC :: cp_result_type,&
            cp_result_p_type

! *** Public subroutines ***
  PUBLIC :: cp_result_create,&
            cp_result_release,&
            cp_result_retain,&
            cp_result_clean,&
            cp_result_copy,&
            cp_result_value_create,&
            cp_result_value_copy,&
            cp_result_value_p_reallocate,&
            cp_result_value_init

! *****************************************************************************
!> \brief low level type for storing real informations
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  TYPE cp_result_value_type
     INTEGER                                              :: type_in_use
     LOGICAL, DIMENSION(:), POINTER                       :: logical_type
     INTEGER, DIMENSION(:), POINTER                       :: integer_type
     REAL(KIND=dp), DIMENSION(:), POINTER                 :: real_type
  END TYPE cp_result_value_type

! *****************************************************************************
  TYPE cp_result_value_p_type
     TYPE(cp_result_value_type), POINTER                  :: value
  END TYPE cp_result_value_p_type

! *****************************************************************************
!> \brief contains arbitrary information which need to be stored
!> \note
!>      result_list is a character list, in which everthing can be stored
!>      before passing any variable just name the variable like '[NAME]'
!>      brackets will be used to identify the start of a new set
!> \author fschiff (12.2007)
! *****************************************************************************
  TYPE cp_result_type
     INTEGER                                              :: ref_count
     TYPE(cp_result_value_p_type), POINTER, DIMENSION(:)  :: result_value
     CHARACTER(LEN=default_string_length),DIMENSION(:),&
          POINTER                                         :: result_label
  END TYPE cp_result_type

! *****************************************************************************
  TYPE cp_result_p_type
     TYPE(cp_result_type), POINTER                        :: results
  END TYPE cp_result_p_type

CONTAINS

! *****************************************************************************
!> \brief Allocates and intitializes the cp_result
!> \param results ...
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_result_create(results)
    TYPE(cp_result_type), POINTER            :: results

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    ALLOCATE(results)
    NULLIFY(results%result_value, results%result_label)
    results%ref_count=1
    ALLOCATE(results%result_label(0))
    ALLOCATE(results%result_value(0))
    CALL timestop(handle)
  END SUBROUTINE cp_result_create

! *****************************************************************************
!> \brief Releases cp_result type
!> \param results ...
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_result_release(results)
    TYPE(cp_result_type), POINTER            :: results

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)
    IF(ASSOCIATED(results))THEN
       IF(.NOT.(results%ref_count>0))CALL cp__a("common/cp_result_types.F",128)
       results%ref_count=results%ref_count-1
       IF (results%ref_count==0) THEN
          ! Description
          IF(ASSOCIATED(results%result_label))THEN
             DEALLOCATE(results%result_label)
          END IF
          ! Values
          IF(ASSOCIATED(results%result_value))THEN
             DO i = 1, SIZE(results%result_value)
                CALL cp_result_value_release(results%result_value(i)%value)
             END DO
             DEALLOCATE(results%result_value)
          END IF
          DEALLOCATE(results)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_result_release

! *****************************************************************************
!> \brief Releases cp_result clean
!> \param results ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_clean(results)
    TYPE(cp_result_type), POINTER            :: results

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_clean', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)
    IF(ASSOCIATED(results))THEN
       ! Description
       IF(ASSOCIATED(results%result_label))THEN
          DEALLOCATE(results%result_label)
       END IF
       ! Values
       IF(ASSOCIATED(results%result_value))THEN
          DO i = 1, SIZE(results%result_value)
             CALL cp_result_value_release(results%result_value(i)%value)
          END DO
          DEALLOCATE(results%result_value)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_result_clean

! *****************************************************************************
!> \brief Retains cp_result type
!> \param results ...
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_result_retain(results)
    TYPE(cp_result_type), POINTER            :: results

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_retain', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(results)))CALL cp__a("common/cp_result_types.F",191)
    IF(.NOT.(results%ref_count>0))CALL cp__a("common/cp_result_types.F",192)
    results%ref_count=results%ref_count+1
  END SUBROUTINE cp_result_retain

! *****************************************************************************
!> \brief Allocates and intitializes the cp_result_value type
!> \param value ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_value_create(value)
    TYPE(cp_result_value_type), POINTER      :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_value_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    ALLOCATE(value)
    value%type_in_use = -1
    NULLIFY(value%real_type)
    NULLIFY(value%logical_type)
    NULLIFY(value%integer_type)
    CALL timestop(handle)
  END SUBROUTINE cp_result_value_create

! *****************************************************************************
!> \brief Setup of the cp_result_value type
!> \param value ...
!> \param type_in_use ...
!> \param size_value ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_value_init(value, type_in_use, size_value)
    TYPE(cp_result_value_type), POINTER      :: value
    INTEGER, INTENT(IN)                      :: type_in_use, size_value

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_value_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    IF(.NOT.(ASSOCIATED(value)))CALL cp__a("common/cp_result_types.F",235)
    value%type_in_use = type_in_use
    SELECT CASE(value%type_in_use)
    CASE(result_type_real)
       ALLOCATE(value%real_type(size_value))
    CASE(result_type_integer)
       ALLOCATE(value%integer_type(size_value))
    CASE(result_type_logical)
       ALLOCATE(value%logical_type(size_value))
    CASE DEFAULT
       ! Type not implemented in cp_result_type
       CALL cp__b("common/cp_result_types.F",246,"")
    END SELECT
    CALL timestop(handle)
  END SUBROUTINE cp_result_value_init

! *****************************************************************************
!> \brief Releases the cp_result_value type
!> \param value ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_value_release(value)
    TYPE(cp_result_value_type), POINTER      :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_value_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    IF (ASSOCIATED(value)) THEN
       SELECT CASE(value%type_in_use)
       CASE(result_type_real)
          IF (ASSOCIATED(value%real_type)) THEN
             DEALLOCATE(value%real_type)
          END IF
          IF(.NOT.(.NOT.ASSOCIATED(value%integer_type)))CALL cp__a("common/cp_result_types.F",271)
          IF(.NOT.(.NOT.ASSOCIATED(value%logical_type)))CALL cp__a("common/cp_result_types.F",272)
       CASE(result_type_integer)
          IF (ASSOCIATED(value%integer_type)) THEN
             DEALLOCATE(value%integer_type)
          END IF
          IF(.NOT.(.NOT.ASSOCIATED(value%real_type)))CALL cp__a("common/cp_result_types.F",277)
          IF(.NOT.(.NOT.ASSOCIATED(value%logical_type)))CALL cp__a("common/cp_result_types.F",278)
       CASE(result_type_logical)
          IF (ASSOCIATED(value%logical_type)) THEN
             DEALLOCATE(value%logical_type)
          END IF
          IF(.NOT.(.NOT.ASSOCIATED(value%integer_type)))CALL cp__a("common/cp_result_types.F",283)
          IF(.NOT.(.NOT.ASSOCIATED(value%real_type)))CALL cp__a("common/cp_result_types.F",284)
       CASE DEFAULT
          ! Type not implemented in cp_result_type
          CALL cp__b("common/cp_result_types.F",287,"")
       END SELECT
       DEALLOCATE(value)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_result_value_release

! *****************************************************************************
!> \brief Copies the cp_result type
!> \param results_in ...
!> \param results_out ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_copy(results_in, results_out)
    TYPE(cp_result_type), POINTER            :: results_in, results_out

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ndim
    LOGICAL                                  :: check

    CALL timeset(routineN,handle)
    IF(.NOT.(ASSOCIATED(results_in)))CALL cp__a("common/cp_result_types.F",310)
    IF(.NOT.(ASSOCIATED(results_out)))CALL cp__a("common/cp_result_types.F",311)
    CALL cp_result_clean(results_out)

    check = SIZE(results_in%result_label)==SIZE(results_in%result_value)
    IF(.NOT.(check))CALL cp__a("common/cp_result_types.F",315)
    ndim  = SIZE(results_in%result_value)
    ALLOCATE(results_out%result_label(ndim))
    ALLOCATE(results_out%result_value(ndim))
    DO i = 1, ndim
       results_out%result_label(i) = results_in%result_label(i)
       CALL cp_result_value_create(results_out%result_value(i)%value)
       CALL cp_result_value_copy(results_out%result_value(i)%value,&
            results_in%result_value(i)%value)
    END DO
    CALL timestop(handle)
  END SUBROUTINE cp_result_copy

! *****************************************************************************
!> \brief Copies the cp_result_value type
!> \param value_out ...
!> \param value_in ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_value_copy(value_out, value_in)
    TYPE(cp_result_value_type), POINTER      :: value_out, value_in

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_value_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isize

    CALL timeset(routineN,handle)
    IF(.NOT.(ASSOCIATED(value_in)))CALL cp__a("common/cp_result_types.F",343)
    IF(.NOT.(ASSOCIATED(value_out)))CALL cp__a("common/cp_result_types.F",344)
    value_out%type_in_use = value_in%type_in_use
    SELECT CASE(value_out%type_in_use)
    CASE(result_type_real)
       isize = SIZE(value_in%real_type)
       ALLOCATE(value_out%real_type(isize))
       value_out%real_type = value_in%real_type
    CASE(result_type_integer)
       isize = SIZE(value_in%integer_type)
       ALLOCATE(value_out%integer_type(isize))
       value_out%integer_type = value_in%integer_type
    CASE(result_type_logical)
       isize = SIZE(value_in%logical_type)
       ALLOCATE(value_out%logical_type(isize))
       value_out%logical_type = value_in%logical_type
    CASE DEFAULT
       ! Type not implemented in cp_result_type
       CALL cp__b("common/cp_result_types.F",361,"")
    END SELECT
    CALL timestop(handle)
  END SUBROUTINE cp_result_value_copy

! *****************************************************************************
!> \brief Reallocates the cp_result_value type
!> \param result_value ...
!> \param istart ...
!> \param iend ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE cp_result_value_p_reallocate(result_value, istart, iend)
    TYPE(cp_result_value_p_type), &
      DIMENSION(:), POINTER                  :: result_value
    INTEGER, INTENT(in)                      :: istart, iend

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_value_p_reallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb_size, ub_size
    TYPE(cp_result_value_p_type), &
      DIMENSION(:), POINTER                  :: tmp_value

    CALL timeset(routineN,handle)
    ub_size  = 0
    lb_size  = 0
    IF (ASSOCIATED(result_value)) THEN
       ub_size = UBOUND(result_value,1)
       lb_size = LBOUND(result_value,1)
    END IF
    ! Allocate and copy new values while releases old
    ALLOCATE(tmp_value(istart:iend))
    DO i = istart, iend
       NULLIFY(tmp_value(i)%value)
       CALL cp_result_value_create(tmp_value(i)%value)
       IF ((i<=ub_size).AND.(i>=lb_size)) THEN
          CALL cp_result_value_copy(tmp_value(i)%value, result_value(i)%value)
          CALL cp_result_value_release(result_value(i)%value)
       END IF
    END DO
    DEALLOCATE(result_value)
    result_value => tmp_value
    CALL timestop(handle)
  END SUBROUTINE cp_result_value_p_reallocate

END MODULE cp_result_types
