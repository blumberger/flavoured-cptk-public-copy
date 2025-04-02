# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_fb_trial_fns_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_fb_trial_fns_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_trial_fns_types


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
# 9 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_fb_trial_fns_types.F" 2
  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_trial_fns_obj

! public methods
  PUBLIC :: fb_trial_fns_retain,&
            fb_trial_fns_release,&
            fb_trial_fns_nullify,&
            fb_trial_fns_associate,&
            fb_trial_fns_has_data,&
            fb_trial_fns_create,&
            fb_trial_fns_get,&
            fb_trial_fns_set

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_trial_fns_types'
  INTEGER, PRIVATE, SAVE :: last_fb_trial_fns_id = 0

! *****************************************************************************
!> \brief data containing information on trial functions used by filter
!>        matrix diagonalisation method
!> \param nfunctions : nfunctions(ikind) = number of trial functions for
!>                     atomic kind ikind
!> \param functions  : functions(itrial,ikind) = the index of the
!>                     GTO atomic orbital corresponding to itrial-th trial
!>                     function for kind ikind
!> \param id_nr      : unique id for the object
!> \param ref_count  : reference counter for the object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_trial_fns_data
     INTEGER :: id_nr, ref_count
     INTEGER, DIMENSION(:), POINTER :: nfunctions
     INTEGER, DIMENSION(:,:), POINTER :: functions
  END TYPE fb_trial_fns_data


! *****************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to fb_trial_fns objects
!> \param obj : pointer to the fb_trial_fns object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_trial_fns_obj
     TYPE(fb_trial_fns_data), POINTER, PRIVATE :: obj
  END TYPE fb_trial_fns_obj


CONTAINS


! *****************************************************************************
!> \brief retains given object
!> \brief ...
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_retain(trial_fns)
    ! note INTENT(IN) is okay because the obj pointer contained in the
    ! obj type will not be changed
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_retain', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(trial_fns%obj)))CALL cp__a("qs_fb_trial_fns_types.F",76)
    IF(.NOT.(trial_fns%obj%ref_count>0))CALL cp__a("qs_fb_trial_fns_types.F",77)
    trial_fns%obj%ref_count = trial_fns%obj%ref_count + 1
  END SUBROUTINE fb_trial_fns_retain


! *****************************************************************************
!> \brief releases given object
!> \brief ...
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_release(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(trial_fns%obj)) THEN
       IF(.NOT.(trial_fns%obj%ref_count>0))CALL cp__a("qs_fb_trial_fns_types.F",95)
       trial_fns%obj%ref_count = trial_fns%obj%ref_count - 1
       IF (trial_fns%obj%ref_count == 0) THEN
          trial_fns%obj%ref_count = 1
          IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
             DEALLOCATE(trial_fns%obj%nfunctions)
          END IF
          IF (ASSOCIATED(trial_fns%obj%functions)) THEN
             DEALLOCATE(trial_fns%obj%functions)
          END IF
          trial_fns%obj%ref_count = 0
          DEALLOCATE(trial_fns%obj)
       END IF
    ELSE
       NULLIFY(trial_fns%obj)
    END IF
  END SUBROUTINE fb_trial_fns_release


! *****************************************************************************
!> \brief nullifies the content of given object
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_nullify(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(trial_fns%obj)
  END SUBROUTINE fb_trial_fns_nullify


! *****************************************************************************
!> \brief associates the content of an object to that of another object
!>        of the same type
!> \param a : the output object
!> \param b : the input object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_associate(a, b)
    TYPE(fb_trial_fns_obj), INTENT(OUT)      :: a
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_trial_fns_associate


! *****************************************************************************
!> \brief check if the object has data associated to it
!> \param trial_fns : the fb_trial_fns object in question
!> \retval res : true if trial_fns%obj is associated, false otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_trial_fns_has_data(trial_fns) RESULT(res)
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(trial_fns%obj)
  END FUNCTION fb_trial_fns_has_data


! *****************************************************************************
!> \brief creates an fb_trial_fns object and initialises it
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_create(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(trial_fns%obj)))CALL cp__a("qs_fb_trial_fns_types.F",175)
    ALLOCATE(trial_fns%obj)
    NULLIFY(trial_fns%obj%nfunctions)
    NULLIFY(trial_fns%obj%functions)
    trial_fns%obj%ref_count = 1
    trial_fns%obj%id_nr = last_fb_trial_fns_id + 1
    last_fb_trial_fns_id = trial_fns%obj%id_nr
  END SUBROUTINE fb_trial_fns_create


! *****************************************************************************
!> \brief initialises an fb_trial_fns object
!> \param trial_fns : the fb_trial_fns object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_init(trial_fns)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_init', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(trial_fns%obj)))CALL cp__a("qs_fb_trial_fns_types.F",196)
    ! if halo_atoms are associated, then deallocate and de-associate
    IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
       DEALLOCATE(trial_fns%obj%nfunctions)
    END IF
    IF (ASSOCIATED(trial_fns%obj%functions)) THEN
       DEALLOCATE(trial_fns%obj%functions)
    END IF
  END SUBROUTINE fb_trial_fns_init


! *****************************************************************************
!> \brief get values of the attributes of a fb_trial_fns object
!> \param trial_fns  : the fb_trial_fns object in question
!> \param nfunctions : outputs pointer to trial_fns%obj%nfunctions
!> \param functions  : outputs pointer to trial_fns%obj%functions
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_get(trial_fns, &
                              nfunctions, &
                              functions)
    TYPE(fb_trial_fns_obj), INTENT(IN)       :: trial_fns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nfunctions
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: functions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_get', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(trial_fns%obj)))CALL cp__a("qs_fb_trial_fns_types.F",225)
    IF (PRESENT(nfunctions)) nfunctions => trial_fns%obj%nfunctions
    IF (PRESENT(functions)) functions => trial_fns%obj%functions
  END SUBROUTINE fb_trial_fns_get


! *****************************************************************************
!> \brief sets the attributes of a fb_trial_fns object
!> \param trial_fns  : the fb_trial_fns object in question
!> \param nfunctions : associates trial_fns%obj%nfunctions to this pointer
!> \param functions  : associates trial_fns%obj%nfunctions to this pointer
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_trial_fns_set(trial_fns, &
                              nfunctions, &
                              functions)
    TYPE(fb_trial_fns_obj), INTENT(INOUT)    :: trial_fns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nfunctions
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: functions

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_trial_fns_set', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(trial_fns%obj)))CALL cp__a("qs_fb_trial_fns_types.F",249)
    IF (PRESENT(nfunctions)) THEN
       IF (ASSOCIATED(trial_fns%obj%nfunctions)) THEN
          DEALLOCATE(trial_fns%obj%nfunctions)
       END IF
       trial_fns%obj%nfunctions => nfunctions
    END IF
    IF (PRESENT(functions)) THEN
       IF (ASSOCIATED(trial_fns%obj%functions)) THEN
          DEALLOCATE(trial_fns%obj%functions)
       END IF
       trial_fns%obj%functions => functions
    END IF
  END SUBROUTINE fb_trial_fns_set

END MODULE qs_fb_trial_fns_types
