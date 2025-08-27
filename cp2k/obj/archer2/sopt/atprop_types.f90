# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/subsys/atprop_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/subsys/atprop_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Holds information on atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE atprop_types
  
  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/subsys/../base/base_uses.f90" 1
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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/subsys/atprop_types.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atprop_types'

  PUBLIC :: atprop_type, atprop_create, atprop_release, atprop_init
  PUBLIC :: atprop_array_init, atprop_array_add

! *****************************************************************************
!> \brief type for the atomic properties
!> \author fawzi
! *****************************************************************************
TYPE atprop_type
   LOGICAL                                   :: energy
   LOGICAL                                   :: stress
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atener
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ateb
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atexc
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ateself
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atevdw
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atecc
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ate1c
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atecoul
   REAL(KIND=dp), DIMENSION(:,:,:), POINTER  :: atstress
END TYPE atprop_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atprop_env ...
! *****************************************************************************
  SUBROUTINE atprop_create(atprop_env)
    TYPE(atprop_type), POINTER               :: atprop_env

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_create', &
      routineP = moduleN//':'//routineN

  CALL atprop_release(atprop_env)
  ALLOCATE(atprop_env)
  NULLIFY(atprop_env%atener,atprop_env%atstress)
  NULLIFY(atprop_env%ateb,atprop_env%atevdw,atprop_env%atecc,atprop_env%atecoul)
  NULLIFY(atprop_env%ateself,atprop_env%atexc,atprop_env%ate1c)
  atprop_env%energy = .FALSE.
  atprop_env%stress = .FALSE.

  END SUBROUTINE atprop_create

! *****************************************************************************
!> \brief ...
!> \param atprop_env ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE atprop_init(atprop_env,natom)
    TYPE(atprop_type), POINTER               :: atprop_env
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_init', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(atprop_env)))CALL cp__a("subsys/atprop_types.F",77)

  IF(atprop_env%energy) THEN
    CALL atprop_array_init(atprop_env%atener,natom)
    CALL atprop_array_release(atprop_env%ateb)
    CALL atprop_array_release(atprop_env%atevdw)
    CALL atprop_array_release(atprop_env%atecc)
    CALL atprop_array_release(atprop_env%atecoul)
    CALL atprop_array_release(atprop_env%ateself)
    CALL atprop_array_release(atprop_env%atexc)
    CALL atprop_array_release(atprop_env%ate1c)
  END IF

  IF(atprop_env%stress) THEN
    IF(ASSOCIATED(atprop_env%atstress)) THEN
      IF(.NOT.(SIZE(atprop_env%atstress,3)==natom))CALL cp__a("subsys/atprop_types.F",92)
    ELSE
      ALLOCATE(atprop_env%atstress(3,3,natom))
    END IF
    atprop_env%atstress = 0._dp
  END IF

  END SUBROUTINE atprop_init

! *****************************************************************************
!> \brief ...
!> \param atarray ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE atprop_array_init(atarray,natom)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atarray
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_init', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(atarray)) THEN
      IF(.NOT.(SIZE(atarray)==natom))CALL cp__a("subsys/atprop_types.F",114)
    ELSE
      ALLOCATE(atarray(natom))
    END IF
    atarray = 0._dp

  END SUBROUTINE atprop_array_init

! *****************************************************************************
!> \brief ...
!> \param atarray ...
! *****************************************************************************
  SUBROUTINE atprop_array_release(atarray)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atarray

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(atarray)) THEN
      DEALLOCATE(atarray)
    END IF

  END SUBROUTINE atprop_array_release

! *****************************************************************************
!> \brief ...
!> \param array_a ...
!> \param array_b ...
! *****************************************************************************
  SUBROUTINE atprop_array_add(array_a,array_b)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: array_a, array_b

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_add', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(array_b)) THEN
      IF(.NOT.(ASSOCIATED(array_a)))CALL cp__a("subsys/atprop_types.F",150)
      array_a = array_a + array_b
    END IF

  END SUBROUTINE atprop_array_add

! *****************************************************************************
!> \brief releases the atprop
!> \param atprop_env the object to release
!> \author fawzi
! *****************************************************************************
SUBROUTINE atprop_release(atprop_env)
    TYPE(atprop_type), POINTER               :: atprop_env

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_release', &
      routineP = moduleN//':'//routineN

  IF (ASSOCIATED(atprop_env)) THEN
     ! energy
     CALL atprop_array_release(atprop_env%atener)
     CALL atprop_array_release(atprop_env%ateb)
     CALL atprop_array_release(atprop_env%ateself)
     CALL atprop_array_release(atprop_env%atexc)
     CALL atprop_array_release(atprop_env%atevdw)
     CALL atprop_array_release(atprop_env%atecc)
     CALL atprop_array_release(atprop_env%ate1c)
     CALL atprop_array_release(atprop_env%atecoul)
     ! stress
     IF (ASSOCIATED(atprop_env%atstress)) THEN
        DEALLOCATE(atprop_env%atstress)
     END IF
     ! atprop type
     DEALLOCATE(atprop_env)
  END IF
  NULLIFY(atprop_env)
END SUBROUTINE atprop_release

END MODULE atprop_types
