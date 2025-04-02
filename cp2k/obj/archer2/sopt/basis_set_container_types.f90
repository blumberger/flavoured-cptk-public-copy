# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/aobasis/basis_set_container_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/aobasis/basis_set_container_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - Container to hold basis sets 
!> \author JGH (09.07.2015)
! *****************************************************************************
MODULE basis_set_container_types

  USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                             gto_basis_set_type
  USE kinds,                           ONLY: default_string_length

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/aobasis/../base/base_uses.f90" 1
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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/aobasis/basis_set_container_types.F" 2

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'basis_set_container_types'

! *****************************************************************************
  INTEGER, PARAMETER                       :: unknown_basis        = 100,&
                                              orbital_basis        = 101,&
                                              auxiliary_basis      = 102,&
                                              ri_aux_basis         = 103,&
                                              lri_aux_basis        = 104,&
                                              aux_fit_basis        = 105,&
                                              soft_basis           = 106,&
                                              hard_basis           = 107
! *****************************************************************************
  TYPE basis_set_container_type
     PRIVATE
     CHARACTER(LEN=default_string_length)       :: basis_type = ""
     INTEGER                                    :: basis_type_nr = 0
     TYPE(gto_basis_set_type), POINTER          :: basis_set => NULL()
  END TYPE basis_set_container_type
! *****************************************************************************

  PUBLIC :: basis_set_container_type

  PUBLIC :: remove_basis_set_container,&
            add_basis_set_to_container, get_basis_from_container,&
            remove_basis_from_container, translate_basis_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param basis ...
! *****************************************************************************
  SUBROUTINE remove_basis_set_container(basis)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: basis

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_basis_set_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(basis)
       basis(i)%basis_type = ""
       basis(i)%basis_type_nr = 0
       IF(ASSOCIATED(basis(i)%basis_set)) THEN
          CALL deallocate_gto_basis_set(basis(i)%basis_set)
       END IF
    END DO

  END SUBROUTINE remove_basis_set_container

! *****************************************************************************
!> \brief ...
!> \param basis_set_type ...
!> \retval basis_type_nr ...
! *****************************************************************************
  FUNCTION get_basis_type(basis_set_type) RESULT(basis_type_nr)
    CHARACTER(len=*)                         :: basis_set_type
    INTEGER                                  :: basis_type_nr

    SELECT CASE (basis_set_type)
       CASE ("ORB")
          basis_type_nr = orbital_basis
       CASE ("AUX")
          basis_type_nr = auxiliary_basis
       CASE ("RI_AUX")
          basis_type_nr = ri_aux_basis
       CASE ("LRI")
          basis_type_nr = lri_aux_basis
       CASE ("AUX_FIT")
          basis_type_nr = aux_fit_basis
       CASE ("SOFT")
          basis_type_nr = soft_basis
       CASE ("HARD")
          basis_type_nr = hard_basis
       CASE DEFAULT
          basis_type_nr = unknown_basis
    END SELECT

  END FUNCTION get_basis_type
! *****************************************************************************
!> \brief ...
!> \param basis_set_id ...
!> \retval basis_set_type ...
! *****************************************************************************
  FUNCTION translate_basis_type(basis_set_id) RESULT(basis_set_type)
    INTEGER                                  :: basis_set_id
    CHARACTER(len=default_string_length)     :: basis_set_type

    INTEGER, PARAMETER :: use_aux_basis_set = 3, use_aux_fit_basis_set = 2, &
      use_lri_basis_set = 5, use_orb_basis_set = 1, use_ri_aux_basis_set = 4

! This parameter identifies basis sets

    SELECT CASE (basis_set_id)
       CASE (use_orb_basis_set)
          basis_set_type = "ORB"
       CASE (use_ri_aux_basis_set)
          basis_set_type = "RI_AUX"
       CASE (use_aux_basis_set)
          basis_set_type = "AUX"
       CASE (use_aux_fit_basis_set)
          basis_set_type = "AUX_FIT"
       CASE (use_lri_basis_set)
          basis_set_type = "LRI"
       CASE DEFAULT
          basis_set_type = ""
    END SELECT

  END FUNCTION translate_basis_type

! *****************************************************************************
!> \brief ...
!> \param container ...
!> \param basis_set ...
!> \param basis_set_type ...
! *****************************************************************************
  SUBROUTINE add_basis_set_to_container(container,basis_set,basis_set_type)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    CHARACTER(len=*)                         :: basis_set_type

    CHARACTER(len=*), PARAMETER :: routineN = 'add_basis_set_to_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: success

    success = .FALSE.
    DO i=1,SIZE(container)
       IF(container(i)%basis_type_nr == 0) THEN
          container(i)%basis_type = basis_set_type
          container(i)%basis_set => basis_set
          container(i)%basis_type_nr = get_basis_type(basis_set_type)
          success = .TRUE.
          EXIT
       END IF
    END DO
    IF(.NOT.(success))CALL cp__a("aobasis/basis_set_container_types.F",165)

  END SUBROUTINE add_basis_set_to_container

! *****************************************************************************
!> \brief ...
!> \param container ...
!> \param inum ...
!> \param basis_type ...
! *****************************************************************************
  SUBROUTINE remove_basis_from_container(container,inum,basis_type)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    INTEGER, INTENT(IN), OPTIONAL            :: inum
    CHARACTER(len=*), OPTIONAL               :: basis_type

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_basis_from_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: basis_nr, i, ibas

    IF(PRESENT(inum)) THEN
       IF(.NOT.(inum <= SIZE(container)))CALL cp__a("aobasis/basis_set_container_types.F",187)
       IF(.NOT.(inum >= 1))CALL cp__a("aobasis/basis_set_container_types.F",188)
       ibas = inum
    ELSE IF(PRESENT(basis_type)) THEN
       basis_nr = get_basis_type(basis_type)
       ibas = 0
       DO i=1,SIZE(container)
          IF(container(i)%basis_type_nr == basis_nr) THEN
             ibas = i
             EXIT
          END IF
       END DO
    ELSE
       CALL cp__b("aobasis/basis_set_container_types.F",200,"")
    END IF
    !
    IF(ibas /= 0) THEN
       container(ibas)%basis_type = ""
       container(ibas)%basis_type_nr = 0
       IF(ASSOCIATED(container(ibas)%basis_set)) THEN
          CALL deallocate_gto_basis_set(container(ibas)%basis_set)
       END IF
    END IF
    ! shift other basis sets
    DO i=ibas+1,SIZE(container)
       IF(container(i)%basis_type_nr == 0) CYCLE
       container(i-1)%basis_type = container(i)%basis_type
       container(i-1)%basis_set => container(i)%basis_set
       container(i-1)%basis_type_nr = container(i)%basis_type_nr
    END DO

  END SUBROUTINE remove_basis_from_container

! *****************************************************************************
!> \brief Retrieve a basis set from the container
!> \param container ...
!> \param basis_set ...
!> \param inumbas ...
!> \param basis_type ...
! *****************************************************************************
  SUBROUTINE get_basis_from_container(container,basis_set,inumbas,basis_type)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    INTEGER, OPTIONAL                        :: inumbas
    CHARACTER(len=*), OPTIONAL               :: basis_type

    CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_from_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: basis_nr, i

    IF(PRESENT(inumbas)) THEN
       IF(.NOT.(inumbas <= SIZE(container)))CALL cp__a("aobasis/basis_set_container_types.F",240)
       IF(.NOT.(inumbas >= 1))CALL cp__a("aobasis/basis_set_container_types.F",241)
       basis_set => container(inumbas)%basis_set
       IF(PRESENT(basis_type)) THEN
          basis_type = container(inumbas)%basis_type
       END IF
    ELSE IF(PRESENT(basis_type)) THEN
       NULLIFY(basis_set)
       basis_nr = get_basis_type(basis_type)
       DO i=1,SIZE(container)
          IF(container(i)%basis_type_nr == basis_nr) THEN
             basis_set => container(i)%basis_set
             EXIT
          END IF
       END DO
    ELSE
       CALL cp__b("aobasis/basis_set_container_types.F",256,"")
    END IF

  END SUBROUTINE get_basis_from_container
! *****************************************************************************

END MODULE basis_set_container_types
