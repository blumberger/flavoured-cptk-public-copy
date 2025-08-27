# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmm_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmm_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Basic container type for QM/MM.
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmm_types
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_env_release,&
                                             fist_environment_type
  USE kinds,                           ONLY: dp
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type

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
# 24 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmm_types.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_types'

  PUBLIC :: qmmm_env_type, qmmm_env_retain, qmmm_env_release, qmmm_env_get

  TYPE qmmm_env_type
     INTEGER                                                 :: ref_count = 1
     TYPE(qs_environment_type), POINTER                      :: qs_env => Null()
     TYPE(fist_environment_type), POINTER                    :: fist_env => Null()
     TYPE(qmmm_env_qm_type), POINTER                         :: qm => Null()
  END TYPE qmmm_env_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
! *****************************************************************************
  SUBROUTINE qmmm_env_get(qmmm_env,subsys,potential_energy,kinetic_energy)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_get', &
      routineP = moduleN//':'//routineN

    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(qs_energy_type), POINTER            :: qs_energy

    NULLIFY(qs_energy, thermo)

    IF(.NOT.(ASSOCIATED(qmmm_env)))CALL cp__a("qmmm_types.F",62)
    IF(.NOT.(qmmm_env%ref_count>0))CALL cp__a("qmmm_types.F",63)

    IF (PRESENT(kinetic_energy)) THEN
        CALL fist_env_get(qmmm_env%fist_env,thermo=thermo)
        kinetic_energy = thermo%kin
    END IF
    IF (PRESENT(subsys)) THEN
        CALL fist_env_get(qmmm_env%fist_env,subsys=subsys)
    ENDIF
    IF (PRESENT(potential_energy)) THEN
         ! get the underlying energies from primary subsys.  This is the only subsys
         ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
         CALL fist_env_get(qmmm_env%fist_env, thermo=thermo)
         CALL get_qs_env(qmmm_env%qs_env,energy=qs_energy)
         potential_energy = thermo%pot + qs_energy%total
    ENDIF
  END SUBROUTINE qmmm_env_get


! *****************************************************************************
!> \brief releases the given qmmm_env (see doc/ReferenceCounting.html)
!> \param qmmm_env the object to release
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_release(qmmm_env)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(qmmm_env)) THEN
       IF(.NOT.(qmmm_env%ref_count>0))CALL cp__a("qmmm_types.F",94)
       qmmm_env%ref_count=qmmm_env%ref_count-1
       IF (qmmm_env%ref_count==0) THEN
          CALL qs_env_release(qmmm_env%qs_env)
          CALL qmmm_env_qm_release(qmmm_env%qm)
          CALL fist_env_release(qmmm_env%fist_env)
          DEALLOCATE(qmmm_env)
       END IF
    END IF
    NULLIFY(qmmm_env)
  END SUBROUTINE qmmm_env_release

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_retain(qmmm_env)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_retain', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(qmmm_env)))CALL cp__a("qmmm_types.F",117)
    IF(.NOT.(qmmm_env%ref_count>0))CALL cp__a("qmmm_types.F",118)
    qmmm_env%ref_count=qmmm_env%ref_count+1
  END SUBROUTINE qmmm_env_retain

END MODULE qmmm_types
