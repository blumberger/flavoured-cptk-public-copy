# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Basic container type for QM/MM with force mixing.
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_types
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE kinds,                           ONLY: dp
  USE qmmm_types,                      ONLY: qmmm_env_get,&
                                             qmmm_env_release,&
                                             qmmm_env_type

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qmmmx_types.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_types'

  PUBLIC :: qmmmx_env_type, qmmmx_env_get, qmmmx_env_retain, qmmmx_env_release

  TYPE qmmmx_env_type
     INTEGER                                                 :: ref_count = 1
     TYPE(qmmm_env_type), POINTER                            :: core => Null()
     TYPE(qmmm_env_type), POINTER                            :: ext => Null()
  END TYPE qmmmx_env_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
! *****************************************************************************
  SUBROUTINE qmmmx_env_get(qmmmx_env,subsys,potential_energy,kinetic_energy)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_get', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(qmmmx_env)))CALL cp__a("qmmmx_types.F",49)
    IF(.NOT.(qmmmx_env%ref_count>0))CALL cp__a("qmmmx_types.F",50)

    ! get the underlying energies from primary subsys.  This is the only subsys
    ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
    CALL qmmm_env_get(qmmmx_env%ext,&
                      kinetic_energy=kinetic_energy,&
                      potential_energy=potential_energy,&
                      subsys=subsys)

  END SUBROUTINE qmmmx_env_get


! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_retain(qmmmx_env)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_retain', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(qmmmx_env)))CALL cp__a("qmmmx_types.F",73)
    IF(.NOT.(qmmmx_env%ref_count>0))CALL cp__a("qmmmx_types.F",74)
    qmmmx_env%ref_count = qmmmx_env%ref_count+1
  END SUBROUTINE qmmmx_env_retain


! *****************************************************************************
!> \brief releases the given qmmmx_env (see doc/ReferenceCounting.html)
!> \param qmmmx_env the object to release
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_release(qmmmx_env)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(qmmmx_env)) THEN
       IF(.NOT.(qmmmx_env%ref_count>0))CALL cp__a("qmmmx_types.F",91)
       qmmmx_env%ref_count = qmmmx_env%ref_count-1
       IF (qmmmx_env%ref_count==0) THEN
          CALL qmmm_env_release(qmmmx_env%core)
          CALL qmmm_env_release(qmmmx_env%ext)
          DEALLOCATE(qmmmx_env)
       END IF
    END IF
    NULLIFY(qmmmx_env)
  END SUBROUTINE qmmmx_env_release

END MODULE qmmmx_types
