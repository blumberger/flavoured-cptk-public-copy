# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE manybody_quip

USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE bibliography,                    ONLY: QUIP_ref,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_set,&
                                             fist_nonbond_env_type,&
                                             quip_data_type
  USE kinds,                           ONLY: dp
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             quip_pot_type,&
                                             quip_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom,&
                                             evolt





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
# 31 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F" 2

IMPLICIT NONE

PRIVATE

PUBLIC quip_energy_store_force_virial, quip_add_force_virial

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_quip'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param fist_nonbond_env ...
!> \param pot_quip ...
!> \param para_env ...
! *****************************************************************************
SUBROUTINE quip_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, fist_nonbond_env, &
                                         pot_quip, para_env)
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    REAL(kind=dp)                            :: pot_quip
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'quip_energy_store_force_virial', &
      routineP = moduleN//':'//routineN

# 83 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(particle_set))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(cell))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(atomic_kind_set))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(potparm))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(fist_nonbond_env))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pot_quip))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(para_env))==-1) EXIT ;  END DO ; ENDIF
    CALL cp_abort(cp__l("manybody_quip.F",90),"In order to use QUIP you need to download "//&
         "and install the libAtoms/QUIP library (check CP2K manual for details)")
# 188 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F"
END SUBROUTINE quip_energy_store_force_virial

! *****************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param virial ...
! *****************************************************************************
SUBROUTINE quip_add_force_virial(fist_nonbond_env, force, virial)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    REAL(KIND=dp)                            :: force(:,:), virial(3,3)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'quip_add_force_virial', &
      routineP = moduleN//':'//routineN







    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(fist_nonbond_env))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(force))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(virial))==-1) EXIT ;  END DO ; ENDIF
    RETURN
# 224 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/manybody_quip.F"
END SUBROUTINE quip_add_force_virial

END MODULE manybody_quip
