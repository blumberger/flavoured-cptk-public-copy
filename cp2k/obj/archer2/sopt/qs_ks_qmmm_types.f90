# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_ks_qmmm_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_ks_qmmm_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      05.2004 [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qs_ks_qmmm_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_p_type

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
# 24 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_ks_qmmm_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_qmmm_types'

  PUBLIC :: qs_ks_qmmm_env_type
  PUBLIC :: qs_ks_qmmm_release, qs_ks_qmmm_retain

! *****************************************************************************
!> \brief calculation environement to calculate the ks_qmmm matrix,
!>      holds the QM/MM potential and all the needed variables to
!>      compute the QM/MM electrostatic 1-electron ks matrix
!>      assumes that the core hamiltonian and energy are up to date.
!>      v_metal_rspace is the potential at the metal sites within the image
!>      charge approach
!> \par History
!>      05.2004 created [tlaino]
!>      01.2012 added v_metal_rspace [dgolze]
!> \author Teodoro Laino
! *****************************************************************************
  TYPE qs_ks_qmmm_env_type
     INTEGER :: n_evals, &
                id_nr, ref_count
     REAL(KIND=dp)                               :: pc_ener
     TYPE(pw_env_type), POINTER                  :: pw_env
     TYPE(pw_p_type)                             :: v_qmmm_rspace
     TYPE(pw_p_type)                             :: v_metal_rspace
     TYPE(cube_info_type),DIMENSION(:), POINTER  :: cube_info
     TYPE(cp_dbcsr_p_type), DIMENSION(:), &
          POINTER                                :: matrix_h
  END TYPE qs_ks_qmmm_env_type

! *****************************************************************************
!> \brief type to build arrays of pointers
!> \param ks_env the ks_env pointer
!> \par History
!>      05.2004 [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  TYPE qs_ks_qmmm_env_p_type
     TYPE(qs_ks_qmmm_env_type), POINTER :: ks_env
  END TYPE qs_ks_qmmm_env_p_type
CONTAINS

! *****************************************************************************
!> \brief releases the ks_qmmm_env (see doc/ReferenceCounting.html)
!> \param ks_qmmm_env the ks_qmmm_env to be released
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qs_ks_qmmm_release(ks_qmmm_env)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(pw_pool_type), POINTER              :: pool

    IF (ASSOCIATED(ks_qmmm_env)) THEN
       IF(.NOT.(ks_qmmm_env%ref_count>0))CALL cp__a("qs_ks_qmmm_types.F",87)
       ks_qmmm_env%ref_count=ks_qmmm_env%ref_count-1

       IF (ks_qmmm_env%ref_count<1) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=pool)
          CALL pw_pool_give_back_pw(pool,ks_qmmm_env%v_qmmm_rspace%pw)
          CALL pw_env_release(ks_qmmm_env%pw_env)
          IF (ASSOCIATED(ks_qmmm_env%cube_info))THEN
             DO i=1,SIZE(ks_qmmm_env%cube_info)
                CALL destroy_cube_info(ks_qmmm_env%cube_info(i))
             END DO
             DEALLOCATE(ks_qmmm_env%cube_info)
          END IF
          IF (ASSOCIATED(ks_qmmm_env%matrix_h)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(ks_qmmm_env%matrix_h)
          END IF
          DEALLOCATE(ks_qmmm_env)
       END IF
    END IF
    NULLIFY(ks_qmmm_env)
  END SUBROUTINE qs_ks_qmmm_release

! *****************************************************************************
!> \brief retains the given ks_environment
!> \param ks_qmmm_env the KohnSham QM/MM environment to retain
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
SUBROUTINE qs_ks_qmmm_retain(ks_qmmm_env)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_retain', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(ks_qmmm_env)))CALL cp__a("qs_ks_qmmm_types.F",122)
  IF(.NOT.(ks_qmmm_env%ref_count>0))CALL cp__a("qs_ks_qmmm_types.F",123)
  ks_qmmm_env%ref_count=ks_qmmm_env%ref_count+1
END SUBROUTINE qs_ks_qmmm_retain

END MODULE qs_ks_qmmm_types
