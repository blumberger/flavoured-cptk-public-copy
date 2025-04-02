# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/reftraj_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/reftraj_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialization of the reftraj structure used to analyse
!>     previously generated trajectories
!> \par History
!>      Created 10-07 [MI]
!> \author MI
! *****************************************************************************
MODULE reftraj_types

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/../base/base_uses.f90" 1
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
# 24 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/reftraj_types.F" 2

  IMPLICIT NONE

    PRIVATE
  PUBLIC :: reftraj_type, reftraj_msd_type,&
            create_reftraj, release_reftraj, retain_reftraj

! *****************************************************************************
!> \brief parameters related to the analysis of previously generated trajecorties
!> \author MI
! *****************************************************************************
  TYPE reftraj_info_type
     INTEGER                                  :: first_snapshot
     INTEGER                                  :: last_snapshot
     INTEGER                                  :: stride
     LOGICAL                                  :: eval_ef
     LOGICAL                                  :: variable_volume
     LOGICAL                                  :: msd
     TYPE(cp_parser_type), POINTER            :: traj_parser
     TYPE(cp_parser_type), POINTER            :: cell_parser
  END TYPE reftraj_info_type

! *****************************************************************************
  TYPE reftraj_msd_type
     LOGICAL                                  :: disp_atom, msd_kind, msd_molecule, msd_region
     INTEGER                                  :: num_disp_atom, ref0_unit
     INTEGER, POINTER, DIMENSION(:)           :: disp_atom_index
     REAL(KIND=dp)                            :: disp_atom_tol, drcom(3), ref0_com(3), total_mass
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: disp_atom_dr
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ref0_pos
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ref0_com_molecule
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_kind
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_molecule
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_region
  END TYPE reftraj_msd_type

! *****************************************************************************
  TYPE reftraj_type
     INTEGER                                  :: ref_count
     INTEGER                                  :: itimes
     INTEGER                                  :: itimes0
     INTEGER                                  :: isnap
     INTEGER                                  :: natom
     LOGICAL                                  :: init
     REAL(KIND=dp)                            :: epot, epot0, time, time0
     TYPE(reftraj_info_type), POINTER         :: info
     TYPE(reftraj_msd_type), POINTER          :: msd
  END TYPE reftraj_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_types'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param reftraj ...
!> \param reftraj_section ...
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE create_reftraj(reftraj,reftraj_section,para_env)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'create_reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename

    IF(.NOT.(.NOT. ASSOCIATED(reftraj)))CALL cp__a("motion/reftraj_types.F",94)
    ALLOCATE(reftraj)
    reftraj%ref_count = 1

    NULLIFY(reftraj%info)
    NULLIFY(reftraj%msd)

    ALLOCATE(reftraj%info)
    NULLIFY(reftraj%info%traj_parser)
    NULLIFY(reftraj%info%cell_parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(reftraj_section,"TRAJ_FILE_NAME",c_val=filename)
    CALL parser_create(reftraj%info%traj_parser,filename,para_env=para_env)

    CALL section_vals_val_get(reftraj_section,"VARIABLE_VOLUME",l_val=reftraj%info%variable_volume)
    IF(reftraj%info%variable_volume) THEN
       ! In case requested initialize parser for cell
       CALL section_vals_val_get(reftraj_section,"CELL_FILE_NAME",c_val=filename)
       CALL parser_create(reftraj%info%cell_parser,filename,para_env=para_env)
    END IF

    CALL section_vals_val_get(reftraj_section,"FIRST_SNAPSHOT",i_val=reftraj%info%first_snapshot)
    CALL section_vals_val_get(reftraj_section,"LAST_SNAPSHOT",i_val=reftraj%info%last_snapshot)
    CALL section_vals_val_get(reftraj_section,"STRIDE",i_val=reftraj%info%stride)
    CALL section_vals_val_get(reftraj_section,"EVAL_ENERGY_FORCES",l_val=reftraj%info%eval_ef)

    CALL section_vals_val_get(reftraj_section,"MSD%_SECTION_PARAMETERS_",&
            l_val=reftraj%info%msd)

  END SUBROUTINE create_reftraj

! *****************************************************************************
!> \brief ...
!> \param reftraj ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE retain_reftraj(reftraj)

    TYPE(reftraj_type), POINTER              :: reftraj

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_reftraj', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(reftraj)) THEN
       IF(.NOT.(reftraj%ref_count>0))CALL cp__a("motion/reftraj_types.F",141)
       reftraj%ref_count=reftraj%ref_count+1
    END IF

  END SUBROUTINE retain_reftraj

! *****************************************************************************
!> \brief ...
!> \param reftraj ...
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
  SUBROUTINE release_reftraj(reftraj)

    TYPE(reftraj_type), POINTER              :: reftraj

    CHARACTER(len=*), PARAMETER :: routineN = 'release_reftraj', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(reftraj)) THEN
       IF(.NOT.(reftraj%ref_count>0))CALL cp__a("motion/reftraj_types.F",162)
       reftraj%ref_count=reftraj%ref_count-1
       IF(reftraj%ref_count<1) THEN
          CALL parser_release(reftraj%info%traj_parser)
          CALL parser_release(reftraj%info%cell_parser)
          IF(ASSOCIATED(reftraj%info)) THEN
             DEALLOCATE (reftraj%info)
          END IF
          IF(ASSOCIATED(reftraj%msd)) THEN
             DEALLOCATE(reftraj%msd%ref0_pos)
             IF(reftraj%msd%msd_kind) THEN
                 DEALLOCATE(reftraj%msd%val_msd_kind)
             END IF
             IF(reftraj%msd%msd_molecule) THEN
                 DEALLOCATE(reftraj%msd%val_msd_molecule)
                 DEALLOCATE(reftraj%msd%ref0_com_molecule)
             END IF
             IF(reftraj%msd%disp_atom) THEN
                 DEALLOCATE(reftraj%msd%disp_atom_index)
                 DEALLOCATE(reftraj%msd%disp_atom_dr)
             END IF

             DEALLOCATE (reftraj%msd)
          END IF

          DEALLOCATE( reftraj )
       END IF

    END IF
  END SUBROUTINE release_reftraj

END MODULE reftraj_types
