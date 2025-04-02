# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/csvr_system_init.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/csvr_system_init.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
MODULE csvr_system_init

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE csvr_system_mapping,             ONLY: csvr_to_barostat_mapping,&
                                             csvr_to_particle_mapping,&
                                             csvr_to_shell_mapping
  USE csvr_system_types,               ONLY: csvr_system_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: read_rng_stream,&
                                             rng_record_length
  USE simpar_types,                    ONLY: simpar_type
  USE thermostat_types,                ONLY: thermostat_info_type

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/../../base/base_uses.f90" 1
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
# 29 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/csvr_system_init.F" 2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initialize_csvr_part, initialize_csvr_baro,&
            initialize_csvr_shell

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_init'

CONTAINS

! *****************************************************************************
!> \brief fire up the thermostats, if NPT
!> \param simpar ...
!> \param csvr ...
!> \param csvr_section ...
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE initialize_csvr_baro ( simpar, csvr, csvr_section)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: csvr_section

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_csvr_baro', &
      routineP = moduleN//':'//routineN

    CALL csvr_to_barostat_mapping ( simpar,  csvr)
    CALL restart_csvr( csvr, csvr_section)

  END SUBROUTINE initialize_csvr_baro

! *****************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule ...
!> \param molecule_kind_set ...
!> \param para_env ...
!> \param csvr ...
!> \param csvr_section ...
!> \param gci ...
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE initialize_csvr_part ( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, csvr, csvr_section,&
       gci)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: csvr_section
    TYPE(global_constraint_type), POINTER    :: gci

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_csvr_part', &
      routineP = moduleN//':'//routineN

    CALL csvr_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
         molecule, molecule_kind_set, csvr, para_env, gci)
    CALL restart_csvr( csvr, csvr_section)

  END SUBROUTINE initialize_csvr_part

! *****************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule ...
!> \param molecule_kind_set ...
!> \param para_env ...
!> \param csvr ...
!> \param csvr_section ...
!> \param gci ...
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE initialize_csvr_shell( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, csvr, csvr_section,&
       gci)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: csvr_section
    TYPE(global_constraint_type), POINTER    :: gci

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_csvr_shell', &
      routineP = moduleN//':'//routineN

    CALL csvr_to_shell_mapping(thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, csvr, para_env, gci)
    CALL restart_csvr( csvr, csvr_section)

  END SUBROUTINE  initialize_csvr_shell

! *****************************************************************************
!> \brief ...
!> \param csvr ...
!> \param csvr_section ...
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE restart_csvr(csvr, csvr_section)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: csvr_section

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_csvr', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: i, my_index, n_rep
    LOGICAL                                  :: explicit
    TYPE(section_vals_type), POINTER         :: work_section

! Possibly restart the initial thermostat energy

    work_section => section_vals_get_subs_vals(section_vals=csvr_section,&
         subsection_name="THERMOSTAT_ENERGY")
    CALL section_vals_get(work_section,explicit=explicit)
    IF (explicit) THEN
       CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
            n_rep_val=n_rep)
       IF (n_rep==csvr%glob_num_csvr) THEN
          DO i = 1, csvr%loc_num_csvr
             my_index = csvr%map_info%index(i)
             CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                  i_rep_val=my_index,r_val=csvr%nvt(i)%thermostat_energy)
          END DO
       ELSE
          CALL cp_abort(cp__l("motion/thermostat/csvr_system_init.F",167),&
               'Number pf restartable stream not equal to the number of'//&
               ' total thermostats!')
       END IF
    END IF

    ! Possibly restart the random number generators for the different thermostats
    work_section => section_vals_get_subs_vals(section_vals=csvr_section,&
         subsection_name="RNG_INIT")

    CALL section_vals_get(work_section,explicit=explicit)
    IF (explicit) THEN
       CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
            n_rep_val=n_rep)
       IF (n_rep==csvr%glob_num_csvr) THEN
          DO i = 1, csvr%loc_num_csvr
             my_index = csvr%map_info%index(i)
             CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                  i_rep_val=my_index,c_val=rng_record)
             CALL read_rng_stream(rng_stream=csvr%nvt(i)%gaussian_rng_stream,&
                  rng_record=rng_record)
          END DO
       ELSE
          CALL cp_abort(cp__l("motion/thermostat/csvr_system_init.F",190),&
               'Number pf restartable stream not equal to the number of'//&
               ' total thermostats!')
       END IF
    END IF
  END SUBROUTINE restart_csvr

END MODULE csvr_system_init
