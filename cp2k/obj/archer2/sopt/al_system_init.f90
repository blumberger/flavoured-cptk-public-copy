# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/al_system_init.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/al_system_init.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
MODULE al_system_init

  USE al_system_mapping,               ONLY: al_to_particle_mapping
  USE al_system_types,                 ONLY: al_system_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
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
# 26 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/al_system_init.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_al_part

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'al_system_init'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule ...
!> \param molecule_kind_set ...
!> \param para_env ...
!> \param al ...
!> \param al_section ...
!> \param gci ...
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE initialize_al_part ( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, al, al_section,&
       gci)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(al_system_type), POINTER            :: al
    TYPE(section_vals_type), POINTER         :: al_section
    TYPE(global_constraint_type), POINTER    :: gci

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_al_part', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: restart

    restart=.FALSE.
    CALL al_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
         molecule, molecule_kind_set, al, para_env, gci)

    CALL restart_al( al, al_section, restart)

    IF (.NOT. restart) THEN
      CALL init_al_variables(al)
    ENDIF

  END SUBROUTINE initialize_al_part

! *****************************************************************************
!> \brief ...
!> \param al ...
! *****************************************************************************
  SUBROUTINE init_al_variables(al)
    TYPE(al_system_type), POINTER            :: al

    CHARACTER(len=*), PARAMETER :: routineN = 'init_al_variables', &
      routineP = moduleN//':'//routineN

    al%nvt(:)%mass = al%nvt(:)%nkt * al%tau_nh**2
    al%nvt(:)%chi = 0.0_dp
  END SUBROUTINE init_al_variables

! *****************************************************************************
!> \brief ...
!> \param al ...
!> \param al_section ...
!> \param restart ...
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE restart_al(al, al_section, restart)
    TYPE(al_system_type), POINTER            :: al
    TYPE(section_vals_type), POINTER         :: al_section
    LOGICAL, INTENT(inout)                   :: restart

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_al', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_index, n_rep
    LOGICAL                                  :: explicit
    TYPE(section_vals_type), POINTER         :: work_section

    restart = .FALSE.

    ! Possibly restart the initial thermostat DOF value
    work_section => section_vals_get_subs_vals(section_vals=al_section,&
         subsection_name="CHI")
    CALL section_vals_get(work_section,explicit=explicit)
    restart = explicit
    IF (explicit) THEN
       CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
            n_rep_val=n_rep)
       IF (n_rep==al%glob_num_al) THEN
          DO i = 1, al%loc_num_al
             my_index = al%map_info%index(i)
             CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                  i_rep_val=my_index,r_val=al%nvt(i)%chi)
          END DO
       ELSE
          CALL cp_abort(cp__l("motion/thermostat/al_system_init.F",130),&
               'Number pf restartable stream not equal to the number of'//&
               ' total thermostats!')
       END IF
    END IF

    ! Possibly restart the initial thermostat mass
    work_section => section_vals_get_subs_vals(section_vals=al_section,&
         subsection_name="MASS")
    CALL section_vals_get(work_section,explicit=explicit)
    IF(restart.NEQV.explicit)&
       CALL cp_abort(cp__l("motion/thermostat/al_system_init.F",141),&
            "You need to define both CHI and MASS sections (or none) in the AD_LANGEVIN section")
    restart = restart.and.explicit
    IF (explicit) THEN
       CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
            n_rep_val=n_rep)
       IF (n_rep==al%glob_num_al) THEN
          DO i = 1, al%loc_num_al
             my_index = al%map_info%index(i)
             CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                  i_rep_val=my_index,r_val=al%nvt(i)%mass)
          END DO
       ELSE
          CALL cp_abort(cp__l("motion/thermostat/al_system_init.F",154),&
               'Number pf restartable stream not equal to the number of'//&
               ' total thermostats!')
       END IF
    END IF

  END SUBROUTINE restart_al

END MODULE al_system_init
