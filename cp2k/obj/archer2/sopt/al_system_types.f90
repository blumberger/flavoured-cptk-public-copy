# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/al_system_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/al_system_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type for the canonical sampling through velocity rescaling     
!> \author Teodoro Laino - 09.2007 University of Zurich [tlaino]
! *****************************************************************************
MODULE al_system_types
  USE bibliography,                    ONLY: Jones2011,&
                                             cite_reference
  USE extended_system_types,           ONLY: create_map_info_type,&
                                             map_info_type,&
                                             release_map_info_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE simpar_types,                    ONLY: simpar_type

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
# 21 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/al_system_types.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: al_system_type,&
            al_init,&
            al_dealloc,&
            al_thermo_create

! *****************************************************************************
  TYPE al_thermo_type
     INTEGER                                 :: degrees_of_freedom
     REAL(KIND=dp)                           :: nkt
     REAL(KIND=dp)                           :: chi
     REAL(KIND=dp)                           :: mass
     REAL(KIND=dp)                           :: region_kin_energy
  END TYPE al_thermo_type

! *****************************************************************************
  TYPE al_system_type
     INTEGER                                 :: region, glob_num_al, loc_num_al
     REAL(KIND=dp)                           :: tau_nh, tau_langevin, dt_fact
     REAL(KIND=dp)                           :: dt
     TYPE(al_thermo_type), POINTER           :: nvt(:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE al_system_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'al_system_types'

CONTAINS

! *****************************************************************************
!> \brief Initialize type for Adaptive Langevin (AD_LANGEVIN)
!> \param al ...
!> \param simpar ...
!> \param section ...
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_init(al, simpar, section)
    TYPE(al_system_type), POINTER            :: al
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(section_vals_type), POINTER         :: section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(al%nvt)
    NULLIFY(al%map_info)
    al%loc_num_al=0
    al%glob_num_al=0
    al%dt_fact=1.0_dp
    al%dt=simpar%dt
    CALL cite_reference(Jones2011)
    CALL section_vals_val_get(section,"TIMECON_NH",r_val=al%tau_nh)
    CALL section_vals_val_get(section,"TIMECON_LANGEVIN",r_val=al%tau_langevin)
    CALL create_map_info_type(al%map_info)
    
  END SUBROUTINE al_init

! *****************************************************************************
!> \brief Initialize NVT type for AD_LANGEVIN thermostat
!> \param al ...
!> \author Noam Bernstein [noamb]  02.2012
! *****************************************************************************
  SUBROUTINE al_thermo_create(al)
    TYPE(al_system_type), POINTER            :: al

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_thermo_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: seed

    IF(.NOT.(ASSOCIATED(al)))CALL cp__a("al_system_types.F",96)
    IF(.NOT.(.NOT.ASSOCIATED(al%nvt)))CALL cp__a("al_system_types.F",97)
    
    ALLOCATE ( al%nvt(al%loc_num_al))
    DO i = 1, al%loc_num_al
       al%nvt(i)%chi = 0.0_dp
    END DO
    ! Initialize the gaussian stream random number
    ALLOCATE (seed(3,2,al%glob_num_al))

  END SUBROUTINE al_thermo_create

! *****************************************************************************
!> \brief Deallocate type for AD_LANGEVIN thermostat
!> \param al ...
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_dealloc ( al)
    TYPE(al_system_type), POINTER            :: al

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_dealloc', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(al)) THEN
       CALL al_thermo_dealloc(al%nvt)
       CALL release_map_info_type(al%map_info)
       DEALLOCATE (al)
    ENDIF

  END SUBROUTINE al_dealloc

! *****************************************************************************
!> \brief Deallocate NVT type for AD_LANGEVIN thermostat
!> \param nvt ...
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_thermo_dealloc ( nvt)
    TYPE(al_thermo_type), DIMENSION(:), &
      POINTER                                :: nvt

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_thermo_dealloc', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(nvt)) THEN
       DEALLOCATE (nvt)
    ENDIF
  END SUBROUTINE al_thermo_dealloc
  
END MODULE al_system_types

