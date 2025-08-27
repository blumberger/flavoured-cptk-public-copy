# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/input_cp2k_barostats.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/input_cp2k_barostats.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_barostats
  USE barostat_types,                  ONLY: do_clv_x,&
                                             do_clv_xy,&
                                             do_clv_xyz,&
                                             do_clv_xz,&
                                             do_clv_y,&
                                             do_clv_yz,&
                                             do_clv_z
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_cp2k_thermostats,          ONLY: create_mass_section,&
                                             create_thermostat_section,&
                                             create_velocity_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a

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
# 37 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/thermostat/input_cp2k_barostats.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_barostats'

PUBLIC :: create_barostat_section

!***
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param section will contain the coeff section
!> \author teo
! *****************************************************************************
  SUBROUTINE create_barostat_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_barostat_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, thermo_section

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("motion/thermostat/input_cp2k_barostats.F",63)
    CALL section_create(section,name="barostat",&
         description="Parameters of barostat.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    NULLIFY(keyword,subsection,thermo_section)
    CALL keyword_create(keyword, name="PRESSURE",&
         description="Initial pressure",&
         usage="PRESSURE real",&
         default_r_val=0._dp,unit_str='bar')
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TIMECON",&
         description="Barostat time constant",&
         usage="TIMECON real",&
         default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs"),&
         unit_str='fs')
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMPERATURE",&
         description="Barostat initial temperature. If not set, the ensemble temperature is used instead.",&
         usage="TEMPERATURE real",type_of_var=real_t,&
         unit_str='K')
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMP_TOL",&
         description="Maximum oscillation of the Barostat temperature imposed by rescaling.",&
         usage="TEMP_TOL real",default_r_val=0._dp, &
         unit_str='K')
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VIRIAL",&
         description="For NPT_F only: allows the screening of one or more components of the virial in order"//&
         " to relax the cell only along specific cartesian axis",&
         usage="VIRIAL (XYZ | X | Y | Z | XY| XZ | YZ)",&
         enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
         enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
         default_i_val=do_clv_xyz)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL create_velocity_section(subsection,"BAROSTAT")
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL create_mass_section(subsection,"BAROSTAT")
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL create_thermostat_section(thermo_section, coupled_thermostat=.TRUE.)
    CALL section_add_subsection(section, thermo_section)
    CALL section_release(thermo_section)

    CALL create_print_section(subsection)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

  END SUBROUTINE create_barostat_section

! *****************************************************************************
!> \brief Creates print section for barostat section
!> \param section ...
!> \author teo [tlaino] - University of Zurich - 02.2008
! *****************************************************************************
  SUBROUTINE create_print_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_type), POINTER              :: print_key

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("motion/thermostat/input_cp2k_barostats.F",139)
    NULLIFY(print_key)
    CALL section_create(section,name="PRINT",&
         description="Collects all print_keys for barostat",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    CALL cp_print_key_section_create(print_key,"ENERGY",&
         description="Controls the output of kinetic energy, and potential energy "//&
         " of the defined barostat.", print_level=high_print_level, common_iter_levels=1,&
         filename="")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)
  END SUBROUTINE create_print_section

END MODULE input_cp2k_barostats
