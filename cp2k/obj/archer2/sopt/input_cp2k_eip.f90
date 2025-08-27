# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_eip.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_eip.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Creates the EIP section of the input
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE input_cp2k_eip
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level,&
                                             medium_print_level
  USE input_constants,                 ONLY: use_bazant_eip,&
                                             use_lenosky_eip
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: enum_t
  USE string_utilities,                ONLY: s2a

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
# 29 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_eip.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_eip'

  PUBLIC :: create_eip_section

CONTAINS

! *****************************************************************************
!> \brief Create the input section for EIP
!> \param section the section to create
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE create_eip_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eip_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

!   ------------------------------------------------------------------------

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("input_cp2k_eip.F",58)
    CALL section_create(section,name="EIP", &
         description="This section contains all information to run an "//&
         "Empirical Interatomic Potential (EIP) calculation.", &
         n_keywords=1, n_subsections=1, repeats=.FALSE.)

    NULLIFY(subsection, keyword)

    CALL keyword_create(keyword, name="EIP_MODEL", &
         description="Selects the empirical interaction potential model", &
         usage="EIP_MODEL BAZANT",  type_of_var=enum_t, &
         n_var=1, repeats=.FALSE., variants=(/"EIP-MODEL"/), &
         enum_c_vals=s2a("BAZANT", "EDIP", "LENOSKY"), &
         enum_i_vals=(/use_bazant_eip, use_bazant_eip, use_lenosky_eip/), &
         enum_desc=s2a("Bazant potentials",&
                       "Environment-Dependent Interatomic Potential",&
                       "Lenosky potentials"),&
         default_i_val=use_lenosky_eip)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL create_eip_print_section(subsection)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)


  END SUBROUTINE create_eip_section

! *****************************************************************************
!> \brief Creates the print section for the eip subsection
!> \param section the section to create
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE create_eip_print_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eip_print_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_type), POINTER              :: print_key

!   ------------------------------------------------------------------------

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("input_cp2k_eip.F",103)
    CALL section_create(section, name="PRINT", &
         description="Section of possible print options in EIP code.", &
         n_keywords=0, n_subsections=6, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key, "ENERGIES", &
         description="Controls the printing of the EIP energies.", &
         print_level=medium_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key, "ENERGIES_VAR", &
         description="Controls the printing of the variance of the EIP energies.", &
         print_level=high_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key, "FORCES", &
         description="Controls the printing of the EIP forces.", &
         print_level=medium_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key, "COORD_AVG", &
         description="Controls the printing of the average coordination number.", &
         print_level=high_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key, "COORD_VAR", &
         description="Controls the printing of the variance of the coordination number.", &
         print_level=high_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key, "COUNT", &
         description="Controls the printing of the number of function calls.", &
         print_level=high_print_level, filename="__STD_OUT__")
    CALL section_add_subsection(section, print_key)
    CALL section_release(print_key)

  END SUBROUTINE create_eip_print_section

END MODULE input_cp2k_eip
