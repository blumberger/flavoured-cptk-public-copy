# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_atprop.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_atprop.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE input_cp2k_atprop
  USE bibliography,                    ONLY: Kikuchi2009
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type

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
# 21 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_atprop.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_atprop'

  PUBLIC :: create_atprop_section

CONTAINS

! *****************************************************************************
!> \brief Creates the ATOMIC section
!> \param section the section to create
!> \author JHU
! *****************************************************************************
  SUBROUTINE create_atprop_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atprop_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("input_cp2k_atprop.F",45)
    CALL section_create(section,name="ATOMIC",&
         description="Controls the calculation of atomic properties. "//&
                     "Printing is controled by FORCE_EVAL / PRINT / PROGRAM_RUN_INFO",&
         repeats=.FALSE., &
         citations=(/Kikuchi2009/))

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="ENERGY",&
         description="Calculate atomic energies ",&
         usage="ENERGY {logical}",&
         repeats=.FALSE.,&
         n_var=1,&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PRESSURE",&
         description="Calculate atomic pressure tensors ",&
         usage="PRESSURE {logical}",&
         repeats=.FALSE.,&
         n_var=1,&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

  END SUBROUTINE create_atprop_section

END MODULE input_cp2k_atprop
