# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_mom_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_mom_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief manage control variables for the maximum overlap method
! *****************************************************************************
MODULE qs_mom_types
  USE bibliography,                    ONLY: Gilbert2008
  USE input_constants,                 ONLY: momproj_norm,&
                                             momproj_sum
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t
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
# 22 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_mom_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mom_types'

  PUBLIC :: create_mom_section

CONTAINS

! *****************************************************************************
!> \brief Create CP2K input section for variable occupancy using the Maximum
!>        Overlap Method. Only with diagonalization methods, i.e. not with OT
!> \param section section to create
!> \date       20.06.2013
!> \par History
!>      06.2013 created [MattW]
!>      01.2016 (DE)OCC_ALPHA and (DE)OCC_BETA keywords accept a list of
!>              molecular orbitals. Added two extra keywords: START_ITER
!>              and PROJ_FORMULA [Sergey Chulkov]
!> \author     MattW
!> \version    1.0
! *****************************************************************************
  SUBROUTINE create_mom_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mom_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("qs_mom_types.F",54)

    CALL section_create(section,&
                        name="MOM",&
                        description="Define type and parameters for the maximum overlap method (MOM) "//&
                        "to determine orbital occupancies. "//&
                        "The MOM procedures activated by this section are only active for diagonalization "//&
                        "methods, i.e. not with minimization methods based on OT.",&
                        n_keywords=7, n_subsections=0, repeats=.FALSE., &
                        citations=(/Gilbert2008/))

    NULLIFY (keyword)

    CALL keyword_create(keyword,&
                        name="_SECTION_PARAMETERS_",&
                        description="Controls the activation of the MOM procedure",&
                        usage="MOM ON",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,&
                        name="START_ITER",&
                        description="SCF iteration cycle to start the MOM procedure",&
                        repeats=.FALSE.,&
                        n_var=1,&
                        type_of_var=integer_t,&
                        default_i_val=0,&
                        usage="START_ITER 2")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,&
                        name="DEOCC_ALPHA",&
                        description="Alpha orbitals to be deoccupied",&
                        repeats=.FALSE.,&
                        n_var=-1,&
                        type_of_var=integer_t,&
                        default_i_val=0,&
                        usage="DEOCC_ALPHA 10 11 ...")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,&
                        name="DEOCC_BETA",&
                        description="Beta orbitals to be deoccupied",&
                        repeats=.FALSE.,&
                        n_var=-1,&
                        type_of_var=integer_t,&
                        default_i_val=0,&
                        usage="DEOCC_BETA 10 11 ...")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,&
                        name="OCC_ALPHA",&
                        description="Alpha orbitals to be occupied",&
                        repeats=.FALSE.,&
                        n_var=-1,&
                        type_of_var=integer_t,&
                        default_i_val=0,&
                        usage="OCC_ALPHA 12 15 ...")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,&
                        name="OCC_BETA",&
                        description="Beta orbitals to be occupied",&
                        repeats=.FALSE.,&
                        n_var=-1,&
                        type_of_var=integer_t,&
                        default_i_val=0,&
                        usage="OCC_BETA 12 15 ...")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PROJ_FORMULA",&
         description="Projection formula to be used",&
         usage="SCF_GUESS RESTART", default_i_val=momproj_sum,&
         enum_c_vals=s2a("SUM","NORM"),&
         enum_desc=s2a("The one proposed in the original paper", &
         "The one which ignores the phase of molecular orbitals"), &
         enum_i_vals=(/momproj_sum, momproj_norm/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)


  END SUBROUTINE create_mom_section

END MODULE qs_mom_types
