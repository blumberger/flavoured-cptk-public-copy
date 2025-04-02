# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_distribution.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_distribution.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the distribution section of the input
!> \par History
!>      04.2007 created
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_distribution
  
  USE input_constants,                 ONLY: model_block_count,&
                                             model_block_lmax,&
                                             model_block_surface
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
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
# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/input_cp2k_distribution.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_distribution'

  PUBLIC :: create_distribution_section

CONTAINS

! *****************************************************************************
!> \brief Creates the distribution section
!> \param section the section to create
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_distribution_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_distribution_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    IF(.NOT.(.NOT.ASSOCIATED(section)))CALL cp__a("input_cp2k_distribution.F",49)
    CALL section_create(section,name="DISTRIBUTION",&
         description="can be used used to tune the parallel distribution of the data",&
         n_keywords=2, n_subsections=2, repeats=.FALSE.)

    NULLIFY(keyword)

    CALL keyword_create(keyword,name="COST_MODEL",&
         description="The cost model that needs to be minimized ",&
         usage="COST_MODEL BLOCK_COUNT",&
         enum_c_vals=s2a("BLOCK_COUNT","BLOCK_SURFACE", "BLOCK_LMAX"),&
         enum_i_vals=(/model_block_count, model_block_surface,model_block_lmax/),&
         enum_desc=s2a("the number of blocks",&
                       "the number of blocks weighted by the number elements per block",&
                       "the number of blocks weighted by the sum of the lmax"), &
         default_i_val=model_block_count)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="2D_MOLECULAR_DISTRIBUTION",&
         description="Distribute the atoms so that atoms belonging to a given molecule"//&
                     " are on the same CPU for the 2D distribution. This might give rise to a"//&
                     " worse distribution but reduces memory needs of finding the optimal distribution.",&
         usage="2D_MOLECULAR_DISTRIBUTION TRUE",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SKIP_OPTIMIZATION",&
         description="Do not optimize the distribution, go for something very simple."//&
                     " Might be useful if the optimization, which scales quadratically in system size, is too expensive.",&
         usage="SKIP_OPTIMIZATION TRUE",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="BASIC_OPTIMIZATION",&
         description="Creates a distribution based on a few heuristics using only minimal memory "//&
                     "and CPU time.",&
         usage="BASIC_OPTIMIZATION TRUE",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="BASIC_SPATIAL_OPTIMIZATION",&
         description="Creates a distribution with spatial info, using only minimal memory "//&
                     "and CPU time.",&
         usage="BASIC_SPATIAL_OPTIMIZATION TRUE",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="BASIC_CLUSTER_OPTIMIZATION",&
         description="Creates a distribution with spatial info, using recursively KMEANS clustering. ",&
         usage="BASIC_CLUSTER_OPTIMIZATION TRUE",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SYMMETRIC",&
         description="Take the symmetry of the distribution_2d into account.",&
         usage="SYMMETRIC TRUE",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

  END SUBROUTINE create_distribution_section

END MODULE input_cp2k_distribution
