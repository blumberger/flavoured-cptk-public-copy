# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_init.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_init.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!>      Routine to initalize a real space grid from a given input section
!> \par History
!>      01.2014 moved routine from realspace_grid_types into separate file.
!> \author Ole Schuett
! *****************************************************************************
MODULE cp_realspace_grid_init
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE realspace_grid_types,            ONLY: realspace_grid_input_type,&
                                             rsgrid_automatic,&
                                             rsgrid_replicated

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
# 22 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_init.F" 2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_input_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_realspace_grid_init'

CONTAINS

! *****************************************************************************
!> \brief parses an input section to assign the proper values to the input type
!> \param input_settings ...
!> \param nsmax ...
!> \param rs_grid_section ...
!> \param ilevel ...
!> \param higher_grid_layout the layout of a higher level grid. layouts with
!>       negative or zero values are ignored
!> \par History
!>      01.2008 created [Joost VandeVondele]
!> \note
!>      if rs_grid_section is not present we setup for an replicated setup
! *****************************************************************************
  SUBROUTINE init_input_type(input_settings,nsmax,rs_grid_section,ilevel,higher_grid_layout)
    TYPE(realspace_grid_input_type), &
      INTENT(OUT)                            :: input_settings
    INTEGER, INTENT(IN)                      :: nsmax
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: rs_grid_section
    INTEGER, INTENT(IN)                      :: ilevel
    INTEGER, DIMENSION(3), INTENT(IN)        :: higher_grid_layout

    INTEGER                                  :: isection, &
                                                max_distributed_level, &
                                                nsection
    INTEGER, DIMENSION(:), POINTER           :: tmp

    IF (PRESENT(rs_grid_section)) THEN
       input_settings%nsmax=nsmax
       ! we use the section corresponding to the level, or the largest available one
       ! i.e. the last section defines all following ones
       CALL section_vals_get(rs_grid_section,n_repetition=nsection)
       isection=MAX(1,MIN(ilevel,nsection))
       CALL section_vals_val_get(rs_grid_section,"DISTRIBUTION_TYPE",&
            i_rep_section=isection,&
            i_val=input_settings%distribution_type)
       CALL section_vals_val_get(rs_grid_section,"DISTRIBUTION_LAYOUT",&
            i_rep_section=isection,&
            i_vals=tmp)
       input_settings%distribution_layout=tmp
       CALL section_vals_val_get(rs_grid_section,"MEMORY_FACTOR",&
            i_rep_section=isection,&
            r_val=input_settings%memory_factor)
       CALL section_vals_val_get(rs_grid_section,"HALO_REDUCTION_FACTOR",&
            i_rep_section=isection,&
            r_val=input_settings%halo_reduction_factor)
       CALL section_vals_val_get(rs_grid_section,"LOCK_DISTRIBUTION",&
            i_rep_section=isection,&
            l_val=input_settings%lock_distribution)
       CALL section_vals_val_get(rs_grid_section,"MAX_DISTRIBUTED_LEVEL",&
            i_rep_section=isection,&
            i_val=max_distributed_level)

       ! multigrids that are to coarse are not distributed in the automatic scheme
       IF (input_settings%distribution_type == rsgrid_automatic) THEN
          IF (ilevel>max_distributed_level) THEN
             input_settings%distribution_type=rsgrid_replicated
          ENDIF
       ENDIF
    ELSE
       input_settings%nsmax=-1
       input_settings%distribution_type=rsgrid_replicated
       input_settings%lock_distribution=.FALSE.
       input_settings%halo_reduction_factor=1.0_dp
    ENDIF
    IF (input_settings%lock_distribution) THEN
       IF (ALL(higher_grid_layout>0)) input_settings%distribution_layout=higher_grid_layout
    ENDIF
  END SUBROUTINE init_input_type

END MODULE cp_realspace_grid_init
