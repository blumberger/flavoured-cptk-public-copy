# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/ps_wavelet_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/ps_wavelet_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the ps_wavelet data type.
!> \author Florian Schiffmann (09.2007,fschiff)
! *****************************************************************************
MODULE ps_wavelet_types

  
  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/../base/base_uses.f90" 1
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
# 15 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/ps_wavelet_types.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_types'

  INTEGER, PARAMETER, PUBLIC               :: WAVELET3D = 1400,&
                                              WAVELET2D = 1401,&
                                              WAVELET1D = 1402,&
                                              WAVELET0D = 1403

  PUBLIC :: ps_wavelet_type,&
            ps_wavelet_release

! *****************************************************************************
!> \par History
!>      09.2007 created [Florian Schiffmann]
!> \author fschiff
! *****************************************************************************
  TYPE ps_wavelet_type
     CHARACTER(LEN=1)                                  :: geocode
     CHARACTER(LEN=1)                                  :: datacode
     INTEGER                                           :: itype_scf
     INTEGER                                           :: method, special_dimension
     REAL(kind= dp), POINTER, DIMENSION(:)             :: karray
     REAL (KIND=dp), DIMENSION ( :, :, : ), POINTER    :: rho_z_sliced
     INTEGER,DIMENSION(3)                              :: PS_grid
  END TYPE ps_wavelet_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param wavelet ...
! *****************************************************************************
  SUBROUTINE ps_wavelet_release(wavelet)

    TYPE(ps_wavelet_type), POINTER           :: wavelet

    CHARACTER(len=*), PARAMETER :: routineN = 'ps_wavelet_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(wavelet)) THEN
       IF (ASSOCIATED(wavelet%karray))&
            DEALLOCATE(wavelet%karray)
       IF(ASSOCIATED(wavelet%rho_z_sliced))&
            DEALLOCATE(wavelet%rho_z_sliced)
       DEALLOCATE(wavelet)
    END IF
  END SUBROUTINE ps_wavelet_release

END MODULE ps_wavelet_types
