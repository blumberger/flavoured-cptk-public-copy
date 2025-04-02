# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_charges_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_charges_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief container for information about total charges on the grids
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_charges_types
  
  USE kinds,                           ONLY: dp

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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_charges_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_charges_types'

  PUBLIC :: qs_charges_type
  PUBLIC :: qs_charges_create, qs_charges_retain, qs_charges_release
!***

! *****************************************************************************
!> \brief Container for information about total charges on the grids
!> \param total_rho_core_rspace total charge on the rho_core grid
!> \param total_rho_rspace total charge in the real space
!> \param total_rho_gspace total charge in the g space
!> \note
!>      this type is loosing the reason to exist...
!> \par History
!>      10.2002 created [fawzi]
!>      11.2002 moved total_rho_elec_rspace to qs_rho_type
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_charges_type
     INTEGER :: ref_count
     REAL(KIND = dp) :: total_rho_core_rspace, total_rho_gspace
     REAL(KIND = dp) :: total_rho0_soft_rspace, total_rho0_hard_lebedev
     REAL(KIND = dp) :: total_rho_soft_gspace
     REAL(KIND = dp), DIMENSION(:), POINTER  :: total_rho1_hard,&
                                                total_rho1_soft
     REAL(KIND = dp) :: background
  END TYPE qs_charges_type

CONTAINS

! *****************************************************************************
!> \brief creates a charges object
!> \param qs_charges the charges object to create
!> \param nspins ...
!> \param total_rho_core_rspace ...
!> \param total_rho_gspace ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_create(qs_charges,nspins,total_rho_core_rspace, &
     total_rho_gspace)
    TYPE(qs_charges_type), POINTER           :: qs_charges
    INTEGER, INTENT(in)                      :: nspins
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: total_rho_core_rspace, &
                                                total_rho_gspace

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_create', &
      routineP = moduleN//':'//routineN

  ALLOCATE(qs_charges)
  qs_charges%total_rho_core_rspace=0.0_dp
  IF (PRESENT(total_rho_core_rspace)) &
       qs_charges%total_rho_core_rspace=total_rho_core_rspace
  qs_charges%total_rho_gspace=0.0_dp
  IF (PRESENT(total_rho_gspace)) &
       qs_charges%total_rho_gspace=total_rho_gspace
  qs_charges%total_rho_soft_gspace = 0.0_dp
  qs_charges%total_rho0_hard_lebedev = 0.0_dp
  qs_charges%total_rho_soft_gspace = 0.0_dp
  qs_charges%background = 0.0_dp
  ALLOCATE(qs_charges%total_rho1_hard(nspins))
  qs_charges%total_rho1_hard(:) = 0.0_dp
  ALLOCATE(qs_charges%total_rho1_soft(nspins))
  qs_charges%total_rho1_soft(:) = 0.0_dp
  qs_charges%ref_count=1
END SUBROUTINE qs_charges_create

! *****************************************************************************
!> \brief retains the given qs_charges (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to retain
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_retain(qs_charges)
    TYPE(qs_charges_type), POINTER           :: qs_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_retain', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(qs_charges)))CALL cp__a("qs_charges_types.F",102)
  IF(.NOT.(qs_charges%ref_count>0))CALL cp__a("qs_charges_types.F",103)
  qs_charges%ref_count=qs_charges%ref_count+1
END SUBROUTINE qs_charges_retain

! *****************************************************************************
!> \brief releases the charges object (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to be released
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_release(qs_charges)
    TYPE(qs_charges_type), POINTER           :: qs_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_release', &
      routineP = moduleN//':'//routineN

  IF (ASSOCIATED(qs_charges)) THEN
     IF(.NOT.(qs_charges%ref_count>0))CALL cp__a("qs_charges_types.F",121)
     qs_charges%ref_count=qs_charges%ref_count-1
     IF (qs_charges%ref_count<1) THEN
        DEALLOCATE(qs_charges%total_rho1_hard)
        DEALLOCATE(qs_charges%total_rho1_soft)
        DEALLOCATE(qs_charges)
     END IF
  END IF
  NULLIFY(qs_charges)
END SUBROUTINE qs_charges_release

END MODULE qs_charges_types
