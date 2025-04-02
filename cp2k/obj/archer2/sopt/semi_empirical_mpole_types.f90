# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_mpole_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_mpole_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the semi empirical multipole integral expansions types
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
MODULE semi_empirical_mpole_types
  
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
# 14 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/semi_empirical_mpole_types.F" 2

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_mpole_types'

! *****************************************************************************
!> \brief Semi-empirical integral multipole expansion type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE semi_empirical_mpole_type
     LOGICAL, DIMENSION(3)                    :: task
     INTEGER                                  :: indi, indj
     REAL(KIND=dp)                            :: c
     REAL(KIND=dp), DIMENSION(3)              :: d
     REAL(KIND=dp), DIMENSION(3,3)            :: qc ! quadrupole cartesian
     REAL(KIND=dp), DIMENSION(5)              :: qs ! quadrupole spherical
     ! alternative definition used in GKS integral routines
     REAL(KIND=dp)                            :: cs
     REAL(KIND=dp), DIMENSION(3)              :: ds
     REAL(KIND=dp), DIMENSION(3,3)            :: qq ! quadrupole cartesian
  END TYPE semi_empirical_mpole_type

! *****************************************************************************
!> \brief Semi-empirical integral multipole expansion type - pointer type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE semi_empirical_mpole_p_type
     TYPE(semi_empirical_mpole_type), POINTER :: mpole
  END TYPE semi_empirical_mpole_p_type

! *****************************************************************************
!> \brief Global Multipolar NDDO information type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE nddo_mpole_type
     REAL(KIND=dp), DIMENSION(:), POINTER     :: charge,     efield0
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: dipole,     efield1, efield2
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: quadrupole
  END TYPE nddo_mpole_type


  PUBLIC :: semi_empirical_mpole_type,&
            semi_empirical_mpole_p_type,&
            semi_empirical_mpole_p_create,&
            semi_empirical_mpole_p_release,&
            nddo_mpole_type,&
            nddo_mpole_create,&
            nddo_mpole_release

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \param ndim ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_create(mpole, ndim)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpole
    INTEGER, INTENT(IN)                      :: ndim

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_mpole_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(.NOT.(.NOT.ASSOCIATED(mpole)))CALL cp__a("semi_empirical_mpole_types.F",86)
    ALLOCATE (mpole(ndim))
    DO i = 1, ndim
       NULLIFY(mpole(i)%mpole)
       CALL semi_empirical_mpole_create(mpole(i)%mpole)
    END DO

  END SUBROUTINE semi_empirical_mpole_p_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_release(mpole)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpole

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_mpole_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(mpole)) THEN
       DO i = 1, SIZE(mpole)
          CALL semi_empirical_mpole_release(mpole(i)%mpole)
       END DO
       DEALLOCATE (mpole)
    END IF

  END SUBROUTINE semi_empirical_mpole_p_release

! *****************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_create(mpole)
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(mpole)))CALL cp__a("semi_empirical_mpole_types.F",130)
    ALLOCATE (mpole)
    mpole%task = .FALSE.
    mpole%indi = 0
    mpole%indj = 0
    mpole%c    = HUGE(0.0_dp)
    mpole%d    = HUGE(0.0_dp)
    mpole%qc   = HUGE(0.0_dp)
    mpole%qs   = HUGE(0.0_dp)
    mpole%cs   = HUGE(0.0_dp)
    mpole%ds   = HUGE(0.0_dp)
    mpole%qq   = HUGE(0.0_dp)
  END SUBROUTINE semi_empirical_mpole_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_release(mpole)
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(mpole)) THEN
       DEALLOCATE (mpole)
    END IF

  END SUBROUTINE semi_empirical_mpole_release

! *****************************************************************************
!> \brief Allocate NDDO multipole type
!> \param nddo_mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE nddo_mpole_create(nddo_mpole)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_create', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(.NOT.ASSOCIATED(nddo_mpole)))CALL cp__a("semi_empirical_mpole_types.F",172)
    ALLOCATE (nddo_mpole)
    NULLIFY(nddo_mpole%charge)
    NULLIFY(nddo_mpole%dipole)
    NULLIFY(nddo_mpole%quadrupole)
    NULLIFY(nddo_mpole%efield0)
    NULLIFY(nddo_mpole%efield1)
    NULLIFY(nddo_mpole%efield2)
  END SUBROUTINE nddo_mpole_create

! *****************************************************************************
!> \brief Deallocate NDDO multipole type
!> \param nddo_mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE nddo_mpole_release(nddo_mpole)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(nddo_mpole)) THEN
       IF (ASSOCIATED(nddo_mpole%charge)) THEN
          DEALLOCATE(nddo_mpole%charge)
       END IF
       IF (ASSOCIATED(nddo_mpole%dipole)) THEN
          DEALLOCATE(nddo_mpole%dipole)
       END IF
       IF (ASSOCIATED(nddo_mpole%quadrupole)) THEN
          DEALLOCATE(nddo_mpole%quadrupole)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield0)) THEN
          DEALLOCATE(nddo_mpole%efield0)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield1)) THEN
          DEALLOCATE(nddo_mpole%efield1)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield2)) THEN
          DEALLOCATE(nddo_mpole%efield2)
       END IF
       DEALLOCATE (nddo_mpole)
    END IF

  END SUBROUTINE nddo_mpole_release

END MODULE semi_empirical_mpole_types
