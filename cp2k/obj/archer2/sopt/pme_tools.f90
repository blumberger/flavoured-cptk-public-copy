# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pme_tools.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pme_tools.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Tools common both to PME and SPME
!> \par History
!>      JGH (03-May-2001) : first correctly working version
!>      teo (Feb-2007)    : Merging common routines to spme and pme
!> \author JGH (21-Mar-2001)
! *****************************************************************************
MODULE pme_tools

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             real_to_scaled
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE realspace_grid_types,            ONLY: realspace_grid_type

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
# 23 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pme_tools.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_center, set_list

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pme_tools'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param part ...
!> \param npart ...
!> \param center ...
!> \param p1 ...
!> \param rs ...
!> \param ipart ...
!> \param core_center ...
! *****************************************************************************
  SUBROUTINE set_list ( part, npart, center, p1, rs, ipart, core_center )

    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    INTEGER, INTENT(IN)                      :: npart
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: center
    INTEGER, INTENT(OUT)                     :: p1
    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, INTENT(INOUT)                   :: ipart
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: core_center

    INTEGER                                  :: ndim, npos
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(KIND=dp)                            :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    p1 = 0
    lb = rs % lb_real
    ub = rs % ub_real

    DO
       ipart = ipart + 1
       IF ( ipart > npart ) EXIT
       atomic_kind => part(ipart)%atomic_kind
       CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=charge)
       IF (charge == 0.0_dp.AND.part(ipart)%shell_index==0) CYCLE
       IF (rs % desc % parallel) THEN
          ! check if the rs grid is distributed or not
          IF ( ALL ( rs % desc % group_dim == 1 ) ) THEN
             ndim = rs % desc % group_size
             npos = rs % desc % my_pos
             ! All processors work on the same grid
             IF ( MOD ( ipart, ndim ) == npos ) THEN
                p1 = ipart
                EXIT
             END IF
          ELSE
             ! First check if this atom is on my grid
             IF(part(ipart)%shell_index/=0 .AND. PRESENT(core_center))THEN
                IF ( in_slice(core_center( : , part(ipart)%shell_index ), lb, ub)) THEN
                   p1 = ipart
                ENDIF
             ELSE
                IF ( in_slice ( center ( : , ipart ), lb, ub ) ) THEN
                   p1 = ipart
                   EXIT
                END IF
             END IF
          ENDIF
       ELSE
          p1 = ipart
          EXIT
       END IF
    END DO

  END SUBROUTINE set_list

! *****************************************************************************
!> \brief ...
!> \param pos ...
!> \param lb ...
!> \param ub ...
!> \retval internal ...
! *****************************************************************************
  FUNCTION in_slice ( pos, lb, ub ) RESULT ( internal )

    INTEGER, DIMENSION(3), INTENT(IN)        :: pos, lb, ub
    LOGICAL                                  :: internal

    IF ( ALL ( pos >= lb ) .AND. ALL ( pos <= ub ) ) THEN
       internal = .TRUE.
    ELSE
       internal = .FALSE.
    END IF

  END FUNCTION in_slice

! *****************************************************************************
!> \brief ...
!> \param part ...
!> \param box ...
!> \param center ...
!> \param npts ...
!> \param n ...
! *****************************************************************************
  SUBROUTINE get_center ( part, box, center, npts, n )

    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), POINTER                 :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: center
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: ipart, mp
    REAL(KIND=dp)                            :: rmp
    REAL(KIND=dp), DIMENSION(3)              :: gp, s

    mp = MAXVAL ( npts(:) )
    rmp = REAL ( mp,KIND=dp)
    DO ipart = 1, SIZE ( part )
       ! compute the scaled coordinate of atom ipart
       CALL real_to_scaled(s,part(ipart)%r,box)
       s = s - NINT ( s )
       gp = REAL ( npts,KIND=dp) * s
       ! find the closest grid point (on big grid)
       IF ( MOD ( n, 2 ) == 0 ) THEN
          center ( :, ipart ) = INT ( gp + rmp ) - mp
       ELSE
          center ( :, ipart ) = NINT ( gp )
       END IF
       center ( :, ipart ) = center ( :, ipart ) + npts(:)/2
       center ( :, ipart ) = MODULO ( center ( :, ipart ), npts(:) )
       center ( :, ipart ) = center ( :, ipart ) - npts(:)/2
    END DO

  END SUBROUTINE get_center

END MODULE pme_tools

