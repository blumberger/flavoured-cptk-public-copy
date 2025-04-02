# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_cube.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_cube.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A wrapper around pw_to_cube() which accepts particle_list_type
!> \author Ole Schuett
! *****************************************************************************
MODULE cp_realspace_grid_cube
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_types,                        ONLY: pw_type
  USE realspace_grid_cube,             ONLY: pw_to_cube

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
# 17 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/cp_realspace_grid_cube.F" 2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_pw_to_cube

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_realspace_grid_cube'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles ...
!> \param stride ...
!> \param zero_tails ...
! *****************************************************************************
  SUBROUTINE cp_pw_to_cube ( pw, unit_nr, title, particles, stride, zero_tails)
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: unit_nr
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: stride
    LOGICAL, INTENT(IN), OPTIONAL            :: zero_tails

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_pw_to_cube', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: particles_z
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: particles_r
    TYPE(particle_list_type), POINTER        :: my_particles

    NULLIFY(my_particles)
    my_particles=>particles
    IF (ASSOCIATED(my_particles)) THEN
       n = my_particles%n_els
       ALLOCATE(particles_z(n))
       ALLOCATE(particles_r(3,n))
       DO i=1, n
          CALL get_atomic_kind(my_particles%els(i)%atomic_kind,z=particles_z(i))
          particles_r(:,i) = my_particles%els(i)%r(:)
       END DO

       CALL pw_to_cube(pw=pw, unit_nr=unit_nr, title=title, &
          particles_z=particles_z, particles_r=particles_r,&
          stride=stride, zero_tails=zero_tails)
    ELSE
       CALL pw_to_cube(pw=pw, unit_nr=unit_nr, title=title, &
         stride=stride, zero_tails=zero_tails)
    END IF

  END SUBROUTINE cp_pw_to_cube

END MODULE cp_realspace_grid_cube
