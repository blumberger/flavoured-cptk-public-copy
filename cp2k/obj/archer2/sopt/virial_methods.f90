# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/virial_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/virial_methods.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH [04042007] code refactoring
! *****************************************************************************
MODULE virial_methods

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE virial_types,                    ONLY: virial_type

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
# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/virial_methods.F" 2

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: virial_evaluate, virial_pair_force, virial_update

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_methods'

CONTAINS
! *****************************************************************************
!> \brief Updates the virial given the virial and subsys
!> \param virial ...
!> \param subsys ...
!> \param para_env ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE virial_update(virial, subsys, para_env)
    TYPE(virial_type), INTENT(INOUT)         :: virial
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_update', &
      routineP = moduleN//':'//routineN

    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

    CALL cp_subsys_get(subsys, local_particles=local_particles, atomic_kinds=atomic_kinds,&
         particles=particles)

    CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles,&
         virial, para_env%group)

  END SUBROUTINE virial_update

! *****************************************************************************
!> \brief Computes the kinetic part of the pressure tensor and updates
!>      the full VIRIAL (PV)
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param virial ...
!> \param igroup ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE virial_evaluate ( atomic_kind_set, particle_set, local_particles,&
       virial, igroup)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: igroup

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      nnodes, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: mass, mfl
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    IF ( virial%pv_availability ) THEN
       CALL timeset(routineN,handle)
       NULLIFY(atomic_kind)
       mfl = 0.0_dp
       nparticle_kind = SIZE ( atomic_kind_set )
       virial%pv_kinetic = 0.0_dp
       DO i = 1, 3
          DO j = 1, i
             nnodes = 0
             DO iparticle_kind=1,nparticle_kind
                atomic_kind => atomic_kind_set(iparticle_kind)
                CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass )
                nparticle_local = local_particles%n_el(iparticle_kind)
                DO iparticle_local=1,nparticle_local
                   nnodes = nnodes + 1
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   virial%pv_kinetic(i,j) = virial%pv_kinetic(i,j) + &
                        mass * particle_set(iparticle)%v(i)*particle_set(iparticle)%v(j)
                END DO
             END DO
             virial%pv_kinetic(j,i) = virial%pv_kinetic(i,j)
          END DO
       END DO
       mfl = REAL( 9 * nnodes, KIND=dp) * 2.0_dp * 1.e-6_dp

       CALL mp_sum(virial%pv_kinetic,igroup)

       ! total virial
       virial%pv_total = virial%pv_virial + virial%pv_kinetic + virial%pv_constraint

       CALL timestop(handle)
    ENDIF

  END SUBROUTINE virial_evaluate

! *****************************************************************************
!> \brief Computes the contribution to the stress tensor from two-body
!>      pair-wise forces
!> \param pv_virial ...
!> \param f0 ...
!> \param force ...
!> \param rab ...
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE virial_pair_force ( pv_virial, f0, force, rab)

    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_virial
    REAL(KIND=dp)                            :: f0
    REAL(KIND=dp), DIMENSION(3)              :: force, rab

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_pair_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    DO i=1,3
       DO j=1,3
          pv_virial(i,j) = pv_virial(i,j) + f0 * force(i) * rab(j)
       END DO
    END DO

  END SUBROUTINE virial_pair_force

END MODULE virial_methods

