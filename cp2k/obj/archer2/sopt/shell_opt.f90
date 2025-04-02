# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/shell_opt.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/shell_opt.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
!> \author
! *****************************************************************************
MODULE shell_opt

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cg_optimizer,                    ONLY: geoopt_cg
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_create,&
                                             gopt_f_release,&
                                             gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE input_constants,                 ONLY: default_shellcore_method_id
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE integrator_utils,                ONLY: tmp_variables_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/../base/base_uses.f90" 1
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
# 43 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/shell_opt.F" 2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: optimize_shell_core
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'shell_opt'

CONTAINS

! *****************************************************************************
!> \brief Optimize shell-core positions along an MD run
!> \param force_env ...
!> \param particle_set ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param globenv ...
!> \param tmp ...
!> \param check ...
!> \author
! *****************************************************************************

  SUBROUTINE optimize_shell_core(force_env,particle_set,shell_particle_set,core_particle_set,globenv,tmp,check)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                shell_particle_set, &
                                                core_particle_set
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(tmp_variables_type), OPTIONAL, &
      POINTER                                :: tmp
    LOGICAL, INTENT(IN), OPTIONAL            :: check

    CHARACTER(len=*), PARAMETER :: routineN = 'optimize_shell_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iat, nshell
    LOGICAL                                  :: do_update, explicit, &
                                                my_check, optimize
    REAL(dp), DIMENSION(:), POINTER          :: dvec_sc, dvec_sc_0
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                geo_section, root_section

    NULLIFY(logger)
    logger => cp_get_default_logger()

    IF(.NOT.(ASSOCIATED(force_env)))CALL cp__a("motion/shell_opt.F",94)
    IF(.NOT.(ASSOCIATED(globenv)))CALL cp__a("motion/shell_opt.F",95)

    NULLIFY (gopt_param,force_env_section,gopt_env,dvec_sc,dvec_sc_0,root_section,geo_section)
    root_section      => force_env%root_section
    force_env_section => force_env%force_env_section
    geo_section       => section_vals_get_subs_vals(root_section,"MOTION%SHELL_OPT")

    CALL section_vals_get(geo_section, explicit=explicit)
    IF(.NOT. explicit) RETURN

    CALL timeset(routineN,handle)

    optimize = .FALSE.
    my_check = .FALSE.
    IF(PRESENT(check)) my_check = check
    IF(my_check) THEN
        NULLIFY(subsys, para_env, atomic_kinds, local_particles)
        CALL force_env_get(force_env=force_env, subsys=subsys, para_env=para_env)
        CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds, local_particles=local_particles)
        CALL check_shell_core_distance(atomic_kinds,local_particles,particle_set,shell_particle_set,&
             core_particle_set,para_env,optimize)

        IF(.NOT. optimize) THEN
           CALL timestop(handle)
           RETURN
        END IF
    END IF

    nshell = SIZE(shell_particle_set)
    ALLOCATE(dvec_sc(3*nshell))
    ALLOCATE(dvec_sc_0(3*nshell))
    DO i= 1,nshell
      dvec_sc(1+3*(i-1)) = core_particle_set(i)%r(1)-shell_particle_set(i)%r(1)
      dvec_sc(2+3*(i-1)) = core_particle_set(i)%r(2)-shell_particle_set(i)%r(2)
      dvec_sc(3+3*(i-1)) = core_particle_set(i)%r(3)-shell_particle_set(i)%r(3)
    END DO
    dvec_sc_0 = dvec_sc


    CALL gopt_param_read(gopt_param, geo_section, type_id=default_shellcore_method_id)
    CALL gopt_f_create(gopt_env, gopt_param, force_env=force_env, globenv=globenv,&
         geo_opt_section=geo_section)

    CALL cp_add_iter_level(logger%iter_info,"SHELL_OPT")
    gopt_env%eval_opt_geo = .FALSE.
    CALL geoopt_cg  (force_env,gopt_param,globenv,&
    geo_section, gopt_env, dvec_sc, do_update=do_update)
    IF(.NOT.do_update) THEN
      DO i= 1,nshell
        shell_particle_set(i)%r(1) = -dvec_sc_0(1+3*(i-1)) + core_particle_set(i)%r(1)
        shell_particle_set(i)%r(2) = -dvec_sc_0(2+3*(i-1)) + core_particle_set(i)%r(2)
        shell_particle_set(i)%r(3) = -dvec_sc_0(3+3*(i-1)) + core_particle_set(i)%r(3)
      END DO
    END IF
    CALL cp_rm_iter_level(logger%iter_info,"SHELL_OPT")

    CALL gopt_f_release(gopt_env)
    CALL gopt_param_release(gopt_param)
    DEALLOCATE(dvec_sc)
    DEALLOCATE(dvec_sc_0)

    IF(PRESENT(tmp)) THEN
      DO i=1,nshell
        iat = shell_particle_set(i)%atom_index
       tmp% shell_vel(1:3,i) = tmp%vel(1:3,iat)
       tmp% core_vel(1:3,i) = tmp%vel(1:3,iat)
      END DO
    ELSE
      DO i=1,nshell
        iat = shell_particle_set(i)%atom_index
        shell_particle_set(i)%v(1:3) = particle_set(iat)%v(1:3)
        core_particle_set(i)%v(1:3) = particle_set(iat)%v(1:3)
       END DO
    END IF


    CALL timestop(handle)

  END SUBROUTINE  optimize_shell_core

! *****************************************************************************
!> \brief Check shell_core_distance
!> \param atomic_kinds ...
!> \param local_particles ...
!> \param particle_set ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param para_env ...
!> \param optimize ...
!> \par History
!>      none
!> \author MI (October 2008)
!>     I soliti ignoti
! *****************************************************************************
 SUBROUTINE check_shell_core_distance(atomic_kinds,local_particles,particle_set,&
            shell_particle_set,core_particle_set,para_env,optimize)

    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                shell_particle_set, &
                                                core_particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(INOUT)                   :: optimize

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_shell_core_distance', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iparticle, &
                                                iparticle_local, itest, &
                                                nkind, nparticle_local, &
                                                shell_index
    LOGICAL                                  :: is_shell
    REAL(dp)                                 :: dsc, rc(3), rs(3)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    nkind = atomic_kinds%n_els
    itest = 0
    DO ikind = 1,nkind
      NULLIFY(atomic_kind)
      atomic_kind => atomic_kinds%els(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind, shell_active=is_shell, shell=shell)
      IF(is_shell) THEN
        IF(shell%max_dist > 0.0_dp) THEN
          nparticle_local = local_particles%n_el(ikind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(ikind)%array(iparticle_local)
             shell_index = particle_set(iparticle)%shell_index

             rc(:) = core_particle_set(shell_index)%r(:)
             rs(:) = shell_particle_set(shell_index)%r(:)
             dsc = SQRT((rc(1)-rs(1))**2 + (rc(2)-rs(2))**2 + (rc(3)-rs(3))**2)
             IF ( dsc > shell%max_dist) THEN
                itest = 1
             END IF
          END DO
        END IF
      END IF
    END DO

    CALL mp_sum(itest,para_env%group)
    IF(itest > 0) optimize = .TRUE.

 END SUBROUTINE check_shell_core_distance
END MODULE shell_opt
