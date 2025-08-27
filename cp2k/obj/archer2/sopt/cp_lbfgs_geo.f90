# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cp_lbfgs_geo.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cp_lbfgs_geo.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Main driver for L-BFGS optimizer
!> \par History
!>      none
! *****************************************************************************
MODULE cp_lbfgs_geo
  USE cell_types,                      ONLY: cell_type
  USE cp_external_control,             ONLY: external_control
  USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                             cp_opt_gopt_create,&
                                             cp_opt_gopt_next,&
                                             cp_opt_gopt_release,&
                                             cp_opt_gopt_stop
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: gopt_f_ii,&
                                             gopt_f_io_finalize,&
                                             print_geo_opt_header,&
                                             print_geo_opt_nc
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_ts_method_id
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp

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
# 40 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cp_lbfgs_geo.F" 2

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_lbfgs_geo'

  PUBLIC :: geoopt_lbfgs

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  RECURSIVE SUBROUTINE geoopt_lbfgs(force_env, gopt_param, globenv, geo_section, gopt_env,&
                          x0)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: geo_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0

    CHARACTER(len=*), PARAMETER :: routineN = 'geoopt_lbfgs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iter_nr, its, &
                                                output_unit
    LOGICAL                                  :: converged, should_stop
    REAL(KIND=dp)                            :: trust_radius
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section

    CALL timeset(routineN,handle)

    NULLIFY (optimizer, para_env)
    logger => cp_get_default_logger()
    root_section => force_env%root_section
    IF(.NOT.(ASSOCIATED(force_env)))CALL cp__a("motion/cp_lbfgs_geo.F",89)
    IF(.NOT.(ASSOCIATED(gopt_param)))CALL cp__a("motion/cp_lbfgs_geo.F",90)
    IF(.NOT.(gopt_param%ref_count>0))CALL cp__a("motion/cp_lbfgs_geo.F",91)

    CALL force_env_get(force_env, para_env=para_env, cell=cell)

    ! Geometry optimization starts now
    output_unit = cp_print_key_unit_nr(logger,geo_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".geoLog")
    CALL print_geo_opt_header(gopt_env, output_unit, "L-BFGS")

    ! Stop if not implemented
    IF(gopt_env%type_id == default_ts_method_id) &
       CALL cp__b("motion/cp_lbfgs_geo.F",102,"BFGS method not yet working with DIMER")

    CALL section_vals_val_get(geo_section,"LBFGS%TRUST_RADIUS",r_val=trust_radius)
    CALL cp_opt_gopt_create(optimizer, para_env=para_env, obj_funct=gopt_env,&
         x0=x0, wanted_relative_f_delta=gopt_param%wanted_rel_f_error,&
         wanted_projected_gradient=gopt_param%wanted_proj_gradient, m=gopt_param%max_h_rank,&
         max_f_per_iter=gopt_param%max_f_per_iter,trust_radius=trust_radius)
    CALL cp_iterate(logger%iter_info,increment=0,iter_nr_out=iter_nr)
    converged=.FALSE.

    DO its=iter_nr+1,gopt_param%max_iter
       CALL cp_iterate(logger%iter_info,last=(its==gopt_param%max_iter))
       CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=its)
       CALL gopt_f_ii(its, output_unit)

       ! Real optimization step..
       IF (.NOT.cp_opt_gopt_next(optimizer,geo_section=geo_section,&
            force_env=force_env,gopt_param=gopt_param, converged=converged)) EXIT

       ! Check for an external exit command
       CALL external_control(should_stop,"GEO",globenv=globenv)
       IF (should_stop) THEN
          CALL cp_opt_gopt_stop(optimizer)
          EXIT
       END IF
       IF (its == gopt_param%max_iter) EXIT
    END DO

    IF ((its == gopt_param%max_iter).AND.(.NOT.converged))THEN
       CALL print_geo_opt_nc(gopt_env, output_unit)
    END IF

    ! Write final output information, if converged
    CALL cp_iterate(logger%iter_info,last=.TRUE.,increment=0)
    CALL gopt_f_io_finalize(gopt_env, force_env, optimizer%x, converged, its, root_section,&
         optimizer%para_env, optimizer%master, output_unit)

    CALL cp_opt_gopt_release(optimizer)
    CALL cp_print_key_finished_output(output_unit,logger,geo_section,&
            "PRINT%PROGRAM_RUN_INFO")

    CALL timestop(handle)

  END SUBROUTINE geoopt_lbfgs

END MODULE cp_lbfgs_geo
