# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/dbcsr_performance_driver.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/dbcsr_performance_driver.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Performance tester for DBCSR operations
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
PROGRAM dbcsr_performance_driver
  USE dbcsr_config,                    ONLY: dbcsr_set_default_config
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_lib,                       ONLY: dbcsr_finalize_lib,&
                                             dbcsr_init_lib
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_new,&
                                             dbcsr_mp_release
  USE dbcsr_performance_multiply,      ONLY: dbcsr_perf_multiply
  USE dbcsr_test_methods,              ONLY: dbcsr_test_read_args
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj
  USE kinds,                           ONLY: default_string_length
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/../base/base_uses.f90" 1
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
# 37 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/dbcsr_performance_driver.F" 2

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

  IMPLICIT NONE


  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, io_unit, narg, handle
  INTEGER, DIMENSION(2)                    :: npdims, myploc
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(dbcsr_mp_obj)                       :: mp_env
  CHARACTER(len=default_string_length)     :: args(100)


  CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_multiply'


  !***************************************************************************************

  !
  ! initialize libdbcsr errors
  CALL timeset (routineN, handle)

  !
  ! initialize mpi
  CALL mp_world_init (mp_comm)

  !
  ! setup the mp environment
  npdims(:) = 0
  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)
  CALL mp_environ (numnodes, mynode, group)
  ALLOCATE(pgrid(0:npdims(1)-1,0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        CALL mp_cart_rank (group, (/ prow, pcol /), pgrid(prow, pcol))
     ENDDO
  ENDDO
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))
  DEALLOCATE(pgrid)

  !
  ! set standard output parameters
  io_unit = 0
  IF (mynode.EQ.0) io_unit = default_output_unit

  !
  ! read and distribute input args
  IF (mynode.eq.0) CALL dbcsr_test_read_args (narg, args)
  CALL mp_bcast (narg, 0, group)
  CALL mp_bcast (args, 0, group)
  CALL dbcsr_assert (narg.GE.1, dbcsr_fatal_level, dbcsr_wrong_args_error, &
     routineN, "nargs not correct", 90)

  !
  ! initialize libdbcsr
  CALL dbcsr_init_lib ()
  CALL dbcsr_set_default_config ()

  !
  ! select the operation
  SELECT CASE(args(1))
    CASE('dbcsr_multiply')
       CALL dbcsr_perf_multiply (group, mp_env, npdims, io_unit, narg, args)
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
          routineN, "operation not found", 104)
  END SELECT

  !
  ! finalize libdbcsr
  CALL dbcsr_finalize_lib (mp_comm, io_unit)

  !
  ! clean mp enviroment
  CALL dbcsr_mp_release (mp_env)

  !
  ! finalize mpi
  CALL mp_comm_free(group)
  CALL mp_world_finalize ()

  !
  ! finalize libdbcsr errors
  CALL timestop (handle)

END PROGRAM dbcsr_performance_driver
