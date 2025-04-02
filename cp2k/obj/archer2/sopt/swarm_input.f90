# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_input.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_input.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Declares the input for swarm framework
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_input
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             low_print_level
  USE glbopt_input,                    ONLY: glbopt_declare_input
  USE input_constants,                 ONLY: swarm_do_glbopt
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t
  USE string_utilities,                ONLY: s2a

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/../base/base_uses.f90" 1
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
# 27 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_input.F" 2

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_input'

 PUBLIC :: create_swarm_section

  CONTAINS


! *****************************************************************************
!> \brief Declares the SWARM input section
!> \param swarm_section ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE create_swarm_section(swarm_section)
    TYPE(section_type), POINTER              :: swarm_section

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_section, printkey

    NULLIFY(swarm_section, print_section, printkey, keyword)

    CALL section_create(swarm_section,name="SWARM",&
         description="Section to control swarm runs. "//&
         "The swarm framework provides a common ground for master/worker algorithms.",&
         repeats=.FALSE.)

    CALL keyword_create(keyword, name="BEHAVIOR",&
         description="Which behaviour should control the swarm.",&
         usage="BEHAVIOR <STRING>",&
         default_i_val=swarm_do_glbopt,&
         enum_c_vals=s2a("GLOBAL_OPT"),&
         enum_desc=s2a("Runs global geometry optimisation"),&
         enum_i_vals=(/swarm_do_glbopt/))
    CALL section_add_keyword(swarm_section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_OF_WORKERS",&
        description="Number of workers used for swarm. "//&
        "Of the total number of processors one is used for the master, "//&
        "the remaining processors should be divisible by the number of workers.",&
        type_of_var=integer_t)
    CALL section_add_keyword(swarm_section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="REPLAY_COMMUNICATION_LOG",&
           description="Filename of communication log of previous run. Use this to restart a swarm.",&
           repeats=.FALSE.,&
           usage="REPLAY_COMMUNICATION_LOG <CHARACTER>", default_lc_val="swarm_translog_replay.xyz")
    CALL section_add_keyword(swarm_section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_ITER",&
        description="The maximum number iterations the master should perform",&
        type_of_var=integer_t,default_i_val=HUGE(1))
    CALL section_add_keyword(swarm_section, keyword)
    CALL keyword_release(keyword)

    CALL section_create(print_section,name="PRINT",&
         description="Controls the printing properties during a global optimization run",&
         n_keywords=0, n_subsections=1, repeats=.TRUE.)

    CALL cp_print_key_section_create(printkey,"WORKER_RUN_INFO",&
               description="Controls the printing of the worker's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(print_section,printkey)
    CALL section_release(printkey)

    CALL cp_print_key_section_create(printkey,"MASTER_RUN_INFO",&
               description="Controls the printing of the masters's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(print_section,printkey)
    CALL section_release(printkey)

    CALL cp_print_key_section_create(printkey,"COMMUNICATION_LOG",&
            description="Log all the communication between workers and master. Needed for restart.",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom")
    CALL section_add_subsection(print_section,printkey)
    CALL section_release(printkey)

    CALL section_add_subsection(swarm_section,print_section)
    CALL section_release(print_section)


    CALL glbopt_declare_input(swarm_section)

  END SUBROUTINE create_swarm_section


END MODULE swarm_input

