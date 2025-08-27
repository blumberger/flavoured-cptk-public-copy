# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cell_opt_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cell_opt_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains type used for a Simulation Cell Optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
MODULE cell_opt_types
  USE cell_opt_utils,                  ONLY: get_ut_cell_matrix,&
                                             read_external_press_tensor
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type

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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/motion/cell_opt_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'
  INTEGER, PRIVATE, SAVE :: last_cell_opt_env_id=0

  PUBLIC :: cell_opt_env_type,&
            cell_opt_env_create,&
            cell_opt_env_release

! *****************************************************************************
!> \brief Type containing all informations abour the simulation cell optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  TYPE cell_opt_env_type
     ! Simulation cell optimization parameters
     INTEGER                                    :: ref_count, id_nr
     LOGICAL                                    :: keep_angles,&
                                                   keep_symmetry
     REAL(KIND=dp)                              :: pres_ext, pres_int, pres_tol
     REAL(KIND=dp), DIMENSION(3,3)              :: mtrx
     REAL(KIND=dp), DIMENSION(3,3)              :: rot_matrix
     TYPE(cell_type), POINTER                   :: ref_cell
  END TYPE cell_opt_env_type


CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param force_env ...
!> \param geo_section ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE cell_opt_env_create(cell_env, force_env, geo_section)
    TYPE(cell_opt_env_type), POINTER         :: cell_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: geo_section

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_opt_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp), DIMENSION(3)              :: r
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    IF(.NOT.(.NOT.ASSOCIATED(cell_env)))CALL cp__a("motion/cell_opt_types.F",83)
    ALLOCATE(cell_env)
    NULLIFY(cell_env%ref_cell,cell,subsys,particles)
    cell_env%ref_count=1
    last_cell_opt_env_id=last_cell_opt_env_id+1
    cell_env%id_nr=last_cell_opt_env_id
    CALL force_env_get(force_env, cell=cell, subsys=subsys)
    CALL cell_create(cell_env%ref_cell)
    CALL cell_clone(cell, cell_env%ref_cell)
    CALL section_vals_val_get(geo_section,"KEEP_ANGLES",l_val=cell_env%keep_angles)
    CALL section_vals_val_get(geo_section,"KEEP_SYMMETRY",l_val=cell_env%keep_symmetry)
    CALL section_vals_val_get(geo_section,"PRESSURE_TOLERANCE",r_val=cell_env%pres_tol)

    ! First let's rotate the cell vectors in order to have an upper triangular matrix.
    CALL get_ut_cell_matrix(cell)

    ! Compute the rotation matrix that give the cell vectors in the "canonical" orientation
    cell_env%rot_matrix = MATMUL(cell_env%ref_cell%hmat,cell%h_inv)

    ! Get the external pressure
    CALL read_external_press_tensor(geo_section,cell,cell_env%pres_ext,cell_env%mtrx,&
         cell_env%rot_matrix)

    ! Rotate particles accordingly
    CALL cp_subsys_get(subsys, particles=particles)
    DO ip=1,particles%n_els
       r = MATMUL(TRANSPOSE(cell_env%rot_matrix),particles%els(ip)%r)
       particles%els(ip)%r = r
    END DO
  END SUBROUTINE cell_opt_env_create

! *****************************************************************************
!> \brief ...
!> \param cell_env ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE cell_opt_env_release(cell_env)
    TYPE(cell_opt_env_type), POINTER         :: cell_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_opt_env_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(cell_env)) THEN
       IF(.NOT.(cell_env%ref_count>0))CALL cp__a("motion/cell_opt_types.F",128)
       cell_env%ref_count=cell_env%ref_count-1
       IF (cell_env%ref_count==0) THEN
          CALL cell_release(cell_env%ref_cell)
          DEALLOCATE(cell_env)
       END IF
    END IF
  END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
