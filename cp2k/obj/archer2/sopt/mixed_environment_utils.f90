# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_environment_utils.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_environment_utils.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Util mixed_environment
!> \author Teodoro Laino [tlaino] - 02.2011
! *****************************************************************************
MODULE mixed_environment_utils

  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mixed_energy_types,              ONLY: mixed_force_type
  USE particle_list_types,             ONLY: particle_list_type
  USE virial_types,                    ONLY: virial_p_type,&
                                             virial_type,&
                                             zero_virial

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
# 30 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/mixed_environment_utils.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_environment_utils'

  PUBLIC :: mixed_map_forces,&
            get_subsys_map_index

CONTAINS

! *****************************************************************************
!> \brief Maps forces between the different force_eval sections/environments
!> \param particles_mix ...
!> \param virial_mix ...
!> \param results_mix ...
!> \param global_forces ...
!> \param virials ...
!> \param results ...
!> \param factor ...
!> \param iforce_eval ...
!> \param nforce_eval ...
!> \param map_index ...
!> \param mapping_section ...
!> \param overwrite ...
!> \author Teodoro Laino - University of Zurich [tlaino] - 05.2007
! *****************************************************************************
  SUBROUTINE mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces,&
      virials, results, factor, iforce_eval, nforce_eval, map_index,&
      mapping_section, overwrite)

    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(cp_result_type), POINTER            :: results_mix
    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
    TYPE(virial_p_type), DIMENSION(:), &
      POINTER                                :: virials
    TYPE(cp_result_p_type), DIMENSION(:), &
      POINTER                                :: results
    REAL(KIND=dp), INTENT(IN)                :: factor
    INTEGER, INTENT(IN)                      :: iforce_eval, nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: map_index
    TYPE(section_vals_type), POINTER         :: mapping_section
    LOGICAL, INTENT(IN)                      :: overwrite

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_map_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: iparticle, jparticle, natom, &
                                                nres
    LOGICAL                                  :: dip_exists
    REAL(KIND=dp), DIMENSION(3)              :: dip_mix, dip_tmp

! Get Mapping index array

    natom = SIZE(global_forces(iforce_eval)%forces,2)
    CALL get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, map_index)
    DO iparticle = 1, natom
       jparticle = map_index(iparticle)
       IF (overwrite) THEN
          particles_mix%els(jparticle)%f(:)= factor* global_forces(iforce_eval)%forces(:,iparticle)
       ELSE
          particles_mix%els(jparticle)%f(:)= particles_mix%els(jparticle)%f(:) + &
               factor* global_forces(iforce_eval)%forces(:,iparticle)
       END IF
    END DO
    ! Mixing Virial
    IF (virial_mix%pv_availability) THEN
       IF (overwrite) CALL zero_virial(virial_mix,reset=.FALSE.)
       virial_mix%pv_total      = virial_mix%pv_total + factor*virials(iforce_eval)%virial%pv_total
       virial_mix%pv_kinetic    = virial_mix%pv_kinetic + factor*virials(iforce_eval)%virial%pv_kinetic
       virial_mix%pv_virial     = virial_mix%pv_virial + factor*virials(iforce_eval)%virial%pv_virial
       virial_mix%pv_xc         = virial_mix%pv_xc + factor*virials(iforce_eval)%virial%pv_xc
       virial_mix%pv_fock_4c    = virial_mix%pv_fock_4c + factor*virials(iforce_eval)%virial%pv_fock_4c
       virial_mix%pv_constraint = virial_mix%pv_constraint + factor*virials(iforce_eval)%virial%pv_constraint
    END IF
    ! Deallocate map_index array
    IF (ASSOCIATED(map_index)) THEN
       DEALLOCATE(map_index)
    END IF

    ! Collect Requested Results info
    description = '[DIPOLE]'
    IF (overwrite) CALL cp_results_erase(results_mix)

    dip_exists = test_for_result(results=results(iforce_eval)%results,description=description)
    IF (dip_exists) THEN
      CALL get_results(results=results_mix,description=description,n_rep=nres)
      IF(.NOT.(nres<=1))CALL cp__a("mixed_environment_utils.F",121)
      dip_mix = 0.0_dp
      IF (nres==1) CALL get_results(results=results_mix,description=description,values=dip_mix)
      CALL get_results(results=results(iforce_eval)%results,description=description,n_rep=nres)
      CALL get_results(results=results(iforce_eval)%results,description=description,&
           values=dip_tmp,nval=nres)
      dip_mix = dip_mix + factor*dip_tmp
      CALL cp_results_erase(results=results_mix,description=description)
      CALL put_results(results=results_mix,description=description,values=dip_mix)
    END IF

  END SUBROUTINE mixed_map_forces

! *****************************************************************************
!> \brief performs mapping of the subsystems of different force_eval
!> \param mapping_section ...
!> \param natom ...
!> \param iforce_eval ...
!> \param nforce_eval ...
!> \param map_index ...
!> \author Teodoro Laino - University of Zurich [tlaino] - 05.2007
! *****************************************************************************
  SUBROUTINE get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, map_index)

    TYPE(section_vals_type), POINTER         :: mapping_section
    INTEGER, INTENT(IN)                      :: natom, iforce_eval, &
                                                nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: map_index

    CHARACTER(len=*), PARAMETER :: routineN = 'get_subsys_map_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ival, j, jval, k, &
                                                n_rep, n_rep_loc, n_rep_map, &
                                                n_rep_sys, tmp
    INTEGER, DIMENSION(:), POINTER           :: index_glo, index_loc, list
    LOGICAL                                  :: check, explicit
    TYPE(section_vals_type), POINTER         :: fragments_loc, fragments_sys, &
                                                map_force_ev, map_full_sys

    IF(.NOT.(.NOT.ASSOCIATED(map_index)))CALL cp__a("mixed_environment_utils.F",161)
    ALLOCATE(map_index(natom))
    CALL section_vals_get(mapping_section, explicit=explicit)
    IF (.NOT.explicit) THEN
       ! Standard Mapping.. subsys are assumed to have the same structure
       DO i = 1, natom
          map_index(i) = i
       END DO
    ELSE
       ! Mapping systems with different structures
       map_full_sys => section_vals_get_subs_vals(mapping_section,"FORCE_EVAL_MIXED")
       map_force_ev => section_vals_get_subs_vals(mapping_section,"FORCE_EVAL")
       CALL section_vals_get(map_full_sys, explicit=explicit)
       IF(.NOT.(explicit))CALL cp__a("mixed_environment_utils.F",174)
       CALL section_vals_get(map_force_ev, explicit=explicit, n_repetition=n_rep)
       IF(.NOT.(explicit))CALL cp__a("mixed_environment_utils.F",176)
       IF(.NOT.(n_rep==nforce_eval))CALL cp__a("mixed_environment_utils.F",177)
       DO i = 1, n_rep
          CALL section_vals_val_get(map_force_ev,"_SECTION_PARAMETERS_",i_rep_section=i,i_val=ival)
          IF (ival==iforce_eval) EXIT
       END DO
       IF(.NOT.(i<=nforce_eval))CALL cp__a("mixed_environment_utils.F",182)
       fragments_sys => section_vals_get_subs_vals(map_full_sys,"FRAGMENT")
       fragments_loc => section_vals_get_subs_vals(map_force_ev,"FRAGMENT",i_rep_section=i)
       !Perform few check on the structure of the input mapping section. as provided by the user
       CALL section_vals_get(fragments_loc, n_repetition=n_rep_loc)
       CALL section_vals_get(fragments_sys, explicit=explicit, n_repetition=n_rep_sys)
       IF(.NOT.(explicit))CALL cp__a("mixed_environment_utils.F",188)
       IF(.NOT.(n_rep_sys>=n_rep_loc))CALL cp__a("mixed_environment_utils.F",189)
       IF (n_rep_loc==0) THEN
          NULLIFY(list)
          ! We expect an easier syntax in this case..
          CALL section_vals_val_get(map_force_ev,"DEFINE_FRAGMENTS",i_rep_section=i,n_rep_val=n_rep_map)
          check = (n_rep_map/=0)
          IF(.NOT.(check))CALL cp__a("mixed_environment_utils.F",195)
          CALL section_vals_val_get(map_force_ev,"DEFINE_FRAGMENTS",i_rep_section=i,i_vals=list)
          IF(.NOT.(SIZE(list)>0))CALL cp__a("mixed_environment_utils.F",197)
          iatom = 0
          DO i = 1, SIZE(list)
             jval = list(i)
             DO j = 1, n_rep_sys
                CALL section_vals_val_get(fragments_sys,"_SECTION_PARAMETERS_",i_rep_section=j,i_val=tmp)
                IF (tmp==jval) EXIT
             END DO
             CALL section_vals_val_get(fragments_sys,"_DEFAULT_KEYWORD_",i_rep_section=j,i_vals=index_glo)
             DO k = 0, index_glo(2)-index_glo(1)
                iatom = iatom + 1
 IF(.NOT.(iatom<=natom))CALL cp__a("mixed_environment_utils.F",208)
                map_index(iatom) = index_glo(1)+k
             END DO
          END DO
          check = (iatom==natom)
          IF(.NOT.(check))CALL cp__a("mixed_environment_utils.F",213)
       ELSE
          ! General syntax..
          !Loop over the fragment of the force_eval
          DO i = 1, n_rep_loc
             CALL section_vals_val_get(fragments_loc,"_SECTION_PARAMETERS_",i_rep_section=i,i_val=ival)
             CALL section_vals_val_get(fragments_loc,"MAP",i_rep_section=i,i_val=jval)
             ! Index corresponding to the mixed_force_eval fragment
             DO j = 1, n_rep_sys
                CALL section_vals_val_get(fragments_sys,"_SECTION_PARAMETERS_",i_rep_section=j,i_val=tmp)
                IF (tmp==jval) EXIT
             END DO
             IF(.NOT.(j<=n_rep_sys))CALL cp__a("mixed_environment_utils.F",225)
             CALL section_vals_val_get(fragments_loc,"_DEFAULT_KEYWORD_",i_rep_section=i,i_vals=index_loc)
             CALL section_vals_val_get(fragments_sys,"_DEFAULT_KEYWORD_",i_rep_section=j,i_vals=index_glo)
             check = ((index_loc(2)-index_loc(1))==(index_glo(2)-index_glo(1)))
             IF(.NOT.(check))CALL cp__a("mixed_environment_utils.F",229)
             ! Now let's build the real mapping
             DO k = 0, index_loc(2)-index_loc(1)
                map_index(index_loc(1)+k) = index_glo(1)+k
             END DO
          END DO
       END IF
    END IF

  END SUBROUTINE get_subsys_map_index

END MODULE mixed_environment_utils
