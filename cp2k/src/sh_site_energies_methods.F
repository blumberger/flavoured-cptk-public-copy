
! *****************************************************************************
!> \brief 
!> \par History
!>      SG 
! *****************************************************************************
MODULE sh_site_energies_methods 


     USE kinds,                        ONLY: default_string_length,               &
                                             dp
     USE aom_overlapfinal,             ONLY:  connect_list2!,                      &
     !                                         build_s_tables

     USE aom_utils,                    ONLY: read_listcoeff,                       &
                                             read_listdiabcoeff,                  &
                                             read_list_tresp,                     &
                                             read_list_frz_couplings   

     USE sh_types,                    ONLY: sh_env_type

     USE input_section_types,          ONLY: section_vals_get_subs_vals,          &
                                             section_vals_type,                   &
                                             section_vals_list_get,               &
                                              section_vals_val_get

     USE cp_linked_list_val,           ONLY: &
         cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
         cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
         cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
         cp_sll_val_set_el_at, cp_sll_val_type


     USE cp_log_handling,                 ONLY: cp_get_default_logger,      &
                                                cp_logger_type,             &
                                                cp_to_string

     USE cp_output_handling,           ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_should_output
     USE molecule_kind_types,             ONLY: molecule_kind_type
     USE molecule_types_new,              ONLY: molecule_type,&
                                                get_molecule
     USE cell_types,                      ONLY: cell_type,                  &
                                                init_cell,                  &
                                                pbc
     USE atomic_kind_types,               ONLY: atomic_kind_type,           &
                                                get_atomic_kind
     USE atprop_types,                    ONLY: atprop_type
     USE bibliography,                    ONLY: Ewald1921,&
                                                cite_reference
     USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                                cp_subsys_type

     USE cp_para_types,                   ONLY: cp_para_env_type
     USE distribution_1d_types,           ONLY: distribution_1d_type
     USE dg_rho0_types,                   ONLY: dg_rho0_type
     USE dg_types,                        ONLY: dg_get,&
                                                dg_type
     USE exclusion_types,                 ONLY: exclusion_type
     USE ewald_environment_types,         ONLY: ewald_env_get,              &
                                                ewald_env_set,              &
                                                ewald_environment_type
     USE ewald_pw_methods,                ONLY: ewald_pw_grid_update
     USE ewald_pw_types,                  ONLY: ewald_pw_type,              &
                                                ewald_pw_get
     USE ewalds,                          ONLY: ewald_print
     USE fist_nonbond_force,              ONLY: bonded_correct_gaussian_fssh
     USE fist_environment_types,          ONLY: fist_environment_type,      &
                                                fist_env_get
     USE fist_energy_types,               ONLY: fist_energy_type
     USE fist_neighbor_list_control,      ONLY: list_control
     USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,         &
                                                neighbor_kind_pairs_type
     USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type,      &
                                                pos_type,                   &
                                                fist_nonbond_env_get
     USE kinds,                           ONLY: default_string_length,      &
                                                dp, dp_size, int_size
     USE input_constants,                 ONLY: do_fist_pol_none
     USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                                section_vals_type,          &
                                                section_vals_list_get,      &
                                                section_vals_val_get
     USE machine,                         ONLY: m_memory
     USE mathconstants,                   ONLY: fourpi,                     &
                                                oorootpi,                   &
                                                pi
     USE mathlib,                         ONLY: matvec_3x3
     USE manybody_eam,                    ONLY: density_nonbond
     USE manybody_potential,              ONLY: energy_manybody
     USE message_passing,                 ONLY: mp_sum
     USE mixed_environment_types,         ONLY: get_mixed_env,              &
                                                 mixed_env_retain,          &
                                                 mixed_environment_type
     USE mixed_energy_types,              ONLY: mixed_energy_type,&
                                                 mixed_force_type
     USE molecule_kind_types,             ONLY: molecule_kind_type
     USE molecule_types_new,              ONLY: molecule_type,              &
                                                get_molecule
     USE multipole_types,                 ONLY: multipole_type
     USE sh_types,                        ONLY: sh_env_type
     USE shell_potential_types,           ONLY: get_shell,                  &
                                                shell_kind_type
     USE structure_factor_types,          ONLY: structure_factor_type
     USE structure_factors,               ONLY: structure_factor_allocate,  &
                                                structure_factor_deallocate,&
                                                structure_factor_evaluate
     USE splines_types,                   ONLY: spline_data_p_type,         &
                                                spline_factor_type
     USE termination,                     ONLY: stop_memory
     USE particle_types,                  ONLY: particle_type
     USE pair_potential_coulomb,          ONLY: potential_coulomb
     USE pair_potential_types,            ONLY: pair_potential_pp_type,     &
                                                pair_potential_single_type, &
                                                siepmann_type,              &
                                                tersoff_type

     USE pw_grid_types,                   ONLY: pw_grid_type
     USE pw_poisson_types,                ONLY: do_ewald_none,              &
                                                do_ewald_ewald,             &
                                                do_ewald_dsf,               &
                                                do_ewald_pme,               &
                                                do_ewald_spme,              &
                                                do_ewald_fssh_ewald,        &
                                                do_ewald_fssh_dsf
     USE pw_pool_types,                   ONLY: pw_pool_type
     USE virial_types,                    ONLY: virial_type


#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='sh_site_energies_methods'   

     PUBLIC :: ct_block_energies_forces, &
               add_energies_forces_mem_saver_V2, &
               add_energies_forces, &
               create_ext_press_list, &
               add_energies_forces_mem_saver, &
               multi_block_energies_forces, &
               add_coulomb_interactions, &
               single_phase_energies_forces

CONTAINS

! *****************************************************************************
!>  \brief This creates the logical list of active molecules ext_pres that is
!         used to calculate kinetic energy for the active region in case of
!         decoherence and rescaling. Note in case AOM is used ext_press is
!         created by reading the AOM coeff from external file. 
!         It also creates active_system_info which is an 2D array with info
!         about initial and final atoms of a give molecule. It is useful to
!         consider multi-phase systems where n atoms may vary between different
!         molecules
!>  \author SG 
!> \date 02/02/20
! *****************************************************************************
subroutine create_ext_press_list(molecule_kind_set, local_molecules, molecule_set, natom, my_mol_index, & 
                                ext_pres, active_system_info, full_system_info)

     TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
     TYPE(distribution_1d_type), POINTER                :: local_molecules
     TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
     INTEGER, INTENT(IN)                                :: natom
     INTEGER, DIMENSION(:), POINTER                     :: my_mol_index !this are the decomp list indeces
     LOGICAL, DIMENSION(:), POINTER,  INTENT(INOUT)     :: ext_pres
     INTEGER, DIMENSION(:,:), ALLOCATABLE,               &
                                              OPTIONAL  :: full_system_info
     INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)  :: active_system_info

     CHARACTER(len=*), PARAMETER :: routineN ='create_ext_press_list', &
     routineP = moduleN//':'//routineN

     INTEGER                                            :: imol, nkind,  nmol_per_kind,  ikind, first_atom, &
                                                            last_atom, aoi, my_mol, j
     TYPE(molecule_kind_type), POINTER                  :: molecule_kind
     TYPE(molecule_type), POINTER                       :: molecule
     INTEGER                                            :: total_num_mol !total number of
                                                               !mols in the system 
     INTEGER                                            :: count_
     !> this array is still deallocated 
     ALLOCATE(ext_pres(natom))

     !>active_system_info contains info about the number of initial and final atom for each mol. 
     !active_system_info(1, start index)
     !active_system_info(2, final index)
     ALLOCATE(active_system_info(2,SIZE(my_mol_index)))

     !initialization
     ext_pres(:) = .FALSE.    

     nkind = SIZE ( molecule_kind_set )

     total_num_mol = 0
     count_ = 0
     MOL: DO ikind = 1, nkind !> loop over kind of molecules (the kind depends on the
                               !> different names given to the residues in the
                               !topology file)
           nmol_per_kind = local_molecules % n_el ( ikind )
           DO imol = 1, nmol_per_kind !loop over molecules of a given kind
              my_mol = local_molecules % list ( ikind ) % array ( imol )
              molecule => molecule_set ( my_mol )
              molecule_kind => molecule % molecule_kind  
              total_num_mol = total_num_mol + 1
              DO j = 1, SIZE(my_mol_index)
                 IF ( my_mol_index(j)==my_mol) THEN 
                    count_ = count_ +1
                    CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom )
                    active_system_info(1,count_) =  first_atom
                    active_system_info(2,count_) =  last_atom
                    DO aoi=first_atom, last_atom !loop over atoms of interest
                       ext_pres(aoi) = .TRUE.
                    END DO
                 END IF   
              END DO
           END DO 
     END DO MOL

 

     !>full_system_info contains info about the number of atoms for each mol and
     !whether this mol is active or not. 
     !full_system_info(1, start index)
     !full_system_info(2, final index)
     !full_system_info(3, 1 if active, 0 if not)

     IF (PRESENT(full_system_info)) THEN 
         ALLOCATE(full_system_info(3,total_num_mol))

         full_system_info(:,:) = 0
         count_ = 0
         nkind = SIZE ( molecule_kind_set )
         MOL1: DO ikind = 1, nkind
               nmol_per_kind = local_molecules % n_el ( ikind )
               DO imol = 1, nmol_per_kind !loop over molecules of a given kind
                  my_mol = local_molecules % list ( ikind ) % array ( imol )
                  molecule => molecule_set ( my_mol )
                  molecule_kind => molecule % molecule_kind
                  CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom )
                  !natom_mol = (last_atom-first_atom)+1
                  count_ = count_+1 
                  full_system_info(1,count_) =  first_atom
                  full_system_info(2,count_) =  last_atom
                  DO j = 1, SIZE(my_mol_index)
                     IF ( my_mol_index(j)==my_mol) THEN
                        full_system_info(3,count_) =  1
                     END IF
                  END DO
               END DO
         END DO MOL1
     END IF
              
end subroutine create_ext_press_list

SUBROUTINE single_phase_energies_forces(adiab_section,energies_1_mol_list, energies, &
                                    intra_forces_global, global_forces, my_energies, &
                                    active_system_info, mixed_env, my_itimes, my_mol_index, &
                                    fe_block_keyword, cr_block_keyword) 

   TYPE(section_vals_type), POINTER          :: adiab_section
   TYPE(mixed_environment_type), &
      POINTER                                :: mixed_env
   TYPE(mixed_force_type), DIMENSION(:), &
     POINTER                                 :: global_forces
   REAL(KIND=dp), DIMENSION(:), POINTER      :: energies
   REAL(KIND=dp), DIMENSION(:), POINTER      :: my_energies
   REAL(KIND=dp), ALLOCATABLE                :: intra_forces_global(:,:,:) 
   REAL(KIND=dp), ALLOCATABLE                :: energies_1_mol_list(:,:)
   INTEGER, POINTER                          :: my_itimes !CHANGE_AC
   INTEGER, DIMENSION(:,:), INTENT(IN)       :: active_system_info
   LOGICAL, INTENT(IN)                       :: fe_block_keyword, cr_block_keyword

   CHARACTER(len=*), PARAMETER :: routineN ='single_phase_energies_forces', &
     routineP = moduleN//':'//routineN

   INTEGER :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
   INTEGER :: init_atom_mol_j, final_atom_mol_j, atoms_per_mol_j
   INTEGER                                   :: i,j
   INTEGER                                   :: start_i,end_i,start_j,end_j,atoms
   TYPE(section_vals_type), POINTER          :: site_energy_section
   TYPE(section_vals_type), POINTER          :: ct_block, fect_block  
   CHARACTER(LEN=default_string_length)           :: fe_keyword  
   INTEGER                                   :: total_states, state, nadiab, natom
   INTEGER, DIMENSION(:,:), ALLOCATABLE      :: intra_forces_info
   INTEGER                                   :: count_
   INTEGER, DIMENSION(:), INTENT(IN)         :: my_mol_index !this is the decomp index

   NULLIFY(site_energy_section,ct_block, fect_block)
   NULLIFY( my_energies)
   ! TOTAL number states
   CALL section_vals_val_get(adiab_section,"NUMBER_DIABATIC_STATES",&
                        i_val=nadiab)

   site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
   ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")

      !> this are the output add-sub energies and forces
   ALLOCATE(my_energies(nadiab))
   natom = SIZE(global_forces(1)%forces(:,:), 2)

   !> we do F allocation just at the first step and never again just to NOT
   !allocate and deallocate this big array at every step. NOTE that this array is
   !now never deallocated and just deallocated by the compiler at the very end.
   !One should probably write a subroutine that deallocate all allocatable just
   !at the very end of cp2k
    if(my_itimes .eq. 0) then
       ALLOCATE(mixed_env%F(3,natom,nadiab))
    end if
  
   !> construct a intra_forces_info array with the index of the atoms that are
   !part of the intra_forces_global array that stores only the active atoms.
   !intra_forces_info is therefore an array with the same shape as
   !active_system_info but with sequential numbering of the atoms.
   ALLOCATE(intra_forces_info(2,SIZE(active_system_info,2)))
   atoms = 0
   do i=1,SIZE(active_system_info,2)
         !> find numer of atoms in mol i
         init_atom_mol_i = active_system_info(1,i)
         final_atom_mol_i = active_system_info(2,i)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
         atoms = atoms + atoms_per_mol_i
         intra_forces_info(1,i) = (atoms - atoms_per_mol_i)+1
         intra_forces_info(2,i) = atoms
   end do

   if ((fe_block_keyword) .and. (cr_block_keyword)) then
      !> some checks  N*(N-1)+N+1 is the total number of states
      total_states = SIZE(my_mol_index)*SIZE(my_mol_index) + 1
   else if ((.not. fe_block_keyword) .and. (cr_block_keyword)) then 
      !> some checks  N*(N-1)+1 is the total number of states
      total_states = SIZE(my_mol_index)*SIZE(my_mol_index) - SIZE(my_mol_index) + 1
   else if ((fe_block_keyword) .and. (.not. cr_block_keyword)) then
      !> some checks  N*(N-1)+N is the total number of states
      total_states = SIZE(my_mol_index)*SIZE(my_mol_index)
   else 
      total_states = SIZE(my_mol_index)*SIZE(my_mol_index)
   end if

   IF (total_states .ne. nadiab) THEN 
         CPABORT("NUMBER_DIABATIC_STATES different from total n states")
   END IF


   !> loop over states: note that it does not matter if donor molecules comes
   !after acceptor in the xyz or viceversa, or the molecules are DDDAAA or ADADA
   !as the index of each molecule is referred to DECOMP file
   state = 0
   do i=1,SIZE(my_mol_index)
      do j=1,SIZE(my_mol_index)

         if (i .ne. j) then

         state = state + 1
         !> we create site energies from the 2 force_evals sections: the first
         !neutral, the second charged (with donor and acc Froce-fields).
         my_energies(state) = energies(1) + energies_1_mol_list(2,i) - energies_1_mol_list(1,i) &
            + energies_1_mol_list(3,j) - energies_1_mol_list(1,j)  

         !> we allocate diagonal forces array. All the atoms for a give state i
         !get the value of the neutral force_field that contains: intra +
         !nonbonded interactions
         mixed_env%F(:,:,state) = global_forces(1)%forces
         
         !> find index mol i (which is state i) of the pair
         init_atom_mol_i = active_system_info(1,i)
         final_atom_mol_i = active_system_info(2,i)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

         !> find index mol j (which is state j) of the pair
         init_atom_mol_j = active_system_info(1,j)
         final_atom_mol_j = active_system_info(2,j)
         atoms_per_mol_j = (final_atom_mol_j-init_atom_mol_j)+1

         !> since  mixed_env%F stores all the atoms the numbering is a
         !bit more triky and relies on the "global numbering" (which is an info
         !contained in active_system_info). Whereas intra_forces_global stores
         !just the active atoms and so the numbering should be sequential and
         !this info is included in intra_forces_info

        start_i = intra_forces_info(1,i)
        end_i = intra_forces_info(2,i)

       !> first molecule i of the pair
        mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,1,start_i:end_i) - &
                           intra_forces_global(1,1,start_i:end_i)
        mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,2,start_i:end_i) - &
                           intra_forces_global(1,2,start_i:end_i)
        mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,3,start_i:end_i) - &
                           intra_forces_global(1,3,start_i:end_i)

        start_j = intra_forces_info(1,j)
        end_j = intra_forces_info(2,j)

       !> second mol j of the pair
        mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(3,1,start_j:end_j) - &
                           intra_forces_global(1,1,start_j:end_j)
        mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(3,2,start_j:end_j) - &
                           intra_forces_global(1,2,start_j:end_j)
        mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(3,3,start_j:end_j) - &
                           intra_forces_global(1,3,start_j:end_j)


           end if
      end do
   end do

!                   open(913,file="diagonal_forces_state1.txt",position="append")
!                   do i = 1,size(mixed_env%F(3,:,1))
!                   write(913,*) mixed_env%F(3,i,1)
!                   end do
!                   close(913)
!

  !> THE FE BLOCK (always after the CT Block)
  IF (fe_block_keyword) THEN
      !energies_1_mol_list like in the CT case
  
      !> loop over each (active) state starting from where CT were left
      count_ = 0
      DO i = state+1, total_states !loop over num of diabatic states until the end
         count_ = count_ + 1
         !> get index actual mol excited BUT in the base 1,2,3..
         !ind1 = exc_ind(count_)
   
         !> we create site energies from the 2 force_evals sections: the first
         !neutral, the second charge and the THIRD excited.
         my_energies(i) = energies(1) + energies_1_mol_list(4,count_) - energies_1_mol_list(1,count_)

         !> we allocate diagonal forces array. All the atoms for a give state i
         !assume the value of the neutral force_field that contains: intra +
         !nonbonded interactions
         mixed_env%F(:,:,i) = global_forces(1)%forces
         
         !> find index mol i (which is state i) of the pair
         init_atom_mol_i = active_system_info(1,count_)
         final_atom_mol_i = active_system_info(2,count_)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
   
   
         !> since intra_forces_global stores only active atoms the numbering is a
         !bit more triky and relies on the n of atoms per mol (which is an info
         !contained in active_system_info)
         start_i = intra_forces_info(1,count_)
         end_i = intra_forces_info(2,count_)

         mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(4,1,start_i:end_i) - &
                            intra_forces_global(1,1,start_i:end_i)
         mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(4,2,start_i:end_i) - &
                            intra_forces_global(1,2,start_i:end_i)
         mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(4,3,start_i:end_i) - &
                            intra_forces_global(1,3,start_i:end_i)
      END DO
   END IF

   DEALLOCATE(intra_forces_info)

!open(912,file="site-energies.txt",position="append")
!write(912,*) "1st CT-state, anion opt-config:", my_energies(1)
!write(912,*) "1st XT-state, exciton opt-config:", my_energies(381)
!write(912,*) "total energy, anion opt-config:", energies(1)
!close(912)


END SUBROUTINE single_phase_energies_forces

 ! *****************************************************************************
!> \brief This subroutine does the addition and subtraction scheme for the 
!         CT, FE, CR states  
!> \note This subroutine gives non-degenerate site_energies (degeneracy may be
!        introduce in the fssh specific subroutines)
!> \par History 
!> \author SG 
! *****************************************************************************
SUBROUTINE multi_block_energies_forces(adiab_section,energies_1_mol_list, energies, &
                                    intra_forces_global, global_forces, my_energies, &
                                    active_system_info, mixed_env, my_itimes, my_mol_index, &
                                    fe_block_keyword, cr_block_keyword) 

   TYPE(section_vals_type), POINTER          :: adiab_section
   TYPE(mixed_environment_type), &
      POINTER                                :: mixed_env
   TYPE(mixed_force_type), DIMENSION(:), &
     POINTER                                 :: global_forces
   REAL(KIND=dp), DIMENSION(:), POINTER      :: energies
   REAL(KIND=dp), DIMENSION(:), POINTER      :: my_energies
   REAL(KIND=dp), ALLOCATABLE                :: intra_forces_global(:,:,:) 
   REAL(KIND=dp), ALLOCATABLE                :: energies_1_mol_list(:,:)
   INTEGER, POINTER                          :: my_itimes !CHANGE_AC
   INTEGER, DIMENSION(:,:), INTENT(IN)       :: active_system_info
   LOGICAL, INTENT(IN)                       :: fe_block_keyword, cr_block_keyword

   CHARACTER(len=*), PARAMETER :: routineN ='multi_block_energies_forces', &
     routineP = moduleN//':'//routineN

   INTEGER :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
   INTEGER :: init_atom_mol_j, final_atom_mol_j, atoms_per_mol_j
   INTEGER                                   :: i,j, ind1, ind2!, indfirst,indsec
   INTEGER                                   :: start_i,end_i,start_j,end_j,atoms
   TYPE(section_vals_type), POINTER          :: site_energy_section
   TYPE(section_vals_type), POINTER          :: ct_block, fect_block  
   CHARACTER(LEN=default_string_length)           :: fe_keyword  
   INTEGER, DIMENSION(:), POINTER            :: donor_index, acceptor_index !, exc_index
   INTEGER, DIMENSION(:), allocatable            :: exc_index
   INTEGER                                   :: total_states, state,nadiab, natom
   INTEGER, DIMENSION(:,:), ALLOCATABLE      :: intra_forces_info
   INTEGER, DIMENSION(:), ALLOCATABLE        :: first_ind, second_ind, exc_ind
   INTEGER                                   :: count_
   INTEGER, DIMENSION(:), INTENT(IN)         :: my_mol_index !this is the decomp index

   NULLIFY(site_energy_section,ct_block, fect_block)
   NULLIFY( my_energies)
   ! TOTAL number states
   CALL section_vals_val_get(adiab_section,"NUMBER_DIABATIC_STATES",&
                        i_val=nadiab)

   site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
   ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")

   CALL section_vals_val_get(ct_block,"INDEX_DONOR_MOLS",&
                        i_vals=donor_index)
   CALL section_vals_val_get(ct_block,"INDEX_ACCEPTOR_MOLS",&
                        i_vals=acceptor_index)

   !open(911, file="exciton_opt_geom.txt", position="append")
   !write(911,*) "neutral FF energy is:", energies_1_mol_list(1,1)
   !write(911,*) "anion FF energy is:", energies_1_mol_list(2,1)
   !write(911,*) "cation FF energy is:", energies_1_mol_list(2,1)             
   !write(911,*) "exciton FF energy is:", energies_1_mol_list(4,1)
   !close(911)


   ind1 = donor_index(1) 
   do i=2,SIZE(donor_index)
      if (donor_index(i) .le. ind1) THEN 
        CPABORT("Donor mol must be ascending")
      end if
   end do
   ind1 = acceptor_index(1) 
   do i=2,SIZE(acceptor_index)
      if (acceptor_index(i) .le. ind1) THEN 
        CPABORT("Acceptor mol must be ascending")
      end if
   end do

   !> this are the output add-sub energies and forces
   ALLOCATE(my_energies(nadiab))
   natom = SIZE(global_forces(1)%forces(:,:), 2)

   !> we do F allocation just at the first step and never again just to NOT
   !allocate and deallocate this big array at every step. NOTE that this array is
   !now never deallocated and just deallocated by the compiler at the very end.
   !One should probably write a subroutine that deallocate all allocatable just
   !at the very end of cp2k
    if(my_itimes .eq. 0) then
       ALLOCATE(mixed_env%F(3,natom,nadiab))
    end if
  
   !> construct a intra_forces_info array with the index of the atoms that are
   !part of the intra_forces_global array that stores only the active atoms.
   !intra_forces_info is therefore an array with the same shape as
   !active_system_info but with sequential numbering of the atoms.
   ALLOCATE(intra_forces_info(2,SIZE(active_system_info,2)))
   atoms = 0
   do i=1,SIZE(active_system_info,2)
         !> find numer of atoms in mol i
         init_atom_mol_i = active_system_info(1,i)
         final_atom_mol_i = active_system_info(2,i)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
         atoms = atoms + atoms_per_mol_i
         intra_forces_info(1,i) = (atoms - atoms_per_mol_i)+1
         intra_forces_info(2,i) = atoms
   end do

!print *, "intra_forces_info(1,i)", intra_forces_info(1,:)   
!print *, "intra_forces_info(2,i)", intra_forces_info(2,:)   
   !get the CT index in the 1,2,3... base (active region basis)
   !independently on the actual molecule
   !numbering
   ! FINDLOC function not supported so loops instead
   ALLOCATE(second_ind(SIZE(acceptor_index)))
   do i=1,SIZE(acceptor_index)
     do j=1,SIZE(my_mol_index) 
         if (acceptor_index(i) .eq. my_mol_index(j)) then
            second_ind(i) = j 
            EXIT
         end if 
     end do 
   end do
   
   ALLOCATE(first_ind(SIZE(donor_index)))
   do i=1,SIZE(donor_index)
     do j=1,SIZE(my_mol_index) 
         if (donor_index(i) .eq. my_mol_index(j)) then
            first_ind(i) = j 
            EXIT
         end if 
     end do 
   end do
   
    if (fe_block_keyword) then
        fect_block => section_vals_get_subs_vals(site_energy_section,"FE_CT_BLOCK")
        CALL section_vals_val_get(fect_block,"FE_IS", c_val=fe_keyword)

        SELECT CASE (fe_keyword)
        CASE("DONOR")
            ALLOCATE(exc_index(SIZE(donor_index)))
            exc_index = donor_index
        CASE("ACCEPTOR")
            ALLOCATE(exc_index(SIZE(acceptor_index)))
            exc_index = acceptor_index
        END SELECT
        ind1 = exc_index(1) 
        do i=2,SIZE(exc_index)
            if (exc_index(i) .le. ind1) THEN 
                CPABORT("EXC mol must be ascending")
            end if
        end do
    end if 

   if ((fe_block_keyword) .and. (cr_block_keyword)) then
      !> some checks  N*N+N+1 is the total number of states
      total_states = SIZE(donor_index)*SIZE(acceptor_index) + SIZE(exc_index) +1
   else if ((.not. fe_block_keyword) .and. (cr_block_keyword)) then 
      !> some checks  N*N+1 is the total number of states
      total_states = SIZE(donor_index)*SIZE(acceptor_index) +1
   else if ((fe_block_keyword) .and. (.not. cr_block_keyword)) then
      !> some checks  N*N+N+1 is the total number of states
      total_states = SIZE(donor_index)*SIZE(acceptor_index) + SIZE(exc_index)
   else 
      total_states = SIZE(donor_index)*SIZE(acceptor_index)
   end if

   IF (total_states .ne. nadiab) THEN 
         CPABORT("NUMBER_DIABATIC_STATES different from total n states")
   END IF


   !> loop over states: not that it does not matter if donor molecules comes
   !after acceptor in the xyz or viceversa, or the molecules are DDDAAA or ADADA
   !as the index of each molecule is referred to DECOMP file
   state = 0
   do i=1,SIZE(first_ind)
      ind1 = first_ind(i)
      do j=1,SIZE(second_ind)
         ind2 = second_ind(j)
         state = state + 1
!print *, "STATES", ind1,ind2
!print *, "STATES:", state
         !> we create site energies from the 2 force_evals sections: the first
         !neutral, the second charged (with donor and acc Froce-fields).
         my_energies(state) = energies(1) + &
                  energies_1_mol_list(2,ind1) - energies_1_mol_list(1,ind1) + & 
                  energies_1_mol_list(2,ind2) - energies_1_mol_list(1,ind2)  

         !> we allocate diagonal forces array. All the atoms for a give state i
         !get the value of the neutral force_field that contains: intra +
         !nonbonded interactions
         mixed_env%F(:,:,state) = global_forces(1)%forces
         
         !> find index mol i (which is state i) of the pair
         init_atom_mol_i = active_system_info(1,ind1)
         final_atom_mol_i = active_system_info(2,ind1)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

         !> find index mol j (which is state j) of the pair
         init_atom_mol_j = active_system_info(1,ind2)
         final_atom_mol_j = active_system_info(2,ind2)
         atoms_per_mol_j = (final_atom_mol_j-init_atom_mol_j)+1

         !> since  mixed_env%F stores all the atoms the numbering is a
         !bit more triky and relies on the "global numbering" (which is an info
         !contained in active_system_info). Whereas intra_forces_global stores
         !just the active atoms and so the numbering should be sequential and
         !this info is included in intra_forces_info
        start_i = intra_forces_info(1,ind1)
        end_i = intra_forces_info(2,ind1)
!print *, "index don", init_atom_mol_i, final_atom_mol_i
!print *, "NEW", start_i, end_i
        !> first molecule i of the pair
        mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,1,start_i:end_i) - &
                           intra_forces_global(1,1,start_i:end_i)
        mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,2,start_i:end_i) - &
                           intra_forces_global(1,2,start_i:end_i)
        mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,3,start_i:end_i) - &
                           intra_forces_global(1,3,start_i:end_i)

        start_j = intra_forces_info(1,ind2)
        end_j = intra_forces_info(2,ind2)
!print*, "index acc", init_atom_mol_j, final_atom_mol_j
!print *, "NEW", start_j, end_j
        !> second mol j of the pair
        mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,1,start_j:end_j) - &
                           intra_forces_global(1,1,start_j:end_j)
        mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,2,start_j:end_j) - &
                           intra_forces_global(1,2,start_j:end_j)
        mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,3,start_j:end_j) - &
                           intra_forces_global(1,3,start_j:end_j)

      end do
   end do
  !> THE FE BLOCK (always after the CT Block)
  IF (fe_block_keyword) THEN
      !exc_ind  is in the 1,2,3..  base (active mol base) to be consisent with
      !energies_1_mol_list like in the CT case
      ALLOCATE(exc_ind(SIZE(exc_index)))
      do i=1,SIZE(exc_index)
        do j=1,SIZE(my_mol_index) 
            if (exc_index(i) .eq. my_mol_index(j)) then
               exc_ind(i) = j 
               EXIT
            end if 
        end do 
      end do
   
      if (cr_block_keyword) then
         ! in case CR is present
         total_states = nadiab-1
      else 
         total_states = nadiab
      end if
   
   !print *, "STATE to BEFORE FE", state
      !> loop over each (active) state starting from where CT were left
      count_ = 0
      DO i = state+1, total_states !loop over num of adiabatic states until the end
   !print *, "GLOB_STATE", i
         count_ = count_ + 1
         !> get index actual mol excited BUT in the base 1,2,3..
         ind1 = exc_ind(count_)
   
   !print *, "FE STATE", ind1
         !> we create site energies from the 2 force_evals sections: the first
         !neutral, the second charge and the THIRD excited.
         my_energies(i) = energies(1) + energies_1_mol_list(3,ind1) - energies_1_mol_list(1,ind1) !+ sh_env%e_offset
         !>WTP_4.2: double check sh_env%e_offset term!! 1) default = 0.0 2) input value correctly passed here

         !> we allocate diagonal forces array. All the atoms for a give state i
         !assume the value of the neutral force_field that contains: intra +
         !nonbonded interactions
         mixed_env%F(:,:,i) = global_forces(1)%forces
         
         !> find index mol i (which is state i) of the pair
         init_atom_mol_i = active_system_info(1,ind1)
         final_atom_mol_i = active_system_info(2,ind1)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
   
   
         !> since intra_forces_global stores only active atoms the numbering is a
         !bit more triky and relies on the n of atoms per mol (which is an info
         !contained in active_system_info)
         start_i = intra_forces_info(1,ind1)
         end_i = intra_forces_info(2,ind1)
   !print *, "index exc", init_atom_mol_i, final_atom_mol_i
   !print *, "NEW", start_i, end_i
         mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(3,1,start_i:end_i) - &
                            intra_forces_global(1,1,start_i:end_i)
         mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(3,2,start_i:end_i) - &
                            intra_forces_global(1,2,start_i:end_i)
         mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) = &
                            mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) + &
                            intra_forces_global(3,3,start_i:end_i) - &
                            intra_forces_global(1,3,start_i:end_i)
      END DO
      DEALLOCATE(exc_ind)
   END IF

  !> THE CR BLOCK (always after the last Block)
  !> WTP_4.2: double check...
   if (cr_block_keyword) then
      ! in case CR is present
      total_states = nadiab
      my_energies(total_states) = energies(1)
      mixed_env%F(:,:,total_states) = global_forces(1)%forces
   end if
   

   DEALLOCATE(intra_forces_info)
   DEALLOCATE(first_ind)
   DEALLOCATE(second_ind)
   if (fe_block_keyword) then
       DEALLOCATE(exc_index) 
   end if

END SUBROUTINE multi_block_energies_forces

! *****************************************************************************
!> \brief  This subroutine will apply the electrostatic interactions to energies
!>         and forces.
!> \param  mixed_env => The mixed force eval section
!>         fist_env  => The fist env which needs coulomb interactions
!>                      calculating.
!>         my_energies => The site_energies (should be in sh_env)
!>         active_system_info => The molecules that are active (should be in
!>                               sh_env)
!> \date   10/06/2020
!> \author Matt
!> \note   The important variables in this subroutine are the F array in the
!>         mixed_env and the my_energies. These are the force of state i (i.e.
!>         where mol i is charged all other are neutral) and site-energies or
!>         the energy of state i.
!>
!>         Shape of F array:     (3, natom, nstate)
!>         Shape of my_energies: (nstate)
! *****************************************************************************
SUBROUTINE add_coulomb_interactions (mixed_env, fist_env, my_energies,         &
                                     iforce_eval, my_mol_index)
   TYPE(fist_environment_type), POINTER            :: fist_env
   TYPE(mixed_environment_type), POINTER           :: mixed_env
   REAL(KIND=dp), DIMENSION(:), POINTER            :: my_energies
   INTEGER                                         :: iforce_eval
   INTEGER, DIMENSION(:), INTENT(IN)               :: my_mol_index !index of
                                                    ! active mols

   CHARACTER(len=*), PARAMETER :: routineN ='add_coulomb_interactions',        &
                                  routineP = moduleN//':'//routineN

   TYPE(atomic_kind_type), POINTER                 :: atomic_kind_set(:),    &
                                                      atomic_kind
   !TYPE(atprop_type), POINTER                      :: atprop_env
   TYPE(cell_type), POINTER                        :: cell
   TYPE(cp_logger_type), POINTER                   :: logger
   TYPE(cp_subsys_type), POINTER                   :: subsys
   TYPE(distribution_1d_type), POINTER             :: local_particles,       &
                                                      local_molecules
   TYPE(ewald_environment_type), POINTER           :: ewald_env
   TYPE(ewald_pw_type), POINTER                    :: ewald_pw 
   TYPE(fist_energy_type), POINTER                 :: thermo
   TYPE(fist_nonbond_env_type), POINTER            :: fist_nonbond_env
   TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
   TYPE(molecule_type), POINTER                    :: molecule_set(:)
   TYPE(section_vals_type), POINTER                :: fist_env_section,      &
                                                      mm_section
   TYPE(particle_type), POINTER                    :: particle_set(:)
   TYPE(virial_type), POINTER                      :: virial
   LOGICAL                                         :: shell_present,         &
                                                      use_virial,            &
                                                      do_multipoles,         &
                                                      do_print_charges
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: fR_neut, f_bond_C,     &
                                                      f_bond_n
   !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: fR_corr
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: e_bond_C, e_bond_N
   REAL(KIND=dp), DIMENSION(:), POINTER            :: charge_state_q,        &
                                                      neutral_charges,       &
                                                      ! set only once in sh part
                                                      charged_charges
   REAL(KIND=dp)                                   :: ES_pot, q, vR_coulomb, &
                                                      self_pot_corr,         &
                                                      neut_pot_corr,         &
                                                      neut_charge,           &
                                                      sum_e_bond_n
   ! Should probably put this in the sh part and set only once
   INTEGER, ALLOCATABLE, DIMENSION(:)              :: states
   INTEGER                                         :: nstate, nat, i, j, iw, &
                                                      handle, ist, iat,      &
                                                      ikind, n_of_kind,      &
                                                      ipart, nkind, do_ipol, &
                                                      nmol_per_kind,         &
                                                      nlocal_particles,      &
                                                      nparticle_kind,        &
                                                      start_atom, end_atom,  &
                                                      my_mlen, file_num,     &
                                                      ewald_type, imol

   IF (mixed_env%sh_env%finished_coulomb_step) RETURN

   !------- Decide what type of ES to handle ---------------------------------!
   ewald_type = mixed_env%ewald_types(iforce_eval)
   SELECT CASE(ewald_type)
   CASE (do_ewald_none)
      ! No need to do any coulomb calcs
      RETURN
   CASE (do_ewald_fssh_ewald)
      ! Don't need to do anything here
   CASE (do_ewald_ewald)
      CPABORT("Coding Error! Only do_ewald_fssh_... should get to here!")
   CASE (do_ewald_fssh_dsf)
      ! Don't need to do anything here
   CASE (do_ewald_dsf)
      CPABORT("Coding Error! Only do_ewald_fssh_... should get to here!")
   CASE (do_ewald_pme)
      CPABORT("PME not supported in surface hopping currently.")
   CASE (do_ewald_spme)
      CPABORT("SPME not supported in surface hopping currently.")
   CASE default
      CPABORT("Coding Error! Unrecognised option for ES method in FSSH.")
   END SELECT
   !--------------------------------------------------------------------------!

   PRINT *, "Adding Coulomb Interactions"

   !----------------------- Init Stuff ---------------------------------------!
   CALL timeset(routineN,handle)
   NULLIFY (logger, thermo, atomic_kind_set, atomic_kind, local_molecules,    &
            ewald_env, charge_state_q)
   logger => cp_get_default_logger()

   CALL fist_env_get(fist_env, fist_nonbond_env=fist_nonbond_env,             &
                     particle_set=particle_set,                               &
                     atomic_kind_set=atomic_kind_set,                         &
                     subsys=subsys,                                           &
                     molecule_kind_set=molecule_kind_set,                     &
                     local_molecules=local_molecules,                         &
                     molecule_set=molecule_set,                               &
                     ewald_env=ewald_env,                                     &
                     local_particles=local_particles,                         &
                     shell_model=shell_present,                               &
                     cell=cell,                                               &
                     thermo=thermo,                                           &
                     ewald_pw=ewald_pw)

   CALL ewald_env_get(ewald_env=ewald_env, do_ipol=do_ipol,                   &
                      do_multipoles=do_multipoles)

   CALL cp_subsys_get(subsys, virial=virial)
   use_virial = virial%pv_availability .AND. (.NOT.virial%pv_numer)
   
   ! Get some useful sizes
   nstate = SIZE(my_energies)
   nat = SIZE(particle_set) 
   nkind = SIZE ( molecule_kind_set )
   nmol_per_kind = local_molecules % n_el(1)
   nparticle_kind = SIZE(atomic_kind_set)
   nlocal_particles = 0
   DO ikind = 1, nparticle_kind
       nlocal_particles = nlocal_particles + local_particles%n_el(ikind)
   ENDDO


   NULLIFY(mm_section, fist_env_section)
   CALL fist_env_get(fist_env, input=fist_env_section)
   mm_section => section_vals_get_subs_vals(fist_env_section,"MM")
   iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%DERIVATIVES",&
                             extension=".mmLog")

   ! Some error checking
   CPASSERT(ASSOCIATED(mixed_env%F))
   CPASSERT(SIZE(mixed_env%F, 1) == 3)
   CPASSERT(SIZE(mixed_env%F, 2) == nat)
   CPASSERT(SIZE(SHAPE(mixed_env%F)) == 3)
   IF (nkind /= 1) &
      CPABORT("ES for FSSH accepts only 1 molecule kind not"//cp_to_string(nkind))
   !IF (nmol_per_kind /= SIZE(mixed_env%F, 3)) CPABORT("# mols /= # charge states")

   ! Checks for implementation first
   IF (do_multipoles) THEN
      CPABORT("Mulitpoles not implemented within Surface Hopping")
   END IF
   IF (do_ipol /= do_fist_pol_none) THEN
      CPABORT("Polarisable force-fields not implemented within Surface Hopping")
   END IF
   IF (shell_present) THEN
      CPABORT("Core-Shell model not yet implemented within Ewald sums.")
   END IF
   IF (use_virial) THEN
      CPABORT("Virial not implemented within Surface Hopping")
   END IF

   ! Initialize ewald grids -if required
   IF (ewald_type .eq. do_ewald_fssh_ewald) THEN
       CALL init_cell(cell)
       CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat)
   END IF
   !--------------------------------------------------------------------------!


   !--------------------- Allocate and set charges ---------------------!
   ALLOCATE(neutral_charges(nat))
   ALLOCATE(charged_charges(nat))

   ! Allocate the fist_nonbond_env charge array
   IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
      CPASSERT(SIZE(fist_nonbond_env%charges) == nat)
   ELSE
      ALLOCATE(fist_nonbond_env%charges(nat))
   END IF

   ! Store the neutral charges locally
   DO ikind=1, SIZE(atomic_kind_set)
      ! Get charge and num of particles in a given kind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind, qeff=q)
      n_of_kind = local_particles%n_el(ikind)

      DO ipart=1, n_of_kind
         ! Get index of atom and set the charge
         iat = local_particles%list(ikind)%array(ipart)
         neutral_charges(iat) = q
      END DO
   END DO

   CALL ewald_env_get(ewald_env, charge_state_q=charge_state_q,              &
                      print_charges=do_print_charges)
   !-------------------------------------------------------------------------!
   

   !------------------ Real Space -------------------------------------------!
   ALLOCATE(fR_neut(3, nat))
   ALLOCATE(states(nat))
   states(:) = 0
   fist_nonbond_env%charges = neutral_charges

   ! Do Neutral Real and Self
   fR_neut    = 0.0_dp
   ES_pot     = 0.0_dp
   CALL calc_ewald_real(fist_nonbond_env, ewald_env, particle_set, cell,      &
                        ES_pot, fR_neut, ewald_type=ewald_type)

   ! Self energy fully neutral
   thermo%e_self = 0.0_dp
   IF (ewald_type == do_ewald_fssh_ewald) THEN
       neut_charge = 0.0_dp
       CALL calc_ewald_self(ewald_env, cell, atomic_kind_set, thermo%e_self,      &
                            thermo%e_neut, fist_nonbond_env%charges, neut_charge)
   END IF
   ! =====

   ! Print Neutral Charges
   IF (do_print_charges) THEN
      file_num = 179534
      OPEN(unit=file_num, file="charges.list")
      WRITE(file_num, "(A)") "Neutral State Charges:"
      CALL print_charges(neutral_charges, SIZE(mixed_env%F, 3),               &
                         local_molecules, molecule_set, file_num)
   END IF

   ! Do addition subtraction scheme
   DO ist=1, SIZE(mixed_env%F, 3)
      vR_coulomb = 0.0_dp

      charged_charges = neutral_charges
      ! Get start and end of mol
      imol = my_mol_index(ist)
      CALL get_mol_lims(start_atom, end_atom, imol, local_molecules, molecule_set)

      states(start_atom:end_atom) = ist

      ! Set the charges for the index of the molecule as ind istate
      my_mlen = 1 + end_atom - start_atom
      IF (SIZE(charge_state_q) /= my_mlen) &
         CPABORT("Length of charge array /= atoms per mol.")

      charged_charges(start_atom: end_atom) = charge_state_q

      ! Do some printing if necessary
      IF (do_print_charges) THEN
         WRITE(file_num, "(A,I3,A)") "Charge State", ist, ":"
         CALL print_charges(charged_charges, SIZE(mixed_env%F, 3),                  &
                            local_molecules, molecule_set, file_num)
      END IF

      ! Correct the self energy
      self_pot_corr = 0.0_dp
      neut_pot_corr = 0.0_dp
      IF (ewald_type .eq. do_ewald_fssh_ewald) THEN
          CALL calc_ewald_self(ewald_env, cell, atomic_kind_set, self_pot_corr,         &
                               neut_pot_corr, fist_nonbond_env%charges, neut_charge,    &
                               charge_state_q, start_atom, end_atom)
      END IF
               
      mixed_env%F(:, :, ist) = mixed_env%F(:, :, ist) + fR_neut
      my_energies(ist) = my_energies(ist) + neut_pot_corr +   &
                         self_pot_corr + thermo%e_self + ES_pot
   END DO

   ! Set the charged state charges 
   charged_charges(:) = neutral_charges(:)
   DO ist=1, SIZE(mixed_env%F, 3)
      imol = my_mol_index(ist)
      ! Get start and end of mol
      CALL get_mol_lims(start_atom, end_atom, imol, local_molecules, molecule_set)
      charged_charges(start_atom:end_atom) = charge_state_q
   END DO

   ! Correct the real-space energy and forces
   CALL calc_ewald_real_add_sub_corr2(fist_nonbond_env, ewald_env, particle_set, &
                                     cell, charged_charges, states, mixed_env,   &
                                     my_energies, ewald_type=ewald_type)

   ! I haven't decide how to handle the printing of the mixed energies yet
   mixed_env%energies(iforce_eval) = 0

   DEALLOCATE(fR_neut)
   DEALLOCATE(states)
   !--------------------- End Self/Real Space --------------------------------!



   ! ---------------------- Reciprocal Space ---------------------------------!
   IF (ewald_type == do_ewald_fssh_ewald) THEN
       charged_charges(:) = neutral_charges(:)
       ! Set the charged state charges 
       DO ist=1, SIZE(mixed_env%F, 3)
          imol = my_mol_index(ist)
          ! Get start and end of mol
          CALL get_mol_lims(start_atom, end_atom, imol, local_molecules, molecule_set)
          charged_charges(start_atom:end_atom) = charge_state_q
       END DO
    
       ! Reciprocal part of the Ewald Sum (ammends F and my_energies)
       CALL calc_ewald_reciprocal(ewald_env, ewald_pw, cell, atomic_kind_set,      &
                                  particle_set, local_particles, mixed_env,        &
                                  my_energies, neutral_charges, fist_env,          &
                                  my_mol_index )
   END IF
   ! ------------------ End Reciprocal Space --------------------------------!


   ! Only print charges once per force_eval
   ES_pot = 0.0_dp
   IF (do_print_charges) THEN
      CLOSE(file_num)
      CALL ewald_env_set(ewald_env, print_charges=.FALSE.)
   END IF

   
   ! ------------------ Bonded Term Corrections -----------------------------!
   IF (ewald_type .EQ. do_ewald_fssh_ewald) THEN
        ALLOCATE(f_bond_C(3, nat))
        ALLOCATE(f_bond_N(3, nat))
        ALLOCATE(e_bond_N(nat))
        ALLOCATE(e_bond_C(nat))
        e_bond_C= 0.0_dp
        e_bond_N = 0.0_dp
        f_bond_C= 0.0_dp
        f_bond_N = 0.0_dp

        ! Can be slightly optimised by ignoring inactive molecules. However, 
        !  this is only used in Ewald and the reciprocal space sum is so much
        !  more expensive I didn't bother.

        ! Do calc for all charge state
        fist_nonbond_env%charges  = charged_charges
        CALL bonded_correct_gaussian_fssh(fist_nonbond_env, atomic_kind_set,        &
                                          particle_set, ewald_env,                  &
                                          e_bond_C, cell=cell, force_out=f_bond_C )
        
        ! Do calc for all neutral states
        fist_nonbond_env%charges  = neutral_charges
        CALL bonded_correct_gaussian_fssh(fist_nonbond_env, atomic_kind_set,        &
                                          particle_set, ewald_env,                  &
                                          e_bond_N, cell=cell, force_out=f_bond_N )
        sum_e_bond_n = SUM(e_bond_N)
        
        ! Combine the charge and neutral states for each site energy/force
        DO ist=1, SIZE(mixed_env%F, 3)
           imol = my_mol_index(ist)
           CALL get_mol_lims(start_atom, end_atom, imol, local_molecules, molecule_set)

           ! Save Forces
           mixed_env%F(:, :start_atom-1, ist) = mixed_env%F(:, :start_atom-1, ist) &
                                                + f_bond_N(:, :start_atom-1)
           mixed_env%F(:, start_atom: end_atom, ist) = mixed_env%F(:, start_atom: end_atom, ist) &
                                                       + f_bond_C(:, start_atom: end_atom)
           mixed_env%F(:, end_atom+1:, ist) = mixed_env%F(:, end_atom+1:, ist)     &
                                              + f_bond_N(:, end_atom+1:)

           ! Save Energies
           my_energies(ist) = my_energies(ist) +                                   &
                                           (sum_e_bond_n                           &
                                            - SUM(e_bond_N(start_atom: end_atom))  &
                                            + SUM(e_bond_C(start_atom: end_atom))  )
        END DO


        DEALLOCATE(f_bond_C, f_bond_N)
        DEALLOCATE(e_bond_C, e_bond_N)
        DEALLOCATE(neutral_charges)
        DEALLOCATE(charged_charges)
   END IF
   ! ---------------- End Bonded Term Corrections ---------------------------!

   ! This line is only for printing
   mixed_env%energies(iforce_eval) = mixed_env%energies(iforce_eval) + ES_pot
   !--------------------------------------------------------------------------!


   !------------------- Finalise ---------------------------------------------!
   IF(iw>0)THEN
      CALL fist_env_get(fist_env, thermo=thermo)

      WRITE(iw,'(A)')" FIST::(3)TOTAL FORCES - THE END..."
      WRITE(iw,'(3f15.9)')((particle_set(i)%f(j),j=1,3), i = 1, nat)
      WRITE(iw,'(A,f15.9)')"Energy after FIST calculation.. exiting now ::",thermo%pot
   END IF

   mixed_env%sh_env%finished_coulomb_step = .TRUE.
   CALL timestop(handle)
   !--------------------------------------------------------------------------!

END SUBROUTINE add_coulomb_interactions

! *****************************************************************************
!> \brief Will print out the charges grouped by molecule.
!> \param charges -> The charges to print
!> \param nstate -> The number of charge states
!> \param nstate -> local_molecules ...
!> \param nstate -> molecule_set ...
!> \author Matt
!> \note  Only really used for debugging.
! *****************************************************************************
SUBROUTINE print_charges(charges, nstate, local_molecules, molecule_set, &
                         file_num)
   REAL(KIND=dp), DIMENSION(:), INTENT(IN)          :: charges
   INTEGER, INTENT(IN)                              :: nstate
   TYPE(distribution_1d_type), POINTER, INTENT(IN)  :: local_molecules
   TYPE(molecule_type), POINTER, INTENT(IN)         :: molecule_set(:)
   INTEGER, INTENT(IN)                              :: file_num


   CHARACTER(len=*), PARAMETER               :: routineN="print_charges",     &
                                                routineP=moduleN//':'//routineN

   INTEGER                                          :: ist, iat, start_atom,  &
                                                       end_atom

   DO ist=1, nstate
      CALL get_mol_lims(start_atom, end_atom, ist, local_molecules,           &
                        molecule_set)
      DO iat=start_atom, end_atom
         WRITE(file_num, "(F16.6,2X)", ADVANCE="no") charges(iat)
      END DO
      WRITE(file_num, "(A)") " "
   END DO
END SUBROUTINE print_charges

! ****************************************************************************
!> \brief Will get the atomic indices that act as the limits of the molecule
!>        with index istate.
!> \param fist_nonbond_env => The env type that stores the charges array
!> \param start_atom => Will be set to the index of the first atom in the
!>                      molecule.
!> \param end_atom => Will be set to the index of the last atom in the
!>                      molecule.
!> \param istate   => The charge state index.
!> \param local_molecules => The distribution_1d_type that has info on molecules
!> \param molecule_set => An array of molecules
!> \author Matt
! ****************************************************************************
SUBROUTINE get_mol_lims(start_atom, end_atom, istate, local_molecules,       &
                        molecule_set)
   INTEGER, INTENT(OUT)                      :: start_atom, end_atom
   INTEGER, INTENT(IN)                       :: istate
   TYPE(distribution_1d_type), POINTER       :: local_molecules
   TYPE(molecule_type), POINTER              :: molecule_set(:)

   CHARACTER(len=*), PARAMETER               :: routineN="get_mol_lims",     &
                                                routineP=moduleN//':'//routineN

   TYPE(molecule_type), POINTER              :: molecule
   INTEGER                                   :: my_mol_ind

   ! Set the charge state
   my_mol_ind = local_molecules % list(1) % array(istate)
   molecule => molecule_set ( my_mol_ind )

   ! Now ammend the charges of the mol in the charged state
   CALL get_molecule( molecule, first_atom=start_atom, last_atom=end_atom )

END SUBROUTINE get_mol_lims


! **************************************************************************************************
!> \brief Computes the reciprocal part of the ewald sum. Also ammends the
!>        mixed_env%F and the my_energies array
!>      Ref.: J.-P. Hansen, Enrico Fermi School, 1985
!>      Note: Only the positive G-vectors are used in the sum.
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param mixed_env -> contains F array and num states
!> \param my_energies -> energies to ammend
!> \param charges ...
!> \par History
!>      JGH (21-Feb-2001) : changed name
!>      Matt Ellis, 95ellismle@gmail.com, (2-06-2020) : Moved to SH part of code and modified for surface
!>                         hopping.
!> \note The subroutine is copied from the ewald evaluate from the ewals.F file
!> \author Matt Ellis, 95ellismle@gmail.com, (based on SR by CJM)
! **************************************************************************************************
   SUBROUTINE calc_ewald_reciprocal(ewald_env, ewald_pw, cell, atomic_kind_set, particle_set, &
                                    local_particles, mixed_env, my_energies, neutral_q,     &
                                    fist_env, my_mol_index)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: my_energies
      REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(IN)   :: neutral_q
      TYPE(fist_environment_type), POINTER, INTENT(IN)   :: fist_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_mol_index ! the
                                                          ! index of active mols

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ewald_reciprocal', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: S_neut, S_prime
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sk_neut,  &
                                                            Ak,       &
                                                            Sk_prime
      INTEGER :: gpt, group, handle, iparticle, iparticle_kind, iparticle_local, lp, mp, nnodes, &
         node, np, nparticle_kind, nparticle_local, ist, iat, start_atom, end_atom, node1, imol
      INTEGER, DIMENSION(:, :), POINTER                  :: bds
      REAL(KIND=dp)                                      :: alpha, factor,   &
                                                            four_alpha_sq, vg_coulomb, &
                                                            gauss, pref, e_igdotr(3)
      REAL(KIND=dp), DIMENSION(3)                        :: vec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fg_coulomb!,         &
                                                            !fNg_coulomb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho0
      REAL(KIND=dp), DIMENSION(:), POINTER               :: charge_state_q
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: charge_corrections
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: at_to_node_map,     &
                                                            node_to_at_map
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(dg_rho0_type), POINTER                        :: dg_rho0
      TYPE(dg_type), POINTER                             :: dg
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(structure_factor_type)                        :: exp_igr

      CALL timeset(routineN, handle)
      CALL cite_reference(Ewald1921)

      ! pointing
      NULLIFY(charge_state_q, local_molecules, molecule_set)
      CALL fist_env_get(fist_env=fist_env,                     &
                        local_molecules=local_molecules,       &
                        molecule_set=molecule_set )
      CALL ewald_env_get(ewald_env, alpha=alpha, group=group,  &
                         charge_state_q=charge_state_q)
      CALL ewald_pw_get(ewald_pw, pw_big_pool=pw_pool, dg=dg)
      CALL dg_get(dg, dg_rho0=dg_rho0)
      rho0 => dg_rho0%density%pw%cr3d
      pw_grid => pw_pool%pw_grid
      bds => pw_grid%bounds

      ! allocating
      nparticle_kind = SIZE(atomic_kind_set)
      nnodes = 0
      DO iparticle_kind = 1, nparticle_kind
         nnodes = nnodes + local_particles%n_el(iparticle_kind)
      ENDDO

      CALL structure_factor_allocate(pw_grid%bounds, nnodes, exp_igr)

      ALLOCATE (S_neut(1:pw_grid%ngpts_cut))
      ALLOCATE (S_prime(1:pw_grid%ngpts_cut))

      ALLOCATE (node_to_at_map(1:nnodes))
      ALLOCATE (at_to_node_map(SIZE(mixed_env%F, 2)))
      ALLOCATE (charge_corrections(SIZE(mixed_env%F, 2)))

      ! Initializing vg_coulomb and fg_coulomb
      ! defining four_alpha_sq
      four_alpha_sq = 4.0_dp*alpha**2
      ! zero node count
      node = 0

      !------------ Evaluate Structure Factor -----------------!
      !Loop over all charged sites -calculate structure factor for every R.g
      !Loop over sites with different charges
      !Don't add the charge yet
      DO iparticle_kind = 1, nparticle_kind

         ! Number of each types of particle (e.g. num of elements or num of particles with different charges)
         nparticle_local = local_particles%n_el(iparticle_kind)

         ! loop over all particles get charges from charge array and evaluate
         ! structure factor at every grid point for every particle in x, y, z
         DO iparticle_local = 1, nparticle_local
            node = node + 1
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            node_to_at_map(node) = iparticle
            at_to_node_map(iparticle) = node
            CALL matvec_3x3(vec, cell%h_inv, particle_set(iparticle)%r)
            CALL structure_factor_evaluate(vec, exp_igr%lb,     &
                                           exp_igr%ex(:, node), &
                                           exp_igr%ey(:, node), &
                                           exp_igr%ez(:, node))
         END DO
      END DO

      ! Will loop over charged sites and multiply the structure factor by q_i 
      !summme is the structure factor for each grid point, gpt
      S_neut(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      ! looping over the positive g-vectors (grid points)
      ! Multiply every element of the structure factor by the charge on 
      !  the corresponding site.
      DO gpt = 1, pw_grid%ngpts_cut_local

         lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
         mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
         np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))

         ! lp, mp and np are indexers for the grid spacings.
         lp = lp + bds(1, 1)
         mp = mp + bds(1, 2)
         np = np + bds(1, 3)
         
         ! If the denominator is below a threshold then ignore it.
         IF (pw_grid%gsq(gpt) <= 1.0E-10_dp) CYCLE

         ! initializing sum to be used in the energy and force for each node/particle
         ! Add charge to the structure factor...
         DO node = 1, nnodes
            iat = node_to_at_map(node)
            S_neut(gpt) = S_neut(gpt) + neutral_q(iat) *         &
                                        (  exp_igr%ex(lp, node)  &
                                          * exp_igr%ey(mp, node) &
                                          * exp_igr%ez(np, node) )
         END DO
      END DO
      !----------- End Structure Factor ----------------------!

      ! Pre-calculate the charges
      charge_corrections(:) = 0.0_dp
      DO ist=1, SIZE(mixed_env%F, 3)
         imol = my_mol_index(ist)
         CALL get_mol_lims(start_atom, end_atom, imol, local_molecules,         &
                           molecule_set)
         node1 = 0
         DO iat=start_atom, end_atom
            node1 = node1 + 1

            IF (node1 > SIZE(charge_state_q)) &
               CPABORT("Each mol must have same number of atoms")

            charge_corrections(iat) = charge_state_q(node1) - neutral_q(iat)
         END DO
      END DO

      !--------------- Calculate Potential -------------------------!
      ! The following loop will use the addition-subtraction scheme to calculate
      !  only potentials.
      !
      !Define prefactor
      pref = fourpi/pw_grid%vol

      States: DO ist=1, SIZE(mixed_env%F, 3)
         vg_coulomb = 0.0_dp

         ! Addition Subtraction Scheme
         imol = my_mol_index(ist)
         CALL get_mol_lims(start_atom, end_atom, imol, local_molecules,         &
                           molecule_set)

         S_prime(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

         ! looping over the positive g-vectors (grid points)
         K_Loop: DO gpt = 1, pw_grid%ngpts_cut_local
            ! Get grid points
            lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
            mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
            np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))
   
            lp = lp + bds(1, 1)
            mp = mp + bds(1, 2)
            np = np + bds(1, 3)

            ! If the denominator is below a threshold then ignore it.
            IF (pw_grid%gsq(gpt) <= 1.0E-10_dp) CYCLE

            ! Addition-Subtraction only for potential
            DO iat=start_atom, end_atom
               node = at_to_node_map(iat)

               IF (charge_corrections(iat) == 0) CYCLE

               S_prime(gpt) = S_prime(gpt) + charge_corrections(iat) * &
                                              (  exp_igr%ex(lp, node)    &
                                                * exp_igr%ey(mp, node)   &
                                                * exp_igr%ez(np, node) )
            END DO

            ! computing the potential energy
            !Gaussian/(k^2) prefactor
            gauss = (rho0(lp, mp, np) * pw_grid%vol)**2 / pw_grid%gsq(gpt)
            factor = REAL(  (S_neut(gpt) + S_prime(gpt))       &
                           * CONJG(S_neut(gpt) + S_prime(gpt)) &
                          , KIND=dp) * gauss

            ! Reciprocal Potential
            vg_coulomb = vg_coulomb + factor
         END DO K_loop

         vg_coulomb = vg_coulomb * pref
         my_energies(ist) = my_energies(ist) + vg_coulomb
      END DO States
      
      DEALLOCATE(S_prime)
      
      ALLOCATE (Sk_neut(3, 1:pw_grid%ngpts_cut))
      Sk_neut(1, :) = S_neut
      Sk_neut(2, :) = S_neut
      Sk_neut(3, :) = S_neut
      DEALLOCATE(S_neut)
      ALLOCATE(Sk_prime(3, pw_grid%ngpts_cut))
      ALLOCATE(Ak(3, pw_grid%ngpts_cut))
      ALLOCATE(fg_coulomb(3, SIZE(mixed_env%F, 2)))
      
      ! Pre-calculate some useful values
      Ak = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO gpt=1, pw_grid%ngpts_cut_local
         lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
         mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
         np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))
   
         lp = lp + bds(1, 1)
         mp = mp + bds(1, 2)
         np = np + bds(1, 3)

         IF (pw_grid%gsq(gpt) <= 1.0E-10_dp) CYCLE
         
         gauss = (rho0(lp, mp, np)*pw_grid%vol)**2/pw_grid%gsq(gpt)
         Ak(:, gpt) = pw_grid%g(:, gpt) * gauss
         Sk_neut(:, gpt) = Sk_neut(:, gpt) * pw_grid%g(:, gpt) * gauss
      END DO

      !------------------- Force calculation -------------------------------!
      ! Loop over all states and calculate the force.
      States2: DO ist=1, SIZE(mixed_env%F, 3)
         fg_coulomb(:, :) = 0.0_dp

         ! Addition Subtraction Scheme
         imol = my_mol_index(ist)
         CALL get_mol_lims(start_atom, end_atom, imol, local_molecules,         &
                           molecule_set)


         Sk_prime(:, :) = Sk_neut(:, :)

         ! looping over the positive g-vectors (grid points)
         K_Loop2: DO gpt = 1, pw_grid%ngpts_cut_local
            ! Get grid points
            lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
            mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
            np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))
   
            lp = lp + bds(1, 1)
            mp = mp + bds(1, 2)
            np = np + bds(1, 3)

            ! If the denominator is below a threshold then ignore it.
            IF (pw_grid%gsq(gpt) <= 1.0E-10_dp) CYCLE

            ! Get the charge and neutral structure factors
            node1 = 0
            DO iat=start_atom, end_atom
               node1 = node1 + 1
               IF (charge_corrections(iat) == 0) CYCLE

               node = at_to_node_map(iat)
               Sk_prime(:, gpt) = Sk_prime(:, gpt) + Ak(:, gpt) *        &
                                               charge_corrections(iat) * &
                                               (  exp_igr%ex(lp, node)   &
                                                * exp_igr%ey(mp, node)   &
                                                * exp_igr%ez(np, node) )
            END DO

            ! Computing the force
            At_Loop2: DO node=1, nnodes
               ! e_igdotr = SUM_j [ q_j * sin(2 pi k * Rij) ]
               e_igdotr(:) = AIMAG(Sk_prime(:, gpt) *           &
                                    CONJG(  exp_igr%ex(lp, node) &
                                          * exp_igr%ey(mp, node) &
                                          * exp_igr%ez(np, node)))

               ! Neutral Reciprocal Force
               fg_coulomb(:, node) = fg_coulomb(:, node) + e_igdotr(:)
            END DO At_Loop2
         END DO K_loop2

         ! Add the outer loop charge now
         DO iat=1, start_atom-1
            node = at_to_node_map(iat)
            fg_coulomb(:, node) = neutral_q(iat) * fg_coulomb(:, node)
         END DO
         node1 = 1
         DO iat=start_atom, end_atom
            node = at_to_node_map(iat)
            fg_coulomb(:, node) = charge_state_q(node1) * fg_coulomb(:, node)
            node1 = node1 + 1
         END DO
         DO iat=end_atom+1, SIZE(mixed_env%F, 2)
            node = at_to_node_map(iat)
            fg_coulomb(:, node) = neutral_q(iat) * fg_coulomb(:, node)
         END DO

         fg_coulomb = fg_coulomb * (2.0_dp * pref)

         ! Pass the forces back to the mixed_env type
         DO node=1, nnodes
            iat = node_to_at_map(node)
            mixed_env%F(:, iat, ist) = mixed_env%F(:, iat, ist) + fg_coulomb(:, node)
         END DO

      END DO States2
      !-------------------- End Force Calculation -------------------------!


      DEALLOCATE(node_to_at_map, at_to_node_map, &
                 charge_corrections, Sk_neut, Sk_prime)
      CALL structure_factor_deallocate(exp_igr)

      CALL timestop(handle)

   END SUBROUTINE calc_ewald_reciprocal

! **************************************************************************************************
!> \brief computes the real part of the coulomb potential (for the surface hopping site 
!>                                                         energies/forces)
!> \param fist_nonbond_env ...
!> \param ewald_env ...
!> \param particle_set ...
!> \param cell ...
!> \param pot_coul <REAL(kind=dp)> The real coulomb potential output
!> \param frc_coul <REAL(kind=dp)(3, nat)> The real coulomb force output
!> \param atomic_kind_set ...    
!> \author Matt
! **************************************************************************************************
SUBROUTINE calc_ewald_real(fist_nonbond_env, ewald_env, particle_set, cell,     &
                           pot_coul, frc_coul, ewald_type)
      
    TYPE(fist_nonbond_env_type), POINTER                  :: fist_nonbond_env
    TYPE(ewald_environment_type), POINTER                 :: ewald_env
    TYPE(particle_type), DIMENSION(:), INTENT(IN)         :: particle_set
    TYPE(cell_type), POINTER                              :: cell
    REAL(KIND=dp), INTENT(OUT)                            :: pot_coul 
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)         :: frc_coul
    INTEGER                                               :: ewald_type
    
    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ewald_real', &
      routineP = moduleN//':'//routineN
     
    TYPE(cp_logger_type), POINTER                         :: logger 
    TYPE(fist_neighbor_type), POINTER                     :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER               :: neighbor_kind_pair   
    TYPE(pair_potential_pp_type), POINTER                 :: potparm, potparm14
    TYPE(pair_potential_single_type), POINTER             :: pot
    TYPE(pos_type), DIMENSION(:), POINTER                 :: r_last_update_pbc

    INTEGER                                               :: atom_i, atom_j,      &    
                                                             handle, iend, igrp,  &
                                                             ilist, ipair,        &
                                                             istart, kind_a,      &
                                                             npairs, nkind,       &
                                                             kind_b, count_pairs, &
                                                             max_at, min_at
    INTEGER, DIMENSION(:, :), POINTER                     :: list
    LOGICAL                                               :: full_nl
    REAL(KIND=dp)                                         :: alpha, beta,         &
                                                             fac_ei, fac_kind,    &
                                                             fscalar, rab2_max,   &
                                                             qeff_a, qeff_b, rab2,&
                                                             rab2_com, rcut
    REAL(KIND=dp), DIMENSION(3)                           :: cell_v, cvi, rab_com,&
                                                             rab, fr, rab_cc
    REAL(KIND=dp), DIMENSION(:, :), POINTER               :: ij_kind_full_fac,    &
                                                             ei_interaction_cutoffs_dsf
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER            :: ei_interaction_cutoffs
     
    CALL timeset(routineN,handle)
    NULLIFY (logger)
    logger => cp_get_default_logger()
     
    NULLIFY(pot, ij_kind_full_fac)
    CALL fist_nonbond_env_get(fist_nonbond_env, nonbonded=nonbonded,           &
                              potparm14=potparm14, potparm=potparm,            &
                              r_last_update_pbc=r_last_update_pbc,             &
                              natom_types=nkind,                               &
                              ij_kind_full_fac=ij_kind_full_fac                )
    
    CALL ewald_env_get(ewald_env, alpha=alpha,                                 &
                       interaction_cutoffs=ei_interaction_cutoffs,             &
                       interaction_cutoffs_dsf=ei_interaction_cutoffs_dsf,     &
                       rcut=rcut)
    
    frc_coul = 0.0_dp
    pot_coul = 0.0_dp
    qeff_a   = 0.0_dp
    qeff_b   = 0.0_dp
       
    CPASSERT(ASSOCIATED(fist_nonbond_env%charges))


    ! Starting the force loop
    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs = neighbor_kind_pair%npairs
       count_pairs = 0
      
       neighbor_kind_pair%do_list = .TRUE.
       IF (npairs == 0) THEN
         neighbor_kind_pair%do_list = .FALSE.
         CYCLE
       END IF

       ! Allocate some arrays for optimisation
       IF (ALLOCATED(neighbor_kind_pair%skip_grp)) &
         DEALLOCATE(neighbor_kind_pair%skip_grp)
       ALLOCATE(neighbor_kind_pair%skip_grp(neighbor_kind_pair%ngrp_kind))

       IF (ALLOCATED(neighbor_kind_pair%max_ats)) &
         DEALLOCATE(neighbor_kind_pair%max_ats)
       ALLOCATE(neighbor_kind_pair%max_ats(neighbor_kind_pair%ngrp_kind))

       IF (ALLOCATED(neighbor_kind_pair%min_ats)) &
         DEALLOCATE(neighbor_kind_pair%min_ats)
       ALLOCATE(neighbor_kind_pair%min_ats(neighbor_kind_pair%ngrp_kind))

       neighbor_kind_pair%max_at = 0
       neighbor_kind_pair%min_at = 99999999

       list  => neighbor_kind_pair%list
       cvi   =  neighbor_kind_pair%cell_vector
       CALL matvec_3x3(cell_v,cell%hmat, cvi)
     
       Kind_Group_Loop: DO igrp=1,neighbor_kind_pair%ngrp_kind
          istart = neighbor_kind_pair%grp_kind_start(igrp)
          iend   = neighbor_kind_pair%grp_kind_end(igrp)
      
          neighbor_kind_pair%skip_grp(igrp) = .TRUE.

          ! Save max and min atoms for later add sub scheme
          max_at = 0
          min_at = 99999999

          Pairs: DO ipair=istart,iend
             atom_i = list(1,ipair)
             atom_j = list(2,ipair)

             ! Get the electrostatic parameters for the atoms a and b
             qeff_a = fist_nonbond_env%charges(atom_i)
             qeff_b = fist_nonbond_env%charges(atom_j)
   
             IF (qeff_a == 0 .OR. qeff_b == 0) THEN
               CYCLE
             END IF
     
             ! Get actual atomic kinds, since atom_i is not always of
             ! kind_a and atom_j of kind_b, ie. they might be swapped.
             kind_a = particle_set(atom_i)%atomic_kind%kind_number
             kind_b = particle_set(atom_j)%atomic_kind%kind_number
     
             ! Get scaling factor for kind pair
             fac_kind = ij_kind_full_fac(kind_a,kind_b)
     
             ! take the proper potential
             pot => potparm%pot(kind_a,kind_b)%pot
     
             ! Determine the scaling factors
             fac_ei = fac_kind
             full_nl = ANY(pot%type == tersoff_type).OR.ANY(pot%type == siepmann_type)
             IF ((.NOT.full_nl).AND.(atom_i == atom_j)) THEN
                fac_ei = 0.5_dp*fac_ei
             END IF 
     
             ! Exclude any atoms that are bonded to each other.
             IF (ipair <= neighbor_kind_pair%nscale) &
                fac_ei = fac_ei*neighbor_kind_pair%ei_scale(ipair)
             IF (fac_ei <= 0.0) THEN
               !list(1, ipair) = -999
               CYCLE
             END IF
     
             ! Set beta to 0 (an input to potential_coulomb)
             beta = 0.0_dp
     
             ! not do shell
             rab_cc  = r_last_update_pbc(atom_j)%r - r_last_update_pbc(atom_i)%r
             rab_com = rab_cc
             rab_com  = rab_com + cell_v
             rab2_com = rab_com(1)**2 + rab_com(2)**2 + rab_com(3)**2
     
             ! compute the interactions for the current pair
             rab2_max = pot%rcutsq
             IF (rab2_com <= rab2_max) THEN
                ! Ion-Ion: no shell model coulomb
                rab  = rab_com
                rab2 = rab2_com
                fscalar = 0.0_dp
                ! note that potential_coulomb increments fscalar
                pot_coul = pot_coul + &
                             potential_coulomb(rab2, fscalar, fac_ei*qeff_a*qeff_b,       &
                                               ewald_type, alpha, beta,                   &
                                               ei_interaction_cutoffs(3, kind_a, kind_b), &
                                               dsf_cutoff=ei_interaction_cutoffs_dsf(kind_a,kind_b),&
                                               rc=rcut)
                
                ! Forces
                fr = fscalar*rab
                frc_coul(:, atom_i) = frc_coul(:, atom_i) - fr(:)
                frc_coul(:, atom_j) = frc_coul(:, atom_j) + fr(:)
                
                count_pairs = count_pairs + 1

                max_at = MAX(atom_i, max_at)
                max_at = MAX(atom_j, max_at)
                min_at = MIN(atom_i, min_at)
                min_at = MIN(atom_j, min_at)
             ELSE
               !list(1, ipair) = -999
             END IF

             neighbor_kind_pair%skip_grp(igrp) = .FALSE.
          END DO Pairs

          neighbor_kind_pair%min_ats(igrp) = min_at
          neighbor_kind_pair%max_ats(igrp) = max_at

          neighbor_kind_pair%max_at = MAX(neighbor_kind_pair%max_ats(igrp), &
                                          neighbor_kind_pair%max_at)
          neighbor_kind_pair%min_at = MIN(neighbor_kind_pair%min_ats(igrp), &
                                          neighbor_kind_pair%min_at)

       END DO Kind_Group_Loop


       IF (count_pairs == 0) THEN
          neighbor_kind_pair%do_list = .FALSE.
       END IF

    END DO Lists

    !sample peak memory
    CALL m_memory()

    CALL timestop(handle)

END SUBROUTINE calc_ewald_real

! **************************************************************************************************
!> \brief computes the real part of the addition-subtraction correction to the
!>        potential and forces.
!> \param fist_nonbond_env ...
!> \param ewald_env ...
!> \param particle_set ...
!> \param cell ...
!> \param pot_coul <REAL(kind=dp)> The real coulomb potential output
!> \param frc_coul <REAL(kind=dp)(3, nat)> The real coulomb force output
!> \author Matt
! **************************************************************************************************
SUBROUTINE calc_ewald_real_add_sub_corr2(fist_nonbond_env, ewald_env, particle_set, cell,     &
                                         charged_q, states, mixed_env, my_energies, ewald_type)
      
    TYPE(fist_nonbond_env_type), POINTER, INTENT(IN)   :: fist_nonbond_env
    TYPE(ewald_environment_type), POINTER, INTENT(IN)  :: ewald_env
    TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
    TYPE(cell_type), POINTER, INTENT(IN)               :: cell
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charged_q
    INTEGER, DIMENSION(:)                              :: states
    TYPE(mixed_environment_type), POINTER              :: mixed_env
    REAL(KIND=dp), DIMENSION(:), POINTER               :: my_energies
    INTEGER                                            :: ewald_type
    
    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ewald_real_add_sub_corr2', &
      routineP = moduleN//':'//routineN
     
    TYPE(cp_logger_type), POINTER                      :: logger 
    TYPE(fist_neighbor_type), POINTER                  :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair   
    TYPE(pair_potential_pp_type), POINTER              :: potparm, potparm14
    TYPE(pair_potential_single_type), POINTER          :: pot
    TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc

    INTEGER                                            :: atom_i, atom_j,      &    
                                                          handle, iend, igrp,  &
                                                          ilist, ipair,        &
                                                          istart, kind_a,      &
                                                          nkind, kind_b,       &
                                                          max_at, min_at,      &
                                                          ist_i, ist_j
    INTEGER, DIMENSION(:, :), POINTER                  :: list
    REAL(KIND=dp)                                      :: alpha, beta,         &
                                                          fac_ei, fac_kind,    &
                                                          fscalar, rab2,       &
                                                          rab2_com, q_cn, q_nn,&
                                                          q_cc, qic, qjc, qin, &
                                                          qjn, pot_coul,       &
                                                          qcc_nn, qcn_nn, rcut
    REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi, rab_com,&
                                                          rab, fr
    REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ij_kind_full_fac,    &
                                                          ei_interaction_cutoffs_dsf
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: ei_interaction_cutoffs



    CALL timeset(routineN,handle)
    NULLIFY (logger)
    logger => cp_get_default_logger()
     
    NULLIFY(pot, ij_kind_full_fac)
    CALL fist_nonbond_env_get(fist_nonbond_env, nonbonded=nonbonded,           &
                              potparm14=potparm14, potparm=potparm,            &
                              r_last_update_pbc=r_last_update_pbc,             &
                              natom_types=nkind,                               &
                              ij_kind_full_fac=ij_kind_full_fac                )
    
    CALL ewald_env_get(ewald_env, alpha=alpha,                                 &
                       interaction_cutoffs=ei_interaction_cutoffs,             &
                       interaction_cutoffs_dsf=ei_interaction_cutoffs_dsf,     &
                       rcut=rcut)
    
    CPASSERT(ASSOCIATED(fist_nonbond_env%charges))


    ! Starting the force loop
    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       IF (neighbor_kind_pair%do_list .eqv. .FALSE.) CYCLE

       list  => neighbor_kind_pair%list
       cvi   =  neighbor_kind_pair%cell_vector
       CALL matvec_3x3(cell_v,cell%hmat, cvi)
     
       Kind_Group_Loop: DO igrp=1,neighbor_kind_pair%ngrp_kind
          IF (neighbor_kind_pair%skip_grp(igrp)) CYCLE

          min_at = neighbor_kind_pair%min_ats(igrp)
          max_at = neighbor_kind_pair%max_ats(igrp)

          istart = neighbor_kind_pair%grp_kind_start(igrp)
          iend   = neighbor_kind_pair%grp_kind_end(igrp)

          Pairs: DO ipair=istart,iend
             atom_i = list(1,ipair)
             atom_j = list(2,ipair)
             ist_i = states(atom_i)
             ist_j = states(atom_j)
             IF ((ist_i == 0) .AND. (ist_j == 0)) CONTINUE

             ! Get actual atomic kinds, since atom_i is not always of
             ! kind_a and atom_j of kind_b, ie. they might be swapped.
             kind_a = particle_set(atom_i)%atomic_kind%kind_number
             kind_b = particle_set(atom_j)%atomic_kind%kind_number
     
             ! Get scaling factor for kind pair
             fac_kind = ij_kind_full_fac(kind_a,kind_b)
     
             ! take the proper potential
             pot => potparm%pot(kind_a,kind_b)%pot
     
             ! Determine the scaling factors
             fac_ei = fac_kind
             IF (atom_i == atom_j) THEN
                fac_ei = 0.5_dp*fac_ei
             END IF 
     
             ! Exclude any atoms that are bonded to each other.
             IF (ipair <= neighbor_kind_pair%nscale) &
                fac_ei = fac_ei*neighbor_kind_pair%ei_scale(ipair)
             IF (fac_ei <= 0.0) CYCLE

             ! Set beta to 0 (an input to potential_coulomb)
             beta = 0.0_dp
     
             ! not do shell
             rab_com  = r_last_update_pbc(atom_j)%r - r_last_update_pbc(atom_i)%r + cell_v
             rab2_com = rab_com(1)**2 + rab_com(2)**2 + rab_com(3)**2
     
             ! compute the interactions for the current pair
             IF (rab2_com <= pot%rcutsq) THEN
                ! Ion-Ion: no shell model coulomb
                rab  = rab_com
                rab2 = rab2_com
                fscalar = 0.0_dp
             
                ! Get the charges for atoms i and j
                qin = fist_nonbond_env%charges(atom_i)
                qjn = fist_nonbond_env%charges(atom_j)
                qjc = charged_q(atom_j)
                IF (ist_i == ist_j) THEN 
                     qic = charged_q(atom_i)
                ELSE
                     qic = fist_nonbond_env%charges(atom_i)
                END IF
                qcc_nn = fac_ei*((qic*qjc) - (qin*qjn))
               
                IF (ABS(qcc_nn) > 1e-10 .AND. ist_j /= 0) THEN
                   pot_coul = potential_coulomb(rab2, fscalar, qcc_nn,                     &
                                                ewald_type, alpha, beta,                   &
                                                ei_interaction_cutoffs(3, kind_a, kind_b), &
                                                dsf_cutoff=ei_interaction_cutoffs_dsf(kind_a,kind_b),&
                                                rc=rcut)

                   fr = fscalar*rab
                   my_energies(ist_j) = my_energies(ist_j) + pot_coul
                   mixed_env%F(:, atom_i, ist_j) = mixed_env%F(:, atom_i, ist_j) - fr(:)
                   mixed_env%F(:, atom_j, ist_j) = mixed_env%F(:, atom_j, ist_j) + fr(:)
                END IF

                fscalar = 0.0_dp
                IF (ist_i /= ist_j .AND. ist_i /= 0) THEN
                   qic = charged_q(atom_i)
                   IF (qic == qin) CYCLE 
                   qcn_nn = fac_ei * qjn * (qic - qin)

                   pot_coul = potential_coulomb(rab2, fscalar, qcn_nn,                     &
                                                ewald_type, alpha, beta,                   &
                                                ei_interaction_cutoffs(3, kind_a, kind_b), &
                                                dsf_cutoff=ei_interaction_cutoffs_dsf(kind_a,kind_b),&
                                                rc=rcut)

                   fr = fscalar*rab
                   my_energies(ist_i) = my_energies(ist_i) + pot_coul
                   mixed_env%F(:, atom_i, ist_i) = mixed_env%F(:, atom_i, ist_i) - fr(:)
                   mixed_env%F(:, atom_j, ist_i) = mixed_env%F(:, atom_j, ist_i) + fr(:)
                END IF
             END IF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists

    !sample peak memory
    CALL m_memory()

    CALL timestop(handle)

END SUBROUTINE calc_ewald_real_add_sub_corr2

! **************************************************************************************************
!> \brief Computes the self interaction from g-space and the neutralizing background. 
!> \param ewald_env ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param e_self ...
!> \param e_neut ...
!> \param neut_charges ...
!> \author Matt Ellis, 95ellismle@gmail.com, (based on subroutine by CJM)
! **************************************************************************************************
   SUBROUTINE calc_ewald_self(ewald_env, cell, atomic_kind_set, e_self, e_neut, neut_charges, &
                              neut_charge, chrg_charges, start_atom, end_atom)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      REAL(KIND=dp), INTENT(OUT)                         :: e_self, e_neut
      REAL(KIND=dp), DIMENSION(:), POINTER               :: neut_charges
      REAL(KIND=dp), INTENT(INOUT)                       :: neut_charge
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: chrg_charges
      INTEGER, INTENT(IN), OPTIONAL                      :: start_atom, end_atom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ewald_self', routineP = moduleN//':'//routineN

      INTEGER                                            :: ewald_type, group, iparticle_kind, &
                                                            nparticle_kind, iat, node
      REAL(KIND=dp)                                      :: alpha, mm_radius, q_self, &
                                                            q_sum
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      CALL ewald_env_get(ewald_env, ewald_type=ewald_type, &
                         alpha=alpha, group=group)
      nparticle_kind = SIZE(atomic_kind_set)

      ! For add-sub charges
      IF (PRESENT(chrg_charges)) THEN  
         node = 1
         q_self = 0.0_dp
         q_sum = neut_charge
         DO iat=start_atom, end_atom
            q_self = q_self + (chrg_charges(node)**2) - (neut_charges(iat)**2)
            q_sum = q_sum + chrg_charges(node) - neut_charges(iat)
            node = node + 1
         END DO

         IF (node /= (2 + end_atom - start_atom)) &
            CPABORT("All mols must have the same # of atoms")

      ! For neutral systems
      ELSE
         q_self = DOT_PRODUCT(neut_charges, neut_charges)
         q_sum = SUM(neut_charges)
         neut_charge = q_sum

         ! check and abort..
         DO iparticle_kind = 1, nparticle_kind
            atomic_kind => atomic_kind_set(iparticle_kind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, mm_radius=mm_radius)
            IF (mm_radius > 0.0_dp) THEN
               CPABORT("Array of neut_charges not implemented for mm_radius>0.0 !!")
            END IF
         END DO
      END IF

      e_self = -q_self * alpha * oorootpi
      e_neut = -q_sum *pi/(2.0_dp*cell%deth)*(q_sum/alpha**2)

   END SUBROUTINE calc_ewald_self

! *****************************************************************************
! THE multi_block_energies_forces is more general so this subroutine is
! DEPRECATED. It is kept below for testing purposes if necessary!
!> \brief This subroutine does the addition and subtraction scheme for the CT
!         states  
!> \note This subroutine gives non-degenerate site_energies (degeneracy may be
!        introduce in the fssh specific subroutines)
!> \par History
!> \author SG 
! *****************************************************************************
SUBROUTINE ct_block_energies_forces(adiab_section,energies_1_mol_list, energies, &
                                    intra_forces_global, global_forces, my_energies, &
                                    active_system_info, mixed_env, my_itimes,my_mol_index)

   TYPE(section_vals_type), POINTER          :: adiab_section
   TYPE(mixed_environment_type), &
      POINTER                                :: mixed_env
   TYPE(mixed_force_type), DIMENSION(:), &
     POINTER                                 :: global_forces
   REAL(KIND=dp), DIMENSION(:), POINTER      :: energies
   REAL(KIND=dp), DIMENSION(:), POINTER      :: my_energies
   REAL(KIND=dp), ALLOCATABLE                :: intra_forces_global(:,:,:) 
   REAL(KIND=dp), ALLOCATABLE                :: energies_1_mol_list(:,:)
   INTEGER, POINTER                          :: my_itimes !CHANGE_AC
   INTEGER, DIMENSION(:,:), INTENT(IN)       :: active_system_info
   INTEGER, DIMENSION(:), INTENT(IN)         :: my_mol_index !CHANGE SAM

   CHARACTER(len=*), PARAMETER :: routineN ='add_energies_forces_mem_saver_V2', &
     routineP = moduleN//':'//routineN

   INTEGER :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
   INTEGER :: init_atom_mol_j, final_atom_mol_j, atoms_per_mol_j
   INTEGER                                   :: i,j, ind1, ind2!, indfirst,indsec
   INTEGER                                   :: start_i,end_i,start_j,end_j,atoms
   TYPE(section_vals_type), POINTER          :: site_energy_section
   TYPE(section_vals_type), POINTER          :: ct_block
   INTEGER, DIMENSION(:), POINTER            :: donor_index, acceptor_index
   INTEGER                                   :: total_states, state,nadiab, natom
   INTEGER, DIMENSION(:,:), ALLOCATABLE      :: intra_forces_info
   INTEGER, DIMENSION(:), ALLOCATABLE        :: first_ind, second_ind
 
   NULLIFY(site_energy_section,ct_block)
   NULLIFY( my_energies)
 
   site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
   ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
   CALL section_vals_val_get(ct_block,"INDEX_DONOR_MOLS",&
                        i_vals=donor_index)
   CALL section_vals_val_get(ct_block,"INDEX_ACCEPTOR_MOLS",&
                        i_vals=acceptor_index)
   CALL section_vals_val_get(adiab_section,"NUMBER_DIABATIC_STATES",&
                        i_val=nadiab)

   !> some checks 
   total_states = SIZE(donor_index)*SIZE(acceptor_index)
   IF (total_states .ne. nadiab) THEN 
      CPABORT("NUMBER_DIABATIC_STATES different from total n states")
   END IF
   ind1 = donor_index(1) 
   do i=2,SIZE(donor_index)
      if (donor_index(i) .le. ind1) THEN 
        CPABORT("Donor mol must be ascending")
      end if
   end do
   ind1 = acceptor_index(1) 
   do i=2,SIZE(acceptor_index)
      if (acceptor_index(i) .le. ind1) THEN 
        CPABORT("Acceptor mol must be ascending")
      end if
   end do
 
   !> this are the output add-sub energies and forces
   ALLOCATE(my_energies(nadiab))
   natom = SIZE(global_forces(1)%forces(:,:), 2)

   !> we do F allocation just at the first step and never again just to NOT
   !allocate and deallocate this big array at every step. NOTE that this array is
   !now never deallocated and just deallocated by the compiler at the very end.
   !One should probably write a subroutine that deallocate all allocatable just
   !this at the very end of cp2k
    if(my_itimes .eq. 0) then
       ALLOCATE(mixed_env%F(3,natom,nadiab))
    end if
  
   !> construct a intra_forces_info array with the index of the atoms that are
   !part of the intra_forces_global array that stores only the active atoms.
   !intra_forces_info is therefore an array with the same shape as
   !active_system_info but with sequential numbering of the atoms (starting from
   !1,2,3...).
   ALLOCATE(intra_forces_info(2,SIZE(active_system_info,2)))
   atoms = 0
   do i=1,SIZE(active_system_info,2)
         !> find numer of atoms in mol i
         init_atom_mol_i = active_system_info(1,i)
         final_atom_mol_i = active_system_info(2,i)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
         atoms = atoms + atoms_per_mol_i
         intra_forces_info(1,i) = (atoms - atoms_per_mol_i)+1
         intra_forces_info(2,i) = atoms
   end do

!print *, "intra_forces_info(1,i)", intra_forces_info(1,:)   
!print *, "intra_forces_info(2,i)", intra_forces_info(2,:)   

   !get the CT index in the 1,2,3... base (active region basis)
   !independently on the actual molecule
   !numbering
   ! FINDLOC function not supported so loops instead
   ALLOCATE(second_ind(SIZE(acceptor_index)))
   do i=1,SIZE(acceptor_index)
     do j=1,SIZE(my_mol_index) 
         if (acceptor_index(i) .eq. my_mol_index(j)) then
            second_ind(i) = j 
            EXIT
         end if 
     end do 
   end do
   
   ALLOCATE(first_ind(SIZE(donor_index)))
   do i=1,SIZE(donor_index)
     do j=1,SIZE(my_mol_index) 
         if (donor_index(i) .eq. my_mol_index(j)) then
            first_ind(i) = j 
            EXIT
         end if 
     end do 
   end do
   
!   print *, "DECOMP", my_mol_index
!   print *, "DONOR IND", donor_index
!   print *, "ACCEPTOR IND", acceptor_index
!   print *, "FIRST IND", first_ind
!   print *, "SEC IND", second_ind

   !> loop over states: not that it does not matter if donor molecules comes
   !after acceptor in the xyz or viceversa, or the molecules are DDDAAA or ADADA
   !as the index of each molecule is referred to DECOMP file
   state = 0
   do i=1,SIZE(first_ind)
      ind1 = first_ind(i)
      do j=1,SIZE(second_ind)
         ind2 = second_ind(j)
         state = state + 1
!print *, "STATES", ind1,ind2
         !> we create site energies from the 2 force_evals sections: the first
         !neutral, the second charged (with donor and acc Froce-fields).
         my_energies(state) = energies(1) + &
                  energies_1_mol_list(2,ind1) - energies_1_mol_list(1,ind1) + & 
                  energies_1_mol_list(2,ind2) - energies_1_mol_list(1,ind2)  
!print *, "ENER:", my_energies(state)

         !> we allocate diagonal forces array. All the atoms for a give state i
         !get the value of the neutral force_field that contains: intra +
         !nonbonded interactions
         mixed_env%F(:,:,state) = global_forces(1)%forces
         
         !> find index mol i (which is state i) of the pair
         init_atom_mol_i = active_system_info(1,ind1)
         final_atom_mol_i = active_system_info(2,ind1)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

         !> find index mol j (which is state j) of the pair
         init_atom_mol_j = active_system_info(1,ind2)
         final_atom_mol_j = active_system_info(2,ind2)
         atoms_per_mol_j = (final_atom_mol_j-init_atom_mol_j)+1

         !> since  mixed_env%F stores all the atoms the numbering is a
         !bit more triky and relies on the "global numbering" (which is an info
         !contained in active_system_info). Whereas intra_forces_global stores
         !just the active atoms and so the numbering should be sequential and
         !this info is included in intra_forces_info
        start_i = intra_forces_info(1,ind1)
        end_i = intra_forces_info(2,ind1)
!print *, "index don", init_atom_mol_i, final_atom_mol_i
!print *, "NEW", start_i, end_i
        !> first molecule i of the pair
        mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,1,start_i:end_i) - &
                           intra_forces_global(1,1,start_i:end_i)
        mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,2,start_i:end_i) - &
                           intra_forces_global(1,2,start_i:end_i)
        mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) = &
                           mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,state) + &
                           intra_forces_global(2,3,start_i:end_i) - &
                           intra_forces_global(1,3,start_i:end_i)

        start_j = intra_forces_info(1,ind2)
        end_j = intra_forces_info(2,ind2)
!print*, "index acc", init_atom_mol_j, final_atom_mol_j
!print *, "NEW", start_j, end_j
        !> second mol j of the pair
        mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(1,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,1,start_j:end_j) - &
                           intra_forces_global(1,1,start_j:end_j)
        mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(2,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,2,start_j:end_j) - &
                           intra_forces_global(1,2,start_j:end_j)
        mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) = &
                           mixed_env%F(3,init_atom_mol_j:final_atom_mol_j,state) + &
                           intra_forces_global(2,3,start_j:end_j) - &
                           intra_forces_global(1,3,start_j:end_j)

      end do
   end do

   DEALLOCATE(intra_forces_info)
   DEALLOCATE(first_ind)
   DEALLOCATE(second_ind)

END SUBROUTINE ct_block_energies_forces

! *****************************************************************************
!> \brief We do addition and substraction scheme for the site energy on order to
!         save a lot of memory, as the code would read fewer FORCE_EVALS SECTION from
!         input. To make sure that the code is doing the correct thing all the
!         results can be bechmarked agaist the same system but passing all the
!         FORCE_EVALS by input. 
!> \note This subroutine assumes that each molecule correspond to one state
!        (one-particle picture) and that the first FORCE_EVAL is the neutral one for
!         which intra+LJ are calculated, while the second is the charged one.   
!> \note This whole add-sub scheme makes sense just becuase electrostatic is not
!        considered 
!> \note This subroutine gives non-degenerate site_energies (degeneracy may be
!        introduce in the fssh specific subroutines)
!> \par History
!> \author SG 
! *****************************************************************************
 SUBROUTINE add_energies_forces_mem_saver_V2( energies_1_mol_list, energies, intra_forces_global, &
                                              global_forces, my_energies, & 
                                               active_system_info, & 
                                               mixed_env, my_itimes)

   TYPE(mixed_environment_type),  &
                        POINTER             :: mixed_env
   TYPE(mixed_force_type), DIMENSION(:), &
     POINTER                                :: global_forces
   REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
   REAL(KIND=dp), DIMENSION(:), POINTER     :: my_energies
   REAL(KIND=dp), ALLOCATABLE                 :: intra_forces_global(:,:,:) 
   REAL(KIND=dp), ALLOCATABLE               :: energies_1_mol_list(:,:)
   INTEGER, POINTER                         :: my_itimes !CHANGE_AC
   INTEGER, DIMENSION(:,:), INTENT(IN)      :: active_system_info

   CHARACTER(len=*), PARAMETER :: routineN ='add_energies_forces_mem_saver_V2', &
     routineP = moduleN//':'//routineN

   INTEGER :: i, natom, init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
   INTEGER, DIMENSION(:,:), ALLOCATABLE      :: intra_forces_info
   INTEGER                                   :: start_i,end_i,atoms

   NULLIFY( my_energies)
   !> WARNING: the following allocation is fine as the number of states is the
   !same as the number of active molecules in this case, for which the
   !decomposition was done.
   ALLOCATE(my_energies(SIZE(energies_1_mol_list(1,:))))
   !> WARNING: here we assume that all force evals have the same number of atoms,
   !this is always the case, I do not see a reason why it should be otherwise
   natom = SIZE(global_forces(1)%forces(:,:), 2)

   !> we do F allocation just at the first step and never again just to NOT
   !allocate and deallocate this big array at every step. NOTE that this array is
   !now never deallocated and just deallocated by the compiler at the very end.
   !One should probably write a subroutine that deallocate all allocatable just
   !this at the very end of cp2k
   if(my_itimes .eq. 0) then
      ALLOCATE(mixed_env%F(3,natom, SIZE(energies_1_mol_list(1,:))))
   end if

   !> construct a intra_forces_info array with the index of the atoms that are
   !part of the intra_forces_global array that stores only the active atoms.
   !intra_forces_info is therefore an array with the same shape as
   !active_system_info but with sequential numbering of the atoms.
   ALLOCATE(intra_forces_info(2,SIZE(active_system_info,2)))
   atoms = 0
   do i=1,SIZE(active_system_info,2)
         !> find numer of atoms in mol i
         init_atom_mol_i = active_system_info(1,i)
         final_atom_mol_i = active_system_info(2,i)
         atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1
         atoms = atoms + atoms_per_mol_i
         intra_forces_info(1,i) = (atoms - atoms_per_mol_i)+1
         intra_forces_info(2,i) = atoms
   end do

   !> loop over each (active) state of the system, namely active mol 
   DO i = 1, SIZE(energies_1_mol_list(1,:)) !loop over num of adiabatic states 
      !> we create site energies from the 2 force_evals sections: the first
      !neutral, the second charged.
      my_energies(i) = energies(1) + energies_1_mol_list(2,i) - energies_1_mol_list(1,i)

      !> we allocate diagonal forces array. All the atoms for a give state i
      !assume the value of the neutral force_field that contains: intra +
      !nonbonded interactions
      mixed_env%F(:,:,i) = global_forces(1)%forces
      
      !> find index mol i (which is state i) of the pair
      init_atom_mol_i = active_system_info(1,i)
      final_atom_mol_i = active_system_info(2,i)
      atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1


      !> since intra_forces_global stores only active atoms the numbering is a
      !bit more triky and relies on the n of atoms per mol (which is an info
      !contained in active_system_info)
      start_i = intra_forces_info(1,i)
      end_i = intra_forces_info(2,i)
      mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) = &
                         mixed_env%F(1,init_atom_mol_i:final_atom_mol_i,i) + &
                         intra_forces_global(2,1,start_i:end_i) - &
                         intra_forces_global(1,1,start_i:end_i)
      mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) = &
                         mixed_env%F(2,init_atom_mol_i:final_atom_mol_i,i) + &
                         intra_forces_global(2,2,start_i:end_i) - &
                         intra_forces_global(1,2,start_i:end_i)
      mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) = &
                         mixed_env%F(3,init_atom_mol_i:final_atom_mol_i,i) + &
                         intra_forces_global(2,3,start_i:end_i) - &
                         intra_forces_global(1,3,start_i:end_i)
   END DO
  DEALLOCATE(intra_forces_info)
  END SUBROUTINE add_energies_forces_mem_saver_V2

! *****************************************************************************
!> \brief  OLD SUBROUTINE, rewritten in add_energies_forces_mem_saver_V2 by
!using some info arrays with the system structure info
!> \par History
!>      SG 
! *****************************************************************************
 SUBROUTINE add_energies_forces_mem_saver( energies_1_mol_list, energies, intra_forces_global, &
                                              my_mol_index, &
                                              global_forces, my_energies, & 
                                            molecule_kind_set, local_molecules, molecule_set, mixed_env, my_itimes)

   TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env
   TYPE(mixed_force_type), DIMENSION(:), &
     POINTER                                :: global_forces
   REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
   REAL(KIND=dp), DIMENSION(:), POINTER     :: my_energies
   TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
   TYPE(distribution_1d_type), POINTER                :: local_molecules
   TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
   REAL(KIND=dp), ALLOCATABLE                 :: intra_forces_global(:,:,:) 
   REAL(KIND=dp), ALLOCATABLE               :: energies_1_mol_list(:,:)
   INTEGER, ALLOCATABLE                     :: used_mols(:)
   !INTEGER                                  :: act_atoms
   INTEGER, DIMENSION(:), POINTER           :: my_mol_index !CHANGE SAM
   INTEGER, POINTER                         :: my_itimes !CHANGE_AC

   TYPE(molecule_kind_type), POINTER        :: molecule_kind
   TYPE(molecule_type), POINTER             :: molecule  

   INTEGER :: i, natom, my_mol, count_idx, count2, count3, j
   INTEGER :: imol, nkind,  nmol_per_kind,  ikind, first_atom, &
              last_atom, aoi!, test

   NULLIFY( my_energies)
   ALLOCATE(my_energies(SIZE(energies_1_mol_list(1,:))))
  !> WARNING: here we assume that all force evals have the same number of atoms
  natom = SIZE(global_forces(1)%forces(:,:), 2)

  ALLOCATE(used_mols(SIZE(energies_1_mol_list(1,:))))
 !> we do F allocation just at the first step and never again just to NOT
 !allocate and deallocate this big array at every step. NOTE that this array is
 !now never deallocated and just deallocated by the compiler at the very end.
 !One should probably write a subroutine that deallocate all allocatable just
 !this at the very end of cp2k
  if(my_itimes .eq. 0) then
     ALLOCATE(mixed_env%F(3,natom, SIZE(energies_1_mol_list(1,:))))
  end if

  count2 = 1
  count3 = 0
  used_mols(:) = -1

    count_idx = 0
   DO i = 1, SIZE(energies_1_mol_list(1,:)) !loop over num of adiabatic states 
      !print *, "INTRA+LJ", energies(1), "   ", energies(2)
      !> we create site energies from the 2 force_evals sections
      my_energies(i) = energies(1) + energies_1_mol_list(2,i) - energies_1_mol_list(1,i)

!      ALLOCATE(my_global_forces(i)%forces( 3,natom))
!      my_global_forces(i)%forces = global_forces(1)%forces
 
      !> we allocate diagonal forces array. All the atoms for a give state i
      !assume the value of the neutral force_field that contains: intra +
      !nonbonded interactions
      mixed_env%F(:,:,i) = global_forces(1)%forces
      nkind = SIZE ( molecule_kind_set )
      count3 = 0
      MOL: DO ikind = 1, nkind !> loop over kind of molecules (the kind depends on the
                               !> different names given to the residues in the topology file)
          nmol_per_kind = local_molecules % n_el ( ikind )
          !> nmol_per_kind depends by the topology file
          DO imol = 1, nmol_per_kind !loop over molecules of a given kind
              my_mol = local_molecules % list ( ikind ) % array ( imol )
              molecule => molecule_set ( my_mol )
              molecule_kind => molecule % molecule_kind
              !IF (ANY(my_mol_index==my_mol) .and. count3 .eq. 0) THEN
              DO j = 1, SIZE(my_mol_index)
                   IF ( my_mol_index(j)==my_mol  .and. my_mol_index(j) /= used_mols(j) .and. count3==0 ) THEN
                    CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom )
                    DO aoi=first_atom, last_atom !loop over atoms of interest
                       count_idx = count_idx + 1
                       !print *, "force_charge_to_add_xyz", intra_forces_global(2,:,count_idx)
                       !print *, "force_neutral_to_add_xyz", intra_forces_global(1,:,count_idx)
                       !print *, "force_already_there", global_forces(1)%forces(:,aoi)
                       !my_global_forces(i)%forces(:, aoi) = my_global_forces(i)%forces(:, aoi) + &
                       !           intra_forces_global(2,:,count_idx) - intra_forces_global(1,:,count_idx)
                       !> Additiona and subtruction scheme to use just two
                       !Force_eval section instead of N and save a huge amount
                       !of memory
                       mixed_env%F(:,aoi,i) = mixed_env%F(:,aoi,i) + &
                                  intra_forces_global(2,:,count_idx) - intra_forces_global(1,:,count_idx)
                    END DO
                    used_mols(count2) = my_mol
                    count3 = 1
                    count2 = count2 + 1  
                   END IF
              END DO
          END DO
      END DO MOL
   END DO
 DEALLOCATE(used_mols)

  END SUBROUTINE add_energies_forces_mem_saver


! *****************************************************************************
!> \brief 
!> \par History
!>      SG 
! *****************************************************************************
  SUBROUTINE add_energies_forces(energies, global_forces, my_energies,  do_speedup_lj, do_speedup_intra, &
                                                            mixed_env, my_itimes)

!TYPE(force_env_type), POINTER                :: force_env

    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
   TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
!    TYPE(mixed_force_type), DIMENSION(:), &
!      POINTER                                :: my_global_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_energies
    INTEGER, POINTER                         :: my_itimes 

!    TYPE(force_env_type), POINTER            :: sub_force_envs
!    TYPE(fist_env_type), POINTER             :: charge_fist_env, &
!                                                neutral_fist_env
!    TYPE(fist_energy_type), POINTER          :: charged_thermo, &
!                                                neutral_thermo


    LOGICAL, INTENT(IN) :: do_speedup_lj, do_speedup_intra

    INTEGER :: i, natom, number_state

    
    NULLIFY( my_energies)
    if ((do_speedup_lj).AND.(.NOT.(do_speedup_intra))) then
       !ALLOCATE(my_global_forces(SIZE(energies)-1))
      !> WARNING: here we assume that all force evals have the same number of atoms
       natom = SIZE(global_forces(1)%forces(:,:), 2)
       if(my_itimes .eq. 0) then
           ALLOCATE(mixed_env%F(3,natom, SIZE(energies)-1))
       end if
       ALLOCATE(my_energies(SIZE(energies)-1))
       DO i = 2, SIZE(energies)
          my_energies(i-1) = energies(i) + energies(1)
          natom = SIZE(global_forces(i)%forces(:,:), 2)
          !ALLOCATE(my_global_forces(i-1)%forces( 3,natom))
          mixed_env%F(:,:,i-1) = global_forces(i)%forces + global_forces(1)%forces
!          my_global_forces(i-1)%forces = global_forces(i)%forces + global_forces(1)%forces
       ENDDO

    else if ((do_speedup_lj).AND.(do_speedup_intra)) then
       number_state = (SIZE(energies)-1) /2
       !ALLOCATE(my_global_forces( number_state ))
      !> WARNING: here we assume that all force evals have the same number of atoms
       natom = SIZE(global_forces(1)%forces(:,:), 2)
       if(my_itimes .eq. 0) then
           ALLOCATE(mixed_env%F(3,natom, number_state))
       end if
       ALLOCATE(my_energies(  number_state ))
       DO i = 1, number_state
          my_energies(i) =  energies(1) + energies(2*i) - energies(2*i + 1)
          natom = SIZE(global_forces(i)%forces(:,:), 2)
!          ALLOCATE(my_global_forces(i)%forces( 3,natom))
           mixed_env%F(:,:,i) = global_forces(1)%forces + global_forces(2*i)%forces - global_forces(2*i + 1)%forces
!          my_global_forces(i)%forces =  global_forces(1)%forces + global_forces(2*i)%forces - global_forces(2*i + 1)%forces
       ENDDO

    else if ((.NOT.(do_speedup_lj)).AND.(.NOT.(do_speedup_intra))) then
       !ALLOCATE(my_global_forces(SIZE(energies)))
      !> WARNING: here we assume that all force evals have the same number of atoms
       natom = SIZE(global_forces(1)%forces(:,:), 2)
       if(my_itimes .eq. 0) then
           ALLOCATE(mixed_env%F(3,natom, SIZE(energies)))
       end if
       ALLOCATE(my_energies(SIZE(energies)))
       DO i = 1, SIZE(energies)
          my_energies(i) = energies(i)
          natom = SIZE(global_forces(i)%forces(:,:), 2)
!          ALLOCATE(my_global_forces(i)%forces( 3,natom))
          mixed_env%F(:,:,i) = global_forces(i)%forces
!          my_global_forces(i)%forces = global_forces(i)%forces
        ENDDO
    ELSE
        CPABORT("Cannot do speed-up intra without speedup-lj")
    ENDIF

  END SUBROUTINE add_energies_forces

END MODULE sh_site_energies_methods
