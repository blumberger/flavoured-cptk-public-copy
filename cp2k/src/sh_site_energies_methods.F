
! *****************************************************************************
!> \brief 
!> \par History
!>      SG 
! *****************************************************************************
MODULE sh_site_energies_methods 


     USE kinds,                        ONLY: default_string_length,               &
                                             dp
     USE aom_overlapfinal,             ONLY:  connect_list2,                      &
                                              build_s_tables

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


     USE cp_log_handling,              ONLY: cp_get_default_logger, &
                                             cp_logger_type

     USE cp_output_handling,           ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_should_output
     USE molecule_kind_types,             ONLY: molecule_kind_type
     USE molecule_types_new,              ONLY: molecule_type,&
                                                get_molecule
     USE distribution_1d_types,           ONLY: distribution_1d_type

     USE mixed_environment_types,         ONLY: get_mixed_env,&
                                                 mixed_env_retain,&
                                                 mixed_environment_type
     USE mixed_energy_types,              ONLY: mixed_energy_type,&
                                                 mixed_force_type

#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='sh_site_energies_methods'   

     PUBLIC :: ct_block_energies_forces, &
               add_energies_forces_mem_saver_V2, &
               add_energies_forces, &
               create_ext_press_list, &
               add_energies_forces_mem_saver, &
               multi_block_energies_forces 

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
   TYPE(section_vals_type), POINTER          :: ct_block, fe_block
   INTEGER, DIMENSION(:), POINTER            :: donor_index, acceptor_index, exc_index
   INTEGER                                   :: total_states, state,nadiab, natom
   INTEGER, DIMENSION(:,:), ALLOCATABLE      :: intra_forces_info
   INTEGER, DIMENSION(:), ALLOCATABLE        :: first_ind, second_ind, exc_ind
   INTEGER                                   :: count_
   INTEGER, DIMENSION(:), INTENT(IN)         :: my_mol_index !this is the decomp index
 
   NULLIFY(site_energy_section,ct_block, fe_block)
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

   if (fe_block_keyword) then 
       fe_block => section_vals_get_subs_vals(site_energy_section,"FE_BLOCK")
       CALL section_vals_val_get(fe_block,"INDEX_EXC_MOLS",&
                        i_vals=exc_index)
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
         my_energies(i) = energies(1) + energies_1_mol_list(3,ind1) - energies_1_mol_list(1,ind1)
   
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
   if (cr_block_keyword) then
      ! in case CR is present
      total_states = nadiab
      my_energies(total_states) = energies(1)
      mixed_env%F(:,:,total_states) = global_forces(1)%forces
   end if
   

   DEALLOCATE(intra_forces_info)
   DEALLOCATE(first_ind)
   DEALLOCATE(second_ind)

END SUBROUTINE multi_block_energies_forces

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
!SAM COMMENT: MAYBE FIND A BETTER WAY TO ALLOCATE THIS
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
