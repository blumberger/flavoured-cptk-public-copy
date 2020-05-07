!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2016.01.15: creation, put only ffsh_aom subroutines
! *****************************************************************************

! *****************************************************************************
!>  \brief This module is used to inizialized all the sh_env variables. It
!         also initializes coupling section specific variables
!>  \authors AC and SG 
! *****************************************************************************
MODULE sh_initialization 

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


     USE parallel_rng_types,           ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             rng_stream_type
 
     USE cp_log_handling,              ONLY: cp_get_default_logger, &
                                             cp_logger_type

     USE cp_output_handling,           ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_should_output
     USE molecule_kind_types,             ONLY: molecule_kind_type
     USE molecule_types_new,              ONLY: molecule_type,&
                                                get_molecule
     USE distribution_1d_types,           ONLY: distribution_1d_type
     USE sh_site_energies_methods,        ONLY: create_ext_press_list


#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='sh_initialization'   

     PUBLIC :: create_sh_env

CONTAINS

! *****************************************************************************
!>  \brief This allocates the local sh_env variables. This local structure is
!          then passed to the more general mixed_env in order to update all the
!          variables between different time step. 
!>  \authors AC and SG 
!>  \date revision 02/02/20
! *****************************************************************************
  SUBROUTINE create_sh_env(nadiab, norbitals, natom, &
                                  positions_tp,     &
                                  adiab_section,    &
                                  sh_env, my_mol_index, ext_pres, &
                                  active_system_info, full_system_info)


     INTEGER, INTENT(IN)                            :: nadiab, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(sh_env_type), INTENT(INOUT)               :: sh_env
     !TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
     !TYPE(distribution_1d_type), POINTER                :: local_molecules
     !TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
     !INTEGER, DIMENSION(:), POINTER, OPTIONAL       :: my_mol_index !this are the decomp list indeces

     CHARACTER(len=*), PARAMETER :: routineN ='create_sh_env', &
     routineP = moduleN//':'//routineN

     TYPE(section_vals_type), POINTER               :: optimization_section, site_energy_section
     TYPE(section_vals_type), POINTER               :: ct_block
     TYPE(section_vals_type), POINTER               :: wf_restart_section, decoherence_section
     TYPE(section_vals_type), POINTER               :: frozen_coupling_section,tresp_coupling_section
     integer                                        :: mat, mol
     TYPE(cp_sll_val_type), POINTER                 :: listdiabcoeff
     REAL(KIND=dp), DIMENSION(3, 2)                 :: initial_seed
     INTEGER                                        :: input_seed
     TYPE(rng_stream_type), POINTER                 :: rng_stream
     CHARACTER(LEN=default_string_length)           :: deco_keyword
     CHARACTER(LEN=default_string_length)           ::coupling_keyword
     REAL(KIND=dp)                                  :: MDtimestep

     TYPE(cp_sll_val_type), POINTER                 :: list_tresp
     TYPE(cp_sll_val_type), POINTER                 :: frz_couplings_list
     integer, dimension(:), ALLOCATABLE             :: ext_elements
     LOGICAL, DIMENSION(:), POINTER                 :: active_logical
     INTEGER                                        :: natom_mol, ninteractions
     INTEGER, DIMENSION(:), POINTER                 :: my_mol_index 
     INTEGER, DIMENSION(:,:), ALLOCATABLE           :: full_system_info
     INTEGER, DIMENSION(:,:), ALLOCATABLE           :: active_system_info
     LOGICAL, DIMENSION(:), POINTER                 :: ext_pres

     NULLIFY(list_tresp)
     NULLIFY(frz_couplings_list)

     NULLIFY(rng_stream)
     NULLIFY(decoherence_section)
     NULLIFY(wf_restart_section)
     NULLIFY(optimization_section)
     NULLIFY(site_energy_section)

     CALL section_vals_val_get(adiab_section,"NUCLEAR_TIMESTEP",&
         r_val=MDtimestep)
     sh_env%deltat = MDtimestep !> MD time step

     ! COMMENT_SG: mol and mat would need to be looked into carefully to see if
     ! they make sense for norbitals != 1.
     mol = nadiab ! CHANGE_AC !> number of adiabatic states considered
     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     sh_env%natom = natom ! total number of atoms for a given force eval

     ALLOCATE(sh_env%forces_adiab(3,natom))
     ALLOCATE(sh_env%forces_adiab_new(3,natom))
     ALLOCATE(sh_env%diabatic_nace_tp(mat,mat))
     ALLOCATE(sh_env%diabatic_nace_t(mat,mat))

     sh_env%forces_adiab = 0.0_dp
     sh_env%forces_adiab_new = 0.0_dp
     sh_env%diabatic_nace_tp = 0.0_dp
 
     ALLOCATE(sh_env%velocities_tp(3, natom))
     !> the following arrays are allocated in adiabaticity
     !ALLOCATE(sh_env%C_t(mat,mat))
     !ALLOCATE(sh_env%C_tp(mat,mat))
     !ALLOCATE(sh_env%E_t(mat))
     !ALLOCATE(sh_env%E_tp(mat))
     ALLOCATE(sh_env%H_t(mat,mat))
     ALLOCATE(sh_env%H_tp(mat,mat))
     !> wavefunction in the diabatic basis
     ALLOCATE(sh_env%ctin(mat))
     ALLOCATE(sh_env%ctin_t(mat))
     !> wavefunction associated with decoherence
     ALLOCATE(sh_env%ctin_before_deco(mat))
     ALLOCATE(sh_env%ctin_after_deco(mat))
     ALLOCATE(sh_env%ctin_after_tranf_correction(mat))
     ALLOCATE(sh_env%active_region(mat))
   
     CALL create_adiab_nace_env(mat, sh_env)

     sh_env%futur_decohere = .FALSE.
     sh_env%decohere       = .FALSE.
     sh_env%no_diab_nace   = .FALSE.

     sh_env%new_state      = -1
     sh_env%old_state      = -1
     sh_env%active_state   = -1
     sh_env%attempt_hop    = .FALSE.
     sh_env%successful_hop = .FALSE.
     sh_env%HS_propagation = .FALSE.

     ALLOCATE(sh_env%adiab_overlap(mat, mat))
     ALLOCATE(sh_env%class_forbidden_hop_t(mat))
     allocate(sh_env%adiab_nace_t(mat,mat))
     allocate(sh_env%adiab_nace_tp(mat,mat))

     !> system info passed to sh_env
     !ALLOCATE(sh_env%index_active_mols(SIZE(my_mol_index)))
     ALLOCATE(sh_env%ext_pres(natom))
     !sh_env%index_active_mols = my_mol_index
     sh_env%full_system_info  = full_system_info
     sh_env%active_system_info = active_system_info
     sh_env%ext_pres = ext_pres
 

     sh_env%do_deco = .FALSE.
     decoherence_section => section_vals_get_subs_vals(adiab_section,"DECOHERENCE")
     CALL section_vals_val_get(decoherence_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
     IF (.NOT.(deco_keyword.EQ."NO_DECO_CORR"))  sh_env%do_deco = .TRUE.

     sh_env%do_reorder = .FALSE.
     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
     CALL section_vals_val_get(optimization_section,"REORDERING_STATES_USING_OVERLAP",&
                                                      l_val=sh_env%do_reorder)

     !> INITIALIZE THE RANDOM GENERATOR
     CALL section_vals_val_get(adiab_section,"SEED",i_val=input_seed)
     initial_seed(:,:) = REAL(input_seed,KIND=dp)
     !CALL create_rng_stream(rng_stream=sh_env%rng_stream,&
     CALL create_rng_stream(rng_stream=rng_stream,&
                             name="aom fssh",&
                             distribution_type=UNIFORM, &
                             seed=initial_seed)
     ALLOCATE(sh_env%rng_stream)
     sh_env%rng_stream = rng_stream

     !> Get wavefunction diabatic coefficients from the input file
     !COMMENT_SG: this section should be optional somehow---STILL TO DO USING e.g.
     !STEP_START_EVAL keyword in the input 
     wf_restart_section => section_vals_get_subs_vals(adiab_section,"WAVEFUNCTION_RESTART")
     CALL section_vals_val_get(wf_restart_section,"RESTART_KEY", l_val=sh_env%restart_key)
     print *, "RESTART_KEY", sh_env%restart_key
     if (sh_env%restart_key) then
         CALL section_vals_val_get(wf_restart_section,"ACTIVE_STATE_RESTART", i_val=sh_env%active_state_res)
         CALL section_vals_list_get(wf_restart_section,"_DEFAULT_KEYWORD_",list=listdiabcoeff) !Read them in the input files
         CALL read_listdiabcoeff(listdiabcoeff, norbitals, nadiab, ctin=sh_env%readcoeff) !Transform them in the complex ctin
     end if


!    !> create the list of active molecules if not present in DECOMP and the
!    !array of active atoms: ext_pres, which is use for kinetic energy
!    !calculation of the active subsystem.
!    if (PRESENT(my_mol_index)) then
!       ALLOCATE(sh_env%index_active_mols(mol))
!       sh_env%index_active_mols = my_mol_index
!    else if (.NOT. PRESENT(my_mol_index)) then 
!        !no need to allocate section_vals_val_get does that for us
!        CALL section_vals_val_get(adiab_section,"INDEX_ACTIVE_MOLS",&
!                i_vals=sh_env%index_active_mols)
!    end if
!    print *, "FOB_SH| INDEX ACTIVE MOLS:", sh_env%index_active_mols
!    CALL create_ext_press_list(molecule_kind_set, local_molecules, molecule_set, natom, &
!                                      sh_env%index_active_mols, sh_env%ext_pres, & 
!                                   sh_env%active_system_info, full_system_info=sh_env%full_system_info)
    !> sh_env%ext_pres will be overwritten in AOM case becuase active molecules
    !may still have inactive atoms (this feature is used for further speed-up)

    !sh_env%total_active_atoms = COUNT(sh_env%ext_pres)
    !print *, "FOB_SH| EXT PRESS:", sh_env%ext_pres
    !print *, "TOT ACTIVE ATOMS", COUNT(sh_env%ext_pres)

    !print *, "active sys info", sh_env%active_system_info(1, :)
    !print *, "active sys info", sh_env%active_system_info(2, :)

    !print *, "full sys info", sh_env%full_system_info(1, :)
    !print *, "full sys info", sh_env%full_system_info(2, :)

    ! Tell the program if you want a multi-state Hamiltonian: meaning an
    ! Hamiltonian like the CT one, where the number of states is different from
    ! the number of active molecules. 
     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
    ! COMMENT_SG: this keyword is temporary and in future forces calculation and
    ! adiabatic nacv calculations should be generalized to handle a multi-state
    ! hamiltoninan even without freezing the couplings
    CALL section_vals_val_get(site_energy_section,"MULTI_STATES",l_val=sh_env%multistate)
    if (sh_env%multistate) then 
       ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
       CALL section_vals_val_get(ct_block,"ADD_COULOMB_POTENTIAL", &
                                                             l_val=sh_env%add_coulomb_pot)
    end if 
  
    !sh_env%multistate = .False.
    !initialize couplings
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    sh_env%aom_couplings_keyword = .FALSE.
    sh_env%frozen_couplings_keyword = .FALSE.
    sh_env%tresp_couplings_keyword = .FALSE.
    sh_env%Gaussian_charges = .FALSE.

    SELECT CASE(coupling_keyword)
    !---------------- AOM array ---------------------------------------------------------------
    CASE("AOM")
         sh_env%aom_couplings_keyword = .TRUE.
         if (sh_env%multistate) then
            CPABORT("Multi-state forces not implemented")
         end if   
         CALL set_aom_env(nadiab, norbitals, natom,positions_tp,adiab_section,sh_env)
    !------------------------------------------------------------------------------------------
    CASE("FROZEN_COUPLINGS")
        ! Frozen coupling would have the possibility of using multi-state for
        ! the CT block for example. 
        !sh_env%multistate = .TRUE.
        sh_env%frozen_couplings_keyword = .TRUE.  
        NULLIFY(frozen_coupling_section)
        frozen_coupling_section => section_vals_get_subs_vals(adiab_section, "FROZEN_COUPLINGS")
        !> Get Hab for the first upper and lower diagonal (tridiagonal matrix
        !will be formed)
         
        CALL section_vals_val_get(frozen_coupling_section,"CONNECTIVITY_FRZ_C",&
             l_val=sh_env%connectivity_frz_c)
        IF (sh_env%connectivity_frz_c) THEN 
           CALL section_vals_val_get(frozen_coupling_section,"INTERACTION_NUMBER",&
                i_val=ninteractions)
           CALL section_vals_list_get(frozen_coupling_section,"_DEFAULT_KEYWORD_",list=frz_couplings_list) !Read list tresp from input
           CALL read_list_frz_couplings(ninteractions, frz_couplings_list, sh_env%frz_coup_connectivity, sh_env%frz_couplings) 
!print *, "frz_couplings", sh_env%frz_coup_connectivity(:,1)
!print *, "frz_couplings", sh_env%frz_coup_connectivity(:,2)
!print *, "sh_env%frz_couplings", sh_env%frz_couplings
        ELSE
           CALL section_vals_val_get(frozen_coupling_section,"FIRST_OFF_DIAGONALS_VAL",&
             r_val=sh_env%first_off_d_Hab)
        END IF
    CASE("TRESP_COUPLINGS")
         if (sh_env%multistate) then
            CPABORT("Multi-state forces not implemented")
         end if   
        sh_env%tresp_couplings_keyword = .TRUE.
        !sh_env%frozen_couplings_keyword = .TRUE.  !test
        !sh_env%first_off_d_Hab = -0.0011024789     !test

        !----------- put in a new subroutine -----------------------------
        NULLIFY(tresp_coupling_section)
        tresp_coupling_section => section_vals_get_subs_vals(adiab_section, "TRESP_COUPLINGS")
        !> READ TRESP CHARGES TEST
        CALL section_vals_list_get(tresp_coupling_section,"_DEFAULT_KEYWORD_",list=list_tresp) !Read list tresp from input
        CALL section_vals_val_get(tresp_coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=sh_env%cutoff_sites)
        ! determine if the TRESP are taken from Gaussian code and divide by
        ! sqrt(2) if that's the case (see Multiwfn.pdf for explanation)
        CALL section_vals_val_get(tresp_coupling_section,"TRESP_FROM_GAUSSIAN",&
             l_val=sh_env%Gaussian_charges)
        !### COMMENT_SG: not general for multiphase: but it could be easily
        !charged by giving different charges for different molecules
        CALL section_vals_val_get(tresp_coupling_section,"NUMBER_AOM_ATOMS_PER_SITE",&
             i_val=natom_mol)
        sh_env%atoms_per_site = natom_mol

        !> COMMENT_SG: sh_env%tresp_charges are the tresp for a single molecule, this would need to be generalized
                      !for multiple molecule kind so that the coupling
                      !calculation is general for any num atoms per molecule
        CALL read_list_tresp(sh_env%atoms_per_site, list_tresp, ext_elements, active_logical, sh_env%tresp_charges) 
        !> active_logical is not necessary since inactive atoms will anyway have 0
        !charges that will nullify the term in the coulomb sum 
        
        !c = nadiab*sh_env%atoms_per_site !> Numbers of active atoms per mol 
        allocate(sh_env%crds(3,COUNT(sh_env%ext_pres),3))  !> crds active atoms
        ALLOCATE(sh_env%NN(mat, mat))
        !COMMENT_SG: COUNT(sh_env%ext_pres) is the total number of active atoms,
        !note that I did not use mat*active_atoms_per_mol to be more general for
        !multi-phase
        ALLOCATE(sh_env%grad_tresp(3,mat, COUNT(sh_env%ext_pres)))
!TEST
        !ALLOCATE(sh_env%test_grad_tresp(3,mat, mat, COUNT(sh_env%ext_pres)))
 

    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT

     DEALLOCATE(rng_stream)
  END SUBROUTINE create_sh_env


! *****************************************************************************
!> \brief This allocates aom_env variables. One could in principle put these
!          variables in a separate derive type structure
!> \author  SG 
!> \date 02/02/20
! *****************************************************************************
  SUBROUTINE set_aom_env(nadiab,  norbitals, natom, &
                                  positions_tp,     &
                                  adiab_section,    &
                                  sh_env)


     INTEGER, INTENT(IN)                            :: nadiab, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(sh_env_type), INTENT(INOUT)               :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='set_aom_env', &
     routineP = moduleN//':'//routineN

     TYPE(section_vals_type), POINTER               :: optimization_section
     integer                                        :: i, j, mat, c, coord, mol, start
     REAL(KIND=dp), dimension(:,:), allocatable     :: crds_loc
     REAL(KIND=dp), dimension(:,:), allocatable     :: Atp_loc
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot, scoefftot
     LOGICAL, DIMENSION(:), ALLOCATABLE             :: atomlist_mask, submask
     integer, dimension(:), ALLOCATABLE             :: ext_elements
     TYPE(cp_sll_val_type), POINTER                 :: listcoeff
     LOGICAL                                        :: continue_while
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma, file_name_s_ppi,&
                                                       prop_wf_keyword
     ! OZ: extra strings for C-S and S-S SAB calc
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_C_S,file_name_s_ppi_C_S
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_S_S,file_name_s_ppi_S_S
     TYPE(section_vals_type), POINTER               :: coupling_section
     INTEGER                                        :: aom_atoms

     NULLIFY(listcoeff)
     NULLIFY(optimization_section)
     NULLIFY(coupling_section)

     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")

     CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword) 

     !> Initialize input parameters: 
     !> cutoff to define the neighborough list

     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=sh_env%cutoff_connect)
     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=sh_env%cutoff_sites)
     !> minao is an atomic overlap cutoff
     CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=sh_env%minao)
     !> dR is the increment for the gradient finite difference calculation
     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
             r_val=sh_env%dR)
     !> analytics would probably calculates analytic forces for a 2 state sytems
     !(Spencer 2016)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=sh_env%analytics)
     !> This remove bias coming from the finite difference calculation of the
      !> gradients of the Hamiltonian from AOM
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=sh_env%do_random)
     !> Multiple time step calculation of aom gradient
     CALL section_vals_val_get(coupling_section,"MULTIPLE_TIME_STEP",&
             i_val=sh_env%multiple_time_step)

     !> Create the overlap look-up table for the AOM coupling calculation
     CALL section_vals_val_get(coupling_section,"PSIGMA_FILE_NAME",&
             c_val=file_name_s_psigma)
     CALL section_vals_val_get(coupling_section,"PPI_FILE_NAME",&
             c_val=file_name_s_ppi)
     ! OZ: C-S and S-S
     CALL section_vals_val_get(coupling_section,"PSIGMA_C_S_FILE_NAME",&
             c_val=file_name_s_psigma_C_S)
     CALL section_vals_val_get(coupling_section,"PPI_C_S_FILE_NAME",&
             c_val=file_name_s_ppi_C_S)
     CALL section_vals_val_get(coupling_section,"PSIGMA_S_S_FILE_NAME",&
             c_val=file_name_s_psigma_S_S)
     CALL section_vals_val_get(coupling_section,"PPI_S_S_FILE_NAME",&
             c_val=file_name_s_ppi_S_S)
     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=sh_env%scaling)
     !> cbar factor : Sab = cbar*Sab_bar
     CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=sh_env%cbar)
     !> number of atoms per site (This array is not general for a double phase
     !with molecules with a different number of atoms per site)
     CALL section_vals_val_get(coupling_section,"NUMBER_AOM_ATOMS_PER_SITE",&
         i_val=aom_atoms)


     mol = nadiab ! CHANGE_AC
     sh_env%atoms_per_site = aom_atoms ! CHANGE_AC
     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     c = nadiab*aom_atoms !> Numbers of active atoms


     ALLOCATE(sh_env%elements(c))
     ALLOCATE(sh_env%connlist1(4,aom_atoms))
     allocate(sh_env%aneighbour(0:sh_env%atoms_per_site,0:sh_env%atoms_per_site))
     allocate(sh_env%bneighbour(0:sh_env%atoms_per_site,0:sh_env%atoms_per_site))
     ALLOCATE(sh_env%pcoeffs(aom_atoms,norbitals))
     ALLOCATE(sh_env%scoeffs(aom_atoms))
     !> crds is an AOM arrays (xyz, atoms, tm-t-tp) the third dimension
     !indicates what time the coordinates refer to
     allocate(sh_env%crds(3,c,3))
     ALLOCATE(sh_env%NN(mat, mat))
     allocate(sh_env%pcoeffsrun(sh_env%atoms_per_site,mat,3)) !pcoeffsrun exists because the p-coefficients
                                               !need to be renormalised for each structure
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(sh_env%At(4,sh_env%atoms_per_site))
     allocate(sh_env%Atp(4,sh_env%atoms_per_site))
     allocate(sh_env%Atm(4,sh_env%atoms_per_site))
     allocate(sh_env%Bt(4,sh_env%atoms_per_site))
     allocate(sh_env%Btp(4,sh_env%atoms_per_site))
     allocate(sh_env%Btm(4,sh_env%atoms_per_site))
     allocate(sh_env%siteenergies(3,mat))  !COMMENT_SG: First dimension is to be removed
     ALLOCATE(sh_env%nacv(3,mat, mat*aom_atoms))
     ALLOCATE(sh_env%rm_t(3, natom))
     ALLOCATE(sh_env%rm_tm(3, natom))
     ALLOCATE(sh_env%Stt(mat,mat))
     ALLOCATE(sh_env%Stptp(mat,mat))
     if (prop_wf_keyword .ne. "DIABATIC_RK_FAST") then
     !> in case DIABATIC_RK_FAST diabatic nace are never calculated, neither in
     !> the
     !> propagate_wf nor in the adiabatic nace calculation for the hopping
     !> probability
     !> so that all the Sttp and co matrices are also not used, no point in
     !> allocate them
        print *, "FOB-SH: ALLOCATION OF S MATRICES"
        ALLOCATE(sh_env%Sttp(mat,mat))
        ALLOCATE(sh_env%Stpt(mat,mat))
        ALLOCATE(sh_env%Stmt(mat,mat))
        ALLOCATE(sh_env%Sttm(mat,mat))
        ALLOCATE(sh_env%Stmtm(mat,mat))
     end if

     sh_env%At = 0.0_dp
     sh_env%Atp = 0.0_dp
     sh_env%Atm = 0.0_dp
     sh_env%Bt  = 0.0_dp
     sh_env%Btp = 0.0_dp
     sh_env%Btm = 0.0_dp
     sh_env%siteenergies = 0.0_dp
     sh_env%nacv = 0.0_dp
     !> these arrays are data for the diabatic nacv calculation
     !> they help efficiency by telling the program which atoms are connected
     !> and thus which pairs to carry out the calculation for
     sh_env%aneighbour = 0
     sh_env%bneighbour = 0

     ! OZ: add extra arguments...
     CALL build_s_tables(file_name_s_psigma, file_name_s_ppi,&
                         file_name_s_psigma_C_S, file_name_s_ppi_C_S,&
                         file_name_s_psigma_S_S, file_name_s_ppi_S_S)
     
     
     !> READ THE AOM COEFFICIENT
     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     CALL read_listcoeff(listcoeff, ext_elements, sh_env%ext_pres, &
                         scoefftot, pcoefftot, natom) !Transform them into four arrays

     sh_env%elements = PACK(ext_elements, sh_env%ext_pres) !Only select active atoms
     allocate(atomlist_mask(c))
     atomlist_mask =.NOT.((scoefftot.eq.0).AND.(pcoefftot.eq.0)) !Only select atoms with orbitals
     allocate(submask(aom_atoms))
     i = 1
     continue_while = .True. 
     do while (continue_while)
        if (atomlist_mask(i)) then
           start = i
           continue_while = .False.
         endif
         i = i+1
     end do
     do i=1, aom_atoms 
        submask(i) = atomlist_mask(i + start - 1)
     enddo

     allocate(sh_env%atomlist(COUNT(submask(1:aom_atoms))))

     sh_env%atomlist = PACK( (/ (i, i=1,aom_atoms) /), submask(1:aom_atoms) ) !List atoms with orbitals
                                                                              !namely,
                                                                              !with
                                                                              !coeff
                                                                              !not
                                                                              !zero
     allocate(sh_env%connlist2(4,size(sh_env%atomlist)))
     do i=1,aom_atoms
          do j=1,norbitals
                sh_env%pcoeffs(i,j) = pcoefftot( start -1 + i + (j-1)*norbitals)
          end do
          sh_env%scoeffs(i) = scoefftot(start - 1 + i)
     end do

     !> CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     allocate(crds_loc(3, c))
     ALLOCATE(Atp_loc(4,aom_atoms))

     do coord=1,3
        crds_loc(coord,:) = PACK(positions_tp(coord,:), sh_env%ext_pres)
     enddo

     Atp_loc(2:4,1:aom_atoms) = crds_loc(1:3,1:aom_atoms)
     !output: sh_env%connlist1, sh_env%connlist2
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist, sh_env%connlist1, sh_env%connlist2, sh_env%cutoff_connect)
     

     DEALLOCATE(atomlist_mask) 
     DEALLOCATE(submask) 
     DEALLOCATE(Atp_loc)
     DEALLOCATE(crds_loc)
     DEALLOCATE(ext_elements)
     DEALLOCATE(scoefftot)
     DEALLOCATE(pcoefftot)

  END SUBROUTINE set_aom_env

! *****************************************************************************
!> \brief allocates nace variables
!> \author  MEllis 
! *****************************************************************************
  SUBROUTINE create_adiab_nace_env(mat, sh_env)
   
     INTEGER                       :: mat
     TYPE(sh_env_type)            :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_adiab_nace_env', &
     routineP = moduleN//':'//routineN

      ALLOCATE(sh_env%adiab_nace_env%T1(mat,mat))
      ALLOCATE(sh_env%adiab_nace_env%T2(mat,mat))
      ALLOCATE(sh_env%adiab_nace_env%DT(mat,mat))
      ALLOCATE(sh_env%adiab_nace_env%diab_D(mat,mat))
      ALLOCATE(sh_env%adiab_nace_env%DS(mat,mat))
      ALLOCATE(sh_env%adiab_nace_env%map_new2old(mat))
      ALLOCATE(sh_env%adiab_nace_env%map_old2new(mat))
      ALLOCATE(sh_env%adiab_nace_env%basic(mat))
      ALLOCATE(sh_env%adiab_nace_env%zeros(mat))
      ALLOCATE(sh_env%adiab_nace_env%not_associated(mat))
      ALLOCATE(sh_env%adiab_nace_env%Ssqrt(mat,mat))

  END SUBROUTINE create_adiab_nace_env



END MODULE sh_initialization
