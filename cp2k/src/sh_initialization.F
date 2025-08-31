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
     !                                         build_s_tables, &
                                              STO_initialize_splines

     USE aom_utils,                    ONLY: read_listcoeff,                      &
                                             read_listdiabcoeff,                  &
                                             read_list_tresp,                     &
                                             calc_multiphase_Xmatrix_elements,    &
                                             read_list_frz_couplings,             &
					     read_list_coulomb_barrier   

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
                                             rng_stream_type,  &
                                             next_random_number
 
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
     INTEGER, DIMENSION(:), POINTER                 :: my_mol_index !this are the decomp list indeces 

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
     CHARACTER(LEN=default_string_length)           :: coupling_keyword
     REAL(KIND=dp)                                  :: MDtimestep

     TYPE(cp_sll_val_type), POINTER                 :: list_tresp
     TYPE(cp_sll_val_type), POINTER                 :: frz_couplings_list
     TYPE(cp_sll_val_type), POINTER                 :: coulomb_barrier_list
     REAL(KIND=dp), DIMENSION(:), POINTER           :: coulomb_barrier_contribution

     integer, dimension(:), ALLOCATABLE             :: ext_elements
     LOGICAL, DIMENSION(:), POINTER                 :: active_logical
     INTEGER                                        :: natom_mol, ninteractions, nsites
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

     ! COMMENT: mol and mat would need to be looked into carefully to see if
     ! they make sense for norbitals != 1 (namenly, degeneracy would probably
     ! not work with the current CODE).

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
     sh_env%index_active_mols = my_mol_index
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

     !> Get wavefunction diabatic coefficients from the input file for the
     !restart file option
     !COMMENT_SG: this section should be optional somehow---STILL TO IMPROVE USING e.g.
     !STEP_START_EVAL keyword in the input 
     wf_restart_section => section_vals_get_subs_vals(adiab_section,"WAVEFUNCTION_RESTART")
     CALL section_vals_val_get(wf_restart_section,"RESTART_KEY", l_val=sh_env%restart_key)
     print *, "RESTART_KEY", sh_env%restart_key
     if (sh_env%restart_key) then
         CALL section_vals_val_get(wf_restart_section,"ACTIVE_STATE_RESTART", i_val=sh_env%active_state_res)
         CALL section_vals_list_get(wf_restart_section,"_DEFAULT_KEYWORD_",list=listdiabcoeff) !Read them in the input files
         CALL read_listdiabcoeff(listdiabcoeff, norbitals, nadiab, ctin=sh_env%readcoeff) !Transform them in the complex ctin
     end if

    ! Tell the program if you want a multi-state Hamiltonian: meaning an
    ! Hamiltonian like the CT one, where the number of states is different from
    ! the number of active molecules. 
     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
    ! COMMENT_SG: this keyword is temporary and in future forces calculation and
    ! adiabatic nacv calculations should be generalized to handle a multi-state
    ! hamiltoninan even without freezing the couplings
    CALL section_vals_val_get(site_energy_section,"MULTI_STATES",l_val=sh_env%multistate)

    !if multi-state, this checks if you want to add a static coulomb
    !potential

    if (sh_env%multistate) then 
       ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
       CALL section_vals_val_get(ct_block,"ADD_COULOMB_POTENTIAL", &
                                                             l_val=sh_env%add_coulomb_pot)
       if (sh_env%add_coulomb_pot) then

       CALL section_vals_val_get(ct_block,"SITE_NUMBER",&
                   i_val=nsites)

       sh_env%number_ct_states = nsites

       CALL section_vals_list_get(ct_block,"_DEFAULT_KEYWORD_",list=coulomb_barrier_list)

       CALL read_list_coulomb_barrier(nsites, coulomb_barrier_list, sh_env%coulomb_barrier)

       !added three new variables to the sh_env type: coulomb_barrier = 1D
       !array of coulomb barrier values, order corresponds to order of CT
       !states in X-SH electronic Hamiltonian; number_ct_states = number of
       !ct states in the system, passed in from the SITE_NUMBER keyword,
       !used to update the site energies of the CT states in sh_main.F;
       !energy_offset = ground-state energy difference between interfacial CT-state
       !(set to zero) and XT manifold -- !FI 05/01/24
       
       end if

       !Added the conditional coulomb barrier section in the ct_block
       !section so it's grouped up with the other site energy keywords --
       !FI 05/01/24 

    end if
  
    !sh_env%multistate = .False.
    !initialize couplings
    CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
    sh_env%aom_couplings_keyword = .FALSE.
    sh_env%frozen_couplings_keyword = .FALSE.
    sh_env%tresp_couplings_keyword = .FALSE.
    sh_env%sp_fect_keyword = .FALSE.

    !WTP_for_X-SH: !! .TRUE.
    !!sh_env%Gaussian_charges = .FALSE. !!
    sh_env%Gaussian_charges = .TRUE. !for X-SH should be True..

    sh_env%multiphase_couplings_keyword = .FALSE.
    sh_env%multifect_couplings_keyword = .FALSE.

    SELECT CASE(coupling_keyword)
    !---------------- AOM array ---------------------------------------------------------------
    CASE("AOM")
         sh_env%aom_couplings_keyword = .TRUE.
         if (sh_env%multistate) then
            CPABORT("Multi-state forces not implemented")
         end if   
         CALL set_aom_env(nadiab, norbitals, natom,positions_tp,adiab_section,sh_env)
         !CALL set_aom_env_multip(nadiab, norbitals, natom,positions_tp,adiab_section,sh_env)
    !------------------------------------------------------------------------------------------
    CASE("MULTI-PHASE")
         !sh_env%aom_couplings_keyword = .TRUE.
         sh_env%multiphase_couplings_keyword = .TRUE.  
         
         !WARNING: hardcoded for testing (you can put it in input)
         sh_env%freeze_multip_forces = .FALSE.
         CALL set_aom_env_multip(nadiab, norbitals, natom,positions_tp,adiab_section,sh_env)
    !------------------------------------------------------------------------------------------
    CASE("MULTI-FECT")
         sh_env%multifect_couplings_keyword = .TRUE.

         CALL section_vals_val_get(ct_block,"XT_CT_OFFSET",&
                   r_val=sh_env%energy_offset)

         sh_env%freeze_multip_forces = .FALSE.  !??? do I need this
         !> FE block (TRESP) is also intialized in the following subroutine

         !CALL read_list_frz_couplings(ninteractions, frz_couplings_list, sh_env%frz_coup_connectivity, sh_env%frz_couplings) 
 
         CALL set_aom_env_fect(nadiab, norbitals, natom,positions_tp,adiab_section,sh_env)

    !------------------------------------------------------------------------------------------
    CASE("SINGLE-PHASE-FECT")

         sh_env%sp_fect_keyword = .TRUE.
         sh_env%freeze_multip_forces = .FALSE.

         CALL section_vals_val_get(ct_block,"XT_CT_OFFSET",&
                   r_val=sh_env%energy_offset)

         CALL set_aom_env_sp_fect(nadiab, norbitals, &
         natom,positions_tp,adiab_section,sh_env, my_mol_index)
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
        !> determine if the TRESP are taken from Gaussian code and divide by
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
        ALLOCATE(sh_env%grad_tresp(3,mat, COUNT(sh_env%ext_pres)))
!TEST
        !ALLOCATE(sh_env%test_grad_tresp(3,mat, mat, COUNT(sh_env%ext_pres)))

    CASE DEFAULT
        CPABORT("There is no method to calculate the coupling")
    END SELECT

     DEALLOCATE(rng_stream)
  END SUBROUTINE create_sh_env

! *****************************************************************************
!> \brief This allocates aom_env variables for a Multi-phase Hamiltonian. 
!          One could in principle put these
!          variables in a separate derive type structure instead of sh_env
!  \Note this subroutine will NOT work for orb !=1 (namenly multiple orbitals
!        per molecule) 
!> \author  SG 
!> \date XX/XX/21
! *****************************************************************************
  SUBROUTINE set_aom_env_multip(nadiab,  norbitals, natom, &
                                  positions_tp,     &
                                  adiab_section,    &
                                  sh_env)


     INTEGER, INTENT(IN)                            :: nadiab, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(sh_env_type), INTENT(INOUT)               :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='set_aom_env_multip', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_d, scoefftot_d, &
                                                       pcoefftot_a, scoefftot_a
     TYPE(section_vals_type), POINTER               :: optimization_section
     integer                                        :: mat, act_atoms,i
     integer, dimension(:), ALLOCATABLE             :: ext_elements_d, ext_elements_a
     LOGICAL, dimension(:), POINTER                 :: active_logical
     TYPE(cp_sll_val_type), POINTER                 :: listcoeff
     !LOGICAL                                        :: continue_while
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma, file_name_s_ppi!,&
                                                       !prop_wf_keyword
     ! OZ: extra strings for C-S and S-S SAB calc
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_C_S,file_name_s_ppi_C_S
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_S_S,file_name_s_ppi_S_S
     TYPE(section_vals_type), POINTER               :: coupling_section, site_energy_section, ct_block
     !INTEGER                                        :: aom_atoms, mol
     integer                                        :: atoms_per_donor, atoms_per_acceptor
     INTEGER, DIMENSION(:), POINTER                 :: donor_index, acceptor_index
     REAL(KIND=dp), dimension(:,:), allocatable     :: Atp_loc
     LOGICAL, DIMENSION(:), ALLOCATABLE             :: submask
     INTEGER                                        :: init_atom_mol_d, final_atom_mol_d, atoms_per_mol_d
     INTEGER                                        :: init_atom_mol_a, final_atom_mol_a, atoms_per_mol_a
     CHARACTER(LEN=default_string_length)           :: STO_input_string
     REAL(KIND=dp)                                  :: STO_rmin,STO_rmax,STO_dr
     INTEGER                                        :: STO_debug

     NULLIFY(listcoeff)
     NULLIFY(optimization_section)
     NULLIFY(coupling_section)
     NULLIFY(site_energy_section, ct_block)

     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")

     !get CT block for donor and acceptor
     ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")

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

     !> Create the overlap look-up table for the AOM coupling calculation
     !CALL section_vals_val_get(coupling_section,"PSIGMA_FILE_NAME",&
     !        c_val=file_name_s_psigma)
     !CALL section_vals_val_get(coupling_section,"PPI_FILE_NAME",&
     !        c_val=file_name_s_ppi)
     !! OZ: C-S and S-S
     !CALL section_vals_val_get(coupling_section,"PSIGMA_C_S_FILE_NAME",&
     !        c_val=file_name_s_psigma_C_S)
     !CALL section_vals_val_get(coupling_section,"PPI_C_S_FILE_NAME",&
     !        c_val=file_name_s_ppi_C_S)
     !CALL section_vals_val_get(coupling_section,"PSIGMA_S_S_FILE_NAME",&
     !        c_val=file_name_s_psigma_S_S)
     !CALL section_vals_val_get(coupling_section,"PPI_S_S_FILE_NAME",&
     !        c_val=file_name_s_ppi_S_S)

     CALL section_vals_val_get(coupling_section,"STO_EXPONENTS",&
             c_val=STO_input_string)
     CALL section_vals_val_get(coupling_section,"STO_RMIN",&
             r_val=STO_rmin)
     CALL section_vals_val_get(coupling_section,"STO_RMAX",&
             r_val=STO_rmax)
     CALL section_vals_val_get(coupling_section,"STO_DR",&
             r_val=STO_dr)
     CALL section_vals_val_get(coupling_section,"STO_DEBUG",&
             i_val=STO_debug)


     !!> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_DONOR",&
             r_val=sh_env%scaling_donor)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_ACCEPTOR",&
             r_val=sh_env%scaling_acceptor)
          
     !CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword) 
     !!> Multiple time step calculation of aom gradient
     !CALL section_vals_val_get(coupling_section,"MULTIPLE_TIME_STEP",&
     !        i_val=sh_env%multiple_time_step)
     
     ! passing N_DONOR_ATOMS and N_ACCEPTOR_ATOMS is not strictly useful since
     ! this info is already in full_system_info or active_system_info 
     CALL section_vals_val_get(ct_block,"N_DONOR_ATOMS",&
              i_val=atoms_per_donor) 
     CALL section_vals_val_get(ct_block,"N_ACCEPTOR_ATOMS",&
              i_val=atoms_per_acceptor)

     CALL section_vals_val_get(ct_block,"INDEX_DONOR_MOLS",&
                 i_vals=donor_index)
     CALL section_vals_val_get(ct_block,"INDEX_ACCEPTOR_MOLS",&
                i_vals=acceptor_index)


     !> READ THE AOM COEFFICIENT of DONOR and ACCEPTOR
     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     !Transform them into four arrays. In this case external presence is not read anymore.
     CALL read_listcoeff(listcoeff, ext_elements_d, active_logical, &
                         scoefftot_d, pcoefftot_d, atoms_per_donor) !donor
     CALL read_listcoeff(listcoeff, ext_elements_a, active_logical, &
                         scoefftot_a, pcoefftot_a, atoms_per_acceptor) !acceptor

     !get the first mol of the donor and the first of the acceptor phase 
     init_atom_mol_d = sh_env%full_system_info(1,donor_index(1)) !just get the first donor
     final_atom_mol_d = sh_env%full_system_info(2,donor_index(1))
     atoms_per_mol_d = (final_atom_mol_d-init_atom_mol_d)+1

     init_atom_mol_a = sh_env%full_system_info(1,acceptor_index(1)) !just get the first acceptor
     final_atom_mol_a = sh_env%full_system_info(2,acceptor_index(1))
     atoms_per_mol_a = (final_atom_mol_a-init_atom_mol_a)+1

     ! sanity checks
     IF (atoms_per_mol_a .ne. atoms_per_acceptor) THEN 
        CPABORT("Error in number of atoms acceptor")
     ENDIF 
     IF (atoms_per_mol_d .ne. atoms_per_donor) THEN 
        CPABORT("Error in number of atoms donor")
     ENDIF 

      IF (norbitals .ne. 1) THEN
        CPABORT("DEGENERACY IS NOT IMPLEMENTED FOR MULTI-PHASE")
      ENDIF

     !norbitals must be 1: we have not implemented degeneracy in the ffollowing
     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     act_atoms = COUNT(sh_env%ext_pres) !> Numbers of active atoms

     !_d represents the donor phase and _a the acceptor phase
     !allocate(sh_env%aneighbour_d(0:atoms_per_mol_d,0:atoms_per_mol_d))
     !allocate(sh_env%bneighbour_d(0:atoms_per_mol_d,0:atoms_per_mol_d))
     !allocate(sh_env%aneighbour_a(0:atoms_per_mol_a,0:atoms_per_mol_a))
     !allocate(sh_env%bneighbour_a(0:atoms_per_mol_a,0:atoms_per_mol_a))

     ALLOCATE(sh_env%nacv_dphase(3,mat, mat*atoms_per_mol_d))
     ALLOCATE(sh_env%nacv_aphase(3,mat, mat*atoms_per_mol_a))
     sh_env%nacv_dphase = 0.0_dp
     sh_env%nacv_aphase = 0.0_dp

     !!> crds is an AOM arrays (xyz, atoms, tm-t-tp) the third dimension
     !!indicates what time the coordinates refer to
     !allocate(sh_env%crds(3,act_atoms,3))
     ALLOCATE(sh_env%NN(SIZE(sh_env%index_active_mols), SIZE(sh_env%index_active_mols)))
     allocate(sh_env%siteenergies(3,mat))  !COMMENT_SG: First dimension is not usefull
     sh_env%siteenergies = 0.0_dp

     ! OZ: add extra arguments...
     CALL STO_initialize_splines(STO_input_string,STO_rmin,STO_rmax,STO_dr,&
                                 STO_debug)
     ! OZ: add extra arguments...
     !CALL build_s_tables(file_name_s_psigma, file_name_s_ppi,&
     !                    file_name_s_psigma_C_S, file_name_s_ppi_C_S,&
     !                    file_name_s_psigma_S_S, file_name_s_ppi_S_S)
     
     !> READ THE AOM COEFFICIENT DONOR
     !> CONSTRUCT THE CONNLIST 
     ! SG: strange way to proceed to get the connectivity of one molecule (which
     ! is then assumed to be the same for all molecules) 
     !we fill out the connlist arrays before calculating the nacv
     allocate(submask(atoms_per_donor))   !temp array
     submask =.NOT.((scoefftot_d.eq.0).AND.(pcoefftot_d.eq.0)) !Only select atoms with orbitals

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_d(COUNT(submask(1:atoms_per_donor))))
     sh_env%atomlist_d = PACK( (/ (i, i=1,atoms_per_donor) /), submask(1:atoms_per_donor) ) 
     allocate(sh_env%connlist2_d(4,size(sh_env%atomlist_d)))
     ALLOCATE(sh_env%connlist1_d(4,atoms_per_donor))

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2 for agiven phase
     ALLOCATE(Atp_loc(4,atoms_per_donor)) !temp array
     Atp_loc(1,1:atoms_per_donor) = dble(ext_elements_d)
     Atp_loc(2:4,1:atoms_per_donor) = positions_tp(1:3,init_atom_mol_d:final_atom_mol_d)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_d, sh_env%connlist1_d, sh_env%connlist2_d, sh_env%cutoff_connect)
     
     ! coefficients donor phase
     ALLOCATE(sh_env%pcoeffs_d(atoms_per_donor))
     ALLOCATE (sh_env%scoeffs_d(atoms_per_donor))
     ALLOCATE (sh_env%ext_elements_d(atoms_per_donor))
     sh_env%pcoeffs_d = pcoefftot_d
     sh_env%scoeffs_d = scoefftot_d
     sh_env%ext_elements_d = ext_elements_d
     DEALLOCATE(ext_elements_d)
     DEALLOCATE(scoefftot_d)
     DEALLOCATE(pcoefftot_d)
     DEALLOCATE(submask) 
     DEALLOCATE(Atp_loc)

     !> READ THE AOM COEFFICIENT ACCEPTOR PHASE! 
     !> CONSTRUCT THE CONNLIST 
     allocate(submask(atoms_per_acceptor))   !temp array
     submask =.NOT.((scoefftot_a.eq.0).AND.(pcoefftot_a.eq.0)) !Only select atoms with orbitals

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_a(COUNT(submask(1:atoms_per_acceptor))))
     sh_env%atomlist_a = PACK( (/ (i, i=1,atoms_per_acceptor) /), submask(1:atoms_per_acceptor) ) 
     allocate(sh_env%connlist2_a(4,size(sh_env%atomlist_a)))
     ALLOCATE(sh_env%connlist1_a(4,atoms_per_acceptor))

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2
     ALLOCATE(Atp_loc(4,atoms_per_acceptor)) !temp array
     Atp_loc(1,1:atoms_per_acceptor) = dble(ext_elements_a)
     Atp_loc(2:4,1:atoms_per_acceptor) = positions_tp(1:3,init_atom_mol_a:final_atom_mol_a)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_a, sh_env%connlist1_a, sh_env%connlist2_a, sh_env%cutoff_connect)
     
     ! coefficients acceptor phase
     ALLOCATE(sh_env%pcoeffs_a(atoms_per_acceptor))
     ALLOCATE (sh_env%scoeffs_a(atoms_per_acceptor))
     ALLOCATE (sh_env%ext_elements_a(atoms_per_acceptor))
     sh_env%pcoeffs_a = pcoefftot_a
     sh_env%scoeffs_a = scoefftot_a
     sh_env%ext_elements_a = ext_elements_a
     DEALLOCATE(ext_elements_a)
     DEALLOCATE(submask) 
     DEALLOCATE(Atp_loc)
     DEALLOCATE(scoefftot_a)
     DEALLOCATE(pcoefftot_a)
     

  END SUBROUTINE set_aom_env_multip

  !> for X-SH: WTP_merge may require some changes here...
  SUBROUTINE set_aom_env_fect(nadiab,  norbitals, natom, &
                                  positions_tp,     &
                                  adiab_section,    &
                                  sh_env)


     INTEGER, INTENT(IN)                            :: nadiab, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(sh_env_type), INTENT(INOUT)               :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='set_aom_env_fect', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_d, scoefftot_d, &
                                                       pcoefftot_a, scoefftot_a
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_fe, scoefftot_fe, &
                                                       pcoefftot_ct, scoefftot_ct
!     TYPE(section_vals_type), POINTER               :: optimization_section
     integer                                        :: mat, i !,act_atoms
     integer, dimension(:), ALLOCATABLE             :: ext_elements_d, ext_elements_a
     integer, dimension(:), ALLOCATABLE             :: ext_elements_fe, ext_elements_ct
     integer, dimension(:), ALLOCATABLE             :: ext_elements
     CHARACTER(LEN=default_string_length)           :: fe_keyword
     REAL(KIND=dp)                                  :: e_offset !WTP_4.2

     LOGICAL, dimension(:), POINTER                 :: active_logical
     TYPE(cp_sll_val_type), POINTER                 :: listcoeff, list_tresp
     !LOGICAL                                        :: continue_while
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma, file_name_s_ppi!,&
                                                       !prop_wf_keyword
     ! OZ: extra strings for C-S and S-S SAB calc
!WTP_merge: 
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_C_S,file_name_s_ppi_C_S
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_S_S,file_name_s_ppi_S_S
     TYPE(section_vals_type), POINTER               :: coupling_section, site_energy_section, ct_block
     TYPE(section_vals_type), POINTER               :: fe_block, fect_block 
     !INTEGER                                        :: aom_atoms, mol
     integer                                        :: atoms_per_donor, atoms_per_acceptor
     !>WTP_comment: maybe we don't need these...
     integer                                        :: atoms_per_fe, atoms_per_ct
     TYPE(section_vals_type), POINTER               :: tresp_coupling_section

     INTEGER, DIMENSION(:), POINTER                 :: donor_index, acceptor_index
     REAL(KIND=dp), dimension(:,:), allocatable     :: Atp_loc
     LOGICAL, DIMENSION(:), ALLOCATABLE             :: submask
     INTEGER                                        :: init_atom_mol_d, final_atom_mol_d, atoms_per_mol_d
     INTEGER                                        :: init_atom_mol_a, final_atom_mol_a, atoms_per_mol_a
     INTEGER                                        :: init_atom_mol_fe, final_atom_mol_fe, atoms_per_mol_fe
     INTEGER                                        :: init_atom_mol_ct, final_atom_mol_ct, atoms_per_mol_ct
     CHARACTER(LEN=default_string_length)           :: STO_input_string
     REAL(KIND=dp)                                  :: STO_rmin,STO_rmax,STO_dr
     INTEGER                                        :: STO_debug
     INTEGER                                        :: atom_index, conn_counter

     NULLIFY(listcoeff)
     NULLIFY(coupling_section)
     NULLIFY(site_energy_section, ct_block)
     NULLIFY(site_energy_section, fe_block)
     NULLIFY(site_energy_section, fect_block)

!     open(324,file="set_aom_env_log.txt",position="append")                                                  

     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")

     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
     !get CT block for donor and acceptor
     ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
     fe_block => section_vals_get_subs_vals(site_energy_section,"FE_BLOCK")
     fect_block => section_vals_get_subs_vals(site_energy_section,"FE_CT_BLOCK")

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

     !WTP_comment:: RANDOM_DERIVATIVE card not exists in the input...
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=sh_env%do_random)

!WTP_merge: no longer requires look-up tables...
     !> Create the overlap look-up table for the AOM coupling calculation
!!     CALL section_vals_val_get(coupling_section,"PSIGMA_FILE_NAME",&
!!             c_val=file_name_s_psigma)
!!     CALL section_vals_val_get(coupling_section,"PPI_FILE_NAME",&
!!             c_val=file_name_s_ppi)
     ! OZ: C-S and S-S
!!     CALL section_vals_val_get(coupling_section,"PSIGMA_C_S_FILE_NAME",&
!!             c_val=file_name_s_psigma_C_S)
!!     CALL section_vals_val_get(coupling_section,"PPI_C_S_FILE_NAME",&
!!             c_val=file_name_s_ppi_C_S)
!!     CALL section_vals_val_get(coupling_section,"PSIGMA_S_S_FILE_NAME",&
!!             c_val=file_name_s_psigma_S_S)
!!     CALL section_vals_val_get(coupling_section,"PPI_S_S_FILE_NAME",&
!!             c_val=file_name_s_ppi_S_S)
     CALL section_vals_val_get(coupling_section,"STO_EXPONENTS",&
             c_val=STO_input_string)
!WTP_merge: STO_EXPONENTS is required input; the rest STO_* have default
     CALL section_vals_val_get(coupling_section,"STO_RMIN",&
             r_val=STO_rmin)
     CALL section_vals_val_get(coupling_section,"STO_RMAX",&
             r_val=STO_rmax)
     CALL section_vals_val_get(coupling_section,"STO_DR",&
             r_val=STO_dr)
     CALL section_vals_val_get(coupling_section,"STO_DEBUG",&
             i_val=STO_debug)

     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_DONOR",&
             r_val=sh_env%scaling_donor)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_ACCEPTOR",&
             r_val=sh_env%scaling_acceptor)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_FECT",&
             r_val=sh_env%scaling_fect)

     !CALL section_vals_val_get(optimization_section,"RK_PROPAGATION",
     !c_val=prop_wf_keyword) 
     !!> Multiple time step calculation of aom gradient
     !CALL section_vals_val_get(coupling_section,"MULTIPLE_TIME_STEP",&
     !        i_val=sh_env%multiple_time_step)

     ! passing N_DONOR_ATOMS and N_ACCEPTOR_ATOMS is not strictly useful since
     ! this info is already in full_system_info or active_system_info 
     CALL section_vals_val_get(ct_block,"N_DONOR_ATOMS",&
              i_val=atoms_per_donor)
     CALL section_vals_val_get(ct_block,"N_ACCEPTOR_ATOMS",&
              i_val=atoms_per_acceptor)

     CALL section_vals_val_get(ct_block,"INDEX_DONOR_MOLS",&
                 i_vals=donor_index)
     CALL section_vals_val_get(ct_block,"INDEX_ACCEPTOR_MOLS",&
                i_vals=acceptor_index)

     !> X-SH: FE_CT_BLOCK subsection...
     CALL section_vals_val_get(fect_block,"FE_IS", c_val=fe_keyword)

     SELECT CASE(fe_keyword)
     CASE("DONOR")
         sh_env%fe_donor = .TRUE.
         atoms_per_fe = atoms_per_donor
         atoms_per_ct = atoms_per_acceptor
     CASE("ACCEPTOR")
         sh_env%fe_donor = .FALSE.
         atoms_per_fe = atoms_per_acceptor
         atoms_per_ct = atoms_per_donor
     END SELECT
     !WTP_4.2_Debug:
     print *, "Debug0"
     !>WTP_4.2: Energy offset between FE and CT states
!     CALL section_vals_val_get(fect_block,"E_OFFSET", r_val=sh_env%e_offset) !WTP_4.2: double check

     !> READ THE AOM COEFFICIENT of DONOR and ACCEPTOR
!WTP_comment: everything (for X-SH, 4 AOM include files) be read in listcoeff,
!then read by read_listcoeff subroutine. This makes sanity check for input
!difficult. Thus, need to make sure 4 AOM files in the order: 1. DONOR 2.
!ACCEPTOR 3. FE 4 CT
!> WTP_4.2: and 5 6 for CR!!

     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     !Transform them into four arrays. In this case external presence is not
     !read anymore.
     CALL read_listcoeff(listcoeff, ext_elements_d, active_logical, &
                         scoefftot_d, pcoefftot_d, atoms_per_donor) !donor

     DEALLOCATE(active_logical)

!     write(324,*) "Hole coefficients:"
!     write(324,*) pcoefftot_d

     CALL read_listcoeff(listcoeff, ext_elements_a, active_logical, &
                         scoefftot_a, pcoefftot_a, atoms_per_acceptor) !acceptor

!     write(324,*) "Electron coefficients:"
!     write(324,*) pcoefftot_a


     DEALLOCATE(active_logical)
                 
     !> AOM for X-SH
     CALL read_listcoeff(listcoeff, ext_elements_fe, active_logical, &
                         scoefftot_fe, pcoefftot_fe, atoms_per_fe) !FE:: donor/acceptor

!     write(324,*) "Exciton(h) coefficients:"
!     write(324,*) pcoefftot_fe

     DEALLOCATE(active_logical)
                 
     CALL read_listcoeff(listcoeff, ext_elements_ct, active_logical, &
                         scoefftot_ct, pcoefftot_ct, atoms_per_ct) !CT:: donor/acceptor

     DEALLOCATE(active_logical)
 
     !>WTP_4.2: add 2 orbs. for CR here!!
     !>maybe with a if(CR) statement???
     !> ....

     !get the first mol of the donor and the first of the acceptor phase 
     init_atom_mol_d = sh_env%full_system_info(1,donor_index(1)) !just get the first donor
     final_atom_mol_d = sh_env%full_system_info(2,donor_index(1))
     atoms_per_mol_d = (final_atom_mol_d-init_atom_mol_d)+1

     init_atom_mol_a = sh_env%full_system_info(1,acceptor_index(1)) !just get the first acceptor
     final_atom_mol_a = sh_env%full_system_info(2,acceptor_index(1))
     atoms_per_mol_a = (final_atom_mol_a-init_atom_mol_a)+1

     SELECT CASE(fe_keyword)
     CASE("DONOR")
        init_atom_mol_fe = init_atom_mol_d
        final_atom_mol_fe = final_atom_mol_d
        atoms_per_mol_fe = atoms_per_mol_d
        init_atom_mol_ct = init_atom_mol_a
        final_atom_mol_ct = final_atom_mol_a
        atoms_per_mol_ct = atoms_per_mol_a
     CASE("ACCEPTOR")
        init_atom_mol_fe = init_atom_mol_a
        final_atom_mol_fe = final_atom_mol_a
        atoms_per_mol_fe = atoms_per_mol_a
        init_atom_mol_ct = init_atom_mol_d
        final_atom_mol_ct = final_atom_mol_d
        atoms_per_mol_ct = atoms_per_mol_d
     END SELECT

     ! sanity checks
     IF (atoms_per_mol_a .ne. atoms_per_acceptor) THEN
        CPABORT("Error in number of atoms acceptor")
     ENDIF
     IF (atoms_per_mol_d .ne. atoms_per_donor) THEN
        CPABORT("Error in number of atoms donor")
     ENDIF

      IF (norbitals .ne. 1) THEN
        CPABORT("DEGENERACY IS NOT IMPLEMENTED FOR MULTI-PHASE")
      ENDIF

     !norbitals must be 1: we have not implemented degeneracy in the ffollowing
     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected

     !> TRESP relavant 
     NULLIFY(tresp_coupling_section)
     NULLIFY(list_tresp)
     tresp_coupling_section => section_vals_get_subs_vals(adiab_section, "TRESP_COUPLINGS")
     !> READ TRESP CHARGES TEST
     CALL section_vals_list_get(tresp_coupling_section,"_DEFAULT_KEYWORD_",list=list_tresp)!Read list tresp from input
     CALL section_vals_val_get(tresp_coupling_section,"CUT_OFF_BETWEEN_FE_SITES",&
          r_val=sh_env%cutoff_sites_fe)
     !> determine if the TRESP are taken from Gaussian code and divide by
     ! sqrt(2) if that's the case (see Multiwfn.pdf for explanation)

!>WTP_comment: maybe check gaussian_charges... 
!!     CALL section_vals_val_get(tresp_coupling_section,"TRESP_FROM_GAUSSIAN",&
!!          l_val=sh_env%Gaussian_charges)
     !### COMMENT_SG: not general for multiphase: but it could be easily
     !charged by giving different charges for different molecules

     IF (sh_env%fe_donor) THEN
       sh_env%atoms_per_site = atoms_per_donor
       ext_elements = ext_elements_d 
     ELSE
       sh_env%atoms_per_site = atoms_per_acceptor
       ext_elements = ext_elements_a
     END IF

     !> COMMENT_SG: sh_env%tresp_charges are the tresp for a single molecule,
     !this would need to be generalized
                   !for multiple molecule kind so that the coupling
                   !calculation is general for any num atoms per molecule
     CALL read_list_tresp(sh_env%atoms_per_site, list_tresp, ext_elements, active_logical, sh_env%tresp_charges)
     DEALLOCATE(active_logical)

     !> active_logical is not necessary since inactive atoms will anyway have
     !0
     !charges that will nullify the term in the coulomb sum 

     !c = nadiab*sh_env%atoms_per_site !> Numbers of active atoms per mol 

!WTP_comment: comment out NN/grad_tresp. deal with them differently in X-SH 
!!     allocate(sh_env%crds(3,COUNT(sh_env%ext_pres),3))  !> crds active atoms
!!     ALLOCATE(sh_env%NN(mat, mat))
!!     ALLOCATE(sh_env%grad_tresp(3,mat, COUNT(sh_env%ext_pres)))

     !WTP_comment: nacv_fe=<CT|grad_I XT>, nacv_ct=<XT|grad_I CT> 
     IF (sh_env%fe_donor) THEN
       ALLOCATE(sh_env%NN(SIZE(donor_index),SIZE(donor_index)))
       ALLOCATE(sh_env%grad_tresp(3,SIZE(donor_index), SIZE(donor_index)*atoms_per_donor))
       ALLOCATE(sh_env%nacv_fe(3,SIZE(acceptor_index), SIZE(donor_index)*atoms_per_donor ))
       ALLOCATE(sh_env%nacv_ct(3,SIZE(donor_index), SIZE(acceptor_index)*atoms_per_acceptor ))
     ELSE
       ALLOCATE(sh_env%NN(SIZE(acceptor_index),SIZE(acceptor_index)))
       ALLOCATE(sh_env%grad_tresp(3,SIZE(acceptor_index),SIZE(acceptor_index)*atoms_per_acceptor))
       ALLOCATE(sh_env%nacv_fe(3,SIZE(donor_index), SIZE(acceptor_index)*atoms_per_acceptor ))
       ALLOCATE(sh_env%nacv_ct(3,SIZE(acceptor_index), SIZE(donor_index)*atoms_per_donor ))
     END IF

     ALLOCATE(sh_env%nacv_dphase(3,SIZE(donor_index)*SIZE(acceptor_index), &
                                             SIZE(donor_index)*SIZE(acceptor_index)*atoms_per_mol_d))
     ALLOCATE(sh_env%nacv_aphase(3,SIZE(donor_index)*SIZE(acceptor_index), &
                                             SIZE(donor_index)*SIZE(acceptor_index)*atoms_per_mol_a))
     sh_env%nacv_dphase = 0.0_dp
     sh_env%nacv_aphase = 0.0_dp
     sh_env%nacv_fe = 0.0_dp
     sh_env%nacv_ct = 0.0_dp
     sh_env%grad_tresp = 0.0_dp 

     !!> crds is an AOM arrays (xyz, atoms, tm-t-tp) the third dimension
     !!indicates what time the coordinates refer to
     !allocate(sh_env%crds(3,act_atoms,3))

!WTP_comment: NN seems not used in the CT part, so..
!WTP_comment: --> NN comment out from CT part and only allocate for FE part...
!!     ALLOCATE(sh_env%NN(SIZE(sh_env%index_active_mols), SIZE(sh_env%index_active_mols)))

     allocate(sh_env%siteenergies(3,mat))  !COMMENT_SG: First dimension is not usefull
     sh_env%siteenergies = 0.0_dp

     ! OZ: add extra arguments...
!!     CALL build_s_tables(file_name_s_psigma, file_name_s_ppi,&
!!                         file_name_s_psigma_C_S, file_name_s_ppi_C_S,&
!!                         file_name_s_psigma_S_S, file_name_s_ppi_S_S)
!!WTP_merge: now use spline... 
     ! OZ: add extra arguments...
     CALL STO_initialize_splines(STO_input_string,STO_rmin,STO_rmax,STO_dr,&
                                 STO_debug)

     !> READ THE AOM COEFFICIENT DONOR
     !> CONSTRUCT THE CONNLIST 
     ! SG: strange way to proceed to get the connectivity of one molecule (which
     ! is then assumed to be the same for all molecules) 
     !we fill out the connlist arrays before calculating the nacv
     allocate(submask(atoms_per_donor))   !temp array
     submask =.NOT.((scoefftot_d.eq.0).AND.(pcoefftot_d.eq.0)) !Only select atoms with orbitals

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_d(COUNT(submask(1:atoms_per_donor))))
     sh_env%atomlist_d = PACK( (/ (i, i=1,atoms_per_donor) /), submask(1:atoms_per_donor) )
     allocate(sh_env%connlist2_d(4,size(sh_env%atomlist_d)))
     ALLOCATE(sh_env%connlist1_d(4,atoms_per_donor))

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2 for agiven phase
     ALLOCATE(Atp_loc(4,atoms_per_donor)) !temp array
     Atp_loc(1,1:atoms_per_donor) = dble(ext_elements_d)
     Atp_loc(2:4,1:atoms_per_donor) = positions_tp(1:3,init_atom_mol_d:final_atom_mol_d)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     !WTP_4.2_Debug:
!     print *, "WTP_Debug1: "
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_d, sh_env%connlist1_d, sh_env%connlist2_d, sh_env%cutoff_connect)
!     print *, "WTP_Debug2: "
!     print *, "size(connlist1(1,:)) and size(connlist1(:,1))", size(sh_env%connlist1_d(1,:)), size(sh_env%connlist1_d(:,1))
!     print *, "size(connlist2(1,:)) and size(connlist2(:,1))", size(sh_env%connlist2_d(1,:)), size(sh_env%connlist2_d(:,1))
!     print *, "connlist1: ", sh_env%connlist1_d
!     print *, "connlist2: ", sh_env%connlist2_d 

     open(422, file="donor-connectivity.xyz", position="append")
     write(422,*) "Printing sulphur atoms' nearest neighbours:"

     do conn_counter = 1,size(sh_env%connlist2_d, dim=2)
        atom_index = sh_env%connlist2_d(1,conn_counter)
        if (Atp_loc(1,atom_index) == 16) then
               write(422,*) sh_env%connlist2_d(:,conn_counter)
        end if
     end do

!     close(422)

     ! coefficients donor phase
     ALLOCATE(sh_env%pcoeffs_d(atoms_per_donor))
     ALLOCATE (sh_env%scoeffs_d(atoms_per_donor))
     ALLOCATE (sh_env%ext_elements_d(atoms_per_donor))
     sh_env%pcoeffs_d = pcoefftot_d
     sh_env%scoeffs_d = scoefftot_d
     sh_env%ext_elements_d = ext_elements_d
     DEALLOCATE(ext_elements_d)
     DEALLOCATE(scoefftot_d)
     DEALLOCATE(pcoefftot_d)
     DEALLOCATE(submask)
     DEALLOCATE(Atp_loc)

     !> READ THE AOM COEFFICIENT ACCEPTOR PHASE! 
     !> CONSTRUCT THE CONNLIST 
     allocate(submask(atoms_per_acceptor))   !temp array
     submask =.NOT.((scoefftot_a.eq.0).AND.(pcoefftot_a.eq.0)) !Only select atoms with orbitals

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_a(COUNT(submask(1:atoms_per_acceptor))))
     sh_env%atomlist_a = PACK( (/ (i, i=1,atoms_per_acceptor) /), submask(1:atoms_per_acceptor) )
     allocate(sh_env%connlist2_a(4,size(sh_env%atomlist_a)))
     ALLOCATE(sh_env%connlist1_a(4,atoms_per_acceptor))

!     write(324,*) "AOM_active atoms acceptor:"
!     write(324,*) sh_env%atomlist_a

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2
     ALLOCATE(Atp_loc(4,atoms_per_acceptor)) !temp array
     Atp_loc(1,1:atoms_per_acceptor) = dble(ext_elements_a)
     Atp_loc(2:4,1:atoms_per_acceptor) = positions_tp(1:3,init_atom_mol_a:final_atom_mol_a)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_a, sh_env%connlist1_a, sh_env%connlist2_a, sh_env%cutoff_connect)

     open(333, file="2phase_pdi_connlist.txt", position="append")
     write(333,*) sh_env%connlist2_a
     close(333)

     open(433, file="acceptor-connectivity.xyz", position="append")
     write(433,*) "Print oxygen atoms' nearest neighbours:"

     do conn_counter = 1,size(sh_env%connlist2_a, dim=2)
        atom_index = sh_env%connlist2_a(1,conn_counter)
        if (Atp_loc(1,atom_index) == 8) then
               write(433,*) sh_env%connlist2_a(:,conn_counter)
        end if
     end do

     close(433)


     ! coefficients acceptor phase
     ALLOCATE(sh_env%pcoeffs_a(atoms_per_acceptor))
     ALLOCATE (sh_env%scoeffs_a(atoms_per_acceptor))
     ALLOCATE (sh_env%ext_elements_a(atoms_per_acceptor))
     sh_env%pcoeffs_a = pcoefftot_a
     sh_env%scoeffs_a = scoefftot_a
     sh_env%ext_elements_a = ext_elements_a
     DEALLOCATE(ext_elements_a)
     DEALLOCATE(submask)
     DEALLOCATE(Atp_loc)
     DEALLOCATE(scoefftot_a)
     DEALLOCATE(pcoefftot_a)

     !> READ THE AOM COEFFICIENT ACCEPTOR PHASE! 
     !> CONSTRUCT THE CONNLIST 
     allocate(submask(atoms_per_fe))   !temp array
     submask =.NOT.((scoefftot_fe.eq.0).AND.(pcoefftot_fe.eq.0)) !Only select atoms with orbitals

     !WRITE(*,*) 'submask: ', submask

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_fe(COUNT(submask(1:atoms_per_fe))))
     sh_env%atomlist_fe = PACK( (/ (i, i=1,atoms_per_fe) /), submask(1:atoms_per_fe) )
     allocate(sh_env%connlist2_fe(4,size(sh_env%atomlist_fe)))
     ALLOCATE(sh_env%connlist1_fe(4,atoms_per_fe))

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2
     ALLOCATE(Atp_loc(4,atoms_per_fe)) !temp array
     Atp_loc(1,1:atoms_per_fe) = dble(ext_elements_fe)
!init_atom_mol_fe, final_atom_mol_fe, ... need to be added or modified...
     Atp_loc(2:4,1:atoms_per_fe) = positions_tp(1:3,init_atom_mol_fe:final_atom_mol_fe)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_fe,sh_env%connlist1_fe, sh_env%connlist2_fe, sh_env%cutoff_connect)

     !open(221,file="fe_atomlist.txt",position="append")
     !write(221,*) sh_env%atomlist_fe
     !close(221)

     !WTP_TEST
!     WRITE(*,*) 'connlist1_fe: ', sh_env%connlist1_fe
!     WRITE(*,*) 'connlist1_fe: ', sh_env%connlist1_fe(1,:)
!     WRITE(*,*) 'connlist2_fe: ', sh_env%connlist2_fe

     ! coefficients acceptor phase
     ALLOCATE(sh_env%pcoeffs_fe(atoms_per_fe))
     ALLOCATE (sh_env%scoeffs_fe(atoms_per_fe))
     ALLOCATE (sh_env%ext_elements_fe(atoms_per_fe))
     sh_env%pcoeffs_fe = pcoefftot_fe
     sh_env%scoeffs_fe = scoefftot_fe
     sh_env%ext_elements_fe = ext_elements_fe
     DEALLOCATE(ext_elements_fe)
     DEALLOCATE(submask)
     DEALLOCATE(Atp_loc)
     DEALLOCATE(scoefftot_fe)
     DEALLOCATE(pcoefftot_fe)

     !WTP_comment: then CT part
     allocate(submask(atoms_per_ct))   !temp array
     submask =.NOT.((scoefftot_ct.eq.0).AND.(pcoefftot_ct.eq.0)) !Only select atoms with orbitals

     !create atome of interest (aoi) list
     !List atoms with orbitals namely, with coeff not zero
     allocate(sh_env%atomlist_ct(COUNT(submask(1:atoms_per_ct))))
     sh_env%atomlist_ct = PACK( (/ (i, i=1,atoms_per_ct) /), submask(1:atoms_per_ct) )
     allocate(sh_env%connlist2_ct(4,size(sh_env%atomlist_ct)))
     ALLOCATE(sh_env%connlist1_ct(4,atoms_per_ct))

     !open(223, file="ct_atomlist.txt", position="append")
     !write(223,*) sh_env%atomlist_ct
     !close(223)

     ! SG: this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2
     ALLOCATE(Atp_loc(4,atoms_per_ct)) !temp array
     Atp_loc(1,1:atoms_per_ct) = dble(ext_elements_ct)
!init_atom_mol_ct, final_atom_mol_ct, ... need to be added or modified...
     Atp_loc(2:4,1:atoms_per_ct) = positions_tp(1:3,init_atom_mol_ct:final_atom_mol_ct)
     !output: sh_env%connlist1, sh_env%connlist2, which are the intramoleculat
     !connectivities with and without considering the zeors
     call connect_list2(Atp_loc(2:4,:), sh_env%atomlist_ct, sh_env%connlist1_ct, sh_env%connlist2_ct, sh_env%cutoff_connect)

     !WTP_TEST:
!     WRITE(*,*) 'connlist1_ct: ', sh_env%connlist1_ct
!     WRITE(*,*) 'connlist1_ct: ', sh_env%connlist1_ct(1,:)
!     WRITE(*,*) 'connlist2_ct: ', sh_env%connlist2_ct

!     write(324,*) "Conn list:"
!     write(324,*) sh_env%connlist2_a


     ! coefficients acceptor phase
     ALLOCATE(sh_env%pcoeffs_ct(atoms_per_ct))
     ALLOCATE (sh_env%scoeffs_ct(atoms_per_ct))
     ALLOCATE (sh_env%ext_elements_ct(atoms_per_ct))
     sh_env%pcoeffs_ct = pcoefftot_ct
     sh_env%scoeffs_ct = scoefftot_ct
     sh_env%ext_elements_ct = ext_elements_ct
     DEALLOCATE(ext_elements_ct)
     DEALLOCATE(submask)
     DEALLOCATE(Atp_loc)
     DEALLOCATE(scoefftot_ct)
     DEALLOCATE(pcoefftot_ct)

  END SUBROUTINE set_aom_env_fect

!*************************************************************************
!Subroutine for setting up the sh_env type for the rest of the dynamics
!by passing in parameters from the input file, such as AOM and TRESP,
!etc.

!INPUTS: nadiab: int, total number of states ; norbitals: No. of
!degenerate orbitals per state (ONLY WORKS FOR 1 CURRENTLY) ; natom:
!int, total No. of atoms in whole system ; positions_tp: array,real, 3D
!coords of all atoms ; adiab_section: pointer to part of input file ;
!sh_env: type where we pass key info from input file ; my_mol_index: 1D
!array of ints, DECOMP indices of active mols

!No specific output, just passing important info from all args into the
!sh_env type so we can use them across different subroutines
!*************************************************************************
  SUBROUTINE set_aom_env_sp_fect(nadiab, norbitals, natom, positions_tp, &
  adiab_section, sh_env, my_mol_index)

     INTEGER, INTENT(IN)                            :: nadiab, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(sh_env_type), INTENT(INOUT)               :: sh_env
     INTEGER, DIMENSION(:), POINTER                 :: my_mol_index

     CHARACTER(len=*), PARAMETER :: routineN ='set_aom_env_sp_fect', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_h, scoefftot_h, & !aom s/p slater expansion coeffs
                                                       pcoefftot_e, scoefftot_e
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_fe_h, scoefftot_fe_h, &
                                                       pcoefftot_ct_h, scoefftot_ct_h !AOM coeff files for XT-CT coupling via HOMO-HOMO overlap
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot_fe_e, scoefftot_fe_e, &
                                                       pcoefftot_ct_e, scoefftot_ct_e !AOM coeff files for XT-CT via LUMO-LUMO overlap (e-transfer)         
     integer                                        :: mat, i
     integer, dimension(:), ALLOCATABLE             :: ext_elements !molecule elements (atomic No.)

     LOGICAL, dimension(:), POINTER                 :: active_logical !true for aom-active atoms in molecule
     TYPE(cp_sll_val_type), POINTER                 :: listcoeff, list_tresp
     TYPE(section_vals_type), POINTER               :: coupling_section, site_energy_section, ct_block
     TYPE(section_vals_type), POINTER               :: fe_block, fect_block 
     TYPE(section_vals_type), POINTER               :: tresp_coupling_section

     REAL(KIND=dp), dimension(:,:), allocatable     :: Atp_loc
     LOGICAL, DIMENSION(:), ALLOCATABLE             :: submask
     INTEGER                                        :: initial_atom_mol1, final_atom_mol1, atoms_per_mol, parsed_atoms_per_mol
     !only one set of atom counters needed here as it's a single-phase system
     CHARACTER(LEN=default_string_length)           :: STO_input_string
     REAL(KIND=dp)                                  :: STO_rmin,STO_rmax,STO_dr
     INTEGER                                        :: STO_debug
     INTEGER                                        :: atom_index, conn_counter

     NULLIFY(listcoeff)
     NULLIFY(coupling_section)
     NULLIFY(site_energy_section, ct_block)
     NULLIFY(site_energy_section, fe_block)
     NULLIFY(site_energy_section, fect_block)

!     open(324, file="set_aom_env_log.txt", position="append")

     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")

     site_energy_section => section_vals_get_subs_vals(adiab_section,"METHOD_SITE_ENERGIES")
     !get CT block for donor and acceptor
     ct_block => section_vals_get_subs_vals(site_energy_section,"CT_BLOCK")
     fe_block => section_vals_get_subs_vals(site_energy_section,"FE_BLOCK")

     !> Initialize input parameters: 
     !> Get number of atoms per molecule from input file keyword
     CALL section_vals_val_get(coupling_section,"NUMBER_AOM_ATOMS_PER_SITE",&
             i_val=atoms_per_mol)
     sh_env%atoms_per_site = atoms_per_mol

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

     CALL section_vals_val_get(coupling_section,"STO_EXPONENTS",&
             c_val=STO_input_string)
!WTP_merge: STO_EXPONENTS is required input; the rest STO_* have default
     CALL section_vals_val_get(coupling_section,"STO_RMIN",&
             r_val=STO_rmin)
     CALL section_vals_val_get(coupling_section,"STO_RMAX",&
             r_val=STO_rmax)
     CALL section_vals_val_get(coupling_section,"STO_DR",&
             r_val=STO_dr)
     CALL section_vals_val_get(coupling_section,"STO_DEBUG",&
             i_val=STO_debug)

     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_h-h",&
             r_val=sh_env%scaling_hh)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_e-e",&
             r_val=sh_env%scaling_ee)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_fect_etransfer",&
             r_val=sh_env%scaling_fect_e)
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR_fect_htransfer",&
             r_val=sh_env%scaling_fect_h)

     !now reading in lists of AOM coefficients, the order in which
     !they're passed into the input file needs to match the order in
     !which they're read in

     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     !Transform them into four arrays. In this case external presence is not
     !read anymore.

     !h-h and e-e coupling

     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                         scoefftot_h, pcoefftot_h, sh_env%atoms_per_site)

     ALLOCATE(sh_env%pcoeffs_hh(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_hh(sh_env%atoms_per_site))

     sh_env%pcoeffs_hh = pcoefftot_h
     sh_env%scoeffs_hh = scoefftot_h

     ALLOCATE(sh_env%elements(size(ext_elements)))
     sh_env%elements_sp = ext_elements
     !conjugated atom list taken from HOMO

!     write(324,*) "Hole coefficients:"
!     write(324,*) sh_env%pcoeffs_hh

     DEALLOCATE(pcoefftot_h)
     DEALLOCATE(scoefftot_h)
     DEALLOCATE(active_logical)
     DEALLOCATE(ext_elements)
     !deallocating active_logical and ext_elements every time, as
     !read_listcoeff allocates them every call

     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                         scoefftot_e, pcoefftot_e, sh_env%atoms_per_site)
     
     ALLOCATE(sh_env%pcoeffs_ee(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_ee(sh_env%atoms_per_site))

     sh_env%pcoeffs_ee = pcoefftot_e
     sh_env%scoeffs_ee = scoefftot_e

!     write(324,*) "Electron coefficients:"
!     write(324,*) sh_env%pcoeffs_ee

     DEALLOCATE(pcoefftot_e)
     DEALLOCATE(scoefftot_e)
     DEALLOCATE(ext_elements)
     DEALLOCATE(active_logical)
                 
     !XT-CT coupling with HOMOS of CT and XT (h-transfer)
     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                         scoefftot_fe_h, pcoefftot_fe_h, sh_env%atoms_per_site)

     ALLOCATE(sh_env%pcoeffs_fe_h(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_fe_h(sh_env%atoms_per_site))

     sh_env%pcoeffs_fe_h = pcoefftot_fe_h
     sh_env%scoeffs_fe_h = scoefftot_fe_h

!     write(324,*) "Exciton(h) coefficients:"
!     write(324,*) sh_env%pcoeffs_fe_h

     DEALLOCATE(active_logical)
     DEALLOCATE(ext_elements)
                 
     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                 scoefftot_ct_h, pcoefftot_ct_h, sh_env%atoms_per_site)

     ALLOCATE(sh_env%pcoeffs_ct_h(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_ct_h(sh_env%atoms_per_site))

     sh_env%pcoeffs_ct_h = pcoefftot_ct_h
     sh_env%scoeffs_ct_h = scoefftot_ct_h

     DEALLOCATE(pcoefftot_fe_h)
     DEALLOCATE(scoefftot_fe_h)
     DEALLOCATE(pcoefftot_ct_h)
     DEALLOCATE(scoefftot_ct_h) 
     DEALLOCATE(active_logical)
     DEALLOCATE(ext_elements)
 
     !XT-CT coupling with LUMOS of CT and XT (e-transfer)
     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                         scoefftot_fe_e, pcoefftot_fe_e, sh_env%atoms_per_site)

     ALLOCATE(sh_env%pcoeffs_fe_e(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_fe_e(sh_env%atoms_per_site))

     sh_env%pcoeffs_fe_e = pcoefftot_fe_e
     sh_env%scoeffs_fe_e = scoefftot_fe_e

     DEALLOCATE(active_logical)
     DEALLOCATE(ext_elements)
                 
     CALL read_listcoeff(listcoeff, ext_elements, active_logical, &
                         scoefftot_ct_e, pcoefftot_ct_e, sh_env%atoms_per_site)

     ALLOCATE(sh_env%pcoeffs_ct_e(sh_env%atoms_per_site))
     ALLOCATE(sh_env%scoeffs_ct_e(sh_env%atoms_per_site))

     sh_env%pcoeffs_ct_e = pcoefftot_ct_e
     sh_env%scoeffs_ct_e = scoefftot_ct_e

     DEALLOCATE(pcoefftot_fe_e)
     DEALLOCATE(scoefftot_fe_e)
     DEALLOCATE(pcoefftot_ct_e)
     DEALLOCATE(scoefftot_ct_e)
     DEALLOCATE(active_logical)
     DEALLOCATE(ext_elements)

     !check molecule parsing
     initial_atom_mol1 = sh_env%full_system_info(1, my_mol_index(1))
     final_atom_mol1 = sh_env%full_system_info(2, my_mol_index(1))
     parsed_atoms_per_mol = final_atom_mol1 - initial_atom_mol1 + 1

     IF (atoms_per_mol .ne. parsed_atoms_per_mol) THEN
        CPABORT("NO. ATOMS PER MOLECULE IN INPUT FILE DOESN'T MATCH MOLECULAR STRUCTURE")
     ENDIF

     IF (norbitals .ne. 1) THEN
        CPABORT("Orbital degeneracy not yet implemented in X-SH")
     ENDIF

     !norbitals must be 1: we have not implemented degeneracy
     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected

     !> TRESP relavant 
     NULLIFY(tresp_coupling_section)
     NULLIFY(list_tresp)
     tresp_coupling_section => section_vals_get_subs_vals(adiab_section, "TRESP_COUPLINGS")
     !> READ TRESP CHARGES TEST
     CALL section_vals_list_get(tresp_coupling_section,"_DEFAULT_KEYWORD_",list=list_tresp)!Read list tresp from input
     CALL section_vals_val_get(tresp_coupling_section,"CUT_OFF_BETWEEN_FE_SITES",&
          r_val=sh_env%cutoff_sites_fe)
     !> determine if the TRESP are taken from Gaussian code and divide by
     ! sqrt(2) if that's the case (see Multiwfn.pdf for explanation)

     CALL read_list_tresp(sh_env%atoms_per_site, list_tresp, ext_elements, active_logical, sh_env%tresp_charges)
     DEALLOCATE(active_logical)

     !allocating nacv storage arrays, every atom can feel force from all
     !diabat types, so conditional not needed

     !nearest-neighbour list for tresp + nacv arrays for e-e and h-h coupling
     ALLOCATE(sh_env%NN(size(my_mol_index), size(my_mol_index)))

     ALLOCATE(sh_env%nacv_hh(3, size(my_mol_index), &
size(my_mol_index)*sh_env%atoms_per_site))

     ALLOCATE(sh_env%nacv_ee(3, size(my_mol_index), &
size(my_mol_index)*sh_env%atoms_per_site))

     !nacv arrays for xt-ct coupling via h-transfer 
     ALLOCATE(sh_env%nacv_fe_h(3, size(my_mol_index), sh_env%atoms_per_site*size(my_mol_index)))
     ALLOCATE(sh_env%nacv_ct_h(3, size(my_mol_index), sh_env%atoms_per_site*size(my_mol_index)))

     !nacv arrays for xt-ct coupling via e-transfer 
     ALLOCATE(sh_env%nacv_fe_e(3, size(my_mol_index), sh_env%atoms_per_site*size(my_mol_index)))
     ALLOCATE(sh_env%nacv_ct_e(3, size(my_mol_index), sh_env%atoms_per_site*size(my_mol_index)))

     !nacv tresp array
     ALLOCATE(sh_env%grad_tresp(3, size(my_mol_index), &
size(my_mol_index)*sh_env%atoms_per_site))

     sh_env%nacv_ee = 0.0_dp
     sh_env%nacv_hh = 0.0_dp
     sh_env%grad_tresp = 0.0_dp
     sh_env%nacv_fe_h = 0.0_dp
     sh_env%nacv_ct_h = 0.0_dp
     sh_env%nacv_fe_e = 0.0_dp
     sh_env%nacv_ct_e = 0.0_dp

     ALLOCATE(sh_env%siteenergies(3,mat))
     sh_env%siteenergies = 0.0_dp

     !passing in parameters needed for Mulliken formulas in AOM
     CALL STO_initialize_splines(STO_input_string,STO_rmin,STO_rmax,STO_dr,&
                                 STO_debug)

     !> construct neighbour list for each atom
     ALLOCATE(submask(sh_env%atoms_per_site))   !temp array
     submask =.NOT.( (sh_env%scoeffs_hh.eq.0) .AND. (sh_env%pcoeffs_hh.eq.0) ) !Only select atoms with orbitals

     ALLOCATE(sh_env%atomlist(COUNT(submask(1:sh_env%atoms_per_site))))
     !counting the number of True elements in the submask

     sh_env%atomlist = PACK( (/ (i, i=1,sh_env%atoms_per_site) /), submask(1:sh_env%atoms_per_site) )
     !PACK takes indices of atoms that =TRUE in the mask, so only
     !indices of AOM-active atoms

!     write(324,*) "AOM-active atoms:"
!     write(324,*) sh_env%atomlist

     ALLOCATE(sh_env%connlist2(4,size(sh_env%atomlist)))
     ALLOCATE(sh_env%connlist1(4,sh_env%atoms_per_site))

     !passing coordinates of AOM-active atoms into Atp_loc
     ALLOCATE(Atp_loc(4,sh_env%atoms_per_site))
     Atp_loc(1,1:sh_env%atoms_per_site) = dble(sh_env%elements_sp)
     Atp_loc(2:4,1:sh_env%atoms_per_site) = positions_tp(1:3,initial_atom_mol1:final_atom_mol1)

     CALL connect_list2(Atp_loc(2:4,:), sh_env%atomlist, sh_env%connlist1, sh_env%connlist2, sh_env%cutoff_connect)

!     write(324,*) "Conn list:"
!     write(324,*) sh_env%connlist2


     DEALLOCATE(submask)
     DEALLOCATE(Atp_loc)
     DEALLOCATE(ext_elements)

!     close(324)


  END SUBROUTINE set_aom_env_sp_fect

! *****************************************************************************
!> \brief This allocates aom_env variables. One could in principle put these
!          variables in a separate derive type structure
!> \author Edited SG from JS implementation
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
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_C_S,file_name_s_ppi_C_S
!     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_S_S,file_name_s_ppi_S_S
     TYPE(section_vals_type), POINTER               :: coupling_section
     INTEGER                                        :: aom_atoms
     CHARACTER(LEN=default_string_length)           :: STO_input_string
     REAL(KIND=dp)                                  :: STO_rmin,STO_rmax,STO_dr
     INTEGER                                        :: STO_debug

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
     !CALL section_vals_val_get(coupling_section,"PSIGMA_FILE_NAME",&
     !        c_val=file_name_s_psigma)
     !CALL section_vals_val_get(coupling_section,"PPI_FILE_NAME",&
     !        c_val=file_name_s_ppi)
     ! OZ: C-S and S-S
     !CALL section_vals_val_get(coupling_section,"PSIGMA_C_S_FILE_NAME",&
     !        c_val=file_name_s_psigma_C_S)
     !CALL section_vals_val_get(coupling_section,"PPI_C_S_FILE_NAME",&
     !        c_val=file_name_s_ppi_C_S)
     !CALL section_vals_val_get(coupling_section,"PSIGMA_S_S_FILE_NAME",&
     !        c_val=file_name_s_psigma_S_S)
     !CALL section_vals_val_get(coupling_section,"PPI_S_S_FILE_NAME",&
     !        c_val=file_name_s_ppi_S_S)
     CALL section_vals_val_get(coupling_section,"STO_EXPONENTS",&
             c_val=STO_input_string)
     CALL section_vals_val_get(coupling_section,"STO_RMIN",&
             r_val=STO_rmin)
     CALL section_vals_val_get(coupling_section,"STO_RMAX",&
             r_val=STO_rmax)
     CALL section_vals_val_get(coupling_section,"STO_DR",&
             r_val=STO_dr)
     CALL section_vals_val_get(coupling_section,"STO_DEBUG",&
             i_val=STO_debug)
     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=sh_env%scaling)
     !> cbar factor : Sab = cbar*Sab_bar. cbar is only used for diabatic
     !non-orthogonal nacv calculation (this nacv is negligible) 
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
     allocate(sh_env%siteenergies(3,mat))  !COMMENT_SG: First dimension is not usefull
     ALLOCATE(sh_env%nacv(3,mat, mat*aom_atoms))
     !ALLOCATE(sh_env%nacv_tot(3,mat, mat, c))    !TEST 
     !sh_env%nacv_tot = 0.0d0
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
     CALL STO_initialize_splines(STO_input_string,STO_rmin,STO_rmax,STO_dr,&
                                 STO_debug)
     !CALL build_s_tables(file_name_s_psigma, file_name_s_ppi,&
     !                    file_name_s_psigma_C_S, file_name_s_ppi_C_S,&
     !                    file_name_s_psigma_S_S, file_name_s_ppi_S_S)
      
     
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

     ! coordinate of the active atoms taken from positions_tp
     do coord=1,3
        crds_loc(coord,:) = PACK(positions_tp(coord,:), sh_env%ext_pres)
     enddo

     ! this line again assume same number of atoms for each molecule and it
     ! assumes all molecules have the same connlist1 connlist2
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
