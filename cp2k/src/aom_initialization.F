!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

!!!!!!!!!! CHANGE_AC !!!!!!!!!!!!!

! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2016.01.15: creation, put only ffsh_aom subroutines
! *****************************************************************************

! *****************************************************************************
!>      jacob: comment
! construct_ham_X_SH_full_degen
!
! This program is designed to propagate a single excess electronic charge
! along a molecular dynamics trajectory, either during the trajectory
! (to do surface hopping simulations) or afterwards (classical path approximation).
!
! This program takes molecular coordinates as inputs.
! From these coordinate inputs, it calculates the electronic Hamiltonian
! in a molecular orbital basis, using the overlapFINAL module to calculate
! the off-diagonals and either accepting the diagonals as inputs or calculating
! them in some other way.
!
! With the electronic Hamiltonian, the program then propagates the electron using the RK4 method.
! At the final electronic state, the probabilities for hopping between adiabatic surfaces are calculated.
! The program then decides what adiabatic surface the system is on and calculates
! the appropriate atomic forces.
!
! In the classical path approximation, the above steps are repeated until the
! end of the input molecular trajectory.  For surface-hopping, the data here is fed back into
! the molecular dynamics and the trajectory moves forward another time step.
! *****************************************************************************

MODULE aom_initialization 

     USE kinds,                        ONLY: default_string_length,&
                                             dp
     USE aom_overlapfinal,             ONLY:  connect_list2, &
                                              !calc_sab, &
                                              !sab_screen_new, &
                                              !calc_dRSab, &
                                              !calc_rvecs, &
                                              !calc_sab_pete, &
                                              build_s_tables

     USE aom_utils,                    ONLY: makeX,                                &
                                             adiabaticity,                         &
                                             calc_hop_nacv_aom,                    &
                                             calc_aom_Xmatrix_elements_first_step, & 
                                             make_adiabX,                          &
                                             !surfacehop,     &
                                             calcforces_aom,                           &
                                             read_listcoeff,                       &
                                             read_listdiabcoeff,                   &
                                             calc_coms_and_nlist_aom

     USE sh_types,                    ONLY: sh_env_type, &
                                                sh_print_env_type

     USE mixed_energy_types,           ONLY: mixed_force_type
     USE input_section_types,          ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_list_get,&
                                              section_vals_val_get

     USE cp_linked_list_val,           ONLY: &
         cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
         cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
         cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
         cp_sll_val_set_el_at, cp_sll_val_type

     USE aom_print,                    ONLY: adiab_print_coeff,              &
                                             adiab_print_pop,                &
                                             adiab_print_adiab_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_site_energies,   &
                                             adiab_print_off_diagonals,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian, &
                                             adiab_print_pseudo_hamiltonian,        &
                                             adiab_print_nacv, &
                                             adiab_print_pvecs

     USE atomic_kind_list_types,       ONLY: atomic_kind_list_type

     USE atomic_kind_types,            ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

     USE distribution_1d_types,        ONLY: distribution_1d_type

     USE non_adiabatic_utils,          ONLY: decohere_wf,        &
                                             choose_firstadiabat, &
                                             print_matrix_standard_out,&
                                             diab_to_adiab,      &
                                             translate_com,      &
                                             propagate_wf,       &
                                             propagate_adiab_wf, &
                                             prop_wf_HS_corr!, & 
                                             !classically_forbidden_hops

     USE parallel_rng_types,           ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             rng_stream_type
 
     USE cp_log_handling,              ONLY: cp_get_default_logger, &
                                             cp_logger_type

     USE cp_output_handling,           ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_should_output
     USE mixed_environment_types,      ONLY: mixed_environment_type


#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_initialization'   

     PUBLIC :: create_sh_env, &
               propagate_aom_namd_first_step

CONTAINS

  SUBROUTINE create_sh_env(nadiab, aom_atoms, norbitals, natom, &
                                  positions_tp,     &
                                  adiab_section,    &
                                  sh_env)


     INTEGER, INTENT(IN)                            :: nadiab,  aom_atoms, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: adiab_section, optimization_section
     TYPE(section_vals_type), POINTER               :: wf_restart_section, decoherence_section
     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_sh_env', &
     routineP = moduleN//':'//routineN

     integer                                        :: i, j, mat, c, coord, mol, atoms, start
     REAL(KIND=dp)                                  :: cutoff_connect
     REAL(KIND=dp), dimension(:,:), allocatable     :: crds
     REAL(KIND=dp), dimension(:,:), allocatable     :: Atp
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: pcoefftot, scoefftot
     LOGICAL, DIMENSION(:), ALLOCATABLE             :: atomlist_mask, submask
     integer, dimension(:), ALLOCATABLE             :: ext_elements
     TYPE(cp_sll_val_type), POINTER                 :: listcoeff
     TYPE(cp_sll_val_type), POINTER                 :: listdiabcoeff
!     LOGICAL                                        :: do_random
     REAL(KIND=dp), DIMENSION(3, 2)                 :: initial_seed
     INTEGER                                        :: input_seed
     TYPE(rng_stream_type), POINTER                 :: rng_stream
     LOGICAL                                        :: continue_while
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma, file_name_s_ppi, deco_keyword,&
                                                       prop_wf_keyword
     ! OZ: extra strings for C-S and S-S SAB calc
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_C_S,file_name_s_ppi_C_S
     CHARACTER(LEN=default_string_length)           :: file_name_s_psigma_S_S,file_name_s_ppi_S_S
     TYPE(section_vals_type), POINTER               :: coupling_section
     CHARACTER(LEN=default_string_length)           :: coupling_keyword

     NULLIFY(rng_stream)
     NULLIFY(listcoeff)
     NULLIFY(decoherence_section)
     NULLIFY(wf_restart_section)
     NULLIFY(optimization_section)
     NULLIFY(coupling_section)


     mol = nadiab ! CHANGE_AC
     atoms = aom_atoms ! CHANGE_AC

     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")

     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword) 

     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     c = nadiab*aom_atoms !> Numbers of active atoms


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
   
     CALL create_adiab_nace_env(mol, sh_env)

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
     !ALLOCATE(sh_env%class_forbidden_hop_t(mat))
     allocate(sh_env%adiab_nace_t(mat,mat))
     allocate(sh_env%adiab_nace_tp(mat,mat))

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


!=================== AOM array ===================================
     CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")

     ALLOCATE(sh_env%elements(c))
     ALLOCATE(sh_env%connlist1(4,aom_atoms))
     allocate(sh_env%aneighbour(0:atoms,0:atoms))
     allocate(sh_env%bneighbour(0:atoms,0:atoms))
     ALLOCATE(sh_env%pcoeffs(aom_atoms,norbitals))
     ALLOCATE(sh_env%scoeffs(aom_atoms))
     allocate(sh_env%coms(3,mol))
     !> crds is an AOM arrays (xyz, atoms, tm-t-tp) the third dimension
     !indicates what time this coordinates refer to
     allocate(sh_env%crds(3,c,3))
     allocate(sh_env%dist(mol,mol))
     ALLOCATE(sh_env%NN(mat, mat))
     allocate(sh_env%pcoeffsrun(atoms,mat,3)) !pcoeffsrun exists because the p-coefficients
                                               !need to be renormalised for each structure
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(sh_env%At(4,atoms))
     allocate(sh_env%Atp(4,atoms))
     allocate(sh_env%Atm(4,atoms))
     allocate(sh_env%Bt(4,atoms))
     allocate(sh_env%Btp(4,atoms))
     allocate(sh_env%Btm(4,atoms))
     allocate(sh_env%siteenergies(3,mat))  !COMMENT_SG: First dimension is to be removed
     !> nacv are memory consuming to store however here we are storing them at time
     !t and tp for the multiple time step feature avoing their calculations at
     !each step
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
     ! OZ: add extra arguments...
     CALL build_s_tables(file_name_s_psigma, file_name_s_ppi,&
                         file_name_s_psigma_C_S, file_name_s_ppi_C_S,&
                         file_name_s_psigma_S_S, file_name_s_ppi_S_S)
     
     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=sh_env%scaling)
     !> cbar factor : Sab = cbar*Sab_bar
     CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=sh_env%cbar)

     
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
     allocate(sh_env%connlist2(4,size(sh_env%atomlist)))
     do i=1,aom_atoms
          do j=1,norbitals
                sh_env%pcoeffs(i,j) = pcoefftot( start -1 + i + (j-1)*norbitals)
          end do
          sh_env%scoeffs(i) = scoefftot(start - 1 + i)
     end do

     !> CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     allocate(crds(3, c))
     ALLOCATE(Atp(4,aom_atoms))

     do coord=1,3
        crds(coord,:) = PACK(positions_tp(coord,:), sh_env%ext_pres)
     enddo

     Atp(2:4,1:aom_atoms) = crds(1:3,1:aom_atoms)
     !output: sh_env%connlist1, sh_env%connlist2
     call connect_list2(Atp(2:4,:), sh_env%atomlist, sh_env%connlist1, sh_env%connlist2, cutoff_connect)
     
     DEALLOCATE(atomlist_mask) 
     DEALLOCATE(submask) 
     DEALLOCATE(Atp)
     DEALLOCATE(crds)
!==================================================================================

     DEALLOCATE(rng_stream)
  END SUBROUTINE create_sh_env

  SUBROUTINE create_adiab_nace_env(mol, sh_env)
   
     INTEGER                       :: mol
     TYPE(sh_env_type)            :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_adiab_nace_env', &
     routineP = moduleN//':'//routineN

      ALLOCATE(sh_env%adiab_nace_env%T1(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%T2(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%DT(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%diab_D(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%DS(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%map_new2old(mol))
      ALLOCATE(sh_env%adiab_nace_env%map_old2new(mol))
      ALLOCATE(sh_env%adiab_nace_env%basic(mol))
      ALLOCATE(sh_env%adiab_nace_env%zeros(mol))
      ALLOCATE(sh_env%adiab_nace_env%not_associated(mol))
      ALLOCATE(sh_env%adiab_nace_env%Ttt(mol,mol))
      ALLOCATE(sh_env%adiab_nace_env%Sttsqrt_var(mol,mol))

  END SUBROUTINE create_adiab_nace_env

  SUBROUTINE propagate_aom_namd_first_step(nadiab, aom_atoms, norbitals,  &
                                  positions_tp,  &
                                  velocities_tp,  &
                                  energies_tp,  &
                                  sh_env, itimes, &
                                  adiab_section, my_bo_keyword,  &
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword, method_keyword, diagonal_F)

     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     INTEGER, INTENT(IN)                            :: aom_atoms
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: velocities_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, POINTER                               :: itimes
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section,optimization_section
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     CHARACTER(LEN=default_string_length)           :: method_keyword

     CHARACTER(len=*), PARAMETER :: routineN ='propagate_aom_namd_first_step', &
     routineP = moduleN//':'//routineN

     integer                                        :: coord, mol, atoms, c, & 
                                                       first_diabat, active_state, mat, orb 
     REAL(KIND=dp)                                  :: cutoff_sites, cutoff_connect, & 
                                                       scaling, deltatMD, deltatE
     REAL(KIND=dp)                                  :: cbar, minao
     COMPLEX(KIND=dp), dimension(:), allocatable    :: ctin
     LOGICAL                                        :: decohere, analytics
     INTEGER                                        :: new_state , ELsteps   
     REAL(KIND=dp)                                  :: dR, times 
     INTEGER                                        :: imol,  natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE     :: exact_forces
     LOGICAL                                        :: keyword_translate
     TYPE(sh_print_env_type)                        :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)           :: keyword_init, keyword_basis_repr, &
                                                        prop_wf_keyword
     real(kind=dp), DIMENSION(:,:), allocatable     :: PC_allrvecs
     LOGICAL                                        :: do_random
     REAL(kind=dp), DIMENSION(:), ALLOCATABLE       :: site_energies, couplings
     REAL(KIND=dp), dimension(:), allocatable       :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable    :: adiab_c
     TYPE(cp_logger_type), POINTER                  :: logger
     INTEGER                                        :: unit_nr
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER       :: diagonal_F
     TYPE(section_vals_type), POINTER               :: coupling_section
     CHARACTER(LEN=default_string_length)           :: coupling_keyword

     NULLIFY(positions_com_tp, velocities_com_tp)
     NULLIFY(logger)
     NULLIFY(optimization_section)
     NULLIFY(coupling_section)

     logger => cp_get_default_logger()

     CALL section_vals_val_get(adiab_section,"METHOD_COUPLING",c_val=coupling_keyword)
!------------------  AOM section -------------------------------------
     !> cutoff to define the neighborough list
     coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")

     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=cutoff_sites)
     !> minao is an atomic overlap cutoff
     CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=minao)
     !> dR is the increment for the gradient finite difference calculation
     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
             r_val=dR)
     !> analytics would probably calculates analytic forces for a 2 state sytems
     !(Spencer 2016)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=analytics)
     !> This remove bias coming from the finite difference calculation of the
      !> gradients of the Hamiltonian from AOM
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)
     !> Multiple time step calculation of aom gradient
     CALL section_vals_val_get(coupling_section,"MULTIPLE_TIME_STEP",&
             i_val=sh_env%multiple_time_step)
!-------------------------------------------------------------------

     !> Calculate the orbital overlap in the COM reference sistem of the
     !molecules
     CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS", l_val=keyword_translate)

     CALL section_vals_val_get(adiab_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ELsteps)

     !> Scaling factor Hab = scaling*Sab_bar
     scaling = sh_env%scaling
     cbar = sh_env%cbar
     atoms = sh_env%atoms !Number of atoms per site
     natom = sh_env%natom !Total number of atoms
     deltatMD = sh_env%deltat !Nuclear timestep

     !> Number of sites
     mol = nadiab ! CHANGE_AC
     !> Number of orbitals per atom
     orb = norbitals ! CHANGE_AC
     !> Size of the vectorial space where the wavefunctions are projected
     mat = orb*mol
     !> Numbers of active atoms 
     c = mol*atoms
     !> electronic timestep is a fraction of the molecular timestep
     deltatE = deltatMD/ELsteps !CHANGE_AC

     sh_print_env%adiab_section => adiab_section
     sh_print_env%nadiab = mat
     sh_print_env%norbitals = orb
     sh_print_env%deltatMD = deltatMD

     ALLOCATE(sh_print_env%adiab_energies_t(mat))
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))
     allocate(ctin(mat))      !ctin goes into the subroutine ZGEMV
                              !these are carrier wf coefficient in the diab
                              !basis

     !READ THE COORDINATES, ENERGIES AND STATE
     IF (method_keyword.eq."FROZEN_HAMILTONIAN") THEN
        velocities_com_tp = 0.0D0
     ELSE
        CALL translate_com(velocities_tp, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     ENDIF
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  

!------------------  AOM section -------------------------------------
     do coord=1,3
        sh_env%crds(coord,:,1)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
        !At time 0, we must guess the positions at time -1 with the inital velocities
        sh_env%crds(coord,:,2)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres) -  &
                          PACK(velocities_com_tp(coord,:), sh_env%ext_pres)*deltatMD
        sh_env%crds(coord,:,3)   = PACK(positions_com_tp(coord,:), sh_env%ext_pres)
        sh_env%siteenergies(2,:) = 0.0D0 ! CHANGE_AC
        sh_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo
    
     !> CALCULATE THE CENTER OF MAS AND CONSTRUCT NEIGHBOURS LIST
     CALL calc_coms_and_nlist_aom(sh_env, atoms, mol, cutoff_sites, orb)

     ! Allocate neighbouring list just at the first step as Jacob suggested 
     !only if sh_env%multiple_time_step .eq. 1 otherwise we use Multiple time
     !step feature and to save memory we recalculate bneighbour and aneighbour
     !every N steps. This is done also for the NACV giving a huge speed-up. NOTE
     !that in general is not convenient to recalculate neighboring list at each
     !step because is to expensive to do and not useful.
     if (sh_env%multiple_time_step .eq. 1) then 
       allocate(sh_env%neighbourlist(mol,mol,0:atoms,0:atoms))
       sh_env%neighbourlist = 0
     end if 

     !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
     CALL calc_aom_Xmatrix_elements_first_step(sh_env, atoms, mol, dR, do_random,cbar, &
                      scaling, PC_allrvecs, orb, minao)

     !> CONVERT NACV AND CALCULATE EIGENPROPERTIES 
     !> THESE ASSOCIATIONS ARE NEEDED FOR THE NACV CALC
     sh_env%nacv = cbar*sh_env%nacv
!-----------------------------------------------------------------------------

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian H_tp
     !It gives us back the eigenvalues E_tp and eigenvectors C_tp
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(sh_env%H_tp), sh_env%C_tp, sh_env%E_tp)
     
     !> PRINT
     !> print site energies 
     allocate(site_energies(mat))
     allocate(couplings(mat-1))
     do imol=1,mat
         site_energies(imol) = -AIMAG(sh_env%H_tp(imol,imol))
         if (imol .ne. mat) then
           couplings(imol) = -AIMAG(sh_env%H_tp(imol,imol+1))
         end if
     end do
   
     times = DBLE(itimes)

     !> print info in the output only if we do not restart the file
     sh_print_env%times = times
     sh_print_env%times = DBLE(itimes)
     if (.NOT. sh_env%restart_key) then
        CALL adiab_print_site_energies(adiab_section, site_energies, nadiab,deltatMD, times) 
        CALL adiab_print_off_diagonals(adiab_section, couplings, nadiab-1, deltatMD, times)
        CALL adiab_print_adiabat_energies(adiab_section, sh_env%E_tp, nadiab, deltatMD, times)
        call adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, atoms, nadiab)
        CALL adiab_print_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        CALL adiab_print_pseudo_hamiltonian(adiab_section, -AIMAG(sh_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        CALL adiab_print_nacv(adiab_section, sh_env%nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
     end if

     !> INITIALIZE THE WAVEFUNCTION (WE DO NOT PROPAGATE IT AT THE FIRST STEP) 
     IF  (.NOT. my_test_keyword) THEN
         decohere = .FALSE.
         sh_env%decohere = .FALSE.
         sh_env%attempt_hop = .FALSE.
         !if we want to restart wavefunction then
         if (sh_env%restart_key) then
             keyword_init = "DIABATIC_RESTART"
         else
             CALL section_vals_val_get(adiab_section, "INITIALIZATION", c_val=keyword_init)
         end if

         CALL section_vals_val_get(adiab_section, "REPRESENTATION", c_val=keyword_basis_repr)
         SELECT CASE(keyword_init)
         CASE("DIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_DIABAT", &
                                       i_val=first_diabat)
             CALL choose_firstadiabat(sh_env%C_tp, first_diabat, sh_print_env, active_state, sh_env) 
             !an active state in needed for the nuclear propagation but the wafunction has to
             !start form an eigenvector of the hamiltonian in the diab_basis
             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(first_diabat) = CMPLX(1.0D0, 0.0D0,dp)  
             CASE("ADIABATIC_BASIS")
                 !!!!! THIS SECTION IS DEPRECATED !!!!!!
                 !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
                 CALL nacv_first_ts(sh_env,diagonal_F, natom, mat, velocities_com_tp) 
                 ctin(:) = sh_env%C_tp(first_diabat,:)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("ADIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_ADIABAT",&
                                      i_val=active_state)
             !> print init
             sh_print_env%first_adiabat    = active_state
             ALLOCATE(sh_print_env%pop_init(SIZE(sh_env%C_tp(:,active_state))))
             sh_print_env%pop_init(:)      = 0.0D0
             sh_print_env%pop_init(active_state)  = 1.0D0
             CALL adiab_print_sh_info(sh_print_env, keyword = 'START')
             DEALLOCATE(sh_print_env%pop_init)

             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = sh_env%C_tp(:, active_state)
             CASE("ADIABATIC_BASIS")
                 !!!!! THIS SECTION IS DEPRECATED !!!!!!
                 !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
                 CALL nacv_first_ts(sh_env, diagonal_F, natom, mat, velocities_com_tp) 
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(active_state) = CMPLX(1.0D0, 0.0D0,dp)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("DIABATIC_RESTART")
                 ctin(:) = sh_env%readcoeff(:)
                 active_state = sh_env%active_state_res
                 ! do not print this duplicated step in run-sh.log
                 !sh_print_env%first_adiabat = active_state 
                 !CALL adiab_print_sh_info(sh_print_env, keyword = 'START') 
         CASE DEFAULT
             CPABORT("INITIALIZATION implemented: DIABATIC, ADIABATIC")
         END SELECT

         allocate(adiab_c(mat)) 
         allocate(adiab_pop(mat))

         !> print info in the output only if we do not restart the file
         if (.NOT. sh_env%restart_key) then
             !NB: this subroutine is wasted computational work if adiab_pop is not
             !printed out, same probably for AOM.main

             unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIAB_POPULATIONS",&
               extension=".xyz", middle_name="adiab_pop",log_filename=.FALSE.)
             if (unit_nr > 0) then 
                CALL diab_to_adiab( sh_env%C_tp, ctin, adiab_c, adiab_pop)     
                CALL adiab_print_adiab_pop(adiab_section, adiab_pop, & 
                                  nadiab,norbitals,itimes,deltatMD, times)
             end if 
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, times)
             CALL adiab_print_pop(adiab_section, abs(ctin),nadiab,norbitals,itimes,deltatMD, times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes, &
                                        deltatMD, times)
         end if
     ELSE 
         sh_env%attempt_hop = .TRUE.
         active_state = 1
         new_state = mod(active_state + 1, mat) 
         IF ((new_state).EQ.0) new_state = mat
     ENDIF 

!------------------  AOM section --------------------------------------------------------
     !> second part CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE
     !RESULT OF THE SURFACE HOP ROUTINE IS.
     print*, "FOBSH| ACTIVE STATE = ", active_state
     call calcforces_aom(state=active_state,Hin=-AIMAG(sh_env%H_tp), C=sh_env%C_tp, nacv=sh_env%nacv, atoms=atoms, & 
                     scaling=scaling, cbar=cbar, Fin=diagonal_F, forces_adiab=sh_env%forces_adiab, &
                     ext_pres=sh_env%ext_pres)

     IF (analytics) THEN
        times = DBLE(itimes)
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF
!-----------------------------------------------------------------------------------------------

     !> This selects the most general DIABATIC_RK_FAST that does not include
     !diabatic nace calculation
     optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
     CALL section_vals_val_get(optimization_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
     if (prop_wf_keyword .eq. "DIABATIC_RK_FAST") then 
         ! in this case diabatic nace are never calculated, neither in the
         ! propagate_wf nor in the adiabatic nace calculation for the hopping
         ! probability
         ! so that all the Sttp and co matrices are also not used and so the memory
         ! allocate for thme can be deallocated 
         sh_env%no_diab_nace = .True.
     end if

     
     !> STORE THE COORDINATES IF NOT SUCCESSFUL RUN
!------------------  AOM section --------------------------------------------------------
     sh_env%rm_tm      = positions_com_tp - velocities_com_tp*deltatMD
     sh_env%rm_t       = positions_com_tp
     sh_env%Stt = sh_env%Stptp
!------------------------------------------------------------------------------------- 
     sh_env%sh_print_env = sh_print_env
     sh_env%ctin        = ctin
     sh_env%ctin_t        = ctin
     sh_env%H_t = sh_env%H_tp 
     sh_env%diabatic_nace_t = sh_env%diabatic_nace_tp
     sh_env%E_t = sh_env%E_tp
     sh_env%C_t = sh_env%C_tp
     sh_env%new_state = new_state
     sh_env%active_state = active_state
     !> IMPORTANT FOR THE FIRST TIME STEP
     sh_env%old_state = sh_env%active_state
     sh_env%new_state = -1
     sh_env%successful_hop = .FALSE.


    !> DEALLOCATE
    deallocate(ctin)
    deallocate(exact_forces)
    deallocate(site_energies, couplings, adiab_c, adiab_pop)
    DEALLOCATE(positions_com_tp)
    DEALLOCATE(velocities_com_tp)
    DEALLOCATE(PC_allrvecs)    

  END SUBROUTINE propagate_aom_namd_first_step

  !Seems to change the nace not the nacv 
  SUBROUTINE nacv_first_ts(sh_env, diagonal_F, natom, mat, velocities_com_tp)

     TYPE(sh_env_type), INTENT(INOUT)              :: sh_env
     INTEGER, intent(in)                            :: natom
     integer, intent(in)                            :: mat
     REAL(kind=dp), DIMENSION(:,:), INTENT(IN)      :: velocities_com_tp
     REAL(kind=dp), DIMENSION(:,:,:), POINTER, & 
                                    INTENT(IN)      :: diagonal_F

     CHARACTER(len=*), PARAMETER                    :: routineN ='nacv_first_ts', &
     routineP = moduleN//':'//routineN
     
     integer                                        :: i,j, iparticle, coord
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv

     sh_env%adiab_nace_tp = 0.0D0
     do i =1, mat
         do j = 1, mat
              if (j .ge. i) then
                  CALL calc_hop_nacv_aom(sh_env, diagonal_F, hop_nacv, i, j, fast=.TRUE.)
                  do iparticle=1,natom
                       do coord = 1,3
                            sh_env%adiab_nace_tp(i,j) = sh_env%adiab_nace_tp(i,j) + &
                                                    velocities_com_tp(coord, iparticle)* & 
                                                            hop_nacv(coord, iparticle)
                       end do
                  end do
              else
                  sh_env%adiab_nace_tp(i,j)= -sh_env%adiab_nace_tp(j,i)
              end if 
         end do
     end do
     sh_env%adiab_nace_tp = -sh_env%adiab_nace_tp 
  END SUBROUTINE nacv_first_ts


END MODULE aom_initialization
