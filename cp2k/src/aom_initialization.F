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
                                              calc_sab, &
                                              sab_screen_new, &
                                              calc_dRSab, &
                                              calc_rvecs, &
                                              calc_sab_pete, &
                                              build_s_tables

     USE aom_utils,                    ONLY: makeX,    &
                                             make_adiabX,      &
                                             adiabaticity,   &
                                             !surfacehop,     &
                                             calcforces,     &
                                             read_listcoeff, &
                                             read_listdiabcoeff, &
                                             choose_firstadiabat, &
                                             simple_vel_rescaling, &
                                             translate_com, &
                                              calc_hop_nacv, &
                                              cal_delta_energy_nacv, &
                                              diab_to_adiab, &
                                              print_matrix_standard_out

     USE aom_types,                    ONLY: aom_env_type, &
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

     USE non_adiabatic_utils,            ONLY: decohere_wf,  &
                                               propagate_wf, &
                                               propagate_adiab_wf, &
                                                  prop_wf_HS_corr, & 
                                               classically_forbidden_hops


     USE parallel_rng_types,           ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             rng_stream_type
 
     USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                                cp_logger_type

     USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                                cp_print_key_should_output


#include "./base/base_uses.f90"

     IMPLICIT NONE

     PRIVATE
 
     CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_initialization'   

     PUBLIC :: create_aom_env, &
               propagate_aom_namd_first_step

CONTAINS

  SUBROUTINE create_aom_env(nadiab, aom_atoms, norbitals, natom, &
                                  positions_tp,     &
                                  coupling_section, &
                                  adiab_section,    &
                                  aom_env)


     INTEGER, INTENT(IN)                            :: nadiab,  aom_atoms, norbitals, natom
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     TYPE(section_vals_type), POINTER               :: coupling_section
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(section_vals_type), POINTER               :: wf_restart_section
     TYPE(aom_env_type), INTENT(INOUT)              :: aom_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_aom_env', &
     routineP = moduleN//':'//routineN

     integer :: i, j, mat, c, coord, mol, atoms, start
     REAL(KIND=dp) :: cutoff_connect
     REAL(KIND=dp), dimension(:,:), allocatable :: crds
     REAL(KIND=dp), dimension(:,:), allocatable :: Atp
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pcoefftot, scoefftot
     LOGICAL, DIMENSION(:), ALLOCATABLE :: atomlist_mask, submask
     integer, dimension(:), ALLOCATABLE           :: ext_elements
     TYPE(cp_sll_val_type), POINTER               :: listcoeff
     TYPE(cp_sll_val_type), POINTER               :: listdiabcoeff
!     TYPE(section_vals_type), POINTER             :: input_key
!     TYPE(cp_logger_type), POINTER                :: logger

!     INTEGER                              :: should_output
 !    LOGICAL                              :: my_do_restart = .FALSE.

     LOGICAL    :: do_random
!  test_key
     REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
     INTEGER    :: input_seed
     TYPE(rng_stream_type), POINTER           :: rng_stream
     LOGICAL                                  :: continue_while
     CHARACTER(LEN=default_string_length)     :: file_name_s_psigma, file_name_s_ppi, deco_keyword
!     COMPLEX(KIND=dp), dimension(:), allocatable  :: readcoeff
!     INTEGER :: active_state_res

     NULLIFY(rng_stream)
     NULLIFY(listcoeff)
!     NULLIFY(logger)
!     logger => cp_get_default_logger()


     mol = nadiab ! CHANGE_AC
     atoms = aom_atoms ! CHANGE_AC

     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)
   

     mat = norbitals*nadiab !> Size of the vectorial space where the wavefunctions are projected
     c = nadiab*aom_atoms !> Numbers of active atoms

     ALLOCATE(aom_env%elements(c))
     ALLOCATE(aom_env%neighbourlist(nadiab, nadiab,0:aom_atoms,0:aom_atoms))
     ALLOCATE(aom_env%connlist1(4,aom_atoms))
     ALLOCATE(aom_env%bneighbour(3,0:aom_atoms))
     ALLOCATE(aom_env%pcoeffs(aom_atoms,norbitals))
     ALLOCATE(aom_env%scoeffs(aom_atoms))

     allocate(aom_env%coms(3,mol))
     allocate(aom_env%crds(3,c,3))
     allocate(aom_env%dist(mol,mol))
     ALLOCATE(aom_env%NN(mat, mat))
     allocate(aom_env%pcoeffsrun(atoms,mat,3)) !pcoeffsrun exists because the p-coefficients
                                               !need to be renormalised for each structure
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(aom_env%At(4,atoms))
     allocate(aom_env%Atp(4,atoms))
     allocate(aom_env%Atm(4,atoms))
     allocate(aom_env%Bt(4,atoms))
     allocate(aom_env%Btp(4,atoms))
     allocate(aom_env%Btm(4,atoms))
     allocate(aom_env%siteenergies(3,mat))
     ALLOCATE(aom_env%nacv(3,mat, mat*aom_atoms))
     allocate(aom_env%F(3,natom,mol)) !array of diagonal atomic forces
     ALLOCATE(aom_env%forces_adiab(3,natom))
     ALLOCATE(aom_env%forces_adiab_new(3,natom))
     ALLOCATE(aom_env%diabatic_nace_tp(mat,mat))
     ALLOCATE(aom_env%diabatic_nace_t(mat,mat))

     aom_env%At = 0.0_dp
     aom_env%Atp = 0.0_dp
     aom_env%Atm = 0.0_dp
     aom_env%Bt  = 0.0_dp
     aom_env%Btp = 0.0_dp
     aom_env%Btm = 0.0_dp
     aom_env%siteenergies = 0.0_dp
     aom_env%nacv = 0.0_dp
     aom_env%F = 0.0_dp
     aom_env%forces_adiab = 0.0_dp
     aom_env%forces_adiab_new = 0.0_dp
     aom_env%diabatic_nace_tp = 0.0_dp
 
     ALLOCATE(aom_env%rm_t(3, natom))
     ALLOCATE(aom_env%rm_tm(3, natom))
     ALLOCATE(aom_env%velocities_tp(3, natom))
     ALLOCATE(aom_env%C_t(mat,mat))
     ALLOCATE(aom_env%C_tp(mat,mat))
     ALLOCATE(aom_env%E_t(mat))
     ALLOCATE(aom_env%E_tp(mat))
     ALLOCATE(aom_env%H_t(mat,mat))
     ALLOCATE(aom_env%H_tp(mat,mat))
     ALLOCATE(aom_env%Stt(mat,mat))
     ALLOCATE(aom_env%Stptp(mat,mat))
     ALLOCATE(aom_env%Sttp(mat,mat))
     ALLOCATE(aom_env%Stpt(mat,mat))
     ALLOCATE(aom_env%Stmt(mat,mat))
     ALLOCATE(aom_env%Sttm(mat,mat))
     ALLOCATE(aom_env%Stmtm(mat,mat))
     ALLOCATE(aom_env%ctin(mat))
     ALLOCATE(aom_env%ctin_t(mat))

     ALLOCATE(aom_env%ctin_before_deco(mat))
     ALLOCATE(aom_env%ctin_after_deco(mat))
     ALLOCATE(aom_env%ctin_after_tranf_correction(mat))
     ALLOCATE(aom_env%active_region(mat))
   
     CALL create_adiab_nace_env(mol, aom_env)

     aom_env%futur_decohere = .FALSE.
     aom_env%decohere       = .FALSE.
     aom_env%no_diab_nace       = .FALSE.

     aom_env%new_state   = -1
     aom_env%old_state   = -1
     aom_env%active_state   = -1
     ALLOCATE(aom_env%adiab_overlap(mat, mat))
     aom_env%attempt_hop = .FALSE.
     aom_env%successful_hop = .FALSE.
     
     aom_env%HS_propagation = .FALSE.
     ALLOCATE(aom_env%class_forbidden_hop_t(mat))
     allocate(aom_env%adiab_nace_t(mat,mat))
     allocate(aom_env%adiab_nace_tp(mat,mat))

     !> Create the overlap look-up table for the AOM coupling calculation
     CALL section_vals_val_get(coupling_section,"PSIGMA_FILE_NAME",&
             c_val=file_name_s_psigma)
     CALL section_vals_val_get(coupling_section,"PPI_FILE_NAME",&
             c_val=file_name_s_ppi)
     CALL build_s_tables(file_name_s_psigma, file_name_s_ppi)

     !> Scaling factor Hab = scaling*Sab_bar
     CALL section_vals_val_get(coupling_section,"SCALING_FACTOR",&
             r_val=aom_env%scaling)
     !> cbar factor : Sab = cbar*Sab_bar
     CALL section_vals_val_get(coupling_section,"CBAR",&
             r_val=aom_env%cbar)


     !> INITIALIZE THE RANDOM GENERATOR
     CALL section_vals_val_get(adiab_section,"SEED",i_val=input_seed)
     initial_seed(:,:) = REAL(input_seed,KIND=dp)
     !CALL create_rng_stream(rng_stream=aom_env%rng_stream,&
     CALL create_rng_stream(rng_stream=rng_stream,&
                             name="aom fssh",&
                             distribution_type=UNIFORM, &
                             seed=initial_seed)
     !stop
     ALLOCATE(aom_env%rng_stream)
     aom_env%rng_stream = rng_stream

  
!     input_key => section_vals_get_subs_vals(adiab_section,&
!                                       "WAVEFUNCTION_RESTART")

!     should_output = cp_print_key_should_output(&
!                     iteration_info=logger%iter_info,&
!                     basis_section=input_key)

     
!     IF (should_output .ne. 0) THEN
!        my_do_restart = .TRUE.
!     END IF
!print *, "my_do_restart", my_do_restart

     !> Get wavefunction diabatic coefficients from the input file
     !this section should be optional somehow---STILL TO DO USING e.g.
     !STEP_START_EVAL keyword in the input 
     wf_restart_section => section_vals_get_subs_vals(adiab_section,"WAVEFUNCTION_RESTART")
     CALL section_vals_val_get(wf_restart_section,"RESTART_KEY", l_val=aom_env%restart_key)
     print *, "RESTART_KEY", aom_env%restart_key
     if (aom_env%restart_key) then
         CALL section_vals_val_get(wf_restart_section,"ACTIVE_STATE_RESTART", i_val=aom_env%active_state_res)
         CALL section_vals_list_get(wf_restart_section,"_DEFAULT_KEYWORD_",list=listdiabcoeff) !Read them in the input files
         CALL read_listdiabcoeff(listdiabcoeff, norbitals, nadiab, ctin=aom_env%readcoeff) !Transform them in the complex ctin
     end if
     
     !> READ THE AOM COEFFICIENT
     CALL section_vals_list_get(coupling_section,"_DEFAULT_KEYWORD_",list=listcoeff) !Read them in the input files
     CALL read_listcoeff(listcoeff, ext_elements, aom_env%ext_pres, &
                         scoefftot, pcoefftot, natom) !Transform them into four arrays
     aom_env%elements = PACK(ext_elements, aom_env%ext_pres) !Only select active atoms
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
     allocate(aom_env%atomlist(COUNT(submask(1:aom_atoms))))
     aom_env%atomlist = PACK( (/ (i, i=1,aom_atoms) /), submask(1:aom_atoms) ) !List atoms with orbitals
     allocate(aom_env%connlist2(4,size(aom_env%atomlist)))
     do i=1,aom_atoms
          do j=1,norbitals
                aom_env%pcoeffs(i,j) = pcoefftot( start -1 + i + (j-1)*norbitals)
          end do
          aom_env%scoeffs(i) = scoefftot(start - 1 + i)
     end do


     !> CONSTRUCT THE CONNLIST
     !we fill out the connlist arrays before calculating the nacv
     allocate(crds(3, c))
     do coord=1,3
        crds(coord,:) = PACK(positions_tp(coord,:), aom_env%ext_pres)
     enddo

     ALLOCATE(Atp(4,aom_atoms))
     Atp(2:4,1:aom_atoms) = crds(1:3,1:aom_atoms)
     call connect_list2(Atp(2:4,:), aom_env%atomlist, aom_env%connlist1, aom_env%connlist2, cutoff_connect)

     
     aom_env%do_deco = .FALSE.
     CALL section_vals_val_get(adiab_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
     IF (.NOT.(deco_keyword.EQ."NO_DECO_CORR"))  aom_env%do_deco = .TRUE.

     aom_env%do_reorder = .FALSE.
     CALL section_vals_val_get(adiab_section,"REORDERING_STATES_USING_OVERLAP",l_val=aom_env%do_reorder)

     DEALLOCATE(rng_stream)
  END SUBROUTINE create_aom_env

  SUBROUTINE create_adiab_nace_env(mol, aom_env)
   
     INTEGER                       :: mol
     TYPE(aom_env_type)            :: aom_env

     CHARACTER(len=*), PARAMETER :: routineN ='create_adiab_nace_env', &
     routineP = moduleN//':'//routineN

      ALLOCATE(aom_env%adiab_nace_env%T1(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%T2(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%DT(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%diab_D(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%DS(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%map_new2old(mol))
      ALLOCATE(aom_env%adiab_nace_env%map_old2new(mol))
      ALLOCATE(aom_env%adiab_nace_env%basic(mol))
      ALLOCATE(aom_env%adiab_nace_env%zeros(mol))
      ALLOCATE(aom_env%adiab_nace_env%not_associated(mol))
      ALLOCATE(aom_env%adiab_nace_env%Ttt(mol,mol))
      ALLOCATE(aom_env%adiab_nace_env%Sttsqrt_var(mol,mol))

  END SUBROUTINE create_adiab_nace_env

  SUBROUTINE propagate_aom_namd_first_step(nadiab, aom_atoms, norbitals,  &
                                  positions_tp,  &
                                  velocities_tp,  &
                                  energies_tp,  &
                                  aom_env, itimes, &
                                  global_forces, coupling_section, &
                                  adiab_section, my_bo_keyword,  &
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword, method_keyword)

     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     INTEGER, INTENT(IN)                            :: aom_atoms
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: velocities_tp
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     TYPE(aom_env_type), INTENT(INOUT)              :: aom_env
     INTEGER, POINTER                               :: itimes
     TYPE(mixed_force_type), DIMENSION(:), POINTER  :: global_forces
     TYPE(section_vals_type), POINTER               :: coupling_section
     LOGICAL, INTENT(IN)                            :: my_bo_keyword
     TYPE(section_vals_type), POINTER               :: adiab_section
     TYPE(atomic_kind_list_type), POINTER           :: atomic_kinds
     TYPE(distribution_1d_type), POINTER            :: local_particles
     LOGICAL, INTENT(IN)                            :: my_test_keyword
     CHARACTER(LEN=default_string_length)           :: method_keyword
     !TYPE(atomic_kind_type), DIMENSION(:), POINTER  :: atomic_kind_set

     CHARACTER(len=*), PARAMETER :: routineN ='propagate_aom_namd_first_step', &
     routineP = moduleN//':'//routineN

     integer :: coord 
     integer :: mol, i, j, k, l, m, d, atoms, c,  first_diabat, active_state, mat, orb
     integer, dimension(:,:), allocatable :: NN,  aneighbour, bneighbour  
     REAL(KIND=dp) :: cutoff_sites, cutoff_connect, scaling, deltatMD, deltatE
     REAL(KIND=dp) :: S
     REAL(KIND=dp)    :: cbar, minao
     REAL(KIND=dp), dimension(:,:,:), allocatable :: nacv
     REAL(KIND=dp), dimension(:,:,:), allocatable :: crds, pcoeffsrun
     REAL(KIND=dp), dimension(:,:), allocatable :: Stptp
     REAL(KIND=dp), dimension(:,:), allocatable :: dist
     COMPLEX(KIND=dp), dimension(:), allocatable :: ctin
     REAL(KIND=dp), dimension(:,:), allocatable ::  Atp, Btp, coms, siteenergies
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE ::  scoefftot
     REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: F
     LOGICAL                                      :: decohere, analytics
     INTEGER                                      :: new_state , ELsteps   
     REAL(KIND=dp)                                :: dR, times 
     REAL(KIND=dp)                                :: decohere_criterion
     INTEGER                                      :: imol, iatom, natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: exact_forces
     LOGICAL                                      :: keyword_translate
     TYPE(sh_print_env_type)                      :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER        :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)         :: keyword_init, keyword_basis_repr, &
                                                     prop_wf_keyword
     real(kind=dp), DIMENSION(:,:), allocatable :: PC_allrvecs
     real(kind=dp), allocatable :: PC_rvecs(:,:)
     LOGICAL    :: do_random
     REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: site_energies, couplings
     REAL(KIND=dp), dimension(:), allocatable     :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable               :: adiab_c
     TYPE(cp_logger_type), POINTER                :: logger
     INTEGER                                    :: unit_nr
     !INTEGER                                    :: nparticle_kind
     !INTEGER, dimension(:), allocatable          :: class_forbidden_hop_tp

     NULLIFY(positions_com_tp, velocities_com_tp)
     NULLIFY(logger)

     logger => cp_get_default_logger()


     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_BETWEEN_SITES",&
             r_val=cutoff_sites)
     !> cutoff to define the neighborough list
     CALL section_vals_val_get(coupling_section,"CUT_OFF_CONNECTIVITY",&
             r_val=cutoff_connect)
     !minao is an atomic overlap cutoff
     CALL section_vals_val_get(coupling_section,"ATOMIC_OVERLAP_CUTOFF",&
             r_val=minao)
     CALL section_vals_val_get(coupling_section,"ELECTRONIC_PARTIAL_STEP",&
             i_val=ELsteps)
     CALL section_vals_val_get(coupling_section,"DECOHERENCE_CRITERION",&
             r_val=decohere_criterion)
     CALL section_vals_val_get(coupling_section,"NACV_INCREMENT",&
             r_val=dR)
     CALL section_vals_val_get(coupling_section,"ANALYTICS", l_val=analytics)
     CALL section_vals_val_get(coupling_section,"RANDOM_DERIVATIVE", l_val=do_random)
     CALL section_vals_val_get(coupling_section,"CENTER_OF_MASS", l_val=keyword_translate)

     scaling = aom_env%scaling
     cbar = aom_env%cbar
     atoms = aom_env%atoms !Number of atoms per site
     natom = aom_env%natom !Total number of atoms
     deltatMD = aom_env%deltat !Nuclear timestep

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
     ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))
     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stptp(mat,mat))
     allocate(crds(3,c,3))
     !allocate(NACE_ad_tp(mat,mat))
     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     !allocate(nacv(3,mat,mat,mat*atoms))      !non-adiabatic coupling vectors
     allocate(nacv(3,mat,mat*atoms))      !non-adiabatic coupling vectors
     nacv = 0.0D0
     allocate(siteenergies(3,mat))      !diagonals of the Hamiltonian
     allocate(pcoeffsrun(atoms,mat,3)) !pcoeffsrun exists because the p-coefficients need to be renormalised for each structure
     allocate(scoefftot(c))
     allocate(coms(3,mol))
     allocate(NN(mat,mat))
     allocate(dist(mol,mol))
     !allocate the arrays that will contain atomic coordinates for the overlap routine
     !first column is now atomic number of atom (e.g. 6 for carbon)
     !Naming convention: At = A(t), Atp = A(t+1), Atm = A(t-1), etc
     allocate(Atp(4,atoms))
     allocate(Btp(4,atoms))
     allocate(F(3,natom,mol)) !array of diagonal atomic forces
     !these arrays are data for the non-adiabatic coupling calculation
     !they help efficiency by telling the program which atoms are connected
     !and thus which pairs to carry out the calculation for
     allocate(aneighbour(0:atoms,0:atoms))
     allocate(bneighbour(0:atoms,0:atoms))
     aneighbour = 0
     bneighbour = 0
     !allocate(neighbourlist(mol,mol,0:atoms,0:atoms))
     aom_env%neighbourlist = 0


     !READ THE COORDINATES, ENERGIES AND STATE
     IF (method_keyword.eq."FROZEN_HAMILTONIAN") THEN
        velocities_com_tp = 0.0D0
     ELSE
        CALL translate_com(velocities_tp, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     ENDIF
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     do coord=1,3
        crds(coord,:,1) = PACK(positions_com_tp(coord,:), aom_env%ext_pres)
        !At time 0, we must guess the positions at time -1 with the inital velocities
        crds(coord,:,2) = PACK(positions_com_tp(coord,:), aom_env%ext_pres) -  &
                          PACK(velocities_com_tp(coord,:), aom_env%ext_pres)*deltatMD
        crds(coord,:,3) = PACK(positions_com_tp(coord,:), aom_env%ext_pres)
        siteenergies(2,:) = 0.0D0 ! CHANGE_AC
        siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo


     !CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     coms = 0
     do i=1,mol
          do m=1,atoms
               do d=1,3
                    coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do


     !CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               dist(i,j) = ((coms(1,i) - coms(1,j))*(coms(1,i) - coms(1,j)) + (coms(2,i) - coms(2,j))*&
               &(coms(2,i) - coms(2,j)) + (coms(3,i) - coms(3,j))*(coms(3,i) - coms(3,j)))
               if (dist(i,j) < cutoff_sites*cutoff_sites) then
                    do d=1,orb
                         do m=1,orb
                              NN(orb*(i-1)+d,orb*(j-1)+m) = 1
                         end do
                    end do
               else
                    do d=1,orb
                         do m=1,orb
                              NN(orb*(i-1)+d,orb*(j-1)+m) = 0
                         end do
                    end do
               end if                    
          end do
     end do


     !CALCULATE THE HAMILTONIAN AND THE NACV
     do i=1,mol
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               Atp(1,1:atoms) = dble(aom_env%elements((i-1)*atoms+1:i*atoms))
               do j=1,atoms
                    Atp(2:4,j) = crds(1:3,(i-1)*atoms+j,3)
               end do

               !then for each other molecule
               do k=1,mol
                    !create the array of its atomic coordinates at both t and t+dt
                    Btp(1,1:atoms) = dble(aom_env%elements((k-1)*atoms+1:k*atoms))
                    do l=1,atoms
                        Btp(2:4,l) = crds(1:3,(k-1)*atoms+l,3)
                    end do
                    if (i == 1) then
                         do d=1,orb
                              pcoeffsrun(1:atoms,orb*(k-1)+d,3) = aom_env%pcoeffs(1:atoms,d)
                              call calc_sab(Btp, aom_env%connlist2, pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if  ! (i == 1)
                    do d=1,orb
                         do m=1,orb
                              call sab_screen_new(Atp, aom_env%atomlist, aom_env%connlist1, pcoeffsrun(:,orb*(i-1)+d,3),&
                              & aom_env%atomlist, Btp, aom_env%atomlist, aom_env%connlist1, pcoeffsrun(:,orb*(k-1)+m,3),&
                              & aneighbour(:,1:atoms), bneighbour(:,1:atoms), minao)
                         end do
                    end do
                    aom_env%neighbourlist(i,k,:,1:atoms) = bneighbour(:,1:atoms)
                    aom_env%neighbourlist(k,i,:,1:atoms) = aneighbour(:,1:atoms)
                    if (NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              do d=1,orb
                                   do m=1,orb
                                        !These diagonals are 1 by definition
                                        Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                                        aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, (-siteenergies(3,k)),dp)
                                   end do
                              end do

                         !now the off-diagonals                         
                         else ! ( i == k)
                              !symmetric matrices
                              !no sense calling calc_sab twice as often as you need to
                              if (k < i) then
                                   do m=1,orb
                                        do d=1,orb
                                             aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = aom_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = Stptp(orb*(k-1)+m,orb*(i-1)+d)
                                        end do
                                   end do
                              else
                                   do d=1,orb
                                        do m=1,orb
                                             call calc_sab(Atp, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                                              S, Btp, aom_env%connlist2, pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                                             Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                                             aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, -scaling*S,dp)

                                             call calc_dRSab(Atp, aom_env%atomlist, aom_env%connlist1, &
                                             pcoeffsrun(1:atoms,orb*(i-1)+d,3), aom_env%scoeffs, &
                                             !nacv(1:3,orb*(k-1)+d,orb*(i-1)+m,(i-1)*atoms+1:i*atoms), &
                                             !nacv(1:3,orb*(i-1)+m,orb*(k-1)+d,(k-1)*atoms+1:k*atoms), &
                                             nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                                             nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                                             Btp, aom_env%atomlist, aom_env%connlist1, pcoeffsrun(1:atoms,orb*(k-1)+m,3),&
                                             aom_env%scoeffs, aom_env%neighbourlist(k,i,:,:), aom_env%neighbourlist(i,k,:,:), dR, &
                                             do_random, aom_env)
                                        end do
                                   end do
                              endif 
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                              end do
                         end do
                    end if !(NN(i,k) == 1
               end do !loop over k 
     enddo !loop-over i


     !> CONVERT NACV AND CALCULATE EIGENPROPERTIES 
     !> THESE ASSOCIATIONS ARE NEEDED FOR THE NACV CALC
     nacv = cbar*nacv
!TEST NO NACV 
!nacv = 0.0D0

     CALL adiabaticity(-AIMAG(aom_env%H_tp), aom_env%C_tp, aom_env%E_tp)
     
     
     !> first part: CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
     !> HOP ROUTINE IS. 
     DO imol=1,mol
        DO iatom=1, natom
           F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
        ENDDO
     ENDDO

     !> PRINT
     !> print site energies 
     allocate(site_energies(mat))
     allocate(couplings(mat-1))
     do imol=1,mat
         site_energies(imol) = -AIMAG(aom_env%H_tp(imol,imol))
         if (imol .ne. mat) then
           couplings(imol) = -AIMAG(aom_env%H_tp(imol,imol+1))
         end if
     end do
   
times = DBLE(itimes)


     !> print info in the output only if we do not restart the file
     sh_print_env%times = times
     sh_print_env%times = DBLE(itimes)
     if (.NOT. aom_env%restart_key) then
        CALL adiab_print_site_energies(adiab_section, site_energies, nadiab,deltatMD, times) 
        CALL adiab_print_off_diagonals(adiab_section, couplings, nadiab-1, deltatMD, times)
        CALL adiab_print_adiabat_energies(adiab_section, aom_env%E_tp, nadiab, deltatMD, times)
        call adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, atoms, nadiab)
        CALL adiab_print_hamiltonian(adiab_section, -AIMAG(aom_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        CALL adiab_print_pseudo_hamiltonian(adiab_section, -AIMAG(aom_env%H_tp), nadiab,norbitals,itimes, deltatMD)
        CALL adiab_print_nacv(adiab_section, nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
     end if

     !> INITIALIZE  
     IF  (.NOT. my_test_keyword) THEN
         decohere = .FALSE.
         aom_env%decohere = .FALSE.
         aom_env%attempt_hop = .FALSE.
         !if we want to restart wavefunction then
         if (aom_env%restart_key) then
             keyword_init = "DIABATIC_RESTART"
         else
             CALL section_vals_val_get(adiab_section, "INITIALIZATION", c_val=keyword_init)
         end if
         CALL section_vals_val_get(adiab_section, "REPRESENTATION", c_val=keyword_basis_repr)
         SELECT CASE(keyword_init)
         CASE("DIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_DIABAT", &
                                       i_val=first_diabat)
             CALL choose_firstadiabat(aom_env%C_tp, first_diabat, sh_print_env, active_state, aom_env) 
             !an active state in needed for the nuclear propagation but the wafunction has to
             !start form an eigenvector of the hamiltonian in the diab_basis
             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(first_diabat) = CMPLX(1.0D0, 0.0D0,dp)  
             CASE("ADIABATIC_BASIS")
             !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
             CALL nacv_first_ts(aom_env, natom, mat, velocities_com_tp) 
!ONLY FOR THE TEST
                 ctin(:) = aom_env%C_tp(first_diabat,:)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("ADIABATIC")
             CALL section_vals_val_get(adiab_section,"FIRST_ADIABAT",&
                                       i_val=active_state)
             SELECT CASE(keyword_basis_repr)
             CASE("DIABATIC_BASIS")
                 ctin(:) = aom_env%C_tp(:, active_state)
             CASE("ADIABATIC_BASIS")
                 !> CALCULATE NACV FOR ADIAB PROPAGATION AT t=0
                 CALL nacv_first_ts(aom_env, natom, mat, velocities_com_tp) 
                 ctin(:) = CMPLX(0.0D0, 0.0D0,dp)
                 !> WARNING: THIS LINE HAS TO BE MODIFIED TO INCLUDE MORE ORBITAL PER MOLECULE
                 ctin(active_state) = CMPLX(1.0D0, 0.0D0,dp)
             CASE DEFAULT
                 CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
             END SELECT
         CASE("DIABATIC_RESTART")
                 ctin(:) = aom_env%readcoeff(:)
                 active_state = aom_env%active_state_res
                 ! do not print this duplicated step in run-sh.log
                 !sh_print_env%first_adiabat = active_state 
                 !CALL adiab_print_sh_info(sh_print_env, keyword = 'START') 
         CASE DEFAULT
             CPABORT("INITIALIZATION implemented: DIABATIC, ADIABATIC")
         END SELECT
!         times = 0.0

         allocate(adiab_c(mat)) 
         allocate(adiab_pop(mat))

         !> print info in the output only if we do not restart the file
         if (.NOT. aom_env%restart_key) then
             !NB: this subroutine is wasted computational work if adiab_pop is not
             !printed out, same probably for AOM.main

             unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIAB_POPULATIONS",&
               extension=".xyz", middle_name="adiab_pop",log_filename=.FALSE.)
             if (unit_nr > 0) then 
                CALL diab_to_adiab( aom_env%C_tp, ctin, adiab_c, adiab_pop)     
                CALL adiab_print_adiab_pop(adiab_section, adiab_pop, & 
                                  nadiab,norbitals,itimes,deltatMD, times)
             end if 
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, times)
             CALL adiab_print_pop(adiab_section, abs(ctin),nadiab,norbitals,itimes,deltatMD, times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes, &
                                        deltatMD, times)
         end if
     ELSE 
         aom_env%attempt_hop = .TRUE.
         active_state = 1
         new_state = mod(active_state + 1, mat) 
         IF ((new_state).EQ.0) new_state = mat
     ENDIF 


     !> second part CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE
     !RESULT OF THE SURFACE HOP ROUTINE IS.
     print*, "FOBSH| ACTIVE STATE = ", active_state
     call calcforces(state=active_state,Hin=-AIMAG(aom_env%H_tp), C=aom_env%C_tp, nacv=nacv, atoms=atoms, & 
                     scaling=scaling, cbar=cbar, Fin=F, forces_adiab=aom_env%forces_adiab, &
                     ext_pres=aom_env%ext_pres)

     IF (analytics) THEN
        times = DBLE(itimes)
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF

     CALL section_vals_val_get(adiab_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
     SELECT CASE(prop_wf_keyword)                                 
     CASE("DIABATIC_RK")
     CASE("DIABATIC_RK_FAST")
         aom_env%no_diab_nace = .True.
         DEALLOCATE(aom_env%Sttp)
         DEALLOCATE(aom_env%Stpt)
         DEALLOCATE(aom_env%Stmt)
         DEALLOCATE(aom_env%Sttm)
         DEALLOCATE(aom_env%Stmtm)
     CASE("DIABATIC_RK_DLI")
     CASE("ADIABATIC_RK")
     CASE("HS_CORRECTED_RK")  
          !> CALC ARRAY OF CLASSICALLY FORBIDDEN HOPS
          !CALL classically_forbidden_hops(active_state, aom_env, mat,  &
          !        velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
          !                atomic_kind_set, class_forbidden_hop_tp)
          !aom_env%class_forbidden_hop_t = class_forbidden_hop_tp
     CASE DEFAULT                                           
           CPABORT("There is no method to propagate the wavefunction.") 
     END SELECT                                             
     
print *, "NO_NACE", aom_env%no_diab_nace
!TEST PRINT DIAB NACE
print *, "DIAB NACE TO STORE FROM NOW ON"
CALL print_matrix_standard_out(aom_env%diabatic_nace_tp)


     !> STORE THE COORDINATES IF NOT SUCCESSFUL RUN
     aom_env%rm_tm      = positions_com_tp - velocities_com_tp*deltatMD
     aom_env%rm_t       = positions_com_tp
     aom_env%sh_print_env = sh_print_env
     aom_env%ctin        = ctin
     aom_env%ctin_t        = ctin
     aom_env%H_t = aom_env%H_tp 
     aom_env%diabatic_nace_t = aom_env%diabatic_nace_tp
     aom_env%E_t = aom_env%E_tp
     aom_env%C_t = aom_env%C_tp
     aom_env%nacv        = nacv
     aom_env%new_state = new_state
     aom_env%active_state = active_state
     !> IMPORTANT FOR THE FIRST TIME STEP
     aom_env%Stt = Stptp
     aom_env%old_state = aom_env%active_state
     !aom_env%active_state = -1
     aom_env%new_state = -1
     aom_env%successful_hop = .FALSE.
     !aom_env%Stptp = Stptp
     aom_env%nacv = nacv
     aom_env%NN = NN


    !> DEALLOCATE
    deallocate(scoefftot)
    deallocate(pcoeffsrun)
    deallocate(ctin)
    deallocate(siteenergies)
    deallocate(crds)
    deallocate(coms)
    deallocate(NN)
    deallocate(Atp)
    deallocate(Btp)
    deallocate(Stptp)
    deallocate(nacv, dist)
    deallocate(exact_forces)
    deallocate(F)
    deallocate(site_energies, couplings, adiab_c, adiab_pop)
    DEALLOCATE(positions_com_tp)
    DEALLOCATE(velocities_com_tp)

  END SUBROUTINE propagate_aom_namd_first_step

  !Seems to change the nace not the nacv 
  SUBROUTINE nacv_first_ts(aom_env, natom, mat, velocities_com_tp)

     TYPE(aom_env_type), INTENT(INOUT)              :: aom_env
     INTEGER, intent(in)                            :: natom
     integer, intent(in)                            :: mat
     REAL(kind=dp), DIMENSION(:,:), INTENT(IN)      :: velocities_com_tp

     CHARACTER(len=*), PARAMETER :: routineN ='nacv_first_ts', &
     routineP = moduleN//':'//routineN
     
     integer             :: i,j, iparticle, coord
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv

     aom_env%adiab_nace_tp = 0.0D0
     do i =1, mat
         do j = 1, mat
              if (j .ge. i) then
                  CALL calc_hop_nacv(aom_env, hop_nacv, i, j, fast=.TRUE.)
                  do iparticle=1,natom
                       do coord = 1,3
                            aom_env%adiab_nace_tp(i,j) = aom_env%adiab_nace_tp(i,j) + &
                                                    velocities_com_tp(coord, iparticle)* & 
                                                            hop_nacv(coord, iparticle)
                       end do
                  end do
              else
                  aom_env%adiab_nace_tp(i,j)= -aom_env%adiab_nace_tp(j,i)
              end if 
         end do
     end do
     aom_env%adiab_nace_tp = -aom_env%adiab_nace_tp 
  END SUBROUTINE nacv_first_ts

END MODULE aom_initialization
