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

MODULE aom_main 

     USE kinds,                        ONLY: default_string_length,&
                                             dp


     USE aom_overlapfinal,             ONLY:  connect_list2, &
                                              calc_sab, &
                                              sab_screen_new, &
                                              calc_dRSab, &
                                              calc_rvecs, &
                                              calc_sab_pete

     USE aom_utils,                    ONLY: makeX,    &
                                             make_adiabX,      &
                                             adiabaticity,   &
                                             surf_hop_type,     &
                                             calcforces,     &
                                             read_listcoeff, &
                                             choose_firstadiabat, &
                                             simple_vel_rescaling, &
                                             translate_com, &
                                             calc_adiab_nace, &
                                             calc_adiab_nace_old, &
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
                                             adiab_print_adiab_pop,                &
                                             adiab_print_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_site_energies,   &
                                             adiab_print_off_diagonals,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian, &
                                             adiab_print_pseudo_hamiltonian, &
                                             adiab_print_nacv, &
                                             adiab_print_pvecs



     USE atomic_kind_list_types,       ONLY: atomic_kind_list_type

     USE atomic_kind_types,            ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set

     USE distribution_1d_types,        ONLY: distribution_1d_type

     USE non_adiabatic_utils,            ONLY: decohere_wf,  &
                                               propagate_wf, &
                                               propagate_wf_fast, &
                                               propagate_adiab_wf, &
                                               propagate_wf_diab_nace_interpol, &
                                               prop_wf_HS_corr, & 
                                               classically_forbidden_hops, &
                                               define_active_region, &
                                               calc_psi_after_correction

     USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream

     USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                                cp_logger_type

     USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                                cp_print_key_should_output


#include "./base/base_uses.f90"


     IMPLICIT NONE

        PRIVATE
 
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='aom_main'   

        PUBLIC :: construct_ham_X_SH

CONTAINS

    SUBROUTINE construct_ham_X_SH(nadiab, aom_atoms, norbitals,  &
                                  positions_tp,  &
                                  energies_tp,  &
                                  aom_env, itimes, &
                                  global_forces, coupling_section, &
                                  adiab_section, my_bo_keyword,                 & 
                                  atomic_kinds, local_particles,                &
                                  my_test_keyword)


     INTEGER, INTENT(IN)                            :: nadiab,  norbitals 
     INTEGER, INTENT(IN)                            :: aom_atoms
     REAL(kind=dp), DIMENSION(:,:),POINTER          :: positions_tp
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
     INTEGER                                    :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), POINTER                       :: atomic_kind_set

     CHARACTER(len=*), PARAMETER :: routineN ='construct_ham_X_SH', &
     routineP = moduleN//':'//routineN

     integer :: coord
     integer :: mol, atoms, c,  mat, orb
     REAL(KIND=dp) :: cutoff_sites, cutoff_connect, scaling, deltatMD, deltatE
     REAL(KIND=dp)    :: cbar, minao
     INTEGER, dimension(:), allocatable          :: class_forbidden_hop_tp
     LOGICAL                                      :: analytics, do_nace_speedup
     INTEGER                                      :: ELsteps   
     REAL(KIND=dp)                                :: dR, times 
     REAL(KIND=dp)                                :: decohere_criterion
     INTEGER                                      :: imol, iatom, natom
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: exact_forces
     REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: site_energies, couplings
     LOGICAL                                      :: keyword_translate
     TYPE(sh_print_env_type)                      :: sh_print_env
     REAL(kind=dp), DIMENSION(:,:),POINTER        :: positions_com_tp, velocities_com_tp
     CHARACTER(LEN=default_string_length)         :: prop_wf_keyword 
     CHARACTER(LEN=default_string_length)         :: reversal_keyword
     real(kind=dp), DIMENSION(:,:), allocatable   :: PC_allrvecs
     LOGICAL                                      :: do_random 
     REAL(KIND=dp), dimension(:), allocatable     :: adiab_pop
     COMPLEX(KIND=dp), dimension(:), allocatable               :: adiab_c
     TYPE(cp_logger_type), POINTER                :: logger
     INTEGER                                    :: unit_nr
!TEST
!     INTEGER                                    :: d, i,j

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
     CALL section_vals_val_get(adiab_section, "DO_SPEEDUP_NACE", l_val=do_nace_speedup)
     CALL section_vals_val_get(adiab_section,"MULTIPLE_TIME_STEP",&
             i_val=aom_env%multiple_time_step)

!TO READ FROM EXTERNAL FILE
!     aom_env%multiple_time_step = 10 
!print *, "MULTIPLE TIME STEP", aom_env%multiple_time_step
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

     sh_print_env = aom_env%sh_print_env
     times = DBLE(itimes)
     sh_print_env%times = times
     
     ALLOCATE(positions_com_tp(3, natom))
     ALLOCATE(velocities_com_tp(3, natom))
     allocate(exact_forces(3,natom))

     !> OBTAIN COORDINATES, ENERGIES, STATES, etc. FROM PREVIOUS STEP
     CALL translate_com(positions_tp, keyword_translate, positions_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     CALL translate_com(aom_env%velocities_tp, keyword_translate, velocities_com_tp, local_particles, atomic_kinds = atomic_kinds)  
     aom_env%siteenergies = 0.0_dp
     do coord=1,3
        aom_env%crds(coord,:,1) = PACK(aom_env%rm_tm(coord,:), aom_env%ext_pres)
        aom_env%crds(coord,:,2) = PACK(aom_env%rm_t(coord,:), aom_env%ext_pres)
        aom_env%crds(coord,:,3) = PACK(positions_com_tp(coord,:), aom_env%ext_pres)
        aom_env%siteenergies(3,:) = energies_tp(:) ! CHANGE_AC
     enddo


     !ctin(:) = aom_env%ctin(:)
     !> store adiabatic nace and diabatic non-orthogonal nacv (aom_env%nacv) from the previous
     !step
     aom_env%adiab_nace_t = aom_env%adiab_nace_tp
!     aom_env%nacv_t = aom_env%nacv
!print *, "NACV from prev step", aom_env%nacv_t(1,:,:) 

     !> CALCULATE THE CENTER OF MAS AND CONSTRUCT NEIGHBOURS LIST
     CALL calc_coms_and_nlist(aom_env, atoms, mol, cutoff_sites, orb)

     !> CALCULATE THE HAMILTONIAN, THE STTP MATRIX AND THE NACV
     CALL calc_Xmatrix_elements(aom_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb, minao)

!print *, "Stptp"
!CALL print_matrix_standard_out(aom_env%Stptp)
!print *, "Stt"
!CALL print_matrix_standard_out(aom_env%Stt)

     !> CONVERT THE diabatic NACV WITH THE GOOD UNITS
!print *, "NACV from inside", aom_env%nacv(1,:,:)
    if (mod(itimes-1, aom_env%multiple_time_step) .eq. 0) then
        print *, "FOBSH| *** OFF-DIAGONAL GRADIENT CALCULATION"
        aom_env%nacv = cbar*aom_env%nacv
    end if
!TEST NACV ABSENCE:
!aom_env%nacv = 0.0D0
!print *, "NACV tp", aom_env%nacv(1,:,:)


 !THIS KEYWORD helps treating initial itimes correctly in case of restart file
 !usage 
     aom_env%restart_key = .False.

     !> HAMILTONIAN DIAGONALIZATION
     !The adiabaticity subroutine diagonalises the diabatic Hamiltonian H_tp
     !It gives us back the eigenvalues E_tp and eigenvectors C_tp
     !initialize Adiabatic energies that will be used for trivial hopping and
     !decoherence corrections
     CALL adiabaticity(-AIMAG(aom_env%H_tp), aom_env%C_tp, aom_env%E_tp)
     !Calculate the overlap between C_tp and C_t, reorder state if required,
     !calculate the nace
     IF (do_nace_speedup) THEN
        CALL calc_adiab_nace(aom_env, adiab_section)
     ELSE
        CALL calc_adiab_nace_old(aom_env, adiab_section)
     END IF


!>this if statement is usefull at t=1 becuase aom_env%diabatic_nace_tp
!> is not present and we would need it in the interpolation scheme during the
!> wavefunction propagation, so no we approximate that aom_env%diabatic_nace_t =0
!> and aom_env%diabatic_nace_tp = aom_env%diabatic_nace_t
     if (itimes.eq.1 .or. aom_env%restart_key)  then
       aom_env%diabatic_nace_tp = aom_env%diabatic_nace_t
     end if
     

     !> PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
     !>the CHARGE CARRIER WF is called aom_env%ctin
     aom_env%HS_propagation = .FALSE.
     CALL section_vals_val_get(adiab_section,"RK_PROPAGATION", c_val=prop_wf_keyword)
     SELECT CASE(prop_wf_keyword)                                 
     CASE("DIABATIC_RK")
           CALL propagate_wf(aom_env%ctin, aom_env, nadiab,  norbitals, ELsteps, &
                             itimes, adiab_section, deltatMD, deltatE , times, mat)
     CASE("DIABATIC_RK_FAST")
           CALL propagate_wf_fast(aom_env%ctin, aom_env, nadiab,  norbitals, ELsteps, &
                             itimes, adiab_section, deltatMD, deltatE , times, mat)
     CASE("DIABATIC_RK_DLI")
           CALL propagate_wf_diab_nace_interpol(aom_env%ctin, aom_env, nadiab,  norbitals, ELsteps, &
                             itimes, adiab_section, deltatMD, deltatE , times, mat)
     CASE("ADIABATIC_RK")
           CALL propagate_adiab_wf(aom_env%E_t, aom_env%E_tp, aom_env%ctin, aom_env, nadiab, norbitals, ELsteps, &
                                   itimes, adiab_section, deltatMD, deltatE , times,  mat)
     CASE("HS_CORRECTED_RK")  
           !> ADIABATIC NACV EVALUATION FOR HS CORRECTION
           !> class_forbidden_hop_tp is calculated refferring to tp to be
           !consistent with the hopping probability, which is evaluated at tp.
           !These transition will be set to zero in the hamiltonian evolution
           !between t and tp.
           CALL classically_forbidden_hops(aom_env%active_state, aom_env, mat,  &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, class_forbidden_hop_tp)
           CALL prop_wf_HS_corr(aom_env%ctin, aom_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times, &
                            mat, class_forbidden_hop_tp, aom_env%active_state)
           CALL section_vals_val_get(coupling_section,"METHOD_REVERSAL",c_val=reversal_keyword)
           if (reversal_keyword .eq. "ALWAYS") then
                aom_env%HS_propagation = .FALSE.
           else 
                aom_env%HS_propagation = .TRUE.
           end if
     CASE DEFAULT                                           
           CPABORT("There is no method to propagate the wavefunction.") 
     END SELECT                                             


     !> SURFACE HOPPING TEST
     !The eigenvectors are needed to calculate the surface hopping probabilities
     !The surfacehop subroutine calculates these probabilites and decides on the
     !new (possibly unchanged) surface for simulation
     aom_env%decohere = .FALSE.
     aom_env%attempt_hop = .FALSE.
     IF (my_test_keyword) THEN   
        aom_env%attempt_hop = .TRUE.
        aom_env%active_state = 1
        aom_env%new_state = mod(aom_env%active_state + 1, mat) 
        IF ((aom_env%new_state).EQ.0) aom_env%new_state = mat
     ELSE IF (.NOT.my_bo_keyword) THEN ! THEN IT IS A FSSH RUN
            call surf_hop_type(adiab_section, aom_env%C_t, aom_env%C_tp,   &
                                 deltatMD, &
                                 aom_env%new_state, &
                                 aom_env%attempt_hop, &
                                 decohere_criterion, &
                                 sh_print_env, aom_env)
     ENDIF

     !CALCULATE THE FORCES IN THE ACTIVE STATE, WHATEVER THE RESULT OF THE SURFACE
     !HOP ROUTINE IS. 
     !DO imol=1,mol
     !   DO iatom=1, natom
     !      !aom_env%F(:,iatom,imol) =  global_forces(imol)%forces(:,iatom)
     !      aom_env%F(:,iatom,imol) =  aom_env%diagonal_F(:,iatom,imol)
     !   ENDDO
     !ENDDO
print *, "AOM_ENV_FORCES"
DO imol=1,mol
    DO iatom=1, natom
      print *, aom_env%F(1,iatom,imol)
    ENDDO
ENDDO

     print*, "FOBSH| ACTIVE STATE = ", aom_env%active_state


!TEST PRINTING DIAGONAL GRADIENTS 
!do d=1,3 
!   do j=1,natom
!     do i=1,mat
!print *, "DIAGONAL FORCE", aom_env%F(d,j,i)
!     end do
!   end do
!end do
!print *, "Fdiag", aom_env%F(1,:,:)

!print *, "TOTAL_NUM_ATOMS", natom
!print *, "aom_atoms", aom_atoms
!print *, "mat", mat
!print *, "nacvx", aom_env%nacv(1,:,:)
!print *, "nacvy", aom_env%nacv(2,:,:)
!print *, "nacvz", aom_env%nacv(3,:,:)

     call calcforces(state=aom_env%active_state,Hin=-AIMAG(aom_env%H_tp), C=aom_env%C_tp, nacv=aom_env%nacv, &
                           atoms=atoms, & 
                           scaling=scaling, cbar=cbar, Fin=aom_env%F, forces_adiab=aom_env%forces_adiab, &
                           ext_pres=aom_env%ext_pres)
     IF (aom_env%attempt_hop) THEN
         call calcforces(state=aom_env%new_state,Hin=-AIMAG(aom_env%H_tp), C=aom_env%C_tp, nacv=aom_env%nacv, &
                              atoms=atoms, & 
                              scaling=scaling, cbar=cbar, Fin=aom_env%F, forces_adiab=aom_env%forces_adiab_new, &
                              ext_pres=aom_env%ext_pres)
     ENDIF


     !> PRINT
     CALL adiab_print_hamiltonian(adiab_section, -AIMAG(aom_env%H_tp), nadiab,norbitals,itimes, deltatMD)
     CALL adiab_print_pseudo_hamiltonian(adiab_section, -AIMAG(aom_env%H_tp), nadiab,norbitals,itimes, deltatMD)

     !> print site energies 
     allocate(site_energies(mat))
     allocate(couplings(mat-1))
     allocate(adiab_c(mat))
     allocate(adiab_pop(mat))
 
     do imol=1,mat
         site_energies(imol) = -AIMAG(aom_env%H_tp(imol,imol))
         if (imol .ne. mat) then
           couplings(imol) = -AIMAG(aom_env%H_tp(imol,imol+1))
         end if 
     end do


     unit_nr=cp_print_key_unit_nr(logger,adiab_section,"PRINT%ADIAB_POPULATIONS",&
               extension=".xyz", middle_name="adiab_pop",log_filename=.FALSE.)
     if (unit_nr > 0) then
        !CALL rotation_diab_to_adiab(aom_env%C_tp, aom_env%ctin, mat, adiab_ctin, adiab_pop)
        CALL diab_to_adiab( aom_env%C_tp, aom_env%ctin, adiab_c, adiab_pop)
        CALL adiab_print_adiab_pop(adiab_section, adiab_pop, nadiab,norbitals,itimes, deltatMD, times) 
     end if

     CALL adiab_print_site_energies(adiab_section, site_energies, nadiab, deltatMD, times)
     CALL adiab_print_off_diagonals(adiab_section, couplings, nadiab-1, deltatMD, times)
     CALL adiab_print_nacv(adiab_section, aom_env%nacv, nadiab, norbitals, aom_atoms,itimes, deltatMD)
     CALL adiab_print_adiabat_energies(adiab_section, aom_env%E_tp, nadiab, deltatMD, times)
     CALL adiab_print_pvecs(adiab_section, PC_allrvecs, itimes, deltatMD, times, atoms, nadiab)
     IF (analytics) THEN
        CALL adiab_print_exact_forces(adiab_section, exact_forces, natom, deltatMD, times)
     ENDIF

     !> TRANSLATE THE COORDINATE FOR ONE TIMESTEP
     aom_env%rm_tm      = aom_env%rm_t
     aom_env%rm_t       = positions_com_tp
     aom_env%sh_print_env = sh_print_env

     !> DEALLOCATE
     DEALLOCATE(positions_com_tp) 
     DEALLOCATE(velocities_com_tp)
     deallocate(exact_forces)
     deallocate(site_energies, couplings, adiab_c, adiab_pop)
    END SUBROUTINE construct_ham_X_SH


SUBROUTINE calc_coms_and_nlist(aom_env, atoms, mol, cutoff_sites, orb)
     TYPE(aom_env_type), INTENT(INOUT)     :: aom_env
     INTEGER, INTENT(IN)                   :: mol, atoms, orb
     REAL(KIND=dp), intent (IN)            :: cutoff_sites

     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_coms_and_nlist', &
     routineP = moduleN//':'//routineN

     integer                               :: i, j, m, d
    
     !> CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     aom_env%coms = 0
     do i=1,mol
          do m=1,atoms
               do d=1,3
                    aom_env%coms(d,i) = aom_env%coms(d,i) + (aom_env%crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do

     !> CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               aom_env%dist(i,j) = ((aom_env%coms(1,i) - aom_env%coms(1,j))*(aom_env%coms(1,i) - &
                                 aom_env%coms(1,j)) + (aom_env%coms(2,i) - aom_env%coms(2,j))*&
                      (aom_env%coms(2,i) - aom_env%coms(2,j)) + (aom_env%coms(3,i) - &
                          aom_env%coms(3,j))*(aom_env%coms(3,i) - aom_env%coms(3,j)))
               if (aom_env%dist(i,j) < cutoff_sites*cutoff_sites) then
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 1
                         end do
                    end do
               else
                    do d=1,orb
                         do m=1,orb
                              aom_env%NN(orb*(i-1)+d,orb*(j-1)+m) = 0
                         end do
                    end do
               end if                    
          end do
     end do
END SUBROUTINE calc_coms_and_nlist


SUBROUTINE calc_Xmatrix_elements(aom_env, itimes, atoms, mol, dR, do_random, cbar, &
                      scaling, PC_allrvecs, orb, minao)

     TYPE(aom_env_type)                      :: aom_env
     INTEGER, POINTER                        :: itimes
     integer, intent(IN)                     :: atoms
     integer, intent(in)                     :: mol, orb
     REAL(KIND=dp), intent(in)                   :: dR  
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in)                   :: minao
     REAL(KIND=dp), intent(in) :: scaling
     real(kind=dp), DIMENSION(:,:), allocatable, intent(out) :: PC_allrvecs

     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_Xmatrix_elements', &
     routineP = moduleN//':'//routineN

     REAL(KIND=dp)                   :: S
     integer :: i, j, k, l, m, d, o
     real(kind=dp), allocatable :: PC_rvecs(:,:)
      
     ALLOCATE(PC_allrvecs(5,mol*atoms))
     ALLOCATE(PC_rvecs(3,atoms))
!COMMENT: I guess Atm and Btm are only need for the first timestep (itimes = 1)
!one could duplicate this subroutine and create one only for the first timestep
     !for each molecule i
!print *, "aom_env%NN", aom_env%NN
     do i=1,mol
!print *, "iiiii", i
               !create the array of all its atomic coordinates at t
               !numerical factor converts from Angstroms into Bohr radii
               do j=1,atoms
                    aom_env%At(1,j)  = dble(aom_env%elements(j))
                    aom_env%Atp(1,j) = dble(aom_env%elements(j))
                    aom_env%Atm(1,j) = dble(aom_env%elements(j))
                    do d=2,4
                         !> CHANGE_AC : already coordonates in atomic units
                         aom_env%Atm(d,j) = aom_env%crds(d-1,(i-1)*atoms+j,1)
                         aom_env%At(d,j)  = aom_env%crds(d-1,(i-1)*atoms+j,2)
                         aom_env%Atp(d,j) = aom_env%crds(d-1,(i-1)*atoms+j,3)
                    end do
               end do

               !then for each other molecule
               do k=1,mol
!print *, "k", k
                    !create the array of its atomic coordinates at both t and t+dt
                    do l=1,atoms
                         aom_env%Bt(1,l)  = dble(aom_env%elements(l))
                         aom_env%Btp(1,l) = dble(aom_env%elements(l))
                         aom_env%Btm(1,l) = dble(aom_env%elements(l))
                         do d=2,4
                              !> CHANGE_AC : already coordonates in atomic units
                              aom_env%Btm(d,l) = aom_env%crds(d-1,(k-1)*atoms+l,1)
                              aom_env%Bt(d,l)  = aom_env%crds(d-1,(k-1)*atoms+l,2)
                              aom_env%Btp(d,l) = aom_env%crds(d-1,(k-1)*atoms+l,3)
                         end do
                    end do
                    if (i == 1) then
                         do d=1,orb
!print *, "d", d
                              aom_env%pcoeffsrun(:,orb*(k-1)+d,1) = aom_env%pcoeffs(:,d)
!print *, "call calc_sab 1"
                              call calc_sab(aom_env%Btm, aom_env%connlist2, &
                                                 aom_env%pcoeffsrun(:,orb*(k-1)+d,1))

                              aom_env%pcoeffsrun(:,orb*(k-1)+d,2) = aom_env%pcoeffs(:,d)
!print *, "call calc_sab 2"
                              call calc_sab(aom_env%Bt, aom_env%connlist2, &
                                                 aom_env%pcoeffsrun(:,orb*(k-1)+d,2))

                              aom_env%pcoeffsrun(:,orb*(k-1)+d,3) = aom_env%pcoeffs(:,d)
!print *, "call calc_sab 3"
                              call calc_sab(aom_env%Btp, aom_env%connlist2, & 
                                               aom_env%pcoeffsrun(:,orb*(k-1)+d,3))
                         end do
                    end if
                    !neighbourlist(:,:,:,:) = aom_env%neighbourlist(:,:,:,:)

                    !!recalculate neighbouring lists from scratch
                    !! you could even reduce the times we calculate neighbouring
                    !list even further
                    if (mod(itimes-1, aom_env%multiple_time_step) .eq. 0 & 
                          .and. aom_env%multiple_time_step .ne. 1) then 
                       do d=1,orb
                            do m=1,orb
                                 call sab_screen_new(aom_env%Atp, aom_env%atomlist, &
                                  aom_env%connlist1, aom_env%pcoeffsrun(:,orb*(i-1)+d,3),&
                                 & aom_env%atomlist, aom_env%Btp, aom_env%atomlist, &
                                  aom_env%connlist1, aom_env%pcoeffsrun(:,orb*(k-1)+m,3),&
                                 & aom_env%aneighbour(:,1:atoms), aom_env%bneighbour(:,1:atoms), & 
                                                                                 minao)
                            end do
                       end do
                     end if

                    !only call overlap if neighbours list says so
                    if (aom_env%NN(i,k) == 1) then
                         !first the diagonals
                         if (i == k) then
                              CALL diagonal_elem(i,k,orb,itimes,aom_env)
                         !now the off-diagonals                         
                         else ! ( i == k)
                              CALL off_diagonal_elem(aom_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)
                         end if !(i==k)
                    else !(NN(i,k) == 1
                         do d=1,orb
                              do m=1,orb
                                   !fill in elements of non-overlapping pairs with 0
                                   aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0, 0.0D0,dp)
                                   aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 0.0D0
                                   aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
     ! in case no_diab_nace is true nace in the diabatic basis are not used
     ! neither to
     ! propagate_wf nor in the adiabatic nace calculation for the hopping
     ! probability
     ! so that all the Sttp and co matrices are also not used and not calculated
     ! NB probably the same could be done for aom_env%Stt
                                   if (.NOT. aom_env%no_diab_nace) then
                                       aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 0.0D0
                                       aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                       aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = 0.0D0
                                   end if
                              end do
                         end do
                    end if !(NN(i,k) == 1

               end do !loop over k
               call calc_sab_pete(aom_env%Atp, aom_env%connlist2, PC_rvecs)
               do o=1, size(PC_rvecs,2)
                  PC_allrvecs(1,((i-1)*atoms)+o) = i
                  PC_allrvecs(2,((i-1)*atoms)+o) = o
                  PC_allrvecs (3:5,((i-1)*atoms)+o) = PC_rvecs(:,o)
               end do
     end do !loop over i
END SUBROUTINE calc_Xmatrix_elements

SUBROUTINE diagonal_elem(i,k,orb,itimes,aom_env)

     integer, intent(IN) :: i, k,orb
     INTEGER, POINTER                               :: itimes
     TYPE(aom_env_type)              :: aom_env

     CHARACTER(len=*), PARAMETER :: routineN ='diagonal_elem', &
     routineP = moduleN//':'//routineN

     integer ::  m, d

     do d=1,orb
          do m=1,orb
               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
               aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
     ! in case no_diab_nace is true nace in the diabatic basis are not used
     ! neither to
     ! propagate_wf nor in the adiabatic nace calculation for the hopping
     ! probability
     ! so that all the Sttp and co matrices are also not used and not calculated
     ! NB probably the same could be done for aom_env%Stt
               if (.NOT. aom_env%no_diab_nace) then
                  aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                  aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
               end if
               aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0, &
                                         (-aom_env%siteenergies(3,k)),dp)
          end do
     end do
     if (itimes.eq.1 .or. aom_env%restart_key)  then
          do d=1,orb
               do m=1,orb
               !These diagonals are 1 by definition
               !A bit of math to be convinced for Sttp,
               !Stpt: dSii/dt = 0
                    aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = 1.0D0
                    if (.NOT. aom_env%no_diab_nace) then
                      aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                      aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = 1.0D0
                      aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = 1.0D0
                    end if
               end do
          end do     
     end if ! (.not.inrun)
END SUBROUTINE diagonal_elem

SUBROUTINE off_diagonal_elem(aom_env,i,k,orb,itimes, atoms, &
                                   dR, S, do_random, cbar, scaling)

     TYPE(aom_env_type)                           :: aom_env
     integer, intent(IN) :: i, k,orb, atoms
     INTEGER, POINTER                               :: itimes
     REAL(KIND=dp), intent(in)                   :: dR  
     REAL(KIND=dp), intent(out)                  :: S
     LOGICAL, intent(in)                         :: do_random
     REAL(KIND=dp), intent(in)                   :: cbar
     REAL(KIND=dp), intent(in) :: scaling

     CHARACTER(len=*), PARAMETER :: routineN ='off_diagonal_elem', &
     routineP = moduleN//':'//routineN

     integer ::  m, d
     !symmetric matrices
     !no sense calling calc_sab twice as often as you need to
     if (k < i) then
          do m=1,orb
               do d=1,orb
                    aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = aom_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stt(orb*(i-1)+d,orb*(k-1)+m)   = &
                                             aom_env%Stt(orb*(k-1)+m,orb*(i-1)+d)
                    aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = &
                                            aom_env%Stptp(orb*(k-1)+m,orb*(i-1)+d)
     ! in case no_diab_nace is true nace in the diabatic basis are not used
     ! neither to
     ! propagate_wf nor in the adiabatic nace calculation for the hopping
     ! probability
     ! so that all the Sttp and co matrices are also not used and not calculated
     ! NB probably the same could be done for aom_env%Stt
                    if (.NOT. aom_env%no_diab_nace) then
                       aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = &
                                                aom_env%Stmtm(orb*(k-1)+m,orb*(i-1)+d)
                       !A cross symmetry for Sttp, Stpt,
                       !Sttm, Stmt
                       aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                aom_env%Stpt(orb*(k-1)+m,orb*(i-1)+d)
                       aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m)  = & 
                                                 aom_env%Sttp(orb*(k-1)+m,orb*(i-1)+d)
                       if (itimes.eq.1 .or. aom_env%restart_key)  then
                          aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                 aom_env%Stmt(orb*(k-1)+m,orb*(i-1)+d)
                          aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m)  = &
                                                 aom_env%Sttm(orb*(k-1)+m,orb*(i-1)+d)
                       ENDIF
                    end if
               end do
          end do
     else
          do d=1,orb
               do m=1,orb
               !PC added print_rvecs =yes so p vectors and coords would be printed 
!print *, "call calc_sab 4"
                    call calc_sab(aom_env%Atp, aom_env%connlist2, &
                        aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3),&
                     S, aom_env%Btp, aom_env%connlist2, &
                               aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3),& 
                     print_rvecs='yes', molA_index=i, molB_index=k)
                    aom_env%Stptp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                    aom_env%H_tp(orb*(i-1)+d,orb*(k-1)+m) = cmplx(0.0D0, -scaling*S,dp)

!print *, "call calc_dRsab 1"
!print *, "REST", mod(itimes-1, aom_env%multiple_time_step)
                    !> In the following diabatic non-orthog. nace are
                    !calculated, since calculating them is very time consuming
                    !and expensive we opted for a multiple time step scheme.
                    !while the nacv are within the counter(chosen in the input 
                    !file) they are considered a
                    !constant and taken from the prevous time step.
                    !Calcul of the NACV
                    if (mod(itimes-1, aom_env%multiple_time_step) .eq. 0 &
                          .and. aom_env%multiple_time_step .ne. 1) then 
                        call calc_dRSab(aom_env%Atp, aom_env%atomlist, &
                                      aom_env%connlist1, &
                        aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), aom_env%scoeffs, &
                        aom_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                        aom_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                        aom_env%Btp, aom_env%atomlist, aom_env%connlist1, &
                        aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3), &
                         !aom_env%scoeffs, aom_env%neighbourlist(k,i,:,:), &
                         !         aom_env%neighbourlist(i,k,:,:), dR, &
                         aom_env%scoeffs, aom_env%aneighbour(:,:), &
                                          aom_env%bneighbour(:,:), dR, &
                         do_random, aom_env)
                        !> note that calc_dRSab uses a radom number to randomize
                        !forwards derivative in the NACV (this improves energy
                        !conservation by one/two order of magnitude. To be
                        !consistent with random sequence, one should:
                       !else
                             !call random generator function

                   else if (aom_env%multiple_time_step .eq. 1) then 
                        call calc_dRSab(aom_env%Atp, aom_env%atomlist, &
                                      aom_env%connlist1, &
                        aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), aom_env%scoeffs, &
                        aom_env%nacv(1:3,orb*(k-1)+d,(i-1)*atoms+1:i*atoms), &
                        aom_env%nacv(1:3,orb*(i-1)+m,(k-1)*atoms+1:k*atoms), &
                        aom_env%Btp, aom_env%atomlist, aom_env%connlist1, &
                        aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3), &
                         aom_env%scoeffs, aom_env%neighbourlist(k,i,:,:), &
                                  aom_env%neighbourlist(i,k,:,:), dR, &
                         do_random, aom_env)
                        !> note that calc_dRSab uses a radom number to randomize
                        !forwards derivative in the NACV (this improves energy
                        !conservation by one/two order of magnitude. To be
                        !consistent with random sequence, one should:
                       !else
                             !call random generator function
          
                    end if
               end do
          end do
     ! in case no_diab_nace is true nace in the diabatic basis are not used
     ! neither to
     ! propagate_wf nor in the adiabatic nace calculation for the hopping
     ! probability
     ! so that all the Sttp and co matrices are also not used and not calculated
     ! NB probably the same could be done for aom_env%Stt and aom_env%Stptp
         if (.NOT. aom_env%no_diab_nace) then
             do d=1,orb
                  do m=1,orb
!    print *, "call calc_sab 5"
                       call calc_sab(aom_env%At, aom_env%connlist2, &
                        aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                        aom_env%Btp,&
                       & aom_env%connlist2, aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,3))
                       aom_env%Sttp(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
!    print *, "call calc_sab 6"
                       call calc_sab(aom_env%Atp, aom_env%connlist2, &
                            aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,3), S,& 
                        aom_env%Bt,&
                       & aom_env%connlist2, aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                       aom_env%Stpt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                  end do
             end do
             if (itimes.eq.1 .or. aom_env%restart_key)  then
                  do d=1,orb
                       do m=1,orb
!    print *, "call calc_sab 7"
                            call calc_sab(aom_env%Atm, aom_env%connlist2, &
                            &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                                aom_env%Btm, aom_env%connlist2,&
                            &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                            aom_env%Stmtm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
    
!    print *, "call calc_sab 8"
                            call calc_sab(aom_env%At, aom_env%connlist2, &
                            &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,2), S, &
                                       aom_env%Btm, aom_env%connlist2,&
                            &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,1))
                            aom_env%Sttm(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
    
!    print *, "call calc_sab 9"
                            call calc_sab(aom_env%Atm, aom_env%connlist2, &
                            &aom_env%pcoeffsrun(1:atoms,orb*(i-1)+d,1), S, &
                              aom_env%Bt, aom_env%connlist2, &
                            &aom_env%pcoeffsrun(1:atoms,orb*(k-1)+m,2))
                            aom_env%Stmt(orb*(i-1)+d,orb*(k-1)+m) = cbar*S
                       end do
                  end do
    
             end if ! (.not.inrun)
         end if 
     end if ! ( k < i)
END SUBROUTINE off_diagonal_elem

END MODULE aom_main
