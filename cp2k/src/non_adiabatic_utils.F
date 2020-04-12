!-----------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations           
!   !                                                                           
!   Copyright (C) 2000 - 2015  CP2K developers group                            
!   !                                                                           
!-----------------------------------------------------------------------------! 

!!!!!!!!!!!!!!!!!!!! CHANGE_SG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE non_adiabatic_utils
     
     USE kinds,                        ONLY: default_string_length, &
                                             dp                                    

     USE input_section_types,          ONLY:  section_vals_type, &
                                             section_vals_get, &
                                             section_vals_get_subs_vals, &
                                             section_vals_list_get, &
                                             section_vals_val_get

     USE aom_print,                    ONLY: adiab_print_coeff,              &
                                             adiab_print_pop,                &
                                             adiab_print_adiabat_energies,   &
                                             adiab_print_matrix_density,     &
                                             adiab_print_exact_forces,       &
                                             adiab_print_active_state,       &
                                             adiab_print_sh_info,            &
                                             adiab_print_hamiltonian,        &
                                             adiab_print_nace,               &
                                             adiab_print_nacv,               &
                                             adiab_print_hop_nacv

     USE parallel_rng_types,           ONLY: next_random_number    

     USE sh_types,                    ONLY: sh_env_type, &      
                                                sh_print_env_type 

     USE aom_utils,                    ONLY: makeX,                & 
                                             adiabaticity,         &
                                             calc_hop_nacv_aom,    &
                                             make_adiabX,          &
                                             calc_diabatic_nace,   & 
                                             matinvsqrt

!     USE non_adiabatic_types,           ONLY: sh_env_type

     USE atomic_kind_list_types,          ONLY: atomic_kind_list_type   

     USE atomic_kind_types,               ONLY: atomic_kind_type,&      
                                           get_atomic_kind,&       
                                           get_atomic_kind_set     
                                                                   
     USE distribution_1d_types,           ONLY: distribution_1d_type    

     USE force_env_types,                ONLY: force_env_get, &
                                            force_env_type



#include "./base/base_uses.f90" 

  IMPLICIT NONE
                                                                                        
        PRIVATE
                                                                                                       
        CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN ='non_adiabatic_utils'
                                                                                                       
        PUBLIC :: decohere_wf,                     &
                  calc_coms_and_nlist,             &
                  calc_tresp_coup,                 &
                  !calc_tresp_coup_test,                 &
                  calc_frozen_H,                   &
                  calcforces_smart_grad,           &
                  !calcforces_smart_grad_test,           &
                  !calcforces_frozenH,              &
                  surf_hop_type,                   &
                  check_hop_attempt,               &
                  choose_firstadiabat,             &
                  print_matrix_standard_out,       &
                  diab_to_adiab,                   &
                  cal_delta_energy_nacv,           &
                  translate_com,                   &  
                  simple_vel_rescaling,            &
                  total_momentum,                  &
                  propagate_wf,                    &
                  propagate_wf_fast,               &
                  propagate_adiab_wf,              &
                  prop_wf_HS_corr,                 &
                  calc_kinetic_tot,                &
                  classically_forbidden_hops,      &
                  ida_collapse,                    &
                  define_active_region,            &
                  calc_psi_after_correction,       &
                  apply_decoherence,               & 
                  propagate_wf_diab_nace_interpol
                  

CONTAINS

!**********************************************************************************
!> \brief it calculates the COMs of the molecules and the distance between these
!         COM points. Using an external input, it determines which COMs are within a
!         given cut-off. 
!> \note this is general for degeneracy and also for different num of atoms per mol
!> \author AC [Edited SG]
!**********************************************************************************
SUBROUTINE calc_coms_and_nlist(mol, cutoff_sites, orb, crds, active_system_info, NN)

     INTEGER, INTENT(IN)                             :: mol, orb
     REAL(KIND=dp), intent (IN)                      :: cutoff_sites
     REAL(KIND=dp), dimension(:,:,:), POINTER, &
                                      INTENT(IN)     :: crds !crds active atoms in COM ref
     INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                      INTENT(IN)     :: active_system_info
     INTEGER, DIMENSION(:,:), POINTER, INTENT(OUT)   :: NN

     CHARACTER(len=*), PARAMETER                     :: routineN='calc_coms_and_nlist', &
     routineP = moduleN//':'//routineN

     integer                                         :: i, j, m, d, atoms
     REAL(KIND=dp)                                   :: dist
     REAL(KIND=dp), dimension(:,:), ALLOCATABLE      :: coms
  
     allocate(coms(3,mol))
    
     !CALCULATE THE CENTER OF MASS
     !calculate CoMs at first timestep
     !coms(1 to 3,i) are the x, y and z coordinates of CoM of i
     coms = 0
     do i=1,mol
          !> counts the active atoms for each molecules to make it general
          atoms = (active_system_info(2,i) - active_system_info(1,i)) + 1
          do m=1,atoms
               do d=1,3
                    coms(d,i) = coms(d,i) + (crds(d,(i-1)*atoms + m,3))/atoms
               end do
          end do
     end do

     !CONSTRUCT NEIGHBOURS LIST
     do i=1,mol
          do j=1,mol
               !> this is like a DOT_PRODUCT to calculate distance
               dist = ((coms(1,i) - coms(1,j))*(coms(1,i) &
                - coms(1,j)) + (coms(2,i) - coms(2,j))*&
               &(coms(2,i) - coms(2,j)) + (coms(3,i) & 
                - coms(3,j))*(coms(3,i) - coms(3,j)))
               if (dist < cutoff_sites*cutoff_sites) then
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

    deallocate(coms)
END SUBROUTINE calc_coms_and_nlist

!**********************************************************************************
!> \brief This calculates the adiabatic nacv according to Eq 11 in Carof 17. It
!          is a generalization of calc_hop_nacv_aom for Hamiltonians 
!          different from the AOM one. It uses only the "a" term in Carof 17
!          paper.
!> \note This subroutine has not been tested for degeneracy, but it should be
!        general enough to handle it.
!> \author SG [taken from calc_hop_nacv_aom of AC]
!**********************************************************************************
SUBROUTINE calc_hop_nacv_approx(sh_env, diagonal_F, C, active_system_info, &
                                     orb, hop_nacv, nacv, convert_nacv, frozen_c, &
                                     ext_old_state, ext_new_state)

        TYPE(sh_env_type)                                        :: sh_env
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE  :: hop_nacv
        REAL(KIND=dp), intent(in), dimension(:,:)                :: C
        REAL(KIND=dp), DIMENSION(:,:,:), INTENT(IN)              :: diagonal_F
        INTEGER, INTENT(IN)                                      :: orb
        INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                      INTENT(IN)                 :: active_system_info
        REAL(KIND=dp), dimension(:,:,:), INTENT(IN), OPTIONAL    :: nacv

        INTEGER, OPTIONAL                                        :: ext_old_state, ext_new_state
        LOGICAL, OPTIONAL, INTENT(IN)                            :: frozen_c
        REAL(KIND=dp), intent(in), OPTIONAL                      :: convert_nacv


        CHARACTER(len=*), PARAMETER :: routineN = 'calc_hop_nacv_approx', &
        routineP = moduleN//':'//routineN

        INTEGER                                                  :: num_mol, num_j, natom
        INTEGER                                                  :: old_state, new_state
        REAL(KIND=dp)                                            :: diff_energy
        REAL(KIND=dp)                                            :: E_old_state, E_new_state
        integer                                                  :: i, j, mat, d,m 
        REAL(KIND=dp), dimension(:,:), allocatable               :: F_1, F_11
        LOGICAL, dimension(:), allocatable                       :: ext_pres
        REAL(KIND=dp)                                            :: my_convert_nacv 
        LOGICAL                                                  :: my_frozen_c

        !COMMENT_SG This subrotine is practically a copy of the calcforces and
        !they could be merged into a single one with some more options.      

        ! Based on paper Carof, Giannini, Blumberger, J. Chem. Phys. 147,
        ! 214113 (2017), [Carof17]

        !The following subroutine may NOT work if orb !=1 (degeneracy). It needs
        !testing.


        natom  = SIZE(diagonal_F, 2)
        mat = SIZE(diagonal_F, 3) !this is actually mol

        my_frozen_c = .FALSE.
        IF (PRESENT(frozen_c)) my_frozen_c = frozen_c

        allocate(ext_pres(natom))
        allocate(F_1(3,natom))
        allocate(F_11(3,natom))
        allocate(hop_nacv(3, natom)) 

        
        !nacv = sh_env%nacv
        !Fin = sh_env%F
        ext_pres = sh_env%ext_pres       
        IF (PRESENT(ext_old_state)) THEN
           old_state = ext_old_state
        ELSE
           old_state = sh_env%active_state
        ENDIF
        IF (PRESENT(ext_new_state)) THEN
           new_state = ext_new_state 
        ELSE
           new_state = sh_env%new_state
        ENDIF
        E_old_state = sh_env%E_tp(old_state)
        E_new_state = sh_env%E_tp(new_state)
        IF (E_old_state.EQ.E_new_state) THEN
           diff_energy = 0.0D0
        ELSE
           diff_energy = 1 / (E_new_state - E_old_state)
        ENDIF  


       ! Only calculation of term a_{I, ij} in Eq. (14) of [Carof17], i.e.
       ! Eq. (15)
       num_mol = 0   ! Molecular index 
       do d=1,3      ! Cartesian coordinates
          num_j = 0  ! Active atom index
          do j=1,natom ! Sum over all the atoms
             F_1(d,j) = 0.0D0
             F_11(d,j) = 0.0D0

             IF (ext_pres(j)) THEN
                num_j   = num_j + 1                  ! Update the active atom index
                ! Get the molecule to which j belong, for a general
                ! system. This is a molecular index needed below. Looping
                ! over molecules is necessary here.
                do i=1,SIZE(active_system_info, 2) ! loop over non-deg states
                   IF ((j .GE. active_system_info(1,i)) .AND. (j .LE. active_system_info(2,i))) THEN
                      num_mol = i
                      EXIT
                   END IF
                end do
             ENDIF
             do i=1, SIZE(active_system_info, 2) !non-deg
                do m=1,orb
                   !F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                   ! Diagonal term, it depends on the forces
                   F_1(d,j) = F_1(d,j) + &
                       C(orb*(i-1)+m,orb*(old_state-1)+m)*diagonal_F(d,j,i)*C(orb*(i-1)+m,orb*(new_state-1)+m)*diff_energy

                   ! calculate off diagonal forces if couplings are not frozen
                   IF (.NOT. my_frozen_c) THEN
                     IF (ext_pres(j)) THEN
                         IF (i.ne.num_mol) THEN
                            IF (PRESENT(convert_nacv)) THEN 
                                !grandient in case of AOM formalism
                                my_convert_nacv = convert_nacv
                            ELSE
                                !gradient for a genaral H
                                my_convert_nacv = 1.0D0  ! no scaling factor since we
                                                      ! use gradient directly from other methods
                            END IF
                         ! Of diagonal term, only for active atom j present in
                         ! molecule num_mol
                         !> Note this equation is a generalization of calcforces
                         !one for two different adiabatic states old and new
                         F_11(d,j) = F_11(d,j) + diff_energy * my_convert_nacv*nacv(d, orb*(i-1)+m, num_j) * (&
                                          C(orb*(num_mol-1)+m, orb*(new_state-1)+m) * C(orb*(i-1)+m, orb*(old_state-1)+m) +   &
                                          C(orb*(num_mol-1)+m, orb*(old_state-1)+m) * C(orb*(i-1)+m,orb*(new_state-1)+m) )
                        
                         ENDIF
                     ENDIF   
                   ENDIF
                end do
             end do
          ENDDO
       ENDDO 

       ! minus sign to take the same convention of the paper.
       hop_nacv = -(F_1 + F_11) 
 
       deallocate(F_1)
       deallocate(F_11)

END SUBROUTINE calc_hop_nacv_approx


!**********************************************************************************
!> \brief calculates coupling using TrESP charges from input. 

!> \param energies_tp, site energies
!> \param orb, number of states per mol (for degeneragy)
!> \param NN, nearest neighbour list to apply cut-off radious to coulomb sum 
!> \param positions_com_tp, coordinates of total system from which we choose
!         active ones.
!> \param active_system_info, initial and final index for the active mols
!> \param tresp_charges, initial tresp read from input (to be generalized to different mols)
!> \note BE CAREFULL: the implementation has not been tested for orb !=1 (namely, degeneragy)
!> \output sh_env%H_tp
!> \author SG
!**********************************************************************************
  subroutine calc_tresp_coup(sh_env, energies_tp, orb, NN, positions_com_tp, active_system_info, tresp_charges)

     TYPE(sh_env_type)                               :: sh_env
     REAL(kind=dp), DIMENSION(:),  POINTER           :: energies_tp
     integer, intent(in)                             :: orb
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, &
                                   INTENT(IN)        :: tresp_charges  !tresp
                                                                       !charges taken from input file, they refer 
                                                                       !to a single mol and should be extended
                                                                       !for multi-phase system
     INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)    :: NN
     REAL(kind=dp), DIMENSION(:,:),POINTER, &
                                   INTENT(IN)        :: positions_com_tp
     INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                   INTENT(IN)        :: active_system_info

     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_tresp_coup', &
     routineP = moduleN//':'//routineN

     integer                                         :: i, k, m, d, dim_
     integer                                         :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
     integer                                         :: init_atom_mol_k, final_atom_mol_k, atoms_per_mol_k
     real(kind=dp), allocatable, dimension(:,:)      :: crd_i, crd_k
     real(kind=dp), allocatable, dimension(:)        :: charges1, charges2
     real(kind=dp)                                   :: tresp_coup
     real(kind=dp),  dimension(:,:), allocatable     :: grad1 !gradients acting on atoms of mol1


     !> COMMENT_SG: Tresp charegs array should be generalized for multiple molecules
     !and phases
     !> COMMENT_SG: This is the only NON general bit of the whole subroutine. It
     !assumes just one type of charges
     !#########  WARNING ##############################################
     ALLOCATE(charges1(SIZE(tresp_charges)))
     ALLOCATE(charges2(SIZE(tresp_charges)))
     charges1 = tresp_charges
     charges2 = tresp_charges
     !#################################################################


     sh_env%H_tp(:,:) = cmplx(0.0D0, 0.0D0,dp)
     sh_env%grad_tresp = 0.0D0

     do i=1,SIZE(energies_tp) !loop over all non-degerate states
        do k=1,SIZE(energies_tp)
           ! do calculation only for molecules within cut-off radious. 
           if (NN(i,k) == 1) then
 
              !> diagonal construction
              if (k==i) then 
                 !since you may want to include multiple orbital per state, you
                 !should loop over orbitals
                 do d=1,orb
                    sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0,(-energies_tp(i)),dp) !minus because we propagate -X
                 end do
              else 
                 !loop over upper diagonal
                 if ( k > i ) then
                    do d=1,orb
                       do m=1,orb
                          !> find coordinates mol 1 of the pair
                          init_atom_mol_i = sh_env%active_system_info(1,i)
                          final_atom_mol_i = sh_env%active_system_info(2,i)
                          atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

                          !> find coordinates mol 2 of the pair
                          init_atom_mol_k = sh_env%active_system_info(1,k)
                          final_atom_mol_k = sh_env%active_system_info(2,k)
                          atoms_per_mol_k = (final_atom_mol_k-init_atom_mol_k)+1


                          !> this allocation deallocation is maybe inefficient
                          !but it was the only thing I found to be general for
                          !multi-phases
                          ALLOCATE(crd_i(3,atoms_per_mol_i))
                          ALLOCATE(crd_k(3,atoms_per_mol_k))

                          do dim_=1,3
                             crd_i(dim_,:) = positions_com_tp(dim_, init_atom_mol_i:final_atom_mol_i)
                             crd_k(dim_,:) = positions_com_tp(dim_, init_atom_mol_k:final_atom_mol_k)
                          end do

!print *, "PAIR", i,k
                          ! calc coupling for a given pair of states
                          CALL calc_trsp_ab(crd_i, crd_k, charges1, charges2, tresp_coup)
                          sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0,-tresp_coup,dp) !sign of the off-diag is minus consisently with diagonals
                          !test to recover frozen coup for a chain:
                          !sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0,0.0011024789,dp)

                          !> For the gradient in the case of TrESP (as it was
                          !for nacv in case of AOM) the idea is that the
                          !gradient on atom I are non-zero only if this atom
                          !belong to a molecule of the pair that are
                          !interacting. We save the gradient on atom I belonging
                          !to mol i coming from the other atoms on mol k and
                          !viceversa. See Overleaf note for a better discussion. 
                          
                          ! calc gradient on atoms belonging to mol i, from
                          ! atoms belonging to mol k, and store the former array
!print *, "state", orb*(k-1)+d, "atoms", (i-1)*atoms_per_mol_i+1, i*atoms_per_mol_i 
                          CALL calc_grad_coup( crd_i, crd_k, charges1, charges2, grad1) 

                          !> store gradient in the correct array entry
                          sh_env%grad_tresp(1,orb*(k-1)+d,(i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(1,:)
                          sh_env%grad_tresp(2,orb*(k-1)+d,(i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(2,:)
                          sh_env%grad_tresp(3,orb*(k-1)+d,(i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(3,:)

                          ! calc gradient on atoms belonging to mol k, from
                          ! atoms belonging to mol i, and store the former arrays
!print *, "state", orb*(i-1)+m, "atoms", (k-1)*atoms_per_mol_k+1, k*atoms_per_mol_k 
                          CALL calc_grad_coup( crd_k, crd_i, charges2, charges1, grad1) 
                          
                          !> store gradient in the correct array entry
                          sh_env%grad_tresp(1,orb*(i-1)+m,(k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(1,:)                        
                          sh_env%grad_tresp(2,orb*(i-1)+m,(k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(2,:)                        
                          sh_env%grad_tresp(3,orb*(i-1)+m,(k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(3,:)                        
 
                          ! I should find a way to store grad1 in a more general
                          ! global_grad array that in principle should have
                          ! dimension (3, mat, mat, natoms_TOTAL_active) or in a reduce
                          ! form assuming all atoms_per_mol the same (3, mat, mat*natoms_per_mol))
  
                          DEALLOCATE(crd_i,crd_k)
                      end do
                   end do
                 else
                    do d=1,orb 
                       do m=1,orb 
                          sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m) = sh_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                       end do 
                    end do
                 end if 
              end if 


           end if
        end do 
     end do

!!TEST
!OPEN (UNIT = 21, FILE = "GRAD_formatted.dat")
!
!do i=1,SIZE(sh_env%grad_tresp,3)
!  do k=1, SIZE(sh_env%grad_tresp,2)
!     write(21,*) "atom", i, "state", k, sh_env%grad_tresp(:,k,i)
!  end do
!end do
!
!close (21)
!!stop
  end subroutine  calc_tresp_coup



! COMMENT_SG this subroutine is just a replication of the one below, you should
! merge them together
subroutine calc_grad_coup( crd1, crd2, charges1, charges2, grad1)

   real(kind=dp), allocatable, dimension(:,:), intent(in)  :: crd1, crd2
   real(kind=dp),  dimension(:,:), allocatable, intent(out) :: grad1 !gradients acting on atoms of mol1
   real(kind=dp), allocatable, dimension(:), intent(in) :: charges1, charges2

   real(kind=dp),  dimension(3) :: vect
   real(kind=dp)              :: dist
   real(kind=dp)              :: x=2.0
   logical                    :: Gaussian_charges
   integer :: t,j
   real(kind=dp), allocatable, dimension(:)               :: char1, char2

   !Ang2Bohr=1.889725989
   !AU2eV=27.211399

   Gaussian_charges = .True. !put in input file

   allocate(grad1(3,SIZE(crd1,2)))
   !allocate(grad2(3,SIZE(crd2,2))) 
   allocate(char1(SIZE(charges1)))
   allocate(char2(SIZE(charges2)))

   !> Gaussian charges must be divided by sqrt(2)
   if (Gaussian_charges) then
     char1(:) = (charges1(:)/SQRT(x))
     char2(:) = (charges2(:)/SQRT(x))
   else
     char1(:) = charges1(:)
     char2(:) = charges2(:)
   end if

   ! These the gradient from molecule 2 atoms onto molecule 1 atoms. Basically
   ! each atom on molecule 1 will feel the gradients given by the sum of the
   ! atoms on molecule 2. See Note Overleaf for the equation. 

   grad1(:,:) = 0.0d0 !NB gradient is put to zero here to sum over the atoms of
                      !the other mol
   vect = 0.0d0
   do t=1,SIZE(crd1,2)
      do j=1,SIZE(crd2,2)
         vect(:) = (crd1(:,t))-(crd2(:,j))
         dist = SQRT(DOT_PRODUCT(vect,vect))

         grad1(1,t) =  grad1(1,t) + ((char1(t)*char2(j))/(dist*dist*dist))*vect(1)
         grad1(2,t) =  grad1(2,t) + ((char1(t)*char2(j))/(dist*dist*dist))*vect(2)
         grad1(3,t) =  grad1(3,t) + ((char1(t)*char2(j))/(dist*dist*dist))*vect(3)
      end do
      grad1(:,t) = -grad1(:,t)
!print *, "grad", grad1(:,t)
   end do
!print *, "grad", grad1(:,:)

   deallocate(char1, char2)
end subroutine  calc_grad_coup


!**********************************************************************************
!> \brief given crds of mol1 and mol2 and the respective charges (that might in
!general be different for the two mols), this calculates the tresp_coup. 
!> author SG 
!**********************************************************************************
subroutine calc_trsp_ab(crd1, crd2, charges1, charges2, tresp_coup)

   real(kind=dp), allocatable, dimension(:,:), intent(in)  :: crd1, crd2
   real(kind=dp), intent(out)                              :: tresp_coup
   real(kind=dp), allocatable, dimension(:), intent(in)    :: charges1, charges2

   real(kind=dp),  dimension(3) :: vect
   real(kind=dp)                :: dist !Ang2Bohr, AU2eV
   real(kind=dp)                :: x=2.0
   logical                      :: Gaussian_charges
   integer :: t,j
   real(kind=dp), allocatable, dimension(:)               :: char1, char2

!   Ang2Bohr=1.889725989
!   AU2eV=27.211399

   Gaussian_charges = .True. !put in input file

   allocate(char1(SIZE(charges1)))
   allocate(char2(SIZE(charges2)))

   !> Gaussian charges must be divided by sqrt(2). Have a look at Multiwfn.pdf doc
   if (Gaussian_charges) then
     char1(:) = (charges1(:)/SQRT(x))
     char2(:) = (charges2(:)/SQRT(x))
   else
     char1(:) = charges1(:)
     char2(:) = charges2(:)
   end if

   vect = 0.0d0
   !> coulombic interaction
   tresp_coup = 0.0d0
   do t=1,SIZE(crd1,2)
      do j=1,SIZE(crd2,2)
         vect(:) = (crd1(:,t))-(crd2(:,j))
         dist = SQRT(DOT_PRODUCT(vect,vect))
         tresp_coup = tresp_coup + (char1(t)*char2(j))/dist
      end do
   end do

!print *, "TRESP COUP:", tresp_coup
   deallocate(char1, char2)
end subroutine calc_trsp_ab


!**********************************************************************************
!> \brief for now it calculates a tridiagonal matrix as hamiltonian with
!diagonal elements taken from FF and just off diagonals from input
!> author SG
!**********************************************************************************
  subroutine calc_frozen_H(sh_env, energies_tp, orb)
     TYPE(sh_env_type)                      :: sh_env
     REAL(kind=dp), DIMENSION(:),  POINTER          :: energies_tp
     integer, intent(in)                            :: orb
        
     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_frozen_H', &
     routineP = moduleN//':'//routineN

     integer :: i, k, m, d, frz_ind

     sh_env%H_tp(:,:) = cmplx(0.0D0, 0.0D0,dp)
     ! these loops are inefficients but do the job
     do i=1,SIZE(energies_tp) !loop over all non-generate states 
        do k=1,SIZE(energies_tp)
           !> diagonal construction
           if (k==i) then 
              !since you may want to include multiple orbital per state (for the
              !degenerate case)), you
              !should loop over orbitals
              do d=1,orb
                 sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0,(-energies_tp(i)),dp) !minus because we propagate -X
              end do
           else 
              if (.not. sh_env%connectivity_frz_c) then
                 ! only valid for first off-diag: TO IMPROVE FURTHER for other
                 ! off diagonals
                 if ((k+1 == i) .or.(k-1 == i)) then
                    do d=1,orb
                       do m=1,orb
                          sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0,-sh_env%first_off_d_Hab,dp) !sign of the off-diag 
                                                                                                     !does not matter for 1D chains
                       end do 
                    end do
                 end if
              end if 
          end if 
        end do 
     end do

     ! add frz_ couplings from input and also diagonal if different from zero
     if (sh_env%connectivity_frz_c) then
         do frz_ind=1,SIZE(sh_env%frz_couplings) ! loop over stuff in the input
            i = sh_env%frz_coup_connectivity(frz_ind,1) 
            k = sh_env%frz_coup_connectivity(frz_ind,2)
            do d=1,orb
               ! we use the sum because we initialize H to zero and in this way
               ! we can add the off-diagonals and diagonals if present
               sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = &
                  sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) + cmplx(0.0D0,(-sh_env%frz_couplings(frz_ind)),dp)
            end do
         end do
     end if 
  end subroutine calc_frozen_H

!**********************************************************************************
!> \brief Calculate adiabatic forces according to Spencer 2016 and Carof 2019
!        for each atom in the system on the state passed as input. This
!        subroutines is general for Frozen Coupling case, AOM and TRESP (becuase
!        the off-diagonal gradients are stored in the same way).  
!> \note This subroutine has not been tested for degeneracy, but it should be
!        general enough to handle it.
!>       A huge part of this subroutine has been taken from AC's calcforces_aom
!> \author SG  
!**********************************************************************************
subroutine calcforces_smart_grad(state, C, Fin, forces_adiab, & 
                  ext_pres,  active_system_info, orb, nacv, convert_nacv, frozen_c)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom
        integer, intent(in)                         :: state, orb
        REAL(KIND=dp), intent(in), dimension(:,:)   :: C
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), dimension(:,:), POINTER      :: forces_adiab
        LOGICAL, dimension(:), INTENT(IN)           :: ext_pres
        INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                      INTENT(IN)     :: active_system_info
        REAL(KIND=dp), intent(in), OPTIONAL          :: convert_nacv
        LOGICAL, OPTIONAL, INTENT(IN)                :: frozen_c
        REAL(KIND=dp), intent(in), dimension(:,:,:), &
                                      OPTIONAL       :: nacv

        CHARACTER(len=*), PARAMETER :: routineN ='calcforces_smart_grad', &
                               routineP = moduleN//':'//routineN

        integer                                     :: i, j,  d, m
        INTEGER                                     :: natom
        INTEGER                                     :: num_mol, num_j
        REAL(KIND=dp)                               :: my_convert_nacv 
        LOGICAL                                     :: my_frozen_c

        natom = SIZE(Fin,2)

        !COMMENT_SG: This subroutine should be checked further for degeneracy.
        !Note that Fin(d,j,i) is the array of diagonal forces and its dimension
                            !stops at the number of total force_evals
                            !considered. 
        !The following subroutine may NOT work if orb !=1 (degeneracy). It needs
        !testing.

        ! Based on paper [Spencer16]: Spencer et al., J. Chem. Phys. 145, 064102 (2016)
        ! Use Eq. (21) of paper [Spencer16]

        my_frozen_c = .FALSE.
        IF (PRESENT(frozen_c)) my_frozen_c = frozen_c


        do d=1,3                        ! sum of cartesian 
           num_j = 0               ! active atom index
           do j=1,natom            ! sum over all atoms
              forces_adiab(d,j) = 0.0D0

              IF (ext_pres(j)) THEN   
                 num_j   = num_j + 1                 ! update active atom index
                 ! Get the molecule to which j belong, for a general
                 ! system. This is a molecular index needed below. Looping
                 ! over molecules is necessary here.
                 do i=1,SIZE(active_system_info, 2) ! loop over non-deg states
                    IF ((j .GE.  active_system_info(1,i)) .AND. (j .LE. active_system_info(2,i))) THEN
                       num_mol = i
                       EXIT 
                    END IF
                 end do
              ENDIF

              do i=1,SIZE(active_system_info, 2)
                 do m=1,orb
                   ! Diagonal term, it depends on the forces
                   ! I assume that diagonal forces
                   ! for two degenrate
                   ! states are the same as for site energies
                   forces_adiab(d,j) = forces_adiab(d,j) + (abs(C(orb*(i-1)+m, orb*(state-1)+m))**2)*Fin(d,j,i)
 
                   ! calculate off diagonal forces if
                   ! couplings are not frozen
                   IF (.NOT. my_frozen_c) THEN

                      ! Of diagonal term, only for active atom j present in
                      ! molecule num_mol
                      IF (ext_pres(j)) THEN
                         IF (i.ne.num_mol) THEN
                             IF (PRESENT(convert_nacv)) THEN
                                !> Forces in case of AOM formalism
                                my_convert_nacv = convert_nacv ! to convert the FOB NACV
                             ELSE 
                                !> Forces in case of general gradients (smartly stored as
                                !> nacv in AOM case)
                                my_convert_nacv = 1.0D0 !gradient are already a negative object so no minus sign here
                             END IF

                             forces_adiab(d,j) = forces_adiab(d,j) + &
                                   ((C(orb*(num_mol-1)+m,orb*(state-1)+m))*&
                                     nacv(d, orb*(i-1)+m, num_j)*C(orb*(i-1)+m,orb*(state-1)+m))*my_convert_nacv*2.0D0
                         ENDIF
                      ENDIF
                   END IF
                 end do
              end do 
           ENDDO
        ENDDO

end subroutine calcforces_smart_grad



!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
subroutine surf_hop_type(adiab_section, C1, C2, deltat, &
           new_state, attempt_hop, decohere_criterion, sh_print_env, sh_env)
!Choose whether doing surface hopping with trivial hopping correction or without
!it following the old code 
        
       TYPE(section_vals_type), POINTER               :: adiab_section
       REAL(KIND=dp), dimension(:,:), intent(in) :: C1, C2
       REAL(KIND=dp), intent(in) :: deltat
       INTEGER, INTENT(INOUT)        :: new_state
       LOGICAL, INTENT(OUT)        :: attempt_hop
       REAL(KIND=dp), INTENT(IN) :: decohere_criterion
       TYPE(sh_print_env_type)     :: sh_print_env
       TYPE(sh_env_type)          :: sh_env
       
       CHARACTER(len=*), PARAMETER :: routineN ='surf_hop_type', &
                               routineP = moduleN//':'//routineN

       CHARACTER(LEN=default_string_length)           :: surfhop_keyword
       TYPE(section_vals_type), POINTER               :: optimization_section
       LOGICAL                                        :: overlap_key, do_trivial_corr

       NULLIFY(optimization_section)

       overlap_key = .True.
   
       optimization_section => section_vals_get_subs_vals(adiab_section,"OPTIMIZATION")
       CALL section_vals_val_get(optimization_section,"SURF_HOP_CHOICE", c_val=surfhop_keyword)
       sh_print_env%additional_prints = .FALSE.

       SELECT CASE(surfhop_keyword)
       CASE("TRIVIAL_CORR_WITH_OVERLAP")
            sh_print_env%additional_prints = .TRUE.
            do_trivial_corr = .TRUE.
            CALL modified_surfacehop(C1, C2, deltat, &
              new_state, attempt_hop, decohere_criterion, &
              sh_print_env, sh_env, adiab_section, do_trivial_corr)
       CASE("TRIVIAL_HOP_CORRECT")
            sh_print_env%additional_prints = .TRUE.
            do_trivial_corr = .TRUE.
            CALL modified_surfacehop(C1, C2, deltat, &
              new_state, attempt_hop, decohere_criterion, &
              sh_print_env, sh_env, adiab_section, do_trivial_corr)
       CASE("UNMODIFIED_SURF_HOP")
            do_trivial_corr = .FALSE.
            !CALL  surfacehop(deltat, old_state, C2, &
            !  new_state, attempt_hop, decohere_criterion, sh_print_env, &
            !                                  sh_env, adiab_section)
            CALL modified_surfacehop(C1, C2, deltat, &
              new_state, attempt_hop, decohere_criterion, &
              sh_print_env, sh_env, adiab_section, do_trivial_corr) 
       CASE DEFAULT
            CPABORT("There is no method to add trivial crossing correction")
       END SELECT
END SUBROUTINE surf_hop_type


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
subroutine modified_surfacehop(C1, C2, deltat, &
              new_state, attempt_hop, decohere_criterion, &
              sh_print_env, sh_env, adiab_section,       & 
              do_trivial_corr) 
        !using the method outlined in John Tully's Electronic Transitions paper
        !calculate the probability of the molecular dynamics switching to evolve on a different potential surface
        REAL(KIND=dp), dimension(:,:), intent(in) :: C1, C2
        REAL(KIND=dp), intent(in) :: deltat
        INTEGER, INTENT(INOUT)        :: new_state
        LOGICAL, INTENT(OUT)        :: attempt_hop
        REAL(KIND=dp), INTENT(IN) :: decohere_criterion
        TYPE(sh_print_env_type) :: sh_print_env
        TYPE(sh_env_type), INTENT(INOUT)             :: sh_env
        TYPE(section_vals_type), POINTER              :: adiab_section
        LOGICAL, INTENT(IN)                           :: do_trivial_corr

        CHARACTER(len=*), PARAMETER :: routineN ='modified_surfacehop', &
        routineP = moduleN//':'//routineN

        CHARACTER(LEN=default_string_length)         :: wfpropag_keyword
        integer :: i, j, mol!, k, l
        REAL(KIND=dp) :: hop, sumg, ex_prob, min, g_exa, err
        REAL(KIND=dp), dimension(:), allocatable :: b, g
        !REAL(KIND=dp), dimension(:), allocatable :: g_before, g_after, g_renorm
        COMPLEX(KIND=dp), dimension(:), allocatable :: ct, cttm_ad
        COMPLEX(KIND=dp), dimension(:), allocatable :: ak, akttm
        INTEGER        :: state_min 
        REAL(KIND=dp), dimension(:), ALLOCATABLE :: enerdiff
 
        ! mol should be called "mat" according to the convension used so far,
        ! because "mat" is compatible with degenerate orbitals
        mol = size(sh_env%H_tp,1)

        ! Based on paper [Spencer16]: Spencer et al., J. Chem. Phys. 145, 064102 (2016)
        ! and paper [Wang14] :

        !many allocations
        allocate(ak(mol))
        allocate(b(mol))
        allocate(g(mol))
        allocate(akttm(mol))        
        allocate(enerdiff(mol))
        allocate(ct(mol))
        allocate(cttm_ad(mol)) 
        ALLOCATE(sh_print_env%proba(mol))

        ! Find the state closest in energy according to [Wang14]
        IF (do_trivial_corr) THEN 
             do i=1,mol
                     if (i .ne. sh_env%active_state) then 
                             enerdiff(i) = abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)) 
                     else if (i .eq. sh_env%active_state) then 
                             enerdiff(sh_env%active_state) = 0.0D0
                     end if
             end do
             if (enerdiff(1) == enerdiff(sh_env%active_state)) then
                     min = enerdiff(2)
                     do i= 3,mol
                            if (enerdiff(i) < min .and. i .ne. sh_env%active_state) then
                                    min = enerdiff(i)
                                     state_min = i
                             end if          
                     end do
             else 
                     min = enerdiff(1)
                     do i= 2,mol
                            if (enerdiff(i) < min .and. i .ne. sh_env%active_state) then
                                    min = enerdiff(i)
                                    state_min = i
                             end if          
                     end do
             end if
             do i =1,mol
                     if (enerdiff(i) == min) then
                            state_min = i
                     end if
             end do
             sh_print_env%enermin = min
        ENDIF

        !if the propagation occurs in the diabatic basis, transfor diabatic coefficients
        !into adiabatic coefficients
        CALL section_vals_val_get(adiab_section,"REPRESENTATION", c_val=wfpropag_keyword)
        SELECT CASE(wfpropag_keyword)
        CASE("DIABATIC_BASIS")
               !CALL diab_to_adiab(sh_env, C2, old_state, mol, cstate, ct, &
               !                                    C1, cttm_ad, cstate_tm)
               !     CALL diab_to_adiab(sh_env, C2, old_state, mol, cstate, ct)
               !end if
               CALL diab_to_adiab(C2, sh_env%ctin, ct) 
               CALL diab_to_adiab(C1, sh_env%ctin_t, cttm_ad) 
        CASE("ADIABATIC_BASIS")
               CPABORT("Not Ready")
               !ct(:) = sh_env%ctt(:)
               !cttm_ad(:) = sh_env%cttm(:)
               !cstate = sh_env%ctt(old_state)
               !cstate_tm = sh_env%cttm(old_state)
        CASE DEFAULT
               CPABORT("There is no method to propagate the wavefunction.")
        END SELECT


        do i=1,mol
               !construct state density for t present and t-deltat from the active state to
               !other state
                ak(i) = ct(i)*conjg(ct(sh_env%active_state))
                akttm(i) = cttm_ad(i)*conjg(cttm_ad(sh_env%old_state))
                !>assign 0.0 to b(i) if the hop is forbidden
                if (sh_env%HS_propagation .and. &
                              (i .eq. sh_env%class_forbidden_hop_t(i))) then
                   b(i) = 0.0D0 
                else 
                   b(i) = -2.0D0*REAL(conjg(ak(i))*sh_env%adiab_nace_tp(i,sh_env%active_state))
                end if
        end do

        !calc exact total surface hopping probability (also termed g_tot) 
        !according to Prezhdo and Wang 2014
        !>in this case we need to use the new enforced state to avoid discontinuity
        !in the time derivative
        ex_prob = (abs(akttm(sh_env%old_state))- abs(ak(sh_env%active_state)))/abs(akttm(sh_env%old_state))
        sh_print_env%exprob = ex_prob 

        ! Use Eq. (4) of [Spencer16]
        ! calculate sumg that can also be negative because each individual probability
        ! (flux) could be negative 
        !calculate g and sum(g)
        sumg = 0.0D0
        do i=1,mol
                !g(i) is the probability to hop from the current state (old_state) to any other
                !state (i) 
                g(i) = (deltat*b(i))/abs(ak(sh_env%active_state))
                if (i .ne. sh_env%active_state) then
                        sumg = sumg + g(i)
                end if
        end do
        sh_print_env%sumflux   = sumg
        err = abs(sumg - ex_prob)
        IF (err.gt.0.1) THEN
           print*, "FOBSH| *** WARNING: difference between exact flux and sum of Tully's proba is", err, "***"
        ENDIF

        !EFFECTIVE TRIVIAL CROSSINGS CORRECTION
        ! Use Eq. () of [Wang14]
        g_exa = 0
        if (do_trivial_corr) THEN
           g_exa = ex_prob - sumg + g(state_min)
           g(state_min) = g_exa
        endif 

       !CORRECTION FOR NEGATIVE NUMBER
       do i =1,mol
                if (g(i) < 0.0D0) then
                        g(i) = 0.0D0
                end if
       end do

       !NORMALIZATION
       !due to finite timestep, sometimes the probabilities may sum to greater than 1
       !if this happens, rescale them down to sum = 1
       sumg = 0.0D0
       do i=1,mol
               if (i .ne. sh_env%active_state) then
                       sumg = sumg + g(i)
               end if
       end do
       if (sumg > 1.0D0) then
               print*, "FOBSH| *** WARNING: We have to normalize the probability, the sum was:", sumg, "***"
               do i=1,mol
                       g(i) = g(i)/sumg
                       g(sh_env%active_state) = 0.0D0
               end do
               sh_print_env%renorm = .TRUE.
               sh_print_env%sumprob   = sumg
       else
               sh_print_env%renorm = .FALSE.
               sh_print_env%sumprob   = sumg
               g(sh_env%active_state) = 1.0D0 - sumg
       end if
       sh_print_env%proba     = g

        !do i=1, mol
        !   print *, "AC_TRIVIAL", g_before(i), g_after(i), g_renorm(i)
        !enddo


        hop = next_random_number(sh_env%rng_stream)
        print *, "RANDOM_NUMBER_HOPPING_PROBA", hop
        sh_print_env%hop       = hop
        sh_print_env%old_state = sh_env%active_state
        !and see if this random number causes a surface hop
        sumg = 0.0D0
        do i=1,mol
                sumg = sumg + g(i)
                if (hop <= sumg) then
                        if (i .ne. sh_env%active_state) then
                                new_state = i
                                !write(24,*) i                        
                                !decohere = .FALSE.
                                sh_env%decohere = .FALSE.
                                attempt_hop = .TRUE.
                                print *, "FOBSH| ATTEMPTED_HOP"
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
                                EXIT
                        elseif (i == sh_env%active_state) then
                                new_state = sh_env%active_state 
                                attempt_hop = .FALSE.
                                sh_print_env%attempt_hop    = attempt_hop
                                sh_print_env%new_state = new_state
!                                write(24,*) state
                                
                                !if there has not been any surfacehop, we need to consider the possibility of the decoherence correction
                                !if the non-adiabatic coupling dphi to all other states is below a certain threshold
                                !then a wavefunction collapse will be enforced
                                sh_env%decohere = .TRUE.
                                do j=1,mol
                                        if (abs(sh_env%adiab_nace_tp(sh_env%active_state,j)) & 
                                                     > decohere_criterion .and. sh_env%active_state .ne. j) then
                                               ! decohere = .FALSE.
                                                sh_env%decohere = .FALSE.
                                                EXIT
                                        end if
                                end do
                               ! decohere = .TRUE.
                                EXIT
                        end if
                end if
        end do

        CALL adiab_print_sh_info(sh_print_env, keyword = 'TULLY_TEST') 

        DEALLOCATE(sh_print_env%proba)
        DEALLOCATE(enerdiff)
        deallocate(ak)
        deallocate(b)
        deallocate(G)
        deallocate(akttm)
end subroutine modified_surfacehop

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE check_hop_attempt(velocities, atomic_kind_set, nparticle_kind, local_particles, &
 force_env, rescale_keyword)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        INTEGER, INTENT(IN)                      :: nparticle_kind
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(force_env_type), POINTER            :: force_env
        CHARACTER(LEN=default_string_length), OPTIONAL  :: rescale_keyword

        CHARACTER(len=*), PARAMETER :: routineN = 'check_hop_attempt', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: delta_energy, rescale_factor
        REAL(kind=dp)                            :: kinetic_energy_prime, kinetic_energy_qsys, A, B
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: hop_nacv, hop_nacv_oldold, hop_nacv_newnew
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: hop_nacv_oldnew, hop_nacv_newold
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE  :: F_11, F_1, F_2, F_3
        CHARACTER(LEN=default_string_length)     :: my_rescale_keyword
        CHARACTER(LEN=default_string_length)     :: nacv_keyword, reversal_keyword
        TYPE(section_vals_type), POINTER         :: mixed_section, adiab_section!, coupling_section
        TYPE(section_vals_type), POINTER         :: force_env_section, e_conservation_section
        INTEGER                                  :: old_state, new_state, orb
        TYPE(sh_print_env_type)                  :: sh_print_env
        LOGICAL                                  :: translate_keyword
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities_com
        REAL(KIND=dp)                            :: convert_nacv


        NULLIFY(mixed_section)
        NULLIFY(adiab_section)
        !NULLIFY(coupling_section)
        NULLIFY(e_conservation_section)

        CALL force_env_get(force_env=force_env,&
                       force_env_section=force_env_section)
        mixed_section => section_vals_get_subs_vals(force_env_section,"MIXED")
        adiab_section => section_vals_get_subs_vals(mixed_section, "ADIABATIC")
        !coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        e_conservation_section => section_vals_get_subs_vals(adiab_section, "ENERGY_CONSERVATION")

!!        IF (PRESENT(rescale_keyword)) THEN
!!           my_rescale_keyword = rescale_keyword
!        IF (.NOT. force_env%mixed_env%sh_env%aom_couplings_keyword) THEN
!           CALL section_vals_val_get(e_conservation_section,"METHOD_RESCALING",c_val=my_rescale_keyword)
!            !WARNING:::::: when AOM Hamiltonian is not used the only rescaling
!            !method that will work is SIMPLE rescaling, SIMPLE_QSYS would work
!            !only by defining ext_pres LOGICAL array and NACV FAST rescaling
!            !would need to be generalized for any Hamiltonian
!            if (my_rescale_keyword == "NACV") then
!                CPABORT("Poss. rescaling with H different from AOM are SIMPLE and SIMPLE_QSYS")
!            end if
!            !my_rescale_keyword = "SIMPLE"
!        ELSE
!           CALL section_vals_val_get(e_conservation_section,"METHOD_RESCALING",c_val=my_rescale_keyword)
!        ENDIF

        CALL section_vals_val_get(e_conservation_section,"METHOD_RESCALING",c_val=my_rescale_keyword)
        CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS",l_val=translate_keyword)
        CALL section_vals_val_get(adiab_section,"NUMBER_ORBITALS",i_val=orb)
 
        sh_print_env = force_env%mixed_env%sh_env%sh_print_env

        ALLOCATE(velocities_com(SIZE(velocities, 1), SIZE(velocities, 2)))
        CALL translate_com(velocities, translate_keyword, velocities_com, &
                           local_particles, nparticle_kind = nparticle_kind, atomic_kind_set = atomic_kind_set)

           

        SELECT CASE(my_rescale_keyword)
        CASE("SIMPLE")
           CALL cal_delta_energy(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                 local_particles, delta_energy, kinetic_energy_prime)
           sh_print_env%old_state        = force_env%mixed_env%sh_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%sh_env%new_state
           IF (delta_energy.lt.0) THEN
              force_env%mixed_env%sh_env%successful_hop = .FALSE.
           ELSE 
               !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
              !rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
              !velocities = rescale_factor*velocities
              force_env%mixed_env%sh_env%successful_hop = .TRUE.
              force_env%mixed_env%sh_env%active_state   = force_env%mixed_env%sh_env%new_state
              CALL simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime, rescale_factor)
           ENDIF

           sh_print_env%adiab_energies_t = force_env%mixed_env%sh_env%E_tp
           sh_print_env%delta_energy     = delta_energy
           sh_print_env%successful_hop   = force_env%mixed_env%sh_env%successful_hop
           sh_print_env%kinetic_energy   = kinetic_energy_prime
           IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
           CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')


        CASE("SIMPLE_QSYS")
           CALL cal_delta_ene_QSYS(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                         local_particles, delta_energy, kinetic_energy_qsys)
           sh_print_env%old_state        = force_env%mixed_env%sh_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%sh_env%new_state
           IF (delta_energy.lt.0) THEN
              force_env%mixed_env%sh_env%successful_hop = .FALSE.
           ELSE 
               !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
              !rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
              !velocities = rescale_factor*velocities
              force_env%mixed_env%sh_env%successful_hop = .TRUE.
              !force_env%mixed_env%sh_env%old_state      = force_env%mixed_env%sh_env%active_state
              force_env%mixed_env%sh_env%active_state   = force_env%mixed_env%sh_env%new_state
              CALL simple_vel_resc_qsys(atomic_kind_set, local_particles,  &
                       force_env, nparticle_kind, velocities, delta_energy, kinetic_energy_qsys, &
                         rescale_factor)

           ENDIF

           sh_print_env%adiab_energies_t = force_env%mixed_env%sh_env%E_tp
           sh_print_env%delta_energy     = delta_energy
           sh_print_env%successful_hop   = force_env%mixed_env%sh_env%successful_hop
           sh_print_env%kinetic_energy   = kinetic_energy_qsys
           IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
           CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')


        CASE("NACV")
           CALL section_vals_val_get(e_conservation_section,"METHOD_ADIABATIC_NACV",c_val=nacv_keyword)
           CALL cal_delta_ene_QSYS(nparticle_kind, velocities_com, force_env, atomic_kind_set, &
                                 local_particles, delta_energy, kinetic_energy_prime)
           sh_print_env%kinetic_energy   = kinetic_energy_prime
           sh_print_env%old_state        = force_env%mixed_env%sh_env%active_state
           sh_print_env%new_state        = force_env%mixed_env%sh_env%new_state
           IF ((delta_energy.lt.0).AND.(nacv_keyword.ne."TEST")) THEN
              force_env%mixed_env%sh_env%successful_hop = .FALSE.

              sh_print_env%adiab_energies_t = force_env%mixed_env%sh_env%E_tp
              sh_print_env%delta_energy     = delta_energy
              sh_print_env%successful_hop   = force_env%mixed_env%sh_env%successful_hop
              sh_print_env%old_state        = force_env%mixed_env%sh_env%active_state
              sh_print_env%new_state        = force_env%mixed_env%sh_env%new_state
              IF (sh_print_env%successful_hop) sh_print_env%rescale_factor = rescale_factor
              CALL adiab_print_sh_info(sh_print_env, keyword = 'SIMPLE_RESCALING')
           ELSE 
              !WARNING:::: the only possible NACV rescaling in case of non AOM
              !hamitlonian is FAST. 
              IF (.NOT. force_env%mixed_env%sh_env%aom_couplings_keyword) THEN
                      IF (.NOT. nacv_keyword == "FAST") then
                         CPABORT("NACV rescaling with non AOM H can be only  FAST")
                      END IF
              END IF
              SELECT CASE(nacv_keyword)
              CASE("TEST")
                 old_state = force_env%mixed_env%sh_env%active_state
                 new_state = force_env%mixed_env%sh_env%new_state
                 !sh_env=force_env%mixed_env%sh_env
                 CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env,force_env%mixed_env%F,&
                                                                   hop_nacv_oldold, old_state, old_state, "OLDOLD")
                 CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env, force_env%mixed_env%F,& 
                                                                  hop_nacv_newnew, new_state, new_state, "NEWNEW")
                 CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env, force_env%mixed_env%F,& 
                                                                    hop_nacv_newold, new_state, old_state, "NEWOLD")
                 CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env, force_env%mixed_env%F,& 
                                                                   hop_nacv_oldnew, old_state, new_state, "OLDNEW")
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv_oldold, &
                                  hop_nacv_newnew, hop_nacv_newold, hop_nacv_oldnew) 
                 print*, "ERRORASYM = ", SUM( (hop_nacv_oldnew + hop_nacv_newold ) **2 )
                 hop_nacv = hop_nacv_oldnew
                 print*, "totalnacv =" , SUM(hop_nacv(1,:)),  SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:))
                 print*, "totalnacvasym =" , SUM(hop_nacv_oldnew(1,:) - hop_nacv_newold(1,:)) / 2,  & 
                                             SUM(hop_nacv_oldnew(2,:) - hop_nacv_newold(2,:)) / 2,  &
                                             SUM(hop_nacv_oldnew(3,:) - hop_nacv_newold(3,:)) / 2
              CASE("CONTRIBUTION")
                 CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env, force_env%mixed_env%F,& 
                                                         hop_nacv, exF_11 = F_11, exF_1 = F_1, exF_2 = F_2, exF_3 =F_3)
                 print*, "totalnacv =" , SUM(hop_nacv(1,:)),  SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:))
                 print*, "totalnacvF1 =" , SUM(F_1(1,:)),  SUM(F_1(2,:)), SUM(F_1(3,:))
                 print*, "totalnacvFF1 =" , SUM(F_1(1,:) + F_11(1,:)), &
                                            SUM(F_1(2,:) + F_11(2,:)), &
                                            SUM(F_1(3,:) + F_11(3,:))
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv, F_11, F_1, F_2, F_3) 
                 CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env%mixed_env%sh_env, &
                                              atomic_kind_set, &
                                            local_particles, F_1 + F_11, delta_energy, A, B) 
print*, "DELTAFF1 = ", delta_energy
IF ((delta_energy).lt.0) THEN
     print*, "LOGFF1 = ", 0
ELSE
     print*, "LOGFF1 = ", 1
ENDIF
     
                 CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env%mixed_env%sh_env, &
                                            atomic_kind_set, &
                                            local_particles, F_1, delta_energy, A, B) 
print*, "DELTAF1 = ", delta_energy
IF ((delta_energy).lt.0) THEN
     print*, "LOGF1 = ", 0
ELSE
     print*, "LOGF1 = ", 1
ENDIF
     
              CASE("TOTAL")
                 CALL calc_hop_nacv_aom(force_env%mixed_env%sh_env, force_env%mixed_env%F, hop_nacv)
                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv) 
              CASE("FAST")
                 !> AOM nacv stored
                 IF (force_env%mixed_env%sh_env%aom_couplings_keyword) THEN
                    !COMMENT_SG: calc_hop_nacv_approx is a generalization of calc_hop_nacv_aom in the fast=.TRUE. case
                    !CALL calc_hop_nacv_aom( force_env%mixed_env%sh_env, force_env%mixed_env%F, hop_nacv, fast=.TRUE.)
                    convert_nacv = -force_env%mixed_env%sh_env%scaling/force_env%mixed_env%sh_env%cbar
                    CALL calc_hop_nacv_approx(force_env%mixed_env%sh_env, diagonal_F=force_env%mixed_env%F, &
                                           C=force_env%mixed_env%sh_env%C_tp,&
                                           active_system_info=force_env%mixed_env%sh_env%active_system_info, &
                                           orb=orb, hop_nacv=hop_nacv, &          
                                           nacv=force_env%mixed_env%sh_env%nacv, &
                                           convert_nacv=convert_nacv)
                 !> TRESP gradient stored
                 ELSE IF ( force_env%mixed_env%sh_env%tresp_couplings_keyword) THEN
                    !CALL calc_hop_nacv_frozen_H( force_env%mixed_env%sh_env, force_env%mixed_env%F, hop_nacv)
                    CALL calc_hop_nacv_approx(force_env%mixed_env%sh_env, diagonal_F=force_env%mixed_env%F, &
                                           C=force_env%mixed_env%sh_env%C_tp, &
                                           active_system_info=force_env%mixed_env%sh_env%active_system_info, &
                                           orb=orb, hop_nacv=hop_nacv, &          
                                           nacv=force_env%mixed_env%sh_env%grad_tresp)
                 !> FRZ couplings case
                 ELSE IF ( force_env%mixed_env%sh_env%frozen_couplings_keyword) THEN
                    !CALL calc_hop_nacv_frozen_H( force_env%mixed_env%sh_env, force_env%mixed_env%F, hop_nacv)
                    CALL calc_hop_nacv_approx(force_env%mixed_env%sh_env, diagonal_F=force_env%mixed_env%F, &
                                           C=force_env%mixed_env%sh_env%C_tp,&
                                           active_system_info=force_env%mixed_env%sh_env%active_system_info, &
                                           orb=orb, hop_nacv=hop_nacv, frozen_c=.TRUE.)       
                 ELSE 
                     CPABORT("There is no method to calculate the NACV in FAST approx.")
                 END IF

                 CALL adiab_print_hop_nacv(sh_print_env, hop_nacv) 
              CASE DEFAULT
                CPABORT("There is no method to calculate the NACV.")
              END SELECT 
              CALL cal_delta_energy_nacv(nparticle_kind, velocities_com, force_env%mixed_env%sh_env, &
                                          atomic_kind_set, &
                                          local_particles, hop_nacv, delta_energy, A, B) 
!IF ((delta_energy).lt.0) THEN
!     print*, "LOGNACV = ", 0
!ELSE
!     print*, "LOGNACV = ", 1
!ENDIF
              If (delta_energy.lt.0) THEN
                  force_env%mixed_env%sh_env%successful_hop = .FALSE.
                  CALL section_vals_val_get(e_conservation_section,"METHOD_REVERSAL",c_val=reversal_keyword)
                  CALL hop_nacv_vel_reversing(nparticle_kind, velocities,  & 
                                              force_env%mixed_env%sh_env%forces_adiab, &
                                              force_env%mixed_env%sh_env%forces_adiab_new, &
                                              hop_nacv, atomic_kind_set, local_particles, reversal_keyword, &
                                              sh_print_env,A,B)
              ELSE 
                  force_env%mixed_env%sh_env%successful_hop = .TRUE.
                  force_env%mixed_env%sh_env%active_state   = force_env%mixed_env%sh_env%new_state
                  CALL hop_nacv_vel_rescaling(nparticle_kind, velocities, delta_energy, A, B, hop_nacv, &
                                              atomic_kind_set, local_particles, rescale_factor)
              ENDIF

              sh_print_env%A = A
              sh_print_env%B = B
              sh_print_env%adiab_energies_t = force_env%mixed_env%sh_env%E_tp
              sh_print_env%delta_energy     = delta_energy
              sh_print_env%successful_hop   = force_env%mixed_env%sh_env%successful_hop
              IF (sh_print_env%successful_hop) THEN
                 sh_print_env%rescale_factor = rescale_factor
                 sh_print_env%momentum_shift(:) = rescale_factor *        & 
                          (/ SUM(hop_nacv(1,:)), SUM(hop_nacv(2,:)), SUM(hop_nacv(3,:)) /)
              ENDIF
              CALL adiab_print_sh_info(sh_print_env, keyword = 'NACV_RESCALING')
              !IF (.NOT.(sh_print_env%successful_hop)) THEN
              !   CALL adiab_print_sh_info(sh_print_env, keyword = 'NACV_REVERSING')
              !ENDIF 
           ENDIF 
        CASE DEFAULT
                CPABORT("There is no method to do the rescaling of the velocities.")
        END SELECT

        DEALLOCATE(velocities_com)
   END SUBROUTINE check_hop_attempt




!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE hop_nacv_vel_reversing(nparticle_kind, velocities, forces_old, forces_new, hop_nacv, atomic_kind_set, &
                                     local_particles, reversal_keyword, sh_print_env, A, B)
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: forces_old, forces_new
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: atomic_kind_set
        CHARACTER(LEN=default_string_length)     :: reversal_keyword
        TYPE(sh_print_env_type)                  :: sh_print_env
        REAL(kind=dp), INTENT(IN)                :: A, B

        CHARACTER(len=*), PARAMETER :: routineN = 'hop_nacv_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: product_pd, product_dfo, product_dfn
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        LOGICAL                                  :: reverse
        REAL(KIND=dp)                            :: reversal_factor

        NULLIFY(atomic_kind)
!> TO EXPLAIN  
        IF ((reversal_keyword.ne."ALWAYS").OR.(reversal_keyword.ne."NEVER")) THEN
           product_pd = 0.0D0
           product_dfo = 0.0D0
           product_dfn = 0.0D0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 product_pd = product_pd + mass * (                                   &
                                velocities(1, iparticle) * hop_nacv(1, iparticle) +   &            
                                velocities(2, iparticle) * hop_nacv(2, iparticle) +   &            
                                velocities(3, iparticle) * hop_nacv(3, iparticle) )               
                 product_dfo = product_dfo +  ( &
                                hop_nacv(1, iparticle) * forces_old(1, iparticle) +    &           
                                hop_nacv(2, iparticle) * forces_old(2, iparticle) +    &           
                                hop_nacv(3, iparticle) * forces_old(3, iparticle) )               
                 product_dfn = product_dfn +  ( &
                                hop_nacv(1, iparticle) * forces_new(1, iparticle) +    &           
                                hop_nacv(2, iparticle) * forces_new(2, iparticle) +    &           
                                hop_nacv(3, iparticle) * forces_new(3, iparticle) )               
              ENDDO
           END DO
        ENDIF

        reverse = .FALSE.
        SELECT CASE(reversal_keyword)
        CASE("NEVER")
           reverse = .FALSE.
        CASE("ALWAYS")
           reverse = .TRUE.
        CASE("TRUHLAR")
           IF ((product_pd*product_dfn).lt.0) reverse = .TRUE.
        CASE("SUBOTNIK")
           IF (((product_pd*product_dfn).lt.0) .AND. ((product_dfn*product_dfo).lt.0)) reverse = .TRUE.
        CASE DEFAULT
           CPABORT("There is not method to reverse the velocities after a frustrated hop") 
        END SELECT

        IF (reverse) THEN
           reversal_factor = -B/A
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                      velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                    velocities(:, iparticle ) = velocities(:, iparticle) + (reversal_factor/mass) * hop_nacv(:, iparticle)
                 ENDIF
              ENDDO
           END DO
        ENDIF
        
        sh_print_env%reversal_keyword = reversal_keyword
        sh_print_env%reverse = reverse
        sh_print_env%product_pd = product_pd
        sh_print_env%product_dfn = product_dfn
        sh_print_env%product_dfo = product_dfo

        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE hop_nacv_vel_reversing


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE hop_nacv_vel_rescaling(nparticle_kind, velocities, delta_energy, A, B, hop_nacv, atomic_kind_set, &
                                     local_particles, rescale_factor)
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), INTENT(IN)                :: delta_energy, A, B
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)                :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'hop_nacv_vel_rescaling', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local 
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind

        NULLIFY(atomic_kind)
!> From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular dynamics using semiempirical methods. 
!>     Chem Phys 2008, 349:334–347.
!>  B = -B, compare to the paper

        IF (A == 0) THEN
           rescale_factor = 0
        ELSE
           IF (B < 0) THEN
              rescale_factor = (- B - SQRT(delta_energy))/(2 * A)
           ELSE
              rescale_factor = (- B + SQRT(delta_energy))/(2 * A)
           ENDIF
        ENDIF
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
              IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                   velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                 velocities(:, iparticle ) = velocities(:, iparticle) + (rescale_factor/mass) * hop_nacv(:, iparticle)
              ENDIF
           ENDDO
        END DO
        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE hop_nacv_vel_rescaling

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE simple_vel_resc_qsys(atomic_kind_set,local_particles, force_env, nparticle_kind, & 
                                   velocities, delta_energy, kinetic_energy_qsys, rescale_factor)
        
        INTEGER, INTENT(IN)                      :: nparticle_kind 
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_qsys
        REAL(kind=dp), INTENT(OUT)               :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_resc_qys', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: iparticle_local, iparticle
        INTEGER                                  :: nparticle_local

!> Factor to rescale the velocities (Eq. 49 of PAPER_1)
        rescale_factor = SQRT(  delta_energy /kinetic_energy_qsys)

           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 if (force_env%mixed_env%sh_env%ext_pres(iparticle)) then 
                         velocities(:, iparticle) =rescale_factor*velocities(:, iparticle)
                 end if
              END DO
           END DO
   END SUBROUTINE simple_vel_resc_qsys

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE cal_delta_ene_QSYS(nparticle_kind, velocities, force_env, atomic_kind_set, &
                               local_particles, delta_energy, kinetic_energy_qsys) 
        INTEGER, INTENT(IN)                      :: nparticle_kind
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(force_env_type), POINTER            :: force_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)               :: delta_energy
        REAL(kind=dp), INTENT(OUT)               :: kinetic_energy_qsys

        CHARACTER(len=*), PARAMETER :: routineN = 'cal_delta_ene_QSYS', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind

        NULLIFY(atomic_kind)
           kinetic_energy_qsys = 0.0D0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 if (force_env%mixed_env%sh_env%ext_pres(iparticle)) then 
                      kinetic_energy_qsys = kinetic_energy_qsys +   &
                         mass*( &
                         velocities(1, iparticle)**2 + &
                         velocities(2, iparticle)**2 + &
                         velocities(3, iparticle)**2  )
                 end if
              END DO
           END DO
          kinetic_energy_qsys = 0.5*kinetic_energy_qsys

!print *, "kinetic_energy_qsys", kinetic_energy_qsys
           new_state = force_env%mixed_env%sh_env%new_state
           old_state = force_env%mixed_env%sh_env%active_state
           delta_energy = kinetic_energy_qsys-( &
                     force_env%mixed_env%sh_env%E_tp(new_state) - &
                     force_env%mixed_env%sh_env%E_tp(old_state))
        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE cal_delta_ene_QSYS

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
SUBROUTINE cal_delta_energy(nparticle_kind, velocities, force_env, atomic_kind_set, &
                            local_particles, delta_energy, kinetic_energy_prime) 
     INTEGER, INTENT(IN)                      :: nparticle_kind
     REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
     TYPE(force_env_type), POINTER            :: force_env
     TYPE(distribution_1d_type), POINTER      :: local_particles
     TYPE(atomic_kind_type), DIMENSION(:), &
          POINTER                                :: atomic_kind_set
     REAL(kind=dp), INTENT(OUT)               :: delta_energy
     REAL(kind=dp), INTENT(OUT)               :: kinetic_energy_prime

     CHARACTER(len=*), PARAMETER :: routineN = 'cal_delta_energy', &
     routineP = moduleN//':'//routineN

     REAL(kind=dp)                            :: mass
     INTEGER                                  :: iparticle_kind, new_state, old_state
     INTEGER                                  :: nparticle_local
     INTEGER                                  :: iparticle_local, iparticle
     TYPE(atomic_kind_type), POINTER          :: atomic_kind

     NULLIFY(atomic_kind)
        kinetic_energy_prime = 0.0D0
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
              kinetic_energy_prime = kinetic_energy_prime +   &
                 mass*( &
                 velocities(1, iparticle)**2 + &
                 velocities(2, iparticle)**2 + &
                 velocities(3, iparticle)**2  )
           END DO
        END DO
        kinetic_energy_prime = 0.5*kinetic_energy_prime

        new_state = force_env%mixed_env%sh_env%new_state
        old_state = force_env%mixed_env%sh_env%active_state
        !print *, "OLD_STATE_ISOTROP", old_state
        delta_energy = kinetic_energy_prime-( &
                  force_env%mixed_env%sh_env%E_tp(new_state) - &
                  force_env%mixed_env%sh_env%E_tp(old_state))
     IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
END SUBROUTINE cal_delta_energy

!**********************************************************************************
!> \brief Get first adiabatic state when initialization is DIABATIC
!**********************************************************************************
   SUBROUTINE choose_firstadiabat(C, init, sh_print_env, first_adiabat, sh_env)

        REAL(KIND=dp), dimension(:,:), intent(in) :: C
        INTEGER, INTENT(IN)                       :: init
        TYPE(sh_print_env_type), INTENT(INOUT)    :: sh_print_env
        INTEGER, INTENT(OUT)                      :: first_adiabat
        TYPE(sh_env_type)                        :: sh_env

        CHARACTER(len=*), PARAMETER :: routineN = 'choose_firstadiabat', &
        routineP = moduleN//':'//routineN
        
        INTEGER                                   :: i
        REAL( KIND =dp)                           :: random, sumc

        ALLOCATE(sh_print_env%pop_init(SIZE(C(init,:))))

        random = next_random_number(sh_env%rng_stream)
        sh_print_env%initrandom = random
        sumc = 0.0D0
        DO i=1, SIZE(C(init,:))
           sh_print_env%pop_init(i) = C(init,i)*C(init,i)
        ENDDO
        DO i=1, SIZE(C(init,:))
           sumc = sumc + C(init,i)*C(init,i)
           IF (random.le.sumc) THEN
              first_adiabat = i
              sh_print_env%first_adiabat = i
              EXIT
           ENDIF
        ENDDO
        CALL adiab_print_sh_info(sh_print_env, keyword = 'START')
        DEALLOCATE(sh_print_env%pop_init)

   END SUBROUTINE choose_firstadiabat

!**********************************************************************************
!> \brief Diabatic to adiabatic transformation
!**********************************************************************************
subroutine print_matrix_standard_out(matrix)
     
      INTEGER                                    :: mol,i,j
      REAL(KIND=dp),  DIMENSION(:,:), INTENT(IN) :: matrix

      mol = SIZE(matrix, 1)

PRINT *, "----------------------------------------------------------"
DO i=1, mol
   DO j=1, mol
      WRITE(*, "(A2, e12.4, A2)", advance='no') "  ", matrix(i,j), "  "
   END DO
   PRINT *, " "
END DO
PRINT *, "----------------------------------------------------------"

end subroutine print_matrix_standard_out

!**********************************************************************************
!> \brief Diabatic to adiabatic transformation
!**********************************************************************************
SUBROUTINE diab_to_adiab( C, diab_c, adiab_c, pop)
        REAL(KIND=dp), dimension(:,:), intent(in) :: C
        COMPLEX(KIND=dp), dimension(:), intent(in) :: diab_c
        COMPLEX(KIND=dp), dimension(:), intent(out) :: adiab_c
        REAL(KIND=dp), dimension(:), intent(out), optional    :: pop        

        CHARACTER(len=*), PARAMETER :: routineN ='diab_to_adiab', &
                               routineP = moduleN//':'//routineN

        integer :: i, j!, mol!, k, l

        do i=1, SIZE(diab_c)
           adiab_c(i) = cmplx(0.0, 0.0, dp)
           do j=1, SIZE(diab_c)
              adiab_c(i) = adiab_c(i) + cmplx(C(j,i), 0.0D0,dp)*diab_c(j)
           enddo
           if(PRESENT(pop)) then
              pop(i) = abs(adiab_c(i)*conjg(adiab_c(i)))
           end if
        enddo
END SUBROUTINE

!**********************************************************************************
!> \brief This does energy conservation test according to Fabiano 2008 to check
!> whether the energy conservation criteria for a hop to be successful is met.
!**********************************************************************************
   SUBROUTINE cal_delta_energy_nacv(nparticle_kind, velocities, sh_env, atomic_kind_set, &
                                    local_particles, hop_nacv, delta_energy, A, B, i, j) 
        INTEGER                                 :: nparticle_kind
        INTEGER, OPTIONAL, INTENT(IN)            :: i, j
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(sh_env_type)                       :: sh_env
        TYPE(distribution_1d_type), POINTER      :: local_particles
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                                :: atomic_kind_set
        REAL(kind=dp), DIMENSION(:,:), INTENT(IN)   :: hop_nacv
        REAL(kind=dp), INTENT(OUT)                  :: delta_energy
        REAL(kind=dp), INTENT(OUT)                  :: A, B

        CHARACTER(len=*), PARAMETER :: routineN = 'cal_delta_energy_nacv', &
        routineP = moduleN//':'//routineN

        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind, new_state, old_state
        INTEGER                                  :: nparticle_local, coord
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind


        NULLIFY(atomic_kind)
!From Fabiano E, Keal TW, Thiel W. Implementation of surface hopping molecular dynamics using semiempirical methods. 
!     Chem Phys 2008, 349:334–347.
           A = 0
           B = 0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
!stop
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
!stop
                 DO coord = 1, 3
                    IF ((velocities(1, iparticle)**2 + velocities(2, iparticle)**2 + &
                       velocities(3, iparticle)**2 ).NE.0.0d0) THEN
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ELSE IF ( SUM(velocities).EQ.0) THEN
                       A = A + hop_nacv(coord, iparticle)**2 / mass
                    ENDIF
                    B = B + velocities(coord, iparticle) * hop_nacv(coord, iparticle)
                 ENDDO
              END DO
           END DO
           A = 0.5*A

           IF (PRESENT(j)) THEN
              new_state = j
           ELSE    
              new_state = sh_env%new_state
           ENDIF
           IF (PRESENT(i)) THEN
              old_state = i
           ELSE
              old_state = sh_env%active_state
!print *, "DELTAE_OLDSTATE", old_state
           ENDIF
           delta_energy = B**2 - 4 * A * (  &
                     sh_env%E_tp(new_state) - &
                     sh_env%E_tp(old_state))
        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE cal_delta_energy_nacv

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE translate_com(array, keyword_translate, translate_array, &
                            local_particles, atomic_kinds, nparticle_kind, atomic_kind_set)
        REAL(kind=dp), DIMENSION(:,:),POINTER                     :: array
        LOGICAL                                                   :: keyword_translate
        REAL(KIND=dp), DIMENSION(:,:), POINTER                    :: translate_array
        TYPE(distribution_1d_type), POINTER                       :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL            :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                             :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                                    :: atomic_kind_set

        CHARACTER(len=*), PARAMETER :: routineN = 'translate_com', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), TARGET, ALLOCATABLE, DIMENSION(:) :: my_atomic_kind_set
        INTEGER                                                   :: my_nparticle_kind
        REAL(kind=dp)                                             :: mass
        INTEGER                                                   :: iparticle_kind
        INTEGER                                                   :: nparticle_local
        INTEGER                                                   :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER                           :: atomic_kind
        REAL(kind=dp), DIMENSION(3)                               :: sum_array
        REAL(KIND=dp)                                             :: sum_mass

        NULLIFY(atomic_kind)
        IF (keyword_translate) THEN
           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
              my_atomic_kind_set = atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
                 my_atomic_kind_set = atomic_kind_set
              ENDIF
           ENDIF
           sum_array = 0.0D0
           sum_mass  = 0.0D0
           DO iparticle_kind=1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 sum_array(:) = sum_array(:) +   &
                    mass*array(:, iparticle)
                 sum_mass = sum_mass + mass
              END DO
           END DO
           sum_array = sum_array / sum_mass

           DO iparticle_kind=1, my_nparticle_kind
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 translate_array(:, iparticle) = array(:, iparticle) - sum_array(:)
              END DO   
           END DO   
        ELSE
           translate_array = array
        ENDIF 
      
        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE translate_com

! COMMENT_SG: I think this subroutines is just to print final momentum for
! checking purposes
!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
   SUBROUTINE total_momentum(velocities, atomic_kind_set, nparticle_kind, &
                               local_particles)!, force_env) 
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER                             :: atomic_kind_set
        INTEGER, INTENT(IN)                      :: nparticle_kind
        TYPE(distribution_1d_type), POINTER      :: local_particles
        !TYPE(force_env_type), POINTER            :: force_env
 
        CHARACTER(len=*), PARAMETER :: routineN = 'total_momentum', &
        routineP = moduleN//':'//routineN
 
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        REAL(kind=dp), DIMENSION(3)              :: momentum

        NULLIFY(atomic_kind)
           momentum = 0.0D0 
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                 iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                 momentum(:) = momentum(:) +   &
                    mass*velocities(:, iparticle)
              END DO
           END DO
           !print*, "MOMENTUM = ", MOMENTUM
        IF (ASSOCIATED(atomic_kind)) NULLIFY(atomic_kind)
   END SUBROUTINE total_momentum

!**********************************************************************************
!> \brief Rescale the velocity isotropically for the full system (this means
!>        this rescaling is size estensive it depends on the num. of atoms
!**********************************************************************************
   SUBROUTINE simple_vel_rescaling(velocities, delta_energy, kinetic_energy_prime, rescale_factor)
        REAL(kind=dp), DIMENSION(:,:),POINTER    :: velocities
        REAL(kind=dp), INTENT(IN)                :: delta_energy
        REAL(kind=dp), INTENT(IN)                :: kinetic_energy_prime
        REAL(kind=dp), INTENT(OUT)               :: rescale_factor

        CHARACTER(len=*), PARAMETER :: routineN = 'simple_vel_rescaling', &
        routineP = moduleN//':'//routineN

        !> Factor to rescale the velocities (Eq. 49 of PAPER_1)
        rescale_factor = SQRT(  delta_energy /kinetic_energy_prime )
        velocities = rescale_factor*velocities
   END SUBROUTINE simple_vel_rescaling


!**********************************************************************************
!> \brief Define the active region when the spurious charge transfer correction
!>        is used
!**********************************************************************************
SUBROUTINE define_active_region(sh_env)
     
     TYPE(sh_env_type)                       :: sh_env   
   
     CHARACTER(len=*), PARAMETER :: routineN = 'define_active_region', &
        routineP = moduleN//':'//routineN
    
     REAL(KIND=dp), dimension(:), allocatable :: pop
     REAL(KIND=dp)                            :: sum
     INTEGER                                  :: i 
     INTEGER                                  :: mat

     
     mat = size(sh_env%H_tp,1) 

     ALLOCATE(pop(mat))

     ! do i=1, mol
     pop(:)  = real(sh_env%ctin_before_deco(:)*conjg(sh_env%ctin_before_deco(:)))
     ! end do
     sh_env%active_region(:) = 0
     sum = 0.0
     do i =1, mat
         sum = sum + maxval(pop)
         if (sum .le. 0.995) then
             pop(maxloc(pop)) = -1
         else
             exit
         end if
     end do

     do i =1,mat
        if (pop(i) == -1) then
            sh_env%active_region(i) = 1
        end if
     end do

    !do i =1, mat
    !    write (*,'(I12,A, F12.10)') sh_env%active_region(i),  "  pop = ", pop(i)
    !end do

     DEALLOCATE(pop)
END SUBROUTINE define_active_region

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
SUBROUTINE calc_psi_after_correction(sh_env)
 
     TYPE(sh_env_type)                      :: sh_env
     !COMPLEX(KIND=dp), dimension(:), INTENT(IN) :: ctin, ctin_after_deco
     !COMPLEX(KIND=dp), dimension(:), allocatable, INTENT(out) :: ctin_after_tranf_correction
     !INTEGER, dimension(:), INTENT(in)          :: active_region
 
     CHARACTER(len=*), PARAMETER :: routineN = 'calc_psi_after_correction', &
        routineP = moduleN//':'//routineN
   
     REAL(KIND=dp)             :: sum_ctin, sum_ctin_after_deco, lambda_scaling 
     INTEGER                  :: i
     INTEGER                        :: mat

     mat = size(sh_env%H_tp,1) 
   
     !allocate(ctin_after_tranf_correction(mat))
     sum_ctin = 0.0D0
     sum_ctin_after_deco = 0.0D0
     do i =1, mat
         if (sh_env%active_region(i) == 0) then 
             sum_ctin = sum_ctin + abs(sh_env%ctin_before_deco(i)*&
                                        conjg(sh_env%ctin_before_deco(i)))
         else
             sum_ctin_after_deco = sum_ctin_after_deco + &
                      abs(sh_env%ctin_after_deco(i)*conjg(sh_env%ctin_after_deco(i)))
         end if
     end do

    if (sum_ctin_after_deco == 0.0) then
        lambda_scaling = 1.0D0
    else 
        lambda_scaling = sqrt((1-sum_ctin)/sum_ctin_after_deco)
    end if
    IF (ABS(lambda_scaling-1).GT.1E-03) THEN
       PRINT*, "FOBSH| *** WARNING: Projection to avoid spurious transfer get lambda = ", lambda_scaling, "***"
    ENDIF
    sh_env%ctin_after_tranf_correction(:) = sh_env%ctin_before_deco(:)
    do i =1, mat
       if (sh_env%active_region(i) == 1) then
           sh_env%ctin_after_tranf_correction(i) = sh_env%ctin_after_deco(i)*lambda_scaling
       end if
    end do

END SUBROUTINE calc_psi_after_correction

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
SUBROUTINE classically_forbidden_hops(active_state, diagonal_F, sh_env, mat,  &
                  velocities_com_tp, local_particles, atomic_kinds, nparticle_kind, &
                          atomic_kind_set, class_forbidden_hop)

     INTEGER, INTENT(IN)                                :: active_state, mat
     TYPE(sh_env_type), INTENT(IN)                     :: sh_env
     REAL(kind=dp), DIMENSION(:,:), POINTER             :: velocities_com_tp
     REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN), POINTER   :: diagonal_F
     TYPE(distribution_1d_type), POINTER                :: local_particles
     TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
     INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
     TYPE(atomic_kind_type), DIMENSION(:), &
                   POINTER, OPTIONAL                        :: atomic_kind_set
     INTEGER, dimension(:), allocatable, INTENT(OUT)          :: class_forbidden_hop

     INTEGER :: j
     REAL(kind=dp), DIMENSION(:,:), allocatable     :: hop_nacv
     REAL(kind=dp)                               :: delta_energy
     REAL(kind=dp)                               :: A, B
     REAL(kind=dp)                               :: K_tot, E_cutoff  
     REAL(KIND=dp), dimension(:), allocatable    :: energy_conservation
     
     ALLOCATE(energy_conservation(mat)) 
     ALLOCATE(class_forbidden_hop(mat))

     CALL calc_kinetic_tot(sh_env, velocities_com_tp, local_particles, atomic_kinds, &
                  nparticle_kind, atomic_kind_set, K_tot)
     E_cutoff = 0.0D0
     E_cutoff = sh_env%E_tp(active_state) + K_tot  !THIS CUTOFF IS TO LARGE ALL THE
                                            !STATE ARE BELOW IT
     do j=1, mat
           if ((j .gt. active_state) .and. (sh_env%E_tp(j) .le. E_cutoff)) then
               
              CALL calc_hop_nacv_aom(sh_env, diagonal_F, hop_nacv, active_state, j, fast=.TRUE.)
              CALL cal_delta_energy_nacv(atomic_kinds%n_els, velocities_com_tp, sh_env, atomic_kinds%els, &
                                local_particles, hop_nacv, delta_energy, A, B, active_state, j)
              energy_conservation(j) = delta_energy
              if (energy_conservation(j) .lt. 0.0D0) then 
                   class_forbidden_hop(j) = j
              else 
                   class_forbidden_hop(j) = -1
              end if
           else if ((j .gt. active_state) .and. (sh_env%E_tp(j) .gt. E_cutoff)) then
              class_forbidden_hop(j) = j
              energy_conservation(j) = -2  
           else if (j .le. active_state) then
              energy_conservation(j) = 0.0D0
              class_forbidden_hop(j) = -1
           end if  
     end do
     deallocate(energy_conservation) 
END SUBROUTINE classically_forbidden_hops


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE apply_decoherence(itimes, velocities, force_env, local_particles, nparticle_kind, atomic_kind_set)

       REAL(kind=dp), DIMENSION(:,:),POINTER     :: velocities
       TYPE(atomic_kind_type), DIMENSION(:), &
          POINTER                                :: atomic_kind_set
       INTEGER, INTENT(IN)                       :: nparticle_kind
       TYPE(distribution_1d_type), POINTER       :: local_particles
       TYPE(force_env_type), POINTER             :: force_env
       INTEGER, POINTER                          :: itimes 
       
       CHARACTER(len=*), PARAMETER :: routineN = 'apply_decoherence', &
            routineP = moduleN//':'//routineN 

       REAL(kind=dp), DIMENSION(:,:),POINTER     :: velocities_com_tp
       LOGICAL                                   :: translate_keyword, correct_transfer_keyword 
       TYPE(section_vals_type), POINTER          :: mixed_section, adiab_section!, coupling_section
       TYPE(section_vals_type), POINTER          :: force_env_section, decoherence_section
       TYPE(sh_print_env_type)                   :: sh_print_env

       NULLIFY(mixed_section)      
       NULLIFY(adiab_section)      
       !NULLIFY(coupling_section)      
       NULLIFY(decoherence_section)      

       CALL force_env_get(force_env=force_env,&
                       force_env_section=force_env_section)
        mixed_section => section_vals_get_subs_vals(force_env_section,"MIXED")
        adiab_section => section_vals_get_subs_vals(mixed_section, "ADIABATIC")
        !coupling_section => section_vals_get_subs_vals(adiab_section, "AOM")
        decoherence_section => section_vals_get_subs_vals(adiab_section,"DECOHERENCE") 

        CALL section_vals_val_get(adiab_section,"CENTER_OF_MASS",l_val=translate_keyword)
        sh_print_env = force_env%mixed_env%sh_env%sh_print_env
 
        ALLOCATE(velocities_com_tp(SIZE(velocities, 1), SIZE(velocities, 2)))
        CALL translate_com(velocities, translate_keyword, velocities_com_tp, &
              local_particles, nparticle_kind = nparticle_kind, atomic_kind_set = atomic_kind_set)
   
        force_env%mixed_env%sh_env%ctin_before_deco = force_env%mixed_env%sh_env%ctin
        CALL decohere_wf(adiab_section, force_env%mixed_env%sh_env, force_env%mixed_env%F, sh_print_env, itimes, &
                         velocities_com_tp, local_particles, nparticle_kind,atomic_kind_set)
        CALL section_vals_val_get(decoherence_section,"SPURIOUS_TRANSFER_CORR", l_val=correct_transfer_keyword)
        
        if (correct_transfer_keyword) then
              !> this subroutine define the active regione for local decoherence, it
              !takes sh_env%ctin_before_deco  as argument
              CALL define_active_region(force_env%mixed_env%sh_env)
              !> the new sh_env%ctin after deco is given to sh_env%ctin_after_deco
              force_env%mixed_env%sh_env%ctin_after_deco=force_env%mixed_env%sh_env%ctin
              CALL calc_psi_after_correction(force_env%mixed_env%sh_env)

              !> the value of sh_env%ctin is updated again using the spurious charge
              !transfer correction
              force_env%mixed_env%sh_env%ctin = force_env%mixed_env%sh_env%ctin_after_tranf_correction                        
        end if

        DEALLOCATE(velocities_com_tp)

  END SUBROUTINE

 
!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE decohere_wf(adiab_section,sh_env, diagonal_F, sh_print_env, itimes, &
                         velocities_com_tp, local_particles, nparticle_kind, &
                         atomic_kind_set, atomic_kinds)

        TYPE(section_vals_type), POINTER               :: adiab_section
        TYPE(sh_print_env_type)                        :: sh_print_env
        TYPE(sh_env_type)                             :: sh_env
 !       TYPE(sh_env_type)                              :: sh_env
        INTEGER, POINTER                               :: itimes
        REAL(kind=dp), DIMENSION(:,:), POINTER             :: velocities_com_tp
        REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN), POINTER             :: diagonal_F
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
                  POINTER, OPTIONAL                        :: atomic_kind_set

!        LOGICAL, dimension(:), POINTER, OPTIONAL         :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN ='decohere_wf', &         
                                       routineP = moduleN//':'//routineN

        TYPE(section_vals_type), POINTER               :: decoherence_section
!       logical                                            :: edc_deco_time
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: tau!, tau_forces
        CHARACTER(LEN=default_string_length)           :: wfprop_keyword
        CHARACTER(LEN=default_string_length)           :: deco_keyword
        REAL(kind=dp)                                 :: K_tot
        INTEGER                                        :: mol
        REAL(KIND=dp)                                 :: deltatMD
!        REAL(KIND=dp)                                  :: C_param, E_0_param
  

        NULLIFY(decoherence_section)

        CALL section_vals_val_get(adiab_section,"REPRESENTATION", &
                                                     c_val=wfprop_keyword)
        decoherence_section => section_vals_get_subs_vals(adiab_section,"DECOHERENCE")
        CALL section_vals_val_get(decoherence_section,"DECOHERENCE_CORRECTIONS",c_val=deco_keyword)
        
       
        mol = SIZE(sh_env%H_tp, 1) 
        deltatMD = sh_env%deltat 
        !sh_print_env = sh_env%sh_print_env    

        SELECT CASE(wfprop_keyword)
        CASE("DIABATIC_BASIS")   
            SELECT CASE(deco_keyword)                                                          
            CASE("INSTANT_COLLAPSE")
              CALL wf_collapse(sh_env, sh_print_env, itimes)
            CASE("TRESH_ONLY_COLLAPSE")  
              CALL  wf_tresh_collaps(sh_env, sh_print_env, itimes)
            CASE("DAMPING")
              CALL calc_kinetic_tot(sh_env, velocities_com_tp, local_particles, atomic_kinds,&
                                    nparticle_kind, atomic_kind_set, K_tot, sh_env%ext_pres)
              CALL eval_deco_tau(sh_env,diagonal_F, tau, adiab_section, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set) 
                   
              CALL  wf_damping(sh_env, tau, mol, deltatMD, sh_print_env)  
            CASE("COLLAPSE_WITHIN_DECOTIME")
              CALL calc_kinetic_tot(sh_env, velocities_com_tp, local_particles, atomic_kinds,&
                                    nparticle_kind, atomic_kind_set, K_tot, sh_env%ext_pres)
              CALL eval_deco_tau(sh_env,diagonal_F, tau, adiab_section, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set) 
                   
              CALL  collapse_within_deco_time(tau, mol, deltatMD, sh_env)
            CASE("IDA")
               CALL ida_collapse(sh_env)
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE("ADIABATIC_BASIS")
            SELECT CASE(deco_keyword)                                                          
            CASE("INSTANT_COLLAPSE")
              CALL wf_collapse_adiab_prop(sh_env,sh_print_env, itimes)
            CASE("TRESH_ONLY_COLLAPSE")  
              CPABORT("NOT IMPLEMENTED")
            CASE("DAMPING")
              CALL calc_kinetic_tot(sh_env, velocities_com_tp, local_particles, atomic_kinds,&
                                    nparticle_kind, atomic_kind_set, K_tot, sh_env%ext_pres)
              CALL eval_deco_tau(sh_env,diagonal_F, tau, adiab_section, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds, nparticle_kind, atomic_kind_set) 
              CALL wf_damping_adiab_prop(sh_env, tau, mol, deltatMD, sh_print_env)                    
            CASE DEFAULT                                           
              CPABORT("There is no method to add the decoherence.") 
            END SELECT
        CASE DEFAULT
           CPABORT("Choose REPRESENTATION, only: DIABATIC_BASIS or ADIABATIC_BASIS")
        END SELECT
  END SUBROUTINE decohere_wf


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE wf_collapse_adiab_prop(sh_env,sh_print_env, itimes)
    TYPE(sh_env_type)                             :: sh_env
    TYPE(sh_print_env_type)                        :: sh_print_env
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse_adiab_prop', &         
    routineP = moduleN//':'//routineN                                                                 

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (sh_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
!The wafunction is collapsed in a pure state
                     sh_env%ctin(:) =  CMPLX(0.0D0, 0.0D0,dp)
                     sh_env%ctin(sh_env%active_state) = CMPLX(1.0D0, 0.0D0,dp)
                  sh_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_collapse_adiab_prop


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE wf_collapse(sh_env,sh_print_env, itimes)
    TYPE(sh_print_env_type)                        :: sh_print_env
    TYPE(sh_env_type)                             :: sh_env
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
    

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (sh_env%futur_decohere .and. sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
                  sh_env%futur_decohere = .FALSE.                         
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_collapse 


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE ida_collapse(sh_env)
    TYPE(sh_env_type)                             :: sh_env
    CHARACTER(len=*), PARAMETER :: routineN ='ida_collapse', &         
    routineP = moduleN//':'//routineN                                                                 
        
       if (sh_env%attempt_hop) then
              sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
      end if      
  END SUBROUTINE ida_collapse

 
!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE wf_tresh_collaps(sh_env, sh_print_env, itimes)
    TYPE(sh_env_type)                             :: sh_env
    TYPE(sh_print_env_type)                        :: sh_print_env
    INTEGER, POINTER                               :: itimes

    CHARACTER(len=*), PARAMETER :: routineN ='wf_collapse', &         
    routineP = moduleN//':'//routineN                                                                 

 sh_print_env%decohere = .FALSE.                                      
      IF (itimes.ge.1) THEN                                                
         !Finally, collapse the wavefunction if called for                 
         if (sh_env%decohere) then                                  
!> Check if the surface_hop subroutine asks us to do the decoherence       
                  sh_print_env%decohere = .TRUE.
                     sh_env%ctin(:) = sh_env%C_tp(:, sh_env%active_state)
         end if                                                            
         CALL adiab_print_sh_info(sh_print_env, keyword = 'DECOHERENCE')   
      ENDIF
  END SUBROUTINE wf_tresh_collaps


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE  wf_damping_adiab_prop(sh_env, tau, mol, deltatMD, sh_print_env)                    
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping_adiab_prop', &         
    routineP = moduleN//':'//routineN   
                                                        
    INTEGER  :: i
    REAL(KIND=dp)      :: sum 
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST

    sum = 0.0D0


    !>transfor the diabatic into the adiabatic coefficients and calculate tau
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
!               tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
!                                        (C_param +(E_0_param/K_tot))

    !>Damp the coefficients of the states different from the active state    
                sh_env%ctin(i) = sh_env%ctin(i)*exp(-(deltatMD/tau(i)))
                sum = sum + abs(sh_env%ctin(i))*abs(sh_env%ctin(i))
           end if
    end do

    !>transfer population to the active state for keeping the detailed balance   
    if (sh_env%ctin(sh_env%active_state)== 0) then
           sh_env%ctin(sh_env%active_state) = sqrt(1-sum)
    else 
           sh_env%ctin(sh_env%active_state)= sh_env%ctin(sh_env%active_state)* &
              sqrt((1-sum)/(abs(sh_env%ctin(sh_env%active_state))*abs(sh_env%ctin(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = sh_env%ctin(sh_env%active_state)*conjg(sh_env%ctin(sh_env%active_state))
    deallocate(pop)
  END SUBROUTINE wf_damping_adiab_prop 

  
!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE  eval_deco_tau(sh_env, diagonal_F, tau, adiab_section, mol, K_tot, sh_print_env, &
                          local_particles, atomic_kinds,nparticle_kind, atomic_kind_set)                    

    INTEGER, INTENT(IN)                                :: mol
    TYPE(sh_env_type)                                  :: sh_env
    REAL(kind=dp)                                      :: K_tot
    TYPE(sh_print_env_type)                            :: sh_print_env
    TYPE(distribution_1d_type), POINTER                :: local_particles
    TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
                POINTER, OPTIONAL                      :: atomic_kind_set
    INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: tau
    TYPE(section_vals_type), POINTER                   :: adiab_section
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER           :: diagonal_F



    CHARACTER(len=*), PARAMETER :: routineN ='eval_deco_tau', &         
    routineP = moduleN//':'//routineN   

    CHARACTER(LEN=default_string_length)           :: tau_param
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: tau_forces, tau_pure_deph
    REAL(KIND=dp)                                  :: C_param, E_0_param
    TYPE(atomic_kind_type), ALLOCATABLE, TARGET, DIMENSION(:)     :: my_atomic_kind_set
    INTEGER                                        :: my_nparticle_kind
    REAL(kind=dp)                                  :: mass
    INTEGER                                        :: iparticle_kind, orb
    INTEGER                                        :: nparticle_local
    REAL(kind=dp)                                  :: mixing_threshold
    INTEGER                                        :: iparticle_local, iparticle
    TYPE(atomic_kind_type), POINTER                :: atomic_kind
    INTEGER  :: i, deg_of_freedom
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE       :: adiab_nace_vector
    REAL(KIND=dp)                                  :: KbT, Temp, alpha, T_inp 
    !REAL(kind=dp), DIMENSION(:,:),POINTER          :: active_state_adforces
    REAL(kind=dp), DIMENSION(:,:),POINTER          :: my_forces_adiab
    LOGICAL                                        :: tau_schwartz
    TYPE(section_vals_type), POINTER               :: decoherence_section
    REAL(KIND=dp)                                  :: convert_nacv

    NULLIFY(decoherence_section)
 
    !ALLOCATE(active_state_adforces(3, sh_env%natom))
    ALLOCATE(my_forces_adiab(3,sh_env%natom))
    allocate(tau(mol)) ! decoherence time (EDC method)
    allocate(tau_forces(mol)) ! decoherence time (EDC method)
    allocate(tau_pure_deph(mol)) ! decoherence time (EDC method)
    allocate(sh_print_env%tau_deco(mol))
    allocate(adiab_nace_vector(mol))

    decoherence_section => section_vals_get_subs_vals(adiab_section,"DECOHERENCE")

    adiab_nace_vector(:) = sh_env%adiab_nace_tp(sh_env%active_state,:)
    !>number of degrees of freedom intensive EDC
    deg_of_freedom = (COUNT(sh_env%ext_pres)*3-6)

    CALL section_vals_val_get(adiab_section,"NUMBER_ORBITALS",&
         i_val=orb)


    IF (PRESENT(atomic_kinds)) THEN
       my_nparticle_kind = atomic_kinds%n_els
       ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
       my_atomic_kind_set = atomic_kinds%els
    ELSE
       IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
          my_nparticle_kind = nparticle_kind
          ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
          my_atomic_kind_set = atomic_kind_set
       END IF
    END IF

!TEST TAU_SCHARTZ
!tau_schwartz = .TRUE.
!END TEST
   ! COMMENT_SG: HARDCODED stuff just for testing purposes
   tau_schwartz = .FALSE.
!   sh_print_env%pure_deph(:) = .FALSE.


   !> transfor the diabatic into the adiabatic coefficients and calculate tau
   CALL section_vals_val_get(decoherence_section,"DECO_TIME",c_val=tau_param)
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
              SELECT CASE(tau_param)
              CASE("EDC")
                    CALL section_vals_val_get(decoherence_section,"EDC_C",r_val=C_param)
                    CALL section_vals_val_get(decoherence_section,"EDC_E0",r_val=E_0_param)
                    !TEST WITHOUT KINETIC: useful if statment when E_0=0 and FROZEN HAMILTONIAN are
                    !used
                    if (E_0_param == 0.0D0) then  
                         K_tot = 1.0D0
                    end if
                    tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))* &
                                        (C_param+((E_0_param/K_tot)*deg_of_freedom))
              CASE("FORCES_BASED")
                    !Boltzmann factor to use in alpha parameter
                     KbT = (2*K_tot)/(3*COUNT(sh_env%ext_pres))
                    !instead of using the exact boltzmann factor, at the moment the code uses Temp =
                    !298 K, maybe this point should be modified in future
                    
                    CALL section_vals_val_get(decoherence_section,"THRESHOLD_TAU_FORCES",r_val=mixing_threshold)
                    CALL section_vals_val_get(decoherence_section,"TEMPERATURE_FG_WIDTH",r_val=T_inp)
                    
                    Temp = 3.1693191275167785e-06*T_inp !conver in KbT
                      !Temp = 0.0009444571  !KbT at 298 K
                      ! mixing_threshold = 2.0D0
                    if ((abs(adiab_nace_vector(i))/(abs(sh_env%E_tp(i)-sh_env%E_tp(sh_env%active_state)))) & 
                                         .ge. mixing_threshold) then
                         !if there is mixing then we use force otherwise we use !pure dephasing
                         if (sh_env%aom_couplings_keyword) then
                              ! forces at the moment can be calculate with AOM Hamiltonian only
                              convert_nacv = -sh_env%scaling/sh_env%cbar
                              call calcforces_smart_grad(i , sh_env%C_tp, &
                                    diagonal_F, my_forces_adiab, &
                                    sh_env%ext_pres, sh_env%active_system_info, orb, &
                                       nacv=sh_env%nacv, convert_nacv=convert_nacv)
                         else if (sh_env%tresp_couplings_keyword) then
                              call calcforces_smart_grad(i , sh_env%C_tp, &
                                    diagonal_F, my_forces_adiab, &
                                    sh_env%ext_pres, sh_env%active_system_info, orb, &
                                       nacv=sh_env%grad_tresp)
                         else if (sh_env%frozen_couplings_keyword) then
                              call calcforces_smart_grad(i , sh_env%C_tp, &
                                    diagonal_F, my_forces_adiab, &
                                    sh_env%ext_pres, sh_env%active_system_info, orb, &
                                       frozen_c=.TRUE.)
                         else
                            CPABORT("There is no method to calculate adiabatic forces for this H")
                         end if
                         tau(i) = 0.0D0
                         DO iparticle_kind= 1, my_nparticle_kind
                                atomic_kind => my_atomic_kind_set(iparticle_kind)
                                CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
                                nparticle_local = local_particles%n_el(iparticle_kind)
                                DO iparticle_local= 1, nparticle_local
                                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                                   alpha = 6*mass*Temp                  
                                   if (tau_schwartz) then    
                                        tau(i) = tau(i) + &
                                        (dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle))))/(4*alpha)
                                   else
                                        tau(i) = tau(i) + &
                                        (sqrt(dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)))))/(2*sqrt(alpha))
                                   end if
                                END DO
                         END DO
                         if (tau_schwartz) then
                             tau(i) = sqrt(1/tau(i))
                         else
                             tau(i) = 1/tau(i)
                         end if
                    else
!                         sh_print_env%pure_deph(i) = .TRUE.
                         tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))
                    end if

              CASE("HYBRID_TAU")
                    !Boltzmann factor to use in alpha parameter
                     KbT = (2*K_tot)/(3*COUNT(sh_env%ext_pres))
                    !instead of using the exact boltzmann factor, at the moment the code uses Temp =
                    !298 K, maybe this point should be modified in future
                    
                    CALL section_vals_val_get(decoherence_section,"THRESHOLD_TAU_FORCES",r_val=mixing_threshold)
                    CALL section_vals_val_get(decoherence_section,"TEMPERATURE_FG_WIDTH",r_val=T_inp)
                    
                    Temp = 3.1693191275167785e-06*T_inp !conver in KbT
                      !Temp = 0.0009444571  !KbT at 298 K
                      ! mixing_threshold = 2.0D0
                    if ((abs(adiab_nace_vector(i))/(abs(sh_env%E_tp(i)-sh_env%E_tp(sh_env%active_state)))) & 
                                         .ge. mixing_threshold) then
                         !if there is mixing then we use force otherwise we use !pure dephasing
                         if (sh_env%aom_couplings_keyword) then
                            ! forces at the moment can be calculate with AOM Hamiltonian only
                             convert_nacv = -sh_env%scaling/sh_env%cbar
                             call calcforces_smart_grad(i, sh_env%C_tp, &
                                  diagonal_F, my_forces_adiab, &
                                  sh_env%ext_pres, sh_env%active_system_info, orb, &
                                  nacv=sh_env%nacv, convert_nacv=convert_nacv)
                         else if (sh_env%tresp_couplings_keyword) then
                              call calcforces_smart_grad(i , sh_env%C_tp, &
                                    diagonal_F, my_forces_adiab, &
                                    sh_env%ext_pres, sh_env%active_system_info, orb, &
                                       nacv=sh_env%grad_tresp)
                         else if (sh_env%frozen_couplings_keyword) then
                              call calcforces_smart_grad(i , sh_env%C_tp, &
                                    diagonal_F, my_forces_adiab, &
                                    sh_env%ext_pres, sh_env%active_system_info, orb, &
                                       frozen_c=.TRUE.)
                         else
                            CPABORT("There is no method to calculate adiabatic forces for this H")
                         end if
                         tau(i) = 0.0D0
                         tau_forces(i) =  0.0D0
                          tau_pure_deph(i) = 0.0D0
                         DO iparticle_kind= 1, my_nparticle_kind
                                atomic_kind => my_atomic_kind_set(iparticle_kind)
                                CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
                                nparticle_local = local_particles%n_el(iparticle_kind)
                                DO iparticle_local= 1, nparticle_local
                                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                                   alpha = 6*mass*Temp                  
                                   if (tau_schwartz) then    
                                        tau(i) = tau(i) + &
                                        (dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle))))/(4*alpha)
                                   else
                                        tau(i) = tau(i) + &
                                        (sqrt(dot_product((my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)),&
                                        (my_forces_adiab(:,iparticle)-sh_env%forces_adiab(:,iparticle)))))/(2*sqrt(alpha))
                                   end if
                                END DO
                         END DO
                         if (tau_schwartz) then
                             tau(i) = sqrt(1/tau(i))
                         else
                             tau_forces(i) = 1/tau(i)
                             tau_pure_deph(i) =  1/abs(sh_env%E_tp(i)-sh_env%E_tp(sh_env%active_state))
                            if(tau_forces(i) .ge. tau_pure_deph(i)) then
                                 tau(i) = tau_forces(i)
                            else
                                 tau(i) =  tau_pure_deph(i)
                            end if
                         end if
                    else
!                         sh_print_env%pure_deph(i) = .TRUE.
                         tau(i)= (1/abs(sh_env%E_tp(i) - sh_env%E_tp(sh_env%active_state)))
                    end if
              
              CASE DEFAULT
                  CPABORT("There is no method to choose decoherence time.")
              END SELECT
           else
                tau(i) = 0.0D0
           end if
    end do

    sh_print_env%tau_deco = tau

    
    CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   

    deallocate(tau_forces, tau_pure_deph)
    deallocate(adiab_nace_vector)
    deallocate(sh_print_env%tau_deco)
    DEALLOCATE(my_forces_adiab)
    !DEALLOCATE(active_state_adforces)
!    deallocate(sh_print_env%pure_deph)
  END SUBROUTINE eval_deco_tau


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE collapse_within_deco_time(tau, mol, deltatMD, sh_env)

    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(sh_env_type)                         :: sh_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces

    CHARACTER(len=*), PARAMETER :: routineN ='collapse_within_deco_time', &
    routineP = moduleN//':'//routineN

    INTEGER  :: i, j
    REAL(KIND=dp) :: random_num, collapse_probability
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop
    REAL(KIND=dp)      :: sum  


    allocate(pop(mol)) !TEST
    allocate(ct(mol))  ! adiabatic states

    sum = 0.0D0
    do i=1,mol
         if (i .ne. sh_env%active_state) then
             pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
             ct(i) = cmplx(0.0D0, 0.0D0,dp) 
             !> Tranform coefficients in the adiabatic basis 
             do j=1,mol
                    ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
             end do
             pop(i) = ct(i)*conjg(ct(i))
             !print *, "POPU_BEFORE", pop(i)
        
             !>evaluate which states will be collpased
             collapse_probability = deltatMD/tau(i)
             random_num = next_random_number(sh_env%rng_stream)
             if (random_num .gt. collapse_probability) then 
                    ct(i) = ct(i)  !nothing happens
             else
                    ct(i) = 0.0D0
             end if
             sum = sum + abs(ct(i))*abs(ct(i))
         end if
    end do
 
    !>transfor the active diabatic state into the active adiabatic state
    pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
    ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    do i=1,mol                                                       
           ct(sh_env%active_state) = ct(sh_env%active_state) + &
             cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    end do  
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
    
    !>transfer population to the active state for keeping the detailed balance   
    if (ct(sh_env%active_state)== 0) then
           ct(sh_env%active_state) = sqrt(1-sum)
    else 
           ct(sh_env%active_state)= ct(sh_env%active_state)* &
              sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))
    end if

    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))

    !>transfor back the adiabatic into the diabatic which, then, follow the
      !electronic prpagation    
    do j=1,mol
          sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
          pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
          do i=1,mol 
              sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
          end do
         pop(j)=sh_env%ctin(j)*conjg(sh_env%ctin(j))
    end do

    deallocate(pop) !TEST
    deallocate(ct)  ! adiabatic states

  END SUBROUTINE collapse_within_deco_time


!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE  wf_damping(sh_env, tau, mol, deltatMD, sh_print_env)                    
    REAL(KIND=dp), DIMENSION(:), INTENT(in)    :: tau!, tau_forces
    TYPE(sh_env_type)                          :: sh_env
    INTEGER, INTENT(IN)                        :: mol
    REAL(KIND=dp), INTENT(IN)                  :: deltatMD
    TYPE(sh_print_env_type)                    :: sh_print_env

    CHARACTER(len=*), PARAMETER :: routineN ='wf_damping', &         
    routineP = moduleN//':'//routineN   

    INTEGER  :: j,i
    REAL(KIND=dp)      :: sum  
    COMPLEX(KIND=dp), dimension(:), allocatable :: ct
    COMPLEX(KIND=dp), dimension(:), allocatable :: pop

    allocate(pop(mol)) !TEST
    allocate(ct(mol))  ! adiabatic states

    sum = 0.0D0
    do i=1,mol 
           if (i .ne. sh_env%active_state) then
               pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
               ct(i) = cmplx(0.0D0, 0.0D0,dp) 
               !> if statement to choose the decoherence time

               !> Tranform coefficients in the adiabatic basis 
                do j=1,mol
                        ct(i) = ct(i) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*sh_env%ctin(j)
                end do
                !>damp the coefficients of the states different from the active state    
                pop(i) = ct(i)*conjg(ct(i))
                !ct(i) = ct(i)
                ct(i) = ct(i)*exp(-(deltatMD/tau(i)))
                pop(i) = ct(i)*conjg(ct(i))
                sum = sum + abs(ct(i))*abs(ct(i))
           end if
    end do

    !>transfor the active diabatic state into the active adiabatic state
    pop(sh_env%active_state)= cmplx(0.0D0, 0.0D0, dp)
    ct(sh_env%active_state) = cmplx(0.0D0, 0.0D0, dp)
    do i=1,mol                                                       
           ct(sh_env%active_state) = ct(sh_env%active_state) + &
             cmplx(sh_env%C_tp(i,sh_env%active_state), 0.0D0,dp)*sh_env%ctin(i)
    end do  
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))
    
    !>transfer population to the active state for keeping the detailed balance   
    if (ct(sh_env%active_state)== 0) then
           ct(sh_env%active_state) = sqrt(1-sum)
    else 
           ct(sh_env%active_state)= ct(sh_env%active_state)* &
              sqrt((1-sum)/(abs(ct(sh_env%active_state))*abs(ct(sh_env%active_state))))
    end if
    
    pop(sh_env%active_state) = ct(sh_env%active_state)*conjg(ct(sh_env%active_state))

    !>transfor back the adiabatic into the diabatic which, then, follow the
      !electronic prpagation    
    do j=1,mol
          sh_env%ctin(j) = cmplx(0.0D0, 0.0D0,dp)                           
          pop(j) = cmplx(0.0D0, 0.0D0,dp) !TEST
          do i=1,mol 
              sh_env%ctin(j) = sh_env%ctin(j) + cmplx(sh_env%C_tp(j,i), 0.0D0,dp)*ct(i)
          end do
         pop(j)= sh_env%ctin(j)*conjg(sh_env%ctin(j))
    end do
    
!   CALL adiab_print_sh_info(sh_print_env, keyword = 'DAMPING')   
    
    deallocate(ct) 
    deallocate(pop)
  END SUBROUTINE wf_damping 
 

!**********************************************************************************
!> \brief TOADD
!**********************************************************************************
  SUBROUTINE calc_kinetic_tot(sh_env, array, local_particles, atomic_kinds, &
                                              nparticle_kind, atomic_kind_set, K_tot, ext_pres)
      
        TYPE(sh_env_type)                                 :: sh_env
        REAL(kind=dp), DIMENSION(:,:), POINTER             :: array
        TYPE(distribution_1d_type), POINTER                :: local_particles
        TYPE(atomic_kind_list_type), POINTER, OPTIONAL     :: atomic_kinds
        INTEGER, INTENT(IN), OPTIONAL                      :: nparticle_kind
        TYPE(atomic_kind_type), DIMENSION(:), &
             POINTER, OPTIONAL                             :: atomic_kind_set
        REAL(kind=dp), INTENT(OUT)                         :: K_tot
        LOGICAL, dimension(:), INTENT(IN), OPTIONAL           :: ext_pres

        CHARACTER(len=*), PARAMETER :: routineN = 'calc_kinetic_tot', &
        routineP = moduleN//':'//routineN

        TYPE(atomic_kind_type), ALLOCATABLE, TARGET, DIMENSION(:)     :: my_atomic_kind_set
        INTEGER                                  :: my_nparticle_kind
        REAL(kind=dp)                            :: mass
        INTEGER                                  :: iparticle_kind
        INTEGER                                  :: nparticle_local
        INTEGER                                  :: iparticle_local, iparticle
        TYPE(atomic_kind_type), POINTER          :: atomic_kind
        logical                                  :: qsys_only

           IF (PRESENT(atomic_kinds)) THEN
              my_nparticle_kind = atomic_kinds%n_els
              ALLOCATE(my_atomic_kind_set(SIZE(atomic_kinds%els)))
              my_atomic_kind_set = atomic_kinds%els
           ELSE
              IF ((PRESENT(nparticle_kind)).AND.(PRESENT(atomic_kind_set))) THEN
                 my_nparticle_kind = nparticle_kind
                 ALLOCATE(my_atomic_kind_set(SIZE(atomic_kind_set)))
                 my_atomic_kind_set = atomic_kind_set
              END IF
           END IF

           qsys_only = .FALSE.
           if (PRESENT(ext_pres)) then
               !ext_pres =  sh_env%ext_pres 
               qsys_only = .TRUE.
           else 
               qsys_only = .FALSE.
           end if
     
           K_tot = 0.0D0
           DO iparticle_kind= 1, my_nparticle_kind
              NULLIFY(atomic_kind)
              atomic_kind => my_atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
!calculate kinetic energy according to the COM since the array is velocities_COM
              DO iparticle_local= 1, nparticle_local
                 if (qsys_only) then 
                      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                      if (ext_pres(iparticle)) then
                        K_tot = K_tot + mass*( &
                                   array(1, iparticle)**2 + &
                                   array(2, iparticle)**2 + &
                                   array(3, iparticle)**2  ) 
                      end if
                 else 
!print *, "NOT_PRESENT"
                     iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                     K_tot = K_tot + mass*( &
                                array(1, iparticle)**2 + &
                                array(2, iparticle)**2 + &
                                array(3, iparticle)**2  )
                 end if 
              END DO
           END DO
           K_tot = 0.5*K_tot
   END SUBROUTINE calc_kinetic_tot
 
!**********************************************************************************
!> \brief SG this subroutine just replicates the propagate_wf subroutine, but
!it does not calculate diabatic nace at each electronic step from S inverted
!matrix. However it interpolates diabatic nace between nuclear time steps in the
!same way as for the H
!**********************************************************************************
   SUBROUTINE propagate_wf_diab_nace_interpol(ctin, sh_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)

     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     TYPE(sh_env_type)                              :: sh_env
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp), INTENT(IN)                       :: times 
     integer                                         :: mat 


     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf_nace_interpol', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctotal
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: my_H_t, my_H_tp
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_2, X1, X2, X12
     REAL(KIND=dp), dimension(:,:), allocatable       :: Ddiff, D_0, D_1, D_2

     LOGICAL                                         :: renorm, offset
     REAL(KIND=dp)                                   :: ctmag
     COMPLEX(KIND=dp)                                ::  ctcomp
     integer                                         :: t2,i


     !Various Allocation
     allocate(ctotal(mat))

     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     allocate(Ddiff(mat,mat))
     allocate(D_0(mat,mat))
     allocate(D_1(mat,mat))
     allocate(D_2(mat,mat))
     allocate(Hdiff(mat,mat))
     allocate(H_1(mat,mat))
     allocate(H_2(mat,mat))
     allocate(X1(mat,mat))
     allocate(X2(mat,mat))
     allocate(X12(mat,mat))
     allocate(my_H_t(mat, mat))
     allocate(my_H_tp(mat, mat))

      ! Control the offset of hamiltonian and renormalisation of wf
      ! See Hammes-Schiffer, Tully 1994 (Proton transfer in solution) for offset
      offset = .TRUE.
      renorm = .TRUE.
     
    
     my_H_t = sh_env%H_t
     my_H_tp = sh_env%H_tp
     if (offset) then
        do i = 1, mat
           my_H_t(i,i) = my_H_t(i,i) - sh_env%H_t(1,1)
           my_H_tp(i,i) = my_H_tp(i,i) - sh_env%H_tp(1,1)
        enddo
     endif
      

!COMMENT: INDENTATION
            
             Hdiff = (deltatE/deltatMD)*(my_H_tp - my_H_t)
             Ddiff = (deltatE/deltatMD)*(sh_env%diabatic_nace_tp -sh_env%diabatic_nace_t)
             D_0  = sh_env%diabatic_nace_t

             X1 = my_H_t - D_0

             do t2=1,int(deltatMD/deltatE)
                    !create timestepped equivalents
                    !> interpolate hamiltonian and diabatic NACE (that are in
                    !any case very smallcompare to H)
                    H_1 = my_H_t + 0.5D0*Hdiff
                    D_1 = D_0 + 0.5D0*Ddiff
                    X12 = H_1 - D_1
     
                    H_2 = H_1 + 0.5D0*Hdiff
                    D_2 = D_1 + 0.5D0*Ddiff
                    X2 = H_2 - D_2
     
                    ! Call new (vector) RK4 method
                    CALL RK4(ctin, X1, X12, X2, deltatE, mat)
                    !ctout = ctin(:)

                    ctmag = 0.0D0
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)
                    IF ( ABS( ctmag - 1).GE.1E-5 ) THEN
                       print*, "FOBSH| *** WARNING: Norm is", ctmag
                       CPABORT("Problem in the electronic propagation, norm is not conserved.")
                    ENDIF
  
                    ! Check WF for norm conservation
                    if (renorm) then
                       !rescale to conserve norm = 1
                       do i=1,mat
                            ctin(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctin(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctin(i)),dp)
                       end do
                    endif
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                    end do

                    !update
                     my_H_t = H_2
                     D_0 = D_2
                     X1  = X2
             end do

             CALL adiab_print_nace(adiab_section, D_2, nadiab, norbitals, itimes, deltatMD, times)
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                    times)
             CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                    times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                    times)
       
!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP

       deallocate( Hdiff, H_1, H_2, X1, X2, X12)
       deallocate( my_H_t, my_H_tp, ctotal)
       deallocate(D_0, Ddiff, D_1, D_2)

  END SUBROUTINE propagate_wf_diab_nace_interpol 


  !************************************************************************
  !> \brief The RK4 algorithm for propagating coeffs. This is implemented
  !>        the way it is described on Wikipedia 
  !>        (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).
  !>        The coeffs are multiplied by the matrices to calculate K to
  !>        avoid mat-mat multiplications.
  !> \author Matt
  !************************************************************************   
  SUBROUTINE RK4(ctin, X1, X12, X2, deltaTe, mat) 
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)  :: ctin 
      REAL(KIND=dp)                                :: deltaTe
      INTEGER                                      :: mat
      COMPLEX(KIND=dp), DIMENSION(:,:)              :: X1, X2, X12
 
      CHARACTER(len=*), PARAMETER                  :: routineN='RK4', &
                                                      routineP=routineN//';'//moduleN

      CHARACTER(len=1)                             :: zgemm_type
      COMPLEX(KIND=dp)                             :: alpha, beta 
      INTEGER                                      :: inc
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: K1, K2, K3, K4, Ktmp                

      ALLOCATE(Ktmp(mat), K1(mat), K2(mat), K3(mat), K4(mat))

      inc = 1
      alpha = CMPLX(1, 0, dp)
      beta = CMPLX(0, 0, dp)
      zgemm_type = 'N'

      ! Calculate the K coefficients in the RK4 algorithm
      ! Each K is used to determine the next K, e.g K1 determines K2 etc...
      
      ! ZGEMV is BLAS subroutine for complex matrix-vector multiplication
      ! Z = Complex, GE = general, M = matrix, V = vector
      CALL ZGEMV(zgemm_type, mat, mat, alpha, X1, mat, ctin, inc, beta, K1, inc) 
      K1 = K1 * deltaTe

      Ktmp = ctin + K1/2.0
      CALL ZGEMV(zgemm_type, mat, mat, alpha, X12, mat, Ktmp,&
                 inc, beta, K2, inc) 
      K2 = K2 * deltaTe
     
      Ktmp = ctin + K2/2.0
      CALL ZGEMV(zgemm_type, mat, mat, alpha, X12, mat, Ktmp,&
                 inc, beta, K3, inc) 
      K3 = K3 * deltaTe
     
      Ktmp = ctin + K3
      CALL ZGEMV(zgemm_type, mat, mat, alpha, X2, mat, Ktmp,&
                 inc, beta, K4, inc) 
      K4 = K4 * deltaTe
     
      ! Ktmp = the final K value to add to the coefficient at the end
      Ktmp = (K1 + (2.0 * K2) + (2.0 * K3) + K4) / 6.0
      
      ctin = ctin + Ktmp

      DEALLOCATE(Ktmp, K1, K2, K3, K4)
  END SUBROUTINE RK4


  !**********************************************************************************
  !> \brief Similar to the propagate wf subroutine. However, the diabatic NACE
  !         isn't calculated/used. The diabatic NACE is very small compared to
  !         the coupling for situations in which SH is sensible to use (i.e.
  !         outside of the hopping regime) so this subroutine is applicable to
  !         the majority of sims.
  ! \author Sam (Matt edited RK4)
  ! \date   23/05/19
  !**********************************************************************************
  SUBROUTINE propagate_wf_fast(ctin, sh_env, nadiab,  norbitals, ELsteps, &
                          itimes, adiab_section, deltatMD, deltatE , times,      &
                           mat)
 
      COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
      TYPE(sh_env_type)                              :: sh_env
      INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
      INTEGER                                         :: ELsteps   
      INTEGER, POINTER                                :: itimes
      TYPE(section_vals_type), POINTER                :: adiab_section
      REAL(KIND=dp)                                   :: deltatMD, deltatE
      REAL(KIND=dp), INTENT(IN)                       :: times 
      integer                                         :: mat 
 
 
      CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf_fast', &
      routineP = moduleN//':'//routineN
 
      COMPLEX(KIND=dp), dimension(:), allocatable     :: ctotal
      COMPLEX(KIND=dp), dimension(:,:), allocatable   :: my_H_t, my_H_tp
      COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Hdiff, X1, X2, X12
 
      LOGICAL                                         :: renorm, offset
      REAL(KIND=dp)                                   :: ctmag
      !COMPLEX(KIND=dp)                                :: alpha, beta, ctcomp
      COMPLEX(KIND=dp)                                :: ctcomp
      integer                                         :: t2
      integer                                         :: i
 
      ALLOCATE(ctotal(mat),       Hdiff(mat,mat),  &
               X1(mat,mat),       X2(mat,mat),     &
               X12(mat,mat),      my_H_t(mat, mat),&
               my_H_tp(mat, mat))
      
      ! Control the offset of hamiltonian and renormalisation of wf
      ! See Hammes-Schiffer, Tully 1994 (Proton transfer in solution) for offset
      offset = .TRUE.
      renorm = .TRUE. 
 
      ! Offset the Hamiltonian (I've left some unnecessary arrays here for
      !                         clarity)
      my_H_t = sh_env%H_t
      my_H_tp = sh_env%H_tp
      if (offset) then
         do i = 1, mat
            my_H_t(i,i) = my_H_t(i,i) - sh_env%H_t(1,1)
            my_H_tp(i,i) = my_H_tp(i,i) - sh_env%H_tp(1,1)
         enddo
      endif
 
      ! Difference for linear interpolation
      Hdiff = (deltatE/deltatMD)*(my_H_tp - my_H_t)
 
      ! Now make the X matrices for each electronic half-timestep
      ! X is the matrix used in the propagation equation du/dt = Xu
      !  where u is the coeff, d/dt is an ord. time deriviative.
      ! (see Spencer, 17)
 
      ! Make X for intial timestep
      X1 = my_H_t
      do t2=1,int(deltatMD/deltatE)
          ! Make X for half timestep
          X12 = X1 + 0.5D0*Hdiff ! linearly interpolate hamiltonian
          
          ! Make X for full timestep
          X2 = X12 + 0.5D0*Hdiff ! lin. interp. ham.
      
          ! Do the RK4 algorithm (see wiki)
          CALL RK4(ctin, X1, X12, X2, deltaTe, mat) 
          
          ! Normalisation stuff
          ctmag = 0.0D0
          ctcomp = CMPLX(0.0D0, 0.0D0,dp)
          do i=1,mat
               ctotal(i) = CONJG(ctin(i))*ctin(i)
               ctcomp = ctcomp + ctotal(i)
          end do
   
          ctmag = abs(ctcomp) 
          IF ( ABS( ctmag - 1).GE.1E-5 ) THEN
             PRINT *, "FOBSH| *** WARNING: Norm is", ctmag
             CPABORT("Problem in the electronic propagation, norm is not conserved.")
          ENDIF
           
          ! Check WF for norm conservation
          if (renorm) then
             !rescale to conserve norm = 1
             do i=1,mat
                  ctin(i) = CMPLX((1.0D0/SQRT(ctmag))*REAL(ctin(i)), (1.0D0/SQRT(ctmag))*AIMAG(ctin(i)),dp)
             end do
          endif

          do i=1,mat
               ctotal(i) = conjg(ctin(i))*ctin(i)
          end do

 
          !update
          X1  = X2
      end do
 
 !     CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
      CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                             times)
      CALL adiab_print_pop(adiab_section, ABS(ctotal) ,nadiab,norbitals,itimes, deltatMD, &
                             times)
      CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                             times)
 
      DEALLOCATE(ctotal, Hdiff, X1, X2, X12, my_H_t, my_H_tp)
 
   END SUBROUTINE propagate_wf_fast   
 

  !**********************************************************************************
  !> \brief This subroutine propagates the WF interpolating the overlap matrices
  !         S at each electronic steps and the Hamiltonian at each nuclear
  !         timestep. The diabatic NACE are reconstructed from Overlap matrices
  !         S at each electronic step, making this unusefully accurate and by
  !         far too expensive in term of computer time.
  ! \author Jacob (Sam editted)
  ! \date   23/05/19
  !**********************************************************************************
   SUBROUTINE propagate_wf(ctin, sh_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)

     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     TYPE(sh_env_type)                              :: sh_env
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp), INTENT(IN)                       :: times 
     integer                                         :: mat 


     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_wf', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: my_H_t, my_H_tp
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     REAL(KIND=dp), dimension(:,:), allocatable       :: Stt_int, Sttp_int, Stptp_int, &
                                                          Sttm_int, Stmtm_int

     LOGICAL                                         :: renorm, offset
     REAL(KIND=dp)                                   :: ctmag
     COMPLEX(KIND=dp)                                :: ctcomp
     integer                                         :: t2, i
     character(len=1)                                :: zgemm_type


     !Various Allocation

     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))
     allocate(Stt_int(mat,mat)) 
     !allocate(Stpt(mat,mat))
     allocate(Sttp_int(mat,mat))
     allocate(Stptp_int(mat,mat))
     !allocate(Stmt(mat,mat))
     allocate(Sttm_int(mat,mat))
     allocate(Stmtm_int(mat,mat))

     !allocate(ctout(mat))
     allocate(ctotal(mat))

     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     allocate(NACE(mat,mat))
     allocate(Hdiff(mat,mat))
     allocate(H_0(mat,mat))
     allocate(H_1(mat,mat))
     allocate(H_2(mat,mat))
     allocate(X1(mat,mat))
     allocate(X2(mat,mat))
     allocate(X12(mat,mat))
     allocate(my_H_t(mat, mat))
     allocate(my_H_tp(mat, mat))
     
     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

    
     offset = .TRUE.
     my_H_t = sh_env%H_t
     my_H_tp = sh_env%H_tp
     if (offset) then
        do i = 1, mat
           my_H_t(i,i) = my_H_t(i,i) - sh_env%H_t(1,1)
           my_H_tp(i,i) = my_H_tp(i,i) - sh_env%H_tp(1,1)
        enddo
     endif
      

!COMMENT: INDENTATION
             Hdiff = (deltatE/deltatMD)*(my_H_tp - my_H_t)
             Sttdiff = (deltatE/deltatMD)*(sh_env%Stptp - sh_env%Stt)
             Sttpdiff = (deltatE/deltatMD)*(sh_env%Sttp - sh_env%Stmt)
             Sttmdiff = (deltatE/deltatMD)*(sh_env%Stpt - sh_env%Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(sh_env%Stt - sh_env%Stmtm)
             H_0  = my_H_t

             !Assignment internal variabls used in the RK
             
             Stt_int = sh_env%Stt 
             Sttp_int = sh_env%Sttp
             Stptp_int = sh_env%Stptp
             Sttm_int = sh_env%Sttm
             Stmtm_int = sh_env%Stmtm

             call makeX(mat, H_0, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
                                             Stmtm_int, X1, deltatMD, NACE)

             do t2=1,int(deltatMD/deltatE)
                    !create timestepped equivalents
                    H_1 = H_0 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    call makeX(mat, H_1, Stt_int, Sttp_int, Sttm_int, Stptp_int, & 
                                                       Stmtm_int, X12, deltatMD)
     
                    H_2 = H_1 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    call makeX(mat, H_2, Stt_int, Sttp_int, Sttm_int, Stptp_int, Stmtm_int, &
                                                 X2, deltatMD, NACE)
     
                    !K2 = cmplx(0.0D0, 0.0D0,dp)
                    !K3 = cmplx(0.0D0, 0.0D0,dp)
                    !K4 = cmplx(0.0D0, 0.0D0,dp)
          
                    !!now that X has been found, use 4th order Runge-Kutta numerical integration method
                    !!aiming to find operator such that c(t+1) = (1 + operator)*c
          
                    !!there are four terms in X: each is an increasing multiple. so first term is proportional to X(t),
                    !!second term is as X(t+deltat/2)*X(t), third term is as X(t+deltat/2)*X(t+deltat/2)*X(t),
                    !!fourth term is as X(t+deltat)*X(t+deltat/2)*X(t+deltat/2)*X(t)
                    !!this means each term is some X matrix times the previous term,
                    !!so only three more calls to ZGEMM are needed
                    !!the K matrices are the second, third and fourth Runge-Kutta terms
     
                    !alpha = cmplx(deltatE/2.0D0, 0.0D0,dp)
                    !beta = cmplx(0.0D0, 0.0D0,dp)
     
                    !!second term
                    !call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, X1, &
                    !&mat, beta, K2, mat)
                    !K2 = K2 + X12
                    !
                    !!third term
                    !call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X12, mat, K2, &
                    !&mat, beta, K3, mat)
                    !K3 = K3 + X12
     
                    !!fourth term
                    !alpha = cmplx(deltatE, 0.0D0,dp)
                    !call ZGEMM(zgemm_type, zgemm_type, mat, mat, mat, alpha, X2, mat, K3, &
                    !&mat, beta, K4, mat)
                    !K4 = K4 + X2
 
                    !!put everything together
                    !Ktot = Imat + (deltatE/6.0D0)*(X1 + 2*K2 + 2*K3 + K4)
     
                    !timestep electron propagation and print out
                    !alpha = dcmplx(1.0D0, 0.0D0)
                    !alpha = cmplx(1.0D0, 0.0D0,dp)
                    !inc = 1
                    !call ZGEMV(zgemm_type, mat, mat, alpha, Ktot, mat, ctin, inc, beta, ctout, inc)

                    ! Call new (vector) RK4 method
                    CALL RK4(ctin, X1, X12, X2, deltatE, mat)
                    !ctout = ctin(:)

                    ctmag = 0.0D0
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)
                    IF ( ABS( ctmag - 1).GE.1E-5 ) THEN
                       print*, "FOBSH| *** WARNING: Norm is", ctmag
                       CPABORT("Problem in the electronic propagation, norm is not conserved.")
                    ENDIF
                    renorm = .TRUE.
                    if (renorm) then
                       !rescale to conserve norm = 1
                       do i=1,mat
                            ctin(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctin(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctin(i)),dp)
                       end do
                    endif
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                    end do

                    !update
                     H_0 = H_2
                     X1  = X2

             end do
             CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                    times)
             CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                    times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                    times)
       
!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP

!DEALLOCATE
       DEALLOCATE(NACE)
       deallocate(H_0, Hdiff, H_1, H_2, X1, X2, X12)
       deallocate(my_H_t, my_H_tp, ctotal)
       deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff)
       deallocate(Stt_int, Sttp_int, Stptp_int, Sttm_int, Stmtm_int)
!END DEALLOCATE

  END SUBROUTINE propagate_wf    


  !**********************************************************************************
  !> \brief This subroutine propagates the WF using Elimination of forbidden hop
  !         equations see manuscript A.Carof. 
  ! \author Sam
  ! \date   23/05/19
  !**********************************************************************************
  SUBROUTINE prop_wf_HS_corr(ctin, sh_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat, class_forbidden_hop, active_state)

     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     INTEGER, dimension(:), INTENT(IN)               :: class_forbidden_hop
     integer, INTENT(IN)                             :: active_state, mat 

     CHARACTER(len=*), PARAMETER                     :: routineN ='prop_wf_HS_corr', &
     routineP = moduleN//':'//routineN

     COMPLEX(KIND=dp), dimension(:), allocatable     ::  ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: Sttdiff, Sttmdiff, Sttpdiff, Stmtmdiff, Stptpdiff
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACE , modified_diab_nace
      COMPLEX(KIND=dp), dimension(:,:), allocatable   :: my_H_t, my_H_tp
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12!, Ytot
     REAL(KIND=dp), dimension(:,:), allocatable       :: Stt_int, Sttp_int, Stptp_int, &
                                                          Sttm_int, Stmtm_int
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE       :: U_t, U_t1, U_t2, U_tp, U_t0 
     COMPLEX(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_hamil
     REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: interm_umatrix, udot_matrices 
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE         :: Eig_val_t, Eig_val_t1, Eig_val_t2, Eig_val_t0,& 
                                                          Eig_val_tp
     LOGICAL                                         :: renorm, offset

     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     TYPE(sh_env_type)                              :: sh_env
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp)                                   :: ctmag
     REAL(KIND=dp)                                   :: times 
     !COMPLEX(KIND=dp)                                :: alpha, beta, ctcomp
     COMPLEX(KIND=dp)                                :: ctcomp
     integer                                         :: t2,i, dim_
     !integer                                         :: i,j, inc 
     character(len=1)                                :: zgemm_type

!TEST (KEEP FOR PRINTING PURPOSES)
!     COMPLEX(KIND=dp), dimension(:), allocatable     :: pop
!     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctin_adiab
!     REAL(KIND=dp), dimension(:,:), allocatable      :: Sterm1, UU, UDU, Sterm2, UDU_old_nace
!     integer :: k
!END TEST
     !Various Allocation


     !!Imat is an identity matrix of the dimensions of the other square matrices
     !allocate(Imat(mat,mat))
     !!define once, then use whenever
     !Imat = 0.0D0
     !FORALL(j = 1:mat) Imat(j,j) = 1     
     
     !need to allocate these square matrices
     !first the real, S-related ones
     allocate(interm_hamil(3,int(deltatMD/deltatE),mat,mat))
     allocate(interm_umatrix(3,int(deltatMD/deltatE),mat,mat))
     allocate(udot_matrices(3,int(deltatMD/deltatE),mat,mat))
     allocate(Stptpdiff(mat,mat))
     allocate(Sttpdiff(mat,mat))
     allocate(Sttdiff(mat,mat))
     allocate(Sttmdiff(mat,mat))
     allocate(Stmtmdiff(mat,mat))
     allocate(Stt_int(mat,mat)) 
     !allocate(Stpt(mat,mat))
     allocate(Sttp_int(mat,mat))
     allocate(Stptp_int(mat,mat))
     !allocate(Stmt(mat,mat))
     allocate(Sttm_int(mat,mat))
     allocate(Stmtm_int(mat,mat))

     !allocate(ctout(mat))
     allocate(ctotal(mat))

     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     allocate(NACE(mat,mat))
     allocate(Hdiff(mat,mat))
     allocate(H_0(mat,mat))
     allocate(H_1(mat,mat))
     allocate(H_2(mat,mat))
     allocate(X1(mat,mat))
     allocate(X2(mat,mat))
     allocate(X12(mat,mat))

!TEST
!     allocate(ctin_adiab(mat))
!     allocate(pop(mat))
!    ! allocate(Ytot(mat,mat))
!     allocate(Sterm1(mat,mat)) 
!     allocate(UU(mat,mat)) 
!     allocate(UDU(mat,mat)) 
!     allocate(UDU_old_nace(mat,mat)) 
!     allocate(Sterm2(mat,mat)) 
!END TEST
    
     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'

     allocate(my_H_t(mat, mat))
     allocate(my_H_tp(mat, mat))

     offset = .TRUE.
     my_H_t = sh_env%H_t
     my_H_tp = sh_env%H_tp
     if (offset) then
        do i = 1, mat
           my_H_t(i,i) = my_H_t(i,i) - sh_env%H_t(1,1)
           my_H_tp(i,i) = my_H_tp(i,i) - sh_env%H_tp(1,1)
        enddo
     endif

    
      
!COMMENT: INDENTATION
             !Hdiff = (deltatE/deltatMD)*(sh_env%H_tp - sh_env%H_t)
             Hdiff = (deltatE/deltatMD)*(my_H_tp - my_H_t)
             Sttdiff = (deltatE/deltatMD)*(sh_env%Stptp - sh_env%Stt)
             Sttpdiff = (deltatE/deltatMD)*(sh_env%Sttp - sh_env%Stmt)
             Sttmdiff = (deltatE/deltatMD)*(sh_env%Stpt - sh_env%Sttm)
             Stptpdiff = Sttdiff
             Stmtmdiff = (deltatE/deltatMD)*(sh_env%Stt - sh_env%Stmtm)
             !H_0  = sh_env%H_t
             H_0  = my_H_t
             !H_2 = sh_env%H_tp
             H_2 = my_H_tp


            !> LOOP TO BUILD INTERMEDIATE HAMILTONIANS AND CALCULATE 
            !> RESPECTIVE U MATRICES
!COMMENT: U_t, U_tp, Eig_val_t, Eig_val_tp are not already known?
!WITH FIVE ELECTRONIC TIMESTEP, ONE NEEDS THE FOLLOWING HAMILTONIANS AND U:
!H_t, H_t+deltaEL/2, H_t+deltaEL, H_t+1.5deltaEL, H_t+2deltaEL...
!AT THE END: NINE MATRICES, IN THE LOOP BELOW ONE CALCULATES 15 MATRICES. 
!INDEED: interm_hamil(3,1,:,:) = interm_hamil(1,2,:,:), etc.
!PROBABLY interm_hamil and intemr_umatrix could be only: 
!interm( 9, :, :)
            CALL adiabaticity(-AIMAG(my_H_t), U_t, Eig_val_t)
            CALL adiabaticity(-AIMAG(my_H_tp), U_tp, Eig_val_tp)
            do t2=1,int(deltatMD/deltatE)

                   !create 3 others intermediate hamiltonians usefull to build the X
                   !matrix in the following 
                   interm_hamil(1,t2,:,:) = H_0(:,:)
                   CALL adiabaticity(-AIMAG(H_0), U_t0, Eig_val_t0)
                   interm_umatrix(1,t2,:,:) = U_t0(:,:)

                   H_1 = H_0 + 0.5D0*Hdiff
                   interm_hamil(2,t2,:,:) = H_1(:,:)
                   CALL adiabaticity(-AIMAG(H_1), U_t1, Eig_val_t1)
                   interm_umatrix(2,t2,:,:) = U_t1(:,:)

                   H_2 = H_1 + 0.5D0*Hdiff
                   interm_hamil(3,t2,:,:) = H_2(:,:)
                   CALL adiabaticity(-AIMAG(H_2), U_t2, Eig_val_t2)
                   interm_umatrix(3,t2,:,:) = U_t2(:,:)

                   H_0 = H_2
            end do

     
             !> CONSTRUCT Udot MATRICES
             
             !> 1st and 3rd arraies of the super-array:
             udot_matrices(1,1,:,:) = (interm_umatrix(2,1,:,:) - & 
                                    interm_umatrix(1,1,:,:))/(deltatE/2.0) !forwards
            ! udot_matrices(2,1,:,:) = (interm_umatrix(1,3,:,:) - &
             !                             interm_umatrix(1,1,:,:))/(deltatE) !centering
             do i =1,(int(deltatMD/deltatE)-1)
                   udot_matrices(3,i,:,:) =  (interm_umatrix(2,i+1,:,:) - & 
                                                     interm_umatrix(2,i,:,:))/(deltatE)
                   udot_matrices(1,i+1,:,:) = udot_matrices(3,i,:,:)
             end do 
                                                
             udot_matrices(3,int(deltatMD/deltatE),:,:) = (interm_umatrix(3,int(deltatMD/deltatE),:,:) - & 
                                       interm_umatrix(2,int(deltatMD/deltatE),:,:))/(deltatE/2.0) !backwards
             !> 2nd array of the super-array:
             do i =1,int(deltatMD/deltatE)
                   udot_matrices(2,i,:,:) = (interm_umatrix(3,i,:,:) - & 
                                         interm_umatrix(1,i,:,:))/(deltatE) !all centering
             end do



             !H_0  = sh_env%H_t
             !H_2 = sh_env%H_tp

!Assignment internal variabls used in the RK
             
             Stt_int = sh_env%Stt 
             Sttp_int = sh_env%Sttp
             Stptp_int = sh_env%Stptp
             !Stpt_int = Stpt!
             !Stmt_int = Stmt!
             Sttm_int = sh_env%Sttm
             Stmtm_int = sh_env%Stmtm

             do t2=1,int(deltatMD/deltatE)
                 
                    dim_ = 1
                    !CALL adiabaticity(-AIMAG(sh_env%H_t), mat, U_t, Eig_val_t, .false.)
                    CALL calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, &
                                                        Stptp_int, Stmtm_int, deltatMD, NACE)
                    !CALL calc_diabatic_nace(mat, sh_env%H_t, Stt_int, Sttp_int, Sttm_int, &
                     !                                   Stptp_int, Stmtm_int, X1, deltatMD, NACE)
                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)
                    
                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X1)


!                    call makeX(mat, sh_env%H_t, Stt_int, Sttp_int, Sttm_int, Stptp_int, &
!                                             Stmtm_int, X1, deltatMD, NACE)
     
                    !create timestepped equivalents
                    !H_1 = H_0 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    
                    dim_ = 2 
!                    CALL adiabaticity(-AIMAG(H_1), mat, U_t1, Eig_val_t1, .false.)
                    CALL calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, & 
                                                       Stmtm_int, deltatMD, NACE)

                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)
                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X12)

                   ! H_2 = H_1 + 0.5D0*Hdiff
                    Stt_int = Stt_int + 0.5D0*Sttdiff
                    Sttp_int = Sttp_int + 0.5D0*Sttpdiff
                    Sttm_int = Sttm_int + 0.5D0*Sttmdiff
                    Stptp_int = Stptp_int + 0.5D0*Stptpdiff
                    Stmtm_int = Stmtm_int + 0.5D0*Stmtmdiff
                    
                    dim_ = 3  
!                    CALL adiabaticity(-AIMAG(H_2), mat, U_t2, Eig_val_t2, .false.)
                    call calc_diabatic_nace(mat, interm_hamil(dim_,t2,:,:), Stt_int, Sttp_int, Sttm_int, Stptp_int, &
                                                     Stmtm_int, deltatMD, NACE)
                    CALL calc_modified_diab_nace(class_forbidden_hop, NACE, & 
                          udot_matrices, interm_umatrix, modified_diab_nace, mat, active_state, t2, dim_)

                    CALL hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, X2)

                    ! Call new (vector) RK4 method
                    CALL RK4(ctin, X1, X12, X2, deltatE, mat)
                    !ctout = ctin(:)

                    ctmag = 0.0D0
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)

                    IF ( ABS( ctmag - 1).GE.1E-1 ) THEN
                       print*, "Norm is", ctmag
                       CPABORT("Problem in the electronic propagation, norm is not conserved.")
                    ENDIF
                    renorm = .TRUE.
                    if (renorm) then
                       !rescale to conserve norm = 1
                       do i=1,mat
                            ctin(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctin(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctin(i)),dp)
                       end do
                    end if
                    do i=1,mat
                         ctotal(i) = conjg(ctin(i))*ctin(i)
                    end do

                    !update
                    H_0 = H_2

                    !ctin = ctout

                    !ctout = cmplx(0.0D0, 0.0D0,dp)

                    !times = DBLE(itimes) + (t2-1)*1/DBLE(ELsteps)


!!TEST> TO COMPARE OLD AND NEW NACE!

!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, modified_diab_nace, mat, 0.0D0, Sterm1, mat)
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, Sterm1, mat, interm_umatrix(dim_,t2,:,:), mat, 0.0D0, UDU, mat)

!print *, "UDU", UDU
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, NACE, mat, 0.0D0, Sterm2, mat)
!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, Sterm2, mat, interm_umatrix(dim_,t2,:,:), mat, 0.0D0, UDU_old_nace, mat)
!print *, "UDU_old_nace", UDU_old_nace

!call DGEMM('N', 'N', mat, mat, mat, 1.0D0, interm_umatrix(dim_,t2,:,:), mat, & 
!                                         udot_matrices(dim_,t2,:,:), mat, 0.0D0, UU, mat)

!print *, "UU", UU
!print *, "NEW_NACE_IN_ADIAB", UDU+UU
!print *, "OLD_NACE_IN_ADIAB", UDU_old_nace+UU

!!do i=1,mat
!!  do j=1,mat
!!    Sterm2(i,j)=0
!!    do k=1,mat
!!    Sterm2(i,j) = Sterm2(i,j) + interm_umatrix(dim_,t2,i,k)*modified_diab_nace(k,j)
!!    end do
!!  end do
!!end do

!!print *, "Sterm2", Sterm2

!do i = 1,mat
!   pop(i) = cmplx(0.0D0, 0.0D0,dp) !TEST
!   ctin_adiab(i) = cmplx(0.0D0, 0.0D0,dp)
!   do j = 1, mat

!       ctin_adiab(i) = ctin_adiab(i) + cmplx(interm_umatrix(dim_,t2,j,i), 0.0D0,dp)*ctin(j)
!   end do
!   pop(i) = ctin_adiab(i)*conjg(ctin_adiab(i))
!end do

!print *, "ct_ADIAB", ctin_adiab
!print *, "ct_SECOND_STATE", ctin_adiab(2)
!print *, "POP_SECOND", pop(2)
!print *, "POP_SECOND_REAL", REAL(pop(2))

!END_TEST


             end do

             CALL adiab_print_nace(adiab_section, NACE, nadiab, norbitals, itimes, deltatMD, times)
             CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                    times)
             CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                  times)
             CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                                    times)

!END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP

!DEALLOCATE

       DEALLOCATE(NACE)
       deallocate(interm_hamil,interm_umatrix, udot_matrices)
       deallocate(H_0, Hdiff, H_1, H_2, X1, X2, X12)!, Ktot, K2, K3, K4)
       deallocate( my_H_t, my_H_tp, ctotal)
       deallocate(Sttdiff, Stptpdiff, Stmtmdiff, Sttpdiff, Sttmdiff)
       deallocate(Stt_int, Sttp_int, Stptp_int, Sttm_int, Stmtm_int)
       !deallocate(Imat)
!END DEALLOCATE

  END SUBROUTINE prop_wf_HS_corr                                                                

!COMMENT: IF THIS SUBROUTINE IS ONLY TWO LINES, WHY NOT WRITE IT DIRECTLY IN THE
!MAIN ONE
  SUBROUTINE hs_corrected_makeX(interm_hamil, modified_diab_nace, t2, dim_, Ytot)
      
     COMPLEX(KIND=dp), DIMENSION(:,:,:,:), INTENT(IN)     :: interm_hamil
     integer, INTENT(IN)                                  :: t2,dim_ 
     !COMPLEX(KIND=dp), dimension(:,:), allocatable   :: H_1, Hdiff, H_0, H_2, X1, X2, X12
     COMPLEX(KIND=dp), dimension(:,:), intent(out) :: Ytot
     
     REAL(KIND=dp), dimension(:,:)                 :: modified_diab_nace

     modified_diab_nace = (modified_diab_nace - TRANSPOSE(modified_diab_nace)) / 2

!print *, "modified_diab_nace_MOD", modified_diab_nace

     Ytot = interm_hamil(dim_,t2,:,:) - modified_diab_nace


  END SUBROUTINE hs_corrected_makeX
                                                            
  SUBROUTINE calc_modified_diab_nace(class_forbidden_hop, NACE, udot_matrices, interm_umatrix, & 
                                            modified_diab_nace, mat, active_state, t2, dim_)


     REAL(KIND=dp), DIMENSION(:,:,:,:), INTENT(IN)     :: interm_umatrix, udot_matrices 
     INTEGER, dimension(:), INTENT(IN)          :: class_forbidden_hop
     REAL(KIND=dp), dimension(:,:), INTENT(IN)      :: NACE
     REAL(KIND=dp), dimension(:,:), allocatable, INTENT(OUT)     :: modified_diab_nace
     integer, INTENT(IN)                                         :: active_state,mat,t2,dim_ 
     

     integer                                         :: i,j, k, l
     REAL(KIND=dp), dimension(:,:), allocatable      :: modified_ad_nace
     REAL(KIND=dp), dimension(:,:), allocatable      :: rotated_ad_nace
     
     

     allocate(modified_ad_nace(mat,mat))
     allocate(rotated_ad_nace(mat,mat))
     allocate(modified_diab_nace(mat,mat))


     !> CALCULATE MODIFIED ADIABATIC NACE 
     modified_ad_nace = 0.0D0
!COMMENT: THE DO LOOP IS UNNECESSARY, BECAUSE ONLY i=active_state IS POSSIBLE
!     do i=1,mat
!         if (i .eq. active_state) then
             do j=1,mat
                  if((j .ge. active_state) .and. (j .eq. class_forbidden_hop(j))) then
                      modified_ad_nace(active_state,j) = 0.0D0
                         do k=1,mat
                             modified_ad_nace(active_state,j) = & 
                                modified_ad_nace(active_state,j) + &
                                           interm_umatrix(dim_,t2,k,active_state)*udot_matrices(dim_,t2,k,j)
                             do l=1,mat
                                modified_ad_nace(active_state,j) = & 
                                   modified_ad_nace(active_state,j) + interm_umatrix(dim_,t2,k,active_state)* &
                                        NACE(k,l)*interm_umatrix(dim_,t2,l,j)
!print  *, "NACE(k,l)", NACE(k,l)
 
                             end do
                         end do
                  !else 
                      !modified_ad_nace(i,j) = -modified_ad_nace(j,i)
                  end if 
             end do
         !else 
          !   modified_ad_nace(i,j) = 0.0D0
!         end if
!     end do

     !> SIMMETRIZATION OF modified_ad_nace
     do i=1,mat
         do j=1,mat
              if(j .lt. i) then 
                 modified_ad_nace(i,j) = -modified_ad_nace(j,i)
              end if
         end do
     end do    
  
!print *, "modified_ad_nace", modified_ad_nace
      
     !rotated_ad_nace = modified_ad_nace

     !> ROTATE BACK modified_ad_nace INTO DIABATIC BASIS
!COMMENT: AS modified_ad_nace(k,l) = 0 IF k != active_state, A SHORTCUT FORMULA
!CAN BE USED. THE LOOP OVER k IS PROBABLY NOT USEFUL. 
   do i=1,mat
       do j=1,mat
            rotated_ad_nace(i,j) = 0.0D0
               do k=1,mat
                    do l=1,mat
                         rotated_ad_nace(i,j) =  rotated_ad_nace(i,j) + &
                   interm_umatrix(dim_,t2,i,k)*modified_ad_nace(k,l)*interm_umatrix(dim_,t2,j,l)
!print *, "rotated_ad_nace", rotated_ad_nace
                    end do
               end do
       end do
   end do

!print *, "ROTATED_ARRAY", rotated_ad_nace
     !> CALCULATION OF THE NEW MODIFIED DIABATIC NACV
!NB:THERE WAS A PROBLEM WITH THE SIGN IN THE CORRECTION I JUST CHANGE THIS SIGN
!HERE IN THE modified_diab_nace, BUT I STILL HAVE TO FIND WHERE THE EFFECTIVE
!ERROR IN SIGN IS
!     modified_diab_nace = rotated_ad_nace + NACE 
     modified_diab_nace = NACE - rotated_ad_nace

!print *, "NACE", NACE
!print *, "modified_diab_nace_NEW", modified_diab_nace
   
 
 
     deallocate(modified_ad_nace, rotated_ad_nace)

  END SUBROUTINE calc_modified_diab_nace


! WARNING:::::THIS SUBROUTINE HAS NEVER BEEN USED AND MAY NOT WORK!!!
  SUBROUTINE propagate_adiab_wf(E_t, E_tp, ctin, sh_env, nadiab,  norbitals, ELsteps, &
                           itimes, adiab_section, deltatMD, deltatE , times,      &
                            mat)
     
     REAL(KIND=dp), DIMENSION(:), INTENT(IN)         :: E_t, E_tp 
     COMPLEX(KIND=dp), dimension(:), INTENT(INOUT)   :: ctin
     INTEGER, INTENT(IN)                             :: nadiab,  norbitals 
     TYPE(sh_env_type), intent(in)                  :: sh_env
 
 
     CHARACTER(len=*), PARAMETER                     :: routineN ='propagate_adiab_wf', &
     routineP = moduleN//':'//routineN
 
     COMPLEX(KIND=dp), dimension(:), allocatable     :: ctout, ctotal
     REAL(KIND=dp), dimension(:,:), allocatable      :: NACEdiff
     REAL(KIND=dp), dimension(:,:), allocatable     ::  NACEint
     COMPLEX(KIND=dp), dimension(:,:), allocatable   :: X1, X2, X12
     REAL(KIND=dp), DIMENSION(:), allocatable        :: E_1, Ediff, E_0, E_2
 
     !COMPLEX(KIND=dp), dimension(:,:), allocatable   :: Ktot, K2, K3, K4
     !REAL(KIND=dp), dimension(:,:), allocatable      :: Imat 
 
 
     INTEGER                                         :: ELsteps   
     INTEGER, POINTER                                :: itimes
     TYPE(section_vals_type), POINTER                :: adiab_section
     REAL(KIND=dp)                                   :: deltatMD, deltatE
     REAL(KIND=dp)                                   :: ctmag
     REAL(KIND=dp)                                   :: times 
     !COMPLEX(KIND=dp)                                :: alpha, beta, ctcomp
     COMPLEX(KIND=dp)                                :: ctcomp
     integer                                         :: t2,i 
     !integer                                         :: i,j, inc, mat 
     integer                                         :: mat 
     character(len=1)                                :: zgemm_type
 
 
     !Various Allocation
 
 
     !!Imat is an identity matrix of the dimensions of the other square matrices
     !allocate(Imat(mat,mat))
     !!define once, then use whenever
     !Imat = 0.0D0
     !FORALL(j = 1:mat) Imat(j,j) = 1     
     
     allocate(NACEdiff(mat,mat))
     allocate(NACEint(mat,mat)) 
 
     allocate(ctout(mat))
     allocate(ctotal(mat))
 
     !then the complex ones i.e. H matrices; combinations of H and S
     !Note that in the solution to the Schroedinger equation, the H term is -i*H
     !For ease of calculation therefore, in this program,
     !H_t etc are always -i*H.
     !This is why they're complex variables
     !allocate(NACE(mat,mat))
     allocate(Ediff(mat))
     allocate(E_0(mat))
     allocate(E_1(mat))
     allocate(E_2(mat))
     allocate(X1(mat,mat))
     allocate(X2(mat,mat))
     allocate(X12(mat,mat))
     
    
     !!finally the K matrices, which are the terms in the Runge-Kutta integration
     !!note K1 missing, because it's just a constant times X
     !!Ktot is the final constant by which c(t) is multiplied
     !allocate(Ktot(mat,mat))
     !allocate(K2(mat,mat))
     !allocate(K3(mat,mat))
     !allocate(K4(mat,mat))
 
     !used to give correct commands to ZGEMX subroutines later
     zgemm_type = 'N'
 
 
             Ediff = (deltatE/deltatMD)*(E_tp - E_t)
             NACEdiff = (deltatE/deltatMD)*(sh_env%adiab_nace_tp - &
                                             sh_env%adiab_nace_t)
             
 
             E_0  = E_t
             E_2 = E_tp

             NACEint = sh_env%adiab_nace_t
 
             do t2=1,int(deltatMD/deltatE)
                    call make_adiabX(mat, E_0, NACEint, X1)
 
 
                    !create timestepped equivalents
                    E_1 = E_0 + 0.5D0*Ediff
                    NACEint = NACEint + 0.5D0*NACEdiff
                    
                    call make_adiabX(mat, E_1, NACEint, X12)

     
                    E_2 = E_1 + 0.5D0*Ediff
                    NACEint = NACEint + 0.5D0*NACEdiff

                    call make_adiabX(mat, E_2, NACEint, X2)
     

                    ! Call new (vector) RK4 method
                    CALL RK4(ctin, X1, X12, X2, deltatE, mat)
                    ctout = ctin(:)

                    ctmag = 0.0D0
                    ctcomp = cmplx(0.0D0, 0.0D0,dp)
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                         ctcomp = ctcomp + ctotal(i)
                    end do
                    ctmag = abs(ctcomp)
                    !rescale to conserve norm = 1
                    do i=1,mat
                         ctout(i) = cmplx((1.0D0/sqrt(ctmag))*REAL(ctout(i)), (1.0D0/sqrt(ctmag))*AIMAG(ctout(i)),dp)
                    end do
                    do i=1,mat
                         ctotal(i) = conjg(ctout(i))*ctout(i)
                    end do
 
                    !update
                    E_0 = E_2
 
                    ctin = ctout
 
                    ctout = cmplx(0.0D0, 0.0D0,dp)
 
                    times = DBLE(itimes) + (t2-1)*1/DBLE(ELsteps)

                    CALL adiab_print_coeff(adiab_section, ctin,nadiab, norbitals,itimes, deltatMD, &
                                           times)
                    CALL adiab_print_pop(adiab_section, abs(ctotal),nadiab,norbitals,itimes, deltatMD, &
                                         times)
                    CALL adiab_print_matrix_density(adiab_section, ctin, nadiab,norbitals,itimes,deltatMD, &
                                                    times)
 
             end do
 
 !END PROPAGATE THE ELECTRON FOR ONE MOLECULAR TIMESTEP
 
 !DEALLOCATE
       !DEALLOCATE(NACE)
       deallocate(E_0, Ediff, E_1, E_2, X1, X2, X12)!, Ktot, K2, K3, K4)
       !deallocate(Imat)
 !END DEALLOCATE
 
  END SUBROUTINE propagate_adiab_wf                                                                



!*********************************************************************************
!*********************************************************************************
!!!!!!!!!!!!!! DEPRECATED SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!
!*********************************************************************************
!*********************************************************************************


!**********************************************************************************
!> \brief DEPRECATED! Although this subroutine is good for testing frozen
!coupling case
!**********************************************************************************
subroutine calcforces_frozenH(Hin, C, state, Fin, forces_adiab)
       
        integer, intent(in)                         :: state
        REAL(KIND=dp), intent(in), dimension(:,:)   :: Hin, C
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), dimension(:,:), POINTER      :: forces_adiab


        CHARACTER(len=*), PARAMETER :: routineN ='calcforces_frozenH', &
                               routineP = moduleN//':'//routineN

        integer                                     :: i, j, mat, d
        INTEGER                                     :: natom
                
        natom = SIZE(Fin,2)
        mat = size(Hin,1) !might have denerate states

        !loop over 3 cartesians 
        do d=1,3 
           do j=1,natom
              forces_adiab(d,j) = 0.0D0
              ! if you wanted to be able to account for degeracy this loop
              ! should be more complicated, as the one for diagonal energies in
              ! for example AOM Hamiltonian. 
              do i=1,mat
                 ! Diagonal term, it depends on the forces. The minus is already
                 ! part of diagonal forces calculation from FIST.
                 forces_adiab(d,j) = forces_adiab(d,j) +(abs(C(i,state))**2)*Fin(d,j,i)
              end do 
           end do
        end do
        ! OFF Diagonals grandient calculation put to zero in case of frozen
        ! couplings 

end subroutine calcforces_frozenH

!**********************************************************************************
!> \brief DEPRECATED, Although good for testing Frozen coupling Hamiltonians 
!**********************************************************************************
subroutine calc_hop_nacv_frozen_H(sh_env, diagonal_F, hop_nacv, ext_old_state, ext_new_state)

        TYPE(sh_env_type)                                        :: sh_env
        REAL(kind=dp), DIMENSION(:,:), INTENT(OUT), ALLOCATABLE  :: hop_nacv
        INTEGER, OPTIONAL                                        :: ext_old_state, ext_new_state
        REAL(KIND=dp), DIMENSION(:,:,:), POINTER                 :: diagonal_F


        CHARACTER(len=*), PARAMETER :: routineN = 'calc_hop_nacv_frozen_H', &
        routineP = moduleN//':'//routineN

        !INTEGER                                      :: num_mol,  num_j
        REAL(KIND=dp), dimension(:,:), allocatable   :: C
        INTEGER                                      :: natom
        INTEGER                                      :: old_state, new_state
        REAL(KIND=dp)                                :: diff_energy
        REAL(KIND=dp)                                :: E_old_state, E_new_state
        integer                                      :: i, d, j, mat
        REAL(KIND=dp), dimension(:,:), allocatable   :: F_1!, F_11, F_2, F_22, F_3

        mat    = SIZE(sh_env%H_tp, 1)
        natom  = SIZE(diagonal_F, 2)
        !nadiab = SIZE(diagonal_F, 3)
        !atoms  = SIZE(sh_env%connlist1, 2)


        allocate(C(mat,mat))
        allocate(F_1(3,natom))
        !allocate(F_11(3,natom))
        allocate(hop_nacv(3, natom)) 

        C = sh_env%C_tp

        IF (PRESENT(ext_old_state)) THEN
           old_state = ext_old_state
        ELSE
           old_state = sh_env%active_state
        ENDIF
        IF (PRESENT(ext_new_state)) THEN
           new_state = ext_new_state 
        ELSE
           new_state = sh_env%new_state
        ENDIF
        E_old_state = sh_env%E_tp(old_state)
        E_new_state = sh_env%E_tp(new_state)
        IF (E_old_state.EQ.E_new_state) THEN
           diff_energy = 0.0D0
        ELSE
           diff_energy = 1 / (E_new_state - E_old_state)
        ENDIF  

        ! Only calculation of term a_{I, ij} in Eq. (14) of [Carof17], i.e.
        ! Eq. (15)

        ! COMMENT_SG: the commented part is just to remind myself that this
        ! subroutine could be extended to off-diagonals as well if we had the
        ! nacv for the given Hamiltonian.

        !num_mol = 0   ! Molecular index 
        do d=1,3      ! Cartesian coordinates
           !num_j = 0  ! Active atom index
           do j=1,natom ! Sum over all the atoms
              F_1(d,j) = 0.0D0
              !F_11(d,j) = 0.0D0
              !IF (ext_pres(j)) THEN
              !   num_j   = num_j + 1                  ! Update the active atom index
              !   num_mol = int((num_j-1)/atoms) + 1   ! Update the molecular index
              !   ! WARNING: we suppose that all molecules have the same
              !   ! number of atoms
              !ENDIF
              do i=1,mat
                 !F_1(d,j) = F_1(d,j) + C(i,old_state)*DH(i,i,j,d)*C(i,new_state)*diff_energy
                 ! Diagonal term, it depends on the forces
                 F_1(d,j) = F_1(d,j) + C(i,old_state)*diagonal_F(d,j,i)*C(i,new_state)*diff_energy
                 !IF (ext_pres(j)) THEN
                 !   IF (i.ne.num_mol) THEN
                 !     ! Of diagonal term, only for active atom j present in
                 !     ! molecule num_mol
                 !     F_11(d,j) = F_11(d,j) + diff_energy * convert_nacv*sh_env%nacv(d, i, num_j) * (&
                 !                             C(num_mol, new_state) * C(i, old_state) +   &
                 !                             C(num_mol, old_state) * C(i,new_state) )
                 !   ENDIF
                 !ENDIF
              ENDDO
           ENDDO
        ENDDO 

        ! minus sign to take the same convention of the paper.
        !hop_nacv = -(F_1 + F_11) 
        hop_nacv = -(F_1) ! we consider oncly diagonal elements 
 
        deallocate(F_1)
        deallocate(C)

end subroutine calc_hop_nacv_frozen_H

!**********************************************************************************
!> \brief just a test subroutine for full Gradient storage in case of tresp

!> \author SG
!**********************************************************************************
  subroutine calc_tresp_coup_test(sh_env, energies_tp, orb, NN, positions_com_tp, active_system_info, tresp_charges)

     TYPE(sh_env_type)                               :: sh_env
     REAL(kind=dp), DIMENSION(:),  POINTER           :: energies_tp
     integer, intent(in)                             :: orb
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, &
                                   INTENT(IN)        :: tresp_charges  !tresp
                                                                       !charges taken from input file, they refer 
                                                                       !to a single mol and should be extended
                                                                       !for multi-phase system
     INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)    :: NN
     REAL(kind=dp), DIMENSION(:,:),POINTER, &
                                   INTENT(IN)        :: positions_com_tp
     INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                   INTENT(IN)        :: active_system_info

     CHARACTER(len=*), PARAMETER                     :: routineN ='calc_tresp_coup', &
     routineP = moduleN//':'//routineN

     integer                                         :: i, k, m, d, dim_, l
     integer                                         :: init_atom_mol_i, final_atom_mol_i, atoms_per_mol_i
     integer                                         :: init_atom_mol_k, final_atom_mol_k, atoms_per_mol_k
     real(kind=dp), allocatable, dimension(:,:)      :: crd_i, crd_k
     real(kind=dp), allocatable, dimension(:)        :: charges1, charges2
     real(kind=dp)                                   :: tresp_coup
     real(kind=dp),  dimension(:,:), allocatable     :: grad1 !gradients acting on atoms of mol1


     !> COMMENT_SG: Tresp charegs array should be generalized for multiple molecules
     !and phases
     !> COMMENT_SG: This is the only NON general bit of the whole subroutine. It
     !assumes just one type of charges
     !#########  WARNING ##############################################
     ALLOCATE(charges1(SIZE(tresp_charges)))
     ALLOCATE(charges2(SIZE(tresp_charges)))
     charges1 = tresp_charges
     charges2 = tresp_charges
     !#################################################################


     sh_env%H_tp(:,:) = cmplx(0.0D0, 0.0D0,dp)
     sh_env%test_grad_tresp = 0.0D0

     do i=1,SIZE(energies_tp) !loop over all non-degerate states
        do k=1,SIZE(energies_tp)
           ! do calculation only for molecules within cut-off radious. 
           if (NN(i,k) == 1) then
 
              !> diagonal construction
              if (k==i) then 
                 !since you may want to include multiple orbital per state, you
                 !should loop over orbitals
                 do d=1,orb
                    sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+d) = cmplx(0.0D0,(-energies_tp(i)),dp) !minus because we propagate -X
                 end do
              else 
                 !loop over upper diagonal
                 if ( k > i ) then
                    do d=1,orb
                       do m=1,orb
                          !> find coordinates mol 1 of the pair
                          init_atom_mol_i = sh_env%active_system_info(1,i)
                          final_atom_mol_i = sh_env%active_system_info(2,i)
                          atoms_per_mol_i = (final_atom_mol_i-init_atom_mol_i)+1

                          !> find coordinates mol 2 of the pair
                          init_atom_mol_k = sh_env%active_system_info(1,k)
                          final_atom_mol_k = sh_env%active_system_info(2,k)
                          atoms_per_mol_k = (final_atom_mol_k-init_atom_mol_k)+1


                          !> this allocation deallocation is maybe inefficient
                          !but it was the only thing I found to be general for
                          !multi-phases
                          ALLOCATE(crd_i(3,atoms_per_mol_i))
                          ALLOCATE(crd_k(3,atoms_per_mol_k))

                          do dim_=1,3
                             crd_i(dim_,:) = positions_com_tp(dim_, init_atom_mol_i:final_atom_mol_i)
                             crd_k(dim_,:) = positions_com_tp(dim_, init_atom_mol_k:final_atom_mol_k)
                          end do

print *, "PAIR", i,k
                          ! calc coupling for a given pair of states
                          CALL calc_trsp_ab(crd_i, crd_k, charges1, charges2, tresp_coup)
                          sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0,-tresp_coup,dp) !sign of the off-diag is minus consisently with diagonals
                          !test to recover frozen coup for a chain:
                          !sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m)  = cmplx(0.0D0,0.0011024789,dp)

                          !> For the gradient in the case of TrESP (as it was
                          !for nacv in case of AOM) the idea is that the
                          !gradient on atom I are non-zero only if this atom
                          !belong to a molecule of the pair that are
                          !interacting. We save the gradient on atom I belonging
                          !to mol i coming from the other atoms on mol k and
                          !viceversa. See Overleaf note for a better discussion. 
                          
                          ! calc gradient on atoms belonging to mol i, from
                          ! atoms belonging to mol k, and store the former array
!print *, "state", orb*(k-1)+d, "atoms", (i-1)*atoms_per_mol_i+1, i*atoms_per_mol_i 
                          CALL calc_grad_coup( crd_i, crd_k, charges1, charges2, grad1) 


                          !> store gradient in the correct array entry
                          sh_env%test_grad_tresp(1, orb*(i-1)+d, orb*(k-1)+d,&
                                                   (i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(1,:)
                          sh_env%test_grad_tresp(2, orb*(i-1)+d, orb*(k-1)+d,&
                                                   (i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(2,:)
                          sh_env%test_grad_tresp(3, orb*(i-1)+d, orb*(k-1)+d, &
                                                   (i-1)*atoms_per_mol_i+1:i*atoms_per_mol_i) = grad1(3,:)

                          ! calc gradient on atoms belonging to mol k, from
                          ! atoms belonging to mol i, and store the former arrays
!print *, "state", orb*(i-1)+m, "atoms", (k-1)*atoms_per_mol_k+1, k*atoms_per_mol_k 
                          CALL calc_grad_coup( crd_k, crd_i, charges2, charges1, grad1) 
                          
                      
                          !> store gradient in the correct array entry
                          sh_env%test_grad_tresp(1, orb*(i-1)+m, orb*(k-1)+d, &
                                               (k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(1,:)                        
                          sh_env%test_grad_tresp(2, orb*(i-1)+m, orb*(k-1)+d, &
                                               (k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(2,:)                        
                          sh_env%test_grad_tresp(3, orb*(i-1)+m, orb*(k-1)+d, &
                                               (k-1)*atoms_per_mol_k+1:k*atoms_per_mol_k) = grad1(3,:)                        
 
                          ! I should find a way to store grad1 in a more general
                          ! global_grad array that in principle should have
                          ! dimension (3, mat, mat, natoms_TOTAL_active) or in a reduce
                          ! form assuming all atoms_per_mol the same (3, mat, mat*natoms_per_mol))
  
                          DEALLOCATE(crd_i,crd_k)
                      end do
                   end do
                 else
                    do d=1,orb 
                       do m=1,orb 
                          sh_env%H_tp(orb*(i-1)+d,orb*(k-1)+m) = sh_env%H_tp(orb*(k-1)+m,orb*(i-1)+d)
                          sh_env%test_grad_tresp(1,orb*(i-1)+d,orb*(k-1)+m,:) = &
                                         sh_env%test_grad_tresp(1,orb*(k-1)+d,orb*(i-1)+m,:)
                          sh_env%test_grad_tresp(2,orb*(i-1)+d,orb*(k-1)+m,:) = &
                                         sh_env%test_grad_tresp(2,orb*(k-1)+d,orb*(i-1)+m,:)
                          sh_env%test_grad_tresp(3,orb*(i-1)+d,orb*(k-1)+m,:) = &
                                         sh_env%test_grad_tresp(3,orb*(k-1)+d,orb*(i-1)+m,:)
                       end do 
                    end do
                 end if 
              end if 


           end if
        end do 
     end do

!!TEST
OPEN (UNIT = 21, FILE = "GRAD_formatted_FULL.dat")

do i=1,SIZE(sh_env%test_grad_tresp,4) !atoms
  do k=1, SIZE(sh_env%test_grad_tresp,2)   !mol
     do l=1, SIZE(sh_env%test_grad_tresp,2)   !mol
     write(21,*) "atom", i, "state", k, l, sh_env%test_grad_tresp(:,k,l,i)
     end do
  end do
end do

close (21)
  end subroutine  calc_tresp_coup_test



subroutine calcforces_smart_grad_test(state, C, Fin, forces_adiab, & 
                  ext_pres,  active_system_info, orb, full_grad, convert_nacv, frozen_c)

        !per our full force expression derivation, this calculates the adiabatic forces
        !on each atom
        integer, intent(in)                         :: state, orb
        REAL(KIND=dp), intent(in), dimension(:,:)   :: C
        REAL(KIND=dp), intent(in), dimension(:,:,:) :: Fin
        REAL(KIND=dp), dimension(:,:), POINTER      :: forces_adiab
        LOGICAL, dimension(:), INTENT(IN)           :: ext_pres
        INTEGER, DIMENSION(:,:), ALLOCATABLE, &
                                      INTENT(IN)     :: active_system_info
        REAL(KIND=dp), intent(in), OPTIONAL          :: convert_nacv
        LOGICAL, OPTIONAL, INTENT(IN)                :: frozen_c
        REAL(KIND=dp), intent(in), dimension(:,:,:,:), &
                                      OPTIONAL       :: full_grad

        CHARACTER(len=*), PARAMETER :: routineN ='calcforces_smart_grad', &
                               routineP = moduleN//':'//routineN

        integer                                     :: i, j,  d, l
        INTEGER                                     :: natom
        INTEGER                                     :: num_mol, num_j
        REAL(KIND=dp)                               :: my_convert_nacv 
        LOGICAL                                     :: my_frozen_c

        natom = SIZE(Fin,2)

        !COMMENT_SG: This subroutine should be checked further for degeneracy.
        !Note that Fin(d,j,i) is the array of diagonal forces and its dimension
                            !stops at the number of total force_evals
                            !considered. 
        !The following subroutine may NOT work if orb !=1 (degeneracy). It needs
        !testing.

        ! Based on paper [Spencer16]: Spencer et al., J. Chem. Phys. 145, 064102 (2016)
        ! Use Eq. (21) of paper [Spencer16]

        my_frozen_c = .FALSE.
        IF (PRESENT(frozen_c)) my_frozen_c = frozen_c


        do d=1,3                        ! sum of cartesian 
           num_j = 0               ! active atom index
           do j=1,natom            ! sum over all atoms
              forces_adiab(d,j) = 0.0D0

              IF (ext_pres(j)) THEN   
                 num_j   = num_j + 1                 ! update active atom index
                 ! Get the molecule to which j belong, for a general
                 ! system. This is a molecular index needed below. Looping
                 ! over molecules is necessary here.
                 do i=1,SIZE(active_system_info, 2) ! loop over non-deg states
                    IF ((j .GE.  active_system_info(1,i)) .AND. (j .LE. active_system_info(2,i))) THEN
                       num_mol = i
                       EXIT 
                    END IF
                 end do
              ENDIF

              do i=1,SIZE(active_system_info, 2)
                   ! Diagonal term, it depends on the forces
                   ! I assume that diagonal forces
                   ! for two degenrate
                   ! states are the same as for site energies
                   forces_adiab(d,j) = forces_adiab(d,j) + (abs(C(i, state))**2)*Fin(d,j,i)

                   
                 IF (ext_pres(j)) THEN
                   do l=1, SIZE(active_system_info, 2)
                      my_convert_nacv = 1.0D0
                      if (l .ne. i) THEN 
                         forces_adiab(d,j) = forces_adiab(d,j) + &
                            (C(i, state)*&
                              full_grad(d, i,l, num_j)*C(l,state))
                      end if
                   end do 
                END IF
              end do 
                      
 
            !       ! calculate off diagonal forces if
            !       ! couplings are not frozen
            !       IF (.NOT. my_frozen_c) THEN

            !          ! Of diagonal term, only for active atom j present in
            !          ! molecule num_mol
            !          IF (ext_pres(j)) THEN
            !             IF (i.ne.num_mol) THEN
            !                 IF (PRESENT(convert_nacv)) THEN
            !                    !> Forces in case of AOM formalism
            !                    my_convert_nacv = convert_nacv ! to convert the FOB NACV
            !                 ELSE 
            !                    !> Forces in case of general gradients (smartly stored as
            !                    !> nacv in AOM case)
            !                    my_convert_nacv = 1.0D0 !gradient are already a negative object so no minus sign here
            !                 END IF

            !                 forces_adiab(d,j) = forces_adiab(d,j) + &
            !                       ((C(orb*(num_mol-1)+m,orb*(state-1)+m))*&
            !                         nacv(d, orb*(i-1)+m, num_j)*C(orb*(i-1)+m,orb*(state-1)+m))*my_convert_nacv*2.0D0
            !             ENDIF
            !          ENDIF
            !       END IF
            !     end do
            !  end do 
           ENDDO
        ENDDO

end subroutine calcforces_smart_grad_test


!!!!!!! DEPRECATED !!!!!
! Could be removed
SUBROUTINE adiab_overalp_calc(active_state, sh_env, C_tp, mol, itimes, & 
                                   attempt_hop, enforced_active_state)

     INTEGER, intent(INOUT)                         :: active_state, mol, itimes
     INTEGER, intent(OUT)                         :: enforced_active_state 
     LOGICAL, INTENT(INOUT)                           :: attempt_hop
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: C_tp
     TYPE(sh_env_type), INTENT(IN)              :: sh_env

     CHARACTER(len=*), PARAMETER :: routineN = 'adiab_overalp_calc', &
              routineP = moduleN//':'//routineN
     REAL(KIND=dp), dimension(:), allocatable    :: overlaps
     REAL(KIND=dp)                                :: sum1!, sum
     INTEGER                                     :: i, loc(1) 

     allocate(overlaps(4))

     overlaps = 0.0D0

!     sum = 0.0D0
     sum1 = 0.0D0
     
     do i=1, mol
!         sum = sum + (abs(sh_env%C_t(i, active_state)) *abs(C_tp(i, active_state)))
         sum1 = sum1 + sh_env%C_t(i, active_state)*C_tp(i, active_state)
     end do
     
     
!     print *, "OVERLAP_TRIVIAL", sum
     print *, "PHASE_OVERALAP", abs(sum1)
     
     if (abs(sum1) .lt. 0.7 .and. .not. attempt_hop) then
        print *, "POSSIBLE_TRIVIAL_CROSSING"
        do i =1,mol
            if (active_state .gt. 2) then
                overlaps(1) = overlaps(1) + sh_env%C_t(i, active_state) *C_tp(i,active_state+2)
                overlaps(2) = overlaps(2) + sh_env%C_t(i, active_state) *C_tp(i,active_state+1)
                overlaps(3) = overlaps(3) + sh_env%C_t(i, active_state) *C_tp(i,active_state-1)
                overlaps(4) = overlaps(4) + sh_env%C_t(i, active_state) *C_tp(i,active_state-2)
            else if (active_state .gt. 1) then
                overlaps(1) = overlaps(1) + sh_env%C_t(i, active_state)*C_tp(i,active_state+2)
                overlaps(2) = overlaps(2) + sh_env%C_t(i, active_state)*C_tp(i,active_state+1)
                overlaps(3) = overlaps(3) + sh_env%C_t(i, active_state)*C_tp(i,active_state-1)
            else 
                overlaps(1) = overlaps(1) + sh_env%C_t(i,active_state)*C_tp(i,active_state+2)
                overlaps(2) = overlaps(2) + sh_env%C_t(i,active_state)*C_tp(i,active_state+1)
            end if
     
       end do
       print *, "OVERALPS", abs(overlaps(:))
       print *, "MAX_OVERLAP", itimes, abs(sum1), maxval(abs(overlaps(:)))
print *, "MORE_INFO", itimes, sum1, overlaps(1), overlaps(2), overlaps(3), overlaps(4)
!       print *, "MAXLOC", maxloc(abs(overlaps(:)))
       loc = maxloc(abs(overlaps(:)))
       if (maxval(abs(overlaps(:))) .lt. sum1) then 
           print *, "WARNING: MORE THAN TWO STATES CROSS"
           enforced_active_state = active_state
       else     
           if (loc(1) == 1) then
               enforced_active_state = active_state+2
           else if (loc(1) == 2) then
               enforced_active_state = active_state+1
           else if (loc(1) == 3) then
               enforced_active_state = active_state-1
           else if (loc(1) == 4) then
               enforced_active_state = active_state-2
           end if
       print *, "enforced_active_state", enforced_active_state 
print *, "STATE_TO_GO", active_state, enforced_active_state 
print *, "LOCALIZATION", maxloc(abs(sh_env%C_t(:, active_state)*sh_env%C_t(:, active_state))), &
        maxloc(abs(C_tp(:,enforced_active_state)*C_tp(:,enforced_active_state)))
       end if
      !enforced_active_state = 0
         
     end if
     
     deallocate(overlaps)


END SUBROUTINE adiab_overalp_calc

END MODULE non_adiabatic_utils 
                       
