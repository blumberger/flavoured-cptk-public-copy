!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the MIXED environment
!> \par History
!>      10.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE input_cp2k_mixed
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             low_print_level,&
                                             medium_print_level
  USE input_constants,                 ONLY: mix_coupled,&
                                             mix_generic,&
                                             mix_linear_combination,&
                                             mix_minimum,&
!                                            mix_restrained
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                                             mix_restrained, &
                                             mix_adiab
!> End CHANGE_AC

  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mixed'

  PUBLIC :: create_mix_section

CONTAINS

! *****************************************************************************
!> \brief Create the input section for MIXED.
!> \param section the section to create
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mix_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub2section, sub3section, &
                                                subsection
    TYPE(section_type), POINTER              :: print_key


    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="MIXED",&
         description="This section contains all information to run with a hamiltonian "//&
                     "defined by a mixing of force_evals",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    NULLIFY(keyword, subsection)

    CALL keyword_create(keyword, name="MIXING_TYPE",&
         description="The type of mixing to be employed",&
         usage="MIXING_TYPE LINEAR_COMBINATION",&
         default_i_val=mix_linear_combination,&
         enum_c_vals=s2a("LINEAR_COMBINATION", &
                         "MINIMUM",&
                         "COUPLED",&
                         "RESTRAINT",&
!                        "GENMIX"),&
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                        "GENMIX",&
                        "ADIABATIC"),&
!> End CHANGE_AC
         enum_desc=s2a("Linear combination of force envs (support only 2 force_evals)", &
                       "Use the force env with the minimum energy (support only 2 force_evals)",&
                       "Consider the force envs as a two state system with a given"//&
                       " coupling matrix element (support only 2 force_evals)",&
                       "Use the difference between the energy of the force envs as a"//&
                       " restraint on the first (support only 2 force_evals)",&
!                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)"),&
!         enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic/))
!> CHANGE_AC
!> 15-12-11
!> To create a new way to combine the force
!> To get the adiabatic force 
                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)",&
                 "Combines the diabatic forces to get the adiabatic forces (support for an unlimited number of force_eval)"),&
         enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic, &
                       mix_adiab /))
!> End CHANGE_AC
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="GROUP_PARTITION",&
         description="gives the exact number of processors for each group."//&
         " If not specified processors allocated will be equally distributed for"//&
         " the specified subforce_eval, trying to build a number of groups equal to the"//&
         " number of subforce_eval specified.",&
         usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
         description="Gives the wanted number of groups. If not specified the number"//&
         " of groups is set to the number of subforce_eval defined.",&
         usage="ngroups 4", type_of_var=integer_t)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Double force_eval
    CALL section_create(subsection,name="LINEAR",&
         description="Linear combination between two force_eval:  F= lambda F1 + (1-lambda) F2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="LAMBDA",&
         description="Specify the mixing parameter lambda in the formula.",&
         usage="lambda <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="COUPLING",&
         description="Coupling between two force_eval: E=(E1+E2 - sqrt((E1-E2)**2+4*H12**2))/2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="COUPLING_PARAMETER",&
         description="Coupling parameter H12 used in the coupling",&
         usage="COUPLING_PARAMETER <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="RESTRAINT",&
         description="Restraint between two force_eval: E = E1 + k*(E1-E2-t)**2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="RESTRAINT_TARGET",&
         description="Target value of the restraint (t) ",&
         usage="RESTRAINT_TARGET <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RESTRAINT_STRENGTH",&
         description="Strength of the restraint (k) in "//&
         "k*(E1-E2-t)**2" ,&
         usage="RESTRAINT_STRENGTH <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Adiabatic force eval
!> CHANGE_AC
!> 15-12-11
!> Create the input part for the Non Adiabatic Molecular Dynamics

    !> main ADIABATIC section part of the general MIXED section
    CALL section_create(subsection,name="ADIABATIC",&
         description="Defines the methods to obtain the adiabatic forces.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    !-----------------------------------------------------------------------
    !> start new sub2section OPTIMIZATION, part of ADIABATIC subsection
    NULLIFY(sub2section)
    CALL section_create(sub2section,name="OPTIMIZATION",&
         description="Defines the optimization methods to calculate site,"//&
                      "energies from FIST global energies and forces",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL section_add_subsection(subsection,sub2section)

    CALL keyword_create(keyword, name="DO_SPEEDUP_LJ",&
         description="If DO_SPEEDUP_LJ is TRUE, the first force_eval section"//&
                     " will only calculate LJ interactions (DEPRECATED KEYWORD)",&
         usage="True", default_l_val=.True., &
         lone_keyword_l_val=.True.) 
        ! type_of_var=char_t,&
        ! n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DO_SPEEDUP_INTRA",&
         description="If DO_SPEEDUP_INTRA is TRUE, the first force_eval section will only calculate LJ interactions,"//&
                      " the pair force eval section will calculate intramolecular"//&
                      " interaction for one charged molecule only and the impair for one neutral"//&
                      " molecule only. WARNING: DO_SPEEDUP_LJ should be put to TRUE. ((DEPRECATED KEYWORD)",&
         usage="True", default_l_val=.True., &
         lone_keyword_l_val=.True.) 
        ! type_of_var=char_t,&
        ! n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DO_SPEEDUP_NACE",&
    !     description="If set to TRUE will apply nace optimisations.",&
    !     usage="True", default_l_val=.TRUE., &
    !     lone_keyword_l_val=.True.) 
    !    ! type_of_var=char_t,&
    !    ! n_var=-1)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SURF_HOP_CHOICE",&
         description="Decides whether to use or not the trivial hopping correction (see e.g. Carof17)",&
         usage="SURF_HOP_CHOICE: TRIVIAL_CORR_WITH_OVERLAP or TRIVIAL_HOP_CORRECT or UNMODIFIED_SURF_HOP", &
         type_of_var=char_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RK_PROPAGATION",&
         description="Decides which RK scheme is used and the approximation level (see Giannini20)",&
         usage="DIABATIC_RK or ADIABATIC_RK (deprecated) or HS_CORRECTED_RK"//&
         " or DIABATIC_RK_FAST or DIABATIC_RK_DLI", & 
         type_of_var=char_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="REORDERING_STATES_USING_OVERLAP",&
         description="Choose if you want to reorder the states with a tracking alghoritm"//&
                     " (see e.g. Giannini18, Carof19)",&
         usage="False", default_l_val=.False., &
         lone_keyword_l_val=.True.) 
        ! type_of_var=char_t,&
        ! n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !> releasing OPTIMIZATION subsection
    CALL section_release(sub2section)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !> start new sub2section DECOHERENCE, part of ADIABATIC subsection
   NULLIFY(sub2section)
   CALL section_create(sub2section,name="DECOHERENCE",&
        description="Defines the decoherence methods,"//&
                     " and related parameters to use (see e.g. Giannini18, Carof19, Carof17).",&
        n_keywords=1, n_subsections=0, repeats=.FALSE.)
   CALL section_add_subsection(subsection,sub2section)

    CALL keyword_create(keyword, name="DECOHERENCE_CORRECTIONS",&
         description="Defines the methods to decohere the wavefunction if required.",&
         usage="DECOHERENCE_CORRECTIONS DAMPING", type_of_var=char_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)
    
    CALL keyword_create(keyword, name="DECO_TIME",&
         description="Defines the decoherence time used by the damping (see e.g. Carof19)",&
         usage="DECO_TIME EDC or FORCES_BASED or HYBRID_TAU", type_of_var=char_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SPURIOUS_TRANSFER_CORR",&
         description="Choose if you want to activate the spurious tranfer correction (see e.g. Giannini18).",&
         usage="SPURIOUS_TRANSFER_CORR: True", default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.) 
        ! type_of_var=char_t,&
        ! n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EDC_C",&
         description="C parameter in the Energy-based decoherence correction formula (see e.g. Carof 19)",&
         usage="1.0D0", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EDC_E0",&
         description="E0 parameter in the Energy-based deco. correction formula (see e.g. Carof 19)",&
         usage="0.1D0", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="THRESHOLD_TAU_FORCES",&
         description="Threshold used in the force based calculation to avoid"//&
                     " calculating forces for states high in energy (pure dephasing is used instead).",&
         usage="0.1D0", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMPERATURE_FG_WIDTH",&
         description="Temperature to define FG width in the force-based formula (see Carof 19)",&
         usage="300", default_r_val=300.0_dp, type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DECOHERENCE_CRITERION",&
         description="Threshold criteria related to instantaneous decoherence (see Carof 19)",&
         usage="1.0E-05", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !> releasing DECOHERENCE subsection
    CALL section_release(sub2section)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !> start new sub2section ENERGY_CONSERVATION, part of ADIABATIC subsection
   NULLIFY(sub2section)
   CALL section_create(sub2section,name="ENERGY_CONSERVATION",&
        description="Defines the velocity adj,"//&
                     "and related parameters to use (see Carof 17)",&
        n_keywords=1, n_subsections=0, repeats=.FALSE.)
   CALL section_add_subsection(subsection,sub2section)

   CALL keyword_create(keyword, name="METHOD_RESCALING",&
        description="Defines the methods to rescale the velocities if required.",&
        usage="METHOD_RESCALING SIMPLE or SIMPLE_QSYS or NACV", type_of_var=char_t)
   CALL section_add_keyword(sub2section,keyword)
   CALL keyword_release(keyword)

   CALL keyword_create(keyword, name="METHOD_ADIABATIC_NACV",&
        description="Defines the methods to compute the adiabatic NACV for the rescaling.",&
        usage="METHOD_RESCALING TOTAL or FAST or (others ver used)", type_of_var=char_t)
   CALL section_add_keyword(sub2section,keyword)
   CALL keyword_release(keyword)

   CALL keyword_create(keyword, name="METHOD_REVERSAL",&
        description="Defines the methods to reverse the velocities after a frustrated hop if required.",&
        usage="METHOD_RESCALING NEVER or ALWAYS or TRUHLAR or SUBOTNIK", type_of_var=char_t)
   CALL section_add_keyword(sub2section,keyword)
   CALL keyword_release(keyword)

   !> releasing ENERGY_CONSERVATION subsection
   CALL section_release(sub2section)
    !-----------------------------------------------------------------------

    CALL keyword_create(keyword, name="METHOD_PROPAGATION",&
         description="Defines the methods to obtain the adiabatic forces usually"//&
                     " FSSH or BORN_OPPENHEIMER (some other"//&
                     " options are not used very often and not have been tested in a while, but they are"//& 
                     " usefull for testing purpose when testing e.g. Rabi Oscillations.",&
         usage="METHOD_PROPAGATION : CLASSICAL_PATH, FSSH, BORN_OPPENHEIMER",  type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
 
    CALL keyword_create(keyword, name="PRINT_MORE",&
         description="Choose if you want to print more information in FSSH files",&
         usage="PRINT_MORE: True", default_l_val=.True., &
         lone_keyword_l_val=.TRUE.) 
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INITIALIZATION",&
         description="Decide how to initiate the system for wavefunction and surface populations.",&
         usage="DIABATIC or ADIABATIC", & 
         type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="REPRESENTATION",&
         description="Decide which representation is used: Diab  or Adiab basis (the latter is deprecated)",&
         usage="DIABATIC_BASIS or ADIABATIC_BASIS", & 
         type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ELECTRONIC_PARTIAL_STEP",&
         description="Number of electronic steps per nuclear step in the RK integration",&
         usage="5", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="FIRST_DIABAT",&
         description="Index of the molecule on which to start the wavefunction"//&
                     " in case of DIABATIC initialization (corresponding adiabatic state is"//&
                     " chosen as described e.g. in Carof 19)",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="FIRST_ADIABAT",&
         description="Adiabatic state from which the wavefunction is started in"//&
                     " case of ADIABATIC initialization (corresponding diabatic wavefunction is"//&
                     " retrived by projection on the diabatic basis, see Carof 19) ",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword,name="SEED",&
         description="Initial seed for the global (pseudo)random number "//&
         " generator tocreate a stream of normally Gaussian "//&
         " distributed random numbers.",&
         usage="SEED <INTEGER>",default_i_val=2000)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_DIABATIC_STATES",&
         description="Number of (non degenerate) diabatic states to consider"//&
                     " when constructing the Hamiltonian of the system",&
         usage="2", type_of_var=integer_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_ORBITALS",&
         description="Number of orbitals per sites (e.g. for degenerate states - UNTESTED for !=1 )",&
         usage="1", type_of_var=integer_t, default_i_val=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="NUCLEAR_TIMESTEP",&
         description="Timestep for atomic motion, must be the same as in MD"//&
                     " section (keyword to remove in future becuase replicated)",&
         usage="timestep 0.5", type_of_var=real_t, &
         unit_str = "fs")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CENTER_OF_MASS",&
         description="Decide if do all electronic propagation and velocities "//&
                     "rescaling into the center-of-mass frame.",&
         usage="True", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INDEX_ACTIVE_MOLS",&
         description='Specify the index of mols that are active if DECOMP not present', &
         usage="INDEX_MOL_DECOMP <INTEGER> .. <INTEGER>", &
         type_of_var=integer_t,  n_var=-1, default_i_vals=(/1/))
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    !-----------------------------------------------------------------------
    !> start new sub2section METHOD_SITE_ENERGIES part ADIABATIC subsection
    
   NULLIFY(sub2section,sub3section)
   CALL section_create(sub2section,name="METHOD_SITE_ENERGIES",&
        description="Defines the combination of site energies"//&
                     " for the construction of the Hamiltonian and implement the"//&
                     " addition/subtraction scheme for various blocks of H. ",&
        n_keywords=1, n_subsections=0, repeats=.FALSE.)
   CALL section_add_subsection(subsection,sub2section)

    CALL keyword_create(keyword, name="MULTI_STATES",&
         description="Tell the H and forces calculation if there are more non-degenerate states "//&
                     "than active molecules (e.g. CT block has N_DxN_A states)",&
         usage="false", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.FALSE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="BLOCK_TYPE",&
         description="Determine the right site energy combination for the"//&
                     " generation of a given Hamiltonian (e.g. FE or CT-FE ..)",&
         usage="BLOCK_TYPE CT or FE or CT-FE or CT-FE-CR or CT-CR",  type_of_var=char_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

   CALL section_create(sub3section,name="CT_BLOCK",&
        description="Defines the CT BLOCK,"//&
                     "and its interactions",&
        n_keywords=1, n_subsections=0, repeats=.FALSE.)
   CALL section_add_subsection(sub2section,sub3section)

    CALL keyword_create(keyword, name="ADD_COULOMB_POTENTIAL",&
         description="add screened Coulomb potential in case of CT "//&
                     "states",&
         usage="false", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

!----------------------------------------------------------------------------
    !two keywords added in ct-block section of SITE ENERGIES M,ETHODS
    !SECTION, onSITE NUMBER specifies number of CT states in system, and
    !XT_CT_OFFSET specifies energy offset between interfacial CT-state
    !and XT manifold -- FI 05/01/24

    CALL keyword_create(keyword, name="SITE_NUMBER",&
         description="Number of diabatic states in your full Hamiltonian",&
         usage="positive integer", type_of_var=integer_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="XT_CT_OFFSET",&
         description="energy offset between interfacial CT diabat and XT diabats",&
         usage="real value", type_of_var=real_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

!------------------------------------------------------------------------------

    CALL keyword_create(keyword, name="DIELECTRIC_SCREENING",&
         description="Dielectric epsilon of the screened Coulomb",&
         usage="3.0 WITHOUT UNITS", type_of_var=real_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_DONOR_ATOMS",&
         description="Number of atoms per site",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_ACCEPTOR_ATOMS",&
         description="Number of atoms per site",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="The ground state ESP charges in the format:"//&
         "<p><tt>ATOMIC_KIND  LIST NO TRESP</tt></p>",&
         repeats=.TRUE., usage="{{String} {Int} {Int} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INDEX_DONOR_MOLS",&
         description='Specify the index of donor mols (xyz or topology order)', &
         usage="INDEX_DONOR_MOLS <INTEGER> .. <INTEGER>", &
         type_of_var=integer_t,  n_var=-1, default_i_vals=(/1/))
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INDEX_ACCEPTOR_MOLS",&
         description='Specify the index of acceptor mols (xyz or topology order)', &
         usage="INDEX_ACCEPTOR_MOLS <INTEGER> .. <INTEGER>", &
         type_of_var=integer_t,  n_var=-1, default_i_vals=(/1/))
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

   CALL section_release(sub3section)

   NULLIFY(sub3section)
   CALL section_create(sub3section,name="FE_BLOCK",&
        description="Defines the FE Block,"//&
                     "and the interactions",&
        n_keywords=1, n_subsections=0, repeats=.FALSE.)
   CALL section_add_subsection(sub2section,sub3section)

    CALL keyword_create(keyword, name="INDEX_EXC_MOLS",&
         description='Specify the index of excitonic mols (xyz or topology order)', &
         usage="INDEX_DONOR_MOLS <INTEGER> .. <INTEGER>", &
         type_of_var=integer_t,  n_var=-1, default_i_vals=(/1/))
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_release(sub3section)

    NULLIFY(sub3section)
    CALL section_create(sub3section,name="FE_CT_BLOCK",&
         description="Defines the FE_CT Block,"//&
                     "and the interactions",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL section_add_subsection(sub2section,sub3section)

    CALL keyword_create(keyword, name="FE_IS",&
         description="Defines the FE phase as donor or acceptor.",&
         usage="FE_IS DONOR or ACCEPTOR", type_of_var=char_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

!WTP_4.2
    CALL keyword_create(keyword, name="E_OFFSET",&
         description="E_FE(R_gs) - E_CT(R_gs)",&
         usage="E_OFFSET <REAL>", type_of_var=char_t)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

!WTP_comment: maybe not useful anymore?
!    CALL keyword_create(keyword, name="N_CT_ATOMS",&
!         description="Number of atoms per site",&
!         usage="1", type_of_var=integer_t)
!    CALL section_add_keyword(sub3section,keyword)
!    CALL keyword_release(keyword)

!    CALL keyword_create(keyword, name="N_FE_ATOMS",&
!         description="Number of atoms per site",&
!         usage="1", type_of_var=integer_t)
!    CALL section_add_keyword(sub3section,keyword)
!    CALL keyword_release(keyword)

   CALL section_release(sub3section)


   !> releasing METHOD_SITE_ENERGIES subsection
   CALL section_release(sub2section)

!********************************************************************************

    !> start new sub2section METHOD_COUPLING, part of ADIABATIC subsection
    CALL keyword_create(keyword, name="METHOD_COUPLING",&
         description="Defines the method to calculate off-diagonal matrix elements.",&
         usage="METHOD_COUPLING AOM or TRESP or FROZEN_COUPLINGS",  type_of_var=char_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
!********************************************************************************

    !-----------------------------------------------------------------------
    !> start new sub3section FROZEN_COUPLINGS part of METHOD_COUPLING sub2section
    NULLIFY(sub2section) 
    CALL section_create(sub2section,name="FROZEN_COUPLINGS",&
         description="Couplings are passed as input and kept frozen during dynamics"//&
                     " (only site energies fluctuate).",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL section_add_subsection(subsection,sub2section)

    CALL keyword_create(keyword, name="CONNECTIVITY_FRZ_C",&
         description="Decide if we use the connectivity from input",&
         usage="True", & 
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="INTERACTION_NUMBER",&
         description="Number of interactions included in the Hamiltonian",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="The FRZ COUPLINGS from connectivity in the format:"//&
         "<p><tt>ATOMIC_KIND  LIST NO TRESP</tt></p>",&
         repeats=.TRUE., usage="{{Int} {Int} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="FIRST_OFF_DIAGONALS_VAL",&
         description="Hab for the first upper and lower diagonals or the Hamitlonian (tridiagonal matrix)",&
         usage="001102479, HARTREE", type_of_var=real_t,&
         default_r_val=0.000_dp,    &
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !> releasing FROZEN_COUPLINGS subsection
    CALL section_release(sub2section)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !> start new sub3section TRESP_COUPLINGS part of METHOD_COUPLING sub2section
    NULLIFY(sub2section) 
    CALL section_create(sub2section,name="TRESP_COUPLINGS",&
         description="Uses transition ESP charges to calculate excitonic couplings.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL section_add_subsection(subsection,sub2section)


    CALL keyword_create(keyword, name="NUMBER_AOM_ATOMS_PER_SITE",&
         description="Number of atoms per site",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TRESP_FROM_GAUSSIAN",&
         description="Determines if TRESP are taken from a calculation with Gaussian software."//& 
                     " If so, charges must be divided by sqrt(2). ",&
         usage="True", & 
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_BETWEEN_SITES",&
         description="Cut-off distance over which excitonic couplings are set to zero.",&
         usage="20.0D0", type_of_var=real_t, &
         unit_str="angstrom")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_BETWEEN_FE_SITES",&
         description="Cut-off distance over which excitonic couplings are set to zero.",&
         usage="20.0D0", type_of_var=real_t, &
         unit_str="angstrom")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="The TRESP charges in the format:"//&
         "<p><tt>ATOMIC_KIND  LIST NO TRESP</tt></p>",&
         repeats=.TRUE., usage="{{String} {Int} {Int} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !> releasing FROZEN_COUPLINGS subsection
    CALL section_release(sub2section)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !> start new sub3section AOM part of METHOD_COUPLING sub2section
    NULLIFY(sub2section)
    CALL section_create(sub2section,name="AOM",&
         description="Uses AOM to calculate couplings.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL section_add_subsection(subsection,sub2section)

    CALL keyword_create(keyword, name="NUMBER_AOM_ATOMS_PER_SITE",&
         description="Number of atoms per site",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_BETWEEN_SITES",&
         description="Cut-off distance over which electronic couplings are set to zero.",&
         usage="6.0D0", type_of_var=real_t, &
         unit_str="angstrom")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CUT_OFF_CONNECTIVITY",&
         description="Distance cut-off for H_ab calculation",&
         usage="3.50D0", type_of_var=real_t, &
         unit_str="bohr")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SCALING_FACTOR",&
         description="Hab = C*Sab. This is the scaling C parametrized prior simulation (see Gajdos 14)",&
         usage="0.50D0*0.065190D0, HARTREE", type_of_var=real_t,&
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SCALING_FACTOR_DONOR",&
         description="Hab = C*Sab. This is the scaling C for the donor phase",&
         usage="0.50D0*0.065190D0, HARTREE", type_of_var=real_t,&
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SCALING_FACTOR_ACCEPTOR",&
         description="Hab = C*Sab. This is the scaling C for the acceptor phase",&
         usage="0.50D0*0.065190D0, HARTREE", type_of_var=real_t,&
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="SCALING_FACTOR_FECT",&
         description="Hab = C*Sab. This is the scaling C for the FE-CT coupling",&
         usage="0.50D0*0.065190D0, HARTREE", type_of_var=real_t,&
         unit_str="hartree")
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CBAR",&
         description="Sab constant factor (Not used for fast porpagation)",&
         usage="0.5082 WITHOUT UNITS", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ATOMIC_OVERLAP_CUTOFF",&
         description="An atomic overlap cutoff",&
         usage="1.0D-17", type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)


    !CALL keyword_create(keyword, name="COLLAPSE",&
    !     description="Decide if we use the decoherence procedure",&
    !     usage="True", & 
    !     default_l_val=.FALSE., &
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RANDOM_DERIVATIVE",&
         description="Decide if we use forward and backward derivative "//&
                      "in a random way (True) or only forward (False). (NOTE: keep it True).",&
         usage="True", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ANALYTICS",&
         description="Decide if we use the analytical equation rather than "//&
                      "the numerical expression. Only for 2x2 systems and useful for testing purpose only.",&
         usage="False", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.FALSE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MULTIPLE_TIME_STEP",&
         description="Optimization for code speed up. Number of steps to skip"//&
                      " before the next gradient calculation. See Giannini 20 for details.",&
         usage="1", type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)


    CALL keyword_create(keyword, name="NACV_INCREMENT",&
         description="Spatial Increment to calculate the NACV with Finite Difference. (see Spencer 16).",&
         usage="1E-3*1.889725989D0 (bohr)", & 
         default_r_val=0.00188972598_dp,    &
         unit_str="bohr",                   &
         type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="PPI_FILE_NAME",&
    !     description="Specifies the filename that contains the lookup table",&
    !     usage="PPI_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="PSIGMA_FILE_NAME",&
    !     description="Specifies the filename that contains the lookup table",&
    !     usage="PPI_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    ! OZ: enable C-S and S-S lookup tables keywords

    !CALL keyword_create(keyword, name="PPI_C_S_FILE_NAME",&
    !     description="Pi-pi C-S AOM lookup table",&
    !     usage="PPI_C_S_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="PSIGMA_C_S_FILE_NAME",&
    !     description="Sigma-pi C-S lookup table",&
    !     usage="PPI_C_S_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="PPI_S_S_FILE_NAME",&
    !     description="Pi-pi S-S AOM lookup table",&
    !     usage="PPI_S_S_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="PSIGMA_S_S_FILE_NAME",&
    !     description="Sigma-pi S-S lookup table",&
    !     usage="PPI_S_S_FILE_NAME <FILENAME>",type_of_var=lchar_t)
    !CALL section_add_keyword(sub2section,keyword)
    !CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="The HOMO coefficient in the format:"//&
         "<p><tt>ATOMIC_KIND  LIST NO HOMOS HOMOP</tt></p>",&
         repeats=.TRUE., usage="{{String} {Int} {Int} {Real} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STO_EXPONENTS",&
         description="Slater-type orbital mu exponential coeffs",&
         usage="STO_EXPONENTS <Z:value,Z:value,...>",type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STO_RMIN",&
         description="Spline rmin",&
         usage="STO_RMIN <value>", &
         default_r_val=0.005_dp,    &
         unit_str="bohr",                   &
         type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STO_RMAX",&
         description="Spline rmax",&
         usage="STO_RMAX <value>", &
         default_r_val=19.08622_dp,    &
         unit_str="bohr",                   &
         type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STO_DR",&
         description="Spline dr",&
         usage="STO_DR <value>", &
         default_r_val=0.4_dp,    &
         unit_str="bohr",                   &
         type_of_var=real_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STO_DEBUG",&
         description="Write spline debug info",&
         usage="STO_DEBUG <value>",&
         default_i_val=0,    &
         type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    !> releasing AOM subsection
    CALL section_release(sub2section)
    !-----------------------------------------------------------------------

    !> WAVEFUNCTION_RESTART subsection part of ADIABATIC section
    NULLIFY(sub2section)
    CALL section_create(sub2section,name="WAVEFUNCTION_RESTART",&
         description="WF restart section (in case one wants to use the restart file capability)",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Diabatic WF cofficients in the format:"//&
         "molecule orbital real imaginary",&
         repeats=.TRUE., usage="{{Int} {Int} {Real} {Real}}",&
         type_of_var=lchar_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ACTIVE_STATE_RESTART",&
         description="Active state for the restart (must be the same as the one where the run stopped)",&
         usage="1",  default_i_val=-1, type_of_var=integer_t)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RESTART_KEY",&
         description="Decides whether to use the restart file capability",&
         usage="True", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

!    CALL section_add_subsection(sub2section,sub3section)
!    CALL section_release(sub3section)
!    CALL section_add_subsection(subsection,sub2section)
!    CALL section_release(sub2section)

    CALL section_create(sub2section,name="PRINT",&
         description="Section of possible print options in ADIABATIC env.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of information during the evaluation of "//&
         "the mixed environment. ",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key_section=print_key,&
          name="COEFFICIENTS",&
          description="Controls the output of the wavefunction coefficients (diabatic basis).",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)
    !----------PC
     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="PVECS",&
          description="Controls the output of the p-vectors for PCs visualisation script.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)
     !-------------------

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="POPULATIONS",&
          description="Controls the output of the diabatic populations.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="ADIAB_POPULATIONS",&
          description="Controls the output of the Adiabatic populations.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="ADIABAT_ENERGIES",&
          description="Controls the output of the adiabatic energies.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="SITE_ENERGIES",&
          description="Controls the output of the site energies.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="OFF_DIAGONALS",&
          description="Controls the output of the couplings.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="HAMILTONIAN",&
          description="Controls the output of the Hamiltonian.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="PSEUDO_HAMILTONIAN",&
          description="Controls the output of the Hamiltonian in the pseudo format (to save storage space).",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="ACTIVE_STATE_INFO",&
          description="Print the info relative to the active state (state, energy, wavefunction.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "hartree")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="SH_INFO",&
          description="Prints the info relative to the surface hopping event:"//&
          "Probability, Frustrated Hop, Decoherence, etc.", &
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)


     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="NACE",&
          description="Controls the output of the diabatic NACE.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "fs^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="NACV",&
          description="Controls the output of the diabatic NACV.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "bohr^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="HOP_NACV",&
          description="Controls the output of the adiabatic NACV.",&
          common_iter_levels=1,&
          filename="",         &
          unit_str = "bohr^-1")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="DENSITY_MATRIX",&
          description="Controls the output of the density matrix (diabatic basis).",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="STATE",&
          description="Controls the output of the present state.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

     CALL cp_print_key_section_create(print_key_section=print_key,&
          name="EXACT_FORCES",&
          description="Controls the output of the exact forces in case of ANALYTICS.",&
          common_iter_levels=1,&
          filename="")
     CALL section_add_subsection(sub2section,print_key)
     CALL section_release(print_key)

    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)
!> End CHANGE_AC


    ! Multiple force_eval
    CALL section_create(subsection,name="GENERIC",&
         description="User driven coupling between two or more force_eval.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="MIXING_FUNCTION",&
         description="Specifies the mixing functional form in mathematical notation.",&
         usage="MIXING_FUNCTION (E1+E2-LOG(E1/E2))",  type_of_var=lchar_t,&
         n_var=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VARIABLES",&
         description="Defines the variables of the functional form. To allow an efficient"//&
         " mapping the order of the energy variables will be considered identical to the"//&
         " order of the force_eval in the force_eval_order list.",&
         usage="VARIABLES x",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PARAMETERS",&
         description="Defines the parameters of the functional form",&
         usage="PARAMETERS a b D",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VALUES",&
         description="Defines the values of parameter of the functional form",&
         usage="VALUES ",  type_of_var=real_t,&
         n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="UNITS",&
         description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                     "It is assumed that the corresponding parameter value is specified in this unit.",&
         usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DX",&
         description="Parameter used for computing the derivative with the Ridders' method.",&
         usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ERROR_LIMIT",&
         description="Checks that the error in computing the derivative is not larger than "//&
         "the value set; in case error is larger a warning message is printed.",&
         usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Mapping of atoms
    NULLIFY(sub2section, sub3section)
    CALL section_create(subsection,name="MAPPING",&
         description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."//&
         " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"//&
         " geometrical structure). The mapping is based on defining fragments and the mapping the "//&
         " fragments between the several force_eval and the mixed force_eval",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    ! Mixed force_eval
    CALL section_create(sub2section,name="FORCE_EVAL_MIXED",&
         description="Defines the fragments for the mixed force_eval (reference)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.TRUE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    ! All other force_eval
    CALL section_create(sub2section,name="FORCE_EVAL",&
         description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "//&
         "needs to be provided as parameter)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="DEFINE_FRAGMENTS",&
         description="Specify the fragments definition of the force_eval through the fragments of the"//&
         " force_eval_mixed. This avoids the pedantic definition of the fragments for the force_eval,"//&
         " assuming the order of the fragments for the specified force_eval is the same as the sequence "//&
         " of integers provided. Easier to USE should be preferred to the specification of the single fragments.",&
         usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the force_eval for which fragments and mappings are provided",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAP",&
         description="Provides the index of the fragment of the MIXED force_eval mapped on the"//&
         " locally defined fragment.",&
         usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL create_print_mix_section(subsection)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)
  END SUBROUTINE create_mix_section

! *****************************************************************************
!> \brief Create the print section for mixed
!> \param section the section to create
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_mix_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_type), POINTER              :: print_key

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="print",&
         description="Section of possible print options in MIXED env.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of information during the evaluation of "//&
         "the mixed environment. ",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"DIPOLE",&
         description="Controls the printing of dipole information. "//&
         "Requires the DIPOLE calculation be active for all subforce_eval.", &
         print_level=medium_print_level,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)


  END SUBROUTINE create_print_mix_section

END MODULE input_cp2k_mixed
