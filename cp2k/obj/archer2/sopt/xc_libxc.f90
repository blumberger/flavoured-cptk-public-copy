# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates a functional from libxc and its derivatives
!> \note
!>      LibXC:
!>      (Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!>      For subsequent versions of libxc, the following should be updated if
!>      necessary:
!>      1) The list of functionals for which it is possible to provide input
!>         parameters (in 'xc_libxc_wrap_functional_set_params'). For more
!>         information on the parameters, see subroutines xc_f90_xxx_set_par
!>         in libxc.f90 of the libxc package or xc_f90_lib_m.F.
!>         only checked for functionals up to 2.0.1
!>      2) Reactivate the functionals which are working correctly
!>         (in 'xc_libxc_wrap_functional_buggy').
!>         only checked for functionals up to 3.x.x
!>
!>      WARNING: In the subroutine libxc_lsd_calc, it could be that the
!>      ordering for the 1st index of v2lapltau, v2rholapl, v2rhotau,
!>      v2sigmalapl and v2sigmatau is not correct. For the moment it does not
!>      matter since the calculation of the 2nd derivatives for meta-GGA
!>      functionals is not implemented in CP2K.
!>
!> \par History
!>      01.2013 created [F. Tran]
!>      07.2014 updates to versions 2.1 [JGH]
!>      08.2015 refactoring [A. Gloess (agloess)]
!> \author F. Tran
! *****************************************************************************
MODULE xc_libxc



  USE bibliography,                    ONLY: Marques2012,&
                                             cite_reference
  
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
# 102 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/../base/base_uses.f90" 1
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
# 104 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc'

  PUBLIC :: libxc_lda_info, libxc_lda_eval, libxc_lsd_info, libxc_lsd_eval, &
            libxc_version_info

CONTAINS

! *****************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_info', &
      routineP = moduleN//':'//routineN

# 211 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",211,"In order to use libxc you need to download and install it")


  END SUBROUTINE libxc_lda_info

! *****************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_info', &
      routineP = moduleN//':'//routineN

# 314 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",314,"In order to use libxc you need to download and install it")


  END SUBROUTINE libxc_lsd_info

! *****************************************************************************
!> \brief info about the LibXC version
!> \param version ...
!> \author A. Gloess (agloess)
! *****************************************************************************
  SUBROUTINE libxc_version_info(version)
    CHARACTER(LEN=*), INTENT(OUT)      :: version ! the string that is output

    CHARACTER(LEN=*), PARAMETER :: routineN = 'libxc_version_info', &
      routineP = moduleN//':'//routineN




    CALL cp__b("xc/xc_libxc.F",333,"In order to use libxc you need to download and install it")


  END SUBROUTINE libxc_version_info

! *****************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_eval', &
      routineP = moduleN//':'//routineN

# 580 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",580,"In order to use libxc you need to download and install it")

  END SUBROUTINE libxc_lda_eval

! *****************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_eval', &
      routineP = moduleN//':'//routineN

# 1121 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",1121,"In order to use libxc you need to download and install it")

  END SUBROUTINE libxc_lsd_eval

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rho density
!> \param norm_drho norm of the gradient of the density
!> \param laplace_rho laplacian of the density
!> \param tau kinetic-energy density
!> \param e_0 energy density
!> \param e_rho derivative of the energy density with respect to rho
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_laplace_rho derivative of the energy density with respect to laplace_rho
!> \param e_tau derivative of the energy density with respect to tau
!> \param e_rho_rho derivative of the energy density with respect to rho_rho
!> \param e_ndrho_rho derivative of the energy density with respect to ndrho_rho
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_rho_laplace_rho derivative of the energy density with respect to rho_laplace_rho
!> \param e_rho_tau derivative of the energy density with respect to rho_tau
!> \param e_ndrho_laplace_rho derivative of the energy density with respect to ndrho_laplace_rho
!> \param e_ndrho_tau derivative of the energy density with respect to ndrho_tau
!> \param e_laplace_rho_laplace_rho derivative of the energy density with respect to laplace_rho_laplace_rho
!> \param e_laplace_rho_tau derivative of the energy density with respect to laplace_rho_tau
!> \param e_tau_tau derivative of the energy density with respect to tau_tau
!> \param e_rho_rho_rho derivative of the energy density with respect to rho_rho_rho
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param epsilon_rho ...
!> \param epsilon_tau ...
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_calc(rho,norm_drho,laplace_rho,tau,&
          e_0,e_rho,e_ndrho,e_laplace_rho,e_tau,e_rho_rho,e_ndrho_rho,&
          e_ndrho_ndrho,e_rho_laplace_rho,e_rho_tau,e_ndrho_laplace_rho,&
          e_ndrho_tau,e_laplace_rho_laplace_rho,e_laplace_rho_tau,&
          e_tau_tau,e_rho_rho_rho,&
          grad_deriv,npoints,epsilon_rho,&
          epsilon_tau,func_name,sc,params)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho, laplace_rho, &
                                                tau
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho, &
      e_laplace_rho, e_tau, e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, &
      e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, e_ndrho_tau, &
      e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_tau_tau, e_rho_rho_rho
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_calc', &
      routineP = moduleN//':'//routineN

# 1478 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",1478,"In order to use libxc you need to download and install it")


  END SUBROUTINE libxc_lda_calc

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rhoa alpha density
!> \param rhob beta density
!> \param norm_drho ...
!> \param norm_drhoa norm of the gradient of the alpha density
!> \param norm_drhob norm of the gradient of the beta density
!> \param laplace_rhoa laplacian of the alpha density
!> \param laplace_rhob laplacian of the beta density
!> \param tau_a alpha kinetic-energy density
!> \param tau_b beta kinetic-energy density
!> \param e_0 energy density
!> \param e_rhoa derivative of the energy density with respect to rhoa
!> \param e_rhob derivative of the energy density with respect to rhob
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_ndrhoa derivative of the energy density with respect to ndrhoa
!> \param e_ndrhob derivative of the energy density with respect to ndrhob
!> \param e_laplace_rhoa derivative of the energy density with respect to laplace_rhoa
!> \param e_laplace_rhob derivative of the energy density with respect to laplace_rhob
!> \param e_tau_a derivative of the energy density with respect to tau_a
!> \param e_tau_b derivative of the energy density with respect to tau_b
!> \param e_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa
!> \param e_rhoa_rhob derivative of the energy density with respect to rhoa_rhob
!> \param e_rhob_rhob derivative of the energy density with respect to rhob_rhob
!> \param e_ndrho_rhoa derivative of the energy density with respect to ndrho_rhoa
!> \param e_ndrho_rhob derivative of the energy density with respect to ndrho_rhob
!> \param e_ndrhoa_rhoa derivative of the energy density with respect to ndrhoa_rhoa
!> \param e_ndrhoa_rhob derivative of the energy density with respect to ndrhoa_rhob
!> \param e_ndrhob_rhoa derivative of the energy density with respect to ndrhob_rhoa
!> \param e_ndrhob_rhob derivative of the energy density with respect to ndrhob_rhob
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_ndrho_ndrhoa derivative of the energy density with respect to ndrho_ndrhoa
!> \param e_ndrho_ndrhob derivative of the energy density with respect to ndrho_ndrhob
!> \param e_ndrhoa_ndrhoa derivative of the energy density with respect to ndrhoa_ndrhoa
!> \param e_ndrhoa_ndrhob derivative of the energy density with respect to ndrhoa_ndrhob
!> \param e_ndrhob_ndrhob derivative of the energy density with respect to ndrhob_ndrhob
!> \param e_rhoa_laplace_rhoa derivative of the energy density with respect to rhoa_laplace_rhoa
!> \param e_rhoa_laplace_rhob derivative of the energy density with respect to rhoa_laplace_rhob
!> \param e_rhob_laplace_rhoa derivative of the energy density with respect to rhob_laplace_rhoa
!> \param e_rhob_laplace_rhob derivative of the energy density with respect to rhob_laplace_rhob
!> \param e_rhoa_tau_a derivative of the energy density with respect to rhoa_tau_a
!> \param e_rhoa_tau_b derivative of the energy density with respect to rhoa_tau_b
!> \param e_rhob_tau_a derivative of the energy density with respect to rhob_tau_a
!> \param e_rhob_tau_b derivative of the energy density with respect to rhob_tau_b
!> \param e_ndrho_laplace_rhoa derivative of the energy density with respect to ndrho_laplace_rhoa
!> \param e_ndrho_laplace_rhob derivative of the energy density with respect to ndrho_laplace_rhob
!> \param e_ndrhoa_laplace_rhoa derivative of the energy density with respect to ndrhoa_laplace_rhoa
!> \param e_ndrhoa_laplace_rhob derivative of the energy density with respect to ndrhoa_laplace_rhob
!> \param e_ndrhob_laplace_rhoa derivative of the energy density with respect to ndrhob_laplace_rhoa
!> \param e_ndrhob_laplace_rhob derivative of the energy density with respect to ndrhob_laplace_rhob
!> \param e_ndrho_tau_a derivative of the energy density with respect to ndrho_tau_a
!> \param e_ndrho_tau_b derivative of the energy density with respect to ndrho_tau_b
!> \param e_ndrhoa_tau_a derivative of the energy density with respect to ndrhoa_tau_a
!> \param e_ndrhoa_tau_b derivative of the energy density with respect to ndrhoa_tau_b
!> \param e_ndrhob_tau_a derivative of the energy density with respect to ndrhob_tau_a
!> \param e_ndrhob_tau_b derivative of the energy density with respect to ndrhob_tau_b
!> \param e_laplace_rhoa_laplace_rhoa derivative of the energy density with respect to laplace_rhoa_laplace_rhoa
!> \param e_laplace_rhoa_laplace_rhob derivative of the energy density with respect to laplace_rhoa_laplace_rhob
!> \param e_laplace_rhob_laplace_rhob derivative of the energy density with respect to laplace_rhob_laplace_rhob
!> \param e_laplace_rhoa_tau_a derivative of the energy density with respect to laplace_rhoa_tau_a
!> \param e_laplace_rhoa_tau_b derivative of the energy density with respect to laplace_rhoa_tau_b
!> \param e_laplace_rhob_tau_a derivative of the energy density with respect to laplace_rhob_tau_a
!> \param e_laplace_rhob_tau_b derivative of the energy density with respect to laplace_rhob_tau_b
!> \param e_tau_a_tau_a derivative of the energy density with respect to tau_a_tau_a
!> \param e_tau_a_tau_b derivative of the energy density with respect to tau_a_tau_b
!> \param e_tau_b_tau_b derivative of the energy density with respect to tau_b_tau_b
!> \param e_rhoa_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa_rhoa
!> \param e_rhoa_rhoa_rhob derivative of the energy density with respect to rhoa_rhoa_rhob
!> \param e_rhoa_rhob_rhob derivative of the energy density with respect to rhoa_rhob_rhob
!> \param e_rhob_rhob_rhob derivative of the energy density with respect to rhob_rhob_rhob
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param epsilon_rho ...
!> \param epsilon_tau ...
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_calc(rhoa,rhob,norm_drho,norm_drhoa,&
          norm_drhob,laplace_rhoa,laplace_rhob,tau_a,tau_b,&
          e_0,e_rhoa,e_rhob,e_ndrho,e_ndrhoa,e_ndrhob,&
          e_laplace_rhoa,e_laplace_rhob,e_tau_a,e_tau_b,&
          e_rhoa_rhoa,e_rhoa_rhob,e_rhob_rhob,&
          e_ndrho_rhoa,e_ndrho_rhob,e_ndrhoa_rhoa,&
          e_ndrhoa_rhob,e_ndrhob_rhoa,e_ndrhob_rhob,&
          e_ndrho_ndrho,e_ndrho_ndrhoa,e_ndrho_ndrhob,&
          e_ndrhoa_ndrhoa,e_ndrhoa_ndrhob,e_ndrhob_ndrhob,&
          e_rhoa_laplace_rhoa,e_rhoa_laplace_rhob,&
          e_rhob_laplace_rhoa,e_rhob_laplace_rhob,&
          e_rhoa_tau_a,e_rhoa_tau_b,e_rhob_tau_a,e_rhob_tau_b,&
          e_ndrho_laplace_rhoa,e_ndrho_laplace_rhob,&
          e_ndrhoa_laplace_rhoa,e_ndrhoa_laplace_rhob,&
          e_ndrhob_laplace_rhoa,e_ndrhob_laplace_rhob,&
          e_ndrho_tau_a,e_ndrho_tau_b,&
          e_ndrhoa_tau_a,e_ndrhoa_tau_b,&
          e_ndrhob_tau_a,e_ndrhob_tau_b,&
          e_laplace_rhoa_laplace_rhoa,&
          e_laplace_rhoa_laplace_rhob,&
          e_laplace_rhob_laplace_rhob,&
          e_laplace_rhoa_tau_a,e_laplace_rhoa_tau_b,&
          e_laplace_rhob_tau_a,e_laplace_rhob_tau_b,&
          e_tau_a_tau_a,e_tau_a_tau_b,e_tau_b_tau_b,&
          e_rhoa_rhoa_rhoa,e_rhoa_rhoa_rhob,&
          e_rhoa_rhob_rhob,e_rhob_rhob_rhob,&
          grad_deriv,npoints,epsilon_rho,&
          epsilon_tau,func_name,sc,params)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob, &
                                                laplace_rhoa, laplace_rhob, &
                                                tau_a, tau_b
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rhoa, e_rhob, &
      e_ndrho, e_ndrhoa, e_ndrhob, e_laplace_rhoa, e_laplace_rhob, e_tau_a, &
      e_tau_b, e_rhoa_rhoa, e_rhoa_rhob, e_rhob_rhob, e_ndrho_rhoa, &
      e_ndrho_rhob, e_ndrhoa_rhoa, e_ndrhoa_rhob, e_ndrhob_rhoa, &
      e_ndrhob_rhob, e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, &
      e_ndrhoa_ndrhoa, e_ndrhoa_ndrhob, e_ndrhob_ndrhob, e_rhoa_laplace_rhoa, &
      e_rhoa_laplace_rhob, e_rhob_laplace_rhoa, e_rhob_laplace_rhob, &
      e_rhoa_tau_a, e_rhoa_tau_b, e_rhob_tau_a, e_rhob_tau_b, &
      e_ndrho_laplace_rhoa, e_ndrho_laplace_rhob, e_ndrhoa_laplace_rhoa
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_ndrhoa_laplace_rhob, &
      e_ndrhob_laplace_rhoa, e_ndrhob_laplace_rhob, e_ndrho_tau_a, &
      e_ndrho_tau_b, e_ndrhoa_tau_a, e_ndrhoa_tau_b, e_ndrhob_tau_a, &
      e_ndrhob_tau_b, e_laplace_rhoa_laplace_rhoa, &
      e_laplace_rhoa_laplace_rhob, e_laplace_rhob_laplace_rhob, &
      e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, e_laplace_rhob_tau_a, &
      e_laplace_rhob_tau_b, e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b_tau_b, &
      e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, e_rhoa_rhob_rhob, e_rhob_rhob_rhob
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_calc', &
      routineP = moduleN//':'//routineN

# 2308 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc.F"
    CALL cp__b("xc/xc_libxc.F",2308,"In order to use libxc you need to download and install it")


  END SUBROUTINE libxc_lsd_calc

END MODULE xc_libxc
