# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc_wrap.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc_wrap.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Includes all necessary routines, functions and parameters from
!>        libxc. Provides CP2K routines/functions where the LibXC calling list
!>        is version dependent. The naming convention for such 
!>        routines/functions is xc_f90_XXX --> 'xc_libxc_wrap_XXX'. All version
!>        independent routines/functions are just bypassed to higher level
!>        module file 'xc_libxc'.
!>        
!> \note For LibXC versions 2.2.2 and above.
!>       Once the LibXC-API is stable, remove all 'xc_libxc_wrap_XXX'
!>       routines/functions and use 'xc_f90_lib_m' directly in 'xc_libxc'.
!>       Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!> \par History
!>      08.2015 created [A. Gloess]
!> \author A. Gloessa (agloess)
! *****************************************************************************
MODULE xc_libxc_wrap
# 648 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/xc/xc_libxc_wrap.F"
END MODULE xc_libxc_wrap
