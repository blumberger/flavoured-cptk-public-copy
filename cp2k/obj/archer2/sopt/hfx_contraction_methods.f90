# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contraction_methods.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contraction_methods.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains routines for contraction without dgemms. PLEASE DO NOT MODIFY.
!> \notes Contains specific routines for contraction. The compiler flag
!>        -D__MAX_CONTR defines the maximum angular momentum up to which
!>        specialized code will be used. Default setting is d-functions.
!>        Increasing -D__MAX_CONTR produces faster code but might overburden
!>        the optimization capabilities of some poor compilers.
!>        This file contains specific code up to g-functions. If you need more
!>        look at cp2k/tools/hfx_tools/contraction/
!> \par History
!>      07.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

MODULE hfx_contraction_methods

!** This defines the default behaviour




  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/../base/base_uses.f90" 1
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
# 29 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contraction_methods.F" 2

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: contract


  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param work ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param primitives ...
!> \param buffer1 ...
!> \param buffer2 ...
! *****************************************************************************
  SUBROUTINE  contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                       n_a, n_b, n_c, n_d,nl_a, nl_b, nl_c, nl_d, work,&
                       sphi_a, sphi_b, sphi_c, sphi_d,&
                       primitives, &
                       buffer1, buffer2)

    INTEGER, INTENT(IN)         :: ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod,&
                                   n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d
    REAL(dp), DIMENSION(ncoa*ncob* ncoc* ncod), INTENT(IN) :: work
    REAL(dp), DIMENSION(ncoa,nsoa*nl_a), INTENT(IN)   :: sphi_a
    REAL(dp), DIMENSION(ncob,nsob*nl_b), INTENT(IN)   :: sphi_b
    REAL(dp), DIMENSION(ncoc,nsoc*nl_c), INTENT(IN)   :: sphi_c
    REAL(dp), DIMENSION(ncod,nsod*nl_d), INTENT(IN)   :: sphi_d

    REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b,nsoc*nl_c,nsod*nl_d) :: primitives
    REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)  :: buffer1, buffer2


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ncoa))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ncob))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ncoc))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ncod))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nsoa))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nsob))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nsoc))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nsod))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_a))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_b))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_c))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_d))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nl_a))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nl_b))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nl_c))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nl_d))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(sphi_a))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(sphi_b))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(sphi_c))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(sphi_d))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(work))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(primitives))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(buffer1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(buffer2))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("hfxbase/hfx_contraction_methods.F",108,"libint not compiled in")
# 12607 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contraction_methods.F"
  END SUBROUTINE  contract
# 107603 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contraction_methods.F"
END MODULE hfx_contraction_methods
