# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contract_block.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contract_block.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief routines to contract density matrix blocks with the for center
!>        integrals to yield the Kohn-Sham matrix. The specialized routines
!>        are about 1.2-2.0 as fast as the default one.
!> \par History
!>      10.2009 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE hfx_contract_block
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
# 16 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contract_block.F" 2

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_contract_block'
  PUBLIC :: contract_block
CONTAINS
! *****************************************************************************
!> \brief ...
!> \param ma_max ...
!> \param mb_max ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE contract_block(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: ma_max, mb_max, mc_max, md_max
    REAL(KIND=dp) :: kbd(mb_max*md_max), kbc(mb_max*mc_max), &
      kad(ma_max*md_max), kac(ma_max*mc_max), pbd(mb_max*md_max), &
      pbc(mb_max*mc_max), pad(ma_max*md_max), pac(ma_max*mc_max), &
      prim(ma_max*mb_max*mc_max*md_max), scale


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ma_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(mb_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(mc_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(md_max))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(kbd))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(kbc))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(kad))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(kac))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pbd))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pbc))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pad))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pac))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(prim))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("hfxbase/hfx_contract_block.F",61,"libint not compiled in")
# 4015 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contract_block.F"
  END SUBROUTINE contract_block

# 24687 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfxbase/hfx_contract_block.F"
END MODULE hfx_contract_block
