# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfx_libint_wrapper.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfx_libint_wrapper.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_libint_wrapper

  USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                             C_F_POINTER,&
                                             C_F_PROCPOINTER,&
                                             C_INT,&
                                             C_LOC,&
                                             C_PTR,&
                                             c_funptr
  USE hfx_libint_wrapper_types,        ONLY: build_deriv1_eri_size,&
                                             build_eri_size,&
                                             lib_deriv,&
                                             lib_int,&
                                             libderiv_max_am1,&
                                             libint_max_am,&
                                             prim_data
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: nco

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/./base/base_uses.f90" 1
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
# 32 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfx_libint_wrapper.F" 2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_libint, terminate_libint,&
         initialize_libderiv, get_eris, get_derivs, terminate_libderiv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

# 264 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/hfx_libint_wrapper.F"

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE TO PROVIDE STUB ROUTINES IN ABSENCE OF __LIBINT **!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param lib ...
!> \param max_am ...
! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(lib))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(max_am))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE initialize_libint

! *****************************************************************************
!> \brief ...
!> \param deriv ...
!> \param max_am ...
! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(deriv))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(max_am))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE initialize_libderiv

! *****************************************************************************
!> \brief ...
!> \param lib ...
! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(lib))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE terminate_libint

! *****************************************************************************
!> \brief ...
!> \param deriv ...
! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(deriv))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE terminate_libderiv

! *****************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param prim ...
!> \param p_work ...
!> \param a_mysize ...
! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_a))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_b))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_c))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_d))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(lib))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(prim))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(p_work))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(a_mysize))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE get_eris

! *****************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param deriv ...
!> \param prim ...
!> \param work_forces ...
!> \param a_mysize ...
! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_a))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_b))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_c))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n_d))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(deriv))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(prim))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(work_forces))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(a_mysize))==-1) EXIT ;  END DO ; ENDIF

  END SUBROUTINE get_derivs


END MODULE hfx_libint_wrapper
