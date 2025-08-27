# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_device.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_device.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE acc_device




# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/../base/base_uses.f90" 1
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
# 11 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_device.F" 2

  IMPLICIT NONE

  PUBLIC :: acc_get_ndevices, acc_set_active_device

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_device'

# 40 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_device.F"

CONTAINS

! *****************************************************************************
!> \brief Get number of accelerator devices
!> \retval n number of accelerator devices
! *****************************************************************************
  FUNCTION acc_get_ndevices() RESULT(n)
    INTEGER                                  :: n





     n = 0





  END FUNCTION acc_get_ndevices


! *****************************************************************************
!> \brief Set active accelerator device
!> \param dev_id device ID
! *****************************************************************************
  SUBROUTINE acc_set_active_device(dev_id)
    INTEGER :: dev_id

# 80 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/acc/acc_device.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(dev_id))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("acc/acc_device.F",81,"__ACC not compiled in")

  END SUBROUTINE acc_set_active_device

END MODULE acc_device
