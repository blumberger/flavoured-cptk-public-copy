# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/optbas_opt_utils.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/optbas_opt_utils.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optbas_opt_utils
  USE admm_types,                      ONLY: admm_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_trace,&
                                             cp_fm_upper_to_full
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type

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
# 22 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/optbas_opt_utils.F" 2

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: evaluate_fval, evaluate_energy

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_opt_utils'

CONTAINS

   
! *****************************************************************************
!> \brief ...
!> \param mos ...
!> \param matrix_ks ...
!> \param S_inv_orb ...
!> \param Q ...
!> \param tmp1 ...
!> \param energy ...
! *****************************************************************************
  SUBROUTINE evaluate_energy(mos,matrix_ks,S_inv_orb,Q,tmp1,energy)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_type), POINTER                :: S_inv_orb, Q, tmp1
    REAL(KIND=dp)                            :: energy

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, naux, nmo, norb
    REAL(KIND=dp)                            :: tmp_energy
    TYPE(cp_fm_type), POINTER                :: mo_coeff, QS_inv, tmp, tmp2, &
                                                work, work_orb

    CALL cp_fm_create(QS_inv,matrix_struct=Q%matrix_struct)
    CALL cp_fm_create(tmp,matrix_struct=Q%matrix_struct) 
    CALL cp_fm_create(tmp2,matrix_struct=tmp1%matrix_struct)
    CALL cp_fm_create(work,matrix_struct=S_inv_orb%matrix_struct)
    CALL cp_fm_create(work_orb,matrix_struct=S_inv_orb%matrix_struct)
    CALL cp_fm_get_info(Q,nrow_global=naux,ncol_global=norb)
    CALL cp_gemm('N','N',naux,norb,norb,1.0_dp,Q,S_inv_orb,0.0_dp,QS_inv)
    energy=0.0_dp 
    DO ispin=1,SIZE(matrix_ks)
       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,work)
       CALL cp_fm_upper_to_full(work,work_orb)

       CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo,mo_coeff=mo_coeff)
       CALL cp_gemm('N','N',naux,norb,norb,1.0_dp,QS_inv,work,0.0_dp,tmp)
       CALL cp_gemm('N','T',naux,naux,norb,1.0_dp,tmp,QS_inv,0.0_dp,tmp1)
       CALL cp_gemm('N','T',naux,naux,nmo,1.0_dp,mo_coeff,mo_coeff,0.0_dp,tmp2)
       CALL cp_fm_trace(tmp1,tmp2,tmp_energy)
       energy=energy+tmp_energy*(3.0_dp-REAL(SIZE(matrix_ks),dp))
       
    END DO
    
    CALL cp_fm_release(work_orb)
    CALL cp_fm_release(QS_inv)
    CALL cp_fm_release(tmp)
    CALL cp_fm_release(tmp2)
    CALL cp_fm_release(work)
    CALL cp_fm_release(work_orb)

  END SUBROUTINE evaluate_energy

! *****************************************************************************
!> \brief ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param Q ...
!> \param Snew ...
!> \param admm_env ...
!> \param fval ...
!> \param S_cond_number ...
! *****************************************************************************
  SUBROUTINE evaluate_fval(mos,mos_aux_fit,Q,Snew,admm_env,fval,S_cond_number)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: Q, Snew
    TYPE(admm_type), POINTER                 :: admm_env
    REAL(KIND=dp)                            :: fval, S_cond_number

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_fval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nao_aux_fit, nao_orb, &
                                                nmo, nspins
    REAL(KIND=dp)                            :: trace
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    CALL copy_dbcsr_to_fm(Q,admm_env%Q)
    fval=0.0_dp
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,-2.0_dp,admm_env%Q,mo_coeff,&
                      0.0_dp,admm_env%work_aux_nmo(ispin)%matrix)
      CALL cp_fm_trace(mo_coeff_aux_fit,admm_env%work_aux_nmo(ispin)%matrix,trace)
      fval=fval+trace+2.0_dp*nmo
    END DO

    ALLOCATE(eigenvalues(nao_aux_fit))
    CALL copy_dbcsr_to_fm(Snew,admm_env%work_aux_aux)
    CALL cp_fm_syevd(admm_env%work_aux_aux,admm_env%work_aux_aux2,eigenvalues)
    S_cond_number=MAXVAL(ABS(eigenvalues))/MAX(MINVAL(ABS(eigenvalues)),EPSILON(0.0_dp))
    DEALLOCATE(eigenvalues)

  END SUBROUTINE evaluate_fval

END MODULE optbas_opt_utils
