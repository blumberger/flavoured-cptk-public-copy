# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/molden_utils.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/molden_utils.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Functions handling the MOLDEN format. Split from mode_selective.
!> \author Teodoro Laino, 03.2009
! *****************************************************************************
MODULE molden_utils
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_log_handling,                 ONLY: cp_logger_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/molden_utils.F" 2

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molden_utils'
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  PUBLIC :: molden_out

CONTAINS
! *****************************************************************************
!> \brief writes the output for vibrational analysis in MOLDEN format
!> \param input ...
!> \param particles ...
!> \param freq ...
!> \param eigen_vec ...
!> \param intensities ...
!> \param calc_intens ...
!> \param dump_only_positive ...
!> \param logger ...
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE molden_out(input,particles,freq,eigen_vec,intensities,calc_intens,&
             dump_only_positive,logger)

    TYPE(section_vals_type), POINTER         :: input
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:)              :: freq
    REAL(KIND=dp), DIMENSION(:, :)           :: eigen_vec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: intensities
    LOGICAL                                  :: calc_intens, &
                                                dump_only_positive
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'molden_out', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: i, iw, j, k, l

    iw=cp_print_key_unit_nr(logger,input,"VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB",&
            extension=".mol",file_status='REPLACE')

    IF(iw.GT.0)THEN
       IF(.NOT.(MOD(SIZE(eigen_vec,1),3)==0))CALL cp__a("molden_utils.F",64)
       IF(.NOT.(SIZE(particles)==SIZE(eigen_vec,1)/3))CALL cp__a("molden_utils.F",65)
       IF(.NOT.(SIZE(freq,1)==SIZE(eigen_vec,2)))CALL cp__a("molden_utils.F",66)
       WRITE(iw,'(T2,A)')"[Molden Format]"
       WRITE(iw,'(T2,A)')"[FREQ]"
       DO i=1,SIZE(freq,1)
          IF((.NOT.dump_only_positive).OR.(freq(i)>=0._dp))WRITE(iw,'(T5,F12.6)') freq(i)
       END DO
       WRITE(iw,'(T2,A)')"[FR-COORD]"
       DO i=1,SIZE(particles)
          CALL get_atomic_kind(atomic_kind=particles(i)%atomic_kind,&
               element_symbol=element_symbol)
               WRITE(iw,'(T2,A2,3X,3(F12.6,3X))')&
                    element_symbol, particles((i))%r(:)
       END DO
       WRITE(iw,'(T2,A)')"[FR-NORM-COORD]"
       l=0
       DO i=1,SIZE(eigen_vec,2)
          IF ((.NOT.dump_only_positive).OR.(freq(i)>=0._dp)) THEN
             l=l+1
             WRITE(iw,'(T2,A,1X,I6)')"vibration",l
             DO j=1,SIZE(eigen_vec,1)/3
                k=(j-1)*3
                WRITE(iw,'(T2,3(F12.6,3X))')eigen_vec(k+1,i),eigen_vec(k+2,i),eigen_vec(k+3,i)
             END DO
          END IF
       END DO
       IF(calc_intens)THEN
          WRITE(iw,'(T2,A)')"[INT]"
          DO i=1,SIZE(intensities)
             IF((.NOT.dump_only_positive).OR.(freq(i)>=0._dp))WRITE(iw,'(3X,F18.6)')intensities(i)
          END DO
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,input,"VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB")
  END SUBROUTINE molden_out

END MODULE molden_utils
