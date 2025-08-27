# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_diis_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_diis_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief buffer for the diis of the scf
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author Matthias Krack
! *****************************************************************************
MODULE qs_diis_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_release
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE kinds,                           ONLY: dp

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
# 19 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_diis_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_diis_types'

  PUBLIC :: qs_diis_buffer_type
  PUBLIC :: qs_diis_b_release
  PUBLIC :: qs_diis_buffer_type_sparse, &
            qs_diis_b_release_sparse

! *****************************************************************************
!> \brief keeps a buffer with the previous values of s,p,k
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author Matthias Krack
! *****************************************************************************
  TYPE qs_diis_buffer_type
    INTEGER                                          :: nbuffer,ncall,&
         id_nr,ref_count
    TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: error,PARAMETER
    REAL(KIND = dp), DIMENSION(:,:), POINTER                :: b_matrix
  END TYPE qs_diis_buffer_type

! *****************************************************************************
!> \brief build array of pointers to diis buffers
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE qs_diis_buffer_p_type
     TYPE(qs_diis_buffer_type), POINTER :: diis_buffer
  END TYPE qs_diis_buffer_p_type

! *****************************************************************************
!> \brief build array of pointers to diis buffers for sparse matrix case
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      10.2014 Modified from non-sparse case by Fredy W. Aquino
!> \author fwaq
! *****************************************************************************
  TYPE qs_diis_buffer_type_sparse
    INTEGER                                          :: nbuffer,ncall,&
                                                        id_nr,ref_count
    TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER   :: error,PARAMETER
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: b_matrix
  END TYPE qs_diis_buffer_type_sparse

  TYPE qs_diis_buffer_p_type_sparse
     TYPE(qs_diis_buffer_type_sparse), POINTER       :: diis_buffer
  END TYPE qs_diis_buffer_p_type_sparse

CONTAINS

! *****************************************************************************
!> \brief retains a diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to retain
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE qs_diis_b_retain(diis_buffer)
    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_retain', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(diis_buffer)))CALL cp__a("qs_diis_types.F",88)
  IF(.NOT.(diis_buffer%ref_count>0))CALL cp__a("qs_diis_types.F",89)
  diis_buffer%ref_count=diis_buffer%ref_count+1
END SUBROUTINE qs_diis_b_retain

! *****************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE qs_diis_b_release(diis_buffer)
    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

  IF (ASSOCIATED(diis_buffer)) THEN
     IF(.NOT.(diis_buffer%ref_count>0))CALL cp__a("qs_diis_types.F",109)
     diis_buffer%ref_count=diis_buffer%ref_count-1
     IF (diis_buffer%ref_count<1) THEN
        IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
           DEALLOCATE(diis_buffer%b_matrix)
        END IF
        IF (ASSOCIATED(diis_buffer%error)) THEN
           DO j=1,SIZE(diis_buffer%error,2)
              DO i=1,SIZE(diis_buffer%error,1)
                 CALL cp_fm_release(diis_buffer%error(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%error)
        END IF
        IF (ASSOCIATED(diis_buffer%parameter)) THEN
           DO j=1,SIZE(diis_buffer%parameter,2)
              DO i=1,SIZE(diis_buffer%parameter,1)
                 CALL cp_fm_release(diis_buffer%parameter(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%parameter)
        END IF
        DEALLOCATE(diis_buffer)
     END IF
  END IF
END SUBROUTINE qs_diis_b_release

! *****************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_release
!> \author Fredy W. Aquino
! *****************************************************************************
SUBROUTINE qs_diis_b_release_sparse(diis_buffer)

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_release_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

  IF (ASSOCIATED(diis_buffer)) THEN
        IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
           DEALLOCATE(diis_buffer%b_matrix)
        END IF
        IF (ASSOCIATED(diis_buffer%error)) THEN
           DO j=1,SIZE(diis_buffer%error,2)
              DO i=1,SIZE(diis_buffer%error,1)
               CALL cp_dbcsr_release(diis_buffer%error(i,j)%matrix)
               DEALLOCATE(diis_buffer%error(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%error)
        END IF
        IF (ASSOCIATED(diis_buffer%parameter)) THEN
           DO j=1,SIZE(diis_buffer%parameter,2)
              DO i=1,SIZE(diis_buffer%parameter,1)
               CALL cp_dbcsr_release(diis_buffer%parameter(i,j)%matrix)
               DEALLOCATE(diis_buffer%parameter(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%parameter)
        END IF
        DEALLOCATE(diis_buffer)
  END IF
END SUBROUTINE qs_diis_b_release_sparse

END MODULE qs_diis_types
