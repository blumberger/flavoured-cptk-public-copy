# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility routines for the memory handling.
!>        Uses stop_memory to handle exceptions
!> \par History
!>      none
!> \author Matthias Krack (25.06.1999)
! *****************************************************************************
MODULE memory_utilities

  
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_8,&
                                             int_8_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/../base/base_uses.f90" 1
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
# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'memory_utilities'

  PUBLIC :: reallocate

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_c1,reallocate_c2,reallocate_c3,reallocate_c4,&
                     reallocate_i1,reallocate_i2,reallocate_i3,reallocate_i4,&
                     reallocate_r1,reallocate_r2,reallocate_r3,reallocate_r4,&
                     reallocate_r5,reallocate_s1,reallocate_l1,reallocate_8i1,&
                     reallocate_8i2
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief (Re)Allocate a complex vector with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! *****************************************************************************
  SUBROUTINE reallocate_c1(p,lb1_new,ub1_new)

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c1', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: zero = (0.0_dp,0.0_dp)
    INTEGER, PARAMETER                       :: t_size = 2*dp_size

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 67 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_c1

! *****************************************************************************
!> \brief (Re)Allocate a complex array 2D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! *****************************************************************************
  SUBROUTINE reallocate_c2(p,lb1_new,ub1_new,lb2_new,ub2_new)

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c2', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: zero = (0.0_dp,0.0_dp)
    INTEGER, PARAMETER                       :: t_size = 2*dp_size

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_2.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, ub1, ub1_old, ub2, ub2_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,15,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1))

       END IF
       work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,26,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1))
    END IF
    p(:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
       DEALLOCATE (work)
    END IF
# 97 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_c2

! *****************************************************************************
!> \brief (Re)Allocate a complex array 3D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! *****************************************************************************
  SUBROUTINE reallocate_c3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c3', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: zero = (0.0_dp,0.0_dp)
    INTEGER, PARAMETER                       :: t_size = 2*dp_size

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_3.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,19,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3) = p(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,30,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1))
    END IF
    p(:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (work)
    END IF
# 129 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_c3

! *****************************************************************************
!> \brief (Re)Allocate a complex array 4D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! *****************************************************************************
  SUBROUTINE reallocate_c4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :), POINTER         :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c4', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: zero = (0.0_dp,0.0_dp)
    INTEGER, PARAMETER                       :: t_size = 2*dp_size

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_4.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, ub4_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb4_old = LBOUND(p,4)
       ub4_old = UBOUND(p,4)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       lb4 = MAX(lb4_new,lb4_old)
       ub4 = MIN(ub4_new,ub4_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,23,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1)*&
                                         (ub4-lb4+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new),&
              STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,36,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1)*&
                                   (ub4_new-lb4_new+1))
    END IF
    p(:,:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (work)
    END IF
# 165 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_c4

! *****************************************************************************
!> \brief (Re)Allocate an integer vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_i1(p,lb1_new,ub1_new)

    INTEGER, DIMENSION(:), POINTER           :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i1', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = int_size, zero = 0

    INTEGER, ALLOCATABLE, DIMENSION(:)       :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 189 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_i1

! *****************************************************************************
!> \brief (Re)Allocate an integer array with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_i2(p,lb1_new,ub1_new,lb2_new,ub2_new)

    INTEGER, DIMENSION(:, :), POINTER        :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i2', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = int_size, zero = 0

    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_2.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, ub1, ub1_old, ub2, ub2_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,15,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1))

       END IF
       work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,26,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1))
    END IF
    p(:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
       DEALLOCATE (work)
    END IF
# 216 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_i2

! *****************************************************************************
!> \brief (Re)Allocate an integer array 3D with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_i3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

    INTEGER, DIMENSION(:, :, :), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i3', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = int_size, zero = 0

    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_3.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,19,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3) = p(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,30,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1))
    END IF
    p(:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (work)
    END IF
# 245 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_i3

! *****************************************************************************
!> \brief (Re)Allocate an integer array 3D with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_i4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

    INTEGER, DIMENSION(:, :, :, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i4', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = int_size, zero = 0

    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_4.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, ub4_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb4_old = LBOUND(p,4)
       ub4_old = UBOUND(p,4)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       lb4 = MAX(lb4_new,lb4_old)
       ub4 = MIN(ub4_new,ub4_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,23,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1)*&
                                         (ub4-lb4+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new),&
              STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,36,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1)*&
                                   (ub4_new-lb4_new+1))
    END IF
    p(:,:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (work)
    END IF
# 279 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_i4

! *****************************************************************************
!> \brief (Re)Allocate an integer (int_8) vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_8i1(p,lb1_new,ub1_new)

    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_8i1', &
      routineP = moduleN//':'//routineN
    INTEGER(KIND=int_8), PARAMETER           :: zero = 0
    INTEGER, PARAMETER                       :: t_size = int_8_size

    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 306 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_8i1

! *****************************************************************************
!> \brief (Re)Allocate an integer (int_8) array with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_8i2(p,lb1_new,ub1_new,lb2_new,ub2_new)

    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_8i2', &
      routineP = moduleN//':'//routineN
    INTEGER(KIND=int_8), PARAMETER           :: zero = 0
    INTEGER, PARAMETER                       :: t_size = int_8_size

    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_2.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, ub1, ub1_old, ub2, ub2_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,15,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1))

       END IF
       work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,26,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1))
    END IF
    p(:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
       DEALLOCATE (work)
    END IF
# 336 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_8i2

! *****************************************************************************
!> \brief (Re)Allocate an real vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_r1(p,lb1_new,ub1_new)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r1', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = dp_size
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 361 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_r1

! *****************************************************************************
!> \brief (Re)Allocate an real array with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_r2(p,lb1_new,ub1_new,lb2_new,ub2_new)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r2', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = dp_size
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_2.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, ub1, ub1_old, ub2, ub2_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,15,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1))

       END IF
       work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,26,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1))
    END IF
    p(:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
       DEALLOCATE (work)
    END IF
# 390 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_r2

! *****************************************************************************
!> \brief (Re)Allocate an real array 3D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_r3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r3', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = dp_size
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_3.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,19,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3) = p(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,30,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1))
    END IF
    p(:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
       DEALLOCATE (work)
    END IF
# 422 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_r3

! *****************************************************************************
!> \brief (Re)Allocate an real array 4D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_r4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r4', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = dp_size
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_4.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, ub4_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb4_old = LBOUND(p,4)
       ub4_old = UBOUND(p,4)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       lb4 = MAX(lb4_new,lb4_old)
       ub4 = MIN(ub4_new,ub4_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,23,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1)*&
                                         (ub4-lb4+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new),&
              STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,36,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1)*&
                                   (ub4_new-lb4_new+1))
    END IF
    p(:,:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
       DEALLOCATE (work)
    END IF
# 458 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_r4

! *****************************************************************************
!> \brief (Re)Allocate an real array 5D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \param lb5_new ...
!> \param ub5_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_r5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new,lb5_new,ub5_new)

    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new, lb5_new, &
                                                ub5_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r5', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = dp_size
    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_5.f90" 1
    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      lb5, lb5_old, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, ub4_old, ub5, ub5_old

    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb3_old = LBOUND(p,3)
       ub3_old = UBOUND(p,3)
       lb4_old = LBOUND(p,4)
       ub4_old = UBOUND(p,4)
       lb5_old = LBOUND(p,5)
       ub5_old = UBOUND(p,5)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)
       lb3 = MAX(lb3_new,lb3_old)
       ub3 = MIN(ub3_new,ub3_old)
       lb4 = MAX(lb4_new,lb4_old)
       ub4 = MIN(ub4_new,ub4_old)
       lb5 = MAX(lb5_new,lb5_old)
       ub5 = MIN(ub5_new,ub5_old)
       ALLOCATE (work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,27,&
                           "work",t_size*(ub1-lb1+1)*&
                                         (ub2-lb2+1)*&
                                         (ub3-lb3+1)*&
                                         (ub4-lb4+1)*&
                                         (ub5-lb5+1))
       END IF
       work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5) = p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new,lb5_new:ub5_new),&
              STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,41,&
                        "p",t_size*(ub1_new-lb1_new+1)*&
                                   (ub2_new-lb2_new+1)*&
                                   (ub3_new-lb3_new+1)*&
                                   (ub4_new-lb4_new+1)*&
                                   (ub5_new-lb5_new+1))
    END IF
    p(:,:,:,:,:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5)
       DEALLOCATE (work)
    END IF
# 497 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_r5

! *****************************************************************************
!> \brief (Re)Allocate a logical vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! *****************************************************************************
  SUBROUTINE reallocate_l1(p,lb1_new,ub1_new)

    LOGICAL, DIMENSION(:), POINTER           :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_l1', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: t_size = int_size
    LOGICAL, PARAMETER                       :: zero = .FALSE.

    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 522 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_l1

! *****************************************************************************
!> \brief (Re)Allocate a vector of string variables
!> \param p_short ...
!> \param lb_new ...
!> \param ub_new ...
!> \param p_long ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! *****************************************************************************
  SUBROUTINE reallocate_s1(p_short,lb_new,ub_new,p_long)

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: p_short
    INTEGER, INTENT(IN)                      :: lb_new, ub_new
    CHARACTER(LEN=default_path_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: p_long

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_s1', &
      routineP = moduleN//':'//routineN

    IF      (PRESENT(p_short)) THEN
       CALL reallocate_ss1(p_short,lb_new,ub_new)
    ELSE IF (PRESENT(p_long))  THEN
       CALL reallocate_ls1(p_long,lb_new,ub_new)
    ELSE
       CALL cp__b("common/memory_utilities.F",554,"At least one of the two optional arguments is required")
    END IF

  END SUBROUTINE reallocate_s1

! *****************************************************************************
!> \brief (Re)Allocate a vector of string variables of default string length.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Thomas Chassaing
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! *****************************************************************************
  SUBROUTINE reallocate_ss1(p,lb1_new,ub1_new)

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_ss1', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=default_string_length), &
      PARAMETER                              :: zero = ""
    INTEGER, PARAMETER :: t_size = int_size*default_string_length

    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 588 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_ss1

! *****************************************************************************
!> \brief (Re)Allocate a vector of string variables of default path length.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! *****************************************************************************
  SUBROUTINE reallocate_ls1(p,lb1_new,ub1_new)

    CHARACTER(LEN=default_path_length), &
      DIMENSION(:), POINTER                  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_ls1', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=default_path_length), &
      PARAMETER                              :: zero = ""
    INTEGER, PARAMETER :: t_size = int_size*default_path_length

    CHARACTER(LEN=default_path_length), &
      ALLOCATABLE, DIMENSION(:)              :: work


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/reallocate_1.f90" 1
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                ub1,  ub1_old
    IF (ASSOCIATED(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,10,&
                           "work",t_size*(ub1-lb1+1))
       END IF
       work(lb1:ub1) = p(lb1:ub1)
       DEALLOCATE (p)
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,19,&
                        "p",t_size*(ub1_new-lb1_new+1))
    END IF
    p(:) = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work)
    END IF
# 619 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/memory_utilities.F" 2

  END SUBROUTINE reallocate_ls1

END MODULE memory_utilities
