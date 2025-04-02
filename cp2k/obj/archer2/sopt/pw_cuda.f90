# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on NVIDIA
!> GPUs using CUDA. It depends at execution time on NVIDIA's CUFFT library.
!> \par History
!>      BGL (06-Mar-2008)  : Created
!>      AG  (18-May-2012)  : Refacturing:
!>                           - added explicit interfaces to C routines
!>                           - enable double precision complex transformations
!>      AG  (11-Sept-2012) : Modifications:
!>                          - use pointers if precision mapping is not required
!>                          - use OMP for mapping
!> \author Benjamin G. Levine
! *****************************************************************************
MODULE pw_cuda
  USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                             C_INT,&
                                             C_LOC,&
                                             C_PTR
  USE fast,                            ONLY: zero_c
  USE fft_tools,                       ONLY: &
       cube_transpose_1, cube_transpose_2, fft_scratch_sizes, &
       fft_scratch_type, get_fft_scratch, release_fft_scratch, x_to_yz, &
       xz_to_yz, yz_to_x, yz_to_xz
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE message_passing,                 ONLY: mp_comm_compare,&
                                             mp_environ,&
                                             mp_rank_compare
  USE pw_grid_types,                   ONLY: FULLSPACE
  USE pw_types,                        ONLY: REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE termination,                     ONLY: stop_memory

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/../base/base_uses.f90" 1
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
# 41 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F" 2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pw_cuda_r3dc1d_3d
  PUBLIC :: pw_cuda_c1dr3d_3d
  PUBLIC :: pw_cuda_r3dc1d_3d_ps
  PUBLIC :: pw_cuda_c1dr3d_3d_ps
  PUBLIC :: pw_cuda_init, pw_cuda_finalize

 ! Explicit interfaces to double precision complex transformation
 ! routines (C+CUDA). For more details, see: ../cuda_tools!
  INTERFACE pw_cuda_cfffg_cu
! *****************************************************************************
!> \brief ...
!> \param din ...
!> \param zout ...
!> \param ghatmap ...
!> \param npts ...
!> \param ngpts ...
!> \param scale ...
! *****************************************************************************
    SUBROUTINE pw_cuda_cfffg_z(din, zout, ghatmap, npts, ngpts, scale)&
      BIND(C, name="pw_cuda_cfffg_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: ngpts
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_cfffg_z
  END INTERFACE

  INTERFACE pw_cuda_sfffc_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param dout ...
!> \param ghatmap ...
!> \param npts ...
!> \param ngpts ...
!> \param nmaps ...
!> \param scale ...
! *****************************************************************************
    SUBROUTINE pw_cuda_sfffc_z(zin, dout, ghatmap, npts, ngpts, nmaps, scale)&
      BIND(C, name="pw_cuda_sfffc_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: ngpts, nmaps
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_sfffc_z
  END INTERFACE

  INTERFACE pw_cuda_cff_cu
! *****************************************************************************
!> \brief ...
!> \param din ...
!> \param zout ...
!> \param npts ...
! *****************************************************************************
    SUBROUTINE pw_cuda_cff_z(din, zout, npts)&
      BIND(C, name="pw_cuda_cff_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_cff_z
  END INTERFACE

  INTERFACE pw_cuda_ffc_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param dout ...
!> \param npts ...
! *****************************************************************************
    SUBROUTINE pw_cuda_ffc_z(zin, dout, npts)&
      BIND(C, name="pw_cuda_ffc_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_ffc_z
  END INTERFACE

  INTERFACE pw_cuda_cf_cu
! *****************************************************************************
!> \brief ...
!> \param din ...
!> \param zout ...
!> \param npts ...
! *****************************************************************************
    SUBROUTINE pw_cuda_cf_z(din, zout, npts)&
      BIND(C, name="pw_cuda_cf_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_cf_z
  END INTERFACE

  INTERFACE pw_cuda_fc_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param dout ...
!> \param npts ...
! *****************************************************************************
    SUBROUTINE pw_cuda_fc_z(zin, dout, npts)&
      BIND(C, name="pw_cuda_fc_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_fc_z
  END INTERFACE

  INTERFACE pw_cuda_f_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param zout ...
!> \param dir ...
!> \param n ...
!> \param m ...
! *****************************************************************************
    SUBROUTINE pw_cuda_f_z(zin, zout, dir, n, m)&
      BIND(C, name="pw_cuda_f_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: dir, n, m

    END SUBROUTINE pw_cuda_f_z
  END INTERFACE

  INTERFACE pw_cuda_fg_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param zout ...
!> \param ghatmap ...
!> \param npts ...
!> \param mmax ...
!> \param ngpts ...
!> \param scale ...
! *****************************************************************************
    SUBROUTINE pw_cuda_fg_z(zin, zout, ghatmap, npts, mmax, ngpts, scale)&
      BIND(C, name="pw_cuda_fg_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: mmax, ngpts
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_fg_z
  END INTERFACE

  INTERFACE pw_cuda_sf_cu
! *****************************************************************************
!> \brief ...
!> \param zin ...
!> \param zout ...
!> \param ghatmap ...
!> \param npts ...
!> \param mmax ...
!> \param ngpts ...
!> \param nmaps ...
!> \param scale ...
! *****************************************************************************
    SUBROUTINE pw_cuda_sf_z(zin, zout, ghatmap, npts, mmax, ngpts, nmaps, scale)&
      BIND(C, name="pw_cuda_sf_z_")
      IMPORT
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: mmax, ngpts, nmaps
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_sf_z
  END INTERFACE

  INTERFACE
    FUNCTION pw_cuda_init_cu() RESULT(istat) BIND(C, name="pw_cuda_init")
      IMPORT
      INTEGER(KIND=C_INT)                    :: istat
    END FUNCTION pw_cuda_init_cu
    SUBROUTINE pw_cuda_finalize_cu() BIND(C, name="pw_cuda_finalize")
    END SUBROUTINE pw_cuda_finalize_cu
  END INTERFACE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods_cuda'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

CONTAINS

! *****************************************************************************
!> \brief Allocates resources on the cuda device for cuda fft acceleration
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE pw_cuda_init()






  END SUBROUTINE pw_cuda_init


! *****************************************************************************
!> \brief Releases resources on the cuda device for cuda fft acceleration
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE pw_cuda_finalize()



END SUBROUTINE pw_cuda_finalize


! *****************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \param scale ...
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE pw_cuda_r3dc1d_3d(pw1, pw2, scale)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_r3dc1d_3d', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 334 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_r3dc1d_3d

! *****************************************************************************
!> \brief perform an scatter followed by a fft on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \param scale ...
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE pw_cuda_c1dr3d_3d(pw1, pw2, scale)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_c1dr3d_3d', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 387 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_c1dr3d_3d

! *****************************************************************************
!> \brief perform an parallel fft followed by a gather on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \param scale ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_r3dc1d_3d_ps(pw1, pw2, scale)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_r3dc1d_3d_ps', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 578 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_r3dc1d_3d_ps

! *****************************************************************************
!> \brief perform an parallel scatter followed by a fft on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \param scale ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_c1dr3d_3d_ps(pw1, pw2, scale)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_c1dr3d_3d_ps', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 774 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_c1dr3d_3d_ps

! *****************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 2D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_cff (pw1, pwbuf)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_cff', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
# 816 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_cff

! *****************************************************************************
!> \brief perform a parallel 2D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_ffc (pwbuf, pw2)
    COMPLEX(KIND=dp), DIMENSION(:,:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_ffc', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
# 858 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_ffc

! *****************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 1D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_cf (pw1, pwbuf)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_cf', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
# 900 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_cf

! *****************************************************************************
!> \brief perform a parallel 1D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_fc (pwbuf, pw2)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_fc', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
# 941 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_fc

! *****************************************************************************
!> \brief perform a parallel 1D-FFT on the gpu
!> \param pwbuf1 ...
!> \param pwbuf2 ...
!> \param dir ...
!> \param n ...
!> \param m ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_f(pwbuf1, pwbuf2, dir, n, m)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf2
    INTEGER, INTENT(IN)                      :: dir
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_f', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(dir))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(n))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(m))==-1) EXIT ;  END DO ; ENDIF
# 988 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_f
! *****************************************************************************
!> \brief perform a parallel 1D-FFT followed by a gather on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \param scale ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_fg (pwbuf, pw2, scale)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp), INTENT(IN)                :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_fg', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw2))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 1038 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_fg

! *****************************************************************************
!> \brief perform a parallel scatter followed by a 1D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \param scale ...
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_sf (pw1, pwbuf, scale)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf
    REAL(KIND=dp), INTENT(IN)                :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_sf', &
      routineP = moduleN//':'//routineN


    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pw1))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pwbuf))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(scale))==-1) EXIT ;  END DO ; ENDIF
# 1091 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pw/pw_cuda.F"
  END SUBROUTINE pw_cuda_sf
END MODULE pw_cuda

