# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_calls.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_calls.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief wrappers for the actual blacs calls.
!>      all functionality needed in the code should actually be provide by cp_blacs_env
!>      these functions should be private members of that module
!> \note
!>      http://www.netlib.org/blacs/BLACS/QRef.html
!> \par History
!>      12.2003 created [Joost]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cp_blacs_calls
  
  USE kinds,                           ONLY: dp

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/../base/base_uses.f90" 1
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
# 20 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_calls.F" 2

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_calls'
  PRIVATE
  ! setup / info calls
  PUBLIC :: cp_blacs_gridinit, cp_blacs_set,&
            cp_blacs_gridexit, cp_blacs_gridinfo
  ! actual message passing
  PUBLIC :: cp_blacs_zgebs2d, cp_blacs_zgebr2d, cp_blacs_dgebs2d,&
            cp_blacs_dgebr2d

!***
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param context ...
!> \param order ...
!> \param nprow ...
!> \param npcol ...
! *****************************************************************************
SUBROUTINE cp_blacs_gridinit(context,order,nprow,npcol)
   INTEGER, INTENT(INOUT) :: context
   CHARACTER(len=1), INTENT(IN):: order
   INTEGER, INTENT(IN)    :: nprow, npcol



   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(context))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(order))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nprow))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(npcol))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE cp_blacs_gridinit

! *****************************************************************************
!> \brief ...
!> \param context ...
! *****************************************************************************
SUBROUTINE cp_blacs_gridexit(context)
   INTEGER, INTENT(IN) :: context



   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(context))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE cp_blacs_gridexit

! *****************************************************************************
!> \brief ...
!> \param context ...
!> \param nprow ...
!> \param npcol ...
!> \param myprow ...
!> \param mypcol ...
! *****************************************************************************
SUBROUTINE cp_blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
   INTEGER, INTENT(IN)  :: context
   INTEGER, INTENT(OUT) :: nprow,npcol,myprow,mypcol



   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(context))==-1) EXIT ;  END DO ; ENDIF
   nprow = 1
   npcol = 1
   myprow = 0
   mypcol = 0

END SUBROUTINE cp_blacs_gridinfo

! *****************************************************************************
!> \brief ...
!> \param context ...
!> \param what :
!>     WHAT = 0 : Handle indicating default system context;  ! DO NOT USE (i.e. use para_env%group)
!>     WHAT = 1 : The BLACS message ID range;
!>     WHAT = 2 : The BLACS debug level the library was compiled with;
!>     WHAT = 10: Handle indicating the system context used to define the BLACS context whose handle is ICONTXT;
!>     WHAT = 11: Number of rings multiring topology is presently using;
!>     WHAT = 12: Number of branches general tree topology is presently using.
!>     WHAT = 15: If non-zero, makes topology choice for repeatable collectives
!> \param val ...
! *****************************************************************************
SUBROUTINE cp_blacs_set(context,what,val)
   INTEGER, INTENT(IN)  :: context,what,val



   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(context))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(what))==-1) EXIT ;  END DO ; ENDIF
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(val))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE cp_blacs_set

! *****************************************************************************
!> \brief ...
!> \param ICONTXT ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
! *****************************************************************************
SUBROUTINE cp_blacs_zgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  COMPLEX(KIND=dp)            :: A



  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ICONTXT))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(SCOPE))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(TOP))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(M))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(N))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(A))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(LDA))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE
! *****************************************************************************
!> \brief ...
!> \param ICONTXT ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
!> \param RSRC ...
!> \param CSRC ...
! *****************************************************************************
SUBROUTINE cp_blacs_zgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  COMPLEX(KIND=dp)            :: A



  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ICONTXT))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(SCOPE))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(TOP))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(M))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(N))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(A))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(LDA))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(RSRC))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(CSRC))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param ICONTXT ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
! *****************************************************************************
SUBROUTINE cp_blacs_dgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  REAL(KIND=dp)               :: A



  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ICONTXT))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(SCOPE))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(TOP))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(M))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(N))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(A))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(LDA))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE
! *****************************************************************************
!> \brief ...
!> \param ICONTXT ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
!> \param RSRC ...
!> \param CSRC ...
! *****************************************************************************
SUBROUTINE cp_blacs_dgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  REAL(KIND=dp)               :: A



  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ICONTXT))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(SCOPE))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(TOP))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(M))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(N))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(A))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(LDA))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(RSRC))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(CSRC))==-1) EXIT ;  END DO ; ENDIF

END SUBROUTINE

END MODULE cp_blacs_calls
