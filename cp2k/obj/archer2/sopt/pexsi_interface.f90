# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the PEXSI library, providing wrappers for all PEXSI 
!>        routines that are called inside CP2K. Tested for PEXSI versions 
!>        v0.7.3 - v0.9.0
!> \par History
!>       2014.12 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
MODULE pexsi_interface

# 25 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
  USE kinds,                           ONLY: int_8,&
                                             real_8

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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pexsi_interface'

  PUBLIC :: cp_pexsi_options, cp_pexsi_plan_initialize, & 
            cp_pexsi_load_real_symmetric_hs_matrix, cp_pexsi_dft_driver, &
            cp_pexsi_retrieve_real_symmetric_dft_matrix, cp_pexsi_plan_finalize,&
            cp_pexsi_set_options, cp_pexsi_get_options, cp_pexsi_set_default_options

  TYPE cp_pexsi_options
      PRIVATE



      INTEGER :: unused = -1

  END TYPE cp_pexsi_options

CONTAINS

! *****************************************************************************
!> \brief Set PEXSI internal options
!> \param pexsi_options ...
!> \param temperature ...
!> \param gap ...
!> \param deltaE ...
!> \param numPole ...
!> \param isInertiaCount ...
!> \param maxPEXSIIter ...
!> \param muMin0 ...
!> \param muMax0 ...
!> \param mu0 ...
!> \param muInertiaTolerance ...
!> \param muInertiaExpansion ...
!> \param muPEXSISafeGuard ...
!> \param numElectronPEXSITolerance ...
!> \param matrixType ...
!> \param isSymbolicFactorize ...
!> \param ordering ...
!> \param npSymbFact ...
!> \param verbosity ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_options(pexsi_options, temperature, gap, deltaE, numPole, & 
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, & 
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity)

    TYPE(cp_pexsi_options), INTENT(INOUT)    :: pexsi_options
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: temperature, gap, deltaE
    INTEGER, INTENT(IN), OPTIONAL            :: numPole, isInertiaCount, &
                                                maxPEXSIIter
    REAL(KIND=real_8), INTENT(IN), OPTIONAL :: muMin0, muMax0, mu0, &
      muInertiaTolerance, muInertiaExpansion, muPEXSISafeGuard, &
      numElectronPEXSITolerance
    INTEGER, INTENT(IN), OPTIONAL            :: matrixType, &
                                                isSymbolicFactorize, &
                                                ordering, npSymbFact, &
                                                verbosity

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_set_options', &
      routineP = moduleN//':'//routineN

# 120 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pexsi_options))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(temperature))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(gap))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(deltaE))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numPole))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(isInertiaCount))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(maxPEXSIIter))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(muMin0))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(muMax0))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(mu0))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(muInertiaTolerance))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(muInertiaExpansion))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(muPEXSISafeGuard))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numElectronPEXSITolerance))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(matrixType))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(isSymbolicFactorize))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(ordering))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(npSymbFact))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(verbosity))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("pexsi_interface.F",139,"Requires linking to the PEXSI library.")

  END SUBROUTINE cp_pexsi_set_options

! *****************************************************************************
!> \brief Access PEXSI internal options
!> \param pexsi_options ...
!> \param temperature ...
!> \param gap ...
!> \param deltaE ...
!> \param numPole ...
!> \param isInertiaCount ...
!> \param maxPEXSIIter ...
!> \param muMin0 ...
!> \param muMax0 ...
!> \param mu0 ...
!> \param muInertiaTolerance ...
!> \param muInertiaExpansion ...
!> \param muPEXSISafeGuard ...
!> \param numElectronPEXSITolerance ...
!> \param matrixType ...
!> \param isSymbolicFactorize ...
!> \param ordering ...
!> \param npSymbFact ...
!> \param verbosity ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_get_options(pexsi_options, temperature, gap, deltaE, numPole, &
                                  isInertiaCount, maxPEXSIIter, muMin0, muMax0, mu0, &
                                  muInertiaTolerance, muInertiaExpansion, &
                                  muPEXSISafeGuard, numElectronPEXSITolerance, &
                                  matrixType, isSymbolicFactorize, ordering, & 
                                  npSymbFact, verbosity)
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: temperature, gap, deltaE
    INTEGER, INTENT(OUT), OPTIONAL           :: numPole, isInertiaCount, &
                                                maxPEXSIIter
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: muMin0, muMax0, mu0, &
      muInertiaTolerance, muInertiaExpansion, muPEXSISafeGuard, &
      numElectronPEXSITolerance
    INTEGER, INTENT(OUT), OPTIONAL           :: matrixType, &
                                                isSymbolicFactorize, &
                                                ordering, npSymbFact, &
                                                verbosity

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_get_options', &
      routineP = moduleN//':'//routineN

# 211 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
   IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pexsi_options))==-1) EXIT ;  END DO ; ENDIF
   ! assign intent-out arguments to silence compiler warnings
   IF(PRESENT(temperature)) temperature = 0.0_real_8
   IF(PRESENT(gap)) gap = 0.0_real_8
   IF(PRESENT(deltaE)) deltaE = 0.0_real_8
   IF(PRESENT(numPole)) numPole = -1
   IF(PRESENT(isInertiaCount)) isInertiaCount = -1
   IF(PRESENT(maxPEXSIIter)) maxPEXSIIter = -1
   IF(PRESENT(muMin0)) muMin0 = 0.0_real_8
   IF(PRESENT(muMax0)) muMax0 = 0.0_real_8
   IF(PRESENT(mu0)) mu0 = 0.0_real_8
   IF(PRESENT(muInertiaTolerance)) muInertiaTolerance = 0.0_real_8
   IF(PRESENT(muInertiaExpansion)) muInertiaExpansion = 0.0_real_8
   IF(PRESENT(muPEXSISafeGuard)) muPEXSISafeGuard = 0.0_real_8
   IF(PRESENT(numElectronPEXSITolerance)) numElectronPEXSITolerance = 0.0_real_8
   IF(PRESENT(matrixType)) matrixType = -1
   IF(PRESENT(isSymbolicFactorize)) isSymbolicFactorize = -1
   IF(PRESENT(ordering)) ordering = -1
   IF(PRESENT(npSymbFact)) npSymbFact = -1
   IF(PRESENT(verbosity)) verbosity = -1
   CALL cp__b("pexsi_interface.F",231,"Requires linking to the PEXSI library.")

  END SUBROUTINE cp_pexsi_get_options

! *****************************************************************************
!> \brief ...
!> \param pexsi_options ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_set_default_options(pexsi_options)
    TYPE(cp_pexsi_options), INTENT(OUT)      :: pexsi_options

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_set_default_options', &
      routineP = moduleN//':'//routineN




    CALL cp__b("pexsi_interface.F",248,"Requires linking to the PEXSI library.")

  END SUBROUTINE cp_pexsi_set_default_options

! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param numProcRow ...
!> \param numProcCol ...
!> \param outputFileIndex ...
!> \retval cp_pexsi_plan_initialize ...
! *****************************************************************************
  FUNCTION cp_pexsi_plan_initialize(comm, numProcRow, numProcCol, outputFileIndex)
    INTEGER, INTENT(IN)                      :: comm, numProcRow, numProcCol, &
                                                outputFileIndex
    INTEGER(KIND=int_8)                      :: cp_pexsi_plan_initialize

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_initialize', &
      routineP = moduleN//':'//routineN

# 277 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(comm))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numProcRow))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numProcCol))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(outputFileIndex))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("pexsi_interface.F",281,"Requires linking to the PEXSI library.")

  END FUNCTION cp_pexsi_plan_initialize

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param pexsi_options ...
!> \param nrows ...
!> \param nnz ...
!> \param nnzLocal ...
!> \param numColLocal ...
!> \param colptrLocal ...
!> \param rowindLocal ...
!> \param HnzvalLocal ...
!> \param isSIdentity ...
!> \param SnzvalLocal ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_load_real_symmetric_hs_matrix(plan,pexsi_options,nrows,nnz, &
                                                    nnzLocal,numColLocal,colptrLocal, &
                                                    rowindLocal,HnzvalLocal,isSIdentity, &
                                                    SnzvalLocal)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    INTEGER, INTENT(IN)                      :: nrows, nnz, nnzLocal, &
                                                numColLocal, colptrLocal(*), &
                                                rowindLocal(*)
    REAL(KIND=real_8), INTENT(IN)            :: HnzvalLocal(*)
    INTEGER, INTENT(IN)                      :: isSIdentity
    REAL(KIND=real_8), INTENT(IN)            :: SnzvalLocal(*)

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_load_real_symmetric_hs_matrix', &
      routineP = moduleN//':'//routineN


# 328 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pexsi_options))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nrows))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nnz))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(nnzLocal))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numColLocal))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(isSIdentity))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("pexsi_interface.F",335,"Requires linking to the PEXSI library.")

    ! MARK_USED macro does not work on assumed shape variables
    IF(.FALSE.) THEN; DO
         IF(colptrLocal(1)>rowindLocal(1).OR.HnzvalLocal(1)>SnzvalLocal(1)) EXIT
    ENDDO; ENDIF

  END SUBROUTINE cp_pexsi_load_real_symmetric_hs_matrix

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param pexsi_options ...
!> \param numElectronExact ...
!> \param muPEXSI ...
!> \param numElectronPEXSI ...
!> \param muMinInertia ...
!> \param muMaxInertia ...
!> \param numTotalInertiaIter ...
!> \param numTotalPEXSIIter ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_dft_driver(plan,pexsi_options,numElectronExact,muPEXSI, &
                                 numElectronPEXSI,muMinInertia,muMaxInertia, &
                                 numTotalInertiaIter,numTotalPEXSIIter)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    TYPE(cp_pexsi_options), INTENT(IN)       :: pexsi_options
    REAL(KIND=real_8), INTENT(IN)            :: numElectronExact
    REAL(KIND=real_8), INTENT(out)           :: muPEXSI, numElectronPEXSI, &
                                                muMinInertia, muMaxInertia
    INTEGER, INTENT(out)                     :: numTotalInertiaIter, &
                                                numTotalPEXSIIter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_dft_driver', &
      routineP = moduleN//':'//routineN

# 381 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(numelectronexact))==-1) EXIT ;  END DO ; ENDIF
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(pexsi_options))==-1) EXIT ;  END DO ; ENDIF
    ! assign intent-out arguments to silence compiler warnings
    muPEXSI = 0.0_real_8
    numElectronPEXSI = 0.0_real_8
    muMinInertia = 0.0_real_8
    muMaxInertia = 0.0_real_8
    numTotalInertiaIter = -1
    numTotalPEXSIIter = -1
    CALL cp__b("pexsi_interface.F",391,"Requires linking to the PEXSI library.")

  END SUBROUTINE cp_pexsi_dft_driver

! *****************************************************************************
!> \brief ...
!> \param plan ...
!> \param DMnzvalLocal ...
!> \param EDMnzvalLocal ...
!> \param FDMnzvalLocal ...
!> \param totalEnergyH ...
!> \param totalEnergyS ...
!> \param totalFreeEnergy ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix(plan,DMnzvalLocal,EDMnzvalLocal, &
                                                         FDMnzvalLocal,totalEnergyH, &
                                                         totalEnergyS,totalFreeEnergy)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan
    REAL(KIND=real_8), INTENT(out) :: DMnzvalLocal(*), EDMnzvalLocal(*), &
      FDMnzvalLocal(*), totalEnergyH, totalEnergyS, totalFreeEnergy

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_pexsi_retrieve_real_symmetric_dft_matrix', &
      routineP = moduleN//':'//routineN


# 428 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
    ! assign intent-out arguments to silence compiler warnings
    DMnzvalLocal(1) = 0.0_real_8
    EDMnzvalLocal(1) = 0.0_real_8
    FDMnzvalLocal(1) = 0.0_real_8
    totalEnergyH = 0.0_real_8
    totalEnergyS = 0.0_real_8
    totalFreeEnergy = 0.0_real_8

    CALL cp__b("pexsi_interface.F",437,"Requires linking to the PEXSI library.")

  END SUBROUTINE cp_pexsi_retrieve_real_symmetric_dft_matrix

! *****************************************************************************
!> \brief ...
!> \param plan ...
! *****************************************************************************
  SUBROUTINE cp_pexsi_plan_finalize(plan)
    INTEGER(KIND=int_8), INTENT(IN)          :: plan

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_pexsi_plan_finalize', &
      routineP = moduleN//':'//routineN


# 460 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pexsi_interface.F"
    IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(plan))==-1) EXIT ;  END DO ; ENDIF
    CALL cp__b("pexsi_interface.F",461,"Requires linking to the PEXSI library.")

  END SUBROUTINE

END MODULE pexsi_interface
