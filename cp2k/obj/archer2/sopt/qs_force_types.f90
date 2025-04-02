# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_force_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_force_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Add CP2K error reporting, new add_force routine [07.2014,JGH]
!> \author MK (03.06.2002)
! *****************************************************************************
MODULE qs_force_types

  !USE cp_control_types,                ONLY: qs_control_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
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
# 18 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/qs_force_types.F" 2

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_force_types'
  PRIVATE

  TYPE qs_force_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: all_potential,&
                                         core_overlap,&
                                         gth_ppl,&
                                         gth_nlcc,&
                                         gth_ppnl,&
                                         kinetic,&
                                         overlap,&
                                         overlap_admm,&
                                         rho_core,&
                                         rho_elec,&
                                         rho_lri_elec,&
                                         vhxc_atom,&
                                         g0s_Vh_elec,&
                                         repulsive,&
                                         dispersion,&
                                         other,&
                                         ch_pulay,&
                                         fock_4c,&
                                         ehrenfest,&
                                         efield,&
                                         eev,&
                                         mp2_sep,&
                                         mp2_non_sep,&
                                         total
  END TYPE qs_force_type

  PUBLIC :: qs_force_type

  PUBLIC :: allocate_qs_force,&
            add_qs_force,&
            deallocate_qs_force,&
            zero_qs_force

CONTAINS

! *****************************************************************************
!> \brief   Allocate a Quickstep force data structure.
!> \param qs_force ...
!> \param natom_of_kind ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_qs_force(qs_force,natom_of_kind)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, n, nkind

    IF (ASSOCIATED(qs_force)) CALL deallocate_qs_force(qs_force)

    nkind = SIZE(natom_of_kind)

    ALLOCATE (qs_force(nkind))

    DO ikind=1,nkind
      n = natom_of_kind(ikind)
      ALLOCATE (qs_force(ikind)%all_potential(3,n))
      ALLOCATE (qs_force(ikind)%core_overlap(3,n))
      ALLOCATE (qs_force(ikind)%gth_ppl(3,n))
      ALLOCATE (qs_force(ikind)%gth_nlcc(3,n))
      ALLOCATE (qs_force(ikind)%gth_ppnl(3,n))
      ALLOCATE (qs_force(ikind)%kinetic(3,n))
      ALLOCATE (qs_force(ikind)%overlap(3,n))
      ALLOCATE (qs_force(ikind)%overlap_admm(3,n))
      ALLOCATE (qs_force(ikind)%rho_core(3,n))
      ALLOCATE (qs_force(ikind)%rho_elec(3,n))
      ALLOCATE (qs_force(ikind)%rho_lri_elec(3,n))
      ALLOCATE (qs_force(ikind)%vhxc_atom(3,n))
      ALLOCATE (qs_force(ikind)%g0s_Vh_elec(3,n))
      ALLOCATE (qs_force(ikind)%repulsive(3,n))
      ALLOCATE (qs_force(ikind)%dispersion(3,n))
      ALLOCATE (qs_force(ikind)%other(3,n))
      ALLOCATE (qs_force(ikind)%ch_pulay(3,n))
      ALLOCATE (qs_force(ikind)%ehrenfest(3,n))
      ALLOCATE (qs_force(ikind)%efield(3,n))
      ALLOCATE (qs_force(ikind)%eev(3,n))
      ! Always initialize ch_pulay to zero..
      qs_force(ikind)%ch_pulay  = 0.0_dp
      ALLOCATE (qs_force(ikind)%fock_4c(3,n))
      ALLOCATE (qs_force(ikind)%mp2_sep(3,n))
      ALLOCATE (qs_force(ikind)%mp2_non_sep(3,n))
      ALLOCATE (qs_force(ikind)%total(3,n))
    END DO

  END SUBROUTINE allocate_qs_force

! *****************************************************************************
!> \brief   Deallocate a Quickstep force data structure.
!> \param qs_force ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_force(qs_force)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind

    IF(.NOT.(ASSOCIATED(qs_force)))CALL cp__a("qs_force_types.F",133)

    nkind = SIZE(qs_force)

    DO ikind=1,nkind

      IF (ASSOCIATED(qs_force(ikind)%all_potential)) THEN
        DEALLOCATE (qs_force(ikind)%all_potential)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DEALLOCATE (qs_force(ikind)%core_overlap)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppl)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_nlcc)) THEN
        DEALLOCATE (qs_force(ikind)%gth_nlcc)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DEALLOCATE (qs_force(ikind)%gth_ppnl)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DEALLOCATE (qs_force(ikind)%kinetic)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DEALLOCATE (qs_force(ikind)%overlap)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%overlap_admm)) THEN
        DEALLOCATE (qs_force(ikind)%overlap_admm)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DEALLOCATE (qs_force(ikind)%rho_core)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_elec)
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_lri_elec)) THEN
        DEALLOCATE (qs_force(ikind)%rho_lri_elec)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%vhxc_atom)) THEN
        DEALLOCATE (qs_force(ikind)%vhxc_atom)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%g0s_Vh_elec)) THEN
        DEALLOCATE (qs_force(ikind)%g0s_Vh_elec)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%repulsive)) THEN
        DEALLOCATE (qs_force(ikind)%repulsive)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%dispersion)) THEN
        DEALLOCATE (qs_force(ikind)%dispersion)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%other)) THEN
        DEALLOCATE (qs_force(ikind)%other)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
        DEALLOCATE (qs_force(ikind)%total)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ch_pulay)) THEN
        DEALLOCATE (qs_force(ikind)%ch_pulay)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%fock_4c)) THEN
        DEALLOCATE (qs_force(ikind)%fock_4c)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_sep)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%mp2_non_sep)) THEN
        DEALLOCATE (qs_force(ikind)%mp2_non_sep)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%ehrenfest)) THEN
        DEALLOCATE (qs_force(ikind)%ehrenfest)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%efield)) THEN
        DEALLOCATE (qs_force(ikind)%efield)
      END IF

      IF (ASSOCIATED(qs_force(ikind)%eev)) THEN
        DEALLOCATE (qs_force(ikind)%eev)
      END IF
    END DO

    DEALLOCATE (qs_force)

  END SUBROUTINE deallocate_qs_force

! *****************************************************************************
!> \brief    Initialize a Quickstep force data structure.
!> \param qs_force ...
!> \date    15.07.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE zero_qs_force(qs_force)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    IF(.NOT.(ASSOCIATED(qs_force)))CALL cp__a("qs_force_types.F",256)

    DO ikind=1,SIZE(qs_force)
      qs_force(ikind)%all_potential(:,:) = 0.0_dp
      qs_force(ikind)%core_overlap(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppl(:,:) = 0.0_dp
      qs_force(ikind)%gth_nlcc(:,:) = 0.0_dp
      qs_force(ikind)%gth_ppnl(:,:) = 0.0_dp
      qs_force(ikind)%kinetic(:,:) = 0.0_dp
      qs_force(ikind)%overlap(:,:) = 0.0_dp
      qs_force(ikind)%overlap_admm(:,:) = 0.0_dp
      qs_force(ikind)%rho_core(:,:) = 0.0_dp
      qs_force(ikind)%rho_elec(:,:) = 0.0_dp
      qs_force(ikind)%rho_lri_elec(:,:) = 0.0_dp
      qs_force(ikind)%vhxc_atom(:,:) = 0.0_dp
      qs_force(ikind)%g0s_Vh_elec(:,:) = 0.0_dp
      qs_force(ikind)%repulsive(:,:) = 0.0_dp
      qs_force(ikind)%dispersion(:,:) = 0.0_dp
      qs_force(ikind)%other(:,:) = 0.0_dp
      qs_force(ikind)%fock_4c(:,:) = 0.0_dp
      qs_force(ikind)%ehrenfest(:,:) = 0.0_dp
      qs_force(ikind)%efield(:,:) = 0.0_dp
      qs_force(ikind)%eev(:,:) = 0.0_dp
      qs_force(ikind)%mp2_non_sep(:,:) = 0.0_dp
      qs_force(ikind)%mp2_sep(:,:) = 0.0_dp
      qs_force(ikind)%total(:,:) = 0.0_dp
    END DO

  END SUBROUTINE zero_qs_force

! *****************************************************************************
!> \brief Add force to a force_type  variable.
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param forcetype ...
!> \param atomic_kind_set ...
!> \par History
!>      07.2014 JGH
!> \author JGH
! *****************************************************************************
  SUBROUTINE add_qs_force(force, qs_force, forcetype, atomic_kind_set)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: force
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    CHARACTER(LEN=*), INTENT(IN)             :: forcetype
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'add_qs_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iatom, ikind, natom_kind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ------------------------------------------------------------------------

    IF(.NOT.(ASSOCIATED(qs_force)))CALL cp__a("qs_force_types.F",314)

    SELECT CASE (forcetype)
       CASE ("overlap_admm")
         DO ikind=1,SIZE(atomic_kind_set,1)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom_kind)
            DO ia=1,natom_kind
              iatom = atomic_kind%atom_list(ia)
              qs_force(ikind)%overlap_admm(:,ia) = qs_force(ikind)%overlap_admm(:,ia) + force(:,iatom)
            END DO
         END DO
       CASE DEFAULT
         CALL cp__b("qs_force_types.F",327,"")
    END SELECT

  END SUBROUTINE add_qs_force

END MODULE qs_force_types
