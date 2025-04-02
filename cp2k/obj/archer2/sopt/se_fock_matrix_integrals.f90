# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_integrals.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_integrals.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides the low level routines to build both the exchange and the
!>        Coulomb Fock matrices.. This routines support d-orbitals and should
!>        be changed only if one knows exactly what he is doing..
!> \author Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
! *****************************************************************************
MODULE se_fock_matrix_integrals
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_arrays,       ONLY: se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type

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
# 29 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_integrals.F" 2

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_integrals'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: fock2_1el, dfock2_1el, fock1_2el, fock2_1el_ew, fock2C_ew,&
            fock2C, dfock2C, fock2E, dfock2E, fock2_1el_r3, dfock2_1el_r3,&
            fock2C_r3, dfock2C_r3, se_coulomb_ij_interaction

CONTAINS

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param ksi_block DIMENSION(sepi%natorb, sepi%natorb)
!> \param ksj_block DIMENSION(sepi%natorb, sepi%natorb)
!> \param pi_block ...
!> \param pj_block ...
!> \param ecore ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \date   04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE fock2_1el (sepi, sepj, rij, ksi_block, ksj_block, pi_block, pj_block,&
       ecore, itype, anag, se_int_control, se_taper, store_int_env)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksi_block, ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(2), &
      INTENT(INOUT)                          :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

! Compute integrals

    CALL rotnuc (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
          ksi_block(j1,i1) = ksi_block(i1,j1)
          ecore(1) = ecore(1) + 2.0_dp * e1b(i2) * pi_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
       ecore(1) = ecore(1) + e1b(i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
          ksj_block(j1,i1) = ksj_block(i1,j1)
          ecore(2) = ecore(2) + 2.0_dp * e2a(i2) * pj_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
       ecore(2) = ecore(2) + e2a(i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE fock2_1el

! *****************************************************************************
!> \brief Derivatives of 2-center 1-electron Fock Matrix
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param pi_block ...
!> \param pj_block ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param force ...
!> \param delta ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dfock2_1el (sepi, sepj, rij, pi_block, pj_block, itype, anag,&
       se_int_control, se_taper, force, delta)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(3, 45)          :: de1b, de2a

! Compute integrals

    CALL drotnuc (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, delta=delta)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pi_block(i1,j1)
          force(1) = force(1) + de1b(1,i2) * tmp
          force(2) = force(2) + de1b(2,i2) * tmp
          force(3) = force(3) + de1b(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de1b(1,i2) * pi_block(i1,j1)
       force(2) = force(2) + de1b(2,i2) * pi_block(i1,j1)
       force(3) = force(3) + de1b(3,i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pj_block(i1,j1)
          force(1) = force(1) + de2a(1,i2) * tmp
          force(2) = force(2) + de2a(2,i2) * tmp
          force(3) = force(3) + de2a(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de2a(1,i2) * pj_block(i1,j1)
       force(2) = force(2) + de2a(2,i2) * pj_block(i1,j1)
       force(3) = force(3) + de2a(3,i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE dfock2_1el

! *****************************************************************************
!> \brief Construction of 1-center 2-electron Fock Matrix
!> \param sep ...
!> \param p_tot ...
!> \param p_mat ...
!> \param f_mat DIMENSION(sep%natorb, sep%natorb)
!> \param factor ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE fock1_2el(sep, p_tot, p_mat, f_mat, factor)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: p_tot
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(IN)                             :: p_mat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor

    CHARACTER(len=*), PARAMETER :: routineN = 'fock1_2el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ijw, ikw, iL, im, j, jL, &
                                                jlw, jm, k, kL, klw, l, lL
    REAL(KIND=dp)                            :: sum

!   One-center coulomb and exchange terms for semiempirical_type sep
!
!  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
!                        -(PA(k,l)        )*<i,k|j,l>), k,l on type sep.
!

    DO iL = 1, sep%natorb
       i = se_orbital_pointer(iL)
       DO jL = 1, iL
          j = se_orbital_pointer(jL)

          !    `J' Address IJ in W
          ijw = (iL*(iL-1))/2 + jL
          sum = 0.0_dp
          DO kL = 1, sep%natorb
             k = se_orbital_pointer(kL)
             DO lL =  1, sep%natorb
                l = se_orbital_pointer(lL)

                !    `J' Address KL in W
                im = MAX(kL,lL)
                jm = MIN(kL,lL)
                klw = (im*(im-1))/2 + jm

                !    `K' Address IK in W
                im = MAX(kL,jL)
                jm = MIN(kL,jL)
                ikw = (im*(im-1))/2 + jm

                !    `K' Address JL in W
                im = MAX(lL,iL)
                jm = MIN(lL,iL)
                jlw = (im*(im-1))/2 + jm

                sum = sum + p_tot(k,l) * sep%w(ijw, klw) - p_mat(k,l) * sep%w(ikw, jlw)
             END DO
          END DO
          f_mat(i,j) = f_mat(i,j) + factor*sum
          f_mat(j,i) = f_mat(i,j)
       END DO
    END DO
  END SUBROUTINE fock1_2el

! *****************************************************************************
!> \brief Construction of 2-center 1-electron Fock Matrix (Ewald self term)
!> \param sep ...
!> \param rij ...
!> \param ks_block DIMENSION(sep%natorb, sep%natorb)
!> \param p_block ...
!> \param ecore ...
!> \param itype ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \date 04.2009 [jgh]
!> \author jgh - University of Zurich
! *****************************************************************************
  SUBROUTINE fock2_1el_ew (sep, rij, ks_block, p_block, ecore, itype, anag, &
                           se_int_control, se_taper, store_int_env)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(IN)                             :: p_block
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el_ew', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L, n
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

! Compute integrals

    CALL rotnuc (sep, sep, rij, e1b=e1b, e2a=e2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env)
    !
    ! Add the electron-nuclear attraction term for atom sep
    ! e1b == e2a
    !
    n = (sep%natorb*(sep%natorb+1))/2
    i2 = 0
    DO i1L = 1, sep%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ks_block(i1,j1) = ks_block(i1,j1) + e1b(i2)
          ks_block(j1,i1) = ks_block(i1,j1)
          ecore = ecore + 2._dp * e1b(i2) * p_block(i1,j1)
       END DO
       ! i1L == j1L
       i2 = i2 + 1
       ks_block(i1,i1) = ks_block(i1,i1) + e1b(i2)
       ecore = ecore + e1b(i2) * p_block(i1,i1)
    END DO

  END SUBROUTINE fock2_1el_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Self Terms (Ewald)
!> \param sep ...
!> \param rij ...
!> \param p_tot ...
!> \param f_mat DIMENSION(sep%natorb, sep%natorb)
!> \param factor ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \date 04.2009 [jgh]
!> \author jgh - University of Zurich
! *****************************************************************************
  SUBROUTINE fock2C_ew(sep, rij, p_tot, f_mat, factor, anag, se_int_control, &
                       se_taper, store_int_env)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: p_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C_ew', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb
    REAL(KIND=dp)                            :: a, aa, bb
    REAL(KIND=dp), DIMENSION(2025)           :: w

! Evaluate integrals

    CALL rotint (sep,sep,rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env)
    kr = 0
    natorb = sep%natorb
    DO iL = 1, natorb
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = 0.5_dp*w(kr)*factor
                ! Coulomb
                f_mat(i,j) = f_mat(i,j) + bb * a * p_tot(k,l)
                f_mat(k,l) = f_mat(k,l) + aa * a * p_tot(i,j)
                f_mat(j,i) = f_mat(i,j)
                f_mat(l,k) = f_mat(k,l)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2C_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param pi_tot ...
!> \param fi_mat DIMENSION(sepi%natorb, sepi%natorb)
!> \param pj_tot DIMENSION(sepj%natorb, sepj%natorb)
!> \param fj_mat ...
!> \param factor ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE fock2C(sepi, sepj, rij, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, anag, se_int_control, se_taper, store_int_env)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    REAL(KIND=dp)                            :: a, aa, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

! Evaluate integrals

    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env)
    END IF
    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = w(kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   fi_mat(i,j) = fi_mat(i,j) + bb * a * pj_tot(k,l)
                   fj_mat(k,l) = fj_mat(k,l) + aa * a * pi_tot(i,j)
                   fi_mat(j,i) = fi_mat(i,j)
                   fj_mat(l,k) = fj_mat(k,l)
                ELSE
                   fj_mat(i,j) = fj_mat(i,j) + bb * a * pi_tot(k,l)
                   fi_mat(k,l) = fi_mat(k,l) + aa * a * pj_tot(i,j)
                   fj_mat(j,i) = fj_mat(i,j)
                   fi_mat(l,k) = fi_mat(k,l)
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2C

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - Coulomb Terms
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param pi_tot ...
!> \param pj_tot ...
!> \param factor ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param force ...
!> \param delta ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dfock2C(sepi, sepj, rij, switch, pi_tot, pj_tot, factor, anag,&
       se_int_control, se_taper, force, delta)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    REAL(KIND=dp)                            :: aa, bb, tmp
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

! Evaluate integrals' derivatives

    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a(1) = dw(1,kr)*factor
                a(2) = dw(2,kr)*factor
                a(3) = dw(3,kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   tmp = bb * aa * pj_tot(k,l) * pi_tot(i,j)
                ELSE
                   tmp = bb * aa * pi_tot(k,l) * pj_tot(i,j)
                END IF
                force(1) = force(1) + a(1) * tmp
                force(2) = force(2) + a(2) * tmp
                force(3) = force(3) + a(3) * tmp
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2C

! *****************************************************************************
!> \brief Construction of 2-center Fock Matrix - General Driver
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param isize ...
!> \param pi_mat ...
!> \param fi_mat DIMENSION(isize(1), isize(2))
!> \param factor ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param store_int_env ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE fock2E(sepi, sepj, rij, switch, isize, pi_mat, fi_mat, factor,&
       anag, se_int_control, se_taper, store_int_env)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2E', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    REAL(KIND=dp)                            :: a, aa, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

! Evaluate integrals

    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env)
    END IF
    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = w(kr)*factor
                ! Exchange
                a = a * aa * bb * 0.25_dp
                fi_mat(i,k) = fi_mat(i,k) - a * pi_mat(j,l)
                fi_mat(i,l) = fi_mat(i,l) - a * pi_mat(j,k)
                fi_mat(j,k) = fi_mat(j,k) - a * pi_mat(i,l)
                fi_mat(j,l) = fi_mat(j,l) - a * pi_mat(i,k)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2E

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - General Driver
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param isize ...
!> \param pi_mat ...
!> \param factor ...
!> \param anag ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param force ...
!> \param delta ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dfock2E(sepi, sepj, rij, switch, isize,  pi_mat, factor, anag,&
       se_int_control, se_taper, force, delta)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2E', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    REAL(KIND=dp)                            :: aa, bb, tmp, tmp1, tmp2, &
                                                tmp3, tmp4
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

! Evaluate integrals' derivatives

    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr   = kr + 1
                tmp  = factor * aa * bb * 0.25_dp
                a(1) = dw(1,kr)*tmp
                a(2) = dw(2,kr)*tmp
                a(3) = dw(3,kr)*tmp
                ! Exchange
                tmp1 = pi_mat(j,l) * pi_mat(i,k)
                tmp2 = pi_mat(j,k) * pi_mat(i,l)
                tmp3 = pi_mat(i,l) * pi_mat(j,k)
                tmp4 = pi_mat(i,k) * pi_mat(j,l)

                force(1) = force(1) - a(1) * tmp1
                force(1) = force(1) - a(1) * tmp2
                force(1) = force(1) - a(1) * tmp3
                force(1) = force(1) - a(1) * tmp4

                force(2) = force(2) - a(2) * tmp1
                force(2) = force(2) - a(2) * tmp2
                force(2) = force(2) - a(2) * tmp3
                force(2) = force(2) - a(2) * tmp4

                force(3) = force(3) - a(3) * tmp1
                force(3) = force(3) - a(3) * tmp2
                force(3) = force(3) - a(3) * tmp3
                force(3) = force(3) - a(3) * tmp4
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2E

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix for the residual
!>         (1/R^3) integral part
!> \param sepi ...
!> \param sepj ...
!> \param ksi_block DIMENSION(sepi%natorb, sepi%natorb)
!> \param ksj_block DIMENSION(sepj%natorb, sepj%natorb)
!> \param pi_block ...
!> \param pj_block ...
!> \param e1b ...
!> \param e2a ...
!> \param ecore ...
!> \param rp ...
!> \date   12.2008 [tlaino]
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE fock2_1el_r3 (sepi, sepj, ksi_block, ksj_block, pi_block, pj_block,&
       e1b, e2a, ecore, rp)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksi_block, ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a
    REAL(KIND=dp), DIMENSION(2), &
      INTENT(INOUT)                          :: ecore
    REAL(KIND=dp), INTENT(IN)                :: rp

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2

!
! Add the electron-nuclear residual attraction term for atom sepi
!

    i2 = 0
    DO i1L = 1, sepi%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksi_block(i1,i1) = ksi_block(i1,i1) + e1b(i2) * rp
       ecore(1) = ecore(1) + e1b(i2) * rp * pi_block(i1,i1)
    END DO
    !
    ! Add the electron-nuclear residual attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksj_block(i1,i1) = ksj_block(i1,i1) + e2a(i2) * rp
       ecore(2) = ecore(2) + e2a(i2) * rp * pj_block(i1,i1)
    END DO

  END SUBROUTINE fock2_1el_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center 1-electron Fock Matrix residual (1/R^3)
!>         integral part
!> \param sepi ...
!> \param sepj ...
!> \param drp ...
!> \param pi_block ...
!> \param pj_block ...
!> \param force ...
!> \param e1b ...
!> \param e2a ...
!> \date   12.2008 [tlaino]
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2_1el_r3 (sepi, sepj, drp, pi_block, pj_block, force, e1b, e2a)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2
    REAL(KIND=dp)                            :: tmp

!
! Add the electron-nuclear residual attraction term for atom sepi
!

    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e1b(i2) * pi_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e2a(i2) * pj_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO

  END SUBROUTINE dfock2_1el_r3

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part
!> \param sepi ...
!> \param sepj ...
!> \param switch ...
!> \param pi_tot ...
!> \param fi_mat DIMENSION(sepi%natorb, sepi%natorb)
!> \param pj_tot ...
!> \param fj_mat DIMENSION(sepj%natorb, sepj%natorb)
!> \param factor ...
!> \param w ...
!> \param rp ...
!> \date   12.2008 [tlaino]
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE fock2C_r3(sepi, sepj, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, w, rp)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), INTENT(IN)                :: rp

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    REAL(KIND=dp)                            :: a, w_l(81)

    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          a = w_l(kr) * factor * rp
          ! Coulomb
          IF (.NOT.switch) THEN
             fi_mat(i,i) = fi_mat(i,i) + a * pj_tot(k,k)
             fj_mat(k,k) = fj_mat(k,k) + a * pi_tot(i,i)
          ELSE
             fj_mat(i,i) = fj_mat(i,i) + a * pi_tot(k,k)
             fi_mat(k,k) = fi_mat(k,k) + a * pj_tot(i,i)
          END IF
       END DO
    END DO

  END SUBROUTINE fock2C_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part
!> \param sepi ...
!> \param sepj ...
!> \param switch ...
!> \param pi_tot ...
!> \param pj_tot ...
!> \param factor ...
!> \param w ...
!> \param drp ...
!> \param force ...
!> \date   12.2008 [tlaino]
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2C_r3(sepi, sepj, switch, pi_tot, pj_tot, factor,  w, drp,&
       force)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    REAL(KIND=dp)                            :: a(3), tmp, w_l(81)

    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          tmp  = w_l(kr) * factor
          a(1) = tmp * drp(1)
          a(2) = tmp * drp(2)
          a(3) = tmp * drp(3)
          ! Coulomb
          IF (.NOT.switch) THEN
             tmp = pj_tot(k,k) * pi_tot(i,i)
          ELSE
             tmp = pi_tot(k,k) * pj_tot(i,i)
          END IF
          force(1) = force(1) + a(1) * tmp
          force(2) = force(2) + a(2) * tmp
          force(3) = force(3) + a(3) * tmp
       END DO
    END DO

  END SUBROUTINE dfock2C_r3

! *****************************************************************************
!> \brief  Coulomb interaction multipolar correction
!> \param atom_a ...
!> \param atom_b ...
!> \param my_task ...
!> \param do_forces ...
!> \param do_efield ...
!> \param do_stress ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param force_ab ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \param rab2 ...
!> \param rab ...
!> \param integral_value ...
!> \param ptens11 ...
!> \param ptens12 ...
!> \param ptens13 ...
!> \param ptens21 ...
!> \param ptens22 ...
!> \param ptens23 ...
!> \param ptens31 ...
!> \param ptens32 ...
!> \param ptens33 ...
!> \date   05.2009 [tlaino]
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE se_coulomb_ij_interaction (atom_a, atom_b, my_task, do_forces, do_efield,&
       do_stress, charges, dipoles, quadrupoles, force_ab, efield0, efield1, efield2, &
       rab2, rab, integral_value, ptens11, ptens12, ptens13, ptens21, ptens22, ptens23,&
       ptens31, ptens32, ptens33)
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    LOGICAL, DIMENSION(3)                    :: my_task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupoles
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force_ab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1, efield2
    REAL(KIND=dp), INTENT(IN)                :: rab2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: integral_value
    REAL(KIND=dp), INTENT(INOUT)             :: ptens11, ptens12, ptens13, &
                                                ptens21, ptens22, ptens23, &
                                                ptens31, ptens32, ptens33

    CHARACTER(len=*), PARAMETER :: routineN = 'se_coulomb_ij_interaction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, c, d, e, i, j, k
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, force_eval
    LOGICAL, DIMENSION(3)                    :: do_task
    LOGICAL, DIMENSION(3, 3)                 :: task
    REAL(KIND=dp) :: ch_i, ch_j, ef0_i, ef0_j, eloc, energy, fac, fac_ij, ir, &
      irab2, r, tij, tmp, tmp1, tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, &
      tmp23, tmp31, tmp32, tmp33, tmp_ij, tmp_ji
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: dp_i, dp_j, ef1_i, ef1_j, fr, &
                                                tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: ef2_i, ef2_j, qp_i, qp_j, &
                                                tij_ab
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: tij_abc
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3)     :: tij_abcd
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)  :: tij_abcde

    do_task  = my_task
    energy   = 0.0_dp
    DO i = 1, 3
       IF (do_task(i)) THEN
          SELECT CASE(i)
          CASE(1)
             do_task(1) = (charges(atom_a)/=0.0_dp).OR.(charges(atom_b)/=0.0_dp)
          CASE(2)
             do_task(2) = (ANY(dipoles(:,atom_a)/=0.0_dp)).OR.(ANY(dipoles(:,atom_b)/=0.0_dp))
          CASE(3)
             do_task(3) = (ANY(quadrupoles(:,:,atom_a)/=0.0_dp)).OR.(ANY(quadrupoles(:,:,atom_b)/=0.0_dp))
          END SELECT
       END IF
    END DO
    DO i = 1,3
       DO j =i,3
          task(j,i) = do_task(i).AND.do_task(j)
          task(i,j) = task(j,i)
       END DO
    END DO
    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)

    fac_ij = 1.0_dp
    IF (atom_a==atom_b) fac_ij = 0.5_dp


# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90" 1
                ! Compute the Short Range constribution according the task
                IF (debug_this_module) THEN
                   f         = HUGE(0.0_dp)
                   tij       = HUGE(0.0_dp)
                   tij_a     = HUGE(0.0_dp)
                   tij_ab    = HUGE(0.0_dp)
                   tij_abc   = HUGE(0.0_dp)
                   tij_abcd  = HUGE(0.0_dp)
                   tij_abcde = HUGE(0.0_dp)
                END IF
                r     = SQRT(rab2)
                irab2 = 1.0_dp/rab2
                ir    = 1.0_dp/r

                ! Compute the radial function
# 63 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"
                ! code for point multipole without screening
                f(0)  = ir
                DO i = 1, 5
                   f(i) = irab2*f(i-1)
                END DO
# 85 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"
                ! Compute the Tensor components
                force_eval = do_stress
                IF (task(1,1)) THEN
                   tij         = f(0)*fac_ij
                                                 force_eval = do_forces .OR.do_efield1
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield0
                IF (task(1,2).OR.force_eval) THEN
                   force_eval = do_stress
                   tij_a    = - rab*f(1)*fac_ij
                   IF (task(1,2))                force_eval = force_eval.OR. do_forces
                END IF
                IF (task(1,1))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield0
                IF (task(2,2).OR.task(3,1).OR.force_eval) THEN
                   force_eval = do_stress
                   DO b = 1,3
                      DO a = 1,3
                         tmp = rab(a)*rab(b)*fac_ij
                         tij_ab(a,b) = 3.0_dp*tmp*f(2)
                         IF (a==b) tij_ab(a,b) = tij_ab(a,b) - f(1)*fac_ij
                      END DO
                   END DO
                   IF (task(2,2).OR.task(3,1))   force_eval = force_eval.OR. do_forces
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield1
                IF (task(3,2).OR.force_eval) THEN
                   force_eval = do_stress
                   DO c = 1, 3
                      DO b = 1, 3
                         DO a = 1, 3
                            tmp = rab(a)*rab(b)*rab(c)*fac_ij
                            tij_abc(a,b,c) = - 15.0_dp*tmp*f(3)
                            tmp = 3.0_dp*f(2)*fac_ij
                            IF (a==b) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(c)
                            IF (a==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(b)
                            IF (b==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(a)
                         END DO
                      END DO
                   END DO
                   IF (task(3,2))                force_eval = force_eval.OR. do_forces
                END IF
                IF (task(3,3).OR.force_eval) THEN
                   force_eval = do_stress
                   DO d = 1, 3
                      DO c = 1, 3
                         DO b = 1, 3
                            DO a = 1, 3
                               tmp = rab(a)*rab(b)*rab(c)*rab(d)*fac_ij
                               tij_abcd(a,b,c,d) = 105.0_dp*tmp*f(4)
                               tmp1 = 15.0_dp*f(3)*fac_ij
                               tmp2 =  3.0_dp*f(2)*fac_ij
                               IF (a==b) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(c)*rab(d)
                                  IF (c==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==c) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(d)
                                  IF (b==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(c)
                               IF (b==c) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(d)
                                  IF (a==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (b==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(c)
                               IF (c==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(b)
                            END DO
                         END DO
                      END DO
                   END DO
                   IF (task(3,3))                force_eval = force_eval.OR. do_forces
                END IF
                IF (force_eval) THEN
                   force_eval = do_stress
                   DO e = 1, 3
                      DO d = 1, 3
                         DO c = 1, 3
                            DO b = 1, 3
                               DO a = 1, 3
                                  tmp = rab(a)*rab(b)*rab(c)*rab(d)*rab(e)*fac_ij
                                  tij_abcde(a,b,c,d,e) = -945.0_dp*tmp*f(5)
                                  tmp1 = 105.0_dp*f(4)*fac_ij
                                  tmp2 =  15.0_dp*f(3)*fac_ij
                                  IF (a==b) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(c)*rab(d)*rab(e)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                  END IF
                                  IF (a==c) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(d)*rab(e)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (a==d) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(e)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (a==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(d)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (b==c) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(d)*rab(e)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==d) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(d)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (c==d)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(e)
                                  IF (c==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(d)
                                  IF (d==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(c)
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END IF
                eloc  = 0.0_dp
                fr    = 0.0_dp
                ef0_i = 0.0_dp
                ef0_j = 0.0_dp
                ef1_j = 0.0_dp
                ef1_i = 0.0_dp
                ef2_j = 0.0_dp
                ef2_i = 0.0_dp

# 325 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"

                ! Initialize the charge, dipole and quadrupole for atom A and B
                IF (debug_this_module) THEN
                   ch_j  = HUGE(0.0_dp)
                   ch_i  = HUGE(0.0_dp)
                   dp_j  = HUGE(0.0_dp)
                   dp_i  = HUGE(0.0_dp)
                   qp_j  = HUGE(0.0_dp)
                   qp_i  = HUGE(0.0_dp)
                END IF
                IF (ANY(task(1,:))) THEN
                   ch_j  = charges(atom_a)
                   ch_i  = charges(atom_b)
                END IF
                IF (ANY(task(2,:))) THEN
                   dp_j  = dipoles(:,atom_a)
                   dp_i  = dipoles(:,atom_b)
                END IF
                IF (ANY(task(3,:))) THEN
                   qp_j  = quadrupoles(:,:,atom_a)
                   qp_i  = quadrupoles(:,:,atom_b)
                END IF
                IF (task(1,1)) THEN
                   ! Charge - Charge
                   eloc = eloc + ch_i*tij*ch_j
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      fr(1) = fr(1) - ch_j * tij_a(1) * ch_i
                      fr(2) = fr(2) - ch_j * tij_a(2) * ch_i
                      fr(3) = fr(3) - ch_j * tij_a(3) * ch_i
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + tij * ch_j

                         ef0_j = ef0_j + tij * ch_i
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - tij_a(1) * ch_j
                         ef1_i(2) = ef1_i(2) - tij_a(2) * ch_j
                         ef1_i(3) = ef1_i(3) - tij_a(3) * ch_j

                         ef1_j(1) = ef1_j(1) + tij_a(1) * ch_i
                         ef1_j(2) = ef1_j(2) + tij_a(2) * ch_i
                         ef1_j(3) = ef1_j(3) + tij_a(3) * ch_i

# 383 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"

                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) - tij_ab(1,1) * ch_j
                         ef2_i(2,1) = ef2_i(2,1) - tij_ab(2,1) * ch_j
                         ef2_i(3,1) = ef2_i(3,1) - tij_ab(3,1) * ch_j
                         ef2_i(1,2) = ef2_i(1,2) - tij_ab(1,2) * ch_j
                         ef2_i(2,2) = ef2_i(2,2) - tij_ab(2,2) * ch_j
                         ef2_i(3,2) = ef2_i(3,2) - tij_ab(3,2) * ch_j
                         ef2_i(1,3) = ef2_i(1,3) - tij_ab(1,3) * ch_j
                         ef2_i(2,3) = ef2_i(2,3) - tij_ab(2,3) * ch_j
                         ef2_i(3,3) = ef2_i(3,3) - tij_ab(3,3) * ch_j

                         ef2_j(1,1) = ef2_j(1,1) - tij_ab(1,1) * ch_i
                         ef2_j(2,1) = ef2_j(2,1) - tij_ab(2,1) * ch_i
                         ef2_j(3,1) = ef2_j(3,1) - tij_ab(3,1) * ch_i
                         ef2_j(1,2) = ef2_j(1,2) - tij_ab(1,2) * ch_i
                         ef2_j(2,2) = ef2_j(2,2) - tij_ab(2,2) * ch_i
                         ef2_j(3,2) = ef2_j(3,2) - tij_ab(3,2) * ch_i
                         ef2_j(1,3) = ef2_j(1,3) - tij_ab(1,3) * ch_i
                         ef2_j(2,3) = ef2_j(2,3) - tij_ab(2,3) * ch_i
                         ef2_j(3,3) = ef2_j(3,3) - tij_ab(3,3) * ch_i
                      END IF
                   END IF
                END IF
                IF (task(2,2)) THEN
                   ! Dipole - Dipole
                   tmp= - (dp_i(1)*(tij_ab(1,1)*dp_j(1)+&
                                    tij_ab(2,1)*dp_j(2)+&
                                    tij_ab(3,1)*dp_j(3))+&
                           dp_i(2)*(tij_ab(1,2)*dp_j(1)+&
                                    tij_ab(2,2)*dp_j(2)+&
                                    tij_ab(3,2)*dp_j(3))+&
                           dp_i(3)*(tij_ab(1,3)*dp_j(1)+&
                                    tij_ab(2,3)*dp_j(2)+&
                                    tij_ab(3,3)*dp_j(3)))
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) +  dp_i(1)*(tij_abc(1,1,k)*dp_j(1)+&
                                                   tij_abc(2,1,k)*dp_j(2)+&
                                                   tij_abc(3,1,k)*dp_j(3))&
                                       +  dp_i(2)*(tij_abc(1,2,k)*dp_j(1)+&
                                                   tij_abc(2,2,k)*dp_j(2)+&
                                                   tij_abc(3,2,k)*dp_j(3))&
                                       +  dp_i(3)*(tij_abc(1,3,k)*dp_j(1)+&
                                                   tij_abc(2,3,k)*dp_j(2)+&
                                                   tij_abc(3,3,k)*dp_j(3))
                      END DO
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i - (tij_a(1)*dp_j(1)+&
                                          tij_a(2)*dp_j(2)+&
                                          tij_a(3)*dp_j(3))

                         ef0_j = ef0_j + (tij_a(1)*dp_i(1)+&
                                          tij_a(2)*dp_i(2)+&
                                          tij_a(3)*dp_i(3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) + (tij_ab(1,1)*dp_j(1)+&
                                                tij_ab(2,1)*dp_j(2)+&
                                                tij_ab(3,1)*dp_j(3))
                         ef1_i(2) = ef1_i(2) + (tij_ab(1,2)*dp_j(1)+&
                                                tij_ab(2,2)*dp_j(2)+&
                                                tij_ab(3,2)*dp_j(3))
                         ef1_i(3) = ef1_i(3) + (tij_ab(1,3)*dp_j(1)+&
                                                tij_ab(2,3)*dp_j(2)+&
                                                tij_ab(3,3)*dp_j(3))

                         ef1_j(1) = ef1_j(1) + (tij_ab(1,1)*dp_i(1)+&
                                                tij_ab(2,1)*dp_i(2)+&
                                                tij_ab(3,1)*dp_i(3))
                         ef1_j(2) = ef1_j(2) + (tij_ab(1,2)*dp_i(1)+&
                                                tij_ab(2,2)*dp_i(2)+&
                                                tij_ab(3,2)*dp_i(3))
                         ef1_j(3) = ef1_j(3) + (tij_ab(1,3)*dp_i(1)+&
                                                tij_ab(2,3)*dp_i(2)+&
                                                tij_ab(3,3)*dp_i(3))
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) + (tij_abc(1,1,1)*dp_j(1)+&
                                                    tij_abc(2,1,1)*dp_j(2)+&
                                                    tij_abc(3,1,1)*dp_j(3))
                         ef2_i(1,2) = ef2_i(1,2) + (tij_abc(1,1,2)*dp_j(1)+&
                                                    tij_abc(2,1,2)*dp_j(2)+&
                                                    tij_abc(3,1,2)*dp_j(3))
                         ef2_i(1,3) = ef2_i(1,3) + (tij_abc(1,1,3)*dp_j(1)+&
                                                    tij_abc(2,1,3)*dp_j(2)+&
                                                    tij_abc(3,1,3)*dp_j(3))
                         ef2_i(2,1) = ef2_i(2,1) + (tij_abc(1,2,1)*dp_j(1)+&
                                                    tij_abc(2,2,1)*dp_j(2)+&
                                                    tij_abc(3,2,1)*dp_j(3))
                         ef2_i(2,2) = ef2_i(2,2) + (tij_abc(1,2,2)*dp_j(1)+&
                                                    tij_abc(2,2,2)*dp_j(2)+&
                                                    tij_abc(3,2,2)*dp_j(3))
                         ef2_i(2,3) = ef2_i(2,3) + (tij_abc(1,2,3)*dp_j(1)+&
                                                    tij_abc(2,2,3)*dp_j(2)+&
                                                    tij_abc(3,2,3)*dp_j(3))
                         ef2_i(3,1) = ef2_i(3,1) + (tij_abc(1,3,1)*dp_j(1)+&
                                                    tij_abc(2,3,1)*dp_j(2)+&
                                                    tij_abc(3,3,1)*dp_j(3))
                         ef2_i(3,2) = ef2_i(3,2) + (tij_abc(1,3,2)*dp_j(1)+&
                                                    tij_abc(2,3,2)*dp_j(2)+&
                                                    tij_abc(3,3,2)*dp_j(3))
                         ef2_i(3,3) = ef2_i(3,3) + (tij_abc(1,3,3)*dp_j(1)+&
                                                    tij_abc(2,3,3)*dp_j(2)+&
                                                    tij_abc(3,3,3)*dp_j(3))

                         ef2_j(1,1) = ef2_j(1,1) - (tij_abc(1,1,1)*dp_i(1)+&
                                                    tij_abc(2,1,1)*dp_i(2)+&
                                                    tij_abc(3,1,1)*dp_i(3))
                         ef2_j(1,2) = ef2_j(1,2) - (tij_abc(1,1,2)*dp_i(1)+&
                                                    tij_abc(2,1,2)*dp_i(2)+&
                                                    tij_abc(3,1,2)*dp_i(3))
                         ef2_j(1,3) = ef2_j(1,3) - (tij_abc(1,1,3)*dp_i(1)+&
                                                    tij_abc(2,1,3)*dp_i(2)+&
                                                    tij_abc(3,1,3)*dp_i(3))
                         ef2_j(2,1) = ef2_j(2,1) - (tij_abc(1,2,1)*dp_i(1)+&
                                                    tij_abc(2,2,1)*dp_i(2)+&
                                                    tij_abc(3,2,1)*dp_i(3))
                         ef2_j(2,2) = ef2_j(2,2) - (tij_abc(1,2,2)*dp_i(1)+&
                                                    tij_abc(2,2,2)*dp_i(2)+&
                                                    tij_abc(3,2,2)*dp_i(3))
                         ef2_j(2,3) = ef2_j(2,3) - (tij_abc(1,2,3)*dp_i(1)+&
                                                    tij_abc(2,2,3)*dp_i(2)+&
                                                    tij_abc(3,2,3)*dp_i(3))
                         ef2_j(3,1) = ef2_j(3,1) - (tij_abc(1,3,1)*dp_i(1)+&
                                                    tij_abc(2,3,1)*dp_i(2)+&
                                                    tij_abc(3,3,1)*dp_i(3))
                         ef2_j(3,2) = ef2_j(3,2) - (tij_abc(1,3,2)*dp_i(1)+&
                                                    tij_abc(2,3,2)*dp_i(2)+&
                                                    tij_abc(3,3,2)*dp_i(3))
                         ef2_j(3,3) = ef2_j(3,3) - (tij_abc(1,3,3)*dp_i(1)+&
                                                    tij_abc(2,3,3)*dp_i(2)+&
                                                    tij_abc(3,3,3)*dp_i(3))
                      END IF
                   END IF
                END IF
                IF (task(2,1)) THEN
                   ! Dipole - Charge
                   tmp=   ch_j*(tij_a(1)*dp_i(1)+&
                                tij_a(2)*dp_i(2)+&
                                tij_a(3)*dp_i(3))&
                        - ch_i*(tij_a(1)*dp_j(1)+&
                                tij_a(2)*dp_j(2)+&
                                tij_a(3)*dp_j(3))
# 545 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) -  ch_j *(tij_ab(1,k)*dp_i(1)+&
                                                 tij_ab(2,k)*dp_i(2)+&
                                                 tij_ab(3,k)*dp_i(3))&
                                       +  ch_i *(tij_ab(1,k)*dp_j(1)+&
                                                 tij_ab(2,k)*dp_j(2)+&
                                                 tij_ab(3,k)*dp_j(3))
# 563 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"
                      END DO
                   END IF
                END IF
                IF (task(3,3)) THEN
                   ! Quadrupole - Quadrupole
                   fac  = 1.0_dp/9.0_dp
                   tmp11 = qp_i(1,1)*(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,1)*qp_j(3,3))
                   tmp21 = qp_i(2,1)*(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,2)*qp_j(3,3))
                   tmp31 = qp_i(3,1)*(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,3)*qp_j(3,3))
                   tmp22 = qp_i(2,2)*(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,2)*qp_j(3,3))
                   tmp32 = qp_i(3,2)*(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,3)*qp_j(3,3))
                   tmp33 = qp_i(3,3)*(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,3,3)*qp_j(3,3))
                   tmp12 = tmp21
                   tmp13 = tmp31
                   tmp23 = tmp32
                   tmp   = tmp11 + tmp12 + tmp13 + &
                           tmp21 + tmp22 + tmp23 + &
                           tmp31 + tmp32 + tmp33

                   eloc = eloc + fac*tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         tmp11 = qp_i(1,1)*(tij_abcde(1,1,1,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,1,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,1,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,1,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,1,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,1,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,1,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,1,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,1,1,k)*qp_j(3,3))
                         tmp21 = qp_i(2,1)*(tij_abcde(1,1,2,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,1,k)*qp_j(3,3))
                         tmp31 = qp_i(3,1)*(tij_abcde(1,1,3,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,1,k)*qp_j(3,3))
                         tmp22 = qp_i(2,2)*(tij_abcde(1,1,2,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,2,k)*qp_j(3,3))
                         tmp32 = qp_i(3,2)*(tij_abcde(1,1,3,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,2,k)*qp_j(3,3))
                         tmp33 = qp_i(3,3)*(tij_abcde(1,1,3,3,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,3,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,3,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,3,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,3,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,3,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,3,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,3,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,3,k)*qp_j(3,3))
                         tmp12 = tmp21
                         tmp13 = tmp31
                         tmp23 = tmp32
                         fr(k) = fr(k) - fac * ( tmp11 + tmp12 + tmp13 +&
                                                 tmp21 + tmp22 + tmp23 +&
                                                 tmp31 + tmp32 + tmp33  )
                      END DO
                   END IF
                   ! Electric field
                   IF (do_efield) THEN
                      fac = 1.0_dp/3.0_dp
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + fac*(tij_ab(1,1)*qp_j(1,1)+&
                                              tij_ab(2,1)*qp_j(2,1)+&
                                              tij_ab(3,1)*qp_j(3,1)+&
                                              tij_ab(1,2)*qp_j(1,2)+&
                                              tij_ab(2,2)*qp_j(2,2)+&
                                              tij_ab(3,2)*qp_j(3,2)+&
                                              tij_ab(1,3)*qp_j(1,3)+&
                                              tij_ab(2,3)*qp_j(2,3)+&
                                              tij_ab(3,3)*qp_j(3,3))

                         ef0_j = ef0_j + fac*(tij_ab(1,1)*qp_i(1,1)+&
                                              tij_ab(2,1)*qp_i(2,1)+&
                                              tij_ab(3,1)*qp_i(3,1)+&
                                              tij_ab(1,2)*qp_i(1,2)+&
                                              tij_ab(2,2)*qp_i(2,2)+&
                                              tij_ab(3,2)*qp_i(3,2)+&
                                              tij_ab(1,3)*qp_i(1,3)+&
                                              tij_ab(2,3)*qp_i(2,3)+&
                                              tij_ab(3,3)*qp_i(3,3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - fac*(tij_abc(1,1,1)*qp_j(1,1)+&
                                                    tij_abc(2,1,1)*qp_j(2,1)+&
                                                    tij_abc(3,1,1)*qp_j(3,1)+&
                                                    tij_abc(1,2,1)*qp_j(1,2)+&
                                                    tij_abc(2,2,1)*qp_j(2,2)+&
                                                    tij_abc(3,2,1)*qp_j(3,2)+&
                                                    tij_abc(1,3,1)*qp_j(1,3)+&
                                                    tij_abc(2,3,1)*qp_j(2,3)+&
                                                    tij_abc(3,3,1)*qp_j(3,3))
                         ef1_i(2) = ef1_i(2) - fac*(tij_abc(1,1,2)*qp_j(1,1)+&
                                                    tij_abc(2,1,2)*qp_j(2,1)+&
                                                    tij_abc(3,1,2)*qp_j(3,1)+&
                                                    tij_abc(1,2,2)*qp_j(1,2)+&
                                                    tij_abc(2,2,2)*qp_j(2,2)+&
                                                    tij_abc(3,2,2)*qp_j(3,2)+&
                                                    tij_abc(1,3,2)*qp_j(1,3)+&
                                                    tij_abc(2,3,2)*qp_j(2,3)+&
                                                    tij_abc(3,3,2)*qp_j(3,3))
                         ef1_i(3) = ef1_i(3) - fac*(tij_abc(1,1,3)*qp_j(1,1)+&
                                                    tij_abc(2,1,3)*qp_j(2,1)+&
                                                    tij_abc(3,1,3)*qp_j(3,1)+&
                                                    tij_abc(1,2,3)*qp_j(1,2)+&
                                                    tij_abc(2,2,3)*qp_j(2,2)+&
                                                    tij_abc(3,2,3)*qp_j(3,2)+&
                                                    tij_abc(1,3,3)*qp_j(1,3)+&
                                                    tij_abc(2,3,3)*qp_j(2,3)+&
                                                    tij_abc(3,3,3)*qp_j(3,3))

                         ef1_j(1) = ef1_j(1) + fac*(tij_abc(1,1,1)*qp_i(1,1)+&
                                                    tij_abc(2,1,1)*qp_i(2,1)+&
                                                    tij_abc(3,1,1)*qp_i(3,1)+&
                                                    tij_abc(1,2,1)*qp_i(1,2)+&
                                                    tij_abc(2,2,1)*qp_i(2,2)+&
                                                    tij_abc(3,2,1)*qp_i(3,2)+&
                                                    tij_abc(1,3,1)*qp_i(1,3)+&
                                                    tij_abc(2,3,1)*qp_i(2,3)+&
                                                    tij_abc(3,3,1)*qp_i(3,3))
                         ef1_j(2) = ef1_j(2) + fac*(tij_abc(1,1,2)*qp_i(1,1)+&
                                                    tij_abc(2,1,2)*qp_i(2,1)+&
                                                    tij_abc(3,1,2)*qp_i(3,1)+&
                                                    tij_abc(1,2,2)*qp_i(1,2)+&
                                                    tij_abc(2,2,2)*qp_i(2,2)+&
                                                    tij_abc(3,2,2)*qp_i(3,2)+&
                                                    tij_abc(1,3,2)*qp_i(1,3)+&
                                                    tij_abc(2,3,2)*qp_i(2,3)+&
                                                    tij_abc(3,3,2)*qp_i(3,3))
                         ef1_j(3) = ef1_j(3) + fac*(tij_abc(1,1,3)*qp_i(1,1)+&
                                                    tij_abc(2,1,3)*qp_i(2,1)+&
                                                    tij_abc(3,1,3)*qp_i(3,1)+&
                                                    tij_abc(1,2,3)*qp_i(1,2)+&
                                                    tij_abc(2,2,3)*qp_i(2,2)+&
                                                    tij_abc(3,2,3)*qp_i(3,2)+&
                                                    tij_abc(1,3,3)*qp_i(1,3)+&
                                                    tij_abc(2,3,3)*qp_i(2,3)+&
                                                    tij_abc(3,3,3)*qp_i(3,3))
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_j(3,3))
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_j(3,3))
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_j(3,3))
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_j(3,3))
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_j(3,3))
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_j(3,3))

                         ef2_i(1,1) = ef2_i(1,1) - tmp11
                         ef2_i(1,2) = ef2_i(1,2) - tmp12
                         ef2_i(1,3) = ef2_i(1,3) - tmp13
                         ef2_i(2,1) = ef2_i(2,1) - tmp12
                         ef2_i(2,2) = ef2_i(2,2) - tmp22
                         ef2_i(2,3) = ef2_i(2,3) - tmp23
                         ef2_i(3,1) = ef2_i(3,1) - tmp13
                         ef2_i(3,2) = ef2_i(3,2) - tmp23
                         ef2_i(3,3) = ef2_i(3,3) - tmp33

                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_i(3,3))
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_i(3,3))
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_i(3,3))
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_i(3,3))
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_i(3,3))
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_i(3,3))

                         ef2_j(1,1) = ef2_j(1,1) - tmp11
                         ef2_j(1,2) = ef2_j(1,2) - tmp12
                         ef2_j(1,3) = ef2_j(1,3) - tmp13
                         ef2_j(2,1) = ef2_j(2,1) - tmp12
                         ef2_j(2,2) = ef2_j(2,2) - tmp22
                         ef2_j(2,3) = ef2_j(2,3) - tmp23
                         ef2_j(3,1) = ef2_j(3,1) - tmp13
                         ef2_j(3,2) = ef2_j(3,2) - tmp23
                         ef2_j(3,3) = ef2_j(3,3) - tmp33
                      END IF
                   END IF
                END IF
                IF (task(3,2)) THEN
                   ! Quadrupole - Dipole
                   fac = 1.0_dp/3.0_dp
                   ! Dipole i (locally B) - Quadrupole j (locally A)
                   tmp_ij = dp_i(1)*(tij_abc(1,1,1)*qp_j(1,1)+&
                                     tij_abc(2,1,1)*qp_j(2,1)+&
                                     tij_abc(3,1,1)*qp_j(3,1)+&
                                     tij_abc(1,2,1)*qp_j(1,2)+&
                                     tij_abc(2,2,1)*qp_j(2,2)+&
                                     tij_abc(3,2,1)*qp_j(3,2)+&
                                     tij_abc(1,3,1)*qp_j(1,3)+&
                                     tij_abc(2,3,1)*qp_j(2,3)+&
                                     tij_abc(3,3,1)*qp_j(3,3))+&
                            dp_i(2)*(tij_abc(1,1,2)*qp_j(1,1)+&
                                     tij_abc(2,1,2)*qp_j(2,1)+&
                                     tij_abc(3,1,2)*qp_j(3,1)+&
                                     tij_abc(1,2,2)*qp_j(1,2)+&
                                     tij_abc(2,2,2)*qp_j(2,2)+&
                                     tij_abc(3,2,2)*qp_j(3,2)+&
                                     tij_abc(1,3,2)*qp_j(1,3)+&
                                     tij_abc(2,3,2)*qp_j(2,3)+&
                                     tij_abc(3,3,2)*qp_j(3,3))+&
                            dp_i(3)*(tij_abc(1,1,3)*qp_j(1,1)+&
                                     tij_abc(2,1,3)*qp_j(2,1)+&
                                     tij_abc(3,1,3)*qp_j(3,1)+&
                                     tij_abc(1,2,3)*qp_j(1,2)+&
                                     tij_abc(2,2,3)*qp_j(2,2)+&
                                     tij_abc(3,2,3)*qp_j(3,2)+&
                                     tij_abc(1,3,3)*qp_j(1,3)+&
                                     tij_abc(2,3,3)*qp_j(2,3)+&
                                     tij_abc(3,3,3)*qp_j(3,3))

                   ! Dipole j (locally A) - Quadrupole i (locally B)
                   tmp_ji = dp_j(1)*(tij_abc(1,1,1)*qp_i(1,1)+&
                                     tij_abc(2,1,1)*qp_i(2,1)+&
                                     tij_abc(3,1,1)*qp_i(3,1)+&
                                     tij_abc(1,2,1)*qp_i(1,2)+&
                                     tij_abc(2,2,1)*qp_i(2,2)+&
                                     tij_abc(3,2,1)*qp_i(3,2)+&
                                     tij_abc(1,3,1)*qp_i(1,3)+&
                                     tij_abc(2,3,1)*qp_i(2,3)+&
                                     tij_abc(3,3,1)*qp_i(3,3))+&
                            dp_j(2)*(tij_abc(1,1,2)*qp_i(1,1)+&
                                     tij_abc(2,1,2)*qp_i(2,1)+&
                                     tij_abc(3,1,2)*qp_i(3,1)+&
                                     tij_abc(1,2,2)*qp_i(1,2)+&
                                     tij_abc(2,2,2)*qp_i(2,2)+&
                                     tij_abc(3,2,2)*qp_i(3,2)+&
                                     tij_abc(1,3,2)*qp_i(1,3)+&
                                     tij_abc(2,3,2)*qp_i(2,3)+&
                                     tij_abc(3,3,2)*qp_i(3,3))+&
                            dp_j(3)*(tij_abc(1,1,3)*qp_i(1,1)+&
                                     tij_abc(2,1,3)*qp_i(2,1)+&
                                     tij_abc(3,1,3)*qp_i(3,1)+&
                                     tij_abc(1,2,3)*qp_i(1,2)+&
                                     tij_abc(2,2,3)*qp_i(2,2)+&
                                     tij_abc(3,2,3)*qp_i(3,2)+&
                                     tij_abc(1,3,3)*qp_i(1,3)+&
                                     tij_abc(2,3,3)*qp_i(2,3)+&
                                     tij_abc(3,3,3)*qp_i(3,3))

                   tmp= fac * (tmp_ij - tmp_ji)
                   eloc = eloc + tmp
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         ! Dipole i (locally B) - Quadrupole j (locally A)
                         tmp_ij = dp_i(1)*(tij_abcd(1,1,1,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,1,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,1,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,1,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,1,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,1,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,1,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,1,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,1,k)*qp_j(3,3))+&
                                  dp_i(2)*(tij_abcd(1,1,2,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,2,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,2,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,2,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,2,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,2,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,2,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,2,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,2,k)*qp_j(3,3))+&
                                  dp_i(3)*(tij_abcd(1,1,3,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,3,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,3,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,3,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,3,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,3,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,3,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,3,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,3,k)*qp_j(3,3))

                         ! Dipole j (locally A) - Quadrupole i (locally B)
                         tmp_ji = dp_j(1)*(tij_abcd(1,1,1,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,1,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,1,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,1,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,1,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,1,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,1,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,1,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,1,k)*qp_i(3,3))+&
                                  dp_j(2)*(tij_abcd(1,1,2,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,2,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,2,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,2,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,2,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,2,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,2,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,2,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,2,k)*qp_i(3,3))+&
                                  dp_j(3)*(tij_abcd(1,1,3,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,3,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,3,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,3,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,3,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,3,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,3,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,3,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac * (tmp_ij - tmp_ji)
                      END DO
                   END IF
                END IF
                IF (task(3,1)) THEN
                   ! Quadrupole - Charge
                   fac = 1.0_dp/3.0_dp

                   ! Quadrupole j (locally A) - Charge j (locally B)
                   tmp_ij = ch_i * (tij_ab(1,1)*qp_j(1,1)+&
                                    tij_ab(2,1)*qp_j(2,1)+&
                                    tij_ab(3,1)*qp_j(3,1)+&
                                    tij_ab(1,2)*qp_j(1,2)+&
                                    tij_ab(2,2)*qp_j(2,2)+&
                                    tij_ab(3,2)*qp_j(3,2)+&
                                    tij_ab(1,3)*qp_j(1,3)+&
                                    tij_ab(2,3)*qp_j(2,3)+&
                                    tij_ab(3,3)*qp_j(3,3))

                   ! Quadrupole i (locally B) - Charge j (locally A)
                   tmp_ji = ch_j * (tij_ab(1,1)*qp_i(1,1)+&
                                    tij_ab(2,1)*qp_i(2,1)+&
                                    tij_ab(3,1)*qp_i(3,1)+&
                                    tij_ab(1,2)*qp_i(1,2)+&
                                    tij_ab(2,2)*qp_i(2,2)+&
                                    tij_ab(3,2)*qp_i(3,2)+&
                                    tij_ab(1,3)*qp_i(1,3)+&
                                    tij_ab(2,3)*qp_i(2,3)+&
                                    tij_ab(3,3)*qp_i(3,3))

                   eloc = eloc + fac*(tmp_ij+tmp_ji)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         ! Quadrupole j (locally A) - Charge i (locally B)
                         tmp_ij = ch_i * (tij_abc(1,1,k)*qp_j(1,1)+&
                                          tij_abc(2,1,k)*qp_j(2,1)+&
                                          tij_abc(3,1,k)*qp_j(3,1)+&
                                          tij_abc(1,2,k)*qp_j(1,2)+&
                                          tij_abc(2,2,k)*qp_j(2,2)+&
                                          tij_abc(3,2,k)*qp_j(3,2)+&
                                          tij_abc(1,3,k)*qp_j(1,3)+&
                                          tij_abc(2,3,k)*qp_j(2,3)+&
                                          tij_abc(3,3,k)*qp_j(3,3))

                         ! Quadrupole i (locally B) - Charge j (locally A)
                         tmp_ji = ch_j * (tij_abc(1,1,k)*qp_i(1,1)+&
                                          tij_abc(2,1,k)*qp_i(2,1)+&
                                          tij_abc(3,1,k)*qp_i(3,1)+&
                                          tij_abc(1,2,k)*qp_i(1,2)+&
                                          tij_abc(2,2,k)*qp_i(2,2)+&
                                          tij_abc(3,2,k)*qp_i(3,2)+&
                                          tij_abc(1,3,k)*qp_i(1,3)+&
                                          tij_abc(2,3,k)*qp_i(2,3)+&
                                          tij_abc(3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac *(tmp_ij + tmp_ji)
                      END DO
                   END IF
                END IF
# 1108 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/ewalds_multipole_sr_pure.f90"
                ! Electric fields
                IF (do_efield) THEN
                   ! Potential
                   IF (do_efield0) THEN
                      efield0(  atom_a) = efield0(  atom_a) + ef0_j

                      efield0(  atom_b) = efield0(  atom_b) + ef0_i
                   END IF
                   ! Electric field
                   IF (do_efield1) THEN
                      efield1(1,atom_a) = efield1(1,atom_a) + ef1_j(1)
                      efield1(2,atom_a) = efield1(2,atom_a) + ef1_j(2)
                      efield1(3,atom_a) = efield1(3,atom_a) + ef1_j(3)

                      efield1(1,atom_b) = efield1(1,atom_b) + ef1_i(1)
                      efield1(2,atom_b) = efield1(2,atom_b) + ef1_i(2)
                      efield1(3,atom_b) = efield1(3,atom_b) + ef1_i(3)
                   END IF
                   ! Electric field gradient
                   IF (do_efield2) THEN
                      efield2(1,atom_a) = efield2(1,atom_a) + ef2_j(1,1)
                      efield2(2,atom_a) = efield2(2,atom_a) + ef2_j(1,2)
                      efield2(3,atom_a) = efield2(3,atom_a) + ef2_j(1,3)
                      efield2(4,atom_a) = efield2(4,atom_a) + ef2_j(2,1)
                      efield2(5,atom_a) = efield2(5,atom_a) + ef2_j(2,2)
                      efield2(6,atom_a) = efield2(6,atom_a) + ef2_j(2,3)
                      efield2(7,atom_a) = efield2(7,atom_a) + ef2_j(3,1)
                      efield2(8,atom_a) = efield2(8,atom_a) + ef2_j(3,2)
                      efield2(9,atom_a) = efield2(9,atom_a) + ef2_j(3,3)

                      efield2(1,atom_b) = efield2(1,atom_b) + ef2_i(1,1)
                      efield2(2,atom_b) = efield2(2,atom_b) + ef2_i(1,2)
                      efield2(3,atom_b) = efield2(3,atom_b) + ef2_i(1,3)
                      efield2(4,atom_b) = efield2(4,atom_b) + ef2_i(2,1)
                      efield2(5,atom_b) = efield2(5,atom_b) + ef2_i(2,2)
                      efield2(6,atom_b) = efield2(6,atom_b) + ef2_i(2,3)
                      efield2(7,atom_b) = efield2(7,atom_b) + ef2_i(3,1)
                      efield2(8,atom_b) = efield2(8,atom_b) + ef2_i(3,2)
                      efield2(9,atom_b) = efield2(9,atom_b) + ef2_i(3,3)
                   END IF
                END IF
                IF (do_stress) THEN
                   ptens11 = ptens11 + rab(1) * fr(1)
                   ptens21 = ptens21 + rab(2) * fr(1)
                   ptens31 = ptens31 + rab(3) * fr(1)
                   ptens12 = ptens12 + rab(1) * fr(2)
                   ptens22 = ptens22 + rab(2) * fr(2)
                   ptens32 = ptens32 + rab(3) * fr(2)
                   ptens13 = ptens13 + rab(1) * fr(3)
                   ptens23 = ptens23 + rab(2) * fr(3)
                   ptens33 = ptens33 + rab(3) * fr(3)
                END IF
# 1183 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/se_fock_matrix_integrals.F" 2

    IF (PRESENT(integral_value)) THEN
       integral_value = eloc
    END IF
    IF (do_forces) THEN
       force_ab = fr
    END IF
  END SUBROUTINE se_coulomb_ij_interaction

END MODULE se_fock_matrix_integrals

