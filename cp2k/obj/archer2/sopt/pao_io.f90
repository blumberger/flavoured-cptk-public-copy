# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_io.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_io.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for reading and writing restart files.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_io
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_row_block_sizes
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE pao_input,                       ONLY: id2str
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type

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
# 42 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/pao_io.F" 2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_io'

  PUBLIC :: pao_read_restart, pao_write_restart

CONTAINS

! *****************************************************************************
!> \brief Reads restart file
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_read_restart(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: label
    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: i, iatom, natoms, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: explicit, found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, buffer
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input

    CALL get_qs_env(qs_env, input=input, para_env=para_env)

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%READ_RESTART",&
       c_val=filename, explicit=explicit)

    IF(.NOT. explicit) RETURN

    IF(pao%iw>0) WRITE(pao%iw,'(A,A)') " PAO| Reading restart file: ", TRIM(filename)

    row_blk_sizes => cp_dbcsr_row_block_sizes(pao%matrix_X)
    col_blk_sizes => cp_dbcsr_col_block_sizes(pao%matrix_X)
    IF(SIZE(row_blk_sizes) /= SIZE(col_blk_sizes)) CALL cp__b("pao_io.F",85,"matrix_X not squared")
    natoms = SIZE(col_blk_sizes)

    unit_nr = -1
    IF(para_env%mepos == para_env%source) THEN
       CALL open_file(file_name=filename, file_status="OLD", file_form="FORMATTED",&
                      file_action="READ", unit_number=unit_nr)

       CALL read_restart_header(pao, qs_env, unit_nr)
    ENDIF

    !TODO: this is a serial algorithm
    DO i=1, natoms
       IF(unit_nr>0) THEN
          READ(unit_nr,fmt=*) label, iatom
          IF(.NOT.(TRIM(label)=="Xblock"))CALL cp__a("pao_io.F",100)
       ENDIF
       CALL mp_bcast(iatom, para_env%source, para_env%group)
       IF(pao%iw>0)  WRITE(pao%iw,*) "PAO| Restart found iatom: ", iatom
       IF(.NOT.(iatom <= natoms))CALL cp__a("pao_io.F",104)

       ALLOCATE(buffer(row_blk_sizes(iatom), col_blk_sizes(iatom)))
       IF(unit_nr>0) THEN
          BACKSPACE(unit_nr)
          READ(unit_nr,fmt=*) label, iatom, buffer
       ENDIF
       CALL mp_bcast(buffer, para_env%source, para_env%group)
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)
       IF(ASSOCIATED(block_X)) THEN
          block_X = buffer
       ENDIF
       DEALLOCATE(buffer)
    ENDDO

    IF(unit_nr>0) CALL close_file(unit_number=unit_nr)
  END SUBROUTINE pao_read_restart


! *****************************************************************************
!> \brief Read and check header of restart file
!> \param pao ...
!> \param qs_env ...
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE read_restart_header(pao, qs_env, unit_nr)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'read_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kindname, label, str_in
    INTEGER                                  :: i1, i2, i3, iatom, ikind, &
                                                pao_basis_size, pao_pot_maxl, &
                                                pao_pot_neighbors, z
    REAL(dp)                                 :: diff, pao_pot_beta, r1
    REAL(dp), DIMENSION(3)                   :: pos_in
    REAL(dp), DIMENSION(3, 3)                :: hmat_in
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

   CALL get_qs_env(qs_env,&
                   cell=cell,&
                   particle_set=particle_set,&
                   atomic_kind_set=atomic_kind_set,&
                   qs_kind_set=qs_kind_set)

   DO WHILE(.TRUE.)
       READ(unit_nr, fmt=*) label
       BACKSPACE(unit_nr)

       IF(TRIM(label) == "Parametrization") THEN
          READ(unit_nr,fmt=*) label, str_in
          IF(TRIM(str_in) .NE. TRIM(ADJUSTL(id2str(pao%parameterization))))&
             CALL cp__b("pao_io.F",166,"Restart PAO parametrization does not match")

       ELSE IF(TRIM(label) == "Cell") THEN
          READ(unit_nr,fmt=*) label, hmat_in
          diff = MAXVAL(ABS(hmat_in - cell%hmat * angstrom))
          IF(diff > 1e-14)&
             CALL cp__w("pao_io.F",172,"Restaring from different cell dimension")

       ELSE IF(TRIM(label) == "Kind") THEN
          READ(unit_nr, fmt=*) label, ikind, str_in, i1
          CALL get_atomic_kind(atomic_kind_set(ikind), name=kindname, z=z)
          IF(str_in .NE. kindname)&
             CALL cp__b("pao_io.F",178,"Kind names do not match")
          IF(i1 /= z)&
             CALL cp__b("pao_io.F",180,"Atomic numbers do not match")

       ELSE IF(TRIM(label) == "PrimBasis") THEN
          READ(unit_nr, fmt=*) label, ikind, str_in
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
          IF(TRIM(str_in) .NE. TRIM(basis_set%name))&
             CALL cp__b("pao_io.F",186,"Primary Basis set does not match")

       ELSE IF(TRIM(label) == "PaoBasis") THEN
          READ(unit_nr, fmt=*) label, ikind, i1, i2, i3, r1
          CALL get_qs_kind(qs_kind_set(ikind),&
                           pao_basis_size=pao_basis_size,&
                           pao_potential_maxl=pao_pot_maxl,&
                           pao_potential_neighbors=pao_pot_neighbors,&
                           pao_potential_beta=pao_pot_beta)
          IF(i1 /= pao_basis_size)&
             CALL cp__b("pao_io.F",196,"PAO_BASIS_SIZE does not match")
          IF(i2 /= pao_pot_maxl)&
             CALL cp__b("pao_io.F",198,"PAO_POT_MAXL does not match")
          IF(i3 /= pao_pot_neighbors)&
             CALL cp__b("pao_io.F",200,"PAO_POT_NEIGHBORS does not match")
          IF(r1 /= pao_pot_beta)&
             CALL cp__b("pao_io.F",202,"PAO_POT_BETA does not match")

       ELSE IF(TRIM(label) == "Atom") THEN
          READ(unit_nr,fmt=*) label, iatom, str_in, pos_in
          IF(str_in .NE. particle_set(iatom)%atomic_kind%name)&
             CALL cp__b("pao_io.F",207,"Restart atomic kinds do not match.")
          diff = MAXVAL(ABS(pos_in - particle_set(iatom)%r * angstrom))
          IF(diff > 1e-10)&
             CALL cp__w("pao_io.F",210,"Restarting from different atom positions")

       ELSE IF(TRIM(label) == "Xblock") THEN
          EXIT ! end of header

       ELSE
          !CALL cp__w("pao_io.F",216,"Skipping restart header with label: "//TRIM(label))
          READ(unit_nr,fmt=*) label  ! just read again and ignore
       ENDIF
   ENDDO
 END SUBROUTINE read_restart_header


! *****************************************************************************
!> \brief Writes restart file
!> \param pao ...
!> \param qs_env ...
!> \param energy ...
! *****************************************************************************
  SUBROUTINE pao_write_restart(pao, qs_env, energy)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: energy

    CHARACTER(len=*), PARAMETER :: &
      printkey_section = 'DFT%LS_SCF%PAO%PRINT%RESTART', &
      routineN = 'pao_write_restart', routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, natoms, &
                                                unit_max, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, buffer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()

    CALL get_qs_env(qs_env,&
                    input=input,&
                    natom=natoms,&
                    para_env=para_env)

    ! open file
    unit_nr = cp_print_key_unit_nr(logger,&
                              input,&
                              printkey_section,&
                              extension=".pao",&
                              file_action="WRITE",&
                              file_position="REWIND",&
                              file_status="UNKNOWN",&
                              do_backup=.TRUE.)

    ! although just rank-0 writes the trajectory it requires collective MPI calls
    unit_max = unit_nr
    CALL mp_max(unit_max, para_env%group)
    IF(unit_max>0) THEN
       IF(pao%iw>0) WRITE(pao%iw,'(A,A)') " PAO| Writing restart file."
       IF(unit_nr>0)&
          CALL write_restart_header(pao, qs_env, energy, unit_nr)

       !TODO: this is a serial algorithm
       ! write matrix_X
       row_blk_sizes => cp_dbcsr_row_block_sizes(pao%matrix_X)
       col_blk_sizes => cp_dbcsr_col_block_sizes(pao%matrix_X)
       DO iatom=1, natoms
          ALLOCATE(buffer(row_blk_sizes(iatom), col_blk_sizes(iatom)))
          NULLIFY(block_X)
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)
          IF(ASSOCIATED(block_X)) THEN
             IF(SIZE(block_X) > 0) & ! corner-case of zero pao parameters
               buffer(:,:) = block_X(:,:)
          ELSE
             buffer(:,:) = 0.0_dp
          ENDIF
          CALL mp_sum(buffer, para_env%group)
          IF(unit_nr>0) THEN
             WRITE(unit_nr, fmt="(A,I10,1X)",advance='no') "Xblock ", iatom
             WRITE(unit_nr, *) buffer
          ENDIF
          DEALLOCATE(buffer)
       ENDDO

       ! flush
       IF(unit_nr>0) FLUSH(unit_nr)

    ENDIF

    ! close file
    CALL cp_print_key_finished_output(unit_nr, logger, input, printkey_section)

    CALL timestop(handle)
  END SUBROUTINE pao_write_restart


! *****************************************************************************
!> \brief Writes header of restart file
!> \param pao ...
!> \param qs_env ...
!> \param energy ...
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE write_restart_header(pao, qs_env, energy, unit_nr)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: energy
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'write_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kindname
    INTEGER                                  :: iatom, ikind, pao_basis_size, &
                                                pao_pot_maxl, &
                                                pao_pot_neighbors, z
    REAL(dp)                                 :: pao_pot_beta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL get_qs_env(qs_env,&
                    cell=cell,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set)

    ! write kinds
    WRITE(unit_nr, "(A,5X,F20.10)") "Energy", energy
    WRITE(unit_nr, "(A,5X,I0)") "Step", pao%istep
    WRITE(unit_nr,"(A,5X,A10)") "Parametrization", id2str(pao%parameterization)
    DO ikind=1, SIZE(atomic_kind_set)
       CALL get_atomic_kind(atomic_kind_set(ikind), name=kindname, z=z)
       CALL get_qs_kind(qs_kind_set(ikind),&
                        pao_basis_size=pao_basis_size,&
                        pao_potential_maxl=pao_pot_maxl,&
                        pao_potential_neighbors=pao_pot_neighbors,&
                        pao_potential_beta=pao_pot_beta,&
                        basis_set=basis_set)
       WRITE(unit_nr,"(A,5X,I10,1X,A,1X,I3)") "Kind", ikind, TRIM(kindname), z
       WRITE(unit_nr,"(A,5X,I10,1X,A)") "PrimBasis", ikind, TRIM(basis_set%name)
       WRITE(unit_nr,"(A,5X,I10,1X,I3,1X,I3)",advance='no') "PaoBasis", ikind, pao_basis_size, pao_pot_maxl
       WRITE(unit_nr,"(1X,I3,1X,F20.16)") pao_pot_neighbors, pao_pot_beta
    ENDDO

    ! write cell
    WRITE(unit_nr,fmt="(A,5X)",advance='no') "Cell"
    WRITE(unit_nr,*) cell%hmat * angstrom

    ! write particle positions
    DO iatom=1, SIZE(particle_set)
       kindname = particle_set(iatom)%atomic_kind%name
       WRITE(unit_nr,fmt="(A,5X,I10,5X,A,1X)",advance='no') "Atom ", iatom, TRIM(kindname)
       WRITE(unit_nr,*) particle_set(iatom)%r * angstrom
    ENDDO

  END SUBROUTINE write_restart_header

END MODULE pao_io
