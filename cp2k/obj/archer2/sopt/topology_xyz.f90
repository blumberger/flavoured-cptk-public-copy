# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/topology_xyz.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/topology_xyz.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE topology_xyz
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type

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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/topology_xyz.F" 2

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_xyz'

  PRIVATE
  PUBLIC :: read_coordinate_xyz

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE read_coordinate_xyz (topology,para_env,subsys_section)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_xyz', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_default_index, strtmp
    INTEGER                                  :: frame, handle, iw, j, natom
    LOGICAL                                  :: my_end
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    CALL timeset(routineN,handle)

    NULLIFY (parser,logger)
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                              extension=".subsysLog")

    atom_info => topology%atom_info

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "BEGIN of XYZ data read from file "//TRIM(topology%coord_file_name)
    END IF

    CALL parser_create(parser,topology%coord_file_name, para_env=para_env,&
                       parse_white_lines=.TRUE.)

    ! Element is assigned on the basis of the atm_name
    topology%aa_element = .TRUE.

    natom = 0
    frame = 0
    CALL parser_get_next_line(parser,1)
    Frames: DO
       ! Atom numbers
       CALL parser_get_object(parser,natom)
       frame = frame + 1
       IF (frame == 1) THEN
          CALL reallocate(atom_info%id_molname,1,natom)
          CALL reallocate(atom_info%id_resname,1,natom)
          CALL reallocate(atom_info%resid,1,natom)
          CALL reallocate(atom_info%id_atmname,1,natom)
          CALL reallocate(atom_info%r,1,3,1,natom)
          CALL reallocate(atom_info%atm_mass,1,natom)
          CALL reallocate(atom_info%atm_charge,1,natom)
          CALL reallocate(atom_info%occup,1,natom)
          CALL reallocate(atom_info%beta,1,natom)
          CALL reallocate(atom_info%id_element,1,natom)
       ELSE IF (natom > SIZE(atom_info%id_atmname)) THEN
          CALL cp__b("topology_xyz.F",99,"Atom number differs in different frames!")
       END IF
       ! Dummy line
       CALL parser_get_next_line(parser,2)
       DO j=1,natom
          ! Atom coordinates
          READ (parser%input_line,*) strtmp,&
                                     atom_info%r(1,j),&
                                     atom_info%r(2,j),&
                                     atom_info%r(3,j)
          atom_info%id_atmname(j) = str2id(s2s(strtmp))
          ! For default, set atom name to residue name to molecule name
          WRITE (my_default_index,'(I0)') j
          atom_info%id_molname(j) = str2id(s2s(TRIM(id2str(atom_info%id_atmname(j)))//TRIM(my_default_index)))
          atom_info%id_resname(j) = atom_info%id_molname(j)
          atom_info%resid(j)      = 1
          atom_info%id_element(j) = atom_info%id_atmname(j)
          atom_info%atm_mass(j)   =  HUGE(0.0_dp)
          atom_info%atm_charge(j) = -HUGE(0.0_dp)
          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="(T2,A4,3F8.3,2X,A)")&
               TRIM(id2str(atom_info%id_atmname(j))),&
               atom_info%r(1,j),&
               atom_info%r(2,j),&
               atom_info%r(3,j),&
               ADJUSTL(TRIM(id2str(atom_info%id_molname(j))))
          END IF
          atom_info%r(1,j) = cp_unit_to_cp2k(atom_info%r(1,j),"angstrom")
          atom_info%r(2,j) = cp_unit_to_cp2k(atom_info%r(2,j),"angstrom")
          atom_info%r(3,j) = cp_unit_to_cp2k(atom_info%r(3,j),"angstrom")
          ! If there's a white line or end of file exit.. otherwise read other available
          ! snapshots
          CALL parser_get_next_line(parser,1,at_end=my_end)
          my_end = my_end.OR.(LEN_TRIM(parser%input_line) == 0)
          IF (my_end) THEN
             IF(j/=natom)&
                CALL cp_abort(cp__l("topology_xyz.F",135),&
                     "Number of lines in XYZ format not equal to the number of atoms."//&
                     " Error in XYZ format. Very probably the line with title is missing or is empty."//&
                     " Please check the XYZ file and rerun your job!")
             EXIT Frames
          END IF
       END DO
    END DO Frames
    CALL parser_release(parser)

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "END of XYZ frame data read from file "//TRIM(topology%coord_file_name)
    END IF

    topology%natoms = natom
    topology%molname_generated = .TRUE.

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
                                      "PRINT%TOPOLOGY_INFO/XYZ_INFO")

    CALL timestop(handle)

  END SUBROUTINE read_coordinate_xyz

END MODULE topology_xyz
