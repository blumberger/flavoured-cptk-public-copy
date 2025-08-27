# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_env.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_env.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods related to the blacs parallel environment
!> \par History
!>      08.2002 created [fawzi]
!>      02.2004 modified to associate a blacs_env with a given para_env
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_blacs_env
  USE cp_array_utils_i,                ONLY: cp_2d_i_write
  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinfo,&
                                             cp_blacs_gridinit,&
                                             cp_blacs_set
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: gcd
  USE message_passing,                 ONLY: mp_sum

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
# 27 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_env.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_env'

  ! Blacs type of distribution
  INTEGER, PARAMETER, PUBLIC               :: BLACS_GRID_SQUARE = 1,&
                                              BLACS_GRID_ROW    = 2,&
                                              BLACS_GRID_COL    = 3

  PUBLIC :: cp_blacs_env_type ! make it accessible only through cp_para_types?
  PUBLIC :: cp_blacs_env_create, cp_blacs_env_retain, cp_blacs_env_release
  PUBLIC :: cp_blacs_env_write, get_blacs_info

! *****************************************************************************
!> \brief represent a blacs multidimensional parallel environment
!>      (for the mpi corrispective see cp_paratypes/cp_para_cart_type)
!> \param mepos the position of the actual processor (2D)
!> \param group id of the actual group (context, communicator)
!> \param num_pe number of processors in the group in each dimension
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \param my_pid process id of the actual processor
!> \param n_pid number of process ids
!> \param the para env associated (and compatible) with this blacs_env
!> \param blacs 2mpi: maps mepos(1)-mepos(2) of blacs to its mpi rank
!> \param mpi 2blacs(i,rank): maps the mpi rank to the mepos(i)
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_blacs_env_type
     INTEGER, DIMENSION(2) :: mepos, num_pe
     INTEGER :: group, my_pid, n_pid, ref_count
     TYPE(cp_para_env_type), POINTER :: para_env
     INTEGER, DIMENSION(:,:), POINTER :: blacs2mpi
     INTEGER, DIMENSION(:,:), POINTER :: mpi2blacs
     LOGICAL :: repeatable
  END TYPE cp_blacs_env_type

!***
CONTAINS

! *****************************************************************************
!> \brief   Return informations about the specified BLACS context.
!> \param blacs_env ...
!> \param my_process_row ...
!> \param my_process_column ...
!> \param my_process_number ...
!> \param number_of_process_rows ...
!> \param number_of_process_columns ...
!> \param number_of_processes ...
!> \param para_env ...
!> \param blacs2mpi ...
!> \param mpi2blacs ...
!> \date    19.06.2001
!> \par     History
!>          MM.YYYY moved here from qs_blacs (Joost VandeVondele)
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_blacs_info(blacs_env,my_process_row,my_process_column,&
                            my_process_number,number_of_process_rows,&
                            number_of_process_columns,number_of_processes,&
                            para_env, blacs2mpi, mpi2blacs)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    INTEGER, INTENT(OUT), OPTIONAL :: my_process_row, my_process_column, &
      my_process_number, number_of_process_rows, number_of_process_columns, &
      number_of_processes
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: blacs2mpi, mpi2blacs

    CHARACTER(len=*), PARAMETER :: routineN = 'get_blacs_info', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(blacs_env)) THEN
       CALL cp__b("fm/cp_blacs_env.F",107,"No BLACS environment")
    END IF

    IF (PRESENT(my_process_row)) my_process_row = blacs_env%mepos(1)
    IF (PRESENT(my_process_column)) my_process_column = blacs_env%mepos(2)
    IF (PRESENT(my_process_number)) my_process_number = blacs_env%my_pid
    IF (PRESENT(number_of_process_rows)) number_of_process_rows = blacs_env%num_pe(1)
    IF (PRESENT(number_of_process_columns)) number_of_process_columns = blacs_env%num_pe(2)
    IF (PRESENT(number_of_processes)) number_of_processes = blacs_env%n_pid
    IF (PRESENT(para_env)) para_env => blacs_env%para_env
    IF (PRESENT(blacs2mpi)) blacs2mpi => blacs_env%blacs2mpi
    IF (PRESENT(mpi2blacs)) mpi2blacs => blacs_env%mpi2blacs

  END SUBROUTINE get_blacs_info

! *****************************************************************************
!> \brief allocates and initializes a type that represent a blacs context
!> \param blacs_env the type to initialize
!> \param para_env the para_env for which a blacs env should be created
!> \param blacs_grid_layout ...
!> \param blacs_repeatable ...
!> \param row_major ...
!> \param grid_2d ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_blacs_env_create(blacs_env,para_env,blacs_grid_layout,blacs_repeatable,row_major,grid_2d)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: blacs_grid_layout
    LOGICAL, INTENT(IN), OPTIONAL            :: blacs_repeatable, row_major
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: grid_2d

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipcol, iprow, stat

















  IF(.NOT.(.NOT.ASSOCIATED(blacs_env)))CALL cp__a("fm/cp_blacs_env.F",163)

  ALLOCATE(blacs_env)
  blacs_env%group=0
  blacs_env%ref_count=1
  blacs_env%mepos(:)=0
  blacs_env%num_pe(:)=1
  blacs_env%my_pid=0
  blacs_env%n_pid=1
  CALL cp_para_env_retain(para_env)
  blacs_env%para_env => para_env

# 238 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_blacs_env.F"
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(blacs_grid_layout))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(blacs_repeatable))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(grid_2d))==-1) EXIT ;  END DO ; ENDIF
  IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(row_major))==-1) EXIT ;  END DO ; ENDIF


  ! generate the mappings blacs2mpi and mpi2blacs
  ALLOCATE(blacs_env%blacs2mpi(0:blacs_env%num_pe(1)-1,0:blacs_env%num_pe(2)-1),&
       stat=stat)
  IF(.NOT.(stat==0))CALL cp__a("fm/cp_blacs_env.F",247)
  blacs_env%blacs2mpi=0
  blacs_env%blacs2mpi(blacs_env%mepos(1),blacs_env%mepos(2))=para_env%mepos
  CALL mp_sum(blacs_env%blacs2mpi,para_env%group)
  ALLOCATE(blacs_env%mpi2blacs(2,0:para_env%num_pe-1))
  blacs_env%mpi2blacs=-1
  DO ipcol=0,blacs_env%num_pe(2)-1
     DO iprow=0,blacs_env%num_pe(1)-1
        blacs_env%mpi2blacs(1,blacs_env%blacs2mpi(iprow,ipcol))=iprow
        blacs_env%mpi2blacs(2,blacs_env%blacs2mpi(iprow,ipcol))=ipcol
     END DO
  END DO
END SUBROUTINE cp_blacs_env_create

! *****************************************************************************
!> \brief retains the given blacs env
!> \param blacs_env the blacs env to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_blacs_env_retain(blacs_env)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_retain', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(blacs_env)))CALL cp__a("fm/cp_blacs_env.F",274)
  IF(.NOT.(blacs_env%ref_count>0))CALL cp__a("fm/cp_blacs_env.F",275)
  blacs_env%ref_count=blacs_env%ref_count+1
END SUBROUTINE cp_blacs_env_retain

! *****************************************************************************
!> \brief releases the given blacs_env
!> \param blacs_env the blacs env to relase
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_blacs_env_release(blacs_env)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_release', &
      routineP = moduleN//':'//routineN

  IF (ASSOCIATED(blacs_env)) THEN
     IF(.NOT.(blacs_env%ref_count>0))CALL cp__a("fm/cp_blacs_env.F",293)
     blacs_env%ref_count=blacs_env%ref_count-1
     IF (blacs_env%ref_count<1) THEN
        CALL cp_blacs_gridexit(blacs_env%group)
        CALL cp_para_env_release(blacs_env%para_env)
        DEALLOCATE(blacs_env%mpi2blacs)
        DEALLOCATE(blacs_env%blacs2mpi)
        DEALLOCATE(blacs_env)
     END IF
  END IF
  NULLIFY(blacs_env)
END SUBROUTINE cp_blacs_env_release

! *****************************************************************************
!> \brief writes the description of the given blacs env
!> \param blacs_env the blacs environment to write
!> \param unit_nr the unit number where to write the description of the
!>        blacs environment
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_blacs_env_write(blacs_env, unit_nr)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    INTEGER, INTENT(in)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iostat

  IF (ASSOCIATED(blacs_env)) THEN
     WRITE (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) blacs_env%group, blacs_env%ref_count
     IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",327)
     WRITE (unit=unit_nr,fmt="('  mepos=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%mepos(1), blacs_env%mepos(2)
     IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",330)
     WRITE (unit=unit_nr,fmt="('  num_pe=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%num_pe(1), blacs_env%num_pe(2)
     IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",333)
     IF (ASSOCIATED(blacs_env%blacs2mpi)) THEN
        WRITE (unit=unit_nr,fmt="('  blacs2mpi=')",advance="no",iostat=iostat)
        IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",336)
        CALL cp_2d_i_write(blacs_env%blacs2mpi,unit_nr=unit_nr)
     ELSE
        WRITE (unit=unit_nr,fmt="('  blacs2mpi=*null*')",iostat=iostat)
        IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",340)
     END IF
     IF (ASSOCIATED(blacs_env%para_env)) THEN
        WRITE (unit=unit_nr,fmt="('  para_env=<cp_para_env id=',i6,'>,')")&
             blacs_env%para_env%group
     ELSE
        WRITE (unit=unit_nr,fmt="('  para_env=*null*')")
     END IF
     WRITE (unit=unit_nr,fmt="('  my_pid=',i10,', n_pid=',i10,' }')",&
          iostat=iostat) blacs_env%my_pid, blacs_env%n_pid
     IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",350)
  ELSE
     WRITE (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_blacs_env>:*null* '
     IF(.NOT.(iostat==0))CALL cp__a("fm/cp_blacs_env.F",354)
  END IF
  CALL m_flush(unit_nr)
END SUBROUTINE cp_blacs_env_write

END MODULE cp_blacs_env
