# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR block access
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2010-02-18 Moved from dbcsr_operations and dbcsr_methods
!  - 2010-04-22 Added block buffer operations
! *****************************************************************************
MODULE dbcsr_block_access
  USE array_types,                     ONLY: array_data
  USE btree_I8_k_cp2d_v,               ONLY: btree_2d_data_c => cp2d,&
                                             btree_add_c => btree_add,&
                                             btree_get_c => btree_find
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d,&
                                             btree_add_d => btree_add,&
                                             btree_get_d => btree_find
  USE btree_I8_k_sp2d_v,               ONLY: btree_2d_data_s => sp2d,&
                                             btree_add_s => btree_add,&
                                             btree_get_s => btree_find
  USE btree_I8_k_zp2d_v,               ONLY: btree_2d_data_z => zp2d,&
                                             btree_add_z => btree_add,&
                                             btree_get_z => btree_find
  USE dbcsr_block_operations,          ONLY: dbcsr_data_clear,&
                                             dbcsr_data_set
  USE dbcsr_config,                    ONLY: default_resize_factor
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size_referenced,&
                                             dbcsr_data_set_pointer,&
                                             dbcsr_get_data,&
                                             dbcsr_get_data_p
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_local_cols,&
                                             dbcsr_distribution_local_rows,&
                                             dbcsr_distribution_mp
  USE dbcsr_dist_operations,           ONLY: dbcsr_get_block_index,&
                                             dbcsr_get_stored_block_info,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_failure_level,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error,&
                                             dbcsr_wrong_args_error
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_expand_row_index,&
                                             dbcsr_make_dbcsr_index,&
                                             dbcsr_sort_indices,&
                                             merge_index_arrays
  USE dbcsr_methods,                   ONLY: &
       dbcsr_blk_column_size, dbcsr_blk_row_size, dbcsr_distribution, &
       dbcsr_get_data_type, dbcsr_get_num_blocks, dbcsr_mutable_instantiated, &
       dbcsr_mutable_new, dbcsr_nblkrows_total, dbcsr_use_mutable, &
       dbcsr_wm_use_mutable
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_mynode
  USE dbcsr_ptr_util,                  ONLY: pointer_rank_remap2,&
                                             pointer_view
  USE dbcsr_toollib,                   ONLY: make_coordinate_tuple,&
                                             swap
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_obj, dbcsr_scalar_type, dbcsr_slot_blk_p, &
       dbcsr_slot_col_i, dbcsr_slot_nblks, dbcsr_slot_nze, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_real_4, dbcsr_type_real_4_2d, &
       dbcsr_type_real_8, dbcsr_type_real_8_2d
  USE dbcsr_work_operations,           ONLY: add_work_coordinate,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: dp,&
                                             real_4,&
                                             real_8

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/../../base/base_uses.f90" 1
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
# 79 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F" 2

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_block_access'

  PUBLIC :: dbcsr_get_block_p,&
            dbcsr_put_block, dbcsr_remove_block

  PUBLIC :: dbcsr_reserve_block2d,&
            dbcsr_reserve_blocks, dbcsr_reserve_all_blocks, dbcsr_reserve_diag_blocks


  INTERFACE dbcsr_get_block_p
     MODULE PROCEDURE dbcsr_get_block_p_d, dbcsr_get_block_p_s,&
                      dbcsr_get_block_p_z, dbcsr_get_block_p_c
     MODULE PROCEDURE dbcsr_get_2d_block_p_d, dbcsr_get_2d_block_p_s,&
                      dbcsr_get_2d_block_p_z, dbcsr_get_2d_block_p_c
     MODULE PROCEDURE dbcsr_get_block_p_area
  END INTERFACE

  INTERFACE dbcsr_put_block
     MODULE PROCEDURE dbcsr_put_block_area
     MODULE PROCEDURE dbcsr_put_block_d, dbcsr_put_block_s,&
                      dbcsr_put_block_z, dbcsr_put_block_c
     MODULE PROCEDURE dbcsr_put_block2d_d, dbcsr_put_block2d_s,&
                      dbcsr_put_block2d_z, dbcsr_put_block2d_c
  END INTERFACE

  INTERFACE dbcsr_reserve_block2d
     MODULE PROCEDURE dbcsr_reserve_block2d_s, dbcsr_reserve_block2d_d,&
                      dbcsr_reserve_block2d_c, dbcsr_reserve_block2d_z
  END INTERFACE

  INTERFACE dbcsr_set_block_pointer
     MODULE PROCEDURE dbcsr_set_block_pointer_any
     MODULE PROCEDURE dbcsr_set_block_pointer_2d_s,&
                      dbcsr_set_block_pointer_2d_d,&
                      dbcsr_set_block_pointer_2d_c,&
                      dbcsr_set_block_pointer_2d_z
  END INTERFACE



  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: debug_mod = .FALSE.


  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


  LOGICAL, PARAMETER, PRIVATE :: detailed_timing = .FALSE.

  TYPE block_parameters
     LOGICAL :: tr
     INTEGER :: logical_rows, logical_cols
     INTEGER :: offset, nze
  END TYPE block_parameters

  TYPE dgemm_join
     INTEGER :: p_a, p_b, p_c
     INTEGER :: last_k, last_n
     TYPE(dbcsr_scalar_type) :: alpha, beta
  END TYPE dgemm_join

CONTAINS

! *****************************************************************************
!> \brief Marks a block for removal from a DBCSR matrix. Handles
!>        symmetric matrices.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            row of block to remove
!> \param[in]  col            column of block to remove
!> \param block_nze ...
!> \param[in]  block_number   (optional) the block number, if it is known
! *****************************************************************************
  SUBROUTINE dbcsr_remove_block(matrix, row, col, block_nze, block_number)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col, block_nze
    INTEGER, INTENT(IN), OPTIONAL            :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_remove_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, c, error_handle, r
    LOGICAL                                  :: found, tr

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (PRESENT (block_number)) THEN
       b = block_number
       CALL dbcsr_assert (block_number .LE. matrix%m%nblks, dbcsr_failure_level,&
            dbcsr_caller_error, routineN, "Block number too big.",181)
       found = .TRUE.
    ELSE
       CALL dbcsr_get_block_index (matrix, row, col, r, c, tr, found, b)
    ENDIF
    b = ABS (b)
    IF (found .AND. b .GT. 0) THEN
       ! Mark the block for deletion.
       matrix%m%blk_p(b) = 0
       matrix%m%valid = .FALSE.
       ! update nze accordingly
       matrix%m%nze = matrix%m%nze - block_nze
       IF (debug_mod) THEN
          CALL dbcsr_assert (matrix%m%nze, "GE", 0, dbcsr_failure_level,&
             dbcsr_caller_error, routineN, "nze < 0!",195)
       ENDIF
    ELSE
       IF (debug_mod) THEN
          IF(b.EQ.0)&
             CALL cp__w("dbcsr/block/dbcsr_block_access.F",200,"Block does not exist or is already deleted.")
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_remove_block

! *****************************************************************************
!> \brief Gets a block from a dbcsr matrix as a data area
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
!> \par Data area
!>      The pointer encapsulated in the data area points to data stored in the
!>      matrix. It must be 2-dimensional.
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_area(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, error_handle, iw, &
                                                offset, rsize, stored_col, &
                                                stored_row
    LOGICAL                                  :: stored_tr
    TYPE(btree_2d_data_c)                    :: data_block_c
    TYPE(btree_2d_data_d)                    :: data_block_d
    TYPE(btree_2d_data_s)                    :: data_block_s
    TYPE(btree_2d_data_z)                    :: data_block_z

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handle)
    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    CALL dbcsr_data_clear_pointer (block)
    IF(found) THEN
       CALL dbcsr_set_block_pointer (matrix, block,rsize, csize, stored_tr, offset)
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       iw = 1
!$     iw = omp_get_thread_num()+1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",261)
       IF (dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
          SELECT CASE (block%d%data_type)
          CASE (dbcsr_type_real_4_2d)
             CALL btree_get_s (&
                  matrix%m%wms(iw)%mutable%m%btree_s,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_s, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_s%p)
             ENDIF
          CASE (dbcsr_type_real_8_2d)
             CALL btree_get_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_d, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_d%p)
             ENDIF
          CASE (dbcsr_type_complex_4_2d)
             CALL btree_get_c (&
                  matrix%m%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_c, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_c%p)
             ENDIF
          CASE (dbcsr_type_complex_8_2d)
             CALL btree_get_z (&
                  matrix%m%wms(iw)%mutable%m%btree_z,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_z, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_z%p)
             ENDIF
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
                  routineN, "Only 2-D data for block pointers!",298)
          END SELECT
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_block_p_area


! *****************************************************************************
!> \brief
!>          We allow :
!>                  matrix(dp) [+]= [scale(dp)] * block(dp)
!>                  matrix(dp) [+]= [scale(dp)] * block(sp)
!>                  matrix(sp) [+]= [scale(dp)] * block(sp)
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param transposed ...
!> \param summation ...
!> \param scale ...
!> \param[in]
!> \param[out]
!>
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_area(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj)                     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type_m, error_handle
    LOGICAL                                  :: do_scale

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handle)
    data_type_m = dbcsr_get_data_type (matrix)
    do_scale = PRESENT (scale)
    IF (do_scale) THEN
       !CALL dbcsr_assert (data_type_m .EQ. scale%data_type, dbcsr_fatal_level,&
       !     dbcsr_wrong_args_error, routineN, "Incompatible data types matrix="//&
       !     data_type_m//" scale="//scale%data_type)
    ENDIF
    CALL dbcsr_assert (ASSOCIATED (block%d), dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Can only add valid data block!",349)
    SELECT CASE(block%d%data_type)
    CASE (dbcsr_type_real_4)
       IF (do_scale) THEN
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation, scale=scale%r_sp)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation, scale=REAL(scale%r_sp,real_8))
          ENDIF
       ELSE
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation)
          ENDIF
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation, scale=scale%r_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",427)
    END SELECT
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_put_block_area

! *****************************************************************************
!> \brief Inserts all blocks of a dbcsr matrix to make it a full matrix.
!>        Thus obviously not linear scaling.
!> \param[in,out] matrix      Matrix into which blocks should be added.
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_all_blocks(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_all_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_count, col, col_local, &
                                                col_s, error_handle, myrank, &
                                                rank, row, row_local, row_s
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: columns, rows
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows
    LOGICAL                                  :: tr

    CALL timeset (routineN, error_handle)

    myrank = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix)))
    local_rows => dbcsr_distribution_local_rows (dbcsr_distribution (matrix))
    local_cols => dbcsr_distribution_local_cols (dbcsr_distribution (matrix))

    blk_count=0
    ! should be possible to loop only over the local blockrows/blockcols
    DO row_local = 1, SIZE(local_rows)
     DO col_local = 1, SIZE(local_cols)
        tr = .FALSE.
        row = local_rows(row_local)
        col = local_cols(col_local)
        row_s=row ; col_s=col
        CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, rank)
        ! is that the correct condition for symmetric matrices ?
        IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) blk_count=blk_count+1
     ENDDO
    ENDDO

    ALLOCATE(rows(blk_count),columns(blk_count))

    blk_count=0
    DO row_local = 1, SIZE(local_rows)
     DO col_local = 1, SIZE(local_cols)
        tr = .FALSE.
        row = local_rows(row_local)
        col = local_cols(col_local)
        row_s=row ; col_s=col
        CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, rank)
        IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) THEN
           blk_count=blk_count+1 
           rows(blk_count)=row
           columns(blk_count)=col
        ENDIF
     ENDDO
    ENDDO

    CALL dbcsr_reserve_blocks(matrix,rows,columns)

    CALL timestop (error_handle)

  END SUBROUTINE dbcsr_reserve_all_blocks

! *****************************************************************************
!> \brief Inserts diagonal blocks of a dbcsr matrix to make it a matrix with at least all diagonal blocks present
!> \param[in,out] matrix      Matrix into which blocks should be added.
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_diag_blocks(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_diag_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_count, col, col_s, &
                                                myrank, rank, row, row_s
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: columns, rows
    LOGICAL                                  :: tr

    myrank = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix)))

    blk_count=0
    ! should be possible to loop only over the local blockrows/blockcols
    DO row = 1, dbcsr_nblkrows_total(matrix)
       col = row
       tr = .FALSE.
       row_s=row ; col_s=col
       CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, rank)
       IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) blk_count=blk_count+1
    ENDDO

    ALLOCATE(rows(blk_count),columns(blk_count))

    blk_count=0
    DO row = 1, dbcsr_nblkrows_total(matrix)
       col = row
       tr = .FALSE.
       row_s=row ; col_s=col
       CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, rank)
       IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) THEN
          blk_count=blk_count+1
          rows(blk_count)=row
          columns(blk_count)=col
       ENDIF
    ENDDO

    CALL dbcsr_reserve_blocks(matrix,rows,columns)

  END SUBROUTINE dbcsr_reserve_diag_blocks

! *****************************************************************************
!> \brief Inserts block reservations into a matrix, avoiding the work matrix.
!> \param[in,out] matrix      Matrix into which blocks should be added.
!> \param[in] rows            Rows of the blocks to add
!> \param[in] columns         Columns of the blocks to add
!> \param[in] blk_pointers    (optional) block pointers to use for new blocks
!> \par Data
!>      No data can be specified; instead, space is reserved and zeroed. To
!>      add data, call dbcsr_put_block afterwards.
!> \par Reserving existing blocks
!>      Duplicates are not added, but allocations may be greater than
!>      the minimum necessary.
!> \par blk_pointers
!>      When blk_pointers is passed, the newly added blocks use these pointers.
!>      No data is cleared in this case
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_blocks(matrix, rows, columns, blk_pointers)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: rows, columns
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: blk_pointers

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, data_size_new, data_size_old, handle, nblkrows, &
      nblks_actual_added, nblks_added, nblks_new, nblks_old, new_data_sizes, &
      nze
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: add_blkp, add_cols, add_rows, &
                                                added_sizes, new_blk_p, &
                                                new_col_i, new_row_i, &
                                                old_row_i
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: added_blk_info

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    CALL dbcsr_assert (SIZE(rows), "EQ", SIZE(columns), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "Size of rows and columns array must match.", 579)
    IF (PRESENT (blk_pointers)) THEN
       CALL dbcsr_assert (SIZE(rows), "EQ", SIZE(blk_pointers),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Size of rows and block pointecs arrays must match.",&
            584)
       data_size_old = 0
    ELSE
       ! Get current data size
       data_size_old = dbcsr_data_get_size_referenced(matrix%m%data_area)
    ENDIF
    !> Ensures that the rows and columns are sorted.
    nblks_added = SIZE(rows)
    ALLOCATE (add_rows (nblks_added))
    add_rows(:) = rows(:)
    ALLOCATE (add_cols (nblks_added))
    add_cols(:) = columns(:)
    IF (PRESENT(blk_pointers)) THEN
       ALLOCATE (add_blkp (nblks_added))
       add_blkp(:) = blk_pointers(:)
       CALL dbcsr_sort_indices (nblks_added, add_rows, add_cols,&
            blk_p = add_blkp)
    ELSE
       CALL dbcsr_sort_indices (nblks_added, add_rows, add_cols)
    ENDIF
    nblks_old = dbcsr_get_num_blocks (matrix)
    nblkrows = dbcsr_nblkrows_total(matrix)
    CALL dbcsr_assert (SIZE(rows) .GT. 0, "IMP",&
         nblkrows .GT. 0, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN,&
         "Can not add blocks to matrix with no rows.", 609)
    ! Adjust the index.
    ! Get the old row indices
    ALLOCATE (old_row_i (nblks_old))
    CALL dbcsr_expand_row_index(matrix%m%row_p, old_row_i,&
         nblkrows, nblks_old)
    ! Calculate new block pointers. Possibly high estimates.
    new_data_sizes = 0
    blk_p = data_size_old + 1   ! New blocks start at the end of the old
    ALLOCATE (added_blk_info (3, nblks_added))
    ALLOCATE (added_sizes (nblks_added))
    DO blk = 1, nblks_added
       IF (PRESENT (blk_pointers)) THEN
          blk_p = add_blkp(blk)
       ENDIF
       added_blk_info(1:3,blk) = (/ add_rows(blk), add_cols(blk), blk_p /)
       nze = dbcsr_blk_row_size (matrix, add_rows(blk)) &
            * dbcsr_blk_column_size (matrix, add_cols(blk))
       added_sizes(blk) = nze
       blk_p = blk_p + nze
    ENDDO
    DEALLOCATE (add_rows)
    DEALLOCATE (add_cols)
    IF (PRESENT (blk_pointers)) DEALLOCATE (add_blkp)
    !
    nblks_new = nblks_old + nblks_added ! Possibly high estimate
    ALLOCATE (new_row_i (nblks_new))
    ALLOCATE (new_col_i (nblks_new))
    ALLOCATE (new_blk_p (nblks_new))
    ! Merge the two indices
    IF (PRESENT (blk_pointers)) THEN
       CALL merge_index_arrays (new_row_i, new_col_i, new_blk_p, nblks_new,&
            old_row_i, matrix%m%col_i, matrix%m%blk_p, nblks_old,&
            added_blk_info, nblks_added, added_nblks=nblks_actual_added)
       data_size_new = 0
    ELSE
       CALL merge_index_arrays (new_row_i, new_col_i, new_blk_p, nblks_new,&
            old_row_i, matrix%m%col_i, matrix%m%blk_p, nblks_old,&
            added_blk_info, nblks_added, added_nblks=nblks_actual_added,&
            added_sizes=added_sizes, added_size_offset=data_size_old+1,&
            added_size=data_size_new)
    ENDIF
    nblks_new = nblks_actual_added + nblks_old
    ! Free some memory
    DEALLOCATE (added_blk_info)
    DEALLOCATE (added_sizes)
    DEALLOCATE (old_row_i)
    ! We can skip this if no block was actually added.
    IF (nblks_actual_added .GT. 0) THEN
       ! Write the new index
       matrix%m%nblks = nblks_new
       matrix%m%nze = matrix%m%nze + data_size_new
       matrix%m%index(dbcsr_slot_nblks) = matrix%m%nblks
       matrix%m%index(dbcsr_slot_nze) = matrix%m%index(dbcsr_slot_nze)
       CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_col_i)
       CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_blk_p)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
            new_col_i(1:nblks_new),&
            extra=nblks_new)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
            new_blk_p(1:nblks_new))
       CALL dbcsr_make_dbcsr_index (matrix%m%row_p, new_row_i(1:nblks_new),&
            nblkrows, nblks_new)
       IF (.NOT. PRESENT (blk_pointers)) THEN
          ! Resize data area to fit the new blocks.
          CALL dbcsr_data_ensure_size (matrix%m%data_area,&
               data_size = matrix%m%nze)
          ! Zero the new data blocks.
          CALL dbcsr_data_clear (matrix%m%data_area,&
               lb=data_size_old+1, ub=matrix%m%nze)
       ENDIF
    ENDIF
    CALL timestop (handle)
  END SUBROUTINE dbcsr_reserve_blocks


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param[in,out] pointer_any The pointer to set
!> \param[in] rsize Row sizes of block to point to
!> \param[in] csize Column sizes of block to point to
!> \param[in] main_tr          Whether block is transposed in the matrix
!> \param[in] base_offset      The block pointer
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_any (matrix, pointer_any,&
       rsize, csize, main_tr, base_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    LOGICAL, INTENT(IN)                      :: main_tr
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_block_pointer_any', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

       IF (main_tr) THEN
          CALL dbcsr_data_set_pointer (pointer_any, csize, rsize,&
               matrix%m%data_area, source_lb = base_offset)
       ELSE
          CALL dbcsr_data_set_pointer (pointer_any, rsize, csize,&
               matrix%m%data_area, source_lb = base_offset)
       ENDIF
  END SUBROUTINE dbcsr_set_block_pointer_any




# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access_d.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_d', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_d)          :: data_block
    LOGICAL                                  :: stored_tr
    REAL(kind=real_8), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",43)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       IF(nze.eq.0) THEN
          found = .TRUE.
          block => block0(1:0, 1:0)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, rsize, csize, offset)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 71)
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",75)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_d (&
               matrix%m%wms(iw)%mutable%m%btree_d,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_2d_block_p_d


! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",128)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, 0.0_real_8), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",150)
       CALL dbcsr_assert ("NOT", dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.",153)
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_d


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_d(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_d', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_d)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: original_block

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       !$OMP MASTER
       matrix%m%valid = .FALSE.
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    nwms = SIZE(matrix%m%wms)
    iw = 1
!$  CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$     dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "Number of work matrices not equal to number of threads", &
!$     229)
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_d (matrix%m%wms(iw)%mutable%m%btree_d,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
       matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_d

! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr, do_sum

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    IF (PRESENT (scale)) THEN
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_d

! *****************************************************************************
!> \brief Inserts a block in a dbcsr matrix.
!>
!> If the block exists, the current data is overwritten.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the logical row
!> \param[in]  col            the logical column
!> \param[in]  block          the block to put
!> \param[in]  transposed     (optional) the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_d', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_d)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d

!   ---------------------------------------------------------------------------
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col)
    nze = row_size*col_size
    !
    IF (debug_mod) THEN
       CALL dbcsr_assert (SIZE(block), "GE", nze, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Invalid block dimensions",350)
    ENDIF
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL daxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL daxpy (nze, 1.0_real_8, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ELSE
             IF (PRESENT (scale)) THEN
                CALL dcopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL dcopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 400)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=block, source_lb=1)
          ENDIF
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          IF (PRESENT (scale)) THEN
             data_block%p(:,:) = scale*RESHAPE (block, (/row_size, col_size/))
          ELSE
             data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
          ENDIF
          data_block%tr = tr
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                IF(.NOT.ASSOCIATED(data_block2%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",436,"Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                IF(nze > 0) &
                   CALL daxpy (nze, 1.0_real_8, block(1), 1,&
                        data_block2%p(1,1), 1)
                IF(.NOT.ASSOCIATED(data_block%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",449,"Data was not present in block")
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
  END SUBROUTINE dbcsr_put_block_d


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] base_offset      The block pointer
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_d (&
       matrix, pointer_any, rsize, csize, base_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    REAL(kind=real_8), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handler)
    CALL dbcsr_get_data (matrix%m%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL timestop (error_handler)
  END SUBROUTINE dbcsr_set_block_pointer_2d_d
# 719 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access_z.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_z', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_z)          :: data_block
    LOGICAL                                  :: stored_tr
    COMPLEX(kind=real_8), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_complex_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",43)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       IF(nze.eq.0) THEN
          found = .TRUE.
          block => block0(1:0, 1:0)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, rsize, csize, offset)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 71)
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",75)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_z (&
               matrix%m%wms(iw)%mutable%m%btree_z,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_2d_block_p_z


! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_complex_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",128)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",150)
       CALL dbcsr_assert ("NOT", dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.",153)
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_z


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_z(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_z', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_z)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: original_block

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       !$OMP MASTER
       matrix%m%valid = .FALSE.
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    nwms = SIZE(matrix%m%wms)
    iw = 1
!$  CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$     dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "Number of work matrices not equal to number of threads", &
!$     229)
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_z (matrix%m%wms(iw)%mutable%m%btree_z,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
       matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_z

! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_z(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_z', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr, do_sum

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    IF (PRESENT (scale)) THEN
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_z

! *****************************************************************************
!> \brief Inserts a block in a dbcsr matrix.
!>
!> If the block exists, the current data is overwritten.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the logical row
!> \param[in]  col            the logical column
!> \param[in]  block          the block to put
!> \param[in]  transposed     (optional) the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_z(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_z', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_z)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: block_1d

!   ---------------------------------------------------------------------------
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col)
    nze = row_size*col_size
    !
    IF (debug_mod) THEN
       CALL dbcsr_assert (SIZE(block), "GE", nze, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Invalid block dimensions",350)
    ENDIF
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL zaxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL zaxpy (nze, CMPLX(1.0, 0.0, real_8), block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ELSE
             IF (PRESENT (scale)) THEN
                CALL zcopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL zcopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 400)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=block, source_lb=1)
          ENDIF
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          IF (PRESENT (scale)) THEN
             data_block%p(:,:) = scale*RESHAPE (block, (/row_size, col_size/))
          ELSE
             data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
          ENDIF
          data_block%tr = tr
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_z (&
                  matrix%m%wms(iw)%mutable%m%btree_z,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                IF(.NOT.ASSOCIATED(data_block2%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",436,"Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_z (&
                  matrix%m%wms(iw)%mutable%m%btree_z,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                IF(nze > 0) &
                   CALL zaxpy (nze, CMPLX(1.0, 0.0, real_8), block(1), 1,&
                        data_block2%p(1,1), 1)
                IF(.NOT.ASSOCIATED(data_block%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",449,"Data was not present in block")
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
  END SUBROUTINE dbcsr_put_block_z


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] base_offset      The block pointer
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_z (&
       matrix, pointer_any, rsize, csize, base_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handler)
    CALL dbcsr_get_data (matrix%m%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL timestop (error_handler)
  END SUBROUTINE dbcsr_set_block_pointer_2d_z
# 720 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access_s.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_s', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_4), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_s)          :: data_block
    LOGICAL                                  :: stored_tr
    REAL(kind=real_4), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_4,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",43)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       IF(nze.eq.0) THEN
          found = .TRUE.
          block => block0(1:0, 1:0)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_4), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, rsize, csize, offset)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 71)
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",75)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_s (&
               matrix%m%wms(iw)%mutable%m%btree_s,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_2d_block_p_s


! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_4,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",128)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, 0.0_real_4), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",150)
       CALL dbcsr_assert ("NOT", dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.",153)
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_s


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_s(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_s', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_s)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    REAL(kind=real_4), DIMENSION(:,:), POINTER         :: original_block

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       !$OMP MASTER
       matrix%m%valid = .FALSE.
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    nwms = SIZE(matrix%m%wms)
    iw = 1
!$  CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$     dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "Number of work matrices not equal to number of threads", &
!$     229)
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_s (matrix%m%wms(iw)%mutable%m%btree_s,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
       matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_s

! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_s(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_s', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr, do_sum

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    IF (PRESENT (scale)) THEN
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_s

! *****************************************************************************
!> \brief Inserts a block in a dbcsr matrix.
!>
!> If the block exists, the current data is overwritten.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the logical row
!> \param[in]  col            the logical column
!> \param[in]  block          the block to put
!> \param[in]  transposed     (optional) the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_s(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_s', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_s)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    REAL(kind=real_4), DIMENSION(:), POINTER           :: block_1d

!   ---------------------------------------------------------------------------
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col)
    nze = row_size*col_size
    !
    IF (debug_mod) THEN
       CALL dbcsr_assert (SIZE(block), "GE", nze, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Invalid block dimensions",350)
    ENDIF
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, 0.0_real_4), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL saxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL saxpy (nze, 1.0_real_4, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ELSE
             IF (PRESENT (scale)) THEN
                CALL scopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL scopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 400)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=block, source_lb=1)
          ENDIF
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          IF (PRESENT (scale)) THEN
             data_block%p(:,:) = scale*RESHAPE (block, (/row_size, col_size/))
          ELSE
             data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
          ENDIF
          data_block%tr = tr
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_s (&
                  matrix%m%wms(iw)%mutable%m%btree_s,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                IF(.NOT.ASSOCIATED(data_block2%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",436,"Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_s (&
                  matrix%m%wms(iw)%mutable%m%btree_s,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                IF(nze > 0) &
                   CALL saxpy (nze, 1.0_real_4, block(1), 1,&
                        data_block2%p(1,1), 1)
                IF(.NOT.ASSOCIATED(data_block%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",449,"Data was not present in block")
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
  END SUBROUTINE dbcsr_put_block_s


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] base_offset      The block pointer
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_s (&
       matrix, pointer_any, rsize, csize, base_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=real_4), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    REAL(kind=real_4), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handler)
    CALL dbcsr_get_data (matrix%m%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL timestop (error_handler)
  END SUBROUTINE dbcsr_set_block_pointer_2d_s
# 721 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access_c.f90" 1
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_c', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_c)          :: data_block
    LOGICAL                                  :: stored_tr
    COMPLEX(kind=real_4), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_complex_4,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",43)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       IF(nze.eq.0) THEN
          found = .TRUE.
          block => block0(1:0, 1:0)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, rsize, csize, offset)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 71)
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",75)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_c (&
               matrix%m%wms(iw)%mutable%m%btree_c,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_2d_block_p_c


! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_complex_4,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",128)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",150)
       CALL dbcsr_assert ("NOT", dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.",153)
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_c


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_c(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_c', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_c)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: original_block

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       !$OMP MASTER
       matrix%m%valid = .FALSE.
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    nwms = SIZE(matrix%m%wms)
    iw = 1
!$  CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$     dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "Number of work matrices not equal to number of threads", &
!$     229)
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_c (matrix%m%wms(iw)%mutable%m%btree_c,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
       matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_c

! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_c(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_c', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr, do_sum

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    IF (PRESENT (scale)) THEN
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_c

! *****************************************************************************
!> \brief Inserts a block in a dbcsr matrix.
!>
!> If the block exists, the current data is overwritten.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the logical row
!> \param[in]  col            the logical column
!> \param[in]  block          the block to put
!> \param[in]  transposed     (optional) the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_c(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_c', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_c)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block_1d

!   ---------------------------------------------------------------------------
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col)
    nze = row_size*col_size
    !
    IF (debug_mod) THEN
       CALL dbcsr_assert (SIZE(block), "GE", nze, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Invalid block dimensions",350)
    ENDIF
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL caxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL caxpy (nze, CMPLX(1.0, 0.0, real_4), block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ELSE
             IF (PRESENT (scale)) THEN
                CALL ccopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL ccopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", 400)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=block, source_lb=1)
          ENDIF
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          IF (PRESENT (scale)) THEN
             data_block%p(:,:) = scale*RESHAPE (block, (/row_size, col_size/))
          ELSE
             data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
          ENDIF
          data_block%tr = tr
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_c (&
                  matrix%m%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                IF(.NOT.ASSOCIATED(data_block2%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",436,"Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_c (&
                  matrix%m%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                IF(nze > 0) &
                   CALL caxpy (nze, CMPLX(1.0, 0.0, real_4), block(1), 1,&
                        data_block2%p(1,1), 1)
                IF(.NOT.ASSOCIATED(data_block%p))&
                   CALL cp__w("dbcsr/block/dbcsr_block_access.F",449,"Data was not present in block")
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
  END SUBROUTINE dbcsr_put_block_c


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] base_offset      The block pointer
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_c (&
       matrix, pointer_any, rsize, csize, base_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handler)
    CALL dbcsr_get_data (matrix%m%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL timestop (error_handler)
  END SUBROUTINE dbcsr_set_block_pointer_2d_c
# 722 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/dbcsr/block/dbcsr_block_access.F" 2


END MODULE dbcsr_block_access
