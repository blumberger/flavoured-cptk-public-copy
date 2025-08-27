# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_pool_types.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_pool_types.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief pool for for elements that are retained and released
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_pool_types
  USE cp_fm_struct,                    ONLY: cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_linked_list_fm,               ONLY: cp_sll_fm_dealloc,&
                                             cp_sll_fm_get_first_el,&
                                             cp_sll_fm_insert_el,&
                                             cp_sll_fm_next,&
                                             cp_sll_fm_rm_first_el,&
                                             cp_sll_fm_type
  USE cp_log_handling,                 ONLY: cp_to_string

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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/fm/cp_fm_pool_types.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_pool_types'
  INTEGER, SAVE, PRIVATE :: last_fm_pool_id_nr=0

  PUBLIC :: cp_fm_pool_type, cp_fm_pool_p_type
  PUBLIC :: fm_pool_create, fm_pool_retain,&
            fm_pool_release,&
            fm_pool_create_fm, fm_pool_give_back_fm,&
            fm_pool_get_el_struct
  PUBLIC :: fm_pools_dealloc,&
            fm_pools_create_fm_vect,&
            fm_pools_give_back_fm_vect
!***

! *****************************************************************************
!> \brief represent a pool of elements with the same structure
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param el_struct the structure of the elements stored in this pool
!> \param cache linked list with the elements in the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_fm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(cp_fm_struct_type), POINTER :: el_struct
     
     TYPE(cp_sll_fm_type), POINTER :: cache
  END TYPE cp_fm_pool_type

! *****************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_fm_pool_p_type
     TYPE(cp_fm_pool_type), POINTER :: pool
  END TYPE cp_fm_pool_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a pool of elements
!> \param pool the pool to create
!> \param el_struct the structure of the elements that are stored in
!>        this pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_create(pool, el_struct)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_struct_type), POINTER         :: el_struct

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_create', &
      routineP = moduleN//':'//routineN

    ALLOCATE(pool)
    pool%el_struct=> el_struct
    CALL cp_fm_struct_retain(pool%el_struct)
    last_fm_pool_id_nr=last_fm_pool_id_nr+1
    pool%id_nr=last_fm_pool_id_nr
    pool%ref_count=1
    NULLIFY(pool%cache)
    
  END SUBROUTINE fm_pool_create

! *****************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_retain(pool)
    TYPE(cp_fm_pool_type), POINTER           :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_retain', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(pool)))CALL cp__a("fm/cp_fm_pool_types.F",115)
    IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",116)

    pool%ref_count=pool%ref_count+1
  END SUBROUTINE fm_pool_retain

! *****************************************************************************
!> \brief deallocates all the cached elements
!> \param pool the pool to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_flush_cache(pool)
    TYPE(cp_fm_pool_type), POINTER           :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    TYPE(cp_fm_type), POINTER                :: el_att
    TYPE(cp_sll_fm_type), POINTER            :: iterator

    IF(.NOT.(ASSOCIATED(pool)))CALL cp__a("fm/cp_fm_pool_types.F",137)
    IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",138)
    iterator => pool%cache
    DO
       IF (.NOT.cp_sll_fm_next(iterator,el_att=el_att)) EXIT
       CALL cp_fm_release(el_att)
    END DO
    CALL cp_sll_fm_dealloc(pool%cache)
  END SUBROUTINE fm_pool_flush_cache

! *****************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_release(pool)
    TYPE(cp_fm_pool_type), POINTER           :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(pool)) THEN
       IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",161)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count==0) THEN
          pool%ref_count=1
          CALL fm_pool_flush_cache(pool)
          CALL cp_fm_struct_release(pool%el_struct)
          pool%ref_count=0

          DEALLOCATE(pool)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE fm_pool_release

! *****************************************************************************
!> \brief returns an element, allocating it if none is in the pool
!> \param pool the pool from where you get the element
!> \param element will contain the new element
!>\param name the name for the new matrix (optional)
!> \param name ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_create_fm(pool, element,&
       name)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_type), POINTER                :: element
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_create_fm', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(pool)))CALL cp__a("fm/cp_fm_pool_types.F",194)
    IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",195)
    IF (ASSOCIATED(pool%cache)) THEN
       element => cp_sll_fm_get_first_el(pool%cache)
       CALL cp_sll_fm_rm_first_el(pool%cache)
       
    ELSE
       NULLIFY(element)
       CALL cp_fm_create(element,matrix_struct=pool%el_struct)
    END IF
    
 IF (PRESENT(name)) THEN
   element%name=name
   element%print_count=0
 ELSE
   element%name="tmp-"//TRIM(ADJUSTL(cp_to_string(element%id_nr)))
   element%print_count=0
   ! guarantee output unicity?
 END IF
 
    IF(.NOT.(ASSOCIATED(element)))CALL cp__a("fm/cp_fm_pool_types.F",214)
    IF(.NOT.(element%ref_count==1))CALL cp__a("fm/cp_fm_pool_types.F",215)
  END SUBROUTINE fm_pool_create_fm

! *****************************************************************************
!> \brief returns the element to the pool
!> \param pool the pool where to cache the element
!> \param element the element to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      transfers the ownership of the element to the pool
!>      (it is as if you had called cp_fm_release)
!>      Accept give_backs of non associated elements?
! *****************************************************************************
  SUBROUTINE fm_pool_give_back_fm(pool, element)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_type), POINTER                :: element

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_give_back_fm', &
      routineP = moduleN//':'//routineN

    IF(.NOT.(ASSOCIATED(pool)))CALL cp__a("fm/cp_fm_pool_types.F",237)
    IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",238)
    IF(.NOT.(ASSOCIATED(element)))CALL cp__a("fm/cp_fm_pool_types.F",239)
    IF(pool%el_struct%id_nr/=element%matrix_struct%id_nr)&
       CALL cp__b("fm/cp_fm_pool_types.F",241,"pool cannot reuse matrixes with another structure")

    IF(.NOT.(element%ref_count==1))CALL cp__a("fm/cp_fm_pool_types.F",243)
    CALL cp_sll_fm_insert_el(pool%cache, el=element)
    NULLIFY(element)
  END SUBROUTINE fm_pool_give_back_fm

! *****************************************************************************
!> \brief returns the structure of the elements in this pool
!> \param pool the pool you are interested in
!> \retval res ...
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION fm_pool_get_el_struct(pool) RESULT(res)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_struct_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_get_el_struct', &
      routineP = moduleN//':'//routineN

  IF(.NOT.(ASSOCIATED(pool)))CALL cp__a("fm/cp_fm_pool_types.F",263)
  IF(.NOT.(pool%ref_count>0))CALL cp__a("fm/cp_fm_pool_types.F",264)
  res => pool%el_struct
END FUNCTION fm_pool_get_el_struct

!================== pools ================

! *****************************************************************************
!> \brief shallow copy of an array of pools (retains each pool)
!> \param source_pools the pools to copy
!> \param target_pools will contains the new pools
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_copy(source_pools, target_pools)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

  IF(.NOT.(ASSOCIATED(source_pools)))CALL cp__a("fm/cp_fm_pool_types.F",287)
  ALLOCATE(target_pools(SIZE(source_pools)))
  DO i=1,SIZE(source_pools)
     target_pools(i)%pool => source_pools(i)%pool
     CALL fm_pool_retain(source_pools(i)%pool)
  END DO
END SUBROUTINE fm_pools_copy

! *****************************************************************************
!> \brief deallocate an array of pools (releasing each pool)
!> \param pools the pools to release
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_dealloc(pools)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_release(pools(i)%pool)
     END DO
     DEALLOCATE(pools)
  END IF
END SUBROUTINE fm_pools_dealloc

! *****************************************************************************
!> \brief Returns a vector with an element from each pool
!> \param pools the pools to create the elements from
!> \param elements will contain the vector of elements
!> \param name the name for the new matrixes (optional)
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_create_fm_vect(pools,elements,&
     name)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_create_fm_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(cp_fm_pool_type), POINTER           :: pool

  NULLIFY(pool)

  IF(.NOT.(ASSOCIATED(pools)))CALL cp__a("fm/cp_fm_pool_types.F",344)
  ALLOCATE(elements(SIZE(pools)))
  DO i=1,SIZE(pools)
     NULLIFY(elements(i)%matrix)
     pool => pools(i)%pool
     IF (PRESENT(name)) THEN
      CALL fm_pool_create_fm(pool,elements(i)%matrix,&
        name=name//"-"//ADJUSTL(cp_to_string(i)))
   ELSE
      CALL fm_pool_create_fm(pool,elements(i)%matrix)
   END IF

  END DO
  
END SUBROUTINE fm_pools_create_fm_vect

! *****************************************************************************
!> \brief returns a vector to the pools. The vector is deallocated
!>      (like cp_fm_vect_dealloc)
!> \param pools the pool where to give back the vector
!> \param elements the vector of elements to give back
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      accept unassociated vect?
! *****************************************************************************
SUBROUTINE fm_pools_give_back_fm_vect(pools,elements)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_give_back_fm_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

  IF(.NOT.(ASSOCIATED(pools)))CALL cp__a("fm/cp_fm_pool_types.F",382)
  IF(.NOT.(ASSOCIATED(elements)))CALL cp__a("fm/cp_fm_pool_types.F",383)
  IF(.NOT.(SIZE(pools)==SIZE(elements)))CALL cp__a("fm/cp_fm_pool_types.F",384)
  DO i=1,SIZE(pools)
     CALL fm_pool_give_back_fm(pools(i)%pool,&
          elements(i)%matrix)
  END DO
  DEALLOCATE(elements)
END SUBROUTINE fm_pools_give_back_fm_vect

END MODULE cp_fm_pool_types
