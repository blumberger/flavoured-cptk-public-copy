# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/list_timerenv.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/list_timerenv.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief An array-based list which grows on demand.
!>        When the internal array is full, a new array of twice the size will be
!>        allocated and the items are copied over.
!>
!>        This list can also be used as a stack.
!>        Have look at list_push(), list_pop() and list_peek().
!> \note
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify list__valuetype_.template instead
!>
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! *****************************************************************************

MODULE list_timerenv

  USE timings_types,                   ONLY: timer_env_type

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/../base/base_uses.f90" 1
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
# 28 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/common/list_timerenv.F" 2

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'list_timerenv'

  !this is an internal type
  TYPE private_item_type
     PRIVATE
     TYPE(timer_env_type), POINTER                   :: value
  END TYPE private_item_type

  !this is an internal type
  TYPE private_item_p_type
     PRIVATE
     TYPE(private_item_type), POINTER :: p => Null()
  END TYPE private_item_p_type

  ! this is the public type, which holds a list-instance
  TYPE list_timerenv_type
     PRIVATE
     TYPE(private_item_p_type), DIMENSION(:), POINTER   :: arr => Null()
     INTEGER                                       :: size = -1
  END TYPE list_timerenv_type


  PUBLIC :: list_timerenv_type
  PUBLIC :: list_timerenv_init
  PUBLIC :: list_timerenv_push
  PUBLIC :: list_timerenv_pop
  PUBLIC :: list_timerenv_peek
  PUBLIC :: list_timerenv_insert
  PUBLIC :: list_timerenv_set
  PUBLIC :: list_timerenv_get
  PUBLIC :: list_timerenv_del
  PUBLIC :: list_timerenv_clear
  PUBLIC :: list_timerenv_size
  PUBLIC :: list_timerenv_destroy
  PUBLIC :: list_timerenv_isready


  CONTAINS


! *****************************************************************************
!> \brief Test if the given list has been initialized.
!> \param list ...
!> \retval res ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION list_timerenv_isready(list) RESULT(res)
    TYPE(list_timerenv_type), INTENT(in)     :: list
    LOGICAL                                  :: res

    res = ASSOCIATED(list%arr)
END FUNCTION list_timerenv_isready


! *****************************************************************************
!> \brief Allocates the internal data-structures of the given list.
!>        This has to be called before any of the other routines.
!>        For deallocation call list_timerenv_destroy.
!> \param list ...
!> \param initial_capacity The initial size of the internal array (default=11).
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_init(list, initial_capacity)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    INTEGER, INTENT(in), OPTIONAL            :: initial_capacity

    INTEGER                                  :: initial_capacity_, stat

    initial_capacity_ = 11
    IF(PRESENT(initial_capacity)) initial_capacity_ = initial_capacity

    IF(initial_capacity_ < 0) & 
      CALL cp__b("common/list_timerenv.F",109,"list_timerenv_create: initial_capacity < 0")

    IF(ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",112,"list_timerenv_create: list is already initialized.")

    ALLOCATE(list%arr(initial_capacity_), stat=stat)
    IF (stat/=0)&
       CALL cp__b("common/list_timerenv.F",116,"list_timerenv_init: allocation failed")

    list%size = 0
END SUBROUTINE list_timerenv_init

! *****************************************************************************
!> \brief Deallocated the internal data-structures of the given list.
!>        Caution: If the stored values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_destroy(list)
    TYPE(list_timerenv_type), INTENT(inout)  :: list

    INTEGER                                  :: i

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",136,"list_timerenv_destroy: list is not initialized.")

    DO i=1, list%size
       DEALLOCATE(list%arr(i)%p)
    END DO
    DEALLOCATE(list%arr)
    list%size = -1
END SUBROUTINE list_timerenv_destroy


! *****************************************************************************
!> \brief Assings the given value to the given position in the list. 
!>        Thereby, the former value at that position gets overwritten.
!>        If the position is out of bounds, the program stops.
!> \param list ...
!> \param value ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1.
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_set(list, value, pos)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    TYPE(timer_env_type), INTENT(in), &
      POINTER                                :: value
    INTEGER, INTENT(in)                      :: pos

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",164,"list_timerenv_set: list is not initialized.")
    IF(pos < 1)&
       CALL cp__b("common/list_timerenv.F",166,"list_timerenv_set: pos < 1")
    IF(pos > list%size)&
       CALL cp__b("common/list_timerenv.F",168,"list_timerenv_set: pos > size")
    list%arr(pos)%p%value => value
END SUBROUTINE list_timerenv_set


! *****************************************************************************
!> \brief Appends the given value at the end of the list.
!> \param list ...
!> \param value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_push(list, value)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    TYPE(timer_env_type), INTENT(in), &
      POINTER                                :: value

    INTEGER                                  :: stat

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",189,"list_timerenv_push: list is not initialized.")
    IF(list%size == SIZE(list%arr)) &
       CALL change_capacity(list, 2*SIZE(list%arr)+1)

    list%size = list%size + 1
    ALLOCATE(list%arr(list%size)%p, stat=stat)
    IF (stat/=0)&
       CALL cp__b("common/list_timerenv.F",196,"list_timerenv_push: allocation failed")
    list%arr(list%size)%p%value => value
END SUBROUTINE list_timerenv_push


! *****************************************************************************
!> \brief Inserts the given value at the givenn position within the list.
!>        Values which lay behind the insertion-position move one position up.
!> \param list ...
!> \param value ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+2 .
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_insert(list, value, pos)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    TYPE(timer_env_type), INTENT(in), &
      POINTER                                :: value
    INTEGER, INTENT(in)                      :: pos

    INTEGER                                  :: i, stat

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",220,"list_timerenv_insert: list is not initialized.")
    IF(pos < 1)&
       CALL cp__b("common/list_timerenv.F",222,"list_timerenv_insert: pos < 1")
    IF(pos > list%size+1)&
       CALL cp__b("common/list_timerenv.F",224,"list_timerenv_insert: pos > size+1")

    IF(list%size == SIZE(list%arr)) &
       CALL change_capacity(list, 2*SIZE(list%arr)+1)

    list%size = list%size + 1
    DO i=list%size, pos+1, -1
       list%arr(i)%p => list%arr(i-1)%p
    END DO

    ALLOCATE(list%arr(pos)%p, stat=stat)
     IF (stat/=0)&
        CALL cp__b("common/list_timerenv.F",236,"list_timerenv_insert: allocation failed.")
    list%arr(pos)%p%value => value
END SUBROUTINE list_timerenv_insert


! *****************************************************************************
!> \brief Returns the last element in the list.
!>    Is equivalent to: list_timerenv_get(list, list_timerenv_size(list))
!> \param list ...
!> \retval value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION list_timerenv_peek(list) RESULT(value)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    TYPE(timer_env_type), POINTER            :: value

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",255,"list_timerenv_peek: list is not initialized.")
    IF(list%size < 1) &
       CALL cp__b("common/list_timerenv.F",257,"list_timerenv_peek: list is empty.")

    value => list%arr(list%size)%p%value
END FUNCTION list_timerenv_peek


! *****************************************************************************
!> \brief Returns the last element in the list and removes it.
!>        Is equivialent to:
!>        value = list_timerenv_get(list, list_timerenv_size(list))
!>            call list_timerenv_del(list, list_timerenv_size(list))
!>
!> \param list ...
!> \retval value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION list_timerenv_pop(list) RESULT(value)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    TYPE(timer_env_type), POINTER            :: value

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",280,"list_timerenv_pop: list is not initialized.")
    IF(list%size < 1) &
       CALL cp__b("common/list_timerenv.F",282,"list_timerenv_pop: list is empty.")

    value => list%arr(list%size)%p%value
    DEALLOCATE(list%arr(list%size)%p)
    list%size = list%size - 1
END FUNCTION list_timerenv_pop



! *****************************************************************************
!> \brief Removes all values from the list. The list itself is not deallocated.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_clear(list)
    TYPE(list_timerenv_type), INTENT(inout)  :: list

    INTEGER                                  :: i

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",304,"list_timerenv_clear: list is not initialized.")

    DO i=1, list%size
       DEALLOCATE(list%arr(i)%p)
    END DO
    list%size = 0
END SUBROUTINE list_timerenv_clear


!
! *****************************************************************************
!> \brief Returns the value at the given position from the list.
!> \param list ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1 .
!> \retval value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION list_timerenv_get(list, pos) RESULT(value)
    TYPE(list_timerenv_type), INTENT(in)     :: list
    INTEGER, INTENT(in)                      :: pos
    TYPE(timer_env_type), POINTER            :: value

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",329,"list_timerenv_get: list is not initialized.")
    IF(pos < 1)&
       CALL cp__b("common/list_timerenv.F",331,"list_timerenv_get: pos < 1")
    IF(pos > list%size)&
       CALL cp__b("common/list_timerenv.F",333,"list_timerenv_get: pos > size")

    value => list%arr(pos)%p%value

END FUNCTION list_timerenv_get


! *****************************************************************************
!> \brief Removes the value at the given position from the list.
!> \param list ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1 .
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE list_timerenv_del(list, pos)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    INTEGER, INTENT(in)                      :: pos

    INTEGER                                  :: i

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",355,"list_timerenv_del: list is not initialized.")
    IF(pos < 1)&
       CALL cp__b("common/list_timerenv.F",357,"list_timerenv_det: pos < 1")
    IF(pos > list%size)&
       CALL cp__b("common/list_timerenv.F",359,"list_timerenv_det: pos > size")

    DEALLOCATE(list%arr(pos)%p)
    DO i=pos, list%size-1
       list%arr(i)%p => list%arr(i+1)%p 
    END DO

    list%size = list%size - 1

END SUBROUTINE list_timerenv_del

! *****************************************************************************
!> \brief Returns the current size of the list.
!> \param list ...
!> \retval size ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION list_timerenv_size(list) RESULT(size)
    TYPE(list_timerenv_type), INTENT(in)     :: list
    INTEGER                                  :: size

    IF(.NOT. ASSOCIATED(list%arr)) &
       CALL cp__b("common/list_timerenv.F",383,"list_timerenv_size: list is not initialized.")

    size = list%size
END FUNCTION list_timerenv_size


! *****************************************************************************
!> \brief Internal routine for changing the size of the internal array.
!> \param list ...
!> \param new_capacity ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE change_capacity(list, new_capacity)
    TYPE(list_timerenv_type), INTENT(inout)  :: list
    INTEGER, INTENT(in)                      :: new_capacity

    INTEGER                                  :: i, new_cap, stat
    TYPE(private_item_p_type), &
      DIMENSION(:), POINTER                  :: old_arr

   new_cap = new_capacity
   IF(new_cap < 0) &
      CALL cp__b("common/list_timerenv.F",407,"list_timerenv_change_capacity: new_capacity < 0")
   IF(new_cap < list%size) &
      CALL cp__b("common/list_timerenv.F",409,"list_timerenv_change_capacity: new_capacity < size")
   IF(new_cap > HUGE(i)) THEN
      IF(SIZE(list%arr) == HUGE(i)) &
      CALL cp__b("common/list_timerenv.F",412,"list_timerenv_change_capacity: list has reached integer limit.")
      new_cap = HUGE(i) ! grow as far as possible
   END IF

   old_arr => list%arr
   ALLOCATE(list%arr(new_cap), stat=stat)
   IF (stat/=0)&
      CALL cp__b("common/list_timerenv.F",419,"list_timerenv_change_capacity: allocation failed")

   DO i=1, list%size
      ALLOCATE(list%arr(i)%p, stat=stat)
      IF (stat/=0)&
         CALL cp__b("common/list_timerenv.F",424,"list_timerenv_change_capacity: allocation failed")
      list%arr(i)%p%value => old_arr(i)%p%value
      DEALLOCATE(old_arr(i)%p)
   END DO
   DEALLOCATE(old_arr)

END SUBROUTINE change_capacity

END MODULE list_timerenv
