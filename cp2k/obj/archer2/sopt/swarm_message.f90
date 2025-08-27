# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F"
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Swarm-message, a convenient data-container for with build-in serialization.
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_message

  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE kinds,                           ONLY: default_string_length,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_environ,&
                                             mp_recv,&
                                             mp_send

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/../base/base_uses.f90" 1
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
# 24 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_message'

 TYPE swarm_message_type
   PRIVATE
   TYPE(message_entry_type), POINTER :: root => Null()
 END TYPE swarm_message_type

 INTEGER, PARAMETER  :: key_length      = 20

 TYPE message_entry_type
   CHARACTER(LEN=key_length)                        :: key
   TYPE(message_entry_type),            POINTER     :: next        => Null()
   CHARACTER(LEN=default_string_length),POINTER     :: value_str   => Null()
   INTEGER(KIND=int_4),                 POINTER     :: value_i4    => Null()
   INTEGER(KIND=int_8),                 POINTER     :: value_i8    => Null()
   REAL(KIND=real_4),                   POINTER     :: value_r4    => Null()
   REAL(KIND=real_8),                   POINTER     :: value_r8    => Null()
   INTEGER(KIND=int_4), DIMENSION(:),   POINTER     :: value_1d_i4 => Null()
   INTEGER(KIND=int_8), DIMENSION(:),   POINTER     :: value_1d_i8 => Null()
   REAL(KIND=real_4),   DIMENSION(:),   POINTER     :: value_1d_r4 => Null()
   REAL(KIND=real_8),   DIMENSION(:),   POINTER     :: value_1d_r8 => Null()
 END TYPE message_entry_type


! *****************************************************************************
!> \brief Adds an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 INTERFACE swarm_message_add
    MODULE PROCEDURE swarm_message_add_str
    MODULE PROCEDURE swarm_message_add_i4,    swarm_message_add_i8
    MODULE PROCEDURE swarm_message_add_r4,    swarm_message_add_r8
    MODULE PROCEDURE swarm_message_add_1d_i4, swarm_message_add_1d_i8
    MODULE PROCEDURE swarm_message_add_1d_r4, swarm_message_add_1d_r8
 END INTERFACE swarm_message_add


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 INTERFACE swarm_message_get
    MODULE PROCEDURE swarm_message_get_str
    MODULE PROCEDURE swarm_message_get_i4,    swarm_message_get_i8
    MODULE PROCEDURE swarm_message_get_r4,    swarm_message_get_r8
    MODULE PROCEDURE swarm_message_get_1d_i4, swarm_message_get_1d_i8
    MODULE PROCEDURE swarm_message_get_1d_r4, swarm_message_get_1d_r8
 END INTERFACE swarm_message_get


  PUBLIC :: swarm_message_type, swarm_message_add, swarm_message_get
  PUBLIC :: swarm_message_mpi_send, swarm_message_mpi_recv, swarm_message_mpi_bcast
  PUBLIC :: swarm_message_file_write, swarm_message_file_read
  PUBLIC :: swarm_message_haskey, swarm_message_equal
  PUBLIC :: swarm_message_free


 CONTAINS


! *****************************************************************************
!> \brief Returns the number of entries contained in a swarm-message.
!> \param msg ...
!> \retval l ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_length(msg) RESULT(l)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER                                  :: l

    TYPE(message_entry_type), POINTER        :: curr_entry

    l = 0
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      l = l + 1
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_length


! *****************************************************************************
!> \brief Checks if a swarm-message contains an entry with the given key.
!> \param msg ...
!> \param key ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_haskey(msg, key) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    CHARACTER(LEN=*), INTENT(IN)             :: key
    LOGICAL                                  :: res

    TYPE(message_entry_type), POINTER        :: curr_entry

    res = .FALSE.
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         res = .TRUE.
         EXIT
      END IF
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_haskey


! *****************************************************************************
!> \brief Deallocates all entries contained in a swarm-message.
!> \param msg ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_free(msg)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_free', &
      routineP = moduleN//':'//routineN

    TYPE(message_entry_type), POINTER        :: ENTRY, old_entry

     ENTRY => msg%root
     DO WHILE(ASSOCIATED(ENTRY))
       IF(ASSOCIATED(entry%value_str))      DEALLOCATE(entry%value_str)
       IF(ASSOCIATED(entry%value_i4))       DEALLOCATE(entry%value_i4)
       IF(ASSOCIATED(entry%value_i8))       DEALLOCATE(entry%value_i8)
       IF(ASSOCIATED(entry%value_r4))       DEALLOCATE(entry%value_r4)
       IF(ASSOCIATED(entry%value_r8))       DEALLOCATE(entry%value_r8)
       IF(ASSOCIATED(entry%value_1d_i4))    DEALLOCATE(entry%value_1d_i4)
       IF(ASSOCIATED(entry%value_1d_i8))    DEALLOCATE(entry%value_1d_i8)
       IF(ASSOCIATED(entry%value_1d_r4))    DEALLOCATE(entry%value_1d_r4)
       IF(ASSOCIATED(entry%value_1d_r8))    DEALLOCATE(entry%value_1d_r8)
       old_entry => ENTRY
       ENTRY => entry%next
       DEALLOCATE(old_entry)
     END DO

     NULLIFY(msg%root)

     IF(.NOT.(swarm_message_length(msg)==0))CALL cp__a("swarm/swarm_message.F",166)
  END SUBROUTINE swarm_message_free


! *****************************************************************************
!> \brief Checks if two swarm-messages are equal
!> \param msg1 ...
!> \param msg2 ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_equal(msg1, msg2) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg1, msg2
    LOGICAL                                  :: res

     res =  swarm_message_equal_oneway(msg1, msg2) .AND. &
            swarm_message_equal_oneway(msg2, msg1)

  END FUNCTION swarm_message_equal


! *****************************************************************************
!> \brief Sends a swarm message via MPI.
!> \param msg ...
!> \param group ...
!> \param dest ...
!> \param tag ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_send(msg, group, dest, tag)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER, INTENT(IN)                      :: group, dest, tag

    TYPE(message_entry_type), POINTER        :: curr_entry

    CALL mp_send(swarm_message_length(msg), dest, tag, group)
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_mpi_send(curr_entry, group, dest, tag)
      curr_entry => curr_entry%next
    END DO
  END SUBROUTINE swarm_message_mpi_send


! *****************************************************************************
!> \brief Receives a swarm message via MPI.
!> \param msg ...
!> \param group ...
!> \param src ...
!> \param tag ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_recv(msg, group, src, tag)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(INOUT)                   :: src, tag

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_mpi_recv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, length
    TYPE(message_entry_type), POINTER        :: new_entry

    IF(ASSOCIATED(msg%root)) CALL cp__b("swarm/swarm_message.F",229,"message not empty")
    CALL mp_recv(length, src, tag, group)
    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_mpi_recv(new_entry, group, src, tag)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

  END SUBROUTINE swarm_message_mpi_recv


! *****************************************************************************
!> \brief Broadcasts a swarm message via MPI.
!> \param msg ...
!> \param src ...
!> \param group ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_bcast(msg, src, group)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg
    INTEGER, INTENT(IN)                      :: src, group

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_mpi_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, length, mepos, num_pe
    TYPE(message_entry_type), POINTER        :: curr_entry

    CALL mp_environ(num_pe, mepos, group)

    IF(mepos/=src .AND. ASSOCIATED(msg%root)) CALL cp__b("swarm/swarm_message.F",260,"message not empty")
    length = swarm_message_length(msg)
    CALL mp_bcast(length, src, group)

    IF(mepos==src) curr_entry => msg%root

    DO i=1, length
       IF(mepos/=src) ALLOCATE(curr_entry)

       CALL swarm_message_entry_mpi_bcast(curr_entry, src, group, mepos)

       IF(mepos==src) THEN
          curr_entry => curr_entry%next
       ELSE
          curr_entry%next => msg%root
          msg%root => curr_entry
       END IF
    END DO

  END SUBROUTINE swarm_message_mpi_bcast



! *****************************************************************************
!> \brief Write a swarm-message to a given file / unit.
!> \param msg ...
!> \param unit ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_write(msg, unit)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER, INTENT(IN)                      :: unit

    INTEGER                                  :: handle
    TYPE(message_entry_type), POINTER        :: curr_entry

    IF(unit <= 0) RETURN

    CALL timeset("swarm_message_file_write", handle)
    WRITE(unit,"(A)") "BEGIN SWARM_MESSAGE"
    WRITE(unit,"(A,I10)") "msg_length: ", swarm_message_length(msg)

    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_file_write(curr_entry, unit)
      curr_entry => curr_entry%next
    END DO

    WRITE(unit,"(A)") "END SWARM_MESSAGE"
    WRITE(unit,"()")
    CALL timestop(handle)
  END SUBROUTINE swarm_message_file_write


! *****************************************************************************
!> \brief Reads a swarm-message from a given file / unit.
!> \param msg ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_read(msg, parser, at_end)
    TYPE(swarm_message_type), INTENT(OUT)    :: msg
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end

    INTEGER                                  :: handle

    CALL timeset("swarm_message_file_read", handle)
    CALL swarm_message_file_read_low(msg, parser, at_end)
    CALL timestop(handle)
  END SUBROUTINE swarm_message_file_read


! *****************************************************************************
!> \brief Helper routine, does the actual work of swarm_message_file_read().
!> \param msg ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_read_low(msg, parser, at_end)
    TYPE(swarm_message_type), INTENT(OUT)    :: msg
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_file_read_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: label
    INTEGER                                  :: i, length
    TYPE(message_entry_type), POINTER        :: new_entry

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    IF(.NOT.(TRIM(parser%input_line(1:20))=="BEGIN SWARM_MESSAGE"))CALL cp__a("swarm/swarm_message.F",356)

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:40),*) label, length
    IF(.NOT.(TRIM(label)=="msg_length:"))CALL cp__a("swarm/swarm_message.F",362)

    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_file_read(new_entry, parser, at_end)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    IF(.NOT.(TRIM(parser%input_line(1:20))=="END SWARM_MESSAGE"))CALL cp__a("swarm/swarm_message.F",374)

  END SUBROUTINE swarm_message_file_read_low


! *****************************************************************************
!> \brief Helper routine for swarm_message_equal
!> \param msg1 ...
!> \param msg2 ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_equal_oneway(msg1, msg2) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg1, msg2
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_equal_oneway', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: found
    TYPE(message_entry_type), POINTER        :: entry1, entry2

     res = .FALSE.

     !loop over entries of msg1
     entry1 => msg1%root
     DO WHILE(ASSOCIATED(entry1))

        ! finding matching entry in msg2
        entry2 => msg2%root
        found = .FALSE.
        DO WHILE(ASSOCIATED(entry2))
           IF(TRIM(entry2%key) == TRIM(entry1%key)) THEN
              found = .TRUE.
              EXIT
           END IF
           entry2 => entry2%next
        END DO
        IF(.NOT. found) RETURN

        !compare the two entries
        IF(ASSOCIATED(entry1%value_str)) THEN
           IF(.NOT.ASSOCIATED(entry2%value_str)) RETURN
           IF(TRIM(entry1%value_str) /= TRIM(entry2%value_str)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i4)) RETURN
           IF(entry1%value_i4 /= entry2%value_i4) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i8)) RETURN
           IF(entry1%value_i8 /= entry2%value_i8) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r4)) RETURN
           IF(ABS(entry1%value_r4-entry2%value_r4)>1e-5) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r8)) RETURN
           IF(ABS(entry1%value_r8-entry2%value_r8)>1e-10) RETURN

        ELSE IF(ASSOCIATED(entry1%value_1d_i4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_1d_i4)) RETURN
           IF(ANY(entry1%value_1d_i4 /= entry2%value_1d_i4)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_1d_i8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_1d_i8)) RETURN
           IF(ANY(entry1%value_1d_i8 /= entry2%value_1d_i8)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_1d_r4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_1d_r4)) RETURN
           IF(ANY(ABS(entry1%value_1d_r4-entry2%value_1d_r4)>1e-5)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_1d_r8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_1d_r8)) RETURN
           IF(ANY(ABS(entry1%value_1d_r8-entry2%value_1d_r8)>1e-10)) RETURN
        ELSE
           CALL cp__b("swarm/swarm_message.F",451,"no value ASSOCIATED")
        END IF

        entry1 => entry1%next
     END DO

     ! if we reach this point no differences were found
     res = .TRUE.
  END FUNCTION swarm_message_equal_oneway


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_send.
!> \param ENTRY ...
!> \param group ...
!> \param dest ...
!> \param tag ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_send(ENTRY, group, dest, tag)
    TYPE(message_entry_type), INTENT(IN)     :: ENTRY
    INTEGER, INTENT(IN)                      :: group, dest, tag

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_entry_mpi_send', &
      routineP = moduleN//':'//routineN

    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    key_arr = str2iarr(entry%key)
    CALL mp_send(key_arr, dest, tag, group)

    IF(ASSOCIATED(entry%value_i4)) THEN
       CALL mp_send(1, dest, tag, group)
       CALL mp_send(entry%value_i4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       CALL mp_send(2, dest, tag, group)
       CALL mp_send(entry%value_i8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       CALL mp_send(3, dest, tag, group)
       CALL mp_send(entry%value_r4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       CALL mp_send(4, dest, tag, group)
       CALL mp_send(entry%value_r8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_1d_i4)) THEN
       CALL mp_send(5, dest, tag, group)
       CALL mp_send(SIZE(entry%value_1d_i4), dest, tag, group)
       CALL mp_send(entry%value_1d_i4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_1d_i8)) THEN
       CALL mp_send(6, dest, tag, group)
       CALL mp_send(SIZE(entry%value_1d_i8), dest, tag, group)
       CALL mp_send(entry%value_1d_i8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_1d_r4)) THEN
       CALL mp_send(7, dest, tag, group)
       CALL mp_send(SIZE(entry%value_1d_r4), dest, tag, group)
       CALL mp_send(entry%value_1d_r4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_1d_r8)) THEN
       CALL mp_send(8, dest, tag, group)
       CALL mp_send(SIZE(entry%value_1d_r8), dest, tag, group)
       CALL mp_send(entry%value_1d_r8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       CALL mp_send(9, dest, tag, group)
       value_str_arr = str2iarr(entry%value_str)
       CALL mp_send(value_str_arr, dest, tag, group)
    ELSE
       CALL cp__b("swarm/swarm_message.F",525,"no value ASSOCIATED")
    END IF
  END SUBROUTINE swarm_message_entry_mpi_send


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_recv.
!> \param ENTRY ...
!> \param group ...
!> \param src ...
!> \param tag ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_recv(ENTRY, group, src, tag)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(INOUT)                   :: src, tag

    CHARACTER(len=*), PARAMETER :: routineN = 'swarm_message_entry_mpi_recv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: datatype, s
    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    CALL mp_recv(key_arr, src, tag, group)
    entry%key = iarr2str(key_arr)

    CALL mp_recv(datatype, src, tag, group)

    SELECT CASE(datatype)
    CASE(1)
       ALLOCATE(entry%value_i4)
       CALL mp_recv(entry%value_i4, src, tag, group)
    CASE(2)
       ALLOCATE(entry%value_i8)
       CALL mp_recv(entry%value_i8, src, tag, group)
    CASE(3)
       ALLOCATE(entry%value_r4)
       CALL mp_recv(entry%value_r4, src, tag, group)
    CASE(4)
       ALLOCATE(entry%value_r8)
       CALL mp_recv(entry%value_r8, src, tag, group)

    CASE(5)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_1d_i4(s))
       CALL mp_recv(entry%value_1d_i4, src, tag, group)
    CASE(6)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_1d_i8(s))
       CALL mp_recv(entry%value_1d_i8, src, tag, group)
    CASE(7)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_1d_r4(s))
       CALL mp_recv(entry%value_1d_r4, src, tag, group)
    CASE(8)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_1d_r8(s))
       CALL mp_recv(entry%value_1d_r8, src, tag, group)
    CASE(9)
       ALLOCATE(entry%value_str)
       CALL mp_recv(value_str_arr, src, tag, group)
       entry%value_str = iarr2str(value_str_arr)
    CASE DEFAULT
       CALL cp__b("swarm/swarm_message.F",591,"unknown datatype")
    END SELECT
  END SUBROUTINE swarm_message_entry_mpi_recv


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_bcast.
!> \param ENTRY ...
!> \param src ...
!> \param group ...
!> \param mepos ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_bcast(ENTRY, src, group, mepos)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    INTEGER, INTENT(IN)                      :: src, group, mepos

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'swarm_message_entry_mpi_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: datasize, datatype
    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    IF(src==mepos) key_arr = str2iarr(entry%key)
    CALL mp_bcast(key_arr, src, group)
    IF(src/=mepos) entry%key = iarr2str(key_arr)

    IF(src==mepos) THEN
       datasize = 1
       IF(ASSOCIATED(entry%value_i4)) THEN
          datatype = 1
       ELSE IF(ASSOCIATED(entry%value_i8)) THEN
          datatype = 2
       ELSE IF(ASSOCIATED(entry%value_r4)) THEN
          datatype = 3
       ELSE IF(ASSOCIATED(entry%value_r8)) THEN
          datatype = 4
       ELSE IF(ASSOCIATED(entry%value_1d_i4)) THEN
          datatype = 5
          datasize = SIZE(entry%value_1d_i4)
       ELSE IF(ASSOCIATED(entry%value_1d_i8)) THEN
          datatype = 6
          datasize = SIZE(entry%value_1d_i8)
       ELSE IF(ASSOCIATED(entry%value_1d_r4)) THEN
          datatype = 7
          datasize = SIZE(entry%value_1d_r4)
       ELSE IF(ASSOCIATED(entry%value_1d_r8)) THEN
          datatype = 8
          datasize = SIZE(entry%value_1d_r8)
       ELSE IF(ASSOCIATED(entry%value_str)) THEN
          datatype = 9
       ELSE
          CALL cp__b("swarm/swarm_message.F",646,"no value ASSOCIATED")
       END IF
    END IF
    CALL mp_bcast(datatype, src, group)
    CALL mp_bcast(datasize, src, group)

    SELECT CASE(datatype)
    CASE(1)
       IF(src/=mepos) ALLOCATE(entry%value_i4)
       CALL mp_bcast(entry%value_i4, src, group)
    CASE(2)
       IF(src/=mepos) ALLOCATE(entry%value_i8)
       CALL mp_bcast(entry%value_i8, src, group)
    CASE(3)
       IF(src/=mepos) ALLOCATE(entry%value_r4)
       CALL mp_bcast(entry%value_r4, src, group)
    CASE(4)
       IF(src/=mepos) ALLOCATE(entry%value_r8)
       CALL mp_bcast(entry%value_r8, src, group)
    CASE(5)
       IF(src/=mepos) ALLOCATE(entry%value_1d_i4(datasize))
       CALL mp_bcast(entry%value_1d_i4, src, group)
    CASE(6)
       IF(src/=mepos) ALLOCATE(entry%value_1d_i8(datasize))
       CALL mp_bcast(entry%value_1d_i8, src, group)
    CASE(7)
       IF(src/=mepos) ALLOCATE(entry%value_1d_r4(datasize))
       CALL mp_bcast(entry%value_1d_r4, src, group)
    CASE(8)
       IF(src/=mepos) ALLOCATE(entry%value_1d_r8(datasize))
       CALL mp_bcast(entry%value_1d_r8, src, group)
    CASE(9)
       IF(src==mepos) value_str_arr = str2iarr(entry%value_str)
       CALL mp_bcast(value_str_arr, src, group)
       IF(src/=mepos) THEN
          ALLOCATE(entry%value_str)
          entry%value_str = iarr2str(value_str_arr)
       END IF
    CASE DEFAULT
       CALL cp__b("swarm/swarm_message.F",685,"unknown datatype")
    END SELECT

  END SUBROUTINE swarm_message_entry_mpi_bcast


! *****************************************************************************
!> \brief Helper routine for swarm_message_file_write.
!> \param ENTRY ...
!> \param unit ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_write(ENTRY, unit)
    TYPE(message_entry_type), INTENT(IN)     :: ENTRY
    INTEGER, INTENT(IN)                      :: unit

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'swarm_message_entry_file_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    WRITE(unit,"(A,A)") "key: ",  entry%key
    IF(ASSOCIATED(entry%value_i4)) THEN
       WRITE(unit,"(A)") "datatype: i4"
       WRITE(unit,"(A,I10)") "value: ", entry%value_i4

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       WRITE(unit,"(A)") "datatype: i8"
       WRITE(unit,"(A,I20)") "value: ", entry%value_i8

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       WRITE(unit,"(A)") "datatype: r4"
       WRITE(unit,"(A,E30.20)") "value: ", entry%value_r4

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       WRITE(unit,"(A)") "datatype: r8"
       WRITE(unit,"(A,E30.20)") "value: ", entry%value_r8

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       WRITE(unit,"(A)") "datatype: str"
       WRITE(unit,"(A,A)") "value: ", entry%value_str

    ELSE IF(ASSOCIATED(entry%value_1d_i4)) THEN
       WRITE(unit,"(A)") "datatype: 1d_i4"
       WRITE(unit,"(A,I10)") "size: ", SIZE(entry%value_1d_i4)
       DO i=1, SIZE(entry%value_1d_i4)
         WRITE(unit,*) entry%value_1d_i4(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_1d_i8)) THEN
       WRITE(unit,"(A)") "datatype: 1d_i8"
       WRITE(unit,"(A,I20)") "size: ", SIZE(entry%value_1d_i8)
       DO i=1, SIZE(entry%value_1d_i8)
         WRITE(unit,*) entry%value_1d_i8(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_1d_r4)) THEN
       WRITE(unit,"(A)") "datatype: 1d_r4"
       WRITE(unit,"(A,I8)") "size: ", SIZE(entry%value_1d_r4)
       DO i=1, SIZE(entry%value_1d_r4)
         WRITE(unit,"(1X,E30.20)") entry%value_1d_r4(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_1d_r8)) THEN
       WRITE(unit,"(A)") "datatype: 1d_r8"
       WRITE(unit,"(A,I8)") "size: ", SIZE(entry%value_1d_r8)
       DO i=1, SIZE(entry%value_1d_r8)
         WRITE(unit,"(1X,E30.20)") entry%value_1d_r8(i)
       END DO

    ELSE
       CALL cp__b("swarm/swarm_message.F",757,"no value ASSOCIATED")
    END IF
  END SUBROUTINE swarm_message_entry_file_write


! *****************************************************************************
!> \brief Helper routine for swarm_message_file_read.
!> \param ENTRY ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_read(ENTRY, parser, at_end)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'swarm_message_entry_file_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=15)                        :: datatype, label
    INTEGER                                  :: arr_size, i
    LOGICAL                                  :: is_scalar

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:key_length+10),*) label, entry%key
    IF(.NOT.(TRIM(label)=="key:"))CALL cp__a("swarm/swarm_message.F",786)

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:30),*) label, datatype
    IF(.NOT.(TRIM(label)=="datatype:"))CALL cp__a("swarm/swarm_message.F",792)

    CALL parser_get_next_line(parser, 1, at_end)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN

    is_scalar=.TRUE.
    SELECT CASE (TRIM(datatype))
      CASE("i4")
         ALLOCATE(entry%value_i4)
         READ (parser%input_line(1:40),*) label, entry%value_i4
      CASE("i8")
         ALLOCATE(entry%value_i8)
         READ (parser%input_line(1:40),*) label, entry%value_i8
      CASE("r4")
         ALLOCATE(entry%value_r4)
         READ (parser%input_line(1:40),*) label, entry%value_r4
      CASE("r8")
         ALLOCATE(entry%value_r8)
         READ (parser%input_line(1:40),*) label, entry%value_r8
      CASE("str")
         ALLOCATE(entry%value_str)
         READ (parser%input_line(1:40),*) label, entry%value_str
      CASE DEFAULT
         is_scalar = .FALSE.
    END SELECT

    IF(is_scalar) THEN
       IF(.NOT.(TRIM(label)=="value:"))CALL cp__a("swarm/swarm_message.F",820)
       RETURN
    END IF

    ! musst be an array-datatype
    READ (parser%input_line(1:30),*) label, arr_size
    IF(.NOT.(TRIM(label)=="size:"))CALL cp__a("swarm/swarm_message.F",826)

    SELECT CASE (TRIM(datatype))
      CASE("1d_i4")
        ALLOCATE(entry%value_1d_i4(arr_size))
      CASE("1d_i8")
       ALLOCATE(entry%value_1d_i8(arr_size))
      CASE("1d_r4")
        ALLOCATE(entry%value_1d_r4(arr_size))
      CASE("1d_r8")
        ALLOCATE(entry%value_1d_r8(arr_size))
      CASE DEFAULT
       CALL cp__b("swarm/swarm_message.F",838,"unknown datatype")
    END SELECT

    DO i=1, arr_size
      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
      IF(at_end) RETURN

      !Numbers were written with at most 31 characters.
      SELECT CASE (TRIM(datatype))
        CASE("1d_i4")
          READ (parser%input_line(1:31),*) entry%value_1d_i4(i)
        CASE("1d_i8")
          READ (parser%input_line(1:31),*) entry%value_1d_i8(i)
        CASE("1d_r4")
          READ (parser%input_line(1:31),*) entry%value_1d_r4(i)
        CASE("1d_r8")
          READ (parser%input_line(1:31),*) entry%value_1d_r8(i)
        CASE DEFAULT
          CALL cp__b("swarm/swarm_message.F",857,"swarm_message_entry_file_read: unknown datatype")
      END SELECT
    END DO


  END SUBROUTINE swarm_message_entry_file_read


! *****************************************************************************
!> \brief Helper routine, converts a string into an integer-array
!> \param str ...
!> \retval arr ...
!> \author Ole Schuett
! *****************************************************************************
  PURE FUNCTION str2iarr(str) RESULT(arr)
    CHARACTER(LEN=*), INTENT(IN)             :: str
    INTEGER, DIMENSION(LEN(str))             :: arr

    INTEGER                                  :: i

    DO i=1,LEN(str)
      arr(i) = ICHAR(str(i:i))
    ENDDO
  END FUNCTION str2iarr


! *****************************************************************************
!> \brief Helper routine, converts an integer-array into a string
!> \param arr ...
!> \retval str ...
!> \author Ole Schuett
! *****************************************************************************
  PURE FUNCTION iarr2str(arr) RESULT(str)
    INTEGER, DIMENSION(:), INTENT(IN)        :: arr
    CHARACTER(LEN=SIZE(arr))                 :: str

    INTEGER                                  :: i

    DO i=1,SIZE(arr)
      str(i:i) = CHAR(arr(i))
    ENDDO
  END FUNCTION iarr2str



# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_str.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_str(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   CHARACTER(LEN=*), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_str: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_str)

   new_entry%value_str = value

   !WRITE (*,*) "swarm_message_add_str: key=",key, " value=",new_entry%value_str

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_str


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_str(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   CHARACTER(LEN=default_string_length), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_str: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_str)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_str: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_str
         !WRITE (*,*) "swarm_message_get_str: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_str

!EOF
# 902 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_i4.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   INTEGER(KIND=int_4), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_i4: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_i4)

   new_entry%value_i4 = value

   !WRITE (*,*) "swarm_message_add_i4: key=",key, " value=",new_entry%value_i4

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_i4


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   INTEGER(KIND=int_4), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_i4: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_i4)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_i4: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_i4
         !WRITE (*,*) "swarm_message_get_i4: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_i4

!EOF
# 903 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_i8.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   INTEGER(KIND=int_8), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_i8: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_i8)

   new_entry%value_i8 = value

   !WRITE (*,*) "swarm_message_add_i8: key=",key, " value=",new_entry%value_i8

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_i8


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   INTEGER(KIND=int_8), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_i8: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_i8)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_i8: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_i8
         !WRITE (*,*) "swarm_message_get_i8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_i8

!EOF
# 904 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_r4.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_r4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   REAL(KIND=real_4), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_r4: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_r4)

   new_entry%value_r4 = value

   !WRITE (*,*) "swarm_message_add_r4: key=",key, " value=",new_entry%value_r4

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_r4


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_r4(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   REAL(KIND=real_4), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_r4: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_r4)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_r4: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_r4
         !WRITE (*,*) "swarm_message_get_r4: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_r4

!EOF
# 905 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_r8.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   REAL(KIND=real_8), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_r8: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_r8)

   new_entry%value_r8 = value

   !WRITE (*,*) "swarm_message_add_r8: key=",key, " value=",new_entry%value_r8

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_r8


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   REAL(KIND=real_8), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_r8: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_r8)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_r8: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_r8
         !WRITE (*,*) "swarm_message_get_r8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_r8

!EOF
# 906 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_1d_i4.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_1d_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   INTEGER(KIND=int_4), DIMENSION(:) , INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_1d_i4: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_1d_i4(SIZE(value)))

   new_entry%value_1d_i4 = value

   !WRITE (*,*) "swarm_message_add_1d_i4: key=",key, " value=",new_entry%value_1d_i4

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_1d_i4


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_1d_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   INTEGER(KIND=int_4), DIMENSION(:), POINTER                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_1d_i4: key=",key

   IF(ASSOCIATED(value)) CALL cp__b("swarm/swarm_message.F",52,"swarm_message_get_1d_i4: value already associated")

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_1d_i4)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_1d_i4: value not associated key: "//TRIM(key))
         ALLOCATE(value(SIZE(curr_entry%value_1d_i4)))
         value = curr_entry%value_1d_i4
         !WRITE (*,*) "swarm_message_get_1d_i4: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_1d_i4

!EOF
# 907 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_1d_i8.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_1d_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   INTEGER(KIND=int_8), DIMENSION(:) , INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_1d_i8: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_1d_i8(SIZE(value)))

   new_entry%value_1d_i8 = value

   !WRITE (*,*) "swarm_message_add_1d_i8: key=",key, " value=",new_entry%value_1d_i8

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_1d_i8


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_1d_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   INTEGER(KIND=int_8), DIMENSION(:), POINTER                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_1d_i8: key=",key

   IF(ASSOCIATED(value)) CALL cp__b("swarm/swarm_message.F",52,"swarm_message_get_1d_i8: value already associated")

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_1d_i8)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_1d_i8: value not associated key: "//TRIM(key))
         ALLOCATE(value(SIZE(curr_entry%value_1d_i8)))
         value = curr_entry%value_1d_i8
         !WRITE (*,*) "swarm_message_get_1d_i8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_1d_i8

!EOF
# 908 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_1d_r4.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_1d_r4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   REAL(KIND=real_4), DIMENSION(:) , INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_1d_r4: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_1d_r4(SIZE(value)))

   new_entry%value_1d_r4 = value

   !WRITE (*,*) "swarm_message_add_1d_r4: key=",key, " value=",new_entry%value_1d_r4

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_1d_r4


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_1d_r4(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   REAL(KIND=real_4), DIMENSION(:), POINTER                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_1d_r4: key=",key

   IF(ASSOCIATED(value)) CALL cp__b("swarm/swarm_message.F",52,"swarm_message_get_1d_r4: value already associated")

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_1d_r4)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_1d_r4: value not associated key: "//TRIM(key))
         ALLOCATE(value(SIZE(curr_entry%value_1d_r4)))
         value = curr_entry%value_1d_r4
         !WRITE (*,*) "swarm_message_get_1d_r4: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_1d_r4

!EOF
# 909 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2

# 1 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message_1d_r8.f90" 1
! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_1d_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   REAL(KIND=real_8), DIMENSION(:) , INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL cp__b("swarm/swarm_message.F",16,"swarm_message_add_1d_r8: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_1d_r8(SIZE(value)))

   new_entry%value_1d_r8 = value

   !WRITE (*,*) "swarm_message_add_1d_r8: key=",key, " value=",new_entry%value_1d_r8

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_1d_r8


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_1d_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   REAL(KIND=real_8), DIMENSION(:), POINTER                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_1d_r8: key=",key

   IF(ASSOCIATED(value)) CALL cp__b("swarm/swarm_message.F",52,"swarm_message_get_1d_r8: value already associated")

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_1d_r8)) &
            CALL cp__b("swarm/swarm_message.F",58,"swarm_message_get_1d_r8: value not associated key: "//TRIM(key))
         ALLOCATE(value(SIZE(curr_entry%value_1d_r8)))
         value = curr_entry%value_1d_r8
         !WRITE (*,*) "swarm_message_get_1d_r8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL cp__b("swarm/swarm_message.F",66,"swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_1d_r8

!EOF
# 910 "/work/e05/e05/fivanovic/flavoured-cptk-X-SH-coulomb-barrier/cp2k/src/swarm/swarm_message.F" 2


END MODULE swarm_message


