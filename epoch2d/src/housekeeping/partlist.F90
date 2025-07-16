! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE partlist

  USE shared_data
  USE particle_id_hash_mod
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
  USE random_generator
#endif

  IMPLICIT NONE

  SAVE

  INTEGER :: nvar

  TYPE pointer_item
    TYPE(particle), POINTER :: part
    TYPE(pointer_item), POINTER :: next
  END TYPE pointer_item

  TYPE pointer_list
    TYPE(pointer_item), POINTER :: head, tail
  END TYPE pointer_list

  REAL(num), DIMENSION(:), ALLOCATABLE :: packed_particle_data

CONTAINS

  SUBROUTINE set_partlist_size

    nvar = 3 + c_ndims
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    nvar = nvar+1
#endif
#ifdef DELTAF_METHOD
    nvar = nvar+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    nvar = nvar+2
#endif
#ifdef PARTICLE_DEBUG
    nvar = nvar+2
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    nvar = nvar+1
#endif
#ifdef PHOTONS
    nvar = nvar+1
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    nvar = nvar+1
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    nvar = nvar+1
#endif
#ifdef BREMSSTRAHLUNG
    nvar = nvar+1
#endif
#ifdef PROBE_TIME
    nvar = nvar+1
#endif
#ifdef WORK_DONE_INTEGRATED
    nvar = nvar+6
#endif
    ! Persistent IDs
    IF (any_persistent_subset) nvar = nvar+1

  END SUBROUTINE set_partlist_size



  SUBROUTINE setup_partlists

    LOGICAL :: old_any_persistent_subset

    old_any_persistent_subset = any_persistent_subset
    any_persistent_subset = .TRUE.

    CALL set_partlist_size

    any_persistent_subset = old_any_persistent_subset

    ALLOCATE(packed_particle_data(nvar))

  END SUBROUTINE setup_partlists



  SUBROUTINE deallocate_partlists

    INTEGER :: stat

    IF (ALLOCATED(packed_particle_data)) &
        DEALLOCATE(packed_particle_data, STAT=stat)

  END SUBROUTINE deallocate_partlists



  SUBROUTINE create_empty_partlist(partlist, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies

    NULLIFY(partlist%head)
    NULLIFY(partlist%tail)
    partlist%count = 0
    partlist%id_update = 0
    partlist%safe = .TRUE.
    IF (PRESENT(holds_copies)) THEN
      partlist%holds_copies = holds_copies
    ELSE
      partlist%holds_copies = .FALSE.
    END IF

  END SUBROUTINE create_empty_partlist



  SUBROUTINE create_unsafe_partlist(partlist, a_particle, n_elements, &
      holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist, holds_copies)

    partlist%safe = .FALSE.
    current => a_particle
    ipart = 1
    DO WHILE (ASSOCIATED(current) .AND. ipart < n_elements)
      ipart = ipart+1
      current => current%next
    END DO
    partlist%head => a_particle
    partlist%tail => current
    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist



  SUBROUTINE create_unsafe_partlist_by_tail(partlist, head, tail, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: head, tail
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist, holds_copies)

    partlist%safe = .FALSE.
    partlist%head => head
    partlist%tail => tail

    current => head
    ipart = 0
    DO WHILE (ASSOCIATED(current))
      ipart = ipart+1
      current => current%next
      IF (ASSOCIATED(current)) THEN
        IF (ASSOCIATED(current%prev, TARGET=tail)) EXIT
      END IF
    END DO

    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist_by_tail



  SUBROUTINE create_allocated_partlist(partlist, n_elements, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: new_particle
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist, holds_copies)

    DO ipart = 0, n_elements-1
      CALL create_particle(new_particle)
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    END DO

  END SUBROUTINE create_allocated_partlist



  SUBROUTINE create_filled_partlist(partlist, data_in, n_elements, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: data_in
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: new_particle
    INTEGER(i8) :: ipart, cpos = 0

    CALL set_partlist_size
    CALL create_empty_partlist(partlist, holds_copies)

    DO ipart = 0, n_elements-1
      ALLOCATE(new_particle)
      cpos = ipart*nvar+1
      CALL unpack_particle(data_in(cpos:cpos+nvar-1), new_particle)
#ifdef PARTICLE_DEBUG
      new_particle%processor = rank
#endif
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    END DO

  END SUBROUTINE create_filled_partlist



  FUNCTION test_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: current
    INTEGER :: test_partlist
    INTEGER(i8) :: test_ct

    test_partlist = 0
    test_ct = 0

    ! Empty list is OK
    IF (.NOT. ASSOCIATED(partlist%head) &
        .AND. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = 0
      RETURN
    END IF

    ! List with head or tail but not both is broken
    IF (.NOT. ASSOCIATED(partlist%head) &
        .OR. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = -1
      RETURN
    END IF

    ! Having head and tail elements which are not the end of a list are OK for
    ! unsafe partlists
    IF (ASSOCIATED(partlist%head%prev) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 1)
    IF (ASSOCIATED(partlist%tail%next) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 2)

    ! Since we don't KNOW that count is OK (that's what we're checking)
    ! Have to check both for end of list and for having reached the tail item
    current => partlist%head
    DO WHILE (ASSOCIATED(current))
      test_ct = test_ct+1
      current => current%next
      IF (ASSOCIATED(current)) THEN
        ! This tests if we've just jumped to the tail element
        ! Allows testing of unsafe partlists
        IF (ASSOCIATED(current%prev, TARGET=partlist%tail)) EXIT
      END IF
    END DO

    IF (test_ct /= partlist%count) test_partlist = IOR(test_partlist, 4)

  END FUNCTION test_partlist



  SUBROUTINE destroy_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle, next
    INTEGER(i8) :: ipart

    ! Go through list and delete all the particles in the list
    new_particle => partlist%head
    ipart = 0
    DO WHILE (ipart < partlist%count)
      next => new_particle%next
      ! A partlist that holds copies or an unsafe partlist should not cause
      ! unlinking
      CALL destroy_particle(new_particle, &
          partlist%holds_copies .OR. .NOT.partlist%safe)
      new_particle => next
      ipart = ipart+1
    END DO

    CALL create_empty_partlist(partlist)

  END SUBROUTINE destroy_partlist



  SUBROUTINE copy_partlist(partlist1, partlist2)

    TYPE(particle_list), INTENT(INOUT) :: partlist1, partlist2

    partlist2%head => partlist1%head
    partlist2%tail => partlist1%tail
    partlist2%count = partlist1%count
    partlist2%id_update = partlist1%id_update
    partlist2%holds_copies = partlist1%holds_copies

  END SUBROUTINE copy_partlist



  SUBROUTINE append_partlist(head, tail)

    TYPE(particle_list), INTENT(INOUT) :: head, tail

    IF (.NOT. head%safe .OR. .NOT. tail%safe) THEN
      IF (rank == 0) &
          PRINT *, 'Unable to append partlists because one is not safe'
      RETURN
    END IF

    IF (ASSOCIATED(head%tail)) THEN
      head%tail%next => tail%head
    ELSE
      head%head => tail%head
    END IF
    IF (ASSOCIATED(tail%head)) tail%head%prev => head%tail
    IF (ASSOCIATED(tail%tail)) head%tail => tail%tail
    head%count = head%count + tail%count
    head%id_update = head%id_update + tail%id_update

    CALL create_empty_partlist(tail)

  END SUBROUTINE append_partlist



  SUBROUTINE add_particle_to_partlist(partlist, new_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    ! if (!particle) return;
    IF (.NOT. ASSOCIATED(new_particle)) RETURN
    NULLIFY(new_particle%next, new_particle%prev)

    ! Add particle count
    partlist%count = partlist%count + 1
    partlist%id_update = 1
    IF (.NOT. ASSOCIATED(partlist%tail)) THEN
      ! partlist is empty
      partlist%head => new_particle
      partlist%tail => new_particle
      RETURN
    END IF

    partlist%tail%next => new_particle
    new_particle%prev => partlist%tail
    NULLIFY(new_particle%next)
    partlist%tail => new_particle

  END SUBROUTINE add_particle_to_partlist



  SUBROUTINE remove_particle_from_partlist(partlist, a_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    ! Check whether particle is head or tail of list and unlink
    IF (ASSOCIATED(partlist%head, TARGET=a_particle)) &
        partlist%head => a_particle%next
    IF (ASSOCIATED(partlist%tail, TARGET=a_particle)) &
        partlist%tail => a_particle%prev

    ! Link particles on either side together
    IF (ASSOCIATED(a_particle%next)) a_particle%next%prev => a_particle%prev
    IF (ASSOCIATED(a_particle%prev)) a_particle%prev%next => a_particle%next

    NULLIFY(a_particle%next, a_particle%prev)

    ! Decrement counter
    partlist%count = partlist%count-1

  END SUBROUTINE remove_particle_from_partlist



  SUBROUTINE pack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos, temp_i8

    cpos = 1
    array(cpos:cpos+c_ndims-1) = a_particle%part_pos
    cpos = cpos+c_ndims
    array(cpos:cpos+2) = a_particle%part_p
    cpos = cpos+3
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    array(cpos) = a_particle%weight
    cpos = cpos+1
#endif
#ifdef DELTAF_METHOD
    array(cpos) = a_particle%pvol
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    array(cpos) = a_particle%charge
    array(cpos+1) = a_particle%mass
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    array(cpos) = REAL(a_particle%processor, num)
    array(cpos+1) = REAL(a_particle%processor_at_t0, num)
    cpos = cpos+2
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    array(cpos) = REAL(a_particle%id, num)
    cpos = cpos+1
#endif
#ifdef PHOTONS
    array(cpos) = a_particle%optical_depth
    cpos = cpos+1
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    array(cpos) = a_particle%particle_energy
    cpos = cpos+1
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    array(cpos) = a_particle%optical_depth_tri
    cpos = cpos+1
#endif
#ifdef BREMSSTRAHLUNG
    array(cpos) = a_particle%optical_depth_bremsstrahlung
    cpos = cpos+1
#endif
#ifdef PROBE_TIME
    array(cpos) = a_particle%probe_time
    cpos = cpos+1
#endif
#ifdef WORK_DONE_INTEGRATED
    array(cpos) = a_particle%work_x
    array(cpos+1) = a_particle%work_y
    array(cpos+2) = a_particle%work_z
    array(cpos+3) = a_particle%work_x_total
    array(cpos+4) = a_particle%work_y_total
    array(cpos+5) = a_particle%work_z_total
    cpos = cpos+6
#endif
    IF (any_persistent_subset) THEN
      temp_i8 = id_registry%map(a_particle)
      array(cpos) = TRANSFER(temp_i8, 1.0_num)
      cpos = cpos+1
    END IF

  END SUBROUTINE pack_particle



  SUBROUTINE unpack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(IN) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos, temp_i8

    cpos = 1
    a_particle%part_pos = array(cpos:cpos+c_ndims-1)
    cpos = cpos+c_ndims
    a_particle%part_p = array(cpos:cpos+2)
    cpos = cpos+3
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    a_particle%weight = array(cpos)
    cpos = cpos+1
#endif
#ifdef DELTAF_METHOD
    a_particle%pvol = array(cpos)
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    a_particle%charge = array(cpos)
    a_particle%mass = array(cpos+1)
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    a_particle%processor = rank
    a_particle%processor_at_t0 = NINT(array(cpos+1))
    cpos = cpos+2
#endif
#ifdef PARTICLE_ID4
    a_particle%id = NINT(array(cpos))
    cpos = cpos+1
#elif PARTICLE_ID
    a_particle%id = NINT(array(cpos),i8)
    cpos = cpos+1
#endif
#ifdef PHOTONS
    a_particle%optical_depth = array(cpos)
    cpos = cpos+1
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    a_particle%particle_energy = array(cpos)
    cpos = cpos+1
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    a_particle%optical_depth_tri = array(cpos)
    cpos = cpos+1
#endif
#ifdef BREMSSTRAHLUNG
    a_particle%optical_depth_bremsstrahlung = array(cpos)
    cpos = cpos+1
#endif
#ifdef PROBE_TIME
    a_particle%probe_time = array(cpos)
    cpos = cpos+1
#endif
#ifdef WORK_DONE_INTEGRATED
    a_particle%work_x = array(cpos)
    a_particle%work_y = array(cpos+1)
    a_particle%work_z = array(cpos+2)
    a_particle%work_x_total = array(cpos+3)
    a_particle%work_y_total = array(cpos+4)
    a_particle%work_z_total = array(cpos+5)
    cpos = cpos+6
#endif
    IF (any_persistent_subset) THEN
      CALL id_registry%add_with_map(a_particle, TRANSFER(array(cpos), temp_i8))
      cpos = cpos+1
    END IF

  END SUBROUTINE unpack_particle



  SUBROUTINE init_particle(new_particle)

    TYPE(particle), POINTER :: new_particle

    new_particle%part_p = 0.0_num
    new_particle%part_pos = 0.0_num
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    new_particle%weight = 0.0_num
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    new_particle%charge = 0.0_num
    new_particle%mass = 0.0_num
#endif
#ifdef PARTICLE_DEBUG
    new_particle%processor = 0
    new_particle%processor_at_t0 = 0
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    new_particle%id = 0
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    ! This assigns an optical depth to newly created particle
    new_particle%particle_energy = 0.0_num
#endif
#ifdef PHOTONS
    new_particle%optical_depth = LOG(1.0_num / (1.0_num - random()))
#ifdef TRIDENT_PHOTONS
    new_particle%optical_depth_tri = LOG(1.0_num / (1.0_num - random()))
#endif
#endif
#ifdef BREMSSTRAHLUNG
    new_particle%optical_depth_bremsstrahlung = &
        LOG(1.0_num / (1.0_num - random()))
#endif
#ifdef PROBE_TIME
    new_particle%probe_time = 0.0_num
#endif
#ifdef WORK_DONE_INTEGRATED
    new_particle%work_x = 0_num
    new_particle%work_y = 0_num
    new_particle%work_z = 0_num
    new_particle%work_x_total = 0_num
    new_particle%work_y_total = 0_num
    new_particle%work_z_total = 0_num
#endif 

  END SUBROUTINE init_particle



  SUBROUTINE create_particle(new_particle)

    TYPE(particle), POINTER :: new_particle

    ALLOCATE(new_particle)
    CALL init_particle(new_particle)

  END SUBROUTINE create_particle



  SUBROUTINE destroy_particle(part, is_copy)

    ! Routine to delete a particle. This routine is only safe to use on
    ! a particle that is not in a partlist
    TYPE(particle), POINTER :: part
    LOGICAL, INTENT(IN), OPTIONAL :: is_copy

    IF (any_persistent_subset) THEN
      IF (PRESENT(is_copy)) THEN
        IF (.NOT. is_copy) CALL id_registry%delete_all(part)
      ELSE
        CALL id_registry%delete_all(part)
      END IF
    END IF

    DEALLOCATE(part)

  END SUBROUTINE destroy_particle



  SUBROUTINE display_particle(a_particle)

    TYPE(particle), POINTER :: a_particle

    PRINT *, 'Position', a_particle%part_pos
    PRINT *, 'Momentum', a_particle%part_p

  END SUBROUTINE display_particle



  FUNCTION compare_particles(part1, part2)

    TYPE(particle), POINTER :: part1, part2
    LOGICAL :: compare_particles

    compare_particles = .TRUE.
    IF (MAXVAL(ABS(part1%part_pos-part2%part_pos)) > c_tiny) &
        compare_particles = .FALSE.
    IF (MAXVAL(ABS(part1%part_p - part2%part_p)) > c_tiny) &
        compare_particles = .FALSE.

#ifndef PER_SPECIES_WEIGHT
    IF (ABS(part1%weight - part2%weight) > c_tiny) &
        compare_particles = .FALSE.
#endif

#ifdef PER_PARTICLE_CHARGE_MASS
    IF (ABS(part1%charge - part2%charge) > c_tiny) &
        compare_particles = .FALSE.
    IF (ABS(part1%mass - part2%mass) > c_tiny) &
        compare_particles = .FALSE.
#endif

    IF (.NOT. compare_particles) THEN
      CALL display_particle(part1)
      CALL display_particle(part2)
    END IF

  END FUNCTION compare_particles



  FUNCTION test_packed_particles(partlist, array, npart_in_data)

    TYPE(particle_list), INTENT(IN) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npart_in_data
    TYPE(particle), POINTER :: current
    TYPE(particle), POINTER :: a_particle
    LOGICAL :: test_packed_particles
    INTEGER(i8) :: ipart

    CALL set_partlist_size

    test_packed_particles = .FALSE.

    IF (npart_in_data * nvar /= SIZE(array)) THEN
      PRINT *, 'Size of data array does not match specified on', rank, &
          npart_in_data, SIZE(array)
      RETURN
    END IF
    IF (partlist%count /= npart_in_data) THEN
      PRINT *, 'Size of data array does not match partlist on', rank
      RETURN
    END IF

    ALLOCATE(a_particle)

    current => partlist%head
    DO ipart = 0, npart_in_data-1
      CALL unpack_particle(array(ipart*nvar+1:(ipart+1)*nvar), a_particle)
      IF (.NOT. compare_particles(a_particle, current)) THEN
        PRINT *, 'BAD PARTICLE ', ipart, 'on', rank
        RETURN
      END IF
      current => current%next
    END DO

    DEALLOCATE(a_particle) !DO NOT REPLACE WITH CALL TO destroy_particle

    test_packed_particles = .TRUE.

  END FUNCTION test_packed_particles



  SUBROUTINE partlist_send_nocount(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER :: ipart, nsend, cpos
    TYPE(particle), POINTER :: current

    CALL set_partlist_size

    nsend = INT(partlist%count) * nvar
    ALLOCATE(array(nsend))
    array = 0.0_num

    current => partlist%head
    ipart = 0
    cpos = 0
    DO WHILE (ipart < partlist%count)
      cpos = ipart * nvar + 1
      CALL pack_particle(array(cpos:cpos+nvar-1), current)
      ipart = ipart + 1
      current => current%next
    END DO

    CALL MPI_SEND(array, nsend, mpireal, dest, tag, comm, errcode)

    DEALLOCATE(array)

  END SUBROUTINE partlist_send_nocount



  SUBROUTINE partlist_send(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    INTEGER(i8) :: send_buf(2)

    send_buf(1) = partlist%count
    send_buf(2) = partlist%id_update

    CALL MPI_SEND(send_buf, 2, MPI_INTEGER8, dest, tag, comm, errcode)

    CALL partlist_send_nocount(partlist, dest)

  END SUBROUTINE partlist_send



  SUBROUTINE partlist_recv_nocount(partlist, src, count)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: src
    INTEGER(i8), INTENT(IN) :: count
    INTEGER :: nrecv
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    CALL set_partlist_size
    CALL create_empty_partlist(partlist)

    nrecv = INT(count) * nvar
    ALLOCATE(array(nrecv))
    array = 0.0_num

    CALL MPI_RECV(array, nrecv, mpireal, src, tag, comm, status, errcode)
    CALL create_filled_partlist(partlist, array, count)

    DEALLOCATE(array)

  END SUBROUTINE partlist_recv_nocount



  SUBROUTINE partlist_recv(partlist, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: src
    INTEGER(i8) :: count, recv_buf(2)

    recv_buf = 0
    CALL MPI_RECV(recv_buf, 2, MPI_INTEGER8, src, tag, comm, status, errcode)
    count = recv_buf(1)
    partlist%id_update = partlist%id_update + INT(recv_buf(2))

    CALL partlist_recv_nocount(partlist, src, count)

  END SUBROUTINE partlist_recv



  SUBROUTINE partlist_sendrecv(partlist_send, partlist_recv, dest, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist_send, partlist_recv
    INTEGER, INTENT(IN) :: dest, src
    REAL(num), DIMENSION(:), ALLOCATABLE :: data_send, data_recv
    INTEGER(i8) :: cpos = 0, ipart = 0
    INTEGER(i8) :: npart_recv, send_buf(2), recv_buf(2)
    INTEGER :: nsend, nrecv
    TYPE(particle), POINTER :: current

    ! This subroutine doesn't try to use memory efficient buffering, it sends
    ! all the particles at once. This should work for boundary calls, but
    ! don't try it for any other reason

    CALL set_partlist_size

    recv_buf = 0
    send_buf(1) = partlist_send%count
    send_buf(2) = partlist_send%id_update

    CALL MPI_SENDRECV(send_buf, 2, MPI_INTEGER8, dest, tag, recv_buf, 2, &
        MPI_INTEGER8, src, tag, comm, status, errcode)

    npart_recv = recv_buf(1)
    nsend = INT(send_buf(1)) * nvar
    nrecv = INT(npart_recv) * nvar
    partlist_recv%id_update = partlist_recv%id_update + INT(recv_buf(2))

    ! Copy the data for the particles into a buffer
    ALLOCATE(data_send(nsend))
    ALLOCATE(data_recv(nrecv))

    ! Pack particles to send into buffer
    current => partlist_send%head
    ipart = 0
    DO WHILE (ipart < partlist_send%count)
      cpos = ipart * nvar + 1
      CALL pack_particle(packed_particle_data, current)
      data_send(cpos:cpos+nvar-1) = packed_particle_data(1:nvar)
      ipart = ipart + 1
      current => current%next
    END DO

    ! No longer need the sending partlist, so destroy it to save some memory
    CALL destroy_partlist(partlist_send)

    ! Actual MPI commands
    CALL MPI_SENDRECV(data_send, nsend, mpireal, dest, tag, &
        data_recv, nrecv, mpireal, src, tag, comm, status, errcode)

    DEALLOCATE(data_send)
    CALL create_filled_partlist(partlist_recv, data_recv, npart_recv)
    DEALLOCATE(data_recv)

  END SUBROUTINE partlist_sendrecv



  SUBROUTINE add_particle_to_list(part, list)

    TYPE(particle), POINTER :: part
    TYPE(pointer_list) :: list
    TYPE(pointer_item), POINTER :: item

    ALLOCATE(item)
    item%part => part
    NULLIFY(item%next)

    list%tail%next => item
    list%tail => item

  END SUBROUTINE add_particle_to_list



  SUBROUTINE generate_particle_ids(partlist)

    TYPE(particle_list) :: partlist
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    INTEGER(i8), ALLOCATABLE :: nid_all(:)
    INTEGER(i8) :: nid, part_id
    INTEGER :: i, id_update
    TYPE(particle), POINTER :: current
    TYPE(pointer_list) :: idlist
    TYPE(pointer_item), POINTER :: idcurrent, idnext

    id_update = partlist%id_update

    CALL MPI_ALLREDUCE(id_update, partlist%id_update, 1, MPI_INTEGER, &
        MPI_MAX, comm, errcode)

    IF (partlist%id_update == 0) RETURN

    ALLOCATE(idlist%head)
    idlist%tail => idlist%head
    NULLIFY(idlist%head%next)
    NULLIFY(idlist%head%part)

    ! Scan through particle list and identify particles which need
    ! an ID to be assigned.
    nid = 0
    current => partlist%head
    DO WHILE(ASSOCIATED(current))
      IF (current%id == 0) THEN
        nid = nid + 1
        CALL add_particle_to_list(current, idlist)
      END IF
      current => current%next
    END DO

    ALLOCATE(nid_all(nproc))

    CALL MPI_ALLGATHER(nid, 1, MPI_INTEGER8, nid_all, 1, MPI_INTEGER8, &
        comm, errcode)

    ! Count number of particles on ranks zero to rank-1
    nid = 0
    DO i = 1, rank
      nid = nid + nid_all(i)
    END DO
    part_id = particles_max_id + nid

    ! Count remaining particles
    DO i = rank+1, nproc
      nid = nid + nid_all(i)
    END DO

    particles_max_id = particles_max_id + nid

    DEALLOCATE(nid_all)

    ! Number each particle with a unique id
    idcurrent => idlist%head%next
    DO WHILE(ASSOCIATED(idcurrent))
      part_id = part_id + 1
#if PARTICLE_ID
      idcurrent%part%id = part_id
#else
      idcurrent%part%id = INT(part_id,i4)
#endif
      idnext => idcurrent%next
      DEALLOCATE(idcurrent)
      idcurrent => idnext
    END DO

    DEALLOCATE(idlist%head)

    partlist%id_update = 0
#endif

  END SUBROUTINE generate_particle_ids



  SUBROUTINE update_particle_count

    ! This routine ensures that the particle count for the species_list
    ! objects is accurate. This makes some things easier, but increases
    ! communication
    INTEGER :: ispecies
    LOGICAL, SAVE :: update = .TRUE.

    IF (.NOT.update) RETURN

    DO ispecies = 1, n_species
      CALL MPI_ALLREDUCE(species_list(ispecies)%attached_list%count, &
          species_list(ispecies)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      species_list(ispecies)%count_update_step = step
    END DO

    update = use_particle_count_update

  END SUBROUTINE update_particle_count

END MODULE partlist
