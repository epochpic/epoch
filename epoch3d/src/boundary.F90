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

MODULE boundary

  USE partlist
  USE particle_temperature
  USE laser
  USE mpi_subtype_control
  USE utilities
  USE particle_id_hash_mod
  USE injectors

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_boundaries

    INTEGER :: i, ispecies, bc
    LOGICAL :: error
    CHARACTER(LEN=5), DIMENSION(2*c_ndims) :: &
        boundary = (/ 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max' /)
    CHARACTER(LEN=2*c_max_string_length) :: bc_error

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    add_laser(:) = .FALSE.
    any_open = .FALSE.
    cpml_boundaries = .FALSE.
    DO i = 1, 2*c_ndims
      IF (bc_field(i) == c_bc_other) bc_field(i) = c_bc_clamp
      IF (bc_field(i) == c_bc_cpml_laser &
          .OR. bc_field(i) == c_bc_cpml_outflow) cpml_boundaries = .TRUE.
      IF (bc_field(i) == c_bc_simple_laser) THEN
        add_laser(i) = .TRUE.
        any_open = .TRUE.
      END IF

      ! Note: reflecting EM boundaries not yet implemented.
      IF (bc_field(i) == c_bc_reflect) bc_field(i) = c_bc_clamp
      IF (bc_field(i) == c_bc_open) bc_field(i) = c_bc_simple_outflow
      IF (bc_field(i) == c_bc_simple_outflow) any_open = .TRUE.
    END DO

    error = .FALSE.
    DO ispecies = 1, n_species
      DO i = 1, 2*c_ndims
        bc = species_list(ispecies)%bc_particle(i)
        bc_error = 'Unrecognised "' // TRIM(boundary(i)) // '" boundary for ' &
            // 'species "' // TRIM(species_list(ispecies)%name) // '"'
        error = error .OR. setup_particle_boundary(bc, bc_error)
      END DO
    END DO

    IF (error) THEN
      errcode = c_err_bad_value
      CALL abort_code(errcode)
    END IF

  END SUBROUTINE setup_boundaries



  SUBROUTINE setup_domain_dependent_boundaries

    ! Any boundary condition setup that needs the domain to have already been
    ! created should be added here

    INTEGER :: ispecies, i, bc

    DO ispecies = 1, n_species
      DO i = 1, 2*c_ndims
        bc = species_list(ispecies)%bc_particle(i)
        IF (bc == c_bc_heat_bath) THEN
          CALL create_boundary_injector(ispecies, i)
          species_list(ispecies)%bc_particle(i) = c_bc_open
        END IF
      END DO
    END DO

  END SUBROUTINE setup_domain_dependent_boundaries



  FUNCTION setup_particle_boundary(boundary, bc_error) RESULT(error)

    INTEGER, INTENT(INOUT) :: boundary
    CHARACTER(LEN=*), INTENT(IN) :: bc_error
    LOGICAL :: error

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    IF (boundary == c_bc_other .OR. boundary == c_bc_conduct) &
        boundary = c_bc_reflect

    ! Note, for laser bcs to work, the main bcs must be set IN THE CODE to
    ! simple_laser (or outflow) and the field bcs to c_bc_clamp. Particles
    ! can then be set separately. IN THE DECK, laser bcs are chosen either
    ! by seting the main bcs OR by setting the field bcs to simple_laser
    ! (or outflow).

    ! Laser boundaries assume open particles unless otherwise specified.
    IF (boundary == c_bc_simple_laser &
        .OR. boundary == c_bc_simple_outflow &
        .OR. boundary == c_bc_cpml_laser &
        .OR. boundary == c_bc_cpml_outflow) &
            boundary = c_bc_open

    ! Sanity check on particle boundaries
    error = .FALSE.
    IF (boundary == c_bc_periodic &
        .OR. boundary == c_bc_reflect &
        .OR. boundary == c_bc_thermal &
        .OR. boundary == c_bc_heat_bath &
        .OR. boundary == c_bc_open) RETURN

    IF (rank == 0) THEN
      WRITE(*,*)
      WRITE(*,*) '*** ERROR ***'
      WRITE(*,*) TRIM(bc_error)
    END IF
    error = .TRUE.

  END FUNCTION setup_particle_boundary



  ! Exchanges field values at processor boundaries and applies field
  ! boundary conditions
  SUBROUTINE field_bc(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, ng, nx, ny, nz)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths_slice(field, direction, ng, n1_local, &
      n2_local)

    INTEGER, INTENT(IN) :: direction, ng
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: n1_local, n2_local
    INTEGER :: proc1_min, proc1_max
    INTEGER :: proc2_min, proc2_max
    INTEGER, DIMENSION(c_ndims-1) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    IF (direction == c_dir_x) THEN
      IF (.NOT. y_min_boundary .OR. bc_field(c_bd_y_min) == c_bc_periodic) THEN
        proc1_min = proc_y_min
      ELSE
        proc1_min = MPI_PROC_NULL
      END IF
      IF (.NOT. y_max_boundary .OR. bc_field(c_bd_y_max) == c_bc_periodic) THEN
        proc1_max = proc_y_max
      ELSE
        proc1_max = MPI_PROC_NULL
      END IF
      IF (.NOT. z_min_boundary .OR. bc_field(c_bd_z_min) == c_bc_periodic) THEN
        proc2_min = proc_z_min
      ELSE
        proc2_min = MPI_PROC_NULL
      END IF
      IF (.NOT. z_max_boundary .OR. bc_field(c_bd_z_max) == c_bc_periodic) THEN
        proc2_max = proc_z_max
      ELSE
        proc2_max = MPI_PROC_NULL
      END IF
    ELSE IF (direction == c_dir_y) THEN
      IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
        proc1_min = proc_x_min
      ELSE
        proc1_min = MPI_PROC_NULL
      END IF
      IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
        proc1_max = proc_x_max
      ELSE
        proc1_max = MPI_PROC_NULL
      END IF
      IF (.NOT. z_min_boundary .OR. bc_field(c_bd_z_min) == c_bc_periodic) THEN
        proc2_min = proc_z_min
      ELSE
        proc2_min = MPI_PROC_NULL
      END IF
      IF (.NOT. z_max_boundary .OR. bc_field(c_bd_z_max) == c_bc_periodic) THEN
        proc2_max = proc_z_max
      ELSE
        proc2_max = MPI_PROC_NULL
      END IF
    ELSE
      IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
        proc1_min = proc_x_min
      ELSE
        proc1_min = MPI_PROC_NULL
      END IF
      IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
        proc1_max = proc_x_max
      ELSE
        proc1_max = MPI_PROC_NULL
      END IF
      IF (.NOT. y_min_boundary .OR. bc_field(c_bd_y_min) == c_bc_periodic) THEN
        proc2_min = proc_y_min
      ELSE
        proc2_min = MPI_PROC_NULL
      END IF
      IF (.NOT. y_max_boundary .OR. bc_field(c_bd_y_max) == c_bc_periodic) THEN
        proc2_max = proc_y_max
      ELSE
        proc2_max = MPI_PROC_NULL
      END IF
    END IF

    sizes(1) = n1_local + 2 * ng
    sizes(2) = n2_local + 2 * ng
    starts = 1

    szmax = sizes(1) * ng
    sz = sizes(2) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)

    sz = subsizes(1) * subsizes(2)

    subarray = create_2d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng), 1, subarray, proc1_min, &
        tag, temp, sz, basetype, proc1_max, tag, comm, status, errcode)

    IF (proc1_max /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = n1_local+1, subsizes(1)+n1_local
        field(i,j) = temp(n)
        n = n + 1
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(n1_local+1-ng,1-ng), 1, subarray, proc1_max, &
        tag, temp, sz, basetype, proc1_min, tag, comm, status, errcode)

    IF (proc1_min /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng

    sz = subsizes(1) * subsizes(2)

    subarray = create_2d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1), 1, subarray, proc2_min, &
        tag, temp, sz, basetype, proc2_max, tag, comm, status, errcode)

    IF (proc2_max /= MPI_PROC_NULL) THEN
      n = 1
      DO j = n2_local+1, subsizes(2)+n2_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(1-ng,n2_local+1-ng), 1, subarray, proc2_max, &
        tag, temp, sz, basetype, proc2_min, tag, comm, status, errcode)

    IF (proc2_min /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths_slice



  SUBROUTINE do_field_mpi_with_lengths(field, ng, nx_local, ny_local, &
      nz_local)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    sizes(1) = nx_local + 2 * ng
    sizes(2) = ny_local + 2 * ng
    sizes(3) = nz_local + 2 * ng
    starts = 1

    szmax = sizes(1) * sizes(2) * ng
    sz = sizes(1) * sizes(3) * ng
    IF (sz > szmax) szmax = sz
    sz = sizes(2) * sizes(3) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng,1-ng), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(nx_local+1-ng,1-ng,1-ng), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1,1-ng), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (.NOT. y_max_boundary .OR. bc_field(c_bd_y_max) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(1-ng,ny_local+1-ng,1-ng), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (.NOT. y_min_boundary .OR. bc_field(c_bd_y_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1-ng,1), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (.NOT. z_max_boundary .OR. bc_field(c_bd_z_max) == c_bc_periodic) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(1-ng,1-ng,nz_local+1-ng), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (.NOT. z_min_boundary .OR. bc_field(c_bd_z_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE do_field_mpi_with_lengths_r4(field, ng, nx_local, ny_local, &
      nz_local)

    INTEGER, INTENT(IN) :: ng
    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(r4), ALLOCATABLE :: temp(:)

    basetype = MPI_REAL4

    sizes(1) = nx_local + 2 * ng
    sizes(2) = ny_local + 2 * ng
    sizes(3) = nz_local + 2 * ng
    starts = 1

    szmax = sizes(1) * sizes(2) * ng
    sz = sizes(1) * sizes(3) * ng
    IF (sz > szmax) szmax = sz
    sz = sizes(2) * sizes(3) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng,1-ng), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(nx_local+1-ng,1-ng,1-ng), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1,1-ng), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (.NOT. y_max_boundary .OR. bc_field(c_bd_y_max) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(1-ng,ny_local+1-ng,1-ng), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (.NOT. y_min_boundary .OR. bc_field(c_bd_y_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1-ng,1), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (.NOT. z_max_boundary .OR. bc_field(c_bd_z_max) == c_bc_periodic) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_SENDRECV(field(1-ng,1-ng,nz_local+1-ng), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (.NOT. z_min_boundary .OR. bc_field(c_bd_z_min) == c_bc_periodic) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths_r4



  SUBROUTINE field_zero_gradient(field, stagger_type, boundary)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(i-ng,:,:) = field(ng-i,:,:)
        END DO
      ELSE
        DO i = 1, ng
          field(i-ng,:,:) = field(ng+1-i,:,:)
        END DO
      END IF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(nn+i,:,:) = field(nn-i,:,:)
        END DO
      ELSE
        DO i = 1, ng
          field(nn+i,:,:) = field(nn+1-i,:,:)
        END DO
      END IF

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,i-ng,:) = field(:,ng-i,:)
        END DO
      ELSE
        DO i = 1, ng
          field(:,i-ng,:) = field(:,ng+1-i,:)
        END DO
      END IF
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      nn = ny
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,nn+i,:) = field(:,nn-i,:)
        END DO
      ELSE
        DO i = 1, ng
          field(:,nn+i,:) = field(:,nn+1-i,:)
        END DO
      END IF

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,:,i-ng) = field(:,:,ng-i)
        END DO
      ELSE
        DO i = 1, ng
          field(:,:,i-ng) = field(:,:,ng+1-i)
        END DO
      END IF
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      nn = nz
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,:,nn+i) = field(:,:,nn-i)
        END DO
      ELSE
        DO i = 1, ng
          field(:,:,nn+i) = field(:,:,nn+1-i)
        END DO
      END IF
    END IF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, ng, stagger_type, boundary)

    INTEGER, INTENT(IN) :: ng, stagger_type, boundary
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(i-ng,:,:) = -field(ng-i,:,:)
        END DO
        field(0,:,:) = 0.0_num
      ELSE
        DO i = 1, ng
          field(i-ng,:,:) = -field(ng+1-i,:,:)
        END DO
      END IF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nn,:,:) = 0.0_num
        DO i = 1, ng-1
          field(nn+i,:,:) = -field(nn-i,:,:)
        END DO
      ELSE
        DO i = 1, ng
          field(nn+i,:,:) = -field(nn+1-i,:,:)
        END DO
      END IF

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,i-ng,:) = -field(:,ng-i,:)
        END DO
        field(:,0,:) = 0.0_num
      ELSE
        DO i = 1, ng
          field(:,i-ng,:) = -field(:,ng+1-i,:)
        END DO
      END IF
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      nn = ny
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,nn,:) = 0.0_num
        DO i = 1, ng-1
          field(:,nn+i,:) = -field(:,nn-i,:)
        END DO
      ELSE
        DO i = 1, ng
          field(:,nn+i,:) = -field(:,nn+1-i,:)
        END DO
      END IF

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,:,i-ng) = -field(:,:,ng-i)
        END DO
        field(:,:,0) = 0.0_num
      ELSE
        DO i = 1, ng
          field(:,:,i-ng) = -field(:,:,ng+1-i)
        END DO
      END IF
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      nn = nz
      IF (stagger(c_dir_z,stagger_type)) THEN
        field(:,:,nn) = 0.0_num
        DO i = 1, ng-1
          field(:,:,nn+i) = -field(:,:,nn-i)
        END DO
      ELSE
        DO i = 1, ng
          field(:,:,nn+i) = -field(:,:,nn+1-i)
        END DO
      END IF
    END IF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE particle_reflection_bcs(array, ng, flip_direction, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: flip_direction
    INTEGER, INTENT(IN), OPTIONAL :: species
    INTEGER, DIMENSION(c_ndims) :: sizes
    INTEGER :: nn, n, i, flip_dir, bc
    INTEGER, DIMENSION(2*c_ndims) :: bc_species

    flip_dir = 0
    IF (PRESENT(flip_direction)) flip_dir = flip_direction

    bc_species = bc_allspecies
    IF (PRESENT(species)) THEN
      DO i = 1, 2*c_ndims
        IF (bc_species(i) == c_bc_mixed) THEN
          bc_species(i) = species_list(species)%bc_particle(i)
        END IF
      END DO
    END IF

    sizes = SHAPE(array)
    n = 0
    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    bc = bc_species(n)
    IF (x_min_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(i,:,:) = array(i,:,:) - array(-i,:,:)
          array(-i,:,:) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng-1
          array(i,:,:) = array(i,:,:) + array(1-i,:,:)
          array(1-i,:,:) = 0.0_num
        END DO
      END IF
    END IF

    n = n + 1
    bc = bc_species(n)
    IF (x_max_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(nn-i,:,:) = array(nn-i,:,:) - array(nn+i,:,:)
          array(nn+i,:,:) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng
          array(nn+1-i,:,:) = array(nn+1-i,:,:) + array(nn+i,:,:)
          array(nn+i,:,:) = 0.0_num
        END DO
      END IF
    END IF

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    bc = bc_species(n)
    IF (y_min_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(:,i,:) = array(:,i,:) - array(:,-i,:)
          array(:,-i,:) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng-1
          array(:,i,:) = array(:,i,:) + array(:,1-i,:)
          array(:,1-i,:) = 0.0_num
        END DO
      END IF
    END IF

    n = n + 1
    bc = bc_species(n)
    IF (y_max_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(:,nn-i,:) = array(:,nn-i,:) - array(:,nn+i,:)
          array(:,nn+i,:) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng
          array(:,nn+1-i,:) = array(:,nn+1-i,:) + array(:,nn+i,:)
          array(:,nn+i,:) = 0.0_num
        END DO
      END IF
    END IF

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    bc = bc_species(n)
    IF (z_min_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(:,:,i) = array(:,:,i) - array(:,:,-i)
          array(:,:,-i) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng-1
          array(:,:,i) = array(:,:,i) + array(:,:,1-i)
          array(:,:,1-i) = 0.0_num
        END DO
      END IF
    END IF

    n = n + 1
    bc = bc_species(n)
    IF (z_max_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(:,:,nn-i) = array(:,:,nn-i) - array(:,:,nn+i)
          array(:,:,nn+i) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng
          array(:,:,nn+1-i) = array(:,:,nn+1-i) + array(:,:,nn+i)
          array(:,:,nn+i) = 0.0_num
        END DO
      END IF
    END IF

  END SUBROUTINE particle_reflection_bcs



  SUBROUTINE particle_periodic_bcs(array, ng, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: species
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: n, nn, sz, subarray, i
    INTEGER, DIMENSION(-1:1) :: neighbour_local
    INTEGER, DIMENSION(2*c_ndims) :: bc_species

    ! Transmit and sum all boundaries.
    ! Set neighbour to MPI_PROC_NULL if we don't need to transfer anything

    bc_species = bc_allspecies
    IF (PRESENT(species)) THEN
      DO i = 1, 2*c_ndims
        IF (bc_species(i) == c_bc_mixed) THEN
          bc_species(i) = species_list(species)%bc_particle(i)
        END IF
      END DO
    END IF

    sizes = SHAPE(array)
    starts = 1
    n = 0

    subsizes = sizes
    subsizes(n/2+1) = ng
    nn = sizes(n/2+1) - 2 * ng

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    ! Don't bother communicating non-periodic boundaries
    neighbour_local = neighbour(:,0,0)
    n = n + 1
    IF (x_min_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local(-1) = MPI_PROC_NULL
      END IF
    END IF
    n = n + 1
    IF (x_max_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local( 1) = MPI_PROC_NULL
      END IF
    END IF
    n = n - 2

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn+1,1-ng,1-ng), 1, subarray, &
        neighbour_local( 1), tag, temp, sz, mpireal, &
        neighbour_local(-1), tag, comm, status, errcode)

    n = n + 1
    array(1:ng,:,:) = array(1:ng,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour_local(-1), tag, temp, sz, mpireal, &
        neighbour_local( 1), tag, comm, status, errcode)

    n = n + 1
    array(nn+1-ng:nn,:,:) = array(nn+1-ng:nn,:,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes = sizes
    subsizes(n/2+1) = ng
    nn = sizes(n/2+1) - 2 * ng

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    ! Don't bother communicating non-periodic boundaries
    neighbour_local = neighbour(0,:,0)
    n = n + 1
    IF (y_min_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local(-1) = MPI_PROC_NULL
      END IF
    END IF
    n = n + 1
    IF (y_max_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local( 1) = MPI_PROC_NULL
      END IF
    END IF
    n = n - 2

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,nn+1,1-ng), 1, subarray, &
        neighbour_local( 1), tag, temp, sz, mpireal, &
        neighbour_local(-1), tag, comm, status, errcode)

    n = n + 1
    array(:,1:ng,:) = array(:,1:ng,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour_local(-1), tag, temp, sz, mpireal, &
        neighbour_local( 1), tag, comm, status, errcode)

    n = n + 1
    array(:,nn+1-ng:nn,:) = array(:,nn+1-ng:nn,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes = sizes
    subsizes(n/2+1) = ng
    nn = sizes(n/2+1) - 2 * ng

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    ! Don't bother communicating non-periodic boundaries
    neighbour_local = neighbour(0,0,:)
    n = n + 1
    IF (z_min_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local(-1) = MPI_PROC_NULL
      END IF
    END IF
    n = n + 1
    IF (z_max_boundary) THEN
      IF (bc_species(n) /= c_bc_periodic) THEN
        neighbour_local( 1) = MPI_PROC_NULL
      END IF
    END IF
    n = n - 2

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,nn+1), 1, subarray, &
        neighbour_local( 1), tag, temp, sz, mpireal, &
        neighbour_local(-1), tag, comm, status, errcode)

    n = n + 1
    array(:,:,1:ng) = array(:,:,1:ng) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour_local(-1), tag, temp, sz, mpireal, &
        neighbour_local( 1), tag, comm, status, errcode)

    n = n + 1
    array(:,:,nn+1-ng:nn) = array(:,:,nn+1-ng:nn) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    IF (PRESENT(species)) CALL particle_clear_bcs(array, ng)

  END SUBROUTINE particle_periodic_bcs



  SUBROUTINE particle_clear_bcs(array, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    INTEGER, DIMENSION(c_ndims) :: sizes
    INTEGER :: n, nn

    sizes = SHAPE(array)
    n = 0

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    array(:0,:,:) = 0.0_num
    n = n + 1
    array(nn+1:,:,:) = 0.0_num

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    array(:,:0,:) = 0.0_num
    n = n + 1
    array(:,nn+1:,:) = 0.0_num

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    array(:,:,:0) = 0.0_num
    n = n + 1
    array(:,:,nn+1:) = 0.0_num

  END SUBROUTINE particle_clear_bcs



  SUBROUTINE processor_summation_bcs(array, ng, flip_direction, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: flip_direction
    INTEGER, INTENT(IN), OPTIONAL :: species

    IF (PRESENT(species) .AND. .NOT. ANY(bc_allspecies == c_bc_mixed)) THEN
      RETURN
    END IF

    IF (.NOT. PRESENT(species) .AND. ANY(bc_allspecies == c_bc_mixed)) THEN
      RETURN
    END IF

    ! First apply reflecting boundary conditions
    CALL particle_reflection_bcs(array, ng, flip_direction, species)

    ! Next apply periodic and subdomain boundary conditions
    CALL particle_periodic_bcs(array, ng, species)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs

    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(ex, ng)
    CALL field_bc(ey, ng)
    CALL field_bc(ez, ng)

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_zero_gradient(ex, c_stagger_ex, i)
        CALL field_clamp_zero(ey, ng, c_stagger_ey, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      END IF
    END DO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_zero_gradient(ey, c_stagger_ey, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      END IF
    END DO

    DO i = c_bd_z_min, c_bd_z_max, c_bd_z_max - c_bd_z_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_clamp_zero(ey, ng, c_stagger_ey, i)
        CALL field_zero_gradient(ez, c_stagger_ez, i)
      END IF
    END DO

    DO i = 1, 2*c_ndims
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) == c_bc_clamp &
          .OR. bc_field(i) == c_bc_simple_laser &
          .OR. bc_field(i) == c_bc_simple_outflow) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_clamp_zero(ey, ng, c_stagger_ey, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      END IF

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) == c_bc_zero_gradient &
          .OR. bc_field(i) == c_bc_cpml_laser &
          .OR. bc_field(i) == c_bc_cpml_outflow) THEN
        CALL field_zero_gradient(ex, c_stagger_ex, i)
        CALL field_zero_gradient(ey, c_stagger_ey, i)
        CALL field_zero_gradient(ez, c_stagger_ez, i)
      END IF
    END DO

  END SUBROUTINE efield_bcs



  SUBROUTINE bfield_bcs(mpi_only)

    LOGICAL, INTENT(IN) :: mpi_only
    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(bx, ng)
    CALL field_bc(by, ng)
    CALL field_bc(bz, ng)

    IF (mpi_only) RETURN

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(bx, ng, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      END IF
    END DO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_clamp_zero(by, ng, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      END IF
    END DO

    DO i = c_bd_z_min, c_bd_z_max, c_bd_z_max - c_bd_z_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_clamp_zero(bz, ng, c_stagger_bz, i)
      END IF
    END DO

    DO i = 1, 2*c_ndims
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) == c_bc_clamp &
          .OR. bc_field(i) == c_bc_simple_laser &
          .OR. bc_field(i) == c_bc_simple_outflow) THEN
        CALL field_clamp_zero(bx, ng, c_stagger_bx, i)
        CALL field_clamp_zero(by, ng, c_stagger_by, i)
        CALL field_clamp_zero(bz, ng, c_stagger_bz, i)
      END IF

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) == c_bc_zero_gradient &
          .OR. bc_field(i) == c_bc_cpml_laser &
          .OR. bc_field(i) == c_bc_cpml_outflow) THEN
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      END IF
    END DO

  END SUBROUTINE bfield_bcs



  SUBROUTINE bfield_final_bcs

    INTEGER :: i

    CALL update_laser_omegas
    CALL bfield_bcs(.FALSE.)

    IF (x_min_boundary) THEN
      i = c_bd_x_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_x_min
    END IF

    IF (x_max_boundary) THEN
      i = c_bd_x_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_x_max
    END IF

    IF (y_min_boundary) THEN
      i = c_bd_y_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_y_min
    END IF

    IF (y_max_boundary) THEN
      i = c_bd_y_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_y_max
    END IF

    IF (z_min_boundary) THEN
      i = c_bd_z_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_z_min
    END IF

    IF (z_max_boundary) THEN
      i = c_bd_z_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_z_max
    END IF

    CALL bfield_bcs(.TRUE.)

  END SUBROUTINE bfield_final_bcs



  SUBROUTINE setup_bc_lists

    INTEGER(i8) :: ispecies, ipart
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    REAL(num) :: bnd_x_min, bnd_x_max
    REAL(num) :: bnd_y_min, bnd_y_max
    REAL(num) :: bnd_z_min, bnd_z_max
    TYPE(particle), POINTER :: current
    TYPE(particle_pointer_list), POINTER :: bnd_part_last, bnd_part_next

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head

      IF (species_list(ispecies)%attached_list%count == 0) CYCLE

      bc_species = species_list(ispecies)%bc_particle
      IF ((bc_species(c_bd_x_min) == c_bc_thermal &
          .OR. bc_field(c_bd_x_min) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) &
          .AND. x_min_boundary) THEN
        bnd_x_min = x_min_outer
      ELSE
        bnd_x_min = x_min_local
      END IF
      IF ((bc_species(c_bd_x_max) == c_bc_thermal &
          .OR. bc_field(c_bd_x_max) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) &
          .AND. x_max_boundary) THEN
        bnd_x_max = x_max_outer
      ELSE
        bnd_x_max = x_max_local
      END IF
      IF ((bc_species(c_bd_y_min) == c_bc_thermal &
          .OR. bc_field(c_bd_y_min) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_y_min) == c_bc_cpml_outflow) &
          .AND. y_min_boundary) THEN
        bnd_y_min = y_min_outer
      ELSE
        bnd_y_min = y_min_local
      END IF
      IF ((bc_species(c_bd_y_max) == c_bc_thermal &
          .OR. bc_field(c_bd_y_max) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_y_max) == c_bc_cpml_outflow) &
          .AND. y_max_boundary) THEN
        bnd_y_max = y_max_outer
      ELSE
        bnd_y_max = y_max_local
      END IF
      IF ((bc_species(c_bd_z_min) == c_bc_thermal &
          .OR. bc_field(c_bd_z_min) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_z_min) == c_bc_cpml_outflow) &
          .AND. z_min_boundary) THEN
        bnd_z_min = z_min_outer
      ELSE
        bnd_z_min = z_min_local
      END IF
      IF ((bc_species(c_bd_z_max) == c_bc_thermal &
          .OR. bc_field(c_bd_z_max) == c_bc_cpml_laser &
          .OR. bc_field(c_bd_z_max) == c_bc_cpml_outflow) &
          .AND.z_max_boundary) THEN
        bnd_z_max = z_max_outer
      ELSE
        bnd_z_max = z_max_local
      END IF

      ALLOCATE(species_list(ispecies)%boundary_particles)
      NULLIFY(species_list(ispecies)%boundary_particles%particle)
      NULLIFY(species_list(ispecies)%boundary_particles%next)
      NULLIFY(bnd_part_next)
      bnd_part_last => species_list(ispecies)%boundary_particles

      DO ipart = 1, species_list(ispecies)%attached_list%count
        ! Move particle to boundary candidate list
        IF (current%part_pos(1) < bnd_x_min &
            .OR. current%part_pos(1) > bnd_x_max &
            .OR. current%part_pos(2) < bnd_y_min &
            .OR. current%part_pos(2) > bnd_y_max &
            .OR. current%part_pos(3) < bnd_z_min &
            .OR. current%part_pos(3) > bnd_z_max) THEN
          ALLOCATE(bnd_part_next)
          bnd_part_next%particle => current
          bnd_part_last%next => bnd_part_next
          bnd_part_last => bnd_part_next
        END IF
        current => current%next
      END DO

      ! Boundary list head contains no particle
      bnd_part_last => species_list(ispecies)%boundary_particles
      species_list(ispecies)%boundary_particles &
          => species_list(ispecies)%boundary_particles%next
      DEALLOCATE(bnd_part_last)
      ! Final particle should have null 'next' ptr
      IF (ASSOCIATED(bnd_part_next)) NULLIFY(bnd_part_next%next)
    END DO

  END SUBROUTINE setup_bc_lists



  SUBROUTINE particle_bcs

    TYPE(particle_pointer_list), POINTER :: bnd_part, bnd_part_last
    TYPE(particle), POINTER :: cur
    TYPE(particle_list), DIMENSION(-1:1,-1:1,-1:1) :: send, recv
    INTEGER :: xbd, ybd, zbd
    INTEGER(i8) :: ixp, iyp, izp
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    LOGICAL :: out_of_bounds
    INTEGER :: sgn, bc, ispecies, i, ix, iy, iz
    INTEGER :: cell_x, cell_y, cell_z
    REAL(num), DIMENSION(-1:1) :: gx, gy, gz
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    REAL(num) :: cf2, temp(3)
    REAL(num) :: part_pos
    REAL(num) :: x_shift, y_shift, z_shift

    x_shift = length_x + 2.0_num * dx * REAL(cpml_thickness, num)
    y_shift = length_y + 2.0_num * dy * REAL(cpml_thickness, num)
    z_shift = length_z + 2.0_num * dz * REAL(cpml_thickness, num)

    DO ispecies = 1, n_species
      bnd_part => species_list(ispecies)%boundary_particles
      NULLIFY(bnd_part_last)

      bc_species = species_list(ispecies)%bc_particle

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            CALL create_empty_partlist(send(ix, iy, iz))
            CALL create_empty_partlist(recv(ix, iy, iz))
          END DO
        END DO
      END DO

      DO WHILE (ASSOCIATED(bnd_part))
        bnd_part_last => bnd_part

        cur => bnd_part%particle
        bnd_part => bnd_part%next
        DEALLOCATE(bnd_part_last)

        xbd = 0
        ybd = 0
        zbd = 0
        out_of_bounds = .FALSE.

        part_pos = cur%part_pos(1)
        sgn = -1
        IF (bc_field(c_bd_x_min) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) THEN
          IF (x_min_boundary) THEN
            ! Particle has left the system
            IF (part_pos < x_min_outer) THEN
              xbd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos < x_min_local) xbd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos < x_min_local) THEN
            xbd = sgn
            bc = bc_species(c_bd_x_min)
            IF (bc == c_bc_reflect) THEN
              IF (x_min_boundary) THEN
                xbd = 0
                cur%part_pos(1) = 2.0_num * x_min - part_pos
                cur%part_p(1) = -cur%part_p(1)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (x_min_boundary) THEN
                cur%part_pos(1) = part_pos - sgn * x_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos < x_min_outer) THEN
                xbd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
                cell_y = FLOOR(cell_y_r + 0.5_num)
                cell_frac_y = REAL(cell_y, num) - cell_y_r
                cell_y = cell_y + 1

                cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
                cell_z = FLOOR(cell_z_r + 0.5_num)
                cell_frac_z = REAL(cell_z, num) - cell_z_r
                cell_z = cell_z + 1

                cf2 = cell_frac_y**2
                gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
                gy( 0) = 0.75_num - cf2
                gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

                cf2 = cell_frac_z**2
                gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
                gz( 0) = 0.75_num - cf2
                gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iz = -1, 1
                    DO iy = -1, 1
                      temp(i) = temp(i) + gy(iy) * gz(iz) &
                          * species_list(ispecies)&
                          %ext_temp_x_min(cell_y+iy, cell_z+iz, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                ! y-direction
                i = 2
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! z-direction
                i = 3
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                cur%part_pos(1) = 2.0_num * x_min_outer - part_pos

              ELSE IF (x_min_boundary) THEN
                xbd = 0
              END IF
            ELSE
              IF (part_pos < x_min_outer) THEN
                ! Default to open boundary conditions - remove particle
                xbd = 0
                out_of_bounds = .TRUE.
              ELSE IF (x_min_boundary) THEN
                xbd = 0
              END IF
            END IF
          END IF
        END IF

        sgn = 1
        IF (bc_field(c_bd_x_max) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) THEN
          IF (x_max_boundary) THEN
            ! Particle has left the system
            IF (part_pos >= x_max_outer) THEN
              xbd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos >= x_max_local) xbd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos >= x_max_local) THEN
            xbd = sgn
            bc = bc_species(c_bd_x_max)
            IF (bc == c_bc_reflect) THEN
              IF (x_max_boundary) THEN
                xbd = 0
                cur%part_pos(1) = 2.0_num * x_max - part_pos
                cur%part_p(1) = -cur%part_p(1)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (x_max_boundary) THEN
                cur%part_pos(1) = part_pos - sgn * x_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos >= x_max_outer) THEN
                xbd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
                cell_y = FLOOR(cell_y_r + 0.5_num)
                cell_frac_y = REAL(cell_y, num) - cell_y_r
                cell_y = cell_y + 1

                cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
                cell_z = FLOOR(cell_z_r + 0.5_num)
                cell_frac_z = REAL(cell_z, num) - cell_z_r
                cell_z = cell_z + 1

                cf2 = cell_frac_y**2
                gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
                gy( 0) = 0.75_num - cf2
                gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

                cf2 = cell_frac_z**2
                gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
                gz( 0) = 0.75_num - cf2
                gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iz = -1, 1
                    DO iy = -1, 1
                      temp(i) = temp(i) + gy(iy) * gz(iz) &
                          * species_list(ispecies)&
                          %ext_temp_x_max(cell_y+iy, cell_z+iz, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                ! y-direction
                i = 2
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! z-direction
                i = 3
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                cur%part_pos(1) = 2.0_num * x_max_outer - part_pos

              ELSE IF (x_max_boundary) THEN
                xbd = 0
              END IF
            ELSE
              IF (part_pos >= x_max_outer) THEN
                ! Default to open boundary conditions - remove particle
                xbd = 0
                out_of_bounds = .TRUE.
              ELSE IF (x_max_boundary) THEN
                xbd = 0
              END IF
            END IF
          END IF
        END IF

        part_pos = cur%part_pos(2)
        sgn = -1
        IF (bc_field(c_bd_y_min) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_y_min) == c_bc_cpml_outflow) THEN
          IF (y_min_boundary) THEN
            ! Particle has left the system
            IF (part_pos < y_min_outer) THEN
              ybd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos < y_min_local) ybd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos < y_min_local) THEN
            ybd = sgn
            bc = bc_species(c_bd_y_min)
            IF (bc == c_bc_reflect) THEN
              IF (y_min_boundary) THEN
                ybd = 0
                cur%part_pos(2) = 2.0_num * y_min - part_pos
                cur%part_p(2) = -cur%part_p(2)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (y_min_boundary) THEN
                cur%part_pos(2) = part_pos - sgn * y_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos < y_min_outer) THEN
                ybd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
                cell_x = FLOOR(cell_x_r + 0.5_num)
                cell_frac_x = REAL(cell_x, num) - cell_x_r
                cell_x = cell_x + 1

                cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
                cell_z = FLOOR(cell_z_r + 0.5_num)
                cell_frac_z = REAL(cell_z, num) - cell_z_r
                cell_z = cell_z + 1

                cf2 = cell_frac_x**2
                gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
                gx( 0) = 0.75_num - cf2
                gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

                cf2 = cell_frac_z**2
                gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
                gz( 0) = 0.75_num - cf2
                gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iz = -1, 1
                    DO ix = -1, 1
                      temp(i) = temp(i) + gx(ix) * gz(iz) &
                          * species_list(ispecies)&
                          %ext_temp_y_min(cell_x+ix, cell_z+iz, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! y-direction
                i = 2
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                ! z-direction
                i = 3
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                cur%part_pos(2) = 2.0_num * y_min_outer - part_pos

              ELSE IF (y_min_boundary) THEN
                ybd = 0
              END IF
            ELSE
              IF (part_pos < y_min_outer) THEN
                ! Default to open boundary conditions - remove particle
                ybd = 0
                out_of_bounds = .TRUE.
              ELSE IF (y_min_boundary) THEN
                ybd = 0
              END IF
            END IF
          END IF
        END IF

        sgn = 1
        IF (bc_field(c_bd_y_max) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_y_max) == c_bc_cpml_outflow) THEN
          IF (y_max_boundary) THEN
            ! Particle has left the system
            IF (part_pos >= y_max_outer) THEN
              ybd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos >= y_max_local) ybd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos >= y_max_local) THEN
            ybd = sgn
            bc = bc_species(c_bd_y_max)
            IF (bc == c_bc_reflect) THEN
              IF (y_max_boundary) THEN
                ybd = 0
                cur%part_pos(2) = 2.0_num * y_max - part_pos
                cur%part_p(2) = -cur%part_p(2)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (y_max_boundary) THEN
                cur%part_pos(2) = part_pos - sgn * y_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos >= y_max_outer) THEN
                ybd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
                cell_x = FLOOR(cell_x_r + 0.5_num)
                cell_frac_x = REAL(cell_x, num) - cell_x_r
                cell_x = cell_x + 1

                cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
                cell_z = FLOOR(cell_z_r + 0.5_num)
                cell_frac_z = REAL(cell_z, num) - cell_z_r
                cell_z = cell_z + 1

                cf2 = cell_frac_x**2
                gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
                gx( 0) = 0.75_num - cf2
                gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

                cf2 = cell_frac_z**2
                gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
                gz( 0) = 0.75_num - cf2
                gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iz = -1, 1
                    DO ix = -1, 1
                      temp(i) = temp(i) + gx(ix) * gz(iz) &
                          * species_list(ispecies)&
                          %ext_temp_y_max(cell_x+ix, cell_z+iz, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! y-direction
                i = 2
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                ! z-direction
                i = 3
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                cur%part_pos(2) = 2.0_num * y_max_outer - part_pos

              ELSE IF (y_max_boundary) THEN
                ybd = 0
              END IF
            ELSE
              IF (part_pos >= y_max_outer) THEN
                ! Default to open boundary conditions - remove particle
                ybd = 0
                out_of_bounds = .TRUE.
              ELSE IF (y_max_boundary) THEN
                ybd = 0
              END IF
            END IF
          END IF
        END IF

        part_pos = cur%part_pos(3)
        sgn = -1
        IF (bc_field(c_bd_z_min) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_z_min) == c_bc_cpml_outflow) THEN
          IF (z_min_boundary) THEN
            ! Particle has left the system
            IF (part_pos < z_min_outer) THEN
              zbd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos < z_min_local) zbd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos < z_min_local) THEN
            zbd = sgn
            bc = bc_species(c_bd_z_min)
            IF (bc == c_bc_reflect) THEN
              IF (z_min_boundary) THEN
                zbd = 0
                cur%part_pos(3) = 2.0_num * z_min - part_pos
                cur%part_p(3) = -cur%part_p(3)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (z_min_boundary) THEN
                cur%part_pos(3) = part_pos - sgn * z_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos < z_min_outer) THEN
                zbd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
                cell_x = FLOOR(cell_x_r + 0.5_num)
                cell_frac_x = REAL(cell_x, num) - cell_x_r
                cell_x = cell_x + 1

                cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
                cell_y = FLOOR(cell_y_r + 0.5_num)
                cell_frac_y = REAL(cell_y, num) - cell_y_r
                cell_y = cell_y + 1

                cf2 = cell_frac_x**2
                gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
                gx( 0) = 0.75_num - cf2
                gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

                cf2 = cell_frac_y**2
                gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
                gy( 0) = 0.75_num - cf2
                gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iy = -1, 1
                    DO ix = -1, 1
                      temp(i) = temp(i) + gx(ix) * gy(iy) &
                          * species_list(ispecies)&
                          %ext_temp_z_min(cell_x+ix, cell_y+iy, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! y-direction
                i = 2
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! z-direction
                i = 3
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                cur%part_pos(3) = 2.0_num * z_min_outer - part_pos

              ELSE IF (z_min_boundary) THEN
                zbd = 0
              END IF
            ELSE
              IF (part_pos < z_min_outer) THEN
                ! Default to open boundary conditions - remove particle
                zbd = 0
                out_of_bounds = .TRUE.
              ELSE IF (z_min_boundary) THEN
                zbd = 0
              END IF
            END IF
          END IF
        END IF

        sgn = 1
        IF (bc_field(c_bd_z_max) == c_bc_cpml_laser &
            .OR. bc_field(c_bd_z_max) == c_bc_cpml_outflow) THEN
          IF (z_max_boundary) THEN
            ! Particle has left the system
            IF (part_pos >= z_max_outer) THEN
              zbd = 0
              out_of_bounds = .TRUE.
            END IF
          ELSE
            ! Particle has left this processor
            IF (part_pos >= z_max_local) zbd = sgn
          END IF
        ELSE
          ! Particle has left this processor
          IF (part_pos >= z_max_local) THEN
            zbd = sgn
            bc = bc_species(c_bd_z_max)
            IF (bc == c_bc_reflect) THEN
              IF (z_max_boundary) THEN
                zbd = 0
                cur%part_pos(3) = 2.0_num * z_max - part_pos
                cur%part_p(3) = -cur%part_p(3)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (z_max_boundary) THEN
                cur%part_pos(3) = part_pos - sgn * z_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos >= z_max_outer) THEN
                zbd = 0
                ! Always use the triangle particle weighting for simplicity
                cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
                cell_x = FLOOR(cell_x_r + 0.5_num)
                cell_frac_x = REAL(cell_x, num) - cell_x_r
                cell_x = cell_x + 1

                cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
                cell_y = FLOOR(cell_y_r + 0.5_num)
                cell_frac_y = REAL(cell_y, num) - cell_y_r
                cell_y = cell_y + 1

                cf2 = cell_frac_x**2
                gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
                gx( 0) = 0.75_num - cf2
                gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

                cf2 = cell_frac_y**2
                gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
                gy( 0) = 0.75_num - cf2
                gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

                DO i = 1, 3
                  temp(i) = 0.0_num
                  DO iy = -1, 1
                    DO ix = -1, 1
                      temp(i) = temp(i) + gx(ix) * gy(iy) &
                          * species_list(ispecies)&
                          %ext_temp_z_max(cell_x+ix, cell_y+iy, i)
                    END DO
                  END DO
                END DO

                CALL id_registry%delete_all(cur)

                ! x-direction
                i = 1
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! y-direction
                i = 2
                cur%part_p(i) = momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num)

                ! z-direction
                i = 3
                cur%part_p(i) = flux_momentum_from_temperature(&
                    species_list(ispecies)%mass, temp(i), 0.0_num, &
                    -REAL(sgn, num))

                cur%part_pos(3) = 2.0_num * z_max_outer - part_pos

              ELSE IF (z_max_boundary) THEN
                zbd = 0
              END IF
            ELSE
              IF (part_pos >= z_max_outer) THEN
                ! Default to open boundary conditions - remove particle
                zbd = 0
                out_of_bounds = .TRUE.
              ELSE IF (z_max_boundary) THEN
                zbd = 0
              END IF
            END IF
          END IF
        END IF

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          IF (track_ejected_particles) THEN
            CALL add_particle_to_partlist(&
                ejected_list(ispecies)%attached_list, cur)
          ELSE
            CALL destroy_particle(cur)
          END IF
        ELSE IF (ABS(xbd) + ABS(ybd) + ABS(zbd) > 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd, ybd, zbd), cur)
        END IF
      END DO

      ! swap Particles
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            ixp = -ix
            iyp = -iy
            izp = -iz
            CALL partlist_sendrecv(send(ix, iy, iz), recv(ixp, iyp, izp), &
                neighbour(ix, iy, iz), neighbour(ixp, iyp, izp))
            CALL append_partlist(species_list(ispecies)%attached_list, &
                recv(ixp, iyp, izp))
          END DO
        END DO
      END DO

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            CALL destroy_partlist(send(ix, iy, iz))
            CALL destroy_partlist(recv(ix, iy, iz))
          END DO
        END DO
      END DO

      IF (ASSOCIATED(bnd_part_last)) DEALLOCATE(bnd_part_last)
      IF (ASSOCIATED(species_list(ispecies)%boundary_particles)) THEN
        NULLIFY(species_list(ispecies)%boundary_particles)
      END IF
    END DO

  END SUBROUTINE particle_bcs



  SUBROUTINE current_bcs(species)

    INTEGER, INTENT(IN), OPTIONAL :: species

    ! Domain is decomposed. Just add currents at edges
    CALL processor_summation_bcs(jx, jng, c_dir_x, species)
    CALL processor_summation_bcs(jy, jng, c_dir_y, species)
    CALL processor_summation_bcs(jz, jng, c_dir_z, species)

  END SUBROUTINE current_bcs



  SUBROUTINE set_cpml_helpers(nx, nx_global_min, nx_global_max, &
      ny, ny_global_min, ny_global_max, nz, nz_global_min, nz_global_max)

    INTEGER, INTENT(IN) :: nx, nx_global_min, nx_global_max
    INTEGER, INTENT(IN) :: ny, ny_global_min, ny_global_max
    INTEGER, INTENT(IN) :: nz, nz_global_min, nz_global_max
    INTEGER :: i
    INTEGER :: ix, ix_glob
    INTEGER :: iy, iy_glob
    INTEGER :: iz, iz_glob
    INTEGER, PARAMETER :: cpml_m = 3
    INTEGER, PARAMETER :: cpml_ma = 1
    REAL(num) :: x_pos, x_pos_m, x_pos_ma
    REAL(num) :: y_pos, y_pos_m, y_pos_ma
    REAL(num) :: z_pos, z_pos_m, z_pos_ma
    REAL(num) :: cpml_sigma_maxval

    ALLOCATE(cpml_kappa_ex(1-ng:nx+ng), cpml_kappa_bx(1-ng:nx+ng))
    ALLOCATE(cpml_a_ex(1-ng:nx+ng), cpml_a_bx(1-ng:nx+ng))
    ALLOCATE(cpml_sigma_ex(1-ng:nx+ng), cpml_sigma_bx(1-ng:nx+ng))

    ALLOCATE(cpml_kappa_ey(1-ng:ny+ng), cpml_kappa_by(1-ng:ny+ng))
    ALLOCATE(cpml_a_ey(1-ng:ny+ng), cpml_a_by(1-ng:ny+ng))
    ALLOCATE(cpml_sigma_ey(1-ng:ny+ng), cpml_sigma_by(1-ng:ny+ng))

    ALLOCATE(cpml_kappa_ez(1-ng:nz+ng), cpml_kappa_bz(1-ng:nz+ng))
    ALLOCATE(cpml_a_ez(1-ng:nz+ng), cpml_a_bz(1-ng:nz+ng))
    ALLOCATE(cpml_sigma_ez(1-ng:nz+ng), cpml_sigma_bz(1-ng:nz+ng))

    cpml_kappa_ex = 1.0_num
    cpml_kappa_bx = 1.0_num

    cpml_a_ex = 0.0_num
    cpml_sigma_ex = 0.0_num
    cpml_a_bx = 0.0_num
    cpml_sigma_bx = 0.0_num

    cpml_kappa_ey = 1.0_num
    cpml_kappa_by = 1.0_num

    cpml_a_ey = 0.0_num
    cpml_sigma_ey = 0.0_num
    cpml_a_by = 0.0_num
    cpml_sigma_by = 0.0_num

    cpml_kappa_ez = 1.0_num
    cpml_kappa_bz = 1.0_num

    cpml_a_ez = 0.0_num
    cpml_sigma_ez = 0.0_num
    cpml_a_bz = 0.0_num
    cpml_sigma_bz = 0.0_num

    cpml_sigma_maxval = cpml_sigma_max * c * 0.8_num * (cpml_m + 1.0_num) / dx

    ! ============= x_min boundary =============

    i = c_bd_x_min
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_x_min_start = nx+1
      cpml_x_min_end = 0
      cpml_x_min_offset = 0

      IF (nx_global_min <= cpml_thickness) THEN
        cpml_x_min = .TRUE.
        cpml_x_min_start = 1 ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nx_global_max >= cpml_thickness) THEN
          ! in local grid coordinates
          ! global -> local: ixl = ixg - nx_global_min + 1
          cpml_x_min_end = cpml_thickness - nx_global_min + 1
          cpml_x_min_offset = cpml_thickness - nx_global_min + 1
        ELSE
          cpml_x_min_end = nx ! in local grid coordinates
          cpml_x_min_offset = cpml_thickness
        END IF

        DO ix = cpml_x_min_start,cpml_x_min_end
          ! runs from 1 to cpml_thickness in global coordinates
          ! local -> global: ixg = ixl + nx_global_min - 1
          ix_glob = ix + nx_global_min - 1

          ! runs from 1.0 to nearly 0.0 (actually 0.0 at cpml_thickness+1)
          x_pos = 1.0_num - REAL(ix_glob-1,num) / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_ex(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_ex(ix) = cpml_sigma_maxval * x_pos_m
          cpml_a_ex(ix) = cpml_a_max * x_pos_ma

          ! runs from nearly 1.0 to nearly 0.0 on the half intervals
          ! 1.0 at ix_glob=1-1/2 and 0.0 at ix_glob=cpml_thickness+1/2
          x_pos = 1.0_num - (REAL(ix_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_bx(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_bx(ix) = cpml_sigma_maxval * x_pos_m
          cpml_a_bx(ix) = cpml_a_max * x_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nx_global_min <= cpml_thickness + fng + 1 &
          .AND. nx_global_max >= cpml_thickness + fng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_x_min_laser_idx = cpml_thickness + fng + 1 - nx_global_min
      END IF
    END IF

    ! ============= x_max boundary =============

    i = c_bd_x_max
    ! Same as x_min using the transformation ix -> nx_global - ix + 1
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_x_max_start = nx+1
      cpml_x_max_end = 0
      cpml_x_max_offset = 0

      IF (nx_global_max >= nx_global - cpml_thickness + 1) THEN
        cpml_x_max = .TRUE.
        cpml_x_max_end = nx ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nx_global_min <= nx_global - cpml_thickness + 1) THEN
          ! in local grid coordinates
          ! global -> local: ixl = ixg - nx_global_min + 1
          cpml_x_max_start = nx_global - cpml_thickness + 1 - nx_global_min + 1
          cpml_x_max_offset = cpml_thickness - nx_global + nx_global_max
        ELSE
          cpml_x_max_start = 1 ! in local grid coordinates
          cpml_x_max_offset = cpml_thickness
        END IF

        DO ix = cpml_x_max_start,cpml_x_max_end
          ! runs from cpml_thickness to 1 in global coordinates
          ! local -> global: ixg = ixl + nx_global_min - 1
          ix_glob = nx_global - (ix + nx_global_min - 1) + 1

          ! runs from nearly 0.0 (actually 0.0 at cpml_thickness+1) to 1.0
          x_pos = 1.0_num - REAL(ix_glob-1,num) / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_ex(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_ex(ix) = cpml_sigma_maxval * x_pos_m
          cpml_a_ex(ix) = cpml_a_max * x_pos_ma

          ! runs from nearly 0.0 to nearly 1.0 on the half intervals
          ! 0.0 at ix_glob=cpml_thickness+1/2 and 1.0 at ix_glob=1-1/2
          x_pos = 1.0_num - (REAL(ix_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_bx(ix-1) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_bx(ix-1) = cpml_sigma_maxval * x_pos_m
          cpml_a_bx(ix-1) = cpml_a_max * x_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nx_global_min <= nx_global - cpml_thickness - fng + 2 &
          .AND. nx_global_max >= nx_global - cpml_thickness - fng + 2) THEN
        add_laser(i) = .TRUE.
        cpml_x_max_laser_idx = &
            nx_global - cpml_thickness - fng + 2 - nx_global_min
      END IF
    END IF

    ! ============= y_min boundary =============

    i = c_bd_y_min
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_y_min_start = ny+1
      cpml_y_min_end = 0
      cpml_y_min_offset = 0

      IF (ny_global_min <= cpml_thickness) THEN
        cpml_y_min = .TRUE.
        cpml_y_min_start = 1 ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (ny_global_max >= cpml_thickness) THEN
          ! in local grid coordinates
          ! global -> local: iyl = iyg - ny_global_min + 1
          cpml_y_min_end = cpml_thickness - ny_global_min + 1
          cpml_y_min_offset = cpml_thickness - ny_global_min + 1
        ELSE
          cpml_y_min_end = ny ! in local grid coordinates
          cpml_y_min_offset = cpml_thickness
        END IF

        DO iy = cpml_y_min_start,cpml_y_min_end
          ! runs from 1 to cpml_thickness in global coordinates
          ! local -> global: iyg = iyl + ny_global_min - 1
          iy_glob = iy + ny_global_min - 1

          ! runs from 1.0 to nearly 0.0 (actually 0.0 at cpml_thickness+1)
          y_pos = 1.0_num - REAL(iy_glob-1,num) / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_ey(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_ey(iy) = cpml_sigma_maxval * y_pos_m
          cpml_a_ey(iy) = cpml_a_max * y_pos_ma

          ! runs from nearly 1.0 to nearly 0.0 on the half intervals
          ! 1.0 at iy_glob=1-1/2 and 0.0 at iy_glob=cpml_thickness+1/2
          y_pos = 1.0_num - (REAL(iy_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_by(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_by(iy) = cpml_sigma_maxval * y_pos_m
          cpml_a_by(iy) = cpml_a_max * y_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (ny_global_min <= cpml_thickness + fng + 1 &
          .AND. ny_global_max >= cpml_thickness + fng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_y_min_laser_idx = cpml_thickness + fng + 1 - ny_global_min
      END IF
    END IF

    ! ============= y_max boundary =============

    i = c_bd_y_max
    ! Same as y_min using the transformation iy -> ny_global - iy + 1
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_y_max_start = ny+1
      cpml_y_max_end = 0
      cpml_y_max_offset = 0

      IF (ny_global_max >= ny_global - cpml_thickness + 1) THEN
        cpml_y_max = .TRUE.
        cpml_y_max_end = ny ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (ny_global_min <= ny_global - cpml_thickness + 1) THEN
          ! in local grid coordinates
          ! global -> local: iyl = iyg - ny_global_min + 1
          cpml_y_max_start = ny_global - cpml_thickness + 1 - ny_global_min + 1
          cpml_y_max_offset = cpml_thickness - ny_global + ny_global_max
        ELSE
          cpml_y_max_start = 1 ! in local grid coordinates
          cpml_y_max_offset = cpml_thickness
        END IF

        DO iy = cpml_y_max_start,cpml_y_max_end
          ! runs from cpml_thickness to 1 in global coordinates
          ! local -> global: iyg = iyl + ny_global_min - 1
          iy_glob = ny_global - (iy + ny_global_min - 1) + 1

          ! runs from nearly 0.0 (actually 0.0 at cpml_thickness+1) to 1.0
          y_pos = 1.0_num - REAL(iy_glob-1,num) / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_ey(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_ey(iy) = cpml_sigma_maxval * y_pos_m
          cpml_a_ey(iy) = cpml_a_max * y_pos_ma

          ! runs from nearly 0.0 to nearly 1.0 on the half intervals
          ! 0.0 at iy_glob=cpml_thickness+1/2 and 1.0 at iy_glob=1-1/2
          y_pos = 1.0_num - (REAL(iy_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_by(iy-1) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_by(iy-1) = cpml_sigma_maxval * y_pos_m
          cpml_a_by(iy-1) = cpml_a_max * y_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (ny_global_min <= ny_global - cpml_thickness - fng + 2 &
          .AND. ny_global_max >= ny_global - cpml_thickness - fng + 2) THEN
        add_laser(i) = .TRUE.
        cpml_y_max_laser_idx = &
            ny_global - cpml_thickness - fng + 2 - ny_global_min
      END IF
    END IF

    ! ============= z_min boundary =============

    i = c_bd_z_min
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_z_min_start = nz+1
      cpml_z_min_end = 0
      cpml_z_min_offset = 0

      IF (nz_global_min <= cpml_thickness) THEN
        cpml_z_min = .TRUE.
        cpml_z_min_start = 1 ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nz_global_max >= cpml_thickness) THEN
          ! in local grid coordinates
          ! global -> local: izl = izg - nz_global_min + 1
          cpml_z_min_end = cpml_thickness - nz_global_min + 1
          cpml_z_min_offset = cpml_thickness - nz_global_min + 1
        ELSE
          cpml_z_min_end = nz ! in local grid coordinates
          cpml_z_min_offset = cpml_thickness
        END IF

        DO iz = cpml_z_min_start,cpml_z_min_end
          ! runs from 1 to cpml_thickness in global coordinates
          ! local -> global: izg = izl + nz_global_min - 1
          iz_glob = iz + nz_global_min - 1

          ! runs from 1.0 to nearly 0.0 (actually 0.0 at cpml_thickness+1)
          z_pos = 1.0_num - REAL(iz_glob-1,num) / REAL(cpml_thickness,num)
          z_pos_m = z_pos**cpml_m
          z_pos_ma = (1.0_num - z_pos)**cpml_ma

          cpml_kappa_ez(iz) = 1.0_num + (cpml_kappa_max - 1.0_num) * z_pos_m
          cpml_sigma_ez(iz) = cpml_sigma_maxval * z_pos_m
          cpml_a_ez(iz) = cpml_a_max * z_pos_ma

          ! runs from nearly 1.0 to nearly 0.0 on the half intervals
          ! 1.0 at iz_glob=1-1/2 and 0.0 at iz_glob=cpml_thickness+1/2
          z_pos = 1.0_num - (REAL(iz_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          z_pos_m = z_pos**cpml_m
          z_pos_ma = (1.0_num - z_pos)**cpml_ma

          cpml_kappa_bz(iz) = 1.0_num + (cpml_kappa_max - 1.0_num) * z_pos_m
          cpml_sigma_bz(iz) = cpml_sigma_maxval * z_pos_m
          cpml_a_bz(iz) = cpml_a_max * z_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nz_global_min <= cpml_thickness + fng + 1 &
          .AND. nz_global_max >= cpml_thickness + fng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_z_min_laser_idx = cpml_thickness + fng + 1 - nz_global_min
      END IF
    END IF

    ! ============= z_max boundary =============

    i = c_bd_z_max
    ! Same as z_min using the transformation iz -> nz_global - iz + 1
    IF (bc_field(i) == c_bc_cpml_laser &
        .OR. bc_field(i) == c_bc_cpml_outflow) THEN
      cpml_z_max_start = nz+1
      cpml_z_max_end = 0
      cpml_z_max_offset = 0

      IF (nz_global_max >= nz_global - cpml_thickness + 1) THEN
        cpml_z_max = .TRUE.
        cpml_z_max_end = nz ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nz_global_min <= nz_global - cpml_thickness + 1) THEN
          ! in local grid coordinates
          ! global -> local: izl = izg - nz_global_min + 1
          cpml_z_max_start = nz_global - cpml_thickness + 1 - nz_global_min + 1
          cpml_z_max_offset = cpml_thickness - nz_global + nz_global_max
        ELSE
          cpml_z_max_start = 1 ! in local grid coordinates
          cpml_z_max_offset = cpml_thickness
        END IF

        DO iz = cpml_z_max_start,cpml_z_max_end
          ! runs from cpml_thickness to 1 in global coordinates
          ! local -> global: izg = izl + nz_global_min - 1
          iz_glob = nz_global - (iz + nz_global_min - 1) + 1

          ! runs from nearly 0.0 (actually 0.0 at cpml_thickness+1) to 1.0
          z_pos = 1.0_num - REAL(iz_glob-1,num) / REAL(cpml_thickness,num)
          z_pos_m = z_pos**cpml_m
          z_pos_ma = (1.0_num - z_pos)**cpml_ma

          cpml_kappa_ez(iz) = 1.0_num + (cpml_kappa_max - 1.0_num) * z_pos_m
          cpml_sigma_ez(iz) = cpml_sigma_maxval * z_pos_m
          cpml_a_ez(iz) = cpml_a_max * z_pos_ma

          ! runs from nearly 0.0 to nearly 1.0 on the half intervals
          ! 0.0 at iz_glob=cpml_thickness+1/2 and 1.0 at iz_glob=1-1/2
          z_pos = 1.0_num - (REAL(iz_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          z_pos_m = z_pos**cpml_m
          z_pos_ma = (1.0_num - z_pos)**cpml_ma

          cpml_kappa_bz(iz-1) = 1.0_num + (cpml_kappa_max - 1.0_num) * z_pos_m
          cpml_sigma_bz(iz-1) = cpml_sigma_maxval * z_pos_m
          cpml_a_bz(iz-1) = cpml_a_max * z_pos_ma
        END DO
      END IF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nz_global_min <= nz_global - cpml_thickness - fng + 2 &
          .AND. nz_global_max >= nz_global - cpml_thickness - fng + 2) THEN
        add_laser(i) = .TRUE.
        cpml_z_max_laser_idx = &
            nz_global - cpml_thickness - fng + 2 - nz_global_min
      END IF
    END IF

  END SUBROUTINE set_cpml_helpers



  SUBROUTINE allocate_cpml_fields

    ! I will ignore memory consumption issues and, for simplicity,
    ! allocate the boundary fields throughout the whole simulation box.
    ALLOCATE(cpml_psi_eyx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_ezx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_byx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_bzx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ALLOCATE(cpml_psi_exy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_ezy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_bxy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_bzy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ALLOCATE(cpml_psi_exz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_eyz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_bxz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(cpml_psi_byz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    cpml_psi_eyx = 0.0_num
    cpml_psi_ezx = 0.0_num
    cpml_psi_byx = 0.0_num
    cpml_psi_bzx = 0.0_num

    cpml_psi_exy = 0.0_num
    cpml_psi_ezy = 0.0_num
    cpml_psi_bxy = 0.0_num
    cpml_psi_bzy = 0.0_num

    cpml_psi_exz = 0.0_num
    cpml_psi_eyz = 0.0_num
    cpml_psi_bxz = 0.0_num
    cpml_psi_byz = 0.0_num

  END SUBROUTINE allocate_cpml_fields



  SUBROUTINE deallocate_cpml_helpers

    DEALLOCATE(cpml_kappa_ex, cpml_kappa_bx)
    DEALLOCATE(cpml_a_ex, cpml_a_bx)
    DEALLOCATE(cpml_sigma_ex, cpml_sigma_bx)

    DEALLOCATE(cpml_kappa_ey, cpml_kappa_by)
    DEALLOCATE(cpml_a_ey, cpml_a_by)
    DEALLOCATE(cpml_sigma_ey, cpml_sigma_by)

    DEALLOCATE(cpml_kappa_ez, cpml_kappa_bz)
    DEALLOCATE(cpml_a_ez, cpml_a_bz)
    DEALLOCATE(cpml_sigma_ez, cpml_sigma_bz)

  END SUBROUTINE deallocate_cpml_helpers



  SUBROUTINE cpml_advance_e_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix, iy, iz
    REAL(num) :: acoeff, bcoeff, ccoeff_d, fac
    REAL(num) :: kappa, sigma

    fac = tstep * c**2

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO iy = 1,ny
          DO ipos = cpml_x_min_start,cpml_x_min_end
            kappa = cpml_kappa_ex(ipos)
            sigma = cpml_sigma_ex(ipos)
            acoeff = cpml_a_ex(ipos)
            bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
            ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
                / (sigma + kappa * acoeff) / dx

            cpml_psi_eyx(ipos,iy,iz) = bcoeff * cpml_psi_eyx(ipos,iy,iz) &
                + ccoeff_d * (bz(ipos,iy,iz) - bz(ipos-1,iy,iz))
            ey(ipos,iy,iz) = ey(ipos,iy,iz) - fac * cpml_psi_eyx(ipos,iy,iz)

            cpml_psi_ezx(ipos,iy,iz) = bcoeff * cpml_psi_ezx(ipos,iy,iz) &
                + ccoeff_d * (by(ipos,iy,iz) - by(ipos-1,iy,iz))
            ez(ipos,iy,iz) = ez(ipos,iy,iz) + fac * cpml_psi_ezx(ipos,iy,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO iy = 1,ny
          DO ipos = cpml_x_max_start,cpml_x_max_end
            kappa = cpml_kappa_ex(ipos)
            sigma = cpml_sigma_ex(ipos)
            acoeff = cpml_a_ex(ipos)
            bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
            ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
                / (sigma + kappa * acoeff) / dx

            cpml_psi_eyx(ipos,iy,iz) = bcoeff * cpml_psi_eyx(ipos,iy,iz) &
                + ccoeff_d * (bz(ipos,iy,iz) - bz(ipos-1,iy,iz))
            ey(ipos,iy,iz) = ey(ipos,iy,iz) - fac * cpml_psi_eyx(ipos,iy,iz)

            cpml_psi_ezx(ipos,iy,iz) = bcoeff * cpml_psi_ezx(ipos,iy,iz) &
                + ccoeff_d * (by(ipos,iy,iz) - by(ipos-1,iy,iz))
            ez(ipos,iy,iz) = ez(ipos,iy,iz) + fac * cpml_psi_ezx(ipos,iy,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= y_min boundary =============

    IF (bc_field(c_bd_y_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_min) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO ipos = cpml_y_min_start,cpml_y_min_end
          kappa = cpml_kappa_ey(ipos)
          sigma = cpml_sigma_ey(ipos)
          acoeff = cpml_a_ey(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dy

          DO ix = 1,nx
            cpml_psi_exy(ix,ipos,iz) = bcoeff * cpml_psi_exy(ix,ipos,iz) &
                + ccoeff_d * (bz(ix,ipos,iz) - bz(ix,ipos-1,iz))
            ex(ix,ipos,iz) = ex(ix,ipos,iz) + fac * cpml_psi_exy(ix,ipos,iz)

            cpml_psi_ezy(ix,ipos,iz) = bcoeff * cpml_psi_ezy(ix,ipos,iz) &
                + ccoeff_d * (bx(ix,ipos,iz) - bx(ix,ipos-1,iz))
            ez(ix,ipos,iz) = ez(ix,ipos,iz) - fac * cpml_psi_ezy(ix,ipos,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= y_max boundary =============

    IF (bc_field(c_bd_y_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_max) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO ipos = cpml_y_max_start,cpml_y_max_end
          kappa = cpml_kappa_ey(ipos)
          sigma = cpml_sigma_ey(ipos)
          acoeff = cpml_a_ey(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dy

          DO ix = 1,nx
            cpml_psi_exy(ix,ipos,iz) = bcoeff * cpml_psi_exy(ix,ipos,iz) &
                + ccoeff_d * (bz(ix,ipos,iz) - bz(ix,ipos-1,iz))
            ex(ix,ipos,iz) = ex(ix,ipos,iz) + fac * cpml_psi_exy(ix,ipos,iz)

            cpml_psi_ezy(ix,ipos,iz) = bcoeff * cpml_psi_ezy(ix,ipos,iz) &
                + ccoeff_d * (bx(ix,ipos,iz) - bx(ix,ipos-1,iz))
            ez(ix,ipos,iz) = ez(ix,ipos,iz) - fac * cpml_psi_ezy(ix,ipos,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= z_min boundary =============

    IF (bc_field(c_bd_z_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_z_min) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_z_min_start,cpml_z_min_end
        kappa = cpml_kappa_ez(ipos)
        sigma = cpml_sigma_ez(ipos)
        acoeff = cpml_a_ez(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dz

        DO iy = 1,ny
          DO ix = 1,nx
            cpml_psi_exz(ix,iy,ipos) = bcoeff * cpml_psi_exz(ix,iy,ipos) &
                + ccoeff_d * (by(ix,iy,ipos) - by(ix,iy,ipos-1))
            ex(ix,iy,ipos) = ex(ix,iy,ipos) - fac * cpml_psi_exz(ix,iy,ipos)

            cpml_psi_eyz(ix,iy,ipos) = bcoeff * cpml_psi_eyz(ix,iy,ipos) &
                + ccoeff_d * (bx(ix,iy,ipos) - bx(ix,iy,ipos-1))
            ey(ix,iy,ipos) = ey(ix,iy,ipos) + fac * cpml_psi_eyz(ix,iy,ipos)
          END DO
        END DO
      END DO
    END IF

    ! ============= z_max boundary =============

    IF (bc_field(c_bd_z_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_z_max) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_z_max_start,cpml_z_max_end
        kappa = cpml_kappa_ez(ipos)
        sigma = cpml_sigma_ez(ipos)
        acoeff = cpml_a_ez(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dz

        DO iy = 1,ny
          DO ix = 1,nx
            cpml_psi_exz(ix,iy,ipos) = bcoeff * cpml_psi_exz(ix,iy,ipos) &
                + ccoeff_d * (by(ix,iy,ipos) - by(ix,iy,ipos-1))
            ex(ix,iy,ipos) = ex(ix,iy,ipos) - fac * cpml_psi_exz(ix,iy,ipos)

            cpml_psi_eyz(ix,iy,ipos) = bcoeff * cpml_psi_eyz(ix,iy,ipos) &
                + ccoeff_d * (bx(ix,iy,ipos) - bx(ix,iy,ipos-1))
            ey(ix,iy,ipos) = ey(ix,iy,ipos) + fac * cpml_psi_eyz(ix,iy,ipos)
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE cpml_advance_e_currents



  SUBROUTINE cpml_advance_b_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix, iy, iz
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO iy = 1,ny
          DO ipos = cpml_x_min_start,cpml_x_min_end
            kappa = cpml_kappa_bx(ipos)
            sigma = cpml_sigma_bx(ipos)
            acoeff = cpml_a_bx(ipos)
            bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
            ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
                / (sigma + kappa * acoeff) / dx

            cpml_psi_byx(ipos,iy,iz) = bcoeff * cpml_psi_byx(ipos,iy,iz) &
                + ccoeff_d * (ez(ipos+1,iy,iz) - ez(ipos,iy,iz))
            by(ipos,iy,iz) = by(ipos,iy,iz) + tstep * cpml_psi_byx(ipos,iy,iz)

            cpml_psi_bzx(ipos,iy,iz) = bcoeff * cpml_psi_bzx(ipos,iy,iz) &
                + ccoeff_d * (ey(ipos+1,iy,iz) - ey(ipos,iy,iz))
            bz(ipos,iy,iz) = bz(ipos,iy,iz) - tstep * cpml_psi_bzx(ipos,iy,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO iy = 1,ny
          DO ipos = cpml_x_max_start-1,cpml_x_max_end-1
            kappa = cpml_kappa_bx(ipos)
            sigma = cpml_sigma_bx(ipos)
            acoeff = cpml_a_bx(ipos)
            bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
            ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
                / (sigma + kappa * acoeff) / dx

            cpml_psi_byx(ipos,iy,iz) = bcoeff * cpml_psi_byx(ipos,iy,iz) &
                + ccoeff_d * (ez(ipos+1,iy,iz) - ez(ipos,iy,iz))
            by(ipos,iy,iz) = by(ipos,iy,iz) + tstep * cpml_psi_byx(ipos,iy,iz)

            cpml_psi_bzx(ipos,iy,iz) = bcoeff * cpml_psi_bzx(ipos,iy,iz) &
                + ccoeff_d * (ey(ipos+1,iy,iz) - ey(ipos,iy,iz))
            bz(ipos,iy,iz) = bz(ipos,iy,iz) - tstep * cpml_psi_bzx(ipos,iy,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= y_min boundary =============

    IF (bc_field(c_bd_y_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_min) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO ipos = cpml_y_min_start,cpml_y_min_end
          kappa = cpml_kappa_by(ipos)
          sigma = cpml_sigma_by(ipos)
          acoeff = cpml_a_by(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dy

          DO ix = 1,nx
            cpml_psi_bxy(ix,ipos,iz) = bcoeff * cpml_psi_bxy(ix,ipos,iz) &
                + ccoeff_d * (ez(ix,ipos+1,iz) - ez(ix,ipos,iz))
            bx(ix,ipos,iz) = bx(ix,ipos,iz) - tstep * cpml_psi_bxy(ix,ipos,iz)

            cpml_psi_bzy(ix,ipos,iz) = bcoeff * cpml_psi_bzy(ix,ipos,iz) &
                + ccoeff_d * (ex(ix,ipos+1,iz) - ex(ix,ipos,iz))
            bz(ix,ipos,iz) = bz(ix,ipos,iz) + tstep * cpml_psi_bzy(ix,ipos,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= y_max boundary =============

    IF (bc_field(c_bd_y_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_max) == c_bc_cpml_outflow) THEN
      DO iz = 1,nz
        DO ipos = cpml_y_max_start-1,cpml_y_max_end-1
          kappa = cpml_kappa_by(ipos)
          sigma = cpml_sigma_by(ipos)
          acoeff = cpml_a_by(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dy

          DO ix = 1,nx
            cpml_psi_bxy(ix,ipos,iz) = bcoeff * cpml_psi_bxy(ix,ipos,iz) &
                + ccoeff_d * (ez(ix,ipos+1,iz) - ez(ix,ipos,iz))
            bx(ix,ipos,iz) = bx(ix,ipos,iz) - tstep * cpml_psi_bxy(ix,ipos,iz)

            cpml_psi_bzy(ix,ipos,iz) = bcoeff * cpml_psi_bzy(ix,ipos,iz) &
                + ccoeff_d * (ex(ix,ipos+1,iz) - ex(ix,ipos,iz))
            bz(ix,ipos,iz) = bz(ix,ipos,iz) + tstep * cpml_psi_bzy(ix,ipos,iz)
          END DO
        END DO
      END DO
    END IF

    ! ============= z_min boundary =============

    IF (bc_field(c_bd_z_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_z_min) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_z_min_start,cpml_z_min_end
        kappa = cpml_kappa_bz(ipos)
        sigma = cpml_sigma_bz(ipos)
        acoeff = cpml_a_bz(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dz

        DO iy = 1,ny
          DO ix = 1,nx
            cpml_psi_bxz(ix,iy,ipos) = bcoeff * cpml_psi_bxz(ix,iy,ipos) &
                + ccoeff_d * (ey(ix,iy,ipos+1) - ey(ix,iy,ipos))
            bx(ix,iy,ipos) = bx(ix,iy,ipos) + tstep * cpml_psi_bxz(ix,iy,ipos)

            cpml_psi_byz(ix,iy,ipos) = bcoeff * cpml_psi_byz(ix,iy,ipos) &
                + ccoeff_d * (ex(ix,iy,ipos+1) - ex(ix,iy,ipos))
            by(ix,iy,ipos) = by(ix,iy,ipos) - tstep * cpml_psi_byz(ix,iy,ipos)
          END DO
        END DO
      END DO
    END IF

    ! ============= z_max boundary =============

    IF (bc_field(c_bd_z_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_z_max) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_z_max_start-1,cpml_z_max_end-1
        kappa = cpml_kappa_bz(ipos)
        sigma = cpml_sigma_bz(ipos)
        acoeff = cpml_a_bz(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dz

        DO iy = 1,ny
          DO ix = 1,nx
            cpml_psi_bxz(ix,iy,ipos) = bcoeff * cpml_psi_bxz(ix,iy,ipos) &
                + ccoeff_d * (ey(ix,iy,ipos+1) - ey(ix,iy,ipos))
            bx(ix,iy,ipos) = bx(ix,iy,ipos) + tstep * cpml_psi_bxz(ix,iy,ipos)

            cpml_psi_byz(ix,iy,ipos) = bcoeff * cpml_psi_byz(ix,iy,ipos) &
                + ccoeff_d * (ex(ix,iy,ipos+1) - ex(ix,iy,ipos))
            by(ix,iy,ipos) = by(ix,iy,ipos) - tstep * cpml_psi_byz(ix,iy,ipos)
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE cpml_advance_b_currents

END MODULE boundary
