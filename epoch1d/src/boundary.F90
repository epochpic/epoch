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
        boundary = (/ 'x_min', 'x_max' /)
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
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, ng, nx)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths(field, ng, nx_local)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local
    INTEGER :: basetype, i, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    ALLOCATE(temp(ng))

    CALL MPI_SENDRECV(field(1), ng, basetype, proc_x_min, &
        tag, temp, ng, basetype, proc_x_max, tag, comm, status, errcode)

    IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
      n = 1
      DO i = nx_local+1, nx_local+ng
        field(i) = temp(n)
        n = n + 1
      END DO
    END IF

    CALL MPI_SENDRECV(field(nx_local+1-ng), ng, basetype, proc_x_max, &
        tag, temp, ng, basetype, proc_x_min, tag, comm, status, errcode)

    IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
      n = 1
      DO i = 1-ng, 0
        field(i) = temp(n)
        n = n + 1
      END DO
    END IF

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE do_field_mpi_with_lengths_r4(field, ng, nx_local)

    INTEGER, INTENT(IN) :: ng
    REAL(r4), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local
    INTEGER :: basetype, i, n
    REAL(r4), ALLOCATABLE :: temp(:)

    basetype = MPI_REAL4

    ALLOCATE(temp(ng))

    CALL MPI_SENDRECV(field(1), ng, basetype, proc_x_min, &
        tag, temp, ng, basetype, proc_x_max, tag, comm, status, errcode)

    IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max) == c_bc_periodic) THEN
      n = 1
      DO i = nx_local+1, nx_local+ng
        field(i) = temp(n)
        n = n + 1
      END DO
    END IF

    CALL MPI_SENDRECV(field(nx_local+1-ng), ng, basetype, proc_x_max, &
        tag, temp, ng, basetype, proc_x_min, tag, comm, status, errcode)

    IF (.NOT. x_min_boundary .OR. bc_field(c_bd_x_min) == c_bc_periodic) THEN
      n = 1
      DO i = 1-ng, 0
        field(i) = temp(n)
        n = n + 1
      END DO
    END IF

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths_r4



  SUBROUTINE field_zero_gradient(field, stagger_type, boundary)

    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(i-ng) = field(ng-i)
        END DO
      ELSE
        DO i = 1, ng
          field(i-ng) = field(ng+1-i)
        END DO
      END IF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(nn+i) = field(nn-i)
        END DO
      ELSE
        DO i = 1, ng
          field(nn+i) = field(nn+1-i)
        END DO
      END IF
    END IF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, ng, stagger_type, boundary)

    INTEGER, INTENT(IN) :: ng, stagger_type, boundary
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(i-ng) = -field(ng-i)
        END DO
        field(0) = 0.0_num
      ELSE
        DO i = 1, ng
          field(i-ng) = -field(ng+1-i)
        END DO
      END IF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nn) = 0.0_num
        DO i = 1, ng-1
          field(nn+i) = -field(nn-i)
        END DO
      ELSE
        DO i = 1, ng
          field(nn+i) = -field(nn+1-i)
        END DO
      END IF
    END IF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE particle_reflection_bcs(array, ng, flip_direction, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: array
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
          array(i) = array(i) - array(-i)
          array(-i) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng-1
          array(i) = array(i) + array(1-i)
          array(1-i) = 0.0_num
        END DO
      END IF
    END IF

    n = n + 1
    bc = bc_species(n)
    IF (x_max_boundary .AND. bc == c_bc_reflect) THEN
      IF (flip_dir == (n-1)/2 + 1) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(nn-i) = array(nn-i) - array(nn+i)
          array(nn+i) = 0.0_num
        END DO
      ELSE
        DO i = 1, ng
          array(nn+1-i) = array(nn+1-i) + array(nn+i)
          array(nn+i) = 0.0_num
        END DO
      END IF
    END IF

  END SUBROUTINE particle_reflection_bcs



  SUBROUTINE particle_periodic_bcs(array, ng, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: species
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes
    INTEGER :: n, nn, i
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
    n = 0
    nn = sizes(n/2+1) - 2 * ng

    ALLOCATE(temp(ng))

    ! Don't bother communicating non-periodic boundaries
    neighbour_local = neighbour(:)
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
    CALL MPI_SENDRECV(array(nn+1), ng, mpireal, &
        neighbour_local( 1), tag, temp, ng, mpireal, &
        neighbour_local(-1), tag, comm, status, errcode)

    n = n + 1
    array(1:ng) = array(1:ng) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng), ng, mpireal, &
        neighbour_local(-1), tag, temp, ng, mpireal, &
        neighbour_local( 1), tag, comm, status, errcode)

    n = n + 1
    array(nn+1-ng:nn) = array(nn+1-ng:nn) + temp

    DEALLOCATE(temp)

    IF (PRESENT(species)) CALL particle_clear_bcs(array, ng)

  END SUBROUTINE particle_periodic_bcs



  SUBROUTINE particle_clear_bcs(array, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: array
    INTEGER, DIMENSION(c_ndims) :: sizes
    INTEGER :: n, nn

    sizes = SHAPE(array)
    n = 0

    nn = sizes(n/2+1) - 2 * ng

    n = n + 1
    array(:0) = 0.0_num
    n = n + 1
    array(nn+1:) = 0.0_num

  END SUBROUTINE particle_clear_bcs



  SUBROUTINE processor_summation_bcs(array, ng, flip_direction, species)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: array
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

    CALL bfield_bcs(.TRUE.)

  END SUBROUTINE bfield_final_bcs



  SUBROUTINE setup_bc_lists

    INTEGER(i8) :: ispecies, ipart
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    REAL(num) :: bnd_x_min, bnd_x_max
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

      ALLOCATE(species_list(ispecies)%boundary_particles)
      NULLIFY(species_list(ispecies)%boundary_particles%particle)
      NULLIFY(species_list(ispecies)%boundary_particles%next)
      NULLIFY(bnd_part_next)
      bnd_part_last => species_list(ispecies)%boundary_particles

      DO ipart = 1, species_list(ispecies)%attached_list%count
        ! Move particle to boundary candidate list
        IF (current%part_pos < bnd_x_min &
            .OR. current%part_pos > bnd_x_max) THEN
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
    TYPE(particle_list), DIMENSION(-1:1) :: send, recv
    INTEGER :: xbd
    INTEGER(i8) :: ixp
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    LOGICAL :: out_of_bounds
    INTEGER :: sgn, bc, ispecies, i, ix
    REAL(num) :: temp(3)
    REAL(num) :: part_pos
    REAL(num) :: x_shift

    x_shift = length_x + 2.0_num * dx * REAL(cpml_thickness, num)

    DO ispecies = 1, n_species
      bnd_part => species_list(ispecies)%boundary_particles
      NULLIFY(bnd_part_last)

      bc_species = species_list(ispecies)%bc_particle

      DO ix = -1, 1, 2
        CALL create_empty_partlist(send(ix))
        CALL create_empty_partlist(recv(ix))
      END DO

      DO WHILE (ASSOCIATED(bnd_part))
        bnd_part_last => bnd_part

        cur => bnd_part%particle
        bnd_part => bnd_part%next
        DEALLOCATE(bnd_part_last)

        xbd = 0
        out_of_bounds = .FALSE.

        part_pos = cur%part_pos
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
                cur%part_pos = 2.0_num * x_min - part_pos
                cur%part_p(1) = -cur%part_p(1)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (x_min_boundary) THEN
                cur%part_pos = part_pos - sgn * x_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos < x_min_outer) THEN
                xbd = 0
                DO i = 1, 3
                  temp(i) = species_list(ispecies)%ext_temp_x_min(i)
                END DO

                CALL id_registry%delete_all(cur)
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
                cur%id = generate_id()
#endif

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

                cur%part_pos = 2.0_num * x_min_outer - part_pos

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
                cur%part_pos = 2.0_num * x_max - part_pos
                cur%part_p(1) = -cur%part_p(1)
              END IF
            ELSE IF (bc == c_bc_periodic) THEN
              IF (x_max_boundary) THEN
                cur%part_pos = part_pos - sgn * x_shift
              END IF
            ELSE IF (bc == c_bc_thermal) THEN
              IF (part_pos >= x_max_outer) THEN
                xbd = 0
                DO i = 1, 3
                  temp(i) = species_list(ispecies)%ext_temp_x_max(i)
                END DO

                CALL id_registry%delete_all(cur)
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
                cur%id = generate_id()
#endif

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

                cur%part_pos = 2.0_num * x_max_outer - part_pos

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
        ELSE IF (ABS(xbd) > 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd), cur)
        END IF
      END DO

      ! swap Particles
      DO ix = -1, 1, 2
        ixp = -ix
        CALL partlist_sendrecv(send(ix), recv(ixp), &
            neighbour(ix), neighbour(ixp))
        CALL append_partlist(species_list(ispecies)%attached_list, &
            recv(ixp))
      END DO

      DO ix = -1, 1, 2
        CALL destroy_partlist(send(ix))
        CALL destroy_partlist(recv(ix))
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



  SUBROUTINE set_cpml_helpers(nx, nx_global_min, nx_global_max)

    INTEGER, INTENT(IN) :: nx, nx_global_min, nx_global_max
    INTEGER :: i
    INTEGER :: ix, ix_glob
    INTEGER, PARAMETER :: cpml_m = 3
    INTEGER, PARAMETER :: cpml_ma = 1
    REAL(num) :: x_pos, x_pos_m, x_pos_ma
    REAL(num) :: cpml_sigma_maxval

    ALLOCATE(cpml_kappa_ex(1-ng:nx+ng), cpml_kappa_bx(1-ng:nx+ng))
    ALLOCATE(cpml_a_ex(1-ng:nx+ng), cpml_a_bx(1-ng:nx+ng))
    ALLOCATE(cpml_sigma_ex(1-ng:nx+ng), cpml_sigma_bx(1-ng:nx+ng))

    cpml_kappa_ex = 1.0_num
    cpml_kappa_bx = 1.0_num

    cpml_a_ex = 0.0_num
    cpml_sigma_ex = 0.0_num
    cpml_a_bx = 0.0_num
    cpml_sigma_bx = 0.0_num

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

  END SUBROUTINE set_cpml_helpers



  SUBROUTINE allocate_cpml_fields

    ! I will ignore memory consumption issues and, for simplicity,
    ! allocate the boundary fields throughout the whole simulation box.
    ALLOCATE(cpml_psi_eyx(1-ng:nx+ng))
    ALLOCATE(cpml_psi_ezx(1-ng:nx+ng))
    ALLOCATE(cpml_psi_byx(1-ng:nx+ng))
    ALLOCATE(cpml_psi_bzx(1-ng:nx+ng))

    cpml_psi_eyx = 0.0_num
    cpml_psi_ezx = 0.0_num
    cpml_psi_byx = 0.0_num
    cpml_psi_bzx = 0.0_num

  END SUBROUTINE allocate_cpml_fields



  SUBROUTINE deallocate_cpml_helpers

    DEALLOCATE(cpml_kappa_ex, cpml_kappa_bx)
    DEALLOCATE(cpml_a_ex, cpml_a_bx)
    DEALLOCATE(cpml_sigma_ex, cpml_sigma_bx)

  END SUBROUTINE deallocate_cpml_helpers



  SUBROUTINE cpml_advance_e_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos
    REAL(num) :: acoeff, bcoeff, ccoeff_d, fac
    REAL(num) :: kappa, sigma

    fac = tstep * c**2

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_min_start,cpml_x_min_end
        kappa = cpml_kappa_ex(ipos)
        sigma = cpml_sigma_ex(ipos)
        acoeff = cpml_a_ex(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_eyx(ipos) = bcoeff * cpml_psi_eyx(ipos) &
            + ccoeff_d * (bz(ipos) - bz(ipos-1))
        ey(ipos) = ey(ipos) - fac * cpml_psi_eyx(ipos)

        cpml_psi_ezx(ipos) = bcoeff * cpml_psi_ezx(ipos) &
            + ccoeff_d * (by(ipos) - by(ipos-1))
        ez(ipos) = ez(ipos) + fac * cpml_psi_ezx(ipos)
      END DO
    END IF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_max_start,cpml_x_max_end
        kappa = cpml_kappa_ex(ipos)
        sigma = cpml_sigma_ex(ipos)
        acoeff = cpml_a_ex(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_eyx(ipos) = bcoeff * cpml_psi_eyx(ipos) &
            + ccoeff_d * (bz(ipos) - bz(ipos-1))
        ey(ipos) = ey(ipos) - fac * cpml_psi_eyx(ipos)

        cpml_psi_ezx(ipos) = bcoeff * cpml_psi_ezx(ipos) &
            + ccoeff_d * (by(ipos) - by(ipos-1))
        ez(ipos) = ez(ipos) + fac * cpml_psi_ezx(ipos)
      END DO
    END IF

  END SUBROUTINE cpml_advance_e_currents



  SUBROUTINE cpml_advance_b_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_min_start,cpml_x_min_end
        kappa = cpml_kappa_bx(ipos)
        sigma = cpml_sigma_bx(ipos)
        acoeff = cpml_a_bx(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_byx(ipos) = bcoeff * cpml_psi_byx(ipos) &
            + ccoeff_d * (ez(ipos+1) - ez(ipos))
        by(ipos) = by(ipos) + tstep * cpml_psi_byx(ipos)

        cpml_psi_bzx(ipos) = bcoeff * cpml_psi_bzx(ipos) &
            + ccoeff_d * (ey(ipos+1) - ey(ipos))
        bz(ipos) = bz(ipos) - tstep * cpml_psi_bzx(ipos)
      END DO
    END IF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) == c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_max_start-1,cpml_x_max_end-1
        kappa = cpml_kappa_bx(ipos)
        sigma = cpml_sigma_bx(ipos)
        acoeff = cpml_a_bx(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_byx(ipos) = bcoeff * cpml_psi_byx(ipos) &
            + ccoeff_d * (ez(ipos+1) - ez(ipos))
        by(ipos) = by(ipos) + tstep * cpml_psi_byx(ipos)

        cpml_psi_bzx(ipos) = bcoeff * cpml_psi_bzx(ipos) &
            + ccoeff_d * (ey(ipos+1) - ey(ipos))
        bz(ipos) = bz(ipos) - tstep * cpml_psi_bzx(ipos)
      END DO
    END IF

  END SUBROUTINE cpml_advance_b_currents

END MODULE boundary
