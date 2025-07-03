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

MODULE balance

  USE boundary
  USE mpi_subtype_control
  USE redblack_module
  USE timer
  USE utilities

  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_x_min, new_cell_x_max
  LOGICAL :: overriding
  REAL(num) :: load_av
  INTEGER :: old_comm, old_coordinates(c_ndims)

CONTAINS

  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    REAL(num) :: balance_frac, balance_frac_final, balance_improvement
    REAL(num) :: load_local, load_sum, load_max
    INTEGER(i8) :: npart_local
    INTEGER, SAVE :: balance_check_frequency = 1
    INTEGER, SAVE :: last_check = (1 - HUGE(1)) / 2
    INTEGER, SAVE :: last_full_check = (1 - HUGE(1)) / 2
    LOGICAL, SAVE :: first_flag = .TRUE.
    LOGICAL :: first_message, restarting, full_check, attempt_balance
    LOGICAL :: use_redistribute_domain, use_redistribute_particles
#ifdef PARTICLE_DEBUG
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
#endif

    ! On one processor do nothing to save time
    IF (nproc == 1) RETURN

    full_check = over_ride
    IF (step - last_full_check < dlb_force_interval) THEN
      IF (step - last_check < balance_check_frequency) RETURN
    ELSE
      full_check = .TRUE.
    END IF

    restarting = .FALSE.
    use_redistribute_domain = .FALSE.
    use_redistribute_particles = .FALSE.
    attempt_balance = use_balance

    IF (first_flag) THEN
      attempt_balance = balance_first
      IF (use_exact_restart) attempt_balance = .FALSE.
      first_flag = .FALSE.
      first_message = .TRUE.
      IF (ic_from_restart) THEN
        restarting = .TRUE.
        use_redistribute_particles = .TRUE.
      END IF
    ELSE
      first_message = .FALSE.
    END IF
    last_check = step

    ! count particles
    npart_local = get_total_local_particles()
    load_local = REAL(push_per_field * npart_local + nx, num)

    CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, comm, errcode)
    CALL MPI_ALLREDUCE(load_local, load_sum, 1, mpireal, MPI_SUM, comm, errcode)

    load_av = load_sum / nproc

    balance_frac = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))

    ! The over_ride flag allows the code to force a load balancing sweep
    ! at t = 0
    IF (.NOT. full_check .AND. balance_frac > dlb_threshold) THEN
      balance_check_frequency = &
          MIN(balance_check_frequency * 2, dlb_maximum_interval)
      IF (rank == 0) THEN
        PRINT'(''Skipping redistribution. Balance:'', F6.3, &
              &'', threshold:'', F6.3, '', next: '', i9)', &
              balance_frac, dlb_threshold, &
              MIN(step + balance_check_frequency, &
                  last_full_check + dlb_force_interval)
      END IF
      RETURN
    END IF

    IF (timer_collect) CALL timer_start(c_timer_balance)

    IF (attempt_balance) THEN
      overriding = full_check

      ALLOCATE(load_x(nx_global + 2 * ng))
      CALL get_load(load_x, part_load_func)

      ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))

      CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)

      DEALLOCATE(load_x)

      IF (.NOT.restarting) THEN
        CALL create_npart_per_cell_array

        CALL calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                          new_cell_x_min, new_cell_x_max)

        DEALLOCATE(npart_per_cell_array)

        balance_improvement = (balance_frac_final - balance_frac) / balance_frac

        last_full_check = step

        ! Consider load balancing a success if the load imbalance improved by
        ! more than 5 percent
        IF (balance_improvement > 0.05_num) THEN
          use_redistribute_domain = .TRUE.
          use_redistribute_particles = .TRUE.
          first_message = .FALSE.
          balance_check_frequency = 1
          IF (rank == 0) THEN
            IF (use_balance) THEN
              PRINT'(''Redistributing.          Balance:'', F6.3, &
                    &'',     after:'', F6.3, '', next: '', i9)', &
                    balance_frac, balance_frac_final, &
                    MIN(step + balance_check_frequency, &
                        last_full_check + dlb_force_interval)
            ELSE
              PRINT'(''Redistributing.          Balance:'', F6.3, &
                    &'',     after:'', F6.3, ''  (initial setup)'')', &
                    balance_frac, balance_frac_final
            END IF
          END IF
        ELSE
          IF (.NOT.first_message) THEN
            balance_check_frequency = &
                MIN(balance_check_frequency * 2, dlb_maximum_interval)
            IF (rank == 0) THEN
              IF (use_balance) THEN
                PRINT'(''Skipping redistribution. Balance:'', F6.3, &
                      &'',     after:'', F6.3, '', next: '', i9)', &
                      balance_frac, balance_frac_final, &
                      MIN(step + balance_check_frequency, &
                          last_full_check + dlb_force_interval)
              ELSE
                PRINT'(''Skipping redistribution. Balance:'', F6.3, &
                      &'',     after:'', F6.3, ''  (initial setup)'')', &
                      balance_frac, balance_frac_final
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF

    IF (use_redistribute_domain) THEN
      old_comm = comm
      old_coordinates(:) = coordinates(:)
      CALL redistribute_domain
    END IF

    IF (ALLOCATED(new_cell_x_min)) THEN
      DEALLOCATE(new_cell_x_min, new_cell_x_max)
    END IF

    ! Redistribute the particles onto their new processors
    IF (use_redistribute_particles) CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (full_check) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%processor_at_t0 = rank
          current => current%next
        END DO
      END DO
    END IF
#endif

    IF (first_message) THEN
      npart_local = get_total_local_particles()
      load_local = REAL(push_per_field * npart_local + nx, num)

      CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, &
          comm, errcode)

      balance_frac_final = (load_av + SQRT(load_av)) &
          / (load_max + SQRT(load_max))
      balance_check_frequency = 1

      IF (rank == 0) THEN
        IF (use_redistribute_domain .OR. use_redistribute_particles) THEN
          PRINT'(''Redistributing.          Balance:'', F6.3, &
                & 18X, '' (initial setup)'')', &
                balance_frac_final
        ELSE
          PRINT'(''Skipping redistribution. Balance:'', F6.3, &
                & 18X, '' (initial setup)'')', &
                balance_frac_final
        END IF
      END IF
    END IF

    use_exact_restart = .FALSE.

    IF (timer_collect) CALL timer_stop(c_timer_balance)

  END SUBROUTINE balance_workload



  SUBROUTINE pre_balance_workload(old_communicator, old_coords)

    INTEGER, INTENT(IN), OPTIONAL :: old_communicator, old_coords(:)
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    REAL(num) :: balance_frac, balance_frac_final, balance_improvement
    LOGICAL :: use_redistribute_domain

    ! On one processor do nothing to save time
    IF (nproc == 1 .OR. .NOT.use_pre_balance) RETURN

    overriding = .TRUE.

    ALLOCATE(load_x(nx_global + 2 * ng))

    CALL get_load(load_x, array_load_func)

    ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))

    CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)

    DEALLOCATE(load_x)

    CALL calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                      new_cell_x_min, new_cell_x_max, .TRUE.)

    IF (ALLOCATED(npart_per_cell_array)) DEALLOCATE(npart_per_cell_array)

    balance_improvement = (balance_frac_final - balance_frac) / balance_frac

    ! Consider load balancing a success if the load imbalance improved by
    ! more than 5 percent
    IF (balance_improvement > 0.05_num) THEN
      use_redistribute_domain = .TRUE.
      IF (rank == 0) THEN
        PRINT'(''Redistributing.          Balance:'', F6.3, &
              &'',     after:'', F6.3, '' (pre-load balance)'')', &
              balance_frac, balance_frac_final
      END IF
    ELSE
      use_redistribute_domain = .FALSE.
      IF (rank == 0) THEN
        PRINT'(''Skipping redistribution. Balance:'', F6.3, &
              &'',     after:'', F6.3, '' (pre-load balance)'')', &
              balance_frac, balance_frac_final
      END IF
    END IF

    IF (PRESENT(old_communicator)) use_redistribute_domain = .TRUE.

    IF (use_redistribute_domain) THEN
      IF (PRESENT(old_communicator)) THEN
        old_comm = old_communicator
        old_coordinates(:) = old_coords(:)
        DEALLOCATE(x_grid_mins, x_grid_maxs)
        ALLOCATE(x_grid_mins(0:nprocx-1))
        ALLOCATE(x_grid_maxs(0:nprocx-1))
      ELSE
        old_comm = comm
        old_coordinates(:) = coordinates(:)
      END IF
      CALL redistribute_domain
    END IF

    IF (ALLOCATED(new_cell_x_min)) THEN
      DEALLOCATE(new_cell_x_min, new_cell_x_max)
    END IF

  END SUBROUTINE pre_balance_workload



  SUBROUTINE redistribute_domain

    INTEGER, DIMENSION(c_ndims,2) :: domain

    IF (.NOT.ALLOCATED(new_cell_x_min)) RETURN

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor

    domain(1,:) = (/new_cell_x_min(x_coords+1), new_cell_x_max(x_coords+1)/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    DEALLOCATE(cell_x_min, cell_x_max)
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))

    ! Copy the new lengths into the permanent variables
    cell_x_min(:) = new_cell_x_min(:)
    cell_x_max(:) = new_cell_x_max(:)

    ! Set the new nx
    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)
    n_global_min(1) = nx_global_min
    n_global_max(1) = nx_global_max

    nx = nx_global_max - nx_global_min + 1

    ! Do X array separately because we already have global copies
    DEALLOCATE(x)
    ALLOCATE(x(1-ng:nx+ng))
    DEALLOCATE(xb)
    ALLOCATE(xb(1-ng:nx+ng))

    CALL setup_grid_x

  END SUBROUTINE redistribute_domain



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the field variables over the new
    ! processor layout. If using a field of your own then set the
    ! redistribute_field subroutine to implement it.

    INTEGER :: nx_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp_sum
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4temp_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp, temp2
    TYPE(particle_species_migration), POINTER :: mg
    TYPE(particle_species), POINTER :: sp
    TYPE(initial_condition_block), POINTER :: ic
    INTEGER :: i, ispecies, io, id, nspec_local, mask

    nx_new = new_domain(1,2) - new_domain(1,1) + 1

    ! The following code is quite messy and repetitive. Unfortunately, the
    ! F90 standard does not allow the ALLOCATABLE attribute for subroutine
    ! arguments and POINTER arrays are not as fast.

    ! Full domain arrays

    ALLOCATE(temp(1-ng:nx_new+ng))

    ! Current will be recalculated during the particle push, so there
    ! is no need to copy the contents of the old arrays.
    ! If overriding, then we may not be doing a particle push next
    ! so we still have to balance the arrays.
    ! It is done slightly differently since the arrays may be
    ! a different size.

    IF (overriding) THEN
      ALLOCATE(temp2(1-ng:nx+ng))

      temp2 = jx
      CALL remap_field(temp2, temp)
      DEALLOCATE(jx)
      ALLOCATE(jx(1-jng:nx_new+jng))
      jx(0:nx_new+1) = temp(0:nx_new+1)

      temp2 = jy
      CALL remap_field(temp2, temp)
      DEALLOCATE(jy)
      ALLOCATE(jy(1-jng:nx_new+jng))
      jy(0:nx_new+1) = temp(0:nx_new+1)

      temp2 = jz
      CALL remap_field(temp2, temp)
      DEALLOCATE(jz)
      ALLOCATE(jz(1-jng:nx_new+jng))
      jz(0:nx_new+1) = temp(0:nx_new+1)

      DEALLOCATE(temp2)
    ELSE
      DEALLOCATE(jx)
      DEALLOCATE(jy)
      DEALLOCATE(jz)
      ALLOCATE(jx(1-jng:nx_new+jng))
      ALLOCATE(jy(1-jng:nx_new+jng))
      ALLOCATE(jz(1-jng:nx_new+jng))
    END IF

    CALL remap_field(ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(1-ng:nx_new+ng))
    ex = temp

    CALL remap_field(ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(1-ng:nx_new+ng))
    ey = temp

    CALL remap_field(ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(1-ng:nx_new+ng))
    ez = temp

    CALL remap_field(bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(1-ng:nx_new+ng))
    bx = temp

    CALL remap_field(by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(1-ng:nx_new+ng))
    by = temp

    CALL remap_field(bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(1-ng:nx_new+ng))
    bz = temp

    IF (pre_loading) THEN
      IF (ALLOCATED(global_species_density)) DEALLOCATE(global_species_density)
      IF (ALLOCATED(global_species_temp)) DEALLOCATE(global_species_temp)
      IF (ALLOCATED(global_species_drift)) DEALLOCATE(global_species_drift)
    END IF

    DO ispecies = 1, n_species
      sp => species_list(ispecies)
      mg => sp%migrate

      IF (mg%fluid) THEN
        CALL remap_field(mg%fluid_energy, temp)
        DEALLOCATE(mg%fluid_energy)
        ALLOCATE(mg%fluid_energy(1-ng:nx_new+ng))
        mg%fluid_energy = temp

        CALL remap_field(mg%fluid_density, temp)
        DEALLOCATE(mg%fluid_density)
        ALLOCATE(mg%fluid_density(1-ng:nx_new+ng))
        mg%fluid_density = temp
      END IF

      IF (sp%background_species) THEN
        CALL remap_field(sp%background_density, temp)
        DEALLOCATE(sp%background_density)
        ALLOCATE(sp%background_density(1-ng:nx_new+ng))
        sp%background_density = temp
      END IF

      IF (.NOT.pre_loading) CYCLE

      ! When load-balancing before the particles have been loaded, we must
      ! also redistribute the per-species initial conditions arrays.
      ! These are discarded after the initial setup

      ic => species_list(ispecies)%initial_conditions

      IF (ALLOCATED(ic%density)) THEN
        CALL remap_field(ic%density, temp)
        DEALLOCATE(ic%density)
        ALLOCATE(ic%density(1-ng:nx_new+ng))
        ic%density = temp
      END IF

      IF (ALLOCATED(ic%temp)) THEN
        IF (.NOT. ALLOCATED(temp_sum)) &
            ALLOCATE(temp_sum(1-ng:nx_new+ng,3))

        CALL remap_field(ic%temp(:,1), temp_sum(:,1))
        CALL remap_field(ic%temp(:,2), temp_sum(:,2))
        CALL remap_field(ic%temp(:,3), temp_sum(:,3))

        DEALLOCATE(ic%temp)
        ALLOCATE(ic%temp(1-ng:nx_new+ng,3))
        ic%temp = temp_sum
      END IF

      IF (ALLOCATED(ic%drift)) THEN
        IF (.NOT. ALLOCATED(temp_sum)) &
            ALLOCATE(temp_sum(1-ng:nx_new+ng,3))
        CALL remap_field(ic%drift(:,1), temp_sum(:,1))
        CALL remap_field(ic%drift(:,2), temp_sum(:,2))
        CALL remap_field(ic%drift(:,3), temp_sum(:,3))

        DEALLOCATE(ic%drift)
        ALLOCATE(ic%drift(1-ng:nx_new+ng,3))
        ic%drift = temp_sum
      END IF
    END DO

    IF (ALLOCATED(temp_sum)) DEALLOCATE(temp_sum)

    IF (cpml_boundaries) THEN
      CALL remap_field(cpml_psi_eyx, temp)
      DEALLOCATE(cpml_psi_eyx)
      ALLOCATE(cpml_psi_eyx(1-ng:nx_new+ng))
      cpml_psi_eyx = temp

      CALL remap_field(cpml_psi_byx, temp)
      DEALLOCATE(cpml_psi_byx)
      ALLOCATE(cpml_psi_byx(1-ng:nx_new+ng))
      cpml_psi_byx = temp

      CALL remap_field(cpml_psi_ezx, temp)
      DEALLOCATE(cpml_psi_ezx)
      ALLOCATE(cpml_psi_ezx(1-ng:nx_new+ng))
      cpml_psi_ezx = temp

      CALL remap_field(cpml_psi_bzx, temp)
      DEALLOCATE(cpml_psi_bzx)
      ALLOCATE(cpml_psi_bzx(1-ng:nx_new+ng))
      cpml_psi_bzx = temp

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers(nx_new, new_domain(1,1), new_domain(1,2))
    END IF

    DEALLOCATE(temp)

    ! Full domain arrays with an additional index

    DO id = 1, num_vars_to_dump
      io = averaged_var_block(id)
      IF (io == 0) CYCLE

      mask = io_block_list(io)%dumpmask(id)
      nspec_local = 0
      IF (IAND(mask, c_io_no_sum) == 0) &
          nspec_local = 1
      IF (IAND(mask, c_io_species) /= 0) &
          nspec_local = nspec_local + n_species

      IF (nspec_local <= 0) CYCLE
      nspec_local = nspec_local * averaged_var_dims(id)

      IF (io_block_list(io)%averaged_data(id)%dump_single) THEN
        IF (.NOT. ALLOCATED(io_block_list(io)%averaged_data(id)%r4array)) CYCLE

        ALLOCATE(r4temp_sum(1-ng:nx_new+ng, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field_r4(&
              io_block_list(io)%averaged_data(id)%r4array(:,i), &
              r4temp_sum(:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%r4array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %r4array(1-ng:nx_new+ng, nspec_local))

        io_block_list(io)%averaged_data(id)%r4array = r4temp_sum

        DEALLOCATE(r4temp_sum)
      ELSE
        IF (.NOT. ALLOCATED(io_block_list(io)%averaged_data(id)%array)) CYCLE

        ALLOCATE(temp_sum(1-ng:nx_new+ng, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field(&
              io_block_list(io)%averaged_data(id)%array(:,i), &
              temp_sum(:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %array(1-ng:nx_new+ng, nspec_local))

        io_block_list(io)%averaged_data(id)%array = temp_sum

        DEALLOCATE(temp_sum)
      END IF
    END DO

  END SUBROUTINE redistribute_fields



  SUBROUTINE remap_field(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(num), DIMENSION(:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_1d(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths(field_out, ng, n_new(1))

  END SUBROUTINE remap_field



  SUBROUTINE remap_field_r4(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(r4), DIMENSION(:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_1d_r4(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths_r4(field_out, ng, n_new(1))

  END SUBROUTINE remap_field_r4



  SUBROUTINE redistribute_field_1d(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: old_coords, new_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = mpireal

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      old_coords(i) = old_coordinates(cdim(i))
      new_coords(i) = coordinates(cdim(i))
    END DO

    old_min(1) = old_cell_min1(old_coords(1)+1)
    old_max(1) = old_cell_max1(old_coords(1)+1)
    new_min(1) = new_cell_min1(new_coords(1)+1)
    new_max(1) = new_cell_max1(new_coords(1)+1)

    tag = 0
    sendtypes = 0
    recvtypes = 0

    nprocs(1) = SIZE(new_cell_min1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    ! Create array of sendtypes

    DO i = 1, nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    coord = coordinates

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n) - 1
      IF (new_cell_min1(iproc) <= old_min(n) &
          .AND. new_cell_max1(iproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1, nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        END DO
      END IF

      n = 1
      IF (type_max(n) == old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    END DO

    nprocs(1) = SIZE(old_cell_min1)

    ! Create array of recvtypes

    DO i = 1, nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    coord = old_coordinates

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n) - 1
      IF (old_cell_min1(iproc) <= new_min(n) &
          .AND. old_cell_max1(iproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(old_comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1, nd
          new_0(i) = start(i) - ng
        END DO
        DO i = old_0(1), old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        END DO
      END IF

      n = 1
      IF (type_max(n) == new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0, nproc - 1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d



  SUBROUTINE redistribute_field_1d_r4(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(r4), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: old_coords, new_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = MPI_REAL4

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      old_coords(i) = old_coordinates(cdim(i))
      new_coords(i) = coordinates(cdim(i))
    END DO

    old_min(1) = old_cell_min1(old_coords(1)+1)
    old_max(1) = old_cell_max1(old_coords(1)+1)
    new_min(1) = new_cell_min1(new_coords(1)+1)
    new_max(1) = new_cell_max1(new_coords(1)+1)

    tag = 0
    sendtypes = 0
    recvtypes = 0

    nprocs(1) = SIZE(new_cell_min1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    ! Create array of sendtypes

    DO i = 1, nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    coord = coordinates

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n) - 1
      IF (new_cell_min1(iproc) <= old_min(n) &
          .AND. new_cell_max1(iproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1, nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        END DO
      END IF

      n = 1
      IF (type_max(n) == old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    END DO

    nprocs(1) = SIZE(old_cell_min1)

    ! Create array of recvtypes

    DO i = 1, nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    coord = old_coordinates

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n) - 1
      IF (old_cell_min1(iproc) <= new_min(n) &
          .AND. old_cell_max1(iproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(old_comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1, nd
          new_0(i) = start(i) - ng
        END DO
        DO i = old_0(1), old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        END DO
      END IF

      n = 1
      IF (type_max(n) == new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0, nproc - 1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d_r4



  SUBROUTINE get_load(load_x, load_func)

    ! Calculate total load across the X direction

    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load_x
    INTERFACE
      SUBROUTINE load_func(load_x)
        USE constants
        INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load_x
      END SUBROUTINE load_func
    END INTERFACE
    INTEGER :: st

    load_x = 0

    CALL load_func(load_x)

    ! Now have local densities, so add using MPI
    st = SIZE(load_x)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load_x, st, MPI_INTEGER8, MPI_SUM, &
                       comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load_x = push_per_field * load_x
    load_x(ng+1:st-ng) = load_x(ng+1:st-ng) + 1

  END SUBROUTINE get_load



  SUBROUTINE part_load_func(load_x)

    INTEGER(i8), DIMENSION(1-ng:), INTENT(OUT) :: load_x
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x, ispecies

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_grid_min, NOT x_grid_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos - x_grid_min) / dx) + 1
#else
        cell_x = FLOOR((current%part_pos - x_grid_min) / dx + 1.5_num)
#endif
        load_x(cell_x) = load_x(cell_x) + 1

        current => current%next
      END DO
    END DO

  END SUBROUTINE part_load_func



  SUBROUTINE array_load_func(load_x)

    INTEGER(i8), DIMENSION(1-ng:), INTENT(OUT) :: load_x
    INTEGER :: ix, ii

    IF (.NOT.ALLOCATED(npart_per_cell_array)) RETURN

    ii = nx_global_min - 1
    DO ix = 1, nx
      ii = ii + 1
      load_x(ii) = load_x(ii) + npart_per_cell_array(ix)
    END DO

  END SUBROUTINE array_load_func



  SUBROUTINE calculate_breaks(load, nproc, mins, maxs)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(i8), INTENT(IN), DIMENSION(1-ng:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: mins, maxs
    INTEGER :: sz, idim, proc, old, nextra, i, i0, i1, iter, old_maxs, new_maxs
    INTEGER :: ii
    INTEGER(i8) :: total, total_old, load_per_proc_ideal
    INTEGER(i8) :: load_local, load_max, load_min, load_var_best

    sz = SIZE(load) - 2 * ng
    mins = 1
    maxs = sz

    IF (nproc < 2) RETURN

    load_per_proc_ideal = FLOOR(REAL(SUM(load(1:sz)), num) / nproc + 0.5d0, i8)

    proc = 0
    old = 1
    total = 0
    DO idim = 1, sz
      total_old = total
      total = total + load(idim)
      IF (total >= load_per_proc_ideal) THEN
        proc = proc + 1
        IF (load_per_proc_ideal - total_old &
            < total - load_per_proc_ideal) THEN
          maxs(proc) = idim - 1
        ELSE
          maxs(proc) = idim
        END IF
        ! To communicate ghost cell information correctly, each domain must
        ! contain at least ng cells.
        nextra = old - maxs(proc) + ncell_min
        IF (nextra > 0) THEN
          maxs(proc) = maxs(proc) + nextra
        END IF
        IF (proc == nproc - 1) EXIT
        old = maxs(proc)
        total = total - load_per_proc_ideal
      END IF
    END DO

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Backwards
    old = sz
    DO proc = nproc-1, 1, -1
      IF (old - maxs(proc) < ncell_min) THEN
        maxs(proc) = old - ncell_min
      END IF
      old = maxs(proc)
    END DO

    ! Try perturbing the splits by one cell to see if we get a better answer
    load_var_best = HUGE(1)
    DO iter = 1, 1000
      DO i = 1, nproc - 1
        ! Minus
        old_maxs = maxs(i)
        IF (i == 1) THEN
          old = 0
        ELSE
          ii = i - 1
          old = maxs(ii)
        END IF
        new_maxs = old_maxs
        IF (old_maxs - old - 1 >= ng) new_maxs = old_maxs - 1

        IF (new_maxs /= old_maxs) THEN
          maxs(i) = new_maxs
          load_max = -1
          load_min = HUGE(load_min)
          i0 = 1
          DO proc = 1, nproc
            i1 = maxs(proc)
            load_local = SUM(load(i0:i1))
            IF (load_local > load_max) load_max = load_local
            IF (load_local < load_min) load_min = load_local
            i0 = i1 + 1
          END DO
          IF (load_max - load_min < load_var_best) EXIT
          maxs(i) = old_maxs
        END IF

        ! Plus
        old_maxs = maxs(i)
        old = maxs(i+1)
        new_maxs = old_maxs
        IF (old - old_maxs - 1 >= ng) new_maxs = old_maxs + 1

        IF (new_maxs /= old_maxs) THEN
          maxs(i) = new_maxs
          load_max = -1
          load_min = HUGE(load_min)
          i0 = 1
          DO proc = 1, nproc
            i1 = maxs(proc)
            load_local = SUM(load(i0:i1))
            IF (load_local > load_max) load_max = load_local
            IF (load_local < load_min) load_min = load_local
            i0 = i1 + 1
          END DO
          IF (load_max - load_min < load_var_best) EXIT
          maxs(i) = old_maxs
        END IF
      END DO

      IF (load_max - load_min < load_var_best) THEN
        load_var_best = load_max - load_min
      ELSE
        EXIT
      END IF
    END DO

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Backwards
    old = sz
    DO proc = nproc-1, 1, -1
      IF (old - maxs(proc) < ncell_min) THEN
        maxs(proc) = old - ncell_min
      END IF
      old = maxs(proc)
    END DO

    ! Forwards (unnecessary?)
    old = 0
    DO proc = 1, nproc-1
      IF (maxs(proc) - old < ncell_min) THEN
        maxs(proc) = old + ncell_min
      END IF
      old = maxs(proc)
    END DO

    ! Set mins
    mins(1) = 1
    DO proc = 2, nproc
      mins(proc) = maxs(proc-1) + 1
    END DO

  END SUBROUTINE calculate_breaks



  FUNCTION get_particle_processor(part)

    ! This subroutine calculates which processor a given particles resides on

    TYPE(particle), INTENT(IN) :: part
    INTEGER :: get_particle_processor
    INTEGER :: iproc, coords(c_ndims)
    REAL(num) :: minpos, maxpos

    get_particle_processor = -1
    coords = -1

    ! This could be replaced by a bisection method, but for the moment I
    ! just don't care

    DO iproc = 0, nprocx - 1
      IF (iproc == 0) THEN
        minpos = x_grid_mins(iproc) - dx * (0.5_num + png)
      ELSE
        minpos = x_grid_mins(iproc) - dx * 0.5_num
      END IF
      IF (iproc == nprocx - 1) THEN
        maxpos = x_global(cell_x_max(iproc+1)+1) + dx * (png - 0.5_num)
      ELSE
        maxpos = x_global(cell_x_max(iproc+1)+1) - dx * 0.5_num
      END IF
      IF (part%part_pos >= minpos .AND. part%part_pos < maxpos) THEN
        coords(c_ndims) = iproc
        EXIT
      END IF
    END DO

    IF (MINVAL(coords) < 0) THEN
      WRITE(*,*) 'UNLOCATABLE PARTICLE', coords
      RETURN
    END IF
    CALL MPI_CART_RANK(comm, coords, get_particle_processor, errcode)
    ! IF (get_particle_processor /= rank) PRINT *,

  END FUNCTION get_particle_processor



  ! This subroutine is used to rearrange particles over processors
  SUBROUTINE distribute_particles

    ! This subroutine moves particles which are on the wrong processor
    ! to the correct processor.

    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_send
    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_recv
    TYPE(particle), POINTER :: current, next
    INTEGER :: part_proc, iproc, ispecies
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO iproc = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc))
        CALL create_empty_partlist(pointers_recv(iproc))
      END DO

      DO WHILE(ASSOCIATED(current))
        next => current%next
        part_proc = get_particle_processor(current)
        IF (part_proc < 0) THEN
          PRINT *, 'Unlocatable particle on processor', rank, current%part_pos
          CALL abort_code(c_err_bad_value)
          STOP
        END IF
#ifdef PARTICLE_DEBUG
        current%processor = part_proc
#endif
        IF (part_proc /= rank) THEN
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, current)
          CALL add_particle_to_partlist(pointers_send(part_proc), current)
        END IF
        current => next
      END DO

      DO iproc = 0, nproc - 1
        sendcounts(iproc) = pointers_send(iproc)%count
      END DO

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      CALL redblack(pointers_send, pointers_recv, sendcounts, recvcounts)

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      END DO
    END DO

    DEALLOCATE(sendcounts, recvcounts)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles



  SUBROUTINE create_npart_per_cell_array

    INTEGER :: ispecies
    INTEGER :: cell_x
    TYPE(particle), POINTER :: current, next
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    IF (.NOT.ALLOCATED(npart_per_cell_array)) &
        ALLOCATE(npart_per_cell_array(i0:nx+i1))
    npart_per_cell_array(:) = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx) + 1
#else
        cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx + 1.5_num)
#endif
        npart_per_cell_array(cell_x) = &
            npart_per_cell_array(cell_x) + 1

        current => next
      END DO
    END DO

  END SUBROUTINE create_npart_per_cell_array



  SUBROUTINE calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                          load_x_min, load_x_max, &
                                          get_balance)

    REAL(num), INTENT(OUT) :: balance_frac, balance_frac_final
    INTEGER, INTENT(IN) :: load_x_min(:), load_x_max(:)
    LOGICAL, INTENT(IN), OPTIONAL :: get_balance
    REAL(num), ALLOCATABLE :: load_per_cpu(:)
    INTEGER(i8) :: npart_local
    INTEGER :: i, i0, i1, ix, npx
    INTEGER :: ierr
    REAL(num) :: load_local, load_sum, load_max
    LOGICAL :: original_balance

    original_balance = use_injectors
    IF (PRESENT(get_balance)) THEN
      IF (get_balance) original_balance = .TRUE.
    END IF

    IF (original_balance) THEN
      IF (ALLOCATED(npart_per_cell_array)) THEN
        npart_local = SUM(npart_per_cell_array(1:nx))
      ELSE
        npart_local = 0
      END IF
      load_local = REAL(push_per_field * npart_local + nx, num)

      CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, comm, ierr)
      CALL MPI_ALLREDUCE(load_local, load_sum, 1, mpireal, MPI_SUM, comm, ierr)

      load_av = load_sum / nproc

      balance_frac = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))
    END IF

    npx = SIZE(load_x_min)

    ALLOCATE(load_per_cpu(npx))
    load_per_cpu = 0.0_num

    IF (ALLOCATED(npart_per_cell_array)) THEN
      DO i = 1, npx
        i0 = load_x_min(i) - nx_global_min + 1
        i1 = load_x_max(i) - nx_global_min + 1

        IF (i1 < 1 .OR. i0 > nx) CYCLE

        i0 = MAX(i0, 1)
        i1 = MIN(i1, nx)

        DO ix = i0, i1
          load_per_cpu(i) = load_per_cpu(i) &
              + REAL(push_per_field * npart_per_cell_array(ix) + 1, num)
        END DO ! ix
      END DO ! i
    ELSE
      DO i = 1, npx
        i0 = load_x_min(i) - nx_global_min + 1
        i1 = load_x_max(i) - nx_global_min + 1

        IF (i1 < 1 .OR. i0 > nx) CYCLE

        i0 = MAX(i0, 1)
        i1 = MIN(i1, nx)

        DO ix = i0, i1
          load_per_cpu(i) = load_per_cpu(i) + 1.0_num
        END DO ! ix
      END DO ! i
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load_per_cpu, npx, &
                       mpireal, MPI_SUM, comm, ierr)

    load_av = SUM(load_per_cpu) / npx
    load_max = MAXVAL(load_per_cpu)

    balance_frac_final = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))

    DEALLOCATE(load_per_cpu)

  END SUBROUTINE calculate_new_load_imbalance

END MODULE balance
