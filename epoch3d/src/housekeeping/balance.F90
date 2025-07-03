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
  INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_y_min, new_cell_y_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_z_min, new_cell_z_max
  LOGICAL :: overriding
  REAL(num) :: load_av
  INTEGER :: old_comm, old_coordinates(c_ndims)
  INTEGER :: old_slice_coord, new_slice_coord, slice_dir
  LOGICAL :: max_boundary

CONTAINS

  SUBROUTINE get_optimal_layout

    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_y
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_z
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_x_min, p_x_max
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_y_min, p_y_max
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_z_min, p_z_max
    INTEGER :: ii, jj, npx, npy, npz, npyz
    REAL(num) :: dummy, balance_frac, best_balance_frac

    ! On one processor do nothing to save time
    IF (nproc == 1) RETURN
    IF (use_exact_restart) RETURN

    ALLOCATE(load_x(nx_global + 2 * ng))
    ALLOCATE(load_y(ny_global + 2 * ng))
    ALLOCATE(load_z(nz_global + 2 * ng))
    CALL get_load(load_x, load_y, load_z, array_load_func)

    best_balance_frac = -1.0_num

    IF (rank == 0) PRINT*, 'Calculating optimal processor topology'

    DO ii = 1, nproc
      npx = ii
      npyz = nproc / npx
      DO jj = 1, npyz
        npy = jj
        npz = npyz / npy
        IF (npx * npy * npz /= nproc) CYCLE
        IF (nx_global / npx < ncell_min) CYCLE
        IF (ny_global / npy < ncell_min) CYCLE
        IF (nz_global / npz < ncell_min) CYCLE

        ALLOCATE(p_x_min(npx), p_x_max(npx))
        ALLOCATE(p_y_min(npy), p_y_max(npy))
        ALLOCATE(p_z_min(npz), p_z_max(npz))

        CALL calculate_breaks(load_x, npx, p_x_min, p_x_max)
        CALL calculate_breaks(load_y, npy, p_y_min, p_y_max)
        CALL calculate_breaks(load_z, npz, p_z_min, p_z_max)

        CALL calculate_new_load_imbalance(dummy, balance_frac, &
                                          p_x_min, p_x_max, p_y_min, p_y_max, &
                                          p_z_min, p_z_max)

        IF (balance_frac > best_balance_frac) THEN
          best_balance_frac = balance_frac
          nprocx = npx
          nprocy = npy
          nprocz = npz
        END IF

        DEALLOCATE(p_x_min, p_x_max)
        DEALLOCATE(p_y_min, p_y_max)
        DEALLOCATE(p_z_min, p_z_max)
      END DO
    END DO

    DEALLOCATE(load_x, load_y, load_z)

    IF (rank == 0) THEN
      PRINT*, 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
    END IF

  END SUBROUTINE get_optimal_layout



  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_y
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_z
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
    load_local = REAL(push_per_field * npart_local + nx * ny * nz, num)

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
      ALLOCATE(load_y(ny_global + 2 * ng))
      ALLOCATE(load_z(nz_global + 2 * ng))
      CALL get_load(load_x, load_y, load_z, part_load_func)

      ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))
      ALLOCATE(new_cell_y_min(nprocy), new_cell_y_max(nprocy))
      ALLOCATE(new_cell_z_min(nprocz), new_cell_z_max(nprocz))

      CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)
      CALL calculate_breaks(load_y, nprocy, new_cell_y_min, new_cell_y_max)
      CALL calculate_breaks(load_z, nprocz, new_cell_z_min, new_cell_z_max)

      DEALLOCATE(load_x, load_y, load_z)

      IF (.NOT.restarting) THEN
        CALL create_npart_per_cell_array

        CALL calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                          new_cell_x_min, new_cell_x_max, &
                                          new_cell_y_min, new_cell_y_max, &
                                          new_cell_z_min, new_cell_z_max)

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
      DEALLOCATE(new_cell_y_min, new_cell_y_max)
      DEALLOCATE(new_cell_z_min, new_cell_z_max)
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
      load_local = REAL(push_per_field * npart_local + nx * ny * nz, num)

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
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_y
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_z
    REAL(num) :: balance_frac, balance_frac_final, balance_improvement
    LOGICAL :: use_redistribute_domain

    ! On one processor do nothing to save time
    IF (nproc == 1 .OR. .NOT.use_pre_balance) RETURN

    overriding = .TRUE.

    ALLOCATE(load_x(nx_global + 2 * ng))
    ALLOCATE(load_y(ny_global + 2 * ng))
    ALLOCATE(load_z(nz_global + 2 * ng))
    CALL get_load(load_x, load_y, load_z, array_load_func)

    ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))
    ALLOCATE(new_cell_y_min(nprocy), new_cell_y_max(nprocy))
    ALLOCATE(new_cell_z_min(nprocz), new_cell_z_max(nprocz))

    CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)
    CALL calculate_breaks(load_y, nprocy, new_cell_y_min, new_cell_y_max)
    CALL calculate_breaks(load_z, nprocz, new_cell_z_min, new_cell_z_max)

    DEALLOCATE(load_x, load_y, load_z)

    CALL calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                      new_cell_x_min, new_cell_x_max, &
                                      new_cell_y_min, new_cell_y_max, &
                                      new_cell_z_min, new_cell_z_max, .TRUE.)

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
        DEALLOCATE(y_grid_mins, y_grid_maxs)
        DEALLOCATE(z_grid_mins, z_grid_maxs)
        ALLOCATE(x_grid_mins(0:nprocx-1))
        ALLOCATE(x_grid_maxs(0:nprocx-1))
        ALLOCATE(y_grid_mins(0:nprocy-1))
        ALLOCATE(y_grid_maxs(0:nprocy-1))
        ALLOCATE(z_grid_mins(0:nprocz-1))
        ALLOCATE(z_grid_maxs(0:nprocz-1))
      ELSE
        old_comm = comm
        old_coordinates(:) = coordinates(:)
      END IF
      CALL redistribute_domain
    END IF

    IF (ALLOCATED(new_cell_x_min)) THEN
      DEALLOCATE(new_cell_x_min, new_cell_x_max)
      DEALLOCATE(new_cell_y_min, new_cell_y_max)
      DEALLOCATE(new_cell_z_min, new_cell_z_max)
    END IF

  END SUBROUTINE pre_balance_workload



  SUBROUTINE redistribute_domain

    INTEGER, DIMENSION(c_ndims,2) :: domain

    IF (.NOT.ALLOCATED(new_cell_x_min)) RETURN

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor

    domain(1,:) = (/new_cell_x_min(x_coords+1), new_cell_x_max(x_coords+1)/)
    domain(2,:) = (/new_cell_y_min(y_coords+1), new_cell_y_max(y_coords+1)/)
    domain(3,:) = (/new_cell_z_min(z_coords+1), new_cell_z_max(z_coords+1)/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    DEALLOCATE(cell_x_min, cell_x_max)
    DEALLOCATE(cell_y_min, cell_y_max)
    DEALLOCATE(cell_z_min, cell_z_max)
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))
    ALLOCATE(cell_z_min(nprocz), cell_z_max(nprocz))

    ! Copy the new lengths into the permanent variables
    cell_x_min(:) = new_cell_x_min(:)
    cell_x_max(:) = new_cell_x_max(:)
    cell_y_min(:) = new_cell_y_min(:)
    cell_y_max(:) = new_cell_y_max(:)
    cell_z_min(:) = new_cell_z_min(:)
    cell_z_max(:) = new_cell_z_max(:)

    ! Set the new nx, ny, nz
    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)
    n_global_min(1) = nx_global_min
    n_global_max(1) = nx_global_max

    ny_global_min = cell_y_min(y_coords+1)
    ny_global_max = cell_y_max(y_coords+1)
    n_global_min(2) = ny_global_min
    n_global_max(2) = ny_global_max

    nz_global_min = cell_z_min(z_coords+1)
    nz_global_max = cell_z_max(z_coords+1)
    n_global_min(3) = nz_global_min
    n_global_max(3) = nz_global_max

    nx = nx_global_max - nx_global_min + 1
    ny = ny_global_max - ny_global_min + 1
    nz = nz_global_max - nz_global_min + 1

    ! Do X, Y, Z arrays separately because we already have global copies
    DEALLOCATE(x, y, z)
    ALLOCATE(x(1-ng:nx+ng), y(1-ng:ny+ng), z(1-ng:nz+ng))
    DEALLOCATE(xb, yb, zb)
    ALLOCATE(xb(1-ng:nx+ng), yb(1-ng:ny+ng), zb(1-ng:nz+ng))

    CALL setup_grid_x
    CALL setup_grid_y
    CALL setup_grid_z

  END SUBROUTINE redistribute_domain



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the field variables over the new
    ! processor layout. If using a field of your own then set the
    ! redistribute_field subroutine to implement it.

    INTEGER :: nx_new, ny_new, nz_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: temp_sum
    REAL(r4), DIMENSION(:,:,:,:), ALLOCATABLE :: r4temp_sum
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp, temp2
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp_slice
    TYPE(laser_block), POINTER :: laser, laser_next
    TYPE(injector_block), POINTER :: injector, injector_next
    TYPE(particle_species_migration), POINTER :: mg
    TYPE(particle_species), POINTER :: sp
    TYPE(initial_condition_block), POINTER :: ic
    INTEGER :: i, ispecies, io, id, nspec_local, mask

    nx_new = new_domain(1,2) - new_domain(1,1) + 1
    ny_new = new_domain(2,2) - new_domain(2,1) + 1
    nz_new = new_domain(3,2) - new_domain(3,1) + 1

    ! The following code is quite messy and repetitive. Unfortunately, the
    ! F90 standard does not allow the ALLOCATABLE attribute for subroutine
    ! arguments and POINTER arrays are not as fast.

    ! Full domain arrays

    ALLOCATE(temp(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))

    ! Current will be recalculated during the particle push, so there
    ! is no need to copy the contents of the old arrays.
    ! If overriding, then we may not be doing a particle push next
    ! so we still have to balance the arrays.
    ! It is done slightly differently since the arrays may be
    ! a different size.

    IF (overriding) THEN
      ALLOCATE(temp2(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))

      temp2 = jx
      CALL remap_field(temp2, temp)
      DEALLOCATE(jx)
      ALLOCATE(jx(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
      jx(0:nx_new+1, 0:ny_new+1, 0:nz_new+1) = &
          temp(0:nx_new+1, 0:ny_new+1, 0:nz_new+1)

      temp2 = jy
      CALL remap_field(temp2, temp)
      DEALLOCATE(jy)
      ALLOCATE(jy(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
      jy(0:nx_new+1, 0:ny_new+1, 0:nz_new+1) = &
          temp(0:nx_new+1, 0:ny_new+1, 0:nz_new+1)

      temp2 = jz
      CALL remap_field(temp2, temp)
      DEALLOCATE(jz)
      ALLOCATE(jz(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
      jz(0:nx_new+1, 0:ny_new+1, 0:nz_new+1) = &
          temp(0:nx_new+1, 0:ny_new+1, 0:nz_new+1)

      DEALLOCATE(temp2)
    ELSE
      DEALLOCATE(jx)
      DEALLOCATE(jy)
      DEALLOCATE(jz)
      ALLOCATE(jx(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
      ALLOCATE(jy(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
      ALLOCATE(jz(1-jng:nx_new+jng, 1-jng:ny_new+jng, 1-jng:nz_new+jng))
    END IF

    CALL remap_field(ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ex = temp

    CALL remap_field(ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ey = temp

    CALL remap_field(ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ez = temp

    CALL remap_field(bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
    bx = temp

    CALL remap_field(by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
    by = temp

    CALL remap_field(bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
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
        ALLOCATE(mg%fluid_energy(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng))
        mg%fluid_energy = temp

        CALL remap_field(mg%fluid_density, temp)
        DEALLOCATE(mg%fluid_density)
        ALLOCATE(mg%fluid_density(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng))
        mg%fluid_density = temp
      END IF

      IF (sp%background_species) THEN
        CALL remap_field(sp%background_density, temp)
        DEALLOCATE(sp%background_density)
        ALLOCATE(&
            sp%background_density(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng))
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
        ALLOCATE(ic%density(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng))
        ic%density = temp
      END IF

      IF (ALLOCATED(ic%temp)) THEN
        IF (.NOT. ALLOCATED(temp_sum)) &
            ALLOCATE(temp_sum(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng,3))

        CALL remap_field(ic%temp(:,:,:,1), temp_sum(:,:,:,1))
        CALL remap_field(ic%temp(:,:,:,2), temp_sum(:,:,:,2))
        CALL remap_field(ic%temp(:,:,:,3), temp_sum(:,:,:,3))

        DEALLOCATE(ic%temp)
        ALLOCATE(ic%temp(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng,3))
        ic%temp = temp_sum
      END IF

      IF (ALLOCATED(ic%drift)) THEN
        IF (.NOT. ALLOCATED(temp_sum)) &
            ALLOCATE(temp_sum(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng,3))
        CALL remap_field(ic%drift(:,:,:,1), temp_sum(:,:,:,1))
        CALL remap_field(ic%drift(:,:,:,2), temp_sum(:,:,:,2))
        CALL remap_field(ic%drift(:,:,:,3), temp_sum(:,:,:,3))

        DEALLOCATE(ic%drift)
        ALLOCATE(ic%drift(1-ng:nx_new+ng,1-ng:ny_new+ng,1-ng:nz_new+ng,3))
        ic%drift = temp_sum
      END IF
    END DO

    IF (ALLOCATED(temp_sum)) DEALLOCATE(temp_sum)

    IF (cpml_boundaries) THEN
      CALL remap_field(cpml_psi_eyx, temp)
      DEALLOCATE(cpml_psi_eyx)
      ALLOCATE(cpml_psi_eyx(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_eyx = temp

      CALL remap_field(cpml_psi_byx, temp)
      DEALLOCATE(cpml_psi_byx)
      ALLOCATE(cpml_psi_byx(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_byx = temp

      CALL remap_field(cpml_psi_ezx, temp)
      DEALLOCATE(cpml_psi_ezx)
      ALLOCATE(cpml_psi_ezx(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_ezx = temp

      CALL remap_field(cpml_psi_bzx, temp)
      DEALLOCATE(cpml_psi_bzx)
      ALLOCATE(cpml_psi_bzx(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_bzx = temp

      CALL remap_field(cpml_psi_exy, temp)
      DEALLOCATE(cpml_psi_exy)
      ALLOCATE(cpml_psi_exy(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_exy = temp

      CALL remap_field(cpml_psi_bxy, temp)
      DEALLOCATE(cpml_psi_bxy)
      ALLOCATE(cpml_psi_bxy(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_bxy = temp

      CALL remap_field(cpml_psi_ezy, temp)
      DEALLOCATE(cpml_psi_ezy)
      ALLOCATE(cpml_psi_ezy(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_ezy = temp

      CALL remap_field(cpml_psi_bzy, temp)
      DEALLOCATE(cpml_psi_bzy)
      ALLOCATE(cpml_psi_bzy(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_bzy = temp

      CALL remap_field(cpml_psi_exz, temp)
      DEALLOCATE(cpml_psi_exz)
      ALLOCATE(cpml_psi_exz(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_exz = temp

      CALL remap_field(cpml_psi_bxz, temp)
      DEALLOCATE(cpml_psi_bxz)
      ALLOCATE(cpml_psi_bxz(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_bxz = temp

      CALL remap_field(cpml_psi_eyz, temp)
      DEALLOCATE(cpml_psi_eyz)
      ALLOCATE(cpml_psi_eyz(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_eyz = temp

      CALL remap_field(cpml_psi_byz, temp)
      DEALLOCATE(cpml_psi_byz)
      ALLOCATE(cpml_psi_byz(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng))
      cpml_psi_byz = temp

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers(nx_new, new_domain(1,1), new_domain(1,2), &
          ny_new, new_domain(2,1), new_domain(2,2), &
          nz_new, new_domain(3,1), new_domain(3,2))
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

        ALLOCATE(r4temp_sum(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng, &
            nspec_local))

        DO i = 1, nspec_local
          CALL remap_field_r4(&
              io_block_list(io)%averaged_data(id)%r4array(:,:,:,i), &
              r4temp_sum(:,:,:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%r4array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %r4array(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng, &
            nspec_local))

        io_block_list(io)%averaged_data(id)%r4array = r4temp_sum

        DEALLOCATE(r4temp_sum)
      ELSE
        IF (.NOT. ALLOCATED(io_block_list(io)%averaged_data(id)%array)) CYCLE

        ALLOCATE(temp_sum(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng, &
            nspec_local))

        DO i = 1, nspec_local
          CALL remap_field(&
              io_block_list(io)%averaged_data(id)%array(:,:,:,i), &
              temp_sum(:,:,:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %array(1-ng:nx_new+ng, 1-ng:ny_new+ng, 1-ng:nz_new+ng, &
            nspec_local))

        io_block_list(io)%averaged_data(id)%array = temp_sum

        DEALLOCATE(temp_sum)
      END IF
    END DO

    ! Slice in X-direction

    ALLOCATE(temp_slice(1-ng:ny_new+ng, 1-ng:nz_new+ng))

    laser_next => lasers
    DO WHILE(ASSOCIATED(laser_next))
      laser => laser_next
      laser_next => laser%next

      IF (laser%boundary == c_bd_x_min) THEN
        max_boundary = .FALSE.
      ELSE IF (laser%boundary == c_bd_x_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      CALL remap_field_slice(c_dir_x, laser%profile, temp_slice)
      DEALLOCATE(laser%profile)
      ALLOCATE(laser%profile(1-ng:ny_new+ng, 1-ng:nz_new+ng))
      laser%profile = temp_slice

      CALL remap_field_slice(c_dir_x, laser%phase, temp_slice)
      DEALLOCATE(laser%phase)
      ALLOCATE(laser%phase(1-ng:ny_new+ng, 1-ng:nz_new+ng))
      laser%phase = temp_slice
    END DO

    injector_next => injector_list
    DO WHILE(ASSOCIATED(injector_next))
      injector => injector_next
      injector_next => injector%next

      IF (injector%boundary == c_bd_x_min) THEN
        max_boundary = .FALSE.
      ELSE IF (injector%boundary == c_bd_x_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      IF (ALLOCATED(injector%depth)) THEN
        CALL remap_field_slice(c_dir_x, injector%depth, temp_slice)
        DEALLOCATE(injector%depth)
        ALLOCATE(injector%depth(1-ng:ny_new+ng, 1-ng:nz_new+ng))
        injector%depth = temp_slice
      END IF
    END DO

    max_boundary = .FALSE.

    CALL remap_field_slice(c_dir_x, ex_x_min, temp_slice)
    DEALLOCATE(ex_x_min)
    ALLOCATE(ex_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ex_x_min = temp_slice

    CALL remap_field_slice(c_dir_x, ey_x_min, temp_slice)
    DEALLOCATE(ey_x_min)
    ALLOCATE(ey_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ey_x_min = temp_slice

    CALL remap_field_slice(c_dir_x, ez_x_min, temp_slice)
    DEALLOCATE(ez_x_min)
    ALLOCATE(ez_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ez_x_min = temp_slice

    CALL remap_field_slice(c_dir_x, bx_x_min, temp_slice)
    DEALLOCATE(bx_x_min)
    ALLOCATE(bx_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    bx_x_min = temp_slice

    CALL remap_field_slice(c_dir_x, by_x_min, temp_slice)
    DEALLOCATE(by_x_min)
    ALLOCATE(by_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    by_x_min = temp_slice

    CALL remap_field_slice(c_dir_x, bz_x_min, temp_slice)
    DEALLOCATE(bz_x_min)
    ALLOCATE(bz_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    bz_x_min = temp_slice

    max_boundary = .TRUE.

    CALL remap_field_slice(c_dir_x, ex_x_max, temp_slice)
    DEALLOCATE(ex_x_max)
    ALLOCATE(ex_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ex_x_max = temp_slice

    CALL remap_field_slice(c_dir_x, ey_x_max, temp_slice)
    DEALLOCATE(ey_x_max)
    ALLOCATE(ey_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ey_x_max = temp_slice

    CALL remap_field_slice(c_dir_x, ez_x_max, temp_slice)
    DEALLOCATE(ez_x_max)
    ALLOCATE(ez_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    ez_x_max = temp_slice

    CALL remap_field_slice(c_dir_x, bx_x_max, temp_slice)
    DEALLOCATE(bx_x_max)
    ALLOCATE(bx_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    bx_x_max = temp_slice

    CALL remap_field_slice(c_dir_x, by_x_max, temp_slice)
    DEALLOCATE(by_x_max)
    ALLOCATE(by_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    by_x_max = temp_slice

    CALL remap_field_slice(c_dir_x, bz_x_max, temp_slice)
    DEALLOCATE(bz_x_max)
    ALLOCATE(bz_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng))
    bz_x_max = temp_slice

    DEALLOCATE(temp_slice)

    ! Slice in Y-direction

    ALLOCATE(temp_slice(1-ng:nx_new+ng, 1-ng:nz_new+ng))

    laser_next => lasers
    DO WHILE(ASSOCIATED(laser_next))
      laser => laser_next
      laser_next => laser%next

      IF (laser%boundary == c_bd_y_min) THEN
        max_boundary = .FALSE.
      ELSE IF (laser%boundary == c_bd_y_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      CALL remap_field_slice(c_dir_y, laser%profile, temp_slice)
      DEALLOCATE(laser%profile)
      ALLOCATE(laser%profile(1-ng:nx_new+ng, 1-ng:nz_new+ng))
      laser%profile = temp_slice

      CALL remap_field_slice(c_dir_y, laser%phase, temp_slice)
      DEALLOCATE(laser%phase)
      ALLOCATE(laser%phase(1-ng:nx_new+ng, 1-ng:nz_new+ng))
      laser%phase = temp_slice
    END DO

    injector_next => injector_list
    DO WHILE(ASSOCIATED(injector_next))
      injector => injector_next
      injector_next => injector%next

      IF (injector%boundary == c_bd_y_min) THEN
        max_boundary = .FALSE.
      ELSE IF (injector%boundary == c_bd_y_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      IF (ALLOCATED(injector%depth)) THEN
        CALL remap_field_slice(c_dir_y, injector%depth, temp_slice)
        DEALLOCATE(injector%depth)
        ALLOCATE(injector%depth(1-ng:nx_new+ng, 1-ng:nz_new+ng))
        injector%depth = temp_slice
      END IF
    END DO

    max_boundary = .FALSE.

    CALL remap_field_slice(c_dir_y, ex_y_min, temp_slice)
    DEALLOCATE(ex_y_min)
    ALLOCATE(ex_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ex_y_min = temp_slice

    CALL remap_field_slice(c_dir_y, ey_y_min, temp_slice)
    DEALLOCATE(ey_y_min)
    ALLOCATE(ey_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ey_y_min = temp_slice

    CALL remap_field_slice(c_dir_y, ez_y_min, temp_slice)
    DEALLOCATE(ez_y_min)
    ALLOCATE(ez_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ez_y_min = temp_slice

    CALL remap_field_slice(c_dir_y, bx_y_min, temp_slice)
    DEALLOCATE(bx_y_min)
    ALLOCATE(bx_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    bx_y_min = temp_slice

    CALL remap_field_slice(c_dir_y, by_y_min, temp_slice)
    DEALLOCATE(by_y_min)
    ALLOCATE(by_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    by_y_min = temp_slice

    CALL remap_field_slice(c_dir_y, bz_y_min, temp_slice)
    DEALLOCATE(bz_y_min)
    ALLOCATE(bz_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    bz_y_min = temp_slice

    max_boundary = .TRUE.

    CALL remap_field_slice(c_dir_y, ex_y_max, temp_slice)
    DEALLOCATE(ex_y_max)
    ALLOCATE(ex_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ex_y_max = temp_slice

    CALL remap_field_slice(c_dir_y, ey_y_max, temp_slice)
    DEALLOCATE(ey_y_max)
    ALLOCATE(ey_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ey_y_max = temp_slice

    CALL remap_field_slice(c_dir_y, ez_y_max, temp_slice)
    DEALLOCATE(ez_y_max)
    ALLOCATE(ez_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    ez_y_max = temp_slice

    CALL remap_field_slice(c_dir_y, bx_y_max, temp_slice)
    DEALLOCATE(bx_y_max)
    ALLOCATE(bx_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    bx_y_max = temp_slice

    CALL remap_field_slice(c_dir_y, by_y_max, temp_slice)
    DEALLOCATE(by_y_max)
    ALLOCATE(by_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    by_y_max = temp_slice

    CALL remap_field_slice(c_dir_y, bz_y_max, temp_slice)
    DEALLOCATE(bz_y_max)
    ALLOCATE(bz_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng))
    bz_y_max = temp_slice

    DEALLOCATE(temp_slice)

    ! Slice in Z-direction

    ALLOCATE(temp_slice(1-ng:nx_new+ng, 1-ng:ny_new+ng))

    laser_next => lasers
    DO WHILE(ASSOCIATED(laser_next))
      laser => laser_next
      laser_next => laser%next

      IF (laser%boundary == c_bd_z_min) THEN
        max_boundary = .FALSE.
      ELSE IF (laser%boundary == c_bd_z_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      CALL remap_field_slice(c_dir_z, laser%profile, temp_slice)
      DEALLOCATE(laser%profile)
      ALLOCATE(laser%profile(1-ng:nx_new+ng, 1-ng:ny_new+ng))
      laser%profile = temp_slice

      CALL remap_field_slice(c_dir_z, laser%phase, temp_slice)
      DEALLOCATE(laser%phase)
      ALLOCATE(laser%phase(1-ng:nx_new+ng, 1-ng:ny_new+ng))
      laser%phase = temp_slice
    END DO

    injector_next => injector_list
    DO WHILE(ASSOCIATED(injector_next))
      injector => injector_next
      injector_next => injector%next

      IF (injector%boundary == c_bd_z_min) THEN
        max_boundary = .FALSE.
      ELSE IF (injector%boundary == c_bd_z_max) THEN
        max_boundary = .TRUE.
      ELSE
        CYCLE
      END IF

      IF (ALLOCATED(injector%depth)) THEN
        CALL remap_field_slice(c_dir_z, injector%depth, temp_slice)
        DEALLOCATE(injector%depth)
        ALLOCATE(injector%depth(1-ng:nx_new+ng, 1-ng:ny_new+ng))
        injector%depth = temp_slice
      END IF
    END DO

    max_boundary = .FALSE.

    CALL remap_field_slice(c_dir_z, ex_z_min, temp_slice)
    DEALLOCATE(ex_z_min)
    ALLOCATE(ex_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ex_z_min = temp_slice

    CALL remap_field_slice(c_dir_z, ey_z_min, temp_slice)
    DEALLOCATE(ey_z_min)
    ALLOCATE(ey_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ey_z_min = temp_slice

    CALL remap_field_slice(c_dir_z, ez_z_min, temp_slice)
    DEALLOCATE(ez_z_min)
    ALLOCATE(ez_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ez_z_min = temp_slice

    CALL remap_field_slice(c_dir_z, bx_z_min, temp_slice)
    DEALLOCATE(bx_z_min)
    ALLOCATE(bx_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    bx_z_min = temp_slice

    CALL remap_field_slice(c_dir_z, by_z_min, temp_slice)
    DEALLOCATE(by_z_min)
    ALLOCATE(by_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    by_z_min = temp_slice

    CALL remap_field_slice(c_dir_z, bz_z_min, temp_slice)
    DEALLOCATE(bz_z_min)
    ALLOCATE(bz_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    bz_z_min = temp_slice

    max_boundary = .TRUE.

    CALL remap_field_slice(c_dir_z, ex_z_max, temp_slice)
    DEALLOCATE(ex_z_max)
    ALLOCATE(ex_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ex_z_max = temp_slice

    CALL remap_field_slice(c_dir_z, ey_z_max, temp_slice)
    DEALLOCATE(ey_z_max)
    ALLOCATE(ey_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ey_z_max = temp_slice

    CALL remap_field_slice(c_dir_z, ez_z_max, temp_slice)
    DEALLOCATE(ez_z_max)
    ALLOCATE(ez_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    ez_z_max = temp_slice

    CALL remap_field_slice(c_dir_z, bx_z_max, temp_slice)
    DEALLOCATE(bx_z_max)
    ALLOCATE(bx_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    bx_z_max = temp_slice

    CALL remap_field_slice(c_dir_z, by_z_max, temp_slice)
    DEALLOCATE(by_z_max)
    ALLOCATE(by_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    by_z_max = temp_slice

    CALL remap_field_slice(c_dir_z, bz_z_max, temp_slice)
    DEALLOCATE(bz_z_max)
    ALLOCATE(bz_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng))
    bz_z_max = temp_slice

    DEALLOCATE(temp_slice)

    ! Slice in X-direction with an additional index

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:ny_new+ng, 1-ng:nz_new+ng, 3))

        max_boundary = .FALSE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_x, &
              species_list(ispecies)%ext_temp_x_min(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_x_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_x_min(1-ng:ny_new+ng, 1-ng:nz_new+ng, 3))

        species_list(ispecies)%ext_temp_x_min = temp
      END IF

      IF (species_list(ispecies)%bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:ny_new+ng, 1-ng:nz_new+ng, 3))

        max_boundary = .TRUE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_x, &
              species_list(ispecies)%ext_temp_x_max(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_x_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_x_max(1-ng:ny_new+ng, 1-ng:nz_new+ng, 3))

        species_list(ispecies)%ext_temp_x_max = temp
      END IF
    END DO

    IF (ALLOCATED(temp)) DEALLOCATE(temp)

    ! Slice in Y-direction with an additional index

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:nx_new+ng, 1-ng:nz_new+ng, 3))

        max_boundary = .FALSE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_y, &
              species_list(ispecies)%ext_temp_y_min(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_y_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_y_min(1-ng:nx_new+ng, 1-ng:nz_new+ng, 3))

        species_list(ispecies)%ext_temp_y_min = temp
      END IF

      IF (species_list(ispecies)%bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:nx_new+ng, 1-ng:nz_new+ng, 3))

        max_boundary = .TRUE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_y, &
              species_list(ispecies)%ext_temp_y_max(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_y_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_y_max(1-ng:nx_new+ng, 1-ng:nz_new+ng, 3))

        species_list(ispecies)%ext_temp_y_max = temp
      END IF
    END DO

    IF (ALLOCATED(temp)) DEALLOCATE(temp)

    ! Slice in Z-direction with an additional index

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%bc_particle(c_bd_z_min) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:nx_new+ng, 1-ng:ny_new+ng, 3))

        max_boundary = .FALSE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_z, &
              species_list(ispecies)%ext_temp_z_min(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_z_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_z_min(1-ng:nx_new+ng, 1-ng:ny_new+ng, 3))

        species_list(ispecies)%ext_temp_z_min = temp
      END IF

      IF (species_list(ispecies)%bc_particle(c_bd_z_max) == c_bc_thermal) THEN
        IF (.NOT.ALLOCATED(temp)) &
            ALLOCATE(temp(1-ng:nx_new+ng, 1-ng:ny_new+ng, 3))

        max_boundary = .TRUE.

        DO i = 1, 3
          CALL remap_field_slice(c_dir_z, &
              species_list(ispecies)%ext_temp_z_max(:,:,i), temp(:,:,i))
        END DO

        DEALLOCATE(species_list(ispecies)%ext_temp_z_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_z_max(1-ng:nx_new+ng, 1-ng:ny_new+ng, 3))

        species_list(ispecies)%ext_temp_z_max = temp
      END IF
    END DO

    IF (ALLOCATED(temp)) DEALLOCATE(temp)

  END SUBROUTINE redistribute_fields



  SUBROUTINE remap_field_slice(direction, field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(:,:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: field_out
    INTEGER :: i, n
    INTEGER, DIMENSION(c_ndims-1) :: n_new, cdim

    n_new = SHAPE(field_out) - 2 * ng

    n = 1
    DO i = 1, c_ndims
      IF (i == direction) CYCLE
      cdim(n) = c_ndims + 1 - i
      n = n + 1
    END DO

    old_slice_coord = 0
    new_slice_coord = 0

    IF (direction == c_dir_x) THEN
      slice_dir = 3
      IF (max_boundary) THEN
        old_slice_coord = SIZE(cell_x_min) - 1
        new_slice_coord = SIZE(new_cell_x_min) - 1
      END IF
      CALL redistribute_field_2d(field_in, field_out, cdim, &
          cell_y_min, cell_y_max, new_cell_y_min, new_cell_y_max, &
          cell_z_min, cell_z_max, new_cell_z_min, new_cell_z_max)
    ELSE IF (direction == c_dir_y) THEN
      slice_dir = 2
      IF (max_boundary) THEN
        old_slice_coord = SIZE(cell_y_min) - 1
        new_slice_coord = SIZE(new_cell_y_min) - 1
      END IF
      CALL redistribute_field_2d(field_in, field_out, cdim, &
          cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max, &
          cell_z_min, cell_z_max, new_cell_z_min, new_cell_z_max)
    ELSE
      slice_dir = 1
      IF (max_boundary) THEN
        old_slice_coord = SIZE(cell_z_min) - 1
        new_slice_coord = SIZE(new_cell_z_min) - 1
      END IF
      CALL redistribute_field_2d(field_in, field_out, cdim, &
          cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max, &
          cell_y_min, cell_y_max, new_cell_y_min, new_cell_y_max)
    END IF

    CALL do_field_mpi_with_lengths_slice(field_out, direction, ng, n_new(1), &
        n_new(2))

  END SUBROUTINE remap_field_slice



  SUBROUTINE remap_field(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_3d(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max, &
        cell_y_min, cell_y_max, new_cell_y_min, new_cell_y_max, &
        cell_z_min, cell_z_max, new_cell_z_min, new_cell_z_max)

    CALL do_field_mpi_with_lengths(field_out, ng, n_new(1), n_new(2), &
        n_new(3))

  END SUBROUTINE remap_field



  SUBROUTINE remap_field_r4(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(:,:,:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_3d_r4(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max, &
        cell_y_min, cell_y_max, new_cell_y_min, new_cell_y_max, &
        cell_z_min, cell_z_max, new_cell_z_min, new_cell_z_max)

    CALL do_field_mpi_with_lengths_r4(field_out, ng, n_new(1), n_new(2), &
        n_new(3))

  END SUBROUTINE remap_field_r4



  SUBROUTINE redistribute_field_2d(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1, &
      old_cell_min2, old_cell_max2, new_cell_min2, new_cell_max2)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 2
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min2, old_cell_max2
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min2, new_cell_max2
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER :: j, jproc, jnew
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

    old_min(2) = old_cell_min2(old_coords(2)+1)
    old_max(2) = old_cell_max2(old_coords(2)+1)
    new_min(2) = new_cell_min2(new_coords(2)+1)
    new_max(2) = new_cell_max2(new_coords(2)+1)

    tag = 0
    sendtypes = 0
    recvtypes = 0

    nprocs(1) = SIZE(new_cell_min1)
    nprocs(2) = SIZE(new_cell_min2)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))
    nmin(2) = new_cell_min2(1)
    nmax(2) = new_cell_max2(nprocs(2))

    ! Create array of sendtypes

    DO i = 1, nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    coord = coordinates
    coord(slice_dir) = new_slice_coord

    n = 2
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old y_min resides
    ! This could be sped up by using bisection.
    DO jproc = 1, nprocs(n) - 1
      IF (new_cell_min2(jproc) <= old_min(n) &
          .AND. new_cell_max2(jproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      IF (old_coordinates(slice_dir) /= old_slice_coord) EXIT
      coord(cdim(n)) = jproc - 1
      type_max(n) = new_cell_max2(jproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

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
          sendtypes(irank) = create_2d_array_subtype(basetype, n_local, &
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

      n = 2
      IF (type_max(n) == old_max(n)) EXIT
      jproc = jproc + 1
      type_min(n) = new_cell_min2(jproc)
    END DO

    nprocs(1) = SIZE(old_cell_min1)
    nprocs(2) = SIZE(old_cell_min2)

    ! Create array of recvtypes

    DO i = 1, nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    coord = old_coordinates
    coord(slice_dir) = old_slice_coord

    n = 2
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new y_min resides
    ! This could be sped up by using bisection.
    DO jproc = 1, nprocs(n) - 1
      IF (old_cell_min2(jproc) <= new_min(n) &
          .AND. old_cell_max2(jproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      IF (coordinates(slice_dir) /= new_slice_coord) EXIT
      coord(cdim(n)) = jproc - 1
      type_max(n) = old_cell_max2(jproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

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
          recvtypes(irank) = create_2d_array_subtype(basetype, n_local, &
              n_global, start)
        ELSE
          ! New domain is on the same processor as the old domain.
          ! Just copy the region rather than using MPI.
          DO i = 1, nd
            new_0(i) = start(i) - ng
          END DO
          DO j = old_0(2), old_1(2)
            jnew = new_0(2) + j - old_0(2)
            DO i = old_0(1), old_1(1)
              inew = new_0(1) + i - old_0(1)
              field_out(inew,jnew) = field_in(i,j)
            END DO
          END DO
        END IF

        n = 1
        IF (type_max(n) == new_max(n)) EXIT
        iproc = iproc + 1
        type_min(n) = old_cell_min1(iproc)
      END DO

      n = 2
      IF (type_max(n) == new_max(n)) EXIT
      jproc = jproc + 1
      type_min(n) = old_cell_min2(jproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0, nproc - 1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_2d



  SUBROUTINE redistribute_field_3d(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1, &
      old_cell_min2, old_cell_max2, new_cell_min2, new_cell_max2, &
      old_cell_min3, old_cell_max3, new_cell_min3, new_cell_max3)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 3
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min2, old_cell_max2
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min2, new_cell_max2
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min3, old_cell_max3
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min3, new_cell_max3
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER :: j, jproc, jnew
    INTEGER :: k, kproc, knew
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

    old_min(2) = old_cell_min2(old_coords(2)+1)
    old_max(2) = old_cell_max2(old_coords(2)+1)
    new_min(2) = new_cell_min2(new_coords(2)+1)
    new_max(2) = new_cell_max2(new_coords(2)+1)

    old_min(3) = old_cell_min3(old_coords(3)+1)
    old_max(3) = old_cell_max3(old_coords(3)+1)
    new_min(3) = new_cell_min3(new_coords(3)+1)
    new_max(3) = new_cell_max3(new_coords(3)+1)

    tag = 0
    sendtypes = 0
    recvtypes = 0

    nprocs(1) = SIZE(new_cell_min1)
    nprocs(2) = SIZE(new_cell_min2)
    nprocs(3) = SIZE(new_cell_min3)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))
    nmin(2) = new_cell_min2(1)
    nmax(2) = new_cell_max2(nprocs(2))
    nmin(3) = new_cell_min3(1)
    nmax(3) = new_cell_max3(nprocs(3))

    ! Create array of sendtypes

    DO i = 1, nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    coord = coordinates

    n = 3
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old z_min resides
    ! This could be sped up by using bisection.
    DO kproc = 1, nprocs(n) - 1
      IF (new_cell_min3(kproc) <= old_min(n) &
          .AND. new_cell_max3(kproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = kproc - 1
      type_max(n) = new_cell_max3(kproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      n = 2
      type_min(n) = old_min(n)
      type_max(n) = old_min(n)

      ! Find the new processor on which the old y_min resides
      ! This could be sped up by using bisection.
      DO jproc = 1, nprocs(n) - 1
        IF (new_cell_min2(jproc) <= old_min(n) &
            .AND. new_cell_max2(jproc) >= old_min(n)) EXIT
      END DO

      DO WHILE(type_max(n) <= old_max(n))
        coord(cdim(n)) = jproc - 1
        type_max(n) = new_cell_max2(jproc)
        IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

        ng0 = 0
        ng1 = 0
        IF (type_min(n) == nmin(n)) ng0 = ng
        IF (type_max(n) == nmax(n)) ng1 = ng

        n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
        start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

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
            sendtypes(irank) = create_3d_array_subtype(basetype, n_local, &
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

        n = 2
        IF (type_max(n) == old_max(n)) EXIT
        jproc = jproc + 1
        type_min(n) = new_cell_min2(jproc)
      END DO

      n = 3
      IF (type_max(n) == old_max(n)) EXIT
      kproc = kproc + 1
      type_min(n) = new_cell_min3(kproc)
    END DO

    nprocs(1) = SIZE(old_cell_min1)
    nprocs(2) = SIZE(old_cell_min2)
    nprocs(3) = SIZE(old_cell_min3)

    ! Create array of recvtypes

    DO i = 1, nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    coord = old_coordinates

    n = 3
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new z_min resides
    ! This could be sped up by using bisection.
    DO kproc = 1, nprocs(n) - 1
      IF (old_cell_min3(kproc) <= new_min(n) &
          .AND. old_cell_max3(kproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = kproc - 1
      type_max(n) = old_cell_max3(kproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      n = 2
      type_min(n) = new_min(n)
      type_max(n) = new_min(n)

      ! Find the old processor on which the new y_min resides
      ! This could be sped up by using bisection.
      DO jproc = 1, nprocs(n) - 1
        IF (old_cell_min2(jproc) <= new_min(n) &
            .AND. old_cell_max2(jproc) >= new_min(n)) EXIT
      END DO

      DO WHILE(type_max(n) <= new_max(n))
        coord(cdim(n)) = jproc - 1
        type_max(n) = old_cell_max2(jproc)
        IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

        ng0 = 0
        ng1 = 0
        IF (type_min(n) == nmin(n)) ng0 = ng
        IF (type_max(n) == nmax(n)) ng1 = ng

        n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
        start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

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
            recvtypes(irank) = create_3d_array_subtype(basetype, n_local, &
                n_global, start)
          ELSE
            ! New domain is on the same processor as the old domain.
            ! Just copy the region rather than using MPI.
            DO i = 1, nd
              new_0(i) = start(i) - ng
            END DO
            DO k = old_0(3), old_1(3)
              knew = new_0(3) + k - old_0(3)
              DO j = old_0(2), old_1(2)
                jnew = new_0(2) + j - old_0(2)
                DO i = old_0(1), old_1(1)
                  inew = new_0(1) + i - old_0(1)
                  field_out(inew,jnew,knew) = field_in(i,j,k)
                END DO
              END DO
            END DO
          END IF

          n = 1
          IF (type_max(n) == new_max(n)) EXIT
          iproc = iproc + 1
          type_min(n) = old_cell_min1(iproc)
        END DO

        n = 2
        IF (type_max(n) == new_max(n)) EXIT
        jproc = jproc + 1
        type_min(n) = old_cell_min2(jproc)
      END DO

      n = 3
      IF (type_max(n) == new_max(n)) EXIT
      kproc = kproc + 1
      type_min(n) = old_cell_min3(kproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0, nproc - 1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_3d



  SUBROUTINE redistribute_field_3d_r4(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1, &
      old_cell_min2, old_cell_max2, new_cell_min2, new_cell_max2, &
      old_cell_min3, old_cell_max3, new_cell_min3, new_cell_max3)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 3
    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min2, old_cell_max2
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min2, new_cell_max2
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min3, old_cell_max3
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min3, new_cell_max3
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER :: j, jproc, jnew
    INTEGER :: k, kproc, knew
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

    old_min(2) = old_cell_min2(old_coords(2)+1)
    old_max(2) = old_cell_max2(old_coords(2)+1)
    new_min(2) = new_cell_min2(new_coords(2)+1)
    new_max(2) = new_cell_max2(new_coords(2)+1)

    old_min(3) = old_cell_min3(old_coords(3)+1)
    old_max(3) = old_cell_max3(old_coords(3)+1)
    new_min(3) = new_cell_min3(new_coords(3)+1)
    new_max(3) = new_cell_max3(new_coords(3)+1)

    tag = 0
    sendtypes = 0
    recvtypes = 0

    nprocs(1) = SIZE(new_cell_min1)
    nprocs(2) = SIZE(new_cell_min2)
    nprocs(3) = SIZE(new_cell_min3)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))
    nmin(2) = new_cell_min2(1)
    nmax(2) = new_cell_max2(nprocs(2))
    nmin(3) = new_cell_min3(1)
    nmax(3) = new_cell_max3(nprocs(3))

    ! Create array of sendtypes

    DO i = 1, nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    coord = coordinates

    n = 3
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old z_min resides
    ! This could be sped up by using bisection.
    DO kproc = 1, nprocs(n) - 1
      IF (new_cell_min3(kproc) <= old_min(n) &
          .AND. new_cell_max3(kproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = kproc - 1
      type_max(n) = new_cell_max3(kproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      n = 2
      type_min(n) = old_min(n)
      type_max(n) = old_min(n)

      ! Find the new processor on which the old y_min resides
      ! This could be sped up by using bisection.
      DO jproc = 1, nprocs(n) - 1
        IF (new_cell_min2(jproc) <= old_min(n) &
            .AND. new_cell_max2(jproc) >= old_min(n)) EXIT
      END DO

      DO WHILE(type_max(n) <= old_max(n))
        coord(cdim(n)) = jproc - 1
        type_max(n) = new_cell_max2(jproc)
        IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

        ng0 = 0
        ng1 = 0
        IF (type_min(n) == nmin(n)) ng0 = ng
        IF (type_max(n) == nmax(n)) ng1 = ng

        n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
        start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

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
            sendtypes(irank) = create_3d_array_subtype(basetype, n_local, &
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

        n = 2
        IF (type_max(n) == old_max(n)) EXIT
        jproc = jproc + 1
        type_min(n) = new_cell_min2(jproc)
      END DO

      n = 3
      IF (type_max(n) == old_max(n)) EXIT
      kproc = kproc + 1
      type_min(n) = new_cell_min3(kproc)
    END DO

    nprocs(1) = SIZE(old_cell_min1)
    nprocs(2) = SIZE(old_cell_min2)
    nprocs(3) = SIZE(old_cell_min3)

    ! Create array of recvtypes

    DO i = 1, nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    coord = old_coordinates

    n = 3
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new z_min resides
    ! This could be sped up by using bisection.
    DO kproc = 1, nprocs(n) - 1
      IF (old_cell_min3(kproc) <= new_min(n) &
          .AND. old_cell_max3(kproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = kproc - 1
      type_max(n) = old_cell_max3(kproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      n = 2
      type_min(n) = new_min(n)
      type_max(n) = new_min(n)

      ! Find the old processor on which the new y_min resides
      ! This could be sped up by using bisection.
      DO jproc = 1, nprocs(n) - 1
        IF (old_cell_min2(jproc) <= new_min(n) &
            .AND. old_cell_max2(jproc) >= new_min(n)) EXIT
      END DO

      DO WHILE(type_max(n) <= new_max(n))
        coord(cdim(n)) = jproc - 1
        type_max(n) = old_cell_max2(jproc)
        IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

        ng0 = 0
        ng1 = 0
        IF (type_min(n) == nmin(n)) ng0 = ng
        IF (type_max(n) == nmax(n)) ng1 = ng

        n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
        start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

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
            recvtypes(irank) = create_3d_array_subtype(basetype, n_local, &
                n_global, start)
          ELSE
            ! New domain is on the same processor as the old domain.
            ! Just copy the region rather than using MPI.
            DO i = 1, nd
              new_0(i) = start(i) - ng
            END DO
            DO k = old_0(3), old_1(3)
              knew = new_0(3) + k - old_0(3)
              DO j = old_0(2), old_1(2)
                jnew = new_0(2) + j - old_0(2)
                DO i = old_0(1), old_1(1)
                  inew = new_0(1) + i - old_0(1)
                  field_out(inew,jnew,knew) = field_in(i,j,k)
                END DO
              END DO
            END DO
          END IF

          n = 1
          IF (type_max(n) == new_max(n)) EXIT
          iproc = iproc + 1
          type_min(n) = old_cell_min1(iproc)
        END DO

        n = 2
        IF (type_max(n) == new_max(n)) EXIT
        jproc = jproc + 1
        type_min(n) = old_cell_min2(jproc)
      END DO

      n = 3
      IF (type_max(n) == new_max(n)) EXIT
      kproc = kproc + 1
      type_min(n) = old_cell_min3(kproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0, nproc - 1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_3d_r4



  SUBROUTINE get_load_x(load)

    ! Calculate total load across the X direction
    ! Summed in the Y,Z directions

    INTEGER(i8), DIMENSION(:), INTENT(INOUT) :: load
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, st

    load = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_grid_min, NOT x_grid_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(1) - x_grid_min) / dx) + 1 + ng
#else
        cell = FLOOR((current%part_pos(1) - x_grid_min) / dx + 1.5_num) + ng
#endif
        load(cell) = load(cell) + 1

        current => current%next
      END DO
    END DO

    ! Now have local densities, so add using MPI
    st = SIZE(load)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load, st, MPI_INTEGER8, MPI_SUM, &
                       comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * load
    load(ng+1:st-ng) = load(ng+1:st-ng) + ny_global * nz_global

  END SUBROUTINE get_load_x



  SUBROUTINE get_load_y(load)

    ! Calculate total load across the Y direction
    ! Summed in the X,Z directions

    INTEGER(i8), DIMENSION(:), INTENT(INOUT) :: load
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, st

    load = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so y_grid_min, NOT y_grid_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(2) - y_grid_min) / dy) + 1 + ng
#else
        cell = FLOOR((current%part_pos(2) - y_grid_min) / dy + 1.5_num) + ng
#endif
        load(cell) = load(cell) + 1

        current => current%next
      END DO
    END DO

    ! Now have local densities, so add using MPI
    st = SIZE(load)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load, st, MPI_INTEGER8, MPI_SUM, &
                       comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * load
    load(ng+1:st-ng) = load(ng+1:st-ng) + nx_global * nz_global

  END SUBROUTINE get_load_y



  SUBROUTINE get_load_z(load)

    ! Calculate total load across the Z direction
    ! Summed in the X,Y directions

    INTEGER(i8), DIMENSION(:), INTENT(INOUT) :: load
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, st

    load = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so z_grid_min, NOT z_grid_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(3) - z_grid_min) / dz) + 1 + ng
#else
        cell = FLOOR((current%part_pos(3) - z_grid_min) / dz + 1.5_num) + ng
#endif
        load(cell) = load(cell) + 1

        current => current%next
      END DO
    END DO

    ! Now have local densities, so add using MPI
    st = SIZE(load)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load, st, MPI_INTEGER8, MPI_SUM, &
                       comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * load
    load(ng+1:st-ng) = load(ng+1:st-ng) + nx_global * ny_global

  END SUBROUTINE get_load_z



  SUBROUTINE get_load(load_x, load_y, load_z, load_func)

    ! Calculate total load across the X,Y,Z directions

    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load_x, load_y, load_z
    INTERFACE
      SUBROUTINE load_func(load_x, load_y, load_z)
        USE constants
        INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load_x, load_y, load_z
      END SUBROUTINE load_func
    END INTERFACE
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load
    INTEGER :: st, sx, sy, sz
    INTEGER :: i0, i1, j0, j1, k0, k1

    load_x = 0
    load_y = 0
    load_z = 0

    CALL load_func(load_x, load_y, load_z)

    ! Now have local densities, so add using MPI
    sx = SIZE(load_x)
    sy = SIZE(load_y)
    sz = SIZE(load_z)
    st = sx + sy + sz
    ALLOCATE(load(st))

    i0 = 1;       i1 = i0 + sx - 1
    j0 = i0 + i1; j1 = j0 + sy - 1
    k0 = i0 + j1; k1 = k0 + sz - 1

    load(i0:i1) = load_x
    load(j0:j1) = load_y
    load(k0:k1) = load_z

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load, st, MPI_INTEGER8, MPI_SUM, &
                       comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load_x = push_per_field * load(i0:i1)
    load_x(ng+1:sx-ng) = load_x(ng+1:sx-ng) + ny_global * nz_global

    load_y = push_per_field * load(j0:j1)
    load_y(ng+1:sy-ng) = load_y(ng+1:sy-ng) + nx_global * nz_global

    load_z = push_per_field * load(k0:k1)
    load_z(ng+1:sz-ng) = load_z(ng+1:sz-ng) + nx_global * ny_global

    DEALLOCATE(load)

  END SUBROUTINE get_load



  SUBROUTINE part_load_func(load_x, load_y, load_z)

    INTEGER(i8), DIMENSION(1-ng:), INTENT(OUT) :: load_x, load_y, load_z
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x, cell_y, cell_z, ispecies

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_grid_min, NOT x_grid_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos(1) - x_grid_min) / dx) + 1
        cell_y = FLOOR((current%part_pos(2) - y_grid_min) / dy) + 1
        cell_z = FLOOR((current%part_pos(3) - z_grid_min) / dz) + 1
#else
        cell_x = FLOOR((current%part_pos(1) - x_grid_min) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min) / dy + 1.5_num)
        cell_z = FLOOR((current%part_pos(3) - z_grid_min) / dz + 1.5_num)
#endif
        load_x(cell_x) = load_x(cell_x) + 1
        load_y(cell_y) = load_y(cell_y) + 1
        load_z(cell_z) = load_z(cell_z) + 1

        current => current%next
      END DO
    END DO

  END SUBROUTINE part_load_func



  SUBROUTINE array_load_func(load_x, load_y, load_z)

    INTEGER(i8), DIMENSION(1-ng:), INTENT(OUT) :: load_x, load_y, load_z
    INTEGER :: ix, iy, iz, ii, jj, kk

    IF (.NOT.ALLOCATED(npart_per_cell_array)) RETURN

    kk = nz_global_min - 1
    DO iz = 1, nz
      kk = kk + 1
      jj = ny_global_min - 1
      DO iy = 1, ny
        jj = jj + 1
        ii = nx_global_min - 1
        DO ix = 1, nx
          ii = ii + 1
          load_x(ii) = load_x(ii) + npart_per_cell_array(ix,iy,iz)
          load_y(jj) = load_y(jj) + npart_per_cell_array(ix,iy,iz)
          load_z(kk) = load_z(kk) + npart_per_cell_array(ix,iy,iz)
        END DO
      END DO
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
      IF (part%part_pos(1) >= minpos .AND. part%part_pos(1) < maxpos) THEN
        coords(c_ndims) = iproc
        EXIT
      END IF
    END DO

    DO iproc = 0, nprocy - 1
      IF (iproc == 0) THEN
        minpos = y_grid_mins(iproc) - dy * (0.5_num + png)
      ELSE
        minpos = y_grid_mins(iproc) - dy * 0.5_num
      END IF
      IF (iproc == nprocy - 1) THEN
        maxpos = y_global(cell_y_max(iproc+1)+1) + dy * (png - 0.5_num)
      ELSE
        maxpos = y_global(cell_y_max(iproc+1)+1) - dy * 0.5_num
      END IF
      IF (part%part_pos(2) >= minpos .AND. part%part_pos(2) < maxpos) THEN
        coords(c_ndims-1) = iproc
        EXIT
      END IF
    END DO

    DO iproc = 0, nprocz - 1
      IF (iproc == 0) THEN
        minpos = z_grid_mins(iproc) - dz * (0.5_num + png)
      ELSE
        minpos = z_grid_mins(iproc) - dz * 0.5_num
      END IF
      IF (iproc == nprocz - 1) THEN
        maxpos = z_global(cell_z_max(iproc+1)+1) + dz * (png - 0.5_num)
      ELSE
        maxpos = z_global(cell_z_max(iproc+1)+1) - dz * 0.5_num
      END IF
      IF (part%part_pos(3) >= minpos .AND. part%part_pos(3) < maxpos) THEN
        coords(c_ndims-2) = iproc
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
    INTEGER :: cell_x, cell_y, cell_z
    TYPE(particle), POINTER :: current, next
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    IF (.NOT.ALLOCATED(npart_per_cell_array)) &
        ALLOCATE(npart_per_cell_array(i0:nx+i1,i0:ny+i1,i0:nz+i1))
    npart_per_cell_array(:,:,:) = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx) + 1
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy) + 1
        cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz) + 1
#else
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
        cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz + 1.5_num)
#endif
        npart_per_cell_array(cell_x,cell_y,cell_z) = &
            npart_per_cell_array(cell_x,cell_y,cell_z) + 1

        current => next
      END DO
    END DO

  END SUBROUTINE create_npart_per_cell_array



  SUBROUTINE calculate_new_load_imbalance(balance_frac, balance_frac_final, &
                                          load_x_min, load_x_max, &
                                          load_y_min, load_y_max, &
                                          load_z_min, load_z_max, &
                                          get_balance)

    REAL(num), INTENT(OUT) :: balance_frac, balance_frac_final
    INTEGER, INTENT(IN) :: load_x_min(:), load_x_max(:)
    INTEGER, INTENT(IN) :: load_y_min(:), load_y_max(:)
    INTEGER, INTENT(IN) :: load_z_min(:), load_z_max(:)
    LOGICAL, INTENT(IN), OPTIONAL :: get_balance
    REAL(num), ALLOCATABLE :: load_per_cpu(:,:,:)
    INTEGER(i8) :: npart_local
    INTEGER :: i, i0, i1, ix, npx
    INTEGER :: j, j0, j1, iy, npy
    INTEGER :: k, k0, k1, iz, npz
    INTEGER :: ierr
    REAL(num) :: load_local, load_sum, load_max
    LOGICAL :: original_balance

    original_balance = use_injectors
    IF (PRESENT(get_balance)) THEN
      IF (get_balance) original_balance = .TRUE.
    END IF

    IF (original_balance) THEN
      IF (ALLOCATED(npart_per_cell_array)) THEN
        npart_local = SUM(npart_per_cell_array(1:nx,1:ny,1:nz))
      ELSE
        npart_local = 0
      END IF
      load_local = REAL(push_per_field * npart_local + nx * ny * nz, num)

      CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, comm, ierr)
      CALL MPI_ALLREDUCE(load_local, load_sum, 1, mpireal, MPI_SUM, comm, ierr)

      load_av = load_sum / nproc

      balance_frac = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))
    END IF

    npx = SIZE(load_x_min)
    npy = SIZE(load_y_min)
    npz = SIZE(load_z_min)

    ALLOCATE(load_per_cpu(npx,npy,npz))
    load_per_cpu = 0.0_num

    IF (ALLOCATED(npart_per_cell_array)) THEN
      DO k = 1, npz
        k0 = load_z_min(k) - nz_global_min + 1
        k1 = load_z_max(k) - nz_global_min + 1

        IF (k1 < 1 .OR. k0 > nz) CYCLE

        k0 = MAX(k0, 1)
        k1 = MIN(k1, nz)

        DO j = 1, npy
          j0 = load_y_min(j) - ny_global_min + 1
          j1 = load_y_max(j) - ny_global_min + 1

          IF (j1 < 1 .OR. j0 > ny) CYCLE

          j0 = MAX(j0, 1)
          j1 = MIN(j1, ny)

          DO i = 1, npx
            i0 = load_x_min(i) - nx_global_min + 1
            i1 = load_x_max(i) - nx_global_min + 1

            IF (i1 < 1 .OR. i0 > nx) CYCLE

            i0 = MAX(i0, 1)
            i1 = MIN(i1, nx)

            DO iz = k0, k1
            DO iy = j0, j1
            DO ix = i0, i1
              load_per_cpu(i,j,k) = load_per_cpu(i,j,k) &
                  + REAL(push_per_field * npart_per_cell_array(ix,iy,iz) + 1, &
                         num)
            END DO ! ix
            END DO ! iy
            END DO ! iz
          END DO ! i
        END DO ! j
      END DO ! k
    ELSE
      DO k = 1, npz
        k0 = load_z_min(k) - nz_global_min + 1
        k1 = load_z_max(k) - nz_global_min + 1

        IF (k1 < 1 .OR. k0 > nz) CYCLE

        k0 = MAX(k0, 1)
        k1 = MIN(k1, nz)

        DO j = 1, npy
          j0 = load_y_min(j) - ny_global_min + 1
          j1 = load_y_max(j) - ny_global_min + 1

          IF (j1 < 1 .OR. j0 > ny) CYCLE

          j0 = MAX(j0, 1)
          j1 = MIN(j1, ny)

          DO i = 1, npx
            i0 = load_x_min(i) - nx_global_min + 1
            i1 = load_x_max(i) - nx_global_min + 1

            IF (i1 < 1 .OR. i0 > nx) CYCLE

            i0 = MAX(i0, 1)
            i1 = MIN(i1, nx)

            DO iz = k0, k1
            DO iy = j0, j1
            DO ix = i0, i1
              load_per_cpu(i,j,k) = load_per_cpu(i,j,k) + 1.0_num
            END DO ! ix
            END DO ! iy
            END DO ! iz
          END DO ! i
        END DO ! j
      END DO ! k
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load_per_cpu, npx * npy * npz, &
                       mpireal, MPI_SUM, comm, ierr)

    load_av = SUM(load_per_cpu) / (npx * npy * npz)
    load_max = MAXVAL(load_per_cpu)

    balance_frac_final = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))

    DEALLOCATE(load_per_cpu)

  END SUBROUTINE calculate_new_load_imbalance

END MODULE balance
