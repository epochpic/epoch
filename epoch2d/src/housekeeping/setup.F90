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

MODULE setup

  USE fields
  USE mpi_subtype_control
  USE version_data
  USE welcome
  USE shunt
  USE laser
  USE injectors
  USE window
  USE timer
  USE helper
  USE balance
  USE mpi_routines
  USE sdf
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files, flush_stat_file
  PUBLIC :: setup_species, after_deck_last, set_dt
  PUBLIC :: read_cpu_split, pre_load_balance
  PUBLIC :: open_status_file, close_status_file

  TYPE(particle), POINTER, SAVE :: iterator_list
#ifndef NO_IO
  CHARACTER(LEN=c_max_path_length), SAVE :: stat_file
#endif
  LOGICAL :: got_x_grid_min = .FALSE.
  LOGICAL :: got_y_grid_min = .FALSE.
  REAL(num) :: x_grid_min_val
  REAL(num) :: y_grid_min_val

CONTAINS

  SUBROUTINE minimal_init

    INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
    INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
    INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)

    IF (num == r4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      realsize = 4
      mpireal = MPI_REAL4
    ELSE IF (num == r8) THEN
      realsize = 8
      mpireal = MPI_REAL8
    ELSE
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot determine size of real'
      END IF
      CALL abort_code(c_err_terminate)
      STOP
    END IF

    dt_plasma_frequency = 0.0_num
    dt_multiplier = 0.95_num
    stdout_frequency = 0
    cpml_thickness = 6
    cpml_kappa_max = 20.0_num
    cpml_a_max = 0.15_num
    cpml_sigma_max = 0.7_num
    cpml_x_min_offset = 0
    cpml_x_max_offset = 0
    cpml_y_min_offset = 0
    cpml_y_max_offset = 0

    npart_global = -1
    smooth_currents = .FALSE.
    use_balance = .FALSE.
    use_random_seed = .FALSE.
    use_offset_grid = .FALSE.
    use_particle_lists = .FALSE.
    use_multiphoton = .TRUE.
    use_bsi = .TRUE.
    need_random_state = .FALSE.
    force_first_to_be_restartable = .FALSE.
    force_final_to_be_restartable = .FALSE.
    full_dump_every = -1
    restart_dump_every = -1
    nsteps = -1
    t_end = HUGE(1.0_num)
    particles_max_id = 0
    n_zeros = 4

    laser_inject_local = 0.0_num
    laser_absorb_local = 0.0_num
    field_energy_lost_local = 0.0_num
    old_elapsed_time = 0.0_num
    window_offset = 0.0_num

    NULLIFY(lasers)
    NULLIFY(dist_fns)
    NULLIFY(io_block_list)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL timer_init

    ! This array is true if a field component is staggered in the
    ! given direction.
    stagger = .FALSE.
    stagger(c_dir_x,c_stagger_ex) = .TRUE.
    stagger(c_dir_y,c_stagger_ey) = .TRUE.
    stagger(c_dir_z,c_stagger_ez) = .TRUE.

    stagger(c_dir_x:c_dir_z,c_stagger_bx) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_by) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_bz) = .TRUE.
    stagger(c_dir_x,c_stagger_bx) = .FALSE.
    stagger(c_dir_y,c_stagger_by) = .FALSE.
    stagger(c_dir_z,c_stagger_bz) = .FALSE.

    CALL set_tokenizer_stagger(c_stagger_centre)

    ALLOCATE(x(1), y(1))
    x = 0.0_num
    y = 0.0_num

    ALLOCATE(xb(1), yb(1))
    xb = 0.0_num
    yb = 0.0_num

    CALL eval_stack_init

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    CALL setup_grid
    CALL set_initial_values
    CALL setup_domain_dependent_boundaries

  END SUBROUTINE after_control



  SUBROUTINE setup_grid

    INTEGER :: ix, iy
    REAL(num) :: xb_min, yb_min

    length_x = x_max - x_min
    dx = length_x / REAL(nx_global-2*cpml_thickness, num)
    x_grid_min = x_min - dx * cpml_thickness

    length_y = y_max - y_min
    dy = length_y / REAL(ny_global-2*cpml_thickness, num)
    y_grid_min = y_min - dy * cpml_thickness

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_grid_min
    yb_min = y_grid_min
    x_grid_min = x_grid_min + dx / 2.0_num
    y_grid_min = y_grid_min + dy / 2.0_num

    IF (got_x_grid_min) x_grid_min = x_grid_min_val
    IF (got_y_grid_min) y_grid_min = y_grid_min_val

    ! Setup global grid
    DO ix = 1-ng, nx_global + ng
      x_global(ix) = x_grid_min + (ix - 1) * dx
      xb_global(ix) = xb_min + (ix - 1) * dx
      xb_offset_global(ix) = xb_global(ix)
    END DO
    x_grid_max = x_global(nx_global)

    DO iy = 1-ng, ny_global + ng
      y_global(iy) = y_grid_min + (iy - 1) * dy
      yb_global(iy) = yb_min + (iy - 1) * dy
      yb_offset_global(iy) = yb_global(iy)
    END DO
    y_grid_max = y_global(ny_global)

    CALL setup_grid_x
    CALL setup_grid_y

  END SUBROUTINE setup_grid



  SUBROUTINE after_deck_last

    INTEGER :: i

    CALL setup_data_averaging
    CALL setup_field_boundaries

    cpml_x_min = .FALSE.
    cpml_x_max = .FALSE.
    cpml_y_min = .FALSE.
    cpml_y_max = .FALSE.

    IF (cpml_boundaries) THEN
      CALL allocate_cpml_fields
      CALL set_cpml_helpers(nx, nx_global_min, nx_global_max, &
          ny, ny_global_min, ny_global_max)
    ELSE
      cpml_thickness = 0
      cpml_kappa_max = 1.0_num
      cpml_a_max = 0.0_num
      cpml_sigma_max = 0.0_num
      DO i = 1, n_io_blocks
        io_block_list(i)%dumpmask(c_dump_cpml_psi_eyx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_ezx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_byx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bzx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_exy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_ezy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bxy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bzy) = c_io_never
      END DO
    END IF

    IF (dump_field_energy_lost .AND. cpml_boundaries .AND. rank == 0) THEN
      PRINT*, '*** WARNING ***'
      PRINT*, 'The "field_energy_lost" diagnostic only measures energy leaving'
      PRINT*, 'through simple_laser and simple_outflow boundaries. Energy'
      PRINT*, 'absorbed by CPML boundaries is not included.'
    END IF

  END SUBROUTINE after_deck_last



  SUBROUTINE setup_data_averaging()

    INTEGER :: io, ib, nspec_local, nstep_average, mask
    REAL(num) :: dt_average
    TYPE(averaged_data_block), POINTER :: avg

    IF (.NOT. any_average) RETURN

    averaged_var_dims = 1
    averaged_var_dims(c_dump_ekflux) = 6
    averaged_var_dims(c_dump_poynt_flux) = 3

    DO io = 1, n_io_blocks
      IF (io_block_list(io)%any_average) THEN
        dt_min_average = t_end
      ELSE
        dt_min_average = -1.0_num
      END IF
    END DO

    DO io = 1, num_vars_to_dump
      ib = averaged_var_block(io)
      IF (ib /= 0) THEN
        dt_average  = io_block_list(ib)%dt_average
        nstep_average = io_block_list(ib)%nstep_average

        IF (nstep_average > 0 .AND. dt_average > 0) THEN
          io_block_list(ib)%dt_min_average = &
              MIN(io_block_list(ib)%dt_min_average, &
              dt_average / REAL(nstep_average, num))
        END IF

        mask = io_block_list(ib)%dumpmask(io)
        nspec_local = 0
        IF (IAND(mask, c_io_no_sum) == 0) &
            nspec_local = 1
        IF (IAND(mask, c_io_species) /= 0) &
            nspec_local = nspec_local + n_species

        IF (nspec_local <= 0) CYCLE
        nspec_local = nspec_local * averaged_var_dims(io)

        avg => io_block_list(ib)%averaged_data(io)
        IF (avg%dump_single) THEN
          ALLOCATE(avg%r4array(1-ng:nx+ng, 1-ng:ny+ng, nspec_local))
          avg%r4array = 0.0_num
        ELSE
          ALLOCATE(avg%array(1-ng:nx+ng, 1-ng:ny+ng, nspec_local))
          avg%array = 0.0_num
        END IF

        avg%species_sum = 0
        avg%n_species = 0
        IF (IAND(mask, c_io_no_sum) == 0) &
            avg%species_sum = averaged_var_dims(io)
        IF (IAND(mask, c_io_species) /= 0) &
            avg%n_species = n_species * averaged_var_dims(io)
        avg%real_time = 0.0_num
        avg%started = .FALSE.
      END IF
    END DO

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies, n

    ALLOCATE(species_list(n_species))
    ALLOCATE(io_list_data(n_species))
    io_list => species_list

    DO ispecies = 1, n_species
      species_list(ispecies)%name = blank
      species_list(ispecies)%mass = -1.0_num
      species_list(ispecies)%charge = 0.0_num
      species_list(ispecies)%weight = 1.0_num
      species_list(ispecies)%dumpmask = c_io_always
      species_list(ispecies)%count = -1
      species_list(ispecies)%id = 0
      species_list(ispecies)%npart_per_cell = -1
      species_list(ispecies)%global_count = 0
      species_list(ispecies)%count_update_step = 0
      species_list(ispecies)%species_type = c_species_id_generic
      species_list(ispecies)%immobile = .FALSE.
      NULLIFY(species_list(ispecies)%next)
      NULLIFY(species_list(ispecies)%prev)
      NULLIFY(species_list(ispecies)%ext_temp_x_min)
      NULLIFY(species_list(ispecies)%ext_temp_x_max)
      NULLIFY(species_list(ispecies)%ext_temp_y_min)
      NULLIFY(species_list(ispecies)%ext_temp_y_max)
      NULLIFY(species_list(ispecies)%secondary_list)
      NULLIFY(species_list(ispecies)%background_density)
      species_list(ispecies)%bc_particle = c_bc_null
    END DO

    DO ispecies = 1, n_species
      CALL initialise_stack(species_list(ispecies)%density_function)
      CALL set_stack_zero  (species_list(ispecies)%density_function)
      DO n = 1, 3
        CALL initialise_stack(species_list(ispecies)%temperature_function(n))
        CALL set_stack_zero  (species_list(ispecies)%temperature_function(n))
        CALL initialise_stack(species_list(ispecies)%drift_function(n))
        CALL set_stack_zero  (species_list(ispecies)%drift_function(n))
        CALL initialise_stack(species_list(ispecies)%dist_fn_range(n))
        CALL set_stack_zero  (species_list(ispecies)%dist_fn_range(n), &
            n_zeros=2)
      END DO
      species_list(ispecies)%fractional_tail_cutoff = 0.0001_num
      species_list(ispecies)%ic_df_type = c_ic_df_thermal
      species_list(ispecies)%electron = .FALSE.
      species_list(ispecies)%ionise = .FALSE.
      species_list(ispecies)%recombine = .FALSE.
      species_list(ispecies)%ionise_to_species = -1
      species_list(ispecies)%release_species = -1
      species_list(ispecies)%n = 0
      species_list(ispecies)%l = 0
      species_list(ispecies)%ionisation_energy = HUGE(0.0_num)
      species_list(ispecies)%migrate%this_species = .FALSE.
      species_list(ispecies)%migrate%fluid = .FALSE.
      species_list(ispecies)%migrate%done = .FALSE.
      species_list(ispecies)%migrate%promoteable = .FALSE.
      species_list(ispecies)%migrate%demoteable = .FALSE.
      species_list(ispecies)%migrate%promote_to_species = 0
      species_list(ispecies)%migrate%demote_to_species = 0
      species_list(ispecies)%migrate%promotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%demotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%promotion_density = HUGE(1.0_num)
      species_list(ispecies)%migrate%demotion_density = 0.0_num
      species_list(ispecies)%fill_ghosts = .TRUE.
#ifndef NO_TRACER_PARTICLES
      species_list(ispecies)%zero_current = .FALSE.
#endif
#ifndef NO_PARTICLE_PROBES
      NULLIFY(species_list(ispecies)%attached_probes)
#endif
    END DO

  END SUBROUTINE setup_species



  SUBROUTINE setup_field_boundaries

    INTEGER :: nx0, nx1, ny0, ny1

    ALLOCATE(ex_x_min(1-ng:ny+ng), ex_x_max(1-ng:ny+ng))
    ALLOCATE(ey_x_min(1-ng:ny+ng), ey_x_max(1-ng:ny+ng))
    ALLOCATE(ez_x_min(1-ng:ny+ng), ez_x_max(1-ng:ny+ng))
    ALLOCATE(bx_x_min(1-ng:ny+ng), bx_x_max(1-ng:ny+ng))
    ALLOCATE(by_x_min(1-ng:ny+ng), by_x_max(1-ng:ny+ng))
    ALLOCATE(bz_x_min(1-ng:ny+ng), bz_x_max(1-ng:ny+ng))

    ALLOCATE(ex_y_min(1-ng:nx+ng), ex_y_max(1-ng:nx+ng))
    ALLOCATE(ey_y_min(1-ng:nx+ng), ey_y_max(1-ng:nx+ng))
    ALLOCATE(ez_y_min(1-ng:nx+ng), ez_y_max(1-ng:nx+ng))
    ALLOCATE(bx_y_min(1-ng:nx+ng), bx_y_max(1-ng:nx+ng))
    ALLOCATE(by_y_min(1-ng:nx+ng), by_y_max(1-ng:nx+ng))
    ALLOCATE(bz_y_min(1-ng:nx+ng), bz_y_max(1-ng:nx+ng))

    nx0 = 1
    nx1 = nx
    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser) nx0 = cpml_x_min_laser_idx-1
    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser) nx1 = cpml_x_max_laser_idx+1

    ex_x_min = 0.5_num * (ex(nx0,:) + ex(nx0-1,:))
    ey_x_min = ey(nx0,:)
    ez_x_min = ez(nx0,:)
    ex_x_max = 0.5_num * (ex(nx1,:) + ex(nx1-1,:))
    ey_x_max = ey(nx1,:)
    ez_x_max = ez(nx1,:)

    bx_x_min = bx(nx0,:)
    by_x_min = 0.5_num * (by(nx0,:) + by(nx0-1,:))
    bz_x_min = 0.5_num * (bz(nx0,:) + bz(nx0-1,:))
    bx_x_max = bx(nx1,:)
    by_x_max = 0.5_num * (by(nx1,:) + by(nx1-1,:))
    bz_x_max = 0.5_num * (bz(nx1,:) + bz(nx1-1,:))

    ny0 = 1
    ny1 = ny
    IF (bc_field(c_bd_y_min) == c_bc_cpml_laser) ny0 = cpml_y_min_laser_idx-1
    IF (bc_field(c_bd_y_max) == c_bc_cpml_laser) ny1 = cpml_y_max_laser_idx+1

    ex_y_min = ex(:,ny0)
    ey_y_min = 0.5_num * (ey(:,ny0) + ey(:,ny0-1))
    ez_y_min = ez(:,ny0)
    ex_y_max = ex(:,ny1)
    ey_y_max = 0.5_num * (ey(:,ny1) + ey(:,ny1-1))
    ez_y_max = ez(:,ny1)

    bx_y_min = 0.5_num * (bx(:,ny0) + bx(:,ny0-1))
    by_y_min = by(:,ny0)
    bz_y_min = 0.5_num * (bz(:,ny0) + bz(:,ny0-1))
    bx_y_max = 0.5_num * (bx(:,ny1) + bx(:,ny1-1))
    by_y_max = by(:,ny1)
    bz_y_max = 0.5_num * (bz(:,ny1) + bz(:,ny1-1))

  END SUBROUTINE setup_field_boundaries



  SUBROUTINE open_files

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode
    LOGICAL :: exists

    IF (rank == 0) THEN
      WRITE(stat_file, '(a, ''/epoch2d.dat'')') TRIM(data_dir)
      IF (ic_from_restart) THEN
        INQUIRE(file=stat_file, exist=exists)
        IF (exists) THEN
          OPEN(unit=stat_unit, status='OLD', position='APPEND', &
              file=stat_file, iostat=errcode)
        ELSE
          OPEN(unit=stat_unit, status='NEW', file=stat_file, iostat=errcode)
        END IF
      ELSE
        OPEN(unit=stat_unit, status='REPLACE', file=stat_file, iostat=errcode)
      END IF
      IF (errcode /= 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot create "epoch2d.dat" output file. The most common ' &
            // 'cause of this problem '
        PRINT*, 'is that the ouput directory does not exist'
        CALL abort_code(c_err_io_error)
        STOP
      END IF
      IF (ic_from_restart) THEN
        WRITE(stat_unit,*)
        WRITE(stat_unit,*) 'Restarting from ', TRIM(restart_filename)
        WRITE(stat_unit,*) ascii_header
      ELSE
        WRITE(stat_unit,*) ascii_header
        WRITE(stat_unit,*)
      END IF
    END IF
#endif

  END SUBROUTINE open_files



  SUBROUTINE flush_stat_file

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode

    IF (rank == 0) THEN
      CLOSE(unit=stat_unit)
      OPEN(unit=stat_unit, status='OLD', position='APPEND', &
          file=stat_file, iostat=errcode)
    END IF
#endif

  END SUBROUTINE flush_stat_file



  SUBROUTINE close_files

#ifdef NO_IO
    RETURN
#else
    IF (rank == 0) CLOSE(unit=stat_unit)
#endif

  END SUBROUTINE close_files



  SUBROUTINE open_status_file

#ifndef NO_IO
    IF (rank == 0) THEN
      OPEN(unit=du, status='OLD', position='APPEND', file=status_filename, &
           iostat=errcode)
    END IF
#endif

  END SUBROUTINE open_status_file



  SUBROUTINE close_status_file

#ifndef NO_IO
    IF (rank == 0) THEN
      CLOSE(du)
    END IF
#endif

  END SUBROUTINE close_status_file



  SUBROUTINE set_initial_values

    INTEGER :: seed

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    ! Set up random number seed
    seed = 7842432
    IF (use_random_seed) CALL SYSTEM_CLOCK(seed)
    seed = seed + rank

    CALL random_init(seed)

  END SUBROUTINE set_initial_values



  SUBROUTINE set_plasma_frequency_dt

    INTEGER :: ispecies, ix, iy
    REAL(num) :: min_dt, omega2, omega, k_max, fac1, fac2, clipped_dens
    TYPE(initial_condition_block), POINTER :: ic

    IF (ic_from_restart) RETURN

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / MIN(dx, dy)

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
        CALL setup_ic_density(ispecies)
        CALL setup_ic_temp(ispecies)

        ic => species_list(ispecies)%initial_conditions

        fac1 = q0**2 / species_list(ispecies)%mass / epsilon0
        fac2 = 3.0_num * k_max**2 * kb / species_list(ispecies)%mass
        IF (ic%density_max > 0) THEN
          DO iy = 1, ny
          DO ix = 1, nx
            clipped_dens = MIN(species_density(ix,iy), ic%density_max)
            omega2 = fac1 * clipped_dens &
                + fac2 * MAXVAL(species_temp(ix,iy,:))
            IF (omega2 <= c_tiny) CYCLE
            omega = SQRT(omega2)
            IF (2.0_num * pi / omega < min_dt) min_dt = 2.0_num * pi / omega
          END DO ! ix
          END DO ! iy
        ELSE
          DO iy = 1, ny
          DO ix = 1, nx
            omega2 = fac1 * species_density(ix,iy) &
                + fac2 * MAXVAL(species_temp(ix,iy,:))
            IF (omega2 <= c_tiny) CYCLE
            omega = SQRT(omega2)
            IF (2.0_num * pi / omega < min_dt) min_dt = 2.0_num * pi / omega
          END DO ! ix
          END DO ! iy
        END IF
      END IF
    END DO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency / 2.0_num

  END SUBROUTINE set_plasma_frequency_dt



  SUBROUTINE set_dt        ! sets CFL limited step

    INTEGER :: io
    REAL(num) :: dt_solver

    CALL set_plasma_frequency_dt
    CALL set_laser_dt

    IF (maxwell_solver == c_maxwell_solver_yee) THEN
      ! Default maxwell solver with field_order = 2, 4 or 6
      ! cfl is a function of field_order
      dt = cfl * dx * dy / SQRT(dx**2 + dy**2) / c

    ELSE
      ! R. Lehe, PhD Thesis (2014)
      ! A. Pukhov, Journal of Plasma Physics 61, 425-433 (1999)
      dt = MIN(dx, dy) / c
    END IF

    IF (any_open) THEN
      dt_solver = dx * dy / SQRT(dx**2 + dy**2) / c
      dt = MIN(dt, dt_solver)
    END IF

    IF (maxwell_solver == c_maxwell_solver_custom) THEN
      dt = dt_custom
      IF (dt_multiplier < 1.0_num) THEN
        IF (rank == 0) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'Custom maxwell solver is used with custom timestep specified'
          PRINT*, 'in input deck. In this case "dt_multiplier" should be set to'
          PRINT*, 'unity in order to ensure the correct time step.'
          PRINT*, 'Overriding dt_multiplier now.'
        END IF
        dt_multiplier = 1.0_num
      END IF
    END IF

    dt_solver = dt

    IF (dt_plasma_frequency > c_tiny) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser > c_tiny) dt = MIN(dt, dt_laser)

    IF (maxwell_solver /= c_maxwell_solver_yee .AND. dt < dt_solver) THEN
      IF (rank == 0) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'Time step "dt_plasma_frequency" or "dt_laser" is smaller than'
        PRINT*, 'time step given by CFL condition, making steps shorter ', &
            'than intended.'
        PRINT*, 'This may have an adverse effect on dispersion properties!'
        PRINT*, 'Increase grid resolution to fix this.'
      END IF
    END IF

    dt = dt_multiplier * dt

    IF (.NOT. any_average) RETURN

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      io_block_list(io)%average_time = MAX(io_block_list(io)%dt_average, &
          dt * io_block_list(io)%nstep_average)

      IF (io_block_list(io)%dt_min_average > 0 &
          .AND. io_block_list(io)%dt_min_average < dt) THEN
        IF (rank == 0) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'Time step is too small to satisfy "nstep_average"'
          PRINT*, 'Averaging will occur over fewer time steps than specified'
          PRINT*, 'Set "dt_multiplier" less than ', &
              dt_multiplier * io_block_list(io)%dt_min_average / dt, &
              ' to fix this'
        END IF
        io_block_list(io)%dt_min_average = -1
      END IF
    END DO

  END SUBROUTINE set_dt



  SUBROUTINE find_species_by_blockid(specname, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: specname
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies, i1, i2

    CALL strip_species_blockid(specname, i1, i2)

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(specname(i1:i2), species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      END IF
    END DO

  END SUBROUTINE find_species_by_blockid



  SUBROUTINE find_species_by_id(species_id, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: species_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(species_id, species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      END IF
    END DO

  END SUBROUTINE find_species_by_id



  SUBROUTINE copy_string(str_in, str_out)

    CHARACTER(LEN=*), INTENT(IN) :: str_in
    CHARACTER(LEN=*), INTENT(INOUT) :: str_out
    INTEGER :: len1, len2, olen, i

    len1 = LEN_TRIM(str_in)
    len2 = LEN_TRIM(str_out)
    olen = MIN(len1,len2)
    IF (olen > 0) THEN
      str_out(1:olen) = str_in(1:olen)
      DO i = olen+1,len2
        str_out(i:i) = ' '
      END DO
    ELSE
      DO i = 1,len2
        str_out(i:i) = ' '
      END DO
    END IF

  END SUBROUTINE copy_string



  SUBROUTINE find_species_by_id_or_blockid(species_id, block_id, species_number)

    CHARACTER(LEN=*), INTENT(INOUT) :: species_id
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: i1, i2

    IF (str_cmp(species_id, '__unknown__') &
        .OR. ICHAR(species_id(1:1)) == 0) THEN
      CALL strip_species_blockid(block_id, i1, i2)
      IF (i2 <= LEN_TRIM(species_id)) THEN
        CALL copy_string(block_id(i1:i2), species_id)
      ELSE
        CALL copy_string('__unknown__', species_id)
      END IF
    END IF

    CALL find_species_by_id(species_id, species_number)

  END SUBROUTINE find_species_by_id_or_blockid



  SUBROUTINE strip_species_blockid(name, i1, i2)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: i1, i2
    INTEGER :: ii

    i1 = 2
    i2 = LEN_TRIM(name)
    DO ii = 1, i2
      IF (name(ii:ii) == '/') RETURN
      i1 = i1 + 1
    END DO
    i1 = 1

  END SUBROUTINE strip_species_blockid



  SUBROUTINE restart_data(step)

    INTEGER, INTENT(OUT) :: step
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1, str2
    CHARACTER(LEN=c_id_length) :: species_id
    CHARACTER(LEN=c_max_string_length) :: name, len_string
    INTEGER :: blocktype, datatype, code_io_version, string_len, ispecies
    INTEGER :: i, is, iblock, nblocks, ndims, geometry
    INTEGER(i8) :: npart, npart_local
    INTEGER(i8), ALLOCATABLE :: nparts(:), npart_locals(:), npart_proc(:)
    INTEGER, DIMENSION(4) :: dims
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    INTEGER, ALLOCATABLE :: random_states_per_proc_old(:)
    REAL(num), DIMENSION(2*c_ndims) :: extents
    LOGICAL :: restart_flag, got_full
    LOGICAL, ALLOCATABLE :: species_found(:)
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(particle_species), POINTER :: species
    INTEGER, POINTER :: species_subtypes(:)
    INTEGER, POINTER :: species_subtypes_i4(:), species_subtypes_i8(:)
    REAL(num) :: offset_x_min, full_x_min, offset_x_max

    got_full = .FALSE.
    npart_global = 0
    step = -1

    full_x_min = 0.0_num
    offset_x_min = 0.0_num

    IF (rank == 0) THEN
      PRINT*,'Attempting to restart from file: ',TRIM(full_restart_filename)
    END IF

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file ', TRIM(full_restart_filename), &
            ' is not a restart dump. Unable to continue.'
      END IF
      CALL abort_code(c_err_io_error)
      STOP
    END IF

    IF (.NOT.str_cmp(code_name, 'Epoch2d')) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch2d. Unable to ', &
            'continue.'
      END IF
      CALL abort_code(c_err_io_error)
      STOP
    END IF

    IF (string_len > c_max_string_length) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(string_len, len_string)
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'constants.F90 to ','be at least ',TRIM(len_string)
      END IF
      CALL abort_code(c_err_io_error)
      STOP
    END IF

    ! Reset io_block parameters
    DO i = 1, n_io_blocks
      IF (io_block_list(i)%dump_first_after_restart) THEN
        io_block_list(i)%dump_first = .TRUE.
      ELSE
        io_block_list(i)%dump_first = .FALSE.
      END IF
      IF (io_block_list(i)%dt_snapshot > 0.0_num) THEN
        io_block_list(i)%time_prev = time
      ELSE
        io_block_list(i)%time_prev = 0.0_num
      END IF
      IF (io_block_list(i)%nstep_snapshot > 0) THEN
        io_block_list(i)%nstep_prev = step
      ELSE
        io_block_list(i)%nstep_prev = 0
      END IF
      io_block_list(i)%walltime_prev = time
      IF (ASSOCIATED(io_block_list(i)%dump_at_nsteps)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_nsteps)
          IF (step >= io_block_list(i)%dump_at_nsteps(is)) THEN
            io_block_list(i)%dump_at_nsteps(is) = HUGE(1)
          END IF
        END DO
      END IF
      IF (ASSOCIATED(io_block_list(i)%dump_at_times)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_times)
          IF (time >= io_block_list(i)%dump_at_times(is)) THEN
            io_block_list(i)%dump_at_times(is) = HUGE(1.0_num)
          END IF
        END DO
      END IF
    END DO

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      CALL create_ascii_header
#ifndef NO_IO
      WRITE(stat_unit,*) ascii_header
      WRITE(stat_unit,*)
#endif
    END IF

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    dt_from_restart = 0.0_num

    IF (rank == 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    ALLOCATE(nparts(n_species))
    ALLOCATE(npart_locals(n_species))
    ALLOCATE(species_found(n_species))
    nparts = 0
    npart_locals = 0
    species_found = .FALSE.

    ! Scan file for particle species and allocate storage
    ! Both particle grids and CPU split blocks are interrogated. CPU split
    ! blocks take precedence if found and "use_exact_restart" is set.
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype == c_blocktype_cpu_split) THEN
        IF (.NOT.use_exact_restart .OR. datatype /= c_datatype_integer8) CYCLE
      ELSE IF (blocktype /= c_blocktype_point_mesh) THEN
        CYCLE
      END IF

      ispecies = 0
      IF (blocktype == c_blocktype_point_mesh) THEN
        CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)
        CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
      ELSE IF (blocktype == c_blocktype_cpu_split) THEN
        IF (block_id(1:3) == 'cpu') THEN
          species_id = block_id(5:LEN(block_id))
          CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
        END IF
      END IF

      IF (ispecies == 0) THEN
        IF (rank == 0) THEN
          IF (.NOT. str_cmp(species_id(1:6), 'subset')) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Particle species "', TRIM(species_id), &
                '" from restart dump ', 'not found in input deck. Ignoring.'
          END IF
        END IF
        CYCLE
      END IF

      IF (species_found(ispecies)) CYCLE

      species => species_list(ispecies)

      IF (blocktype == c_blocktype_cpu_split) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)

        ALLOCATE(npart_proc(dims(1)))
        CALL sdf_read_srl_cpu_split(sdf_handle, npart_proc)

        npart = 0
        DO i = 1,dims(1)
          npart = npart + npart_proc(i)
        END DO
        npart_local = npart_proc(rank+1)
        DEALLOCATE(npart_proc)

        npart_locals(ispecies) = npart_local
        nparts(ispecies) = npart
        species_found(ispecies) = .TRUE.
      ELSE IF (blocktype == c_blocktype_point_mesh) THEN
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        nparts(ispecies) = npart

        npart_local = npart / nproc
        IF (npart_local * nproc /= npart) THEN
          IF (rank < npart - npart_local * nproc) &
              npart_local = npart_local + 1
        END IF

        npart_locals(ispecies) = npart_local
      END IF
    END DO

    DEALLOCATE(species_found)

    ! Do the species allocation
    DO ispecies = 1, n_species
      species => species_list(ispecies)
      npart_local = npart_locals(ispecies)

      CALL create_allocated_partlist(species%attached_list, npart_local)

      npart_global = npart_global + nparts(ispecies)
      species%count = nparts(ispecies)
    END DO

    DEALLOCATE(nparts, npart_locals)

    CALL create_subtypes_for_load(species_subtypes, species_subtypes_i4, &
        species_subtypes_i8)

    CALL sdf_seek_start(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_array)
        IF (use_exact_restart .AND. need_random_state &
            .AND. str_cmp(block_id, 'random_states')) THEN
          CALL sdf_read_array_info(sdf_handle, dims)
          IF (datatype == c_datatype_integer4 .AND. dims(1) == 4*nproc) THEN
            ! Older form of random_states output
            ! Missing the box_muller_cache entry
            ALLOCATE(random_states_per_proc(4*nproc))
            CALL sdf_read_srl(sdf_handle, random_states_per_proc)
            CALL set_random_state(random_states_per_proc(4*rank+1:4*(rank+1)))
            DEALLOCATE(random_states_per_proc)
          ELSE IF (rank == 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Unrecognised format for random_states block in ', &
                    'the restart file. Ignoring.'
          END IF
        ELSE IF (use_exact_restart .AND. need_random_state &
            .AND. str_cmp(block_id, 'random_states_full')) THEN
          CALL sdf_read_array_info(sdf_handle, dims)
          IF (datatype == c_datatype_integer4 .AND. dims(1) == 5*nproc) THEN
            ALLOCATE(random_states_per_proc(4*nproc))
            ALLOCATE(random_states_per_proc_old(5*nproc))
            CALL sdf_read_srl(sdf_handle, random_states_per_proc_old)
            DO i = 0, nproc - 1
              random_states_per_proc(4*i+1:4*(i+1)) = &
                  random_states_per_proc_old(5*i+1:5*(i+1)-1)
            END DO
            DEALLOCATE(random_states_per_proc_old)
            CALL set_random_state(random_states_per_proc(4*rank+1:4*(rank+1)))
            DEALLOCATE(random_states_per_proc)
          ELSE IF (rank == 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Unrecognised format for random_states_full block in ', &
                    'the restart file. Ignoring.'
          END IF
        ELSE IF (str_cmp(block_id, 'file_numbers')) THEN
          CALL sdf_read_array_info(sdf_handle, dims)
          IF (ndims /= 1 .OR. dims(1) /= SIZE(file_numbers)) THEN
            IF (rank == 0) THEN
              PRINT*, '*** WARNING ***'
              PRINT*, 'Output file numbers do not agree. Ignoring.'
            END IF
          ELSE
            CALL sdf_read_srl(sdf_handle, file_numbers)
          END IF
        END IF

        DO i = 1, 2 * c_ndims
          CALL read_laser_phases(sdf_handle, block_id, ndims, i)
          CALL read_injector_depths(sdf_handle, block_id, ndims, i)
        END DO

      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'dt_plasma_frequency')) THEN
          CALL sdf_read_srl(sdf_handle, dt_plasma_frequency)
        ELSE IF (str_cmp(block_id, 'dt')) THEN
          CALL sdf_read_srl(sdf_handle, dt_from_restart)
        ELSE IF (str_cmp(block_id, 'window_shift_fraction')) THEN
          CALL sdf_read_srl(sdf_handle, window_shift_fraction)
        ELSE IF (str_cmp(block_id, 'x_grid_min')) THEN
          got_x_grid_min = .TRUE.
          CALL sdf_read_srl(sdf_handle, x_grid_min_val)
        ELSE IF (str_cmp(block_id, 'y_grid_min')) THEN
          got_y_grid_min = .TRUE.
          CALL sdf_read_srl(sdf_handle, y_grid_min_val)
        ELSE IF (block_id(1:7) == 'weight/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies == 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%weight)
        ELSE IF (block_id(1:5) == 'nppc/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies == 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%npart_per_cell)
        ELSE IF (block_id(1:10) == 'time_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(11:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%time_prev)
              EXIT
            END IF
          END DO
        ELSE IF (block_id(1:11) == 'nstep_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(12:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%nstep_prev)
              EXIT
            END IF
          END DO
        ELSE IF (block_id(1:14) == 'walltime_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(15:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%walltime_prev)
              EXIT
            END IF
          END DO
        ELSE IF (str_cmp(block_id, 'elapsed_time')) THEN
          CALL sdf_read_srl(sdf_handle, old_elapsed_time)
        ELSE IF (block_id(1:18) == 'field_energy_lost/') THEN
          CALL read_field_energy_lost(sdf_handle, block_id)
        END IF
      CASE(c_blocktype_plain_mesh)
        IF (str_cmp(block_id, 'grid') .OR. str_cmp(block_id, 'grid_full')) THEN
          CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)
          IF (.NOT.got_full) THEN
            x_min = extents(1)
            x_max = extents(c_ndims+1)
            y_min = extents(2)
            y_max = extents(c_ndims+2)

            dx = (x_max - x_min) / nx_global
            x_min = x_min + dx * cpml_thickness
            x_max = x_max - dx * cpml_thickness
            dy = (y_max - y_min) / ny_global
            y_min = y_min + dy * cpml_thickness
            y_max = y_max - dy * cpml_thickness

            IF (str_cmp(block_id, 'grid_full')) THEN
              got_full = .TRUE.
              use_offset_grid = .TRUE.
              full_x_min = extents(1)
            ELSE
              ! Offset grid is offset only in x
              offset_x_min = extents(1)
              offset_x_max = extents(c_ndims+1)
            END IF
          END IF
        END IF
      CASE(c_blocktype_point_mesh)
        CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)

        CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
        IF (ispecies == 0) CYClE
        species => species_list(ispecies)

        npart_local = species%attached_list%count
        iterator_list => species%attached_list%head

        CALL sdf_read_point_mesh(sdf_handle, npart_local, &
            species_subtypes(ispecies), it_part)

      CASE(c_blocktype_plain_variable)
        CALL sdf_read_plain_variable_info(sdf_handle, dims, str1, mesh_id)

        IF (.NOT.str_cmp(mesh_id, 'grid')) CYCLE

        IF (dims(1) /= nx_global .OR. dims(2) /= ny_global) THEN
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Number of gridpoints in restart dump does not match', &
                ' the input deck.'
            CALL integer_as_string(nx_global, str1)
            CALL integer_as_string(ny_global, str2)
            PRINT*, 'Input deck grid: ', TRIM(str1), ',', TRIM(str2)
            CALL integer_as_string(dims(1), str1)
            CALL integer_as_string(dims(2), str2)
            PRINT*, 'Restart dump grid: ', TRIM(str1), ',', TRIM(str2)
          END IF
          CALL abort_code(c_err_bad_setup)
          STOP
        END IF

        IF (str_cmp(block_id, 'ex')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ex, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ey')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ey, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ez')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ez, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'by')) THEN
          CALL sdf_read_plain_variable(sdf_handle, by, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'jx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jx, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jy, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jz, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'cpml_psi_eyx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_eyx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_ezx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_ezx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_byx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_byx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bzx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bzx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_exy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_exy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_ezy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_ezy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bxy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bxy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bzy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bzy, &
              subtype_field, subarray_field)

        END IF

      CASE(c_blocktype_point_variable)
        CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
            str1, species_id)

        CALL find_species_by_id_or_blockid(species_id, mesh_id, ispecies)
        IF (ispecies == 0) CYClE
        species => species_list(ispecies)

        IF (npart /= species%count) THEN
          IF (rank == 0) THEN
            CALL integer_as_string(species%count, str1)
            CALL integer_as_string(npart, str2)
            PRINT*, '*** ERROR ***'
            PRINT*, 'Malformed restart dump.'
            PRINT*, 'Particle grid for species "', TRIM(species%name), &
                '" has ', TRIM(str1), ' particles'
            PRINT*, 'but "', TRIM(block_id), '" has ', TRIM(str2)
          END IF
          CALL abort_code(c_err_io_error)
          STOP
        END IF

        iterator_list => species%attached_list%head
        npart_local = species%attached_list%count

        IF (block_id(1:3) == 'px/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_px)

        ELSE IF (block_id(1:3) == 'py/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_py)

        ELSE IF (block_id(1:3) == 'pz/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_pz)

        ELSE IF (block_id(1:3) == 'id/') THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
          ! Particle IDs may either be 4 or 8-byte integers, depending on the
          ! PARTICLE_ID[4] flag used when writing the file. We must read in
          ! the data using the precision written to file and then convert to
          ! the currently used precision.
          IF (datatype == c_datatype_integer8) THEN
            CALL sdf_read_point_variable(sdf_handle, npart_local, &
                species_subtypes_i8(ispecies), it_id8)
          ELSE
            CALL sdf_read_point_variable(sdf_handle, npart_local, &
                species_subtypes_i4(ispecies), it_id4)
          END IF
#else
          IF (rank == 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Discarding particle IDs.'
            PRINT*, 'To use, please recompile with the -DPARTICLE_ID option.'
          END IF
#endif
        ELSE IF (block_id(1:18) == 'persistent_subset/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes_i8(ispecies), it_persistent_subset)

        ELSE IF (block_id(1:7) == 'weight/') THEN
#ifndef PER_SPECIES_WEIGHT
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_weight)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with per particle weight.'
            PRINT*, 'Please recompile without the -DPER_SPECIES_WEIGHT option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:14) == 'optical depth/') THEN
#ifdef PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with optical depths.'
            PRINT*, 'Please recompile with the -DPHOTONS option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:11) == 'qed energy/') THEN
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_qed_energy)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with QED energies.'
            PRINT*, 'Please recompile with -DPHOTONS or -DBREMSSTRAHLUNG.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:14) == 'trident depth/') THEN
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth_trident)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with Trident optical depths.'
            PRINT*, 'Please recompile with the -DTRIDENT_PHOTONS option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:21) == 'bremsstrahlung depth/') THEN
#ifdef BREMSSTRAHLUNG
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth_bremsstrahlung)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with bremsstrahlung optical depths.'
            PRINT*, 'Please recompile with the -DBREMSSTRAHLUNG option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif
        
        ELSE IF (block_id(1:23) == 'time_integrated_work_x/') THEN
#ifdef WORK_DONE_INTEGRATED
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_work_x_total)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with time integrated work.'
            PRINT*, 'Please recompile with the -DWORK_DONE_INTEGRATED option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif
          ELSE IF (block_id(1:23) == 'time_integrated_work_y/') THEN
#ifdef WORK_DONE_INTEGRATED
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_work_y_total)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with time integrated work.'
            PRINT*, 'Please recompile with the -DWORK_DONE_INTEGRATED option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif
          ELSE IF (block_id(1:23) == 'time_integrated_work_z/') THEN
#ifdef WORK_DONE_INTEGRATED
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_work_z_total)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with time integrated work.'
            PRINT*, 'Please recompile with the -DWORK_DONE_INTEGRATED option.'
          END IF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif
        END IF
      END SELECT
    END DO

    CALL sdf_close(sdf_handle)
    CALL free_subtypes_for_load(species_subtypes, species_subtypes_i4, &
        species_subtypes_i8)

    ! Reset dump_at_walltimes
    DO i = 1, n_io_blocks
      IF (ASSOCIATED(io_block_list(i)%dump_at_walltimes)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_walltimes)
          IF (old_elapsed_time >= io_block_list(i)%dump_at_walltimes(is)) THEN
            io_block_list(i)%dump_at_walltimes(is) = HUGE(1.0_num)
          END IF
        END DO
      END IF
    END DO

    IF (use_offset_grid) THEN
      window_offset = full_x_min - offset_x_min
      CALL shift_particles_to_window
    END IF

    CALL setup_grid

    IF (use_offset_grid) THEN
      CALL create_moved_window(offset_x_min)
    END IF

    CALL set_thermal_bcs_all
    CALL setup_persistent_subsets
    CALL setup_background_species

    IF (rank == 0) PRINT*, 'Load from restart dump OK'

  END SUBROUTINE restart_data



  SUBROUTINE setup_persistent_subsets

    INTEGER :: isub
    TYPE(subset), POINTER :: sub

    DO isub = 1, SIZE(subset_list)
      sub => subset_list(isub)
      IF (sub%persistent) THEN
        IF (time < sub%persist_start_time &
            .AND. step < sub%persist_start_step) CYCLE
        sub%locked = .TRUE.
      END IF
    END DO

  END SUBROUTINE setup_persistent_subsets



  SUBROUTINE read_laser_phases(sdf_handle, block_id_in, ndims, boundary)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: block_id_in
    INTEGER, INTENT(IN) :: ndims, boundary
    REAL(num), DIMENSION(:), ALLOCATABLE :: laser_phases
    INTEGER, DIMENSION(4) :: dims
    CHARACTER(LEN=17) :: block_id_compare
    CHARACTER(LEN=5), DIMENSION(6) :: direction_name = &
        (/'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'/)

    block_id_compare = 'laser_' // direction_name(boundary) // '_phase'

    IF (str_cmp(block_id_in, block_id_compare)) THEN
      CALL sdf_read_array_info(sdf_handle, dims)

      IF (ndims /= 1 .OR. dims(1) /= n_lasers(boundary)) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'Number of laser phases on ', TRIM(direction_name(boundary)), &
            ' does not match number of lasers.'
        PRINT*, 'Lasers will be populated in order, but correct operation ', &
            'is not guaranteed'
      END IF

      ALLOCATE(laser_phases(dims(1)))
      CALL sdf_read_srl(sdf_handle, laser_phases)
      CALL setup_laser_phases(laser_phases)
      DEALLOCATE(laser_phases)
    END IF

  END SUBROUTINE read_laser_phases



  SUBROUTINE read_field_energy_lost(sdf_handle, block_id_in)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: block_id_in
    REAL(num) :: energy_lost
    INTEGER :: i
    CHARACTER(LEN=5), DIMENSION(6) :: direction_name = &
        (/'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'/)

    DO i = 1, 2*c_ndims
      IF (str_cmp(block_id_in, 'field_energy_lost/' // direction_name(i))) THEN
        CALL sdf_read_srl(sdf_handle, energy_lost)
        ! The stored value is a global sum, so only keep it on one rank
        IF (rank == 0) field_energy_lost_local(i) = energy_lost
        EXIT
      END IF
    END DO

  END SUBROUTINE read_field_energy_lost



  ! Read injector depths from restart and initialise
  ! Requires the same injectors defined from the deck
  SUBROUTINE read_injector_depths(sdf_handle, block_id_in, ndims, boundary)

    TYPE(sdf_file_handle), INTENT(INOUT) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: block_id_in
    INTEGER, INTENT(IN) :: ndims, boundary
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: depths
    INTEGER :: injector_count
    INTEGER, DIMENSION(4) :: dims
    INTEGER :: n_els, sz, starts
    CHARACTER(LEN=21) :: block_id_compare
    CHARACTER(LEN=5), DIMENSION(6) :: direction_name = &
        (/'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'/)

    block_id_compare = 'injector_' // direction_name(boundary) // '_depths'

    IF (str_cmp(block_id_in, block_id_compare)) THEN
      CALL sdf_read_array_info(sdf_handle, dims)

      ! In 1-d there is one value, 2-d there is one strip (per bnd),
      ! in 3-d one plane etc
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        n_els = ny
        sz = ny_global
        starts = ny_global_min
      ELSE
        n_els = nx
        sz = nx_global
        starts = nx_global_min
      END IF

      ALLOCATE(depths(n_els, dims(c_ndims)))

      CALL sdf_read_array(sdf_handle, depths, (/sz, dims(c_ndims)/), &
          (/starts, 1/), null_proc=(.NOT. is_boundary(boundary)))

      CALL setup_injector_depths(boundary, depths, injector_count)

      ! Got count back so can now check and message
      IF (ndims /= c_ndims .OR. dims(c_ndims) /= injector_count) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'Number of depths on ', TRIM(block_id_in), &
            ' does not match number of injectors.'
        PRINT*, 'Injectors will be populated in order, but correct operation', &
            ' is not guaranteed'
      END IF

      DEALLOCATE(depths)
    END IF

  END SUBROUTINE read_injector_depths



  SUBROUTINE read_cpu_split

    CHARACTER(LEN=c_id_length) :: code_name, block_id
    CHARACTER(LEN=c_max_string_length) :: name
    INTEGER :: step, code_io_version, string_len, nblocks, ndims
    INTEGER :: blocktype, datatype, geometry, iblock, npx, npy
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file ', TRIM(full_restart_filename), &
            ' is not a restart dump. Unable to continue.'
      END IF
      CALL abort_code(c_err_io_error)
      STOP
    END IF

    IF (.NOT.str_cmp(code_name, 'Epoch2d')) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch2d. Unable to ', &
            'continue.'
      END IF
      CALL abort_code(c_err_io_error)
      STOP
    END IF

    nblocks = sdf_read_nblocks(sdf_handle)

    CALL sdf_read_blocklist(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype == c_blocktype_cpu_split &
          .AND. datatype /= c_datatype_integer8) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)
        npx = dims(1) + 1
        npy = dims(2) + 1
        IF (npx * npy == nproc) THEN
          nprocx = npx
          nprocy = npy
          ALLOCATE(old_x_max(nprocx))
          ALLOCATE(old_y_max(nprocy))
          CALL sdf_read_srl_cpu_split(sdf_handle, old_x_max, old_y_max)
        ELSE
          IF (rank == 0 .AND. use_exact_restart_set) THEN
            PRINT*, '*** WARNING ***'
            PRINT'('' SDF restart file was generated using'', &
                & i4,'' CPUs.'')', npx * npy
            PRINT*, 'Ignoring "use_exact_restart" flag.'
          END IF
          use_exact_restart = .FALSE.
        END IF
        EXIT
      END IF
    END DO

    CALL sdf_close(sdf_handle)

  END SUBROUTINE read_cpu_split



  FUNCTION it_part(array, npart_this_it, start, direction, param)

    REAL(num) :: it_part
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param

    INTEGER(i8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur => iterator_list

    DO ipart = 1, npart_this_it
      cur%part_pos(direction) = array(ipart)
      cur => cur%next
    END DO

    it_part = 0

  END FUNCTION it_part



  FUNCTION it_px(array, npart_this_it, start, param)

    REAL(num) :: it_px
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(1) = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_px = 0

  END FUNCTION it_px



  FUNCTION it_py(array, npart_this_it, start, param)

    REAL(num) :: it_py
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(2) = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_py = 0

  END FUNCTION it_py



  FUNCTION it_pz(array, npart_this_it, start, param)

    REAL(num) :: it_pz
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(3) = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_pz = 0

  END FUNCTION it_pz



#ifndef PER_SPECIES_WEIGHT
  FUNCTION it_weight(array, npart_this_it, start, param)

    REAL(num) :: it_weight
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%weight = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_weight = 0

  END FUNCTION it_weight
#endif



#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
  FUNCTION it_id4(array, npart_this_it, start, param)

    USE constants
    INTEGER(i4) :: it_id4
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#if PARTICLE_ID
      iterator_list%id = INT(array(ipart),i8)
#else
      iterator_list%id = array(ipart)
#endif
      iterator_list => iterator_list%next
    END DO

    it_id4 = 0

  END FUNCTION it_id4



  FUNCTION it_id8(array, npart_this_it, start, param)

    USE constants
    INTEGER(i8) :: it_id8
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#if PARTICLE_ID
      iterator_list%id = array(ipart)
#else
      iterator_list%id = INT(array(ipart),i4)
#endif
      iterator_list => iterator_list%next
    END DO

    it_id8 = 0

  END FUNCTION it_id8
#endif



  FUNCTION it_persistent_subset(array, npart_this_it, start, param)

    USE constants
    USE particle_id_hash_mod
    INTEGER(i8) :: it_persistent_subset
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      CALL id_registry%add_with_map(iterator_list, array(ipart))
      iterator_list => iterator_list%next
    END DO

    it_persistent_subset = 0

  END FUNCTION it_persistent_subset



#ifdef PHOTONS
  FUNCTION it_optical_depth(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_optical_depth = 0

  END FUNCTION it_optical_depth



#ifdef TRIDENT_PHOTONS
  FUNCTION it_optical_depth_trident(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth_trident
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth_tri = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_optical_depth_trident = 0

  END FUNCTION it_optical_depth_trident
#endif
#endif



#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
  FUNCTION it_qed_energy(array, npart_this_it, start, param)

    REAL(num) :: it_qed_energy
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%particle_energy = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_qed_energy = 0

  END FUNCTION it_qed_energy
#endif



#ifdef BREMSSTRAHLUNG
  FUNCTION it_optical_depth_bremsstrahlung(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth_bremsstrahlung
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth_bremsstrahlung = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_optical_depth_bremsstrahlung = 0

  END FUNCTION it_optical_depth_bremsstrahlung
#endif



#ifdef WORK_DONE_INTEGRATED
  FUNCTION it_work_x_total(array, npart_this_it, start, param)

    REAL(num) :: it_work_x_total
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%work_x_total = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_work_x_total = 0

  END FUNCTION it_work_x_total



  FUNCTION it_work_y_total(array, npart_this_it, start, param)

    REAL(num) :: it_work_y_total
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%work_y_total = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_work_y_total = 0

  END FUNCTION it_work_y_total


  
  FUNCTION it_work_z_total(array, npart_this_it, start, param)

    REAL(num) :: it_work_z_total
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%work_z_total = array(ipart)
      iterator_list => iterator_list%next
    END DO

    it_work_z_total = 0

  END FUNCTION it_work_z_total  
#endif



  SUBROUTINE shift_particles_to_window

    TYPE(particle), POINTER :: current
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle_species), POINTER :: species
    INTEGER :: ispecies

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      partlist => species%attached_list
      current => partlist%head

      DO WHILE(ASSOCIATED(current))
        current%part_pos(1) = current%part_pos(1) + window_offset

        current => current%next
      END DO
    END DO

  END SUBROUTINE shift_particles_to_window



  SUBROUTINE create_moved_window(x_min)

    REAL(num), INTENT(IN) :: x_min
    INTEGER :: ix

    DO ix = 1 - ng, nx_global + ng
      xb_offset_global(ix) = xb_offset_global(ix) - window_offset
    END DO

  END SUBROUTINE



  SUBROUTINE pre_load_balance

    INTEGER :: npx, npy, ierr
    INTEGER :: old_comm, old_coords(c_ndims)

    IF (.NOT.use_pre_balance .OR. nproc == 1) RETURN

    npx = nprocx
    npy = nprocy
    pre_loading = .TRUE.

    CALL auto_load

    IF (use_optimal_layout) THEN
      CALL get_optimal_layout

      IF (npx == nprocx .AND. npy == nprocy) THEN
        pre_loading = .FALSE.
        RETURN
      END IF

      old_coords(:) = coordinates(:)
      CALL MPI_COMM_DUP(comm, old_comm, ierr)
      CALL setup_communicator
      CALL pre_balance_workload(old_comm, old_coords)
      CALL MPI_COMM_FREE(old_comm, ierr)
      CALL setup_grid
    ELSE
      CALL pre_balance_workload
    END IF

    pre_loading = .FALSE.

  END SUBROUTINE pre_load_balance

END MODULE setup
