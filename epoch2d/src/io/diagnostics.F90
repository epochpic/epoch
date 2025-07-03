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

MODULE diagnostics

  USE calc_df
  USE sdf
  USE deck
  USE dist_fn
  USE evaluator
  USE epoch_source_info
  USE iterators
  USE probes
  USE version_data
  USE setup
  USE deck_io_block
  USE strings
  USE window
  USE timer
  USE antennae
  USE particle_id_hash_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_routines, create_full_timestring
  PUBLIC :: cleanup_stop_files, check_for_stop_condition
  PUBLIC :: deallocate_file_list, count_n_zeros
  PUBLIC :: build_persistent_subsets

  CHARACTER(LEN=*), PARAMETER :: stop_file = 'STOP'
  CHARACTER(LEN=*), PARAMETER :: stop_file_nodump = 'STOP_NODUMP'
  CHARACTER(LEN=*), PARAMETER :: request_dump_file = 'DUMP'

  TYPE(sdf_file_handle) :: sdf_handle
  INTEGER(i8), ALLOCATABLE :: species_offset(:)
  INTEGER(i8), ALLOCATABLE :: ejected_offset(:)
  LOGICAL :: reset_ejected, done_species_offset_init, done_subset_init
  LOGICAL :: restart_flag, dump_source_code, dump_input_decks
  LOGICAL :: dump_field_grid, skipped_any_set
  LOGICAL :: got_stop_restart = .FALSE.
  LOGICAL :: got_request_dump_name = .FALSE.
  LOGICAL :: got_request_dump_restart = .FALSE.
  LOGICAL :: from_dump_request = .FALSE.
  CHARACTER(LEN=string_length) :: request_dump_name = ''
  LOGICAL, ALLOCATABLE :: dump_point_grid(:)
  LOGICAL, ALLOCATABLE, SAVE :: prefix_first_call(:)
  INTEGER :: isubset
  INTEGER, DIMENSION(num_vars_to_dump) :: iomask
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: iodumpmask, dumped_skip_dir
  INTEGER, SAVE :: last_step = -1
  INTEGER, SAVE :: sdf_max_string_length, max_string_length

  ! Data structures for tracking the list of strings written to '.visit' files
  TYPE string_entry
    CHARACTER(LEN=string_length) :: text
    TYPE(string_entry), POINTER :: next
  END TYPE string_entry

  TYPE string_list
    TYPE(string_entry), POINTER :: head, tail
    INTEGER :: count
  END TYPE string_list

  TYPE(string_list), ALLOCATABLE, TARGET :: file_list(:)

  INTERFACE write_particle_variable
    MODULE PROCEDURE &
#if defined(PARTICLE_ID4) || defined(PARTICLE_DEBUG)
        write_particle_variable_i4, &
        write_particle_variable_i8, &
#endif
        write_particle_variable_i8, &
        write_particle_variable_num
  END INTERFACE write_particle_variable

CONTAINS

  SUBROUTINE test_output

    INTEGER :: iprefix
    INTEGER, SAVE :: nstep_prev = -1
    LOGICAL :: force, print_arrays

    force = .FALSE.
    IF (step == nstep_prev) RETURN

    DO iprefix = 1,SIZE(file_prefixes)
      CALL io_test(iprefix, step, print_arrays, force, prefix_first_call)

      IF (.NOT.print_arrays) CYCLE
      nstep_prev = step

      file_numbers(iprefix) = file_numbers(iprefix) + 1
    END DO

  END SUBROUTINE test_output



  SUBROUTINE count_n_zeros

    INTEGER :: i, step_orig, n_dumps
    REAL(num) :: time_orig
    INTEGER, ALLOCATABLE :: file_numbers_orig(:)
    TYPE(io_block_type), ALLOCATABLE :: io_block_orig(:)
    INTEGER :: ndt
    INTEGER(i8) :: istep, step_interval
    REAL(num) :: time_start, time0, time1, dt_interval
    REAL(num), PARAMETER :: total_time = 30.0_num

#ifdef NO_IO
    RETURN
#endif

    IF (n_io_blocks <= 0) RETURN
    IF (.NOT.use_accurate_n_zeros) RETURN

    ALLOCATE(file_list(n_io_blocks+2))
    ALLOCATE(prefix_first_call(SIZE(file_prefixes)))
    ALLOCATE(file_numbers_orig(SIZE(file_prefixes)))
    ALLOCATE(io_block_orig(n_io_blocks))

    file_numbers_orig(:) = file_numbers(:)
    prefix_first_call = .TRUE.
    DO i = 1,n_io_blocks+2
      file_list(i)%count = 0
    END DO

    DO i = 1,n_io_blocks
      CALL copy_io_block(io_block_list(i), io_block_orig(i))
    END DO

    step_orig = step
    time_orig = time

    time0 = MPI_WTIME()
    time_start = time0
    istep = 0
    step_interval = 100
    ndt = 10

    IF (.NOT.ic_from_restart) CALL test_output

    DO
      step = step + 1
      time = time + dt / 2.0_num

      IF ((step >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end)) EXIT

      CALL test_output
      time = time + dt / 2.0_num

      istep = istep + 1
      IF (istep == step_interval) THEN
        time1 = MPI_WTIME()
        dt_interval = (total_time + time_start - time1) / ndt
        ndt = ndt - 1
        step_interval = INT(step_interval * dt_interval / (time1 - time0), i8)
        IF (step_interval < 0 .OR. time1 - time_start >= total_time) THEN
          EXIT
        END IF
        time0 = time1
        istep = 0
      END IF
    END DO

    CALL test_output

    n_dumps = 0
    step = step_orig
    time = time_orig
    DO i = 1,n_io_blocks
      CALL copy_io_block(io_block_orig(i), io_block_list(i))
    END DO
    DO i = 1,SIZE(file_prefixes)
      IF (file_numbers(i) > n_dumps) n_dumps = file_numbers(i)
      file_numbers(i) = file_numbers_orig(i)
    END DO

    DEALLOCATE(file_list)
    DEALLOCATE(prefix_first_call)
    DEALLOCATE(io_block_orig)
    DEALLOCATE(file_numbers_orig)

    IF (n_dumps > 1) THEN
      n_dumps = n_dumps - 1
      n_zeros = MAX(n_zeros, FLOOR(LOG10(REAL(n_dumps))) + 1)
    END IF

  END SUBROUTINE count_n_zeros



  SUBROUTINE output_routines(step, force_write)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    LOGICAL, INTENT(IN), OPTIONAL :: force_write
    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=c_max_path_length) :: full_filename
    CHARACTER(LEN=c_max_string_length) :: dump_type, temp_name
    CHARACTER(LEN=c_id_length) :: temp_block_id
    REAL(num) :: eta_time, dr, r0
    REAL(num), DIMENSION(:), ALLOCATABLE :: x_reduced, y_reduced
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: array
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    INTEGER :: code, i, io, ispecies, iprefix, mask, rn, dir, dumped, nval, n
    INTEGER :: errcode
    INTEGER :: random_state(4)
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER, SAVE :: nstep_prev = -1
    LOGICAL :: convert, force, any_written, restart_id, print_arrays
    LOGICAL, SAVE :: first_call = .TRUE.
    TYPE(particle_species), POINTER :: species
    TYPE(subset), POINTER :: sub
    CHARACTER(LEN=15) :: timestring, eta_timestring
    CHARACTER(LEN=1), DIMENSION(3) :: dim_tags = (/'x', 'y', 'z'/)
    CHARACTER(LEN=5), DIMENSION(6) :: dir_tags = &
        (/'x_max', 'y_max', 'z_max', 'x_min', 'y_min', 'z_min'/)
    INTEGER, DIMENSION(6) :: fluxdir = &
        (/c_dir_x, c_dir_y, c_dir_z, -c_dir_x, -c_dir_y, -c_dir_z/)

    ! Clean-up any cached RNG state
    CALL random_flush_cache

#ifdef NO_IO
    RETURN
#endif

    CALL build_persistent_subsets

    timer_walltime = -1.0_num
    IF (step /= last_step) THEN
      last_step = step
      IF (rank == 0 .AND. stdout_frequency > 0) THEN
        IF (MOD(step, stdout_frequency) == 0) THEN
          timer_walltime = MPI_WTIME()
          elapsed_time = timer_walltime - walltime_started

          IF (reset_walltime) THEN
            CALL create_timestring(elapsed_time, timestring)
            elapsed_time = elapsed_time + old_elapsed_time
          ELSE
            elapsed_time = elapsed_time + old_elapsed_time
            CALL create_timestring(elapsed_time, timestring)
          END IF

          IF (print_eta_string) THEN
            eta_timestring = ''
            IF (time > c_tiny) THEN
              ! If nsteps is specified and the limit to runtime, use it to
              ! calculate ETA. Otherwise use t_end.
              IF (nsteps /= -1 .AND. nsteps * dt < t_end) THEN
                eta_time = (nsteps * dt - time) * elapsed_time / time
              ELSE
                eta_time = (t_end - time) * elapsed_time / time
              END IF
              CALL create_timestring(eta_time, eta_timestring)
            END IF
            WRITE(*, '(''Time'', g14.6, '', iteration'', i9, '' after'', &
                & a, '', ETA'',a)') time, step, timestring, eta_timestring
          ELSE
            WRITE(*, '(''Time'', g20.12, '' and iteration'', i12, '' after'', &
                & a)') time, step, timestring
          END IF
          IF (skipped_any_set) THEN
            WRITE(*,*) 'One or more subset ranges were empty: ', &
                  'their fields were not output.'
            skipped_any_set = .FALSE.
          END IF
        END IF
      END IF
    END IF

    IF (n_io_blocks <= 0) RETURN

    force = .FALSE.
    IF (PRESENT(force_write)) force = force_write

    IF (step == nstep_prev .AND. .NOT.force) RETURN

    IF (first_call) THEN
      ALLOCATE(dumped_skip_dir(c_ndims,n_subsets))
      ALLOCATE(file_list(n_io_blocks+2))
      ALLOCATE(prefix_first_call(SIZE(file_prefixes)))
      prefix_first_call = first_call
      DO i = 1,n_io_blocks+2
        file_list(i)%count = 0
      END DO
      first_call = .FALSE.
      ! Setting a large output buffer for point data can often make
      ! output much faster.
      ! The default value is set in deck_io_global_block.F90
      CALL sdf_set_point_array_size(sdf_buffer_size)
      sdf_max_string_length = sdf_get_max_string_length()
      max_string_length = MIN(sdf_max_string_length, c_max_string_length)
    END IF

    dims = (/nx_global, ny_global/)

    reset_ejected = .FALSE.
    any_written = .FALSE.

    IF (n_species > 0) THEN
      ALLOCATE(dump_point_grid(n_species))
    ELSE
      ALLOCATE(dump_point_grid(1))
    END IF

    DO iprefix = 1,SIZE(file_prefixes)
      CALL io_test(iprefix, step, print_arrays, force, prefix_first_call)

      IF (.NOT.print_arrays) CYCLE

      IF (.NOT.any_written) THEN
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        CALL create_subtypes
        any_written = .TRUE.
        IF (timer_collect) THEN
          IF (timer_walltime < 0.0_num) THEN
            CALL timer_start(c_timer_io)
          ELSE
            CALL timer_start(c_timer_io, .TRUE.)
          END IF
        END IF
      END IF

      ! Increase n_zeros if needed

      io = file_numbers(iprefix)
      rn = 1
      nval = 1
      DO i = 1, 1000
        IF (rn > io) THEN
          nval = i - 1
          EXIT
        END IF
        rn = rn * 10
      END DO

      IF (nval > n_zeros) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'n_zeros increased to enable further output'
        END IF
        n_zeros = nval
      END IF

      ! Allows a maximum of 10^999 output dumps, should be enough for anyone
      ! (feel free to laugh when this isn't the case)
      WRITE(filename_fmt, '(''(a, i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
          n_zeros, n_zeros
      WRITE(filename, filename_fmt) TRIM(file_prefixes(iprefix)), &
          file_numbers(iprefix)
      full_filename = TRIM(filesystem) &
          // TRIM(data_dir) // '/' // TRIM(filename)

      ! Always dump the variables with the 'Every' attribute
      code = c_io_always

      ! Only dump variables with the 'FULL' attributre on full dump intervals
      IF (MOD(file_numbers(iprefix), full_dump_every) == 0 &
          .AND. full_dump_every > -1) code = IOR(code, c_io_full)
      IF (restart_flag) code = IOR(code, c_io_restartable)
      dump_field_grid = .FALSE.

      nstep_prev = step

      DO isubset = 1, n_subsets
        errcode = 0
        sub => subset_list(isubset)
        IF (.NOT. sub%time_varying) CYCLE
        DO n = 1, c_subset_max
          IF (sub%use_restriction_function(n)) THEN
            sub%restriction(n) = evaluate(sub%restriction_function(n), errcode)
          END IF
        END DO
        IF (sub%space_restrictions) CALL create_subset_subtypes(isubset)
      END DO

      ! open the file
      CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_write)
      CALL sdf_set_string_length(sdf_handle, c_max_string_length)
      CALL sdf_write_header(sdf_handle, 'Epoch2d', 1, step, time, &
          restart_flag, jobid)
      CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
          c_commit_id, epoch_bytes_checksum, c_compile_machine, &
          c_compile_flags, defines, c_compile_date, run_date)
      CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
          cell_x_max, cell_y_max)

      timer_walltime = MPI_WTIME()
      elapsed_time = old_elapsed_time + timer_walltime - walltime_started
      CALL sdf_write_srl(sdf_handle, 'elapsed_time', 'Wall-time', elapsed_time)

      file_numbers(iprefix) = file_numbers(iprefix) + 1

      IF (restart_flag) THEN
        CALL sdf_write_srl(sdf_handle, 'dt', 'Time increment', dt)
        CALL sdf_write_srl(sdf_handle, 'dt_plasma_frequency', &
            'Plasma frequency timestep restriction', dt_plasma_frequency)
        IF (move_window .AND. window_started) THEN
          CALL sdf_write_srl(sdf_handle, 'window_shift_fraction', &
              'Window Shift Fraction', window_shift_fraction)
        END IF
        CALL sdf_write_srl(sdf_handle, 'x_grid_min', &
            'Minimum grid position', x_grid_min)

        DO i = 1, 2 * c_ndims
          CALL write_laser_phases(sdf_handle, i)
          CALL write_injector_depths(sdf_handle, i)
        END DO

        CALL write_antenna_phases(sdf_handle)

        CALL write_id_starts(sdf_handle)

        DO io = 1, n_io_blocks
          CALL sdf_write_srl(sdf_handle, &
              'time_prev/'//TRIM(io_block_list(io)%name), &
              'time_prev/'//TRIM(io_block_list(io)%name), &
              io_block_list(io)%time_prev)
          CALL sdf_write_srl(sdf_handle, &
              'walltime_prev/'//TRIM(io_block_list(io)%name), &
              'walltime_prev/'//TRIM(io_block_list(io)%name), &
              io_block_list(io)%walltime_prev)
          CALL sdf_write_srl(sdf_handle, &
              'nstep_prev/'//TRIM(io_block_list(io)%name), &
              'nstep_prev/'//TRIM(io_block_list(io)%name), &
              io_block_list(io)%nstep_prev)
        END DO

        DO ispecies = 1, n_species
          species => io_list(ispecies)
          CALL sdf_write_srl(sdf_handle, 'nppc/' // TRIM(species%name), &
              'Particles/Particles Per Cell/' // TRIM(species%name), &
              species%npart_per_cell)
        END DO

        CALL sdf_write_srl(sdf_handle, 'file_prefixes', &
            'Output File Stem Names', file_prefixes)

        CALL sdf_write_srl(sdf_handle, 'file_numbers', &
            'Output File Sequence Numbers', file_numbers)

        IF (need_random_state) THEN
          CALL get_random_state(random_state)
          ALLOCATE(random_states_per_proc(4*nproc))
          CALL MPI_GATHER(random_state, 4, MPI_INTEGER, &
              random_states_per_proc, 4, MPI_INTEGER, 0, comm, errcode)
          CALL sdf_write_srl(sdf_handle, 'random_states', &
              'Random States', random_states_per_proc)
          DEALLOCATE(random_states_per_proc)
        END IF
      END IF

      iomask = iodumpmask(1,:)

      CALL write_field(c_dump_ex, code, 'ex', 'Electric Field/Ex', 'V/m', &
          c_stagger_ex, ex)
      CALL write_field(c_dump_ey, code, 'ey', 'Electric Field/Ey', 'V/m', &
          c_stagger_ey, ey)
      CALL write_field(c_dump_ez, code, 'ez', 'Electric Field/Ez', 'V/m', &
          c_stagger_ez, ez)

      CALL write_field(c_dump_bx, code, 'bx', 'Magnetic Field/Bx', 'T', &
          c_stagger_bx, bx)
      CALL write_field(c_dump_by, code, 'by', 'Magnetic Field/By', 'T', &
          c_stagger_by, by)
      CALL write_field(c_dump_bz, code, 'bz', 'Magnetic Field/Bz', 'T', &
          c_stagger_bz, bz)

      CALL write_field(c_dump_jx, code, 'jx', 'Current/Jx', 'A/m^2', &
          c_stagger_jx, jx)
      CALL write_field(c_dump_jy, code, 'jy', 'Current/Jy', 'A/m^2', &
          c_stagger_jy, jy)
      CALL write_field(c_dump_jz, code, 'jz', 'Current/Jz', 'A/m^2', &
          c_stagger_jz, jz)

      IF (cpml_boundaries) THEN
        CALL sdf_write_srl(sdf_handle, 'boundary_thickness', &
            'Boundary thickness', cpml_thickness)

        CALL write_field(c_dump_cpml_psi_eyx, code, 'cpml_psi_eyx', &
            'CPML/Ey_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_eyx)
        CALL write_field(c_dump_cpml_psi_ezx, code, 'cpml_psi_ezx', &
            'CPML/Ez_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_ezx)
        CALL write_field(c_dump_cpml_psi_byx, code, 'cpml_psi_byx', &
            'CPML/By_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_byx)
        CALL write_field(c_dump_cpml_psi_bzx, code, 'cpml_psi_bzx', &
            'CPML/Bz_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_bzx)

        CALL write_field(c_dump_cpml_psi_exy, code, 'cpml_psi_exy', &
            'CPML/Ex_y', 'A/m^2', c_stagger_cell_centre, cpml_psi_exy)
        CALL write_field(c_dump_cpml_psi_ezy, code, 'cpml_psi_ezy', &
            'CPML/Ez_y', 'A/m^2', c_stagger_cell_centre, cpml_psi_ezy)
        CALL write_field(c_dump_cpml_psi_bxy, code, 'cpml_psi_bxy', &
            'CPML/Bx_y', 'A/m^2', c_stagger_cell_centre, cpml_psi_bxy)
        CALL write_field(c_dump_cpml_psi_bzy, code, 'cpml_psi_bzy', &
            'CPML/Bz_y', 'A/m^2', c_stagger_cell_centre, cpml_psi_bzy)
      END IF

      IF (n_subsets > 0) THEN
        DO i = 1, n_species
          CALL create_empty_partlist(io_list_data(i)%attached_list)
        END DO
      END IF

      DO isubset = 1, n_subsets + 1
        done_species_offset_init = .FALSE.
        done_subset_init = .FALSE.
        dump_point_grid = .FALSE.
        IF (isubset > 1) io_list => io_list_data
        iomask = iodumpmask(isubset,:)

#ifndef PER_SPECIES_WEIGHT
        CALL write_particle_variable(c_dump_part_weight, code, 'Weight', '', &
            it_output_real)
#else
        IF (IAND(iomask(c_dump_part_weight), code) /= 0) THEN
          CALL build_species_subset

          DO ispecies = 1, n_species
            species => io_list(ispecies)
            IF (IAND(species%dumpmask, code) /= 0 &
                .OR. IAND(code, c_io_restartable) /= 0) THEN
              CALL sdf_write_srl(sdf_handle, 'weight/' // TRIM(species%name), &
                  'Particles/Weight/' // TRIM(species%name), species%weight)
            END IF
          END DO
        END IF
#endif

#ifdef PER_PARTICLE_CHARGE_MASS
        CALL write_particle_variable(c_dump_part_charge, code, &
            'Q', 'C', it_output_real)
        CALL write_particle_variable(c_dump_part_mass, code, &
            'Mass', 'kg', it_output_real)
#else
        IF (IAND(iomask(c_dump_part_charge), code) /= 0) THEN
          CALL build_species_subset

          DO ispecies = 1, n_species
            species => io_list(ispecies)
            IF (IAND(species%dumpmask, code) /= 0 &
                .OR. IAND(code, c_io_restartable) /= 0) THEN
              CALL sdf_write_srl(sdf_handle, 'charge/' // TRIM(species%name), &
                  'Particles/Charge/' // TRIM(species%name), species%charge)
            END IF
          END DO
        END IF

        IF (IAND(iomask(c_dump_part_mass), code) /= 0) THEN
          CALL build_species_subset

          DO ispecies = 1, n_species
            species => io_list(ispecies)
            IF (IAND(species%dumpmask, code) /= 0 &
                .OR. IAND(code, c_io_restartable) /= 0) THEN
              CALL sdf_write_srl(sdf_handle, 'mass/' // TRIM(species%name), &
                  'Particles/Mass/' // TRIM(species%name), species%mass)
            END IF
          END DO
        END IF
#endif

        mask = iomask(c_dump_total_energy_sum)
        IF (IAND(mask, code) /= 0) THEN
          CALL build_species_subset

          IF (IAND(mask, c_io_species) == 0) THEN
            CALL calc_total_energy_sum(.FALSE.)
          ELSE
            CALL calc_total_energy_sum(.TRUE.)

            DO ispecies = 1, n_species
              species => io_list(ispecies)
              IF (IAND(species%dumpmask, code) == 0) CYCLE

              CALL sdf_write_srl(sdf_handle, &
                  'total_particle_energy/' // TRIM(species%name), &
                  'Total Particle Energy/' // TRIM(species%name) // ' (J)', &
                  total_particle_energy_species(ispecies))
            END DO
          END IF

          IF (isubset == 1) THEN
            IF (IAND(mask, c_io_no_sum) == 0) THEN
              CALL sdf_write_srl(sdf_handle, 'total_particle_energy', &
                  'Total Particle Energy in Simulation (J)', &
                  total_particle_energy)
            END IF
            CALL sdf_write_srl(sdf_handle, 'total_field_energy', &
                'Total Field Energy in Simulation (J)', total_field_energy)
          END IF
        END IF

        CALL write_particle_variable(c_dump_part_px, code, &
            'Px', 'kg.m/s', it_output_real)
        CALL write_particle_variable(c_dump_part_py, code, &
            'Py', 'kg.m/s', it_output_real)
        CALL write_particle_variable(c_dump_part_pz, code, &
            'Pz', 'kg.m/s', it_output_real)

        CALL write_particle_variable(c_dump_part_vx, code, &
            'Vx', 'm/s', it_output_real)
        CALL write_particle_variable(c_dump_part_vy, code, &
            'Vy', 'm/s', it_output_real)
        CALL write_particle_variable(c_dump_part_vz, code, &
            'Vz', 'm/s', it_output_real)

        CALL write_particle_variable(c_dump_part_ek, code, &
            'Ek', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_rel_mass, code, &
            'Relativistic Mass', 'kg', it_output_real)
        CALL write_particle_variable(c_dump_part_gamma, code, &
            'Gamma', '', it_output_real)
#ifdef PARTICLE_DEBUG
        CALL write_particle_variable(c_dump_part_proc, code, &
            'Processor', '', it_output_integer4)
        CALL write_particle_variable(c_dump_part_proc0, code, &
            'Processor_at_t0', '', it_output_integer4)
#endif
#if defined(PARTICLE_ID)
        CALL write_particle_variable(c_dump_part_id, code, &
            'ID', '#', it_output_integer8)
#elif defined(PARTICLE_ID4)
        CALL write_particle_variable(c_dump_part_id, code, &
            'ID', '#', it_output_integer4)
#endif
        IF (id_registry%get_hash_count(step) > 0) THEN
          CALL write_particle_variable(c_dump_persistent_ids, code, &
              'persistent_subset', '#', it_output_integer8)
        END IF
#ifdef PHOTONS
        CALL write_particle_variable(c_dump_part_opdepth, code, &
            'Optical depth', '', it_output_real)
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
        CALL write_particle_variable(c_dump_part_qed_energy, code, &
            'QED energy', 'J', it_output_real)
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
        CALL write_particle_variable(c_dump_part_opdepth_tri, code, &
            'Trident Depth', '', it_output_real)
#endif
#ifdef BREMSSTRAHLUNG
        CALL write_particle_variable(c_dump_part_opdepth_brem, code, &
            'Bremsstrahlung Depth', '', it_output_real)
#endif
#ifdef WORK_DONE_INTEGRATED
        CALL write_particle_variable(c_dump_part_work_x, code, &
            'Work_x_direction', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_work_y, code, &
            'Work_y_direction', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_work_z, code, &
            'Work_z_direction', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_work_x_total, code, &
            'Time_Integrated_Work_x', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_work_y_total, code, &
            'Time_Integrated_Work_y', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_work_z_total, code, &
            'Time_Integrated_Work_z', 'J', it_output_real)
#endif
        CALL write_particle_grid(code)

        ! These are derived variables from the particles
        CALL write_nspecies_field(c_dump_ekbar, code, &
            'ekbar', 'Average_Particle_Energy', 'J', &
            c_stagger_cell_centre, calc_ekbar, array)

        CALL write_nspecies_field(c_dump_mass_density, code, &
            'mass_density', 'Mass_Density', 'kg/m^3', &
            c_stagger_cell_centre, calc_mass_density, array)

        CALL write_nspecies_field(c_dump_charge_density, code, &
            'charge_density', 'Charge_Density', 'C/m^3', &
            c_stagger_cell_centre, calc_charge_density, array)

        CALL write_nspecies_field(c_dump_number_density, code, &
            'number_density', 'Number_Density', '1/m^3', &
            c_stagger_cell_centre, calc_number_density, array)

        CALL write_nspecies_field(c_dump_ppc, code, &
            'ppc', 'Particles_Per_Cell', 'n_particles', &
            c_stagger_cell_centre, calc_ppc, array)

        CALL write_nspecies_field(c_dump_average_weight, code, &
            'average_weight', 'Particles_Average_Weight', 'weight', &
            c_stagger_cell_centre, calc_average_weight, array)

        CALL write_nspecies_field(c_dump_average_px, code, &
            'average_px', 'Particles_Average_Px', 'kg.m/s', &
            c_stagger_cell_centre, calc_average_momentum, array, (/c_dir_x/))

        CALL write_nspecies_field(c_dump_average_py, code, &
            'average_py', 'Particles_Average_Py', 'kg.m/s', &
            c_stagger_cell_centre, calc_average_momentum, array, (/c_dir_y/))

        CALL write_nspecies_field(c_dump_average_pz, code, &
            'average_pz', 'Particles_Average_Pz', 'kg.m/s', &
            c_stagger_cell_centre, calc_average_momentum, array, (/c_dir_z/))

        CALL write_nspecies_field(c_dump_temperature, code, &
            'temperature', 'Temperature', 'K', &
            c_stagger_cell_centre, calc_temperature, array)

        CALL write_nspecies_field(c_dump_temperature_x, code, &
            'temperature_x', 'Temperature_x', 'K', &
            c_stagger_cell_centre, calc_temperature, array, (/c_dir_x/))

        CALL write_nspecies_field(c_dump_temperature_y, code, &
            'temperature_y', 'Temperature_y', 'K', &
            c_stagger_cell_centre, calc_temperature, array, (/c_dir_y/))

        CALL write_nspecies_field(c_dump_temperature_z, code, &
            'temperature_z', 'Temperature_z', 'K', &
            c_stagger_cell_centre, calc_temperature, array, (/c_dir_z/))

        CALL write_nspecies_field(c_dump_jx, code, &
            'jx', 'Jx', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_current, array, (/c_dir_x/))

        CALL write_nspecies_field(c_dump_jy, code, &
            'jy', 'Jy', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_current, array, (/c_dir_y/))

        CALL write_nspecies_field(c_dump_jz, code, &
            'jz', 'Jz', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_current, array, (/c_dir_z/))

        CALL write_nspecies_field(c_dump_ekflux, code, &
            'ekflux', 'Particle_Energy_Flux', 'W/m^2', &
            c_stagger_cell_centre, calc_ekflux, array, fluxdir, dir_tags)

        CALL write_nspecies_field(c_dump_poynt_flux, code, &
            'poynt_flux', 'Poynting Flux', 'W/m^2', &
            c_stagger_cell_centre, calc_poynt_flux, array, fluxdir(1:3), &
            dim_tags)

        CALL write_nspecies_field(c_dump_cou_log, code, &
            'cou_log', 'Coulomb Logarithm', 'dimensionless', &
            c_stagger_cell_centre, calc_cou_log, array)

        IF (isubset /= 1) THEN
          DO i = 1, n_species
            CALL append_partlist(species_list(i)%attached_list, &
                io_list(i)%attached_list)
          END DO
          DO i = 1, n_species
            CALL create_empty_partlist(io_list(i)%attached_list)
          END DO
        END IF
      END DO

      io_list => species_list
      iomask = iodumpmask(1,:)
      mask = iomask(c_dump_grid)
      restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
      convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

      ! Write the cartesian mesh
      IF (IAND(mask, code) /= 0) dump_field_grid = .TRUE.
      IF (IAND(mask, c_io_never) /= 0) dump_field_grid = .FALSE.
      IF (restart_flag) dump_field_grid = .TRUE.

      use_offset_grid = .FALSE.
      DO io = 1, n_io_blocks
        use_offset_grid = use_offset_grid &
           .OR. (io_block_list(io)%dump .AND. io_block_list(io)%use_offset_grid)
      END DO

      IF (dump_field_grid) THEN
        IF (.NOT. use_offset_grid) THEN
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
              xb_global(1:nx_global+1), yb_global(1:ny_global+1), convert)
        ELSE
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
              xb_offset_global(1:nx_global+1), &
              yb_offset_global(1:ny_global+1), convert)
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid_full', &
              'Grid/Grid_Full', xb_global(1:nx_global+1), &
              yb_global(1:ny_global+1), convert)
        END IF
      END IF

      dumped_skip_dir = 0
      dumped = 1

      DO io = 1, n_subsets
        sub => subset_list(io)
        IF (.NOT.sub%dump_field_grid) CYCLE

        IF (.NOT. sub%skip) THEN
          temp_block_id = 'grid/' // TRIM(sub%name)
          temp_name = 'Grid/' // TRIM(sub%name)

          CALL check_name_length('subset', &
              'Grid/' // TRIM(sub%name))
          ranges = cell_global_ranges(sub)

          IF (.NOT. use_offset_grid) THEN
            CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), xb_global(ranges(1,1):ranges(2,1)), &
                yb_global(ranges(1,2):ranges(2,2)), convert)
          ELSE
            CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), xb_offset_global(ranges(1,1):ranges(2,1)), &
                yb_offset_global(ranges(1,2):ranges(2,2)), convert)
          END IF
        ELSE
          DO i = 1, io - 1
            dumped = dumped + SUM(dumped_skip_dir(:,i) - sub%skip_dir)
          END DO
          IF (dumped == 0) CYCLE
          dumped = 0

          dumped_skip_dir(:,io) = sub%skip_dir

          dir = 1
          rn = sub%n_global(dir) + 1
          ALLOCATE(x_reduced(rn))
          dr = sub%skip_dir(dir) * dx
          i = sub%n_start(dir) + 1
          r0 = xb_global(i) + 0.5_num * (dx - dr)

          DO i = 1, rn
            x_reduced(i) = r0 + (i - 1) * dr
          END DO

          dir = 2
          rn = sub%n_global(dir) + 1
          ALLOCATE(y_reduced(rn))
          dr = sub%skip_dir(dir) * dy
          i = sub%n_start(dir) + 1
          r0 = yb_global(i) + 0.5_num * (dy - dr)

          DO i = 1, rn
            y_reduced(i) = r0 + (i - 1) * dr
          END DO

          IF (.NOT. use_offset_grid) THEN
            temp_block_id = 'grid/r_' // TRIM(sub%name)
            temp_name = 'Grid/Reduced_' // TRIM(sub%name)

            CALL check_name_length('subset', &
                'Grid/Reduced_' // TRIM(sub%name))

            CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), x_reduced, y_reduced, convert)
          ELSE
            temp_block_id = 'grid_full/r_' // TRIM(sub%name)
            temp_name = 'Grid_Full/Reduced_' // TRIM(sub%name)

            CALL check_name_length('subset', &
                'Grid_Full/Reduced_' // TRIM(sub%name))

            CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), x_reduced, y_reduced, convert)

            temp_block_id = 'grid/r_' // TRIM(sub%name)
            temp_name = 'Grid/Reduced_' // TRIM(sub%name)

            CALL check_name_length('subset', &
                'Grid/Reduced_' // TRIM(sub%name))

            dir = 1
            rn = sub%n_global(dir) + 1
            dr = sub%skip_dir(dir) * dx
            i = sub%n_start(dir) + 1
            r0 = xb_offset_global(i) + 0.5_num * (dx - dr)

            DO i = 1, rn
              x_reduced(i) = r0 + (i - 1) * dr
            END DO

            dir = 2
            rn = sub%n_global(dir) + 1
            dr = sub%skip_dir(dir) * dy
            i = sub%n_start(dir) + 1
            r0 = yb_offset_global(i) + 0.5_num * (dy - dr)

            DO i = 1, rn
              y_reduced(i) = r0 + (i - 1) * dr
            END DO

            CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), x_reduced, y_reduced, convert)
          END IF

          DEALLOCATE(x_reduced, y_reduced)
        END IF
        sub%dump_field_grid = .FALSE.
      END DO

      IF (IAND(iomask(c_dump_dist_fns), code) /= 0) THEN
        CALL write_dist_fns(sdf_handle, code, iomask(c_dump_dist_fns))
      END IF

#ifndef NO_PARTICLE_PROBES
      IF (IAND(iomask(c_dump_probes), code) /= 0 .OR. restart_flag) THEN
        CALL write_probes(sdf_handle, code, iomask(c_dump_probes))
      END IF
#endif

      IF (dump_input_decks) CALL write_input_decks(sdf_handle)
      IF (dump_source_code .AND. epoch_bytes_len > 0) &
          CALL write_source_info(sdf_handle)

      IF (IAND(iomask(c_dump_absorption), code) /= 0) THEN
        CALL MPI_ALLREDUCE(laser_absorb_local, laser_absorbed, 1, mpireal, &
            MPI_SUM, comm, errcode)
        CALL MPI_ALLREDUCE(laser_inject_local, laser_injected, 1, mpireal, &
            MPI_SUM, comm, errcode)
        IF (laser_injected > 0.0_num) THEN
          laser_absorbed = laser_absorbed / laser_injected
        ELSE
          laser_absorbed = 0.0_num
        END IF
        CALL sdf_write_srl(sdf_handle, 'laser_enTotal', &
            'Absorption/Total Laser Energy Injected (J)', laser_injected)
        CALL sdf_write_srl(sdf_handle, 'abs_frac', &
            'Absorption/Fraction of Laser Energy Absorbed (%)', laser_absorbed)
      END IF

      ! close the file
      CALL sdf_close(sdf_handle)

      IF (rank == 0) THEN
        DO io = 1, n_io_blocks
          IF (io_block_list(io)%dump) THEN
            dump_type = TRIM(io_block_list(io)%name)
            CALL append_filename(dump_type, filename, io)
          END IF
        END DO
        IF (IAND(code, c_io_restartable) /= 0) THEN
          dump_type = 'restart'
          CALL append_filename(dump_type, filename, n_io_blocks+1)
        END IF
        IF (IAND(code, c_io_full) /= 0) THEN
          dump_type = 'full'
          CALL append_filename(dump_type, filename, n_io_blocks+2)
        END IF
        IF (iprefix > 1) dump_type = TRIM(file_prefixes(iprefix))
        IF (from_dump_request) THEN
          WRITE(stat_unit,'(a)') 'Writing DUMP file request'
          from_dump_request = .FALSE.
        END IF
        WRITE(stat_unit, '(''Wrote '', a7, '', '', a18, '' at time'', &
          & g12.4, '' and iteration'', i8)') dump_type, &
          TRIM(filename), time, step
        CALL flush_stat_file()
      END IF

      IF (force) EXIT
    END DO

    DEALLOCATE(dump_point_grid)

    IF (.NOT.any_written) RETURN

    DEALLOCATE(array)
    IF (ALLOCATED(species_offset))  DEALLOCATE(species_offset)
    IF (ALLOCATED(ejected_offset))  DEALLOCATE(ejected_offset)
    CALL free_subtypes()

    IF (reset_ejected) THEN
      DO i = 1, n_species
        CALL destroy_partlist(ejected_list(i)%attached_list)
      END DO
    END IF

    IF (timer_collect) CALL timer_stop(c_timer_io)

  END SUBROUTINE output_routines



  SUBROUTINE write_laser_phases(sdf_handle, boundary)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(:), ALLOCATABLE :: laser_phases
    INTEGER :: ilas
    TYPE(laser_block), POINTER :: current_laser
    CHARACTER(LEN=17) :: block_name
    CHARACTER(LEN=5), DIMENSION(6) :: direction_name = &
        (/'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'/)

    IF (n_lasers(boundary) < 1) RETURN

    block_name = 'laser_' // direction_name(boundary) // '_phase'

    ALLOCATE(laser_phases(n_lasers(boundary)))
    ilas = 1
    current_laser => lasers

    DO WHILE(ASSOCIATED(current_laser))
      IF (current_laser%boundary == boundary) THEN
        laser_phases(ilas) = current_laser%current_integral_phase
        ilas = ilas + 1
      END IF
      current_laser => current_laser%next
    END DO

    CALL sdf_write_srl(sdf_handle, TRIM(block_name), TRIM(block_name), &
        n_lasers(boundary), laser_phases, 0)
    DEALLOCATE(laser_phases)

  END SUBROUTINE write_laser_phases



  SUBROUTINE write_injector_depths(sdf_handle, boundary)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), POINTER :: injector
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: depths
    INTEGER :: inj
    INTEGER :: n_els, sz, starts
    CHARACTER(LEN=21) :: block_name
    CHARACTER(LEN=5), DIMENSION(6) :: direction_name = &
        (/'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'/)

    inj = 0
    injector => injector_list

    DO WHILE(ASSOCIATED(injector))
      IF (injector%boundary == boundary) THEN
        inj = inj + 1
      END IF
      injector => injector%next
    END DO

    IF (inj == 0) RETURN

    block_name = 'injector_' // direction_name(boundary) // '_depths'

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      n_els = ny
      sz = ny_global
      starts = ny_global_min
    ELSE
      n_els = nx
      sz = nx_global
      starts = nx_global_min
    END IF

    ALLOCATE(depths(n_els, inj))

    inj = 0
    injector => injector_list

    DO WHILE(ASSOCIATED(injector))
      IF (injector%boundary == boundary) THEN
        inj = inj + 1
        depths(:,inj) = injector%depth(1:n_els)
      END IF
      injector => injector%next
    END DO

    CALL sdf_write_array(sdf_handle, TRIM(block_name), TRIM(block_name), &
        depths, (/sz, inj/), (/starts, 1/), &
        null_proc=(.NOT. is_boundary(boundary)))

    DEALLOCATE(depths)

  END SUBROUTINE write_injector_depths



  SUBROUTINE write_antenna_phases(sdf_handle)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
    REAL(num), DIMENSION(:), ALLOCATABLE :: antenna_phases
    INTEGER :: iant

    IF (any_antennae) THEN
      ALLOCATE(antenna_phases(SIZE(antenna_list)))
      DO iant = 1, SIZE(antenna_list)
        antenna_phases(iant) = antenna_list(iant)%contents%phase_history
      END DO
      CALL sdf_write_srl(sdf_handle, 'antenna_phases', 'antenna_phases', &
          SIZE(antenna_list), antenna_phases, 0)
      DEALLOCATE(antenna_phases)
    END IF

  END SUBROUTINE write_antenna_phases



  SUBROUTINE write_id_starts(sdf_handle)

    TYPE(sdf_file_handle), INTENT(IN) :: sdf_handle
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: id_starts
    INTEGER(i8) :: start_local
    INTEGER :: ierr

    ALLOCATE(id_starts(1:nproc+1))

    start_local = INT(highest_id, i8)
    CALL MPI_GATHER(start_local, 1, MPI_INTEGER8, id_starts(2), 1, &
        MPI_INTEGER8, 0, comm, ierr)

    id_starts(1) = n_cpu_bits

    CALL sdf_write_srl(sdf_handle, "id_starts", "id_starts", id_starts, 0)

    DEALLOCATE(id_starts)
#endif

  END SUBROUTINE write_id_starts



  SUBROUTINE check_name_length(shorten, string)

    CHARACTER(LEN=*), INTENT(IN) :: shorten, string
    CHARACTER(LEN=c_max_string_length) :: len_string
    INTEGER :: length

    length = LEN_TRIM(string)

    IF (length > max_string_length) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(length, len_string)
        PRINT*, '*** WARNING ***'
        PRINT*, 'Output block name ', TRIM(string), ' is truncated.'
        IF (length > c_max_string_length &
            .AND. length > sdf_max_string_length) THEN
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'in both EPOCH and the ', &
              'SDF library to at least ', TRIM(len_string)
        ELSE IF (length > sdf_max_string_length) THEN
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'in the ', &
              'SDF library to at least ', TRIM(len_string)
        ELSE
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'to at least ', &
              TRIM(len_string)
        END IF
      END IF
    END IF

  END SUBROUTINE check_name_length



  SUBROUTINE append_filename(listname, filename, list_index)
    ! This routine updates a list of each output type (eg. full, dump, normal)
    ! that can be passed to VisIt to filter the output.

    CHARACTER(LEN=*), INTENT(IN) :: listname, filename
    INTEGER, INTENT(IN) :: list_index
    CHARACTER(LEN=c_max_path_length) :: listfile
    TYPE(string_list), POINTER :: list
    TYPE(string_entry), POINTER :: lcur
    INTEGER :: ierr, i
    LOGICAL :: exists

    list => file_list(list_index)

    listfile = TRIM(data_dir) // '/' // TRIM(listname) // '.visit'
    INQUIRE(file=listfile, exist=exists)

    IF (list%count == 0) THEN
      ALLOCATE(list%head)
      list%tail => list%head
      IF (ic_from_restart .AND. exists) CALL setup_file_list(listfile, list)
    ELSE
      ALLOCATE(list%tail%next)
      list%tail => list%tail%next
    END IF

    IF (list%count > 0) THEN
      lcur  => list%head
      IF (TRIM(lcur%text) == TRIM(filename)) RETURN
      DO i = 2,list%count
        lcur  => lcur%next
        IF (TRIM(lcur%text) == TRIM(filename)) RETURN
      END DO
    END IF

    list%count = list%count + 1
    list%tail%text = TRIM(filename)

    IF (exists) THEN
      OPEN(unit=lu, status='OLD', position='APPEND', file=listfile, iostat=ierr)
    ELSE
      OPEN(unit=lu, status='NEW', file=listfile, iostat=errcode)
    END IF

    WRITE(lu,'(a)') TRIM(filename)
    CLOSE(lu)

  END SUBROUTINE append_filename



  SUBROUTINE setup_file_list(listfile, list)
    ! This routine initialises the file list contained in the given
    ! *.visit file which is required when restarting.

    CHARACTER(LEN=*), INTENT(IN) :: listfile
    TYPE(string_list), POINTER :: list
    CHARACTER(LEN=string_length) :: str
    INTEGER :: ierr

    OPEN(unit=lu, status='OLD', file=listfile, iostat=ierr)
    DO
      READ(lu, '(A)', iostat=ierr) str
      IF (ierr < 0) EXIT
      list%tail%text = TRIM(str)
      list%count = list%count + 1
      ALLOCATE(list%tail%next)
      list%tail => list%tail%next
    END DO
    CLOSE(lu)

  END SUBROUTINE setup_file_list



  SUBROUTINE deallocate_file_list

    INTEGER :: i, n, nlist, stat
    TYPE(string_entry), POINTER :: current, next

    IF (ALLOCATED(file_list)) THEN
      DO i = 1, n_io_blocks+2
        nlist = file_list(i)%count
        IF (nlist > 0) THEN
          current => file_list(i)%head
          DO n = 1, nlist
            next => current%next
            IF (ASSOCIATED(current)) DEALLOCATE(current, STAT=stat)
            current => next
          END DO
        END IF
      END DO
      DEALLOCATE(file_list, STAT=stat)
    END IF
    DEALLOCATE(iodumpmask, STAT=stat)
    DEALLOCATE(dumped_skip_dir, STAT=stat)
    DEALLOCATE(prefix_first_call, STAT=stat)

  END SUBROUTINE deallocate_file_list



  SUBROUTINE io_test(iprefix, step, print_arrays, force, first_call)

    INTEGER, INTENT(IN) :: iprefix, step
    LOGICAL, INTENT(OUT) :: print_arrays
    LOGICAL, INTENT(IN) :: force
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: first_call
    INTEGER :: id, io, is, nstep_next = 0, av_block
    REAL(num) :: t0, t1, time_first
    LOGICAL :: last_call, dump, done_time_sync

    IF (.NOT.ALLOCATED(iodumpmask)) &
        ALLOCATE(iodumpmask(n_subsets+1,num_vars_to_dump))

    restart_flag = .FALSE.
    dump_source_code = .FALSE.
    dump_input_decks = .FALSE.
    print_arrays = .FALSE.
    from_dump_request = .FALSE.
    iomask = c_io_none
    iodumpmask = c_io_none

    IF (time >= t_end .OR. step == nsteps) THEN
      last_call = .TRUE.
    ELSE
      last_call = .FALSE.
    END IF

    done_time_sync = .FALSE.

    DO io = 1, n_io_blocks
      io_block_list(io)%dump = .FALSE.

      IF (io_block_list(io)%disabled) CYCLE
      IF (io_block_list(io)%prefix_index /= iprefix) CYCLE

      IF (last_call .AND. io_block_list(io)%dump_last) &
          io_block_list(io)%dump = .TRUE.
      IF (first_call(iprefix) .AND. io_block_list(io)%dump_first) &
          io_block_list(io)%dump = .TRUE.

      IF (force) THEN
        io_block_list(io)%dump = .TRUE.
        restart_flag = .TRUE.
      END IF

      IF (got_request_dump_name) THEN
        IF (str_cmp(request_dump_name, io_block_list(io)%name)) THEN
          io_block_list(io)%dump = .TRUE.
          from_dump_request = .TRUE.
          got_request_dump_name = .FALSE.
        END IF
      END IF

      IF (.NOT. done_time_sync) THEN
        IF (walltime_start > 0.0_num .OR. walltime_stop < HUGE(walltime_stop) &
            .OR. io_block_list(io)%walltime_start > 0.0_num &
            .OR. io_block_list(io)%walltime_stop < HUGE(walltime_stop) &
            .OR. ALLOCATED(io_block_list(io)%dump_at_walltimes)) THEN
          done_time_sync = .TRUE.
          CALL MPI_BCAST(elapsed_time, 1, mpireal, 0, comm, errcode)
        END IF
      END IF

      IF (elapsed_time < walltime_start) CYCLE
      IF (elapsed_time > walltime_stop)  CYCLE
      IF (elapsed_time < io_block_list(io)%walltime_start) CYCLE
      IF (elapsed_time > io_block_list(io)%walltime_stop)  CYCLE

      t0 = io_block_list(io)%walltime_interval
      IF (t0 > 0.0_num) THEN
        IF (elapsed_time - io_block_list(io)%walltime_prev >= t0) THEN
          io_block_list(io)%dump = .TRUE.
          io_block_list(io)%walltime_prev = elapsed_time
        END IF
      END IF

      IF (ALLOCATED(io_block_list(io)%dump_at_nsteps)) THEN
        DO is = 1, SIZE(io_block_list(io)%dump_at_nsteps)
          IF (step >= io_block_list(io)%dump_at_nsteps(is)) THEN
            io_block_list(io)%dump = .TRUE.
            io_block_list(io)%dump_at_nsteps(is) = HUGE(1)
          END IF
        END DO
      END IF

      IF (ALLOCATED(io_block_list(io)%dump_at_times)) THEN
        DO is = 1, SIZE(io_block_list(io)%dump_at_times)
          IF (time >= io_block_list(io)%dump_at_times(is)) THEN
            io_block_list(io)%dump = .TRUE.
            io_block_list(io)%dump_at_times(is) = HUGE(1.0_num)
          END IF
        END DO
      END IF

      IF (ALLOCATED(io_block_list(io)%dump_at_walltimes)) THEN
        DO is = 1, SIZE(io_block_list(io)%dump_at_walltimes)
          IF (elapsed_time >= io_block_list(io)%dump_at_walltimes(is)) THEN
            io_block_list(io)%dump = .TRUE.
            io_block_list(io)%dump_at_walltimes(is) = HUGE(1.0_num)
          END IF
        END DO
      END IF

      ! Work out the time that the next dump will occur based on the
      ! current timestep
      t0 = HUGE(1.0_num)
      t1 = HUGE(1.0_num)
      IF (io_block_list(io)%dt_snapshot >= 0.0_num) &
          t0 = io_block_list(io)%time_prev + io_block_list(io)%dt_snapshot
      IF (io_block_list(io)%nstep_snapshot >= 0) THEN
        nstep_next = io_block_list(io)%nstep_prev &
            + io_block_list(io)%nstep_snapshot
        t1 = time + dt * (nstep_next - step)
      END IF

      IF (t0 < t1) THEN
        ! Next I/O dump based on dt_snapshot
        time_first = t0
        IF (io_block_list(io)%dt_snapshot > 0 .AND. time >= t0) THEN
          ! Store the most recent output time that qualifies
          DO
            t0 = io_block_list(io)%time_prev + io_block_list(io)%dt_snapshot
            IF (t0 > time) EXIT
            io_block_list(io)%time_prev = t0
          END DO
          dump = .TRUE.
          IF (dump .AND. time < io_block_list(io)%time_start)  dump = .FALSE.
          IF (dump .AND. time > io_block_list(io)%time_stop)   dump = .FALSE.
          IF (dump .AND. step < io_block_list(io)%nstep_start) dump = .FALSE.
          IF (dump .AND. step > io_block_list(io)%nstep_stop)  dump = .FALSE.
          IF (dump .AND. time < time_start)  dump = .FALSE.
          IF (dump .AND. time > time_stop)   dump = .FALSE.
          IF (dump .AND. step < nstep_start) dump = .FALSE.
          IF (dump .AND. step > nstep_stop)  dump = .FALSE.
          IF (dump) io_block_list(io)%dump = .TRUE.
        END IF
      ELSE
        ! Next I/O dump based on nstep_snapshot
        time_first = t1
        IF (io_block_list(io)%nstep_snapshot > 0 &
            .AND. step >= nstep_next) THEN
          ! Store the most recent output step that qualifies
          DO
            nstep_next = io_block_list(io)%nstep_prev &
                + io_block_list(io)%nstep_snapshot
            IF (nstep_next > step) EXIT
            io_block_list(io)%nstep_prev = nstep_next
          END DO
          dump = .TRUE.
          IF (dump .AND. time < io_block_list(io)%time_start)  dump = .FALSE.
          IF (dump .AND. time > io_block_list(io)%time_stop)   dump = .FALSE.
          IF (dump .AND. step < io_block_list(io)%nstep_start) dump = .FALSE.
          IF (dump .AND. step > io_block_list(io)%nstep_stop)  dump = .FALSE.
          IF (dump .AND. time < time_start)  dump = .FALSE.
          IF (dump .AND. time > time_stop)   dump = .FALSE.
          IF (dump .AND. step < nstep_start) dump = .FALSE.
          IF (dump .AND. step > nstep_stop)  dump = .FALSE.
          IF (dump) io_block_list(io)%dump = .TRUE.
        END IF
      END IF

      io_block_list(io)%average_time_start = &
          time_first - io_block_list(io)%average_time
    END DO

    DO io = 1, n_io_blocks
      IF (io_block_list(io)%dump) THEN
        print_arrays = .TRUE.
        IF (io_block_list(io)%restart) restart_flag = .TRUE.
        IF (io_block_list(io)%dump_source_code) dump_source_code = .TRUE.
        IF (io_block_list(io)%dump_input_decks) dump_input_decks = .TRUE.
        IF (file_numbers(iprefix) > io_block_list(io)%dump_cycle) &
            file_numbers(iprefix) = io_block_list(io)%dump_cycle_first_index
        iomask = IOR(iomask, io_block_list(io)%dumpmask)
        IF (n_subsets /= 0) THEN
          DO is = 1, n_subsets
            iodumpmask(1+is,:) = &
                IOR(iodumpmask(1+is,:), subset_list(is)%dumpmask(io,:))
          END DO
        END IF
      END IF
    END DO

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      IF (time >= io_block_list(io)%average_time_start) THEN
        DO id = 1, num_vars_to_dump
          av_block = averaged_var_block(id)
          IF (IAND(io_block_list(io)%dumpmask(id), c_io_averaged) /= 0) THEN
            CALL average_field(id, io_block_list(av_block)%averaged_data(id))
          END IF
        END DO
      END IF
    END DO

    IF (got_request_dump_restart .OR. got_stop_restart) THEN
      got_request_dump_restart = .FALSE.
      got_stop_restart = .FALSE.
      restart_flag = .TRUE.
      print_arrays = .TRUE.
      dump_source_code = .TRUE.
      dump_input_decks = .TRUE.
      iomask = IOR(iomask, io_block_list(1)%dumpmask)
    END IF

    IF (MOD(file_numbers(1), restart_dump_every) == 0 &
        .AND. restart_dump_every > -1) restart_flag = .TRUE.
    IF (first_call(iprefix) .AND. force_first_to_be_restartable) &
        restart_flag = .TRUE.
    IF (last_call .AND. force_final_to_be_restartable) restart_flag = .TRUE.
    IF (force) THEN
      restart_flag = .TRUE.
      print_arrays = .TRUE.
    END IF

    IF (.NOT.restart_flag .AND. .NOT.new_style_io_block) THEN
      dump_source_code = .FALSE.
      dump_input_decks = .FALSE.
    END IF

    IF (first_call(iprefix)) first_call(iprefix) = .FALSE.

    IF (force) iomask = IOR(iomask, io_block_list(1)%dumpmask)
    iodumpmask(1,:) = iomask

  END SUBROUTINE io_test



  SUBROUTINE average_field(ioutput, avg)

    INTEGER, INTENT(IN) :: ioutput
    TYPE(averaged_data_block) :: avg
    INTEGER :: n_species_local, ispecies, idir, is, nd
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: array

    avg%real_time = avg%real_time + dt
    avg%started = .TRUE.

    n_species_local = avg%n_species + avg%species_sum

    IF (n_species_local <= 0) RETURN

    IF (avg%dump_single) THEN
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(ex * dt, r4)
      CASE(c_dump_ey)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(ey * dt, r4)
      CASE(c_dump_ez)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(ez * dt, r4)
      CASE(c_dump_bx)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(bx * dt, r4)
      CASE(c_dump_by)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(by * dt, r4)
      CASE(c_dump_bz)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) + REAL(bz * dt, r4)
      CASE(c_dump_jx)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) &
            + REAL(jx(1-ng:nx+ng,1-ng:ny+ng) * dt, r4)
      CASE(c_dump_jy)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) &
            + REAL(jy(1-ng:nx+ng,1-ng:ny+ng) * dt, r4)
      CASE(c_dump_jz)
        avg%r4array(:,:,1) = avg%r4array(:,:,1) &
            + REAL(jz(1-ng:nx+ng,1-ng:ny+ng) * dt, r4)
      CASE(c_dump_ekbar)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_ppc)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_ppc(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_weight)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_weight(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_px)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_x)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_py)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_y)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_pz)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_z)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_x)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_x)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_y)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_y)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_z)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_z)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_ekflux)
        nd = averaged_var_dims(ioutput)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          is = (ispecies - 1) / nd - avg%species_sum / nd + 1
          idir = ispecies - avg%species_sum - (is - 1) * nd
          CALL calc_ekflux(array, is, idir)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_poynt_flux)
        nd = averaged_var_dims(ioutput)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_poynt_flux(array, 0, ispecies)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      CASE(c_dump_cou_log)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_cou_log(array, ispecies-avg%species_sum)
          avg%r4array(:,:,ispecies) = avg%r4array(:,:,ispecies) &
              + REAL(array * dt, r4)
        END DO
        DEALLOCATE(array)
      END SELECT
    ELSE
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%array(:,:,1) = avg%array(:,:,1) + ex * dt
      CASE(c_dump_ey)
        avg%array(:,:,1) = avg%array(:,:,1) + ey * dt
      CASE(c_dump_ez)
        avg%array(:,:,1) = avg%array(:,:,1) + ez * dt
      CASE(c_dump_bx)
        avg%array(:,:,1) = avg%array(:,:,1) + bx * dt
      CASE(c_dump_by)
        avg%array(:,:,1) = avg%array(:,:,1) + by * dt
      CASE(c_dump_bz)
        avg%array(:,:,1) = avg%array(:,:,1) + bz * dt
      CASE(c_dump_jx)
        avg%array(:,:,1) = avg%array(:,:,1) + jx(1-ng:nx+ng,1-ng:ny+ng) * dt
      CASE(c_dump_jy)
        avg%array(:,:,1) = avg%array(:,:,1) + jy(1-ng:nx+ng,1-ng:ny+ng) * dt
      CASE(c_dump_jz)
        avg%array(:,:,1) = avg%array(:,:,1) + jz(1-ng:nx+ng,1-ng:ny+ng) * dt
      CASE(c_dump_ekbar)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_ppc)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_ppc(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_weight)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_weight(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_px)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_x)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_py)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_y)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_average_pz)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_average_momentum(array, ispecies-avg%species_sum, c_dir_z)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_x)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_x)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_y)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_y)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_temperature_z)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum, c_dir_z)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_ekflux)
        nd = averaged_var_dims(ioutput)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          is = (ispecies - 1) / nd - avg%species_sum / nd + 1
          idir = ispecies - avg%species_sum - (is - 1) * nd
          CALL calc_ekflux(array, is, idir)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_poynt_flux)
        nd = averaged_var_dims(ioutput)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_poynt_flux(array, 0, ispecies)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      CASE(c_dump_cou_log)
        ALLOCATE(array(1-ng:nx+ng,1-ng:ny+ng))
        DO ispecies = 1, n_species_local
          CALL calc_cou_log(array, ispecies-avg%species_sum)
          avg%array(:,:,ispecies) = avg%array(:,:,ispecies) + array * dt
        END DO
        DEALLOCATE(array)
      END SELECT
    END IF

  END SUBROUTINE average_field



  SUBROUTINE write_field(id, code, block_id, name, units, stagger, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: reduced
    INTEGER :: io, mask, dumped
    INTEGER :: i, ii, rnx, j, jj, rny
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    INTEGER, DIMENSION(c_ndims) :: dims
    LOGICAL :: convert, dump_skipped, restart_id, normal_id, unaveraged_id
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    TYPE(averaged_data_block), POINTER :: avg
    TYPE(io_block_type), POINTER :: iob
    TYPE(subset), POINTER :: sub

    mask = iomask(id)

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0
    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    ! The variable is either averaged or has snapshot specified
    unaveraged_id = IAND(mask, c_io_averaged) == 0 &
        .OR. IAND(mask, c_io_snapshot) /= 0

    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (.NOT.restart_id .AND. IAND(mask, c_io_never) /= 0) RETURN

    dims = (/nx_global, ny_global/)

    IF (convert) THEN
      subtype  = subtype_field_r4
      IF (id == c_dump_jx .OR. id == c_dump_jy .OR. id == c_dump_jz) THEN
        subarray = subarray_field_big_r4
      ELSE
        subarray = subarray_field_r4
      END IF
    ELSE
      subtype  = subtype_field
      IF (id == c_dump_jx .OR. id == c_dump_jy .OR. id == c_dump_jz) THEN
        subarray = subarray_field_big
      ELSE
        subarray = subarray_field
      END IF
    END IF

    ! Output unaveraged data if:
    !  1. This is a restart dump and a restart variable
    !  2. This is a normal dump and the variable
    !      a) is not averaged
    !      b) is averaged and also has snapshot specified

    dumped_skip_dir = 0
    dumped = 1
    dump_skipped = .FALSE.

    DO io = 1, n_subsets
      IF (IAND(iodumpmask(io+1,id), code) == 0) CYCLE

      sub => subset_list(io)
      IF (.NOT. (sub%skip .OR. sub%space_restrictions)) CYCLE

      IF (.NOT. sub%skip) THEN
        ! Output every subset. Trust user not to do parts twice
        ! Skip empty subsets
        DO i = 1, c_ndims
          IF (sub%n_global(i) <= 0) THEN
            skipped_any_set = .TRUE.
            CYCLE
          END IF
        END DO

        IF (convert) THEN
          rsubtype  = sub%subtype_r4
          rsubarray = sub%subarray_r4
        ELSE
          rsubtype  = sub%subtype
          rsubarray = sub%subarray
        END IF
        temp_grid_id = 'grid/' // TRIM(sub%name)
        CALL check_name_length('subset', TRIM(name) &
            // '/Core_' // TRIM(sub%name))

        temp_block_id = TRIM(block_id)// '/c_' // TRIM(sub%name)
        temp_name = TRIM(name) // '/Core_' // TRIM(sub%name)

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), array, rsubtype, rsubarray, convert)
        sub%dump_field_grid = .TRUE.

      ELSE
        ! This should prevent a reduced variable from being dumped multiple
        ! times in the same output file
        DO i = 1, io - 1
          dumped = dumped + SUM(dumped_skip_dir(:,i) - sub%skip_dir)
        END DO
        IF (dumped == 0) CYCLE
        dumped = 0

        dumped_skip_dir(:,io) = sub%skip_dir

        rnx = sub%n_local(1)
        rny = sub%n_local(2)

        ALLOCATE(reduced(rnx,rny))

        jj = sub%n_start(2) + 1
        DO j = 1, rny
          ii = sub%n_start(1) + 1
          DO i = 1, rnx
            reduced(i,j) = array(ii,jj)
            ii = ii + sub%skip_dir(1)
          END DO
          jj = jj + sub%skip_dir(2)
        END DO

        IF (convert) THEN
          rsubtype  = sub%subtype_r4
          rsubarray = sub%subarray_r4
        ELSE
          rsubtype  = sub%subtype
          rsubarray = sub%subarray
        END IF

        temp_grid_id = 'grid/r_' // TRIM(sub%name)
        temp_block_id = TRIM(block_id) // '/r_' // TRIM(sub%name)
        temp_name = TRIM(name) // '/Reduced_' // TRIM(sub%name)

        CALL check_name_length('subset', &
            TRIM(name) // '/Reduced_' // TRIM(sub%name))

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

        dump_skipped = .TRUE.
        sub%dump_field_grid = .TRUE.
        DEALLOCATE(reduced)
      END IF
    END DO

    IF (IAND(mask, code) == 0) RETURN

    IF (restart_id .OR. (.NOT.dump_skipped .AND. unaveraged_id)) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(block_id), &
          TRIM(name), TRIM(units), dims, stagger, 'grid', array, &
          subtype, subarray, convert)
      dump_field_grid = .TRUE.
    END IF

    ! Dump averages
    DO io = 1, n_io_blocks
      iob => io_block_list(io)
      IF (.NOT.iob%dump) CYCLE

      IF (IAND(mask, c_io_averaged) == 0) CYCLE

      avg => iob%averaged_data(id)
      IF (.NOT.avg%started) CYCLE

      IF (avg%dump_single) THEN
        avg%r4array = avg%r4array / REAL(avg%real_time, r4)

        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
            TRIM(units), dims, stagger, 'grid', &
            avg%r4array(:,:,1), subtype_field_r4, subarray_field_r4)

        avg%r4array = 0.0_num
      ELSE
        avg%array = avg%array / avg%real_time

        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
            TRIM(units), dims, stagger, 'grid', &
            avg%array(:,:,1), subtype_field, subarray_field)

        avg%array = 0.0_num
      END IF

      dump_field_grid = .TRUE.
      avg%real_time = 0.0_num
      avg%started = .FALSE.
    END DO

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN), OPTIONAL :: dir_tags
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: reduced
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies, io, mask, idir, ndirs, iav
    INTEGER :: i, ii, rnx, j, jj, rny
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    LOGICAL :: convert, dump_sum, dump_species, dump_skipped
    LOGICAL :: normal_id, restart_id, unaveraged_id, dump_part
    TYPE(averaged_data_block), POINTER :: avg
    TYPE(io_block_type), POINTER :: iob
    TYPE(subset), POINTER :: sub

    INTERFACE
      SUBROUTINE func(data_array, current_species, direction)
        USE constants
        REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
        INTEGER, INTENT(IN), OPTIONAL :: direction
      END SUBROUTINE func
    END INTERFACE

    mask = iomask(id)
    IF (IAND(mask, code) == 0) RETURN
    IF (IAND(mask, c_io_never) /= 0) RETURN

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0

    IF (.NOT.normal_id) RETURN

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    ! The variable is either averaged or has snapshot specified
    unaveraged_id = IAND(mask, c_io_averaged) == 0 &
        .OR. IAND(mask, c_io_snapshot) /= 0

    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    END IF

    IF (PRESENT(fluxdir)) THEN
      ndirs = SIZE(fluxdir)
    ELSE
      ndirs = 1
    END IF
    dims = (/nx_global, ny_global/)

    dump_sum = unaveraged_id &
        .AND. IAND(mask, c_io_no_sum) == 0 .AND. IAND(mask, c_io_field) == 0
    dump_species = unaveraged_id .AND. IAND(mask, c_io_species) /= 0

    IF (isubset == 1) THEN
      dump_skipped = .FALSE.
      dump_part = .FALSE.
    ELSE
      sub => subset_list(isubset-1)
      dump_skipped = sub%skip
      dump_part = sub%space_restrictions
      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      END IF
    END IF

    IF (dump_sum .OR. dump_species) THEN
      CALL build_species_subset
      ! Skip empty subsets
      IF (dump_part) THEN
        DO i = 1, c_ndims
          IF (sub%n_global(i) <= 0) THEN
            skipped_any_set = .TRUE.
            RETURN
          END IF
        END DO
      END IF
    END IF

    IF (dump_sum) THEN
      DO idir = 1, ndirs
        IF (PRESENT(dir_tags)) THEN
          CALL check_name_length('dir tag', &
              'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)))

          temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir))
          temp_name = &
              'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir))
        ELSE
          temp_block_id = TRIM(block_id)
          temp_name = 'Derived/' // TRIM(name)
        END IF

        IF (isubset /= 1) THEN
          CALL check_name_length('subset', TRIM(temp_name) &
              // '/Subset_' // TRIM(sub%name))

          temp_block_id = TRIM(temp_block_id) &
              // '/s_' // TRIM(sub%name)
          temp_name = TRIM(temp_name) &
              // '/Subset_' // TRIM(sub%name)
        END IF

        IF (PRESENT(fluxdir)) THEN
          CALL func(array, 0, fluxdir(idir))
        ELSE
          CALL func(array, 0)
        END IF

        IF (dump_skipped) THEN
          rnx = sub%n_local(1)
          rny = sub%n_local(2)

          ALLOCATE(reduced(rnx,rny))

          jj = sub%n_start(2) + 1
          DO j = 1, rny
            ii = sub%n_start(1) + 1
            DO i = 1, rnx
              reduced(i,j) = array(ii,jj)
              ii = ii + sub%skip_dir(1)
            END DO
            jj = jj + sub%skip_dir(2)
          END DO

          IF (convert) THEN
            rsubtype  = sub%subtype_r4
            rsubarray = sub%subarray_r4
          ELSE
            rsubtype  = sub%subtype
            rsubarray = sub%subarray
          END IF

          CALL check_name_length('subset', TRIM(temp_name) // '/Reduced')

          temp_grid_id = 'grid/r_' // TRIM(sub%name)
          temp_block_id = TRIM(temp_block_id) // '/r'
          temp_name = TRIM(temp_name) // '/Reduced'

          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
              TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

          sub%dump_field_grid = .TRUE.
        ELSE IF (dump_part) THEN
          temp_grid_id = 'grid/' // TRIM(sub%name)

          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
              temp_grid_id, array, rsubtype, rsubarray, convert)
          sub%dump_field_grid = .TRUE.
        ELSE
          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
              subtype, subarray, convert)
          dump_field_grid = .TRUE.
        END IF
      END DO
    END IF

    IF (dump_species .AND. dump_skipped) THEN
      rnx = sub%n_local(1)
      rny = sub%n_local(2)

      IF (.NOT.ALLOCATED(reduced)) ALLOCATE(reduced(rnx,rny))

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      END IF

      temp_grid_id = 'grid/r_' // TRIM(sub%name)

      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        DO idir = 1, ndirs
          IF (PRESENT(dir_tags)) THEN
            CALL check_name_length('dir tag', &
                'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name)
            temp_name = &
                'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name)
          ELSE
            CALL check_name_length('species', &
                'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) &
                // '/' // TRIM(io_list(ispecies)%name)
            temp_name = &
                'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name)
          END IF

          CALL check_name_length('subset', &
              TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

          temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
          temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

          IF (PRESENT(fluxdir)) THEN
            CALL func(array, ispecies, fluxdir(idir))
          ELSE
            CALL func(array, ispecies)
          END IF

          jj = sub%n_start(2) + 1
          DO j = 1, rny
            ii = sub%n_start(1) + 1
            DO i = 1, rnx
              reduced(i,j) = array(ii,jj)
              ii = ii + sub%skip_dir(1)
            END DO
            jj = jj + sub%skip_dir(2)
          END DO

          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
              TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

          sub%dump_field_grid = .TRUE.
        END DO
      END DO
    ELSE IF (dump_species) THEN
      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        DO idir = 1, ndirs
          IF (PRESENT(dir_tags)) THEN
            CALL check_name_length('dir tag', &
                'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name)
            temp_name = &
                'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                // '/' // TRIM(io_list(ispecies)%name)
          ELSE
            CALL check_name_length('species', &
                'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) &
                // '/' // TRIM(io_list(ispecies)%name)
            temp_name = &
                'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name)
          END IF

          IF (PRESENT(fluxdir)) THEN
            CALL func(array, ispecies, fluxdir(idir))
          ELSE
            CALL func(array, ispecies)
          END IF

          IF (dump_part) THEN
            ! First subset is main dump so there wont be any restrictions
            temp_grid_id = 'grid/' // TRIM(sub%name)

            CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
                temp_grid_id, array, rsubtype, rsubarray, convert)
            sub%dump_field_grid = .TRUE.
          ELSE
            CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
                subtype, subarray, convert)
            dump_field_grid = .TRUE.
          END IF
        END DO
      END DO
    END IF

    IF (ALLOCATED(reduced)) DEALLOCATE(reduced)

    IF (isubset /= 1) RETURN

    ! Write averaged data
    DO io = 1, n_io_blocks
      iob => io_block_list(io)
      IF (.NOT.iob%dump) CYCLE

      mask = iob%dumpmask(id)
      IF (IAND(mask, c_io_averaged) == 0) CYCLE

      avg => iob%averaged_data(id)
      IF (.NOT.avg%started) CYCLE

      IF (avg%dump_single) THEN
        avg%r4array = avg%r4array / REAL(avg%real_time, r4)

        IF (avg%species_sum > 0 .AND. IAND(mask, c_io_field) == 0) THEN
          DO idir = 1, avg%species_sum
            IF (PRESENT(dir_tags)) THEN
              CALL check_name_length('dir tag', &
                  'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged')

              temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged'
              temp_name = &
                  'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged'
            ELSE
              temp_block_id = TRIM(block_id) // '_averaged'
              temp_name = 'Derived/' // TRIM(name) // '_averaged'
            END IF

            CALL sdf_write_plain_variable(sdf_handle, &
                TRIM(temp_block_id), TRIM(temp_name), &
                TRIM(units), dims, stagger, 'grid', &
                avg%r4array(:,:,idir), subtype_field_r4, subarray_field_r4)
          END DO

          dump_field_grid = .TRUE.
        END IF

        IF (avg%n_species > 0) THEN
          iav = avg%species_sum
          DO ispecies = 1, avg%n_species / ndirs
            IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

            DO idir = 1, ndirs
              IF (PRESENT(dir_tags)) THEN
                CALL check_name_length('dir tag', &
                    'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name))

                temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
                temp_name = &
                    'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
              ELSE
                CALL check_name_length('species', &
                    'Derived/' // TRIM(name) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name))

                temp_block_id = TRIM(block_id) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
                temp_name = &
                    'Derived/' // TRIM(name) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
              END IF

              iav = iav + 1
              CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                  TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                  avg%r4array(:,:,iav), subtype_field_r4, subarray_field_r4)
            END DO

            dump_field_grid = .TRUE.
          END DO
        END IF

        avg%r4array = 0.0_num
      ELSE
        avg%array = avg%array / avg%real_time

        IF (avg%species_sum > 0 .AND. IAND(mask, c_io_field) == 0) THEN
          DO idir = 1, avg%species_sum
            IF (PRESENT(dir_tags)) THEN
              CALL check_name_length('dir tag', &
                  'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged')

              temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged'
              temp_name = &
                  'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                  // '_averaged'
            ELSE
              temp_block_id = TRIM(block_id) // '_averaged'
              temp_name = 'Derived/' // TRIM(name) // '_averaged'
            END IF

            CALL sdf_write_plain_variable(sdf_handle, &
                TRIM(temp_block_id), TRIM(temp_name), &
                TRIM(units), dims, stagger, 'grid', &
                avg%array(:,:,idir), subtype_field, subarray_field)
          END DO

          dump_field_grid = .TRUE.
        END IF

        IF (avg%n_species > 0) THEN
          iav = avg%species_sum
          DO ispecies = 1, avg%n_species / ndirs
            IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

            DO idir = 1, ndirs
              IF (PRESENT(dir_tags)) THEN
                CALL check_name_length('dir tag', &
                    'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name))

                temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
                temp_name = &
                    'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
              ELSE
                CALL check_name_length('species', &
                    'Derived/' // TRIM(name) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name))

                temp_block_id = TRIM(block_id) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
                temp_name = &
                    'Derived/' // TRIM(name) &
                    // '_averaged/' // TRIM(io_list(ispecies)%name)
              END IF

              iav = iav + 1
              CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                  TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                  avg%array(:,:,iav), subtype_field, subarray_field)
            END DO

            dump_field_grid = .TRUE.
          END DO
        END IF

        avg%array = 0.0_num
      END IF

      avg%real_time = 0.0_num
      avg%started = .FALSE.
    END DO

  END SUBROUTINE write_nspecies_field



  SUBROUTINE build_species_subset

    INTEGER :: i, l
    TYPE(particle), POINTER :: current, next
    LOGICAL :: use_particle
    REAL(num) :: part_mc
    TYPE(subset), POINTER :: sub
    CLASS(particle_id_hash), POINTER :: current_hash

    IF (done_subset_init) RETURN
    done_subset_init = .TRUE.

    ! Update particle ID values if using
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    DO i = 1, n_species 
      CALL generate_particle_ids(species_list(i)%attached_list)
    END DO
#endif

    IF (isubset == 1) THEN
      io_list => species_list
      RETURN
    END IF

    io_list => io_list_data

    l = isubset - 1
    sub => subset_list(l)
    current_hash => id_registry%get_existing_hash(sub%name)
    DO i = 1, n_species
      io_list(i) = species_list(i)
      io_list(i)%count = 0
      io_list(i)%name = 'subset_' // TRIM(sub%name) // '/' &
          // TRIM(species_list(i)%name)
      CALL create_empty_partlist(io_list(i)%attached_list)

      IF (.NOT. sub%use_species(i)) THEN
        io_list(i)%dumpmask = c_io_never
        CYCLE
      END IF

      IF (sub%persistent) THEN
        IF (time < sub%persist_start_time &
            .AND. step < sub%persist_start_step) THEN
          io_list(i)%dumpmask = c_io_never
          CYCLE
        END IF
      END IF

      IF (sub%persistent) any_persistent_subset = .TRUE.
      io_list(i)%dumpmask = sub%mask

      part_mc = c * species_list(i)%mass

      current => species_list(i)%attached_list%head
      DO WHILE (ASSOCIATED(current))
        next => current%next
        use_particle = .TRUE.

        IF (sub%persistent .AND. sub%locked) THEN
          use_particle = current_hash%holds(current)
        ELSE
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc = c * current%mass
#endif
          use_particle = test_particle(sub, current, part_mc)
        END IF

        IF (use_particle) THEN
          ! Move particle to io_list
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(io_list(i)%attached_list, current)
        END IF
        current => next
      END DO
    END DO

  END SUBROUTINE build_species_subset



  SUBROUTINE build_persistent_subsets

    INTEGER :: isub, ispec
    TYPE(particle), POINTER :: current, next
    LOGICAL :: use_particle
    REAL(num) :: part_mc
    TYPE(subset), POINTER :: sub
    CLASS(particle_id_hash), POINTER :: current_hash

    IF (.NOT. any_persistent_subset) RETURN

    DO isub = 1, SIZE(subset_list)
      sub => subset_list(isub)

      ! Not a persistent subset
      IF (.NOT. sub%persistent) CYCLE
      ! Already locked in
      IF (sub%locked) CYCLE
      ! Not yet time to lock
      IF (time < sub%persist_start_time &
          .AND. step < sub%persist_start_step) CYCLE

      current_hash => id_registry%get_existing_hash(sub%name)
      DO ispec = 1, n_species
        IF (.NOT. sub%use_species(ispec)) THEN
          CYCLE
        END IF

        part_mc = c * species_list(ispec)%mass

        current => species_list(ispec)%attached_list%head
        DO WHILE (ASSOCIATED(current))
          next => current%next
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc = c * current%mass
#endif
          use_particle = test_particle(sub, current, part_mc)

          ! Add particle ID to persistence list
          IF (use_particle) CALL current_hash%add(current)

          current => next
        END DO
      END DO

      sub%locked = .TRUE.
      CALL current_hash%optimise()
    END DO

  END SUBROUTINE build_persistent_subsets



  FUNCTION test_particle(sub, current, part_mc) RESULT(use_particle)

    TYPE(subset), INTENT(IN) :: sub
    TYPE(particle), INTENT(IN) :: current
    REAL(num), INTENT(INOUT) :: part_mc
    LOGICAL :: use_particle
    REAL(num) :: gamma_rel, random_num
    INTEGER :: n

    use_particle = .TRUE.

    IF (sub%use_gamma) THEN
#ifdef PER_PARTICLE_CHARGE_MASS
      part_mc = c * current%mass
#endif
      gamma_rel = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num)

      n = c_subset_gamma_min
      IF (sub%use_restriction(n)) THEN
        IF (gamma_rel < sub%restriction(n)) &
            use_particle = .FALSE.
      END IF

      n = c_subset_gamma_max
      IF (sub%use_restriction(n)) THEN
        IF (gamma_rel > sub%restriction(n)) &
            use_particle = .FALSE.
      END IF
    END IF

    n = c_subset_x_min
    IF (sub%use_restriction(n)) THEN
      IF (current%part_pos(1) < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_x_max
    IF (sub%use_restriction(n)) THEN
      IF (current%part_pos(1) > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_y_min
    IF (sub%use_restriction(n)) THEN
      IF (current%part_pos(2) < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_y_max
    IF (sub%use_restriction(n)) THEN
      IF (current%part_pos(2) > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_px_min
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(1) < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_px_max
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(1) > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_py_min
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(2) < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_py_max
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(2) > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_pz_min
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(3) < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_pz_max
    IF (sub%use_restriction(n)) THEN
      IF (current%part_p(3) > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

#ifndef PER_SPECIES_WEIGHT
    n = c_subset_weight_min
    IF (sub%use_restriction(n)) THEN
      IF (current%weight < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_weight_max
    IF (sub%use_restriction(n)) THEN
      IF (current%weight > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    n = c_subset_charge_min
    IF (sub%use_restriction(n)) THEN
      IF (current%charge < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_charge_max
    IF (sub%use_restriction(n)) THEN
      IF (current%charge > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_mass_min
    IF (sub%use_restriction(n)) THEN
      IF (current%mass < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_mass_max
    IF (sub%use_restriction(n)) THEN
      IF (current%mass > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    n = c_subset_id_min
    IF (sub%use_restriction(n)) THEN
      IF (current%id < sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

    n = c_subset_id_max
    IF (sub%use_restriction(n)) THEN
      IF (current%id > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF
#endif
    n = c_subset_random
    IF (sub%use_restriction(n)) THEN
      random_num = random()
      IF (random_num > sub%restriction(n)) &
          use_particle = .FALSE.
    END IF

  END FUNCTION test_particle



  SUBROUTINE species_offset_init

    INTEGER(i8) :: species_count
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_species_per_proc
    INTEGER :: i, ispecies
    TYPE(particle_species), POINTER :: spec

    IF (done_species_offset_init) RETURN
    done_species_offset_init = .TRUE.

    IF (.NOT.ALLOCATED(species_offset))  ALLOCATE(species_offset(n_species))

    CALL build_species_subset

    ALLOCATE(npart_species_per_proc(nproc))
    species_offset = 0

    npart_global = 0
    DO ispecies = 1, n_species
      spec => io_list(ispecies)

      CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
          npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
      species_count = 0
      DO i = 1, nproc
        IF (rank == i-1) species_offset(ispecies) = species_count
        species_count = species_count + npart_species_per_proc(i)
      END DO
      spec%count = species_count
      spec%count_update_step = step
      npart_global = npart_global + species_count

      CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
          'CPU split/' // TRIM(spec%name), npart_species_per_proc)
    END DO

    IF (track_ejected_particles &
        .AND. .NOT.ALLOCATED(ejected_offset)) THEN
      ALLOCATE(ejected_offset(n_species))
      ejected_offset = 0

      DO ispecies = 1, n_species
        spec => ejected_list(ispecies)

        CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
            npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
        species_count = 0
        DO i = 1, nproc
          IF (rank == i-1) ejected_offset(ispecies) = species_count
          species_count = species_count + npart_species_per_proc(i)
        END DO
        spec%count = species_count
        spec%count_update_step = step

        CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
            'CPU split/' // TRIM(spec%name), npart_species_per_proc)
      END DO
    END IF

    DEALLOCATE(npart_species_per_proc)

  END SUBROUTINE species_offset_init



  SUBROUTINE write_particle_grid(code)

    INTEGER, INTENT(IN) :: code
    INTEGER :: ispecies, id, mask, io
    LOGICAL :: convert, dump_grid, restart_id

    id = c_dump_part_grid
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    use_offset_grid = .FALSE.
    DO io = 1, n_io_blocks
      use_offset_grid = use_offset_grid &
         .OR. (io_block_list(io)%dump .AND. io_block_list(io)%use_offset_grid)
    END DO

    IF (restart_id .OR. (IAND(mask, c_io_never) == 0 &
        .AND. (IAND(mask, code) /= 0 .OR. ANY(dump_point_grid)))) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        mask = current_species%dumpmask
        dump_grid = dump_point_grid(ispecies)

        IF (IAND(mask, code) /= 0) dump_grid = .TRUE.
        IF (IAND(mask, c_io_never) /= 0) dump_grid = .FALSE.
        IF (IAND(code, c_io_restartable) /= 0) dump_grid = .TRUE.

        IF (dump_grid) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          CALL sdf_write_point_mesh(sdf_handle, &
              'grid/' // TRIM(current_species%name), &
              'Grid/Particles/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              io_list(ispecies)%count, c_dimension_2d, &
              it_output_position, species_offset(ispecies), convert)
        END IF
      END DO
    END IF

    IF (isubset /= 1) RETURN

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        CALL sdf_write_point_mesh(sdf_handle, &
            'grid/' // TRIM(current_species%name), &
            'Grid/Particles/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            ejected_list(ispecies)%count, c_dimension_2d, &
            it_output_position, ejected_offset(ispecies), convert)
      END DO
    END IF

  END SUBROUTINE write_particle_grid



  SUBROUTINE write_particle_variable_num(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE constants
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (restart_id &
        .OR. (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0)) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        END IF
      END DO
    END IF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      END DO
    END IF

  END SUBROUTINE write_particle_variable_num



#if defined(PARTICLE_ID4) || defined(PARTICLE_DEBUG)
  SUBROUTINE write_particle_variable_i4(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npart_it, start, param)
        USE constants
        INTEGER(i4) :: iterator
        INTEGER(i4), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (restart_id &
        .OR. (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0)) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        END IF
      END DO
    END IF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      END DO
    END IF

  END SUBROUTINE write_particle_variable_i4
#endif



  SUBROUTINE write_particle_variable_i8(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npart_it, start, param)
        USE constants
        INTEGER(i8) :: iterator
        INTEGER(i8), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (restart_id &
        .OR. (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0)) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        END IF
      END DO
    END IF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      END DO
    END IF

  END SUBROUTINE write_particle_variable_i8



  SUBROUTINE create_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 15
    INTEGER :: days, hours, minutes, seconds, frac_seconds

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    WRITE(timestring, '(i3,'':'',i2.2,'':'',i2.2,'':'',i2.2,''.'',i2.2)') &
        days, hours, minutes, seconds, frac_seconds

  END SUBROUTINE create_timestring



  SUBROUTINE create_full_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 48
    INTEGER :: days, hours, minutes, seconds, frac_seconds, var
    CHARACTER(LEN=8) :: varstring
    CHARACTER(LEN=4) :: intstring, fracstring
    LOGICAL :: string_started

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    timestring = ''
    string_started = .FALSE.

    var = days
    varstring = ' day'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      END IF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    END IF

    var = hours
    varstring = ' hour'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      END IF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    END IF

    var = minutes
    varstring = ' minute'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      END IF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    END IF

    var = seconds
    varstring = ' seconds'
    IF (var > 0 .OR. frac_seconds > 0 .OR. .NOT.string_started) THEN
      CALL integer_as_string(var, intstring)
      WRITE(fracstring, '(i2.2)') frac_seconds
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      END IF
    END IF

  END SUBROUTINE create_full_timestring



  SUBROUTINE cleanup_stop_files()

    INTEGER :: ierr

    IF (rank /= 0) RETURN

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file), iostat=ierr)
    IF (ierr == 0) CLOSE(lu, status='DELETE')

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump), iostat=ierr)
    IF (ierr == 0) CLOSE(lu, status='DELETE')

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(request_dump_file), iostat=ierr)
    IF (ierr == 0) CLOSE(lu, status='DELETE')

  END SUBROUTINE cleanup_stop_files



  SUBROUTINE check_for_stop_condition(halt, force_dump)

    LOGICAL, INTENT(OUT) :: halt, force_dump
    INTEGER :: ierr
    INTEGER, SAVE :: check_counter = 0
    LOGICAL :: buffer(4), got_stop_condition, got_stop_file
    REAL(num) :: walltime

    IF (check_stop_frequency <= 0 .AND. .NOT.check_walltime) RETURN

    walltime = -1.0_num
    IF (check_walltime) &
        CALL check_walltime_auto(walltime, halt)

    IF (halt) THEN
      force_dump = .TRUE.
      RETURN
    END IF

    IF (check_stop_frequency < 0) RETURN

    got_stop_condition = .FALSE.
    force_dump = .FALSE.
    check_counter = check_counter + 1

    IF (check_counter < check_stop_frequency) RETURN
    check_counter = 0

    IF (rank == 0) THEN
      ! Since we're checking for a STOP file, we might as well check the
      ! walltime as well
      IF (check_walltime) THEN
        IF (walltime < 0.0_num) walltime = MPI_WTIME()
        IF (walltime - real_walltime_start >= stop_at_walltime) THEN
          got_stop_condition = .TRUE.
          force_dump = .TRUE.
          PRINT*,'Stopping because "stop_at_walltime" has been exceeded.'
        END IF
      END IF

      ! Next check if stop file exists
      OPEN(unit=lu, status='OLD', iostat=ierr, &
          file=TRIM(data_dir) // '/' // TRIM(stop_file))
      IF (ierr == 0) THEN
        got_stop_file = .TRUE.
        got_stop_condition = .TRUE.
        force_dump = .TRUE.
        CLOSE(lu, status='DELETE')
      ELSE
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump))
        IF (ierr == 0) THEN
          got_stop_file = .TRUE.
          got_stop_condition = .TRUE.
          force_dump = .FALSE.
          CLOSE(lu, status='DELETE')
        ELSE
          got_stop_file = .FALSE.
          ! If no stop files are found, check if a dump file was requested
          OPEN(unit=lu, status='OLD', iostat=ierr, &
              file=TRIM(data_dir) // '/' // TRIM(request_dump_file))
          IF (ierr == 0) THEN
            READ(lu,'(A)',iostat=ierr) request_dump_name
            CLOSE(lu, status='DELETE')
            IF (ierr == 0) THEN
              got_request_dump_name = .TRUE.
              WRITE(stat_unit,'(a)') 'Found DUMP file output request with ' &
                  // 'contents: ' // TRIM(request_dump_name)
            ELSE
              got_request_dump_restart = .TRUE.
              WRITE(stat_unit,'(a)') 'Found DUMP file output request'
            END IF
            CALL flush_stat_file()
          ELSE
            got_request_dump_name = .FALSE.
            got_request_dump_restart = .FALSE.
          END IF
        END IF
      END IF

      IF (got_stop_file) PRINT*,'Stopping because "STOP" file has been found.'

      buffer(1) = got_stop_condition
      buffer(2) = force_dump
      buffer(3) = got_request_dump_name
      buffer(4) = got_request_dump_restart
    END IF

    CALL MPI_BCAST(buffer, 4, MPI_LOGICAL, 0, comm, errcode)
    got_stop_condition = buffer(1)
    force_dump = buffer(2)
    got_request_dump_name = buffer(3)
    got_request_dump_restart = buffer(4)
    got_stop_restart = got_stop_condition

    IF (got_request_dump_name) THEN
      CALL MPI_BCAST(request_dump_name, string_length, MPI_CHARACTER, 0, &
                     comm, errcode)
    END IF

    IF (got_stop_condition) halt = .TRUE.

  END SUBROUTINE check_for_stop_condition



  SUBROUTINE check_walltime_auto(walltime, halt)

    REAL(num), INTENT(INOUT) :: walltime
    LOGICAL, INTENT(OUT) :: halt
    INTEGER, PARAMETER :: tag = 2001
    INTEGER :: msg, request, i
    INTEGER :: status_ignore(MPI_STATUS_SIZE)
    LOGICAL :: flag
    LOGICAL, ALLOCATABLE, SAVE :: completed(:)
    LOGICAL, SAVE :: yet_to_sync = .TRUE.
    LOGICAL, SAVE :: first = .TRUE.
    LOGICAL, SAVE :: all_completed = .FALSE.
    REAL(num), SAVE :: wall0
    REAL(num), PARAMETER :: frac = 1.0_num
    REAL(num) :: timeout

    halt = all_completed
    IF (all_completed) RETURN

    IF (walltime < 0) walltime = MPI_WTIME()
    IF ((walltime + timer_average(c_timer_step) + timer_average(c_timer_io) &
        + timer_average(c_timer_balance) - real_walltime_start) &
        < frac * stop_at_walltime) RETURN

    IF (rank == 0) THEN
      IF (first) THEN
        ALLOCATE(completed(nproc-1))
        completed = .FALSE.
        first = .FALSE.
      END IF
      wall0 = walltime
      timeout = 2.0_num * timer_average(c_timer_step)
      DO
        all_completed = .TRUE.
        DO i = 1,nproc-1
          IF (.NOT.completed(i)) THEN
            ! The Platform-MPI interface for MPI_IPROBE is broken and will not
            ! allow us to use MPI_STATUS_IGNORE
            CALL MPI_IPROBE(i, tag, comm, flag, status_ignore, errcode)
            completed(i) = flag
            IF (flag) THEN
              CALL MPI_RECV(msg, 0, MPI_INTEGER, i, tag, comm, &
                  MPI_STATUS_IGNORE, errcode)
            ELSE
              all_completed = .FALSE.
            END IF
          END IF
        END DO
        IF (all_completed) THEN
          DEALLOCATE(completed)
          msg = -1
          DO i = 1,nproc-1
            CALL MPI_ISEND(msg, 1, MPI_INTEGER, i, tag, comm, request, errcode)
            CALL MPI_REQUEST_FREE(request, errcode)
          END DO
          halt = all_completed
          PRINT*,'Stopping because "stop_at_walltime" has been exceeded.'
          RETURN
        END IF
        walltime = MPI_WTIME()
        IF (walltime - wall0 > timeout) THEN
          msg = 1
          DO i = 1,nproc-1
            IF (completed(i)) THEN
              CALL MPI_ISEND(msg, 1, MPI_INTEGER, i, tag, comm, request, &
                  errcode)
              CALL MPI_REQUEST_FREE(request, errcode)
            END IF
          END DO
          RETURN
        END IF
      END DO
      RETURN
    END IF

    IF (yet_to_sync) THEN
      CALL MPI_ISEND(0, 0, MPI_INTEGER, 0, tag, comm, request, errcode)
      CALL MPI_REQUEST_FREE(request, errcode)
      yet_to_sync = .FALSE.
    END IF
    CALL MPI_RECV(msg, 1, MPI_INTEGER, 0, tag, comm, &
        MPI_STATUS_IGNORE, errcode)
    IF (msg < 0) all_completed = .TRUE.
    halt = all_completed

  END SUBROUTINE check_walltime_auto



  SUBROUTINE epoch_write_source_info(h)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=c_id_length) :: stitched_ids(3)
    CHARACTER(LEN=c_id_length) :: time_string
    CHARACTER(LEN=512) :: string_array(6)
    INTEGER :: n, i

    n = 0

    IF (SIZE(epoch_bytes) > 1 .OR. (TRIM(epoch_bytes_checksum_type) /= '' &
        .AND. ICHAR(epoch_bytes_checksum_type(1:1)) /= 0)) THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'epoch_source/source', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'EPOCH source code', epoch_bytes, &
          epoch_bytes_padding, epoch_bytes_mimetype, &
          epoch_bytes_checksum_type, epoch_bytes_checksum)
    END IF

    IF (SIZE(epoch_bytes) == 1 .AND. SIZE(epoch_diff_bytes) > 1) THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'epoch_source/diff', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'EPOCH repository differences', epoch_diff_bytes, &
          epoch_diff_bytes_padding, epoch_diff_bytes_mimetype, &
          epoch_diff_bytes_checksum_type, epoch_diff_bytes_checksum)
    END IF

    n = n + 1
    CALL sdf_safe_copy_id(h, 'epoch_source/info', stitched_ids(n))
    WRITE(time_string, '(I20)') epoch_bytes_compile_date

    string_array(1) = trim_string(epoch_bytes_git_version)
    string_array(2) = trim_string(epoch_bytes_compile_date_string)
    string_array(3) = trim_string(time_string)
    string_array(4) = trim_string(epoch_bytes_compile_machine_info)
    string_array(5) = trim_string(epoch_bytes_compiler_info)
    string_array(6) = trim_string(epoch_bytes_compiler_flags)

    ! Prevent truncation warning
    DO i = 1, 6
      string_array(i)(h%string_length:512) = ACHAR(0)
    END DO

    CALL sdf_write_namevalue(h, stitched_ids(n), &
        'EPOCH repository information', &
        (/'git_version         ', &
          'compile_date_string ', &
          'compile_date_seconds', &
          'compile_machine_info', &
          'compiler_info       ', &
          'compiler_flags      '/), string_array)

    CALL sdf_write_stitched(h, 'epoch_source', 'EPOCH source', &
        stitched_ids(1), c_stagger_cell_centre, stitched_ids, n)

  END SUBROUTINE epoch_write_source_info



  SUBROUTINE write_source_info(h)

    TYPE(sdf_file_handle) :: h

    CALL sdf_write_source_info(h)
    CALL epoch_write_source_info(h)
    !CALL write_input_decks(h)

  END SUBROUTINE write_source_info

END MODULE diagnostics
