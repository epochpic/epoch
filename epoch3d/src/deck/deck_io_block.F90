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

MODULE deck_io_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: io_deck_initialise, io_deck_finalise
  PUBLIC :: io_block_start, io_block_end
  PUBLIC :: io_block_handle_element, io_block_check, copy_io_block

  INTEGER :: block_number, nfile_prefixes
  INTEGER :: rolling_restart_io_block
  LOGICAL, DIMENSION(num_vars_to_dump) :: io_block_done
  LOGICAL, PRIVATE :: got_name, got_dump_source_code, got_dump_input_decks
  LOGICAL, PRIVATE :: warning_printed, got_dt_average
  LOGICAL, PRIVATE :: dumpmask_warned
  CHARACTER(LEN=c_id_length), ALLOCATABLE :: io_prefixes(:)
  TYPE(io_block_type), POINTER :: io_block

CONTAINS

  SUBROUTINE io_deck_initialise

    block_number = 0
    IF (deck_state /= c_ds_first) RETURN

    any_average = .FALSE.
    warning_printed = .FALSE.
    dumpmask_warned = .FALSE.

    track_ejected_particles = .FALSE.
    dump_absorption = .FALSE.
    averaged_var_block = 0
    new_style_io_block = .FALSE.
    n_io_blocks = 0
    rolling_restart_io_block = 0
    got_dt_average = .FALSE.

  END SUBROUTINE io_deck_initialise



  SUBROUTINE io_deck_finalise

    INTEGER :: i, io, iu, n_zeros_estimate, n_dumps
    REAL(num) :: dumps, dx
#ifndef NO_IO
    CHARACTER(LEN=c_max_path_length) :: list_filename
#endif

    n_io_blocks = block_number
    block_number = 0

    IF (n_io_blocks > 0) THEN
      ! First pass
      IF (deck_state == c_ds_first) THEN
        IF (.NOT.new_style_io_block .AND. n_io_blocks /= 1) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Cannot have multiple unnamed "output" blocks.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF

        ALLOCATE(io_prefixes(n_io_blocks+1))
        nfile_prefixes = 1
        io_prefixes(1) = ''

        ALLOCATE(io_block_list(n_io_blocks))
        DO i = 1, n_io_blocks
          CALL init_io_block(io_block_list(i))
          IF (i == rolling_restart_io_block) THEN
            io_block_list(i)%rolling_restart = .TRUE.
            io_block_list(i)%dump_cycle = 1
            io_block_list(i)%restart = .TRUE.
            io_block_list(i)%prefix_index = 2
            nfile_prefixes = 2
            io_prefixes(2) = 'roll'
          END IF
        END DO
      ! Second pass
      ELSE
        dumps = 1.0_num

        DO i = 1, n_io_blocks
          IF (io_block_list(i)%disabled) CYCLE

          IF (use_offset_grid) io_block%use_offset_grid = use_offset_grid

          IF (io_block_list(i)%dt_average > t_end) THEN
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Averaging time is longer than t_end, will set', &
                    ' averaging time equal'
                WRITE(io,*) 'to t_end.'
              END DO
            END IF
            io_block_list(i)%dt_average = t_end
          END IF

          IF (io_block_list(i)%dump_cycle_first_index &
                > io_block_list(i)%dump_cycle) THEN
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) '"dump_cycle_first_index" cannot be greater ', &
                    'than "dump_cycle"'
                WRITE(io,*) 'Resetting to zero.'
              END DO
            END IF
            io_block_list(i)%dump_cycle_first_index = 0
          END IF

          IF (io_block_list(i)%dt_snapshot > 0.0_num) THEN
            dumps = MAX(dumps, t_end / io_block_list(i)%dt_snapshot)
          END IF

          ! This might fail if Debye length not resolved
          IF (io_block_list(i)%nstep_snapshot > 0) THEN
            dx = (x_max - x_min) / nx_global
            dx = MIN(dx, (y_max - y_min) / ny_global)
            dx = MIN(dx, (z_max - z_min) / nz_global)
            dumps = MAX(dumps, t_end * c / dx / io_block_list(i)%nstep_snapshot)
          END IF
        END DO

        n_dumps = FLOOR(dumps)

        IF (io_block%dump_last .AND. dumps - n_dumps > c_tiny) &
            n_dumps = n_dumps + 1

        IF (.NOT. io_block%dump_first) &
            n_dumps = n_dumps - 1

        n_zeros_estimate = MAX(n_zeros, FLOOR(LOG10(REAL(MAX(n_dumps, 1)))) + 1)

        IF (n_zeros_control > 0) use_accurate_n_zeros = .FALSE.

        IF (.NOT.use_accurate_n_zeros) THEN
          IF (n_zeros_control > 0 &
              .AND. n_zeros_estimate /= n_zeros_control) THEN
            n_zeros = n_zeros_estimate
            IF (n_zeros > n_zeros_control .AND. rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,'(A,I1,A)') ' Estimated value of n_zeros (', n_zeros, &
                    ') has been overidden by input deck'
              END DO
            END IF
            n_zeros = n_zeros_control
          ELSE IF (n_zeros_estimate > n_zeros) THEN
            n_zeros = n_zeros_estimate
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,'(A,I1,A)') ' n_zeros changed to ', n_zeros, &
                    ' to accomodate requested number of snapshots'
              END DO
            END IF
          END IF
        END IF

        ALLOCATE(file_prefixes(nfile_prefixes))
        ALLOCATE(file_numbers(nfile_prefixes))
        DO i = 1,nfile_prefixes
          file_prefixes(i) = TRIM(io_prefixes(i))
          file_numbers(i) = 0
        END DO
        DEALLOCATE(io_prefixes)

#ifndef NO_IO
        ! Remove any left-over VisIt file lists
        IF (.NOT.ic_from_restart .AND. rank == 0) THEN
          list_filename = TRIM(ADJUSTL(data_dir)) // '/full.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
          list_filename = TRIM(ADJUSTL(data_dir)) // '/normal.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
          list_filename = TRIM(ADJUSTL(data_dir)) // '/restart.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
        END IF
#endif
      END IF
    END IF

  END SUBROUTINE io_deck_finalise



  SUBROUTINE io_block_start

    io_block_done = .FALSE.
    got_name = .FALSE.
    got_dump_source_code = .FALSE.
    got_dump_input_decks = .FALSE.
    block_number = block_number + 1
    IF (deck_state /= c_ds_first .AND. block_number > 0) THEN
      io_block => io_block_list(block_number)
    END IF

  END SUBROUTINE io_block_start



  SUBROUTINE io_block_end

    INTEGER :: io, iu, mask
#ifndef NO_IO
    CHARACTER(LEN=c_max_path_length) :: list_filename
#endif

    IF (deck_state == c_ds_first) RETURN

    IF (io_block%disabled) RETURN

    mask = io_block%dumpmask(c_dump_ejected_particles)
    IF (mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) &
        track_ejected_particles = .TRUE.

    IF (io_block%dumpmask(c_dump_absorption) /= c_io_never) THEN
      dump_absorption = .TRUE.
    END IF

    IF (.NOT. got_name) THEN
      IF (new_style_io_block) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot mix old and new style output blocks.'
            WRITE(io,*) 'You can either have multiple, named output blocks ', &
                'or a single unnamed one.'
          END DO
        END IF
        CALL abort_code(c_err_bad_value)
      END IF
      io_block%name = 'normal'
    END IF

#ifndef NO_IO
    ! Delete any existing visit file lists
    IF (.NOT.ic_from_restart .AND. rank == 0) THEN
      list_filename = TRIM(ADJUSTL(data_dir)) // '/' &
          // TRIM(io_block%name) // '.visit'
      OPEN(unit=lu, status='UNKNOWN', file=list_filename)
      CLOSE(unit=lu, status='DELETE')
    END IF
#endif
    io_block%dumpmask(c_dump_jx) = IOR(io_block%dumpmask(c_dump_jx), c_io_field)
    io_block%dumpmask(c_dump_jy) = IOR(io_block%dumpmask(c_dump_jy), c_io_field)
    io_block%dumpmask(c_dump_jz) = IOR(io_block%dumpmask(c_dump_jz), c_io_field)

    IF (.NOT.got_dump_source_code) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_source_code = .TRUE.
    END IF

    IF (.NOT.got_dump_input_decks) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_input_decks = .TRUE.
    END IF

    IF (rank == 0 .AND. new_style_io_block .AND. .NOT.dumpmask_warned) THEN
      IF (ANY(IAND(io_block%dumpmask,c_io_restartable) == c_io_restartable) &
          .OR. ANY(IAND(io_block%dumpmask,c_io_full) == c_io_full)) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'The use of "full" and "restart" as dumpmasks in new-style ', &
                'output blocks is '
        PRINT*, 'deprecated and will be removed in a future version.'
        PRINT*
        dumpmask_warned = .TRUE.
      END IF
    END IF

    CALL set_restart_dumpmasks

  END SUBROUTINE io_block_end



  FUNCTION io_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, style_error
    INTEGER :: elementselected, mask, fullmask = 0, mask_element
    INTEGER :: i, is, subset, n_list, io, iu
    INTEGER, ALLOCATABLE :: subsets(:)
    LOGICAL :: bad, found, got_element
    INTEGER, PARAMETER :: c_err_new_style_ignore = 1
    INTEGER, PARAMETER :: c_err_new_style_global = 2
    INTEGER, PARAMETER :: c_err_old_style_ignore = 3

    errcode = c_err_none
    IF (value == blank) RETURN

    IF (deck_state == c_ds_first) THEN
      IF (str_cmp(element, 'name')) THEN
        new_style_io_block = .TRUE.

      ELSE IF (str_cmp(element, 'rolling_restart')) THEN
        IF (rolling_restart_io_block > 0) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Cannot have multiple "rolling_restart" blocks.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF
        rolling_restart_io_block = block_number
      END IF

      RETURN
    END IF

    IF (io_block%disabled) THEN
      errcode = c_err_none
      RETURN
    END IF

    errcode = c_err_none
    style_error = c_err_none
    got_element = .TRUE.
    elementselected = 0

    IF (str_cmp(element, 'dt_snapshot')) THEN
      io_block%dt_snapshot = as_real_print(value, element, errcode)
      IF (io_block%dt_snapshot < 0.0_num) io_block%dt_snapshot = 0.0_num

    ELSE IF (str_cmp(element, 'full_dump_every')) THEN
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        full_dump_every = as_integer_print(value, element, errcode)
        IF (full_dump_every == 0) full_dump_every = 1
      END IF

    ELSE IF (str_cmp(element, 'restart_dump_every')) THEN
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        restart_dump_every = as_integer_print(value, element, errcode)
        IF (restart_dump_every == 0) restart_dump_every = 1
      END IF

    ELSE IF (str_cmp(element, 'force_first_to_be_restartable')) THEN
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_first_to_be_restartable = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'force_final_to_be_restartable') &
        .OR. str_cmp(element, 'force_last_to_be_restartable')) THEN
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_final_to_be_restartable = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_offset_grid')) THEN
      io_block%use_offset_grid = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'extended_io_file')) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'The "extended_io_file" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        END DO
      END IF
      CALL abort_code(c_err_unknown_element)

    ELSE IF (str_cmp(element, 'dt_average') &
        .OR. str_cmp(element, 'averaging_period')) THEN
      io_block%dt_average = as_real_print(value, element, errcode)
      got_dt_average = .TRUE.

    ELSE IF (str_cmp(element, 'nstep_average') &
        .OR. str_cmp(element, 'min_cycles_per_average')) THEN
      io_block%nstep_average = as_integer_print(value, element, errcode)
      got_dt_average = .TRUE.

    ELSE IF (str_cmp(element, 'nstep_snapshot')) THEN
      io_block%nstep_snapshot = as_integer_print(value, element, errcode)
      IF (io_block%nstep_snapshot < 0) io_block%nstep_snapshot = 0

    ELSE IF (str_cmp(element, 'dump_source_code')) THEN
      io_block%dump_source_code = as_logical_print(value, element, errcode)
      got_dump_source_code = .TRUE.

    ELSE IF (str_cmp(element, 'dump_input_decks')) THEN
      io_block%dump_input_decks = as_logical_print(value, element, errcode)
      got_dump_input_decks = .TRUE.

    ELSE IF (str_cmp(element, 'dump_first')) THEN
      io_block%dump_first = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dump_last') &
        .OR. str_cmp(element, 'dump_final')) THEN
      io_block%dump_last = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'restartable')) THEN
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      io_block%restart = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'name')) THEN
      DO i = 1,block_number
        IF (TRIM(io_block_list(i)%name) == TRIM(value)) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Output block "' // TRIM(value) &
                  // '" already defined.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF
      END DO
      io_block%name = value
      got_name = .TRUE.

    ELSE IF (str_cmp(element, 'time_start')) THEN
      io_block%time_start = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'time_stop')) THEN
      io_block%time_stop = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'nstep_start')) THEN
      io_block%nstep_start = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'nstep_stop')) THEN
      io_block%nstep_stop = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dump_at_nsteps') &
        .OR. str_cmp(element, 'nsteps_dump')) THEN
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_nsteps, errcode)

    ELSE IF (str_cmp(element, 'dump_at_times') &
        .OR. str_cmp(element, 'times_dump')) THEN
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_times, errcode)

    ELSE IF (str_cmp(element, 'dump_cycle')) THEN
      io_block%dump_cycle = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'file_prefix')) THEN
      found = .FALSE.
      DO i = 1,nfile_prefixes
        IF (TRIM(io_prefixes(i)) == TRIM(value)) THEN
          found = .TRUE.
          io_block%prefix_index = i
          EXIT
        END IF
      END DO
      IF (.NOT.found) THEN
        nfile_prefixes = nfile_prefixes + 1
        io_prefixes(nfile_prefixes) = TRIM(value)
        io_block%prefix_index = nfile_prefixes
      END IF

    ELSE IF (str_cmp(element, 'dump_cycle_first_index')) THEN
      io_block%dump_cycle_first_index = &
          as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'filesystem')) THEN
      filesystem = TRIM(value) // ':'

    ELSE IF (str_cmp(element, 'dump_first_after_restart')) THEN
      io_block%dump_first_after_restart = &
          as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dump_at_walltimes') &
        .OR. str_cmp(element, 'walltimes_dump')) THEN
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_walltimes, errcode)

    ELSE IF (str_cmp(element, 'walltime_interval') &
        .OR. str_cmp(element, 'walltime_snapshot')) THEN
      io_block%walltime_interval = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'walltime_start')) THEN
      io_block%walltime_start = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'walltime_stop')) THEN
      io_block%walltime_stop = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'disabled')) THEN
      io_block%disabled = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'rolling_restart')) THEN
      ! Only handled on first parse

    ELSE IF (str_cmp(element, 'particles') &
        .OR. str_cmp(element, 'particle_grid')) THEN
      elementselected = c_dump_part_grid

    ELSE IF (str_cmp(element, 'grid') &
        .OR. str_cmp(element, 'field_grid')) THEN
      elementselected = c_dump_grid

    ELSE IF (str_cmp(element, 'species_id')) THEN
      elementselected = c_dump_part_species

    ELSE IF (str_cmp(element, 'particle_weight') &
        .OR. str_cmp(element, 'weight')) THEN
      elementselected = c_dump_part_weight

    ELSE IF (str_cmp(element, 'px')) THEN
      elementselected = c_dump_part_px

    ELSE IF (str_cmp(element, 'py')) THEN
      elementselected = c_dump_part_py

    ELSE IF (str_cmp(element, 'pz')) THEN
      elementselected = c_dump_part_pz

    ELSE IF (str_cmp(element, 'vx')) THEN
      elementselected = c_dump_part_vx

    ELSE IF (str_cmp(element, 'vy')) THEN
      elementselected = c_dump_part_vy

    ELSE IF (str_cmp(element, 'vz')) THEN
      elementselected = c_dump_part_vz

    ELSE IF (str_cmp(element, 'charge')) THEN
      elementselected = c_dump_part_charge

    ELSE IF (str_cmp(element, 'mass') &
        .OR. str_cmp(element, 'rest_mass')) THEN
      elementselected = c_dump_part_mass

    ELSE IF (str_cmp(element, 'id')) THEN
      elementselected = c_dump_part_id

    ELSE IF (str_cmp(element, 'ek') &
        .OR. str_cmp(element, 'particle_energy')) THEN
      elementselected = c_dump_part_ek

    ELSE IF (str_cmp(element, 'relativistic_mass')) THEN
      elementselected = c_dump_part_rel_mass

    ELSE IF (str_cmp(element, 'gamma')) THEN
      elementselected = c_dump_part_gamma

#ifdef PHOTONS
    ELSE IF (str_cmp(element, 'optical_depth')) THEN
      elementselected = c_dump_part_opdepth
#endif

#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    ELSE IF (str_cmp(element, 'qed_energy')) THEN
      elementselected = c_dump_part_qed_energy
#endif

#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    ELSE IF (str_cmp(element, 'trident_optical_depth')) THEN
      elementselected = c_dump_part_opdepth_tri
#endif

#ifdef BREMSSTRAHLUNG
    ELSE IF (str_cmp(element, 'bremsstrahlung_optical_depth')) THEN
      elementselected = c_dump_part_opdepth_brem
#endif

#ifdef WORK_DONE_INTEGRATED
    ELSE IF (str_cmp(element, 'work_x')) THEN
      elementselected = c_dump_part_work_x

    ELSE IF (str_cmp(element, 'work_y')) THEN
      elementselected = c_dump_part_work_y

    ELSE IF (str_cmp(element, 'work_z')) THEN
      elementselected = c_dump_part_work_z

    ELSE IF (str_cmp(element, 'work_x_total')) THEN
      elementselected = c_dump_part_work_x_total

    ELSE IF (str_cmp(element, 'work_y_total')) THEN
      elementselected = c_dump_part_work_y_total

    ELSE IF (str_cmp(element, 'work_z_total')) THEN
      elementselected = c_dump_part_work_z_total
#endif

#ifdef TRANSITION_RATES
    ELSE IF (str_cmp(element, 'rate_fi')) THEN
      elementselected = c_dump_part_rate_fi

    ELSE IF (str_cmp(element, 'rate_ci')) THEN
      elementselected = c_dump_part_rate_ci

    ELSE IF (str_cmp(element, 'rate_dr')) THEN
      elementselected = c_dump_part_rate_dr

    ELSE IF (str_cmp(element, 'rate_rr')) THEN
      elementselected = c_dump_part_rate_rr

    ELSE IF (str_cmp(element, 'rate_3br')) THEN
      elementselected = c_dump_part_rate_3br
#endif

    ELSE IF (str_cmp(element, 'ex')) THEN
      elementselected = c_dump_ex

    ELSE IF (str_cmp(element, 'ey')) THEN
      elementselected = c_dump_ey

    ELSE IF (str_cmp(element, 'ez')) THEN
      elementselected = c_dump_ez

    ELSE IF (str_cmp(element, 'bx')) THEN
      elementselected = c_dump_bx

    ELSE IF (str_cmp(element, 'by')) THEN
      elementselected = c_dump_by

    ELSE IF (str_cmp(element, 'bz')) THEN
      elementselected = c_dump_bz

    ELSE IF (str_cmp(element, 'jx')) THEN
      elementselected = c_dump_jx

    ELSE IF (str_cmp(element, 'jy')) THEN
      elementselected = c_dump_jy

    ELSE IF (str_cmp(element, 'jz')) THEN
      elementselected = c_dump_jz

    ELSE IF (str_cmp(element, 'ekbar') &
        .OR. str_cmp(element, 'average_particle_energy')) THEN
      elementselected = c_dump_ekbar

    ELSE IF (str_cmp(element, 'mass_density')) THEN
      elementselected = c_dump_mass_density

    ELSE IF (str_cmp(element, 'charge_density')) THEN
      elementselected = c_dump_charge_density

    ELSE IF (str_cmp(element, 'number_density')) THEN
      elementselected = c_dump_number_density

    ELSE IF (str_cmp(element, 'ppc') &
        .OR. str_cmp(element, 'particles_per_cell')) THEN
      elementselected = c_dump_ppc

    ELSE IF (str_cmp(element, 'average_weight')) THEN
      elementselected = c_dump_average_weight

    ELSE IF (str_cmp(element, 'average_px')) THEN
      elementselected = c_dump_average_px

    ELSE IF (str_cmp(element, 'average_py')) THEN
      elementselected = c_dump_average_py

    ELSE IF (str_cmp(element, 'average_pz')) THEN
      elementselected = c_dump_average_pz

    ELSE IF (str_cmp(element, 'temperature')) THEN
      elementselected = c_dump_temperature

    ELSE IF (str_cmp(element, 'tx') &
        .OR. str_cmp(element, 'temperature_x') &
        .OR. str_cmp(element, 'temp_x')) THEN
      elementselected = c_dump_temperature_x

    ELSE IF (str_cmp(element, 'ty') &
        .OR. str_cmp(element, 'temperature_y') &
        .OR. str_cmp(element, 'temp_y')) THEN
      elementselected = c_dump_temperature_y

    ELSE IF (str_cmp(element, 'tz') &
        .OR. str_cmp(element, 'temperature_z') &
        .OR. str_cmp(element, 'temp_z')) THEN
      elementselected = c_dump_temperature_z

    ELSE IF (str_cmp(element, 'distribution_functions')) THEN
      elementselected = c_dump_dist_fns

    ELSE IF (str_cmp(element, 'particle_probes')) THEN
      elementselected = c_dump_probes

    ELSE IF (str_cmp(element, 'ejected_particles')) THEN
      elementselected = c_dump_ejected_particles

    ELSE IF (str_cmp(element, 'ekflux') &
        .OR. str_cmp(element, 'particle_energy_flux')) THEN
      elementselected = c_dump_ekflux

    ELSE IF (str_cmp(element, 'poynt_flux') &
        .OR. str_cmp(element, 'poynting_flux')) THEN
      elementselected = c_dump_poynt_flux

    ELSE IF (str_cmp(element, 'cpml_psi_eyx')) THEN
      elementselected = c_dump_cpml_psi_eyx

    ELSE IF (str_cmp(element, 'cpml_psi_ezx')) THEN
      elementselected = c_dump_cpml_psi_ezx

    ELSE IF (str_cmp(element, 'cpml_psi_byx')) THEN
      elementselected = c_dump_cpml_psi_byx

    ELSE IF (str_cmp(element, 'cpml_psi_bzx')) THEN
      elementselected = c_dump_cpml_psi_bzx

    ELSE IF (str_cmp(element, 'cpml_psi_exy')) THEN
      elementselected = c_dump_cpml_psi_exy

    ELSE IF (str_cmp(element, 'cpml_psi_ezy')) THEN
      elementselected = c_dump_cpml_psi_ezy

    ELSE IF (str_cmp(element, 'cpml_psi_bxy')) THEN
      elementselected = c_dump_cpml_psi_bxy

    ELSE IF (str_cmp(element, 'cpml_psi_bzy')) THEN
      elementselected = c_dump_cpml_psi_bzy

    ELSE IF (str_cmp(element, 'cpml_psi_exz')) THEN
      elementselected = c_dump_cpml_psi_exz

    ELSE IF (str_cmp(element, 'cpml_psi_eyz')) THEN
      elementselected = c_dump_cpml_psi_eyz

    ELSE IF (str_cmp(element, 'cpml_psi_bxz')) THEN
      elementselected = c_dump_cpml_psi_bxz

    ELSE IF (str_cmp(element, 'cpml_psi_byz')) THEN
      elementselected = c_dump_cpml_psi_byz

    ELSE IF (str_cmp(element, 'absorption')) THEN
      elementselected = c_dump_absorption

    ELSE IF (str_cmp(element, 'total_energy_sum')) THEN
      elementselected = c_dump_total_energy_sum

    ELSE IF (str_cmp(element, 'coulomb_logarithm')) THEN
      elementselected = c_dump_cou_log

    ELSE
      got_element = .FALSE.

    END IF

    IF (style_error == c_err_old_style_ignore) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "' // TRIM(ADJUSTL(element)) &
              // '" not ', 'allowed in an unnamed output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        END DO
      END IF
    ELSE IF (style_error == c_err_new_style_ignore) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "' // TRIM(ADJUSTL(element)) &
              // '" not ', 'allowed in a named output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        END DO
      END IF
    ELSE IF (style_error == c_err_new_style_global) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "' // TRIM(ADJUSTL(element)) &
              // '" should be moved to ', 'an "output_global" block.'
          WRITE(io,*) 'Its value will be applied to all output blocks.'
          WRITE(io,*)
        END DO
      END IF
    END IF

    IF (elementselected == 0) THEN
      IF (got_element) RETURN

      errcode = c_err_unknown_element
      RETURN
    END IF

    IF (io_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    END IF
    io_block_done(elementselected) = .TRUE.
    errcode = c_err_none
    style_error = c_err_none

    mask_element = elementselected
    ALLOCATE(subsets(n_subsets+1))
    CALL as_list(value, subsets, n_list, errcode)

    DO is = 1, n_list
      subset = subsets(is)
      IF (is == 1) THEN
        mask = subset
        fullmask = mask
      ELSE
        mask = IOR(subset_list(subset)%mask,fullmask)
      END IF

      ! If setting dumpmask for features which haven't been compiled
      ! in then issue a warning
#ifdef NO_PARTICLE_PROBES
      IF (mask_element == c_dump_probes &
          .AND. mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-NO_DPARTICLE_PROBES'
        mask = c_io_never
      END IF
#endif

#ifdef PARTICLE_ID4
#define PARTICLE_ID
#endif
#ifndef PARTICLE_ID
      IF (mask_element == c_dump_part_id &
          .AND. mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) THEN
        errcode = c_err_pp_options_missing
        extended_error_string = '-DPARTICLE_ID'
        mask = c_io_never
      END IF
#endif

      ! Setting some flags like species and average
      ! wastes memory if the parameters make no sense. Do sanity checking.

      IF (IAND(mask, c_io_species) /= 0 &
          .OR. IAND(mask, c_io_no_sum) == 0) THEN
        bad = .TRUE.
        ! Check for sensible per species variables
        IF (mask_element == c_dump_ekbar) bad = .FALSE.
        IF (mask_element == c_dump_ekflux) bad = .FALSE.
        IF (mask_element == c_dump_mass_density) bad = .FALSE.
        IF (mask_element == c_dump_charge_density) bad = .FALSE.
        IF (mask_element == c_dump_number_density) bad = .FALSE.
        IF (mask_element == c_dump_ppc) bad = .FALSE.
        IF (mask_element == c_dump_average_weight) bad = .FALSE.
        IF (mask_element == c_dump_average_px) bad = .FALSE.
        IF (mask_element == c_dump_average_py) bad = .FALSE.
        IF (mask_element == c_dump_average_pz) bad = .FALSE.
        IF (mask_element == c_dump_temperature) bad = .FALSE.
        IF (mask_element == c_dump_temperature_x) bad = .FALSE.
        IF (mask_element == c_dump_temperature_y) bad = .FALSE.
        IF (mask_element == c_dump_temperature_z) bad = .FALSE.
        IF (mask_element == c_dump_jx) bad = .FALSE.
        IF (mask_element == c_dump_jy) bad = .FALSE.
        IF (mask_element == c_dump_jz) bad = .FALSE.
        IF (mask_element == c_dump_total_energy_sum) bad = .FALSE.
        IF (mask_element == c_dump_cou_log) bad = .FALSE.
        IF (bad) THEN
          IF (rank == 0 .AND. IAND(mask, c_io_species) /= 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Attempting to set per species property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            END DO
          END IF
          mask = IAND(mask, NOT(c_io_species))
          mask = IOR(mask, c_io_no_sum)
          IF (mask_element == c_dump_poynt_flux) &
              mask = IAND(mask, NOT(c_io_no_sum))
        END IF
      END IF

      IF (IAND(mask, c_io_averaged) /= 0) THEN
        bad = .TRUE.
        ! Check for sensible averaged variables
        IF (mask_element == c_dump_ex) bad = .FALSE.
        IF (mask_element == c_dump_ey) bad = .FALSE.
        IF (mask_element == c_dump_ez) bad = .FALSE.
        IF (mask_element == c_dump_bx) bad = .FALSE.
        IF (mask_element == c_dump_by) bad = .FALSE.
        IF (mask_element == c_dump_bz) bad = .FALSE.
        IF (mask_element == c_dump_jx) bad = .FALSE.
        IF (mask_element == c_dump_jy) bad = .FALSE.
        IF (mask_element == c_dump_jz) bad = .FALSE.
        IF (mask_element == c_dump_poynt_flux) bad = .FALSE.

        ! Unset 'no_sum' dumpmask for grid variables
        IF (.NOT.bad) mask = IAND(mask, NOT(c_io_no_sum))

        IF (mask_element == c_dump_ekbar) bad = .FALSE.
        IF (mask_element == c_dump_mass_density) bad = .FALSE.
        IF (mask_element == c_dump_charge_density) bad = .FALSE.
        IF (mask_element == c_dump_number_density) bad = .FALSE.
        IF (mask_element == c_dump_ppc) bad = .FALSE.
        IF (mask_element == c_dump_average_weight) bad = .FALSE.
        IF (mask_element == c_dump_average_px) bad = .FALSE.
        IF (mask_element == c_dump_average_py) bad = .FALSE.
        IF (mask_element == c_dump_average_pz) bad = .FALSE.
        IF (mask_element == c_dump_temperature) bad = .FALSE.
        IF (mask_element == c_dump_temperature_x) bad = .FALSE.
        IF (mask_element == c_dump_temperature_y) bad = .FALSE.
        IF (mask_element == c_dump_temperature_z) bad = .FALSE.
        IF (mask_element == c_dump_ekflux) bad = .FALSE.
        IF (mask_element == c_dump_cou_log) bad = .FALSE.
        IF (bad) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Attempting to set average property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            END DO
          END IF
          mask = IAND(mask, NOT(c_io_averaged))
        ELSE
          any_average = .TRUE.
          io_block%any_average = .TRUE.
          IF (IAND(mask, c_io_average_single) /= 0 .AND. num /= r4) THEN
            io_block%averaged_data(mask_element)%dump_single = .TRUE.
          END IF
          i = averaged_var_block(mask_element)
          IF (i /= 0 .AND. i /= block_number) THEN
            IF (rank == 0 .AND. .NOT.warning_printed) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
                WRITE(io,*) 'Error occurred whilst assigning the averaging ', &
                    'dumpmask of the variable'
                WRITE(io,*) '"' // TRIM(element) &
                    // '" in ', 'output block number ', block_number
                WRITE(io,*) 'Only one average per variable can be computed'
                WRITE(io,*) 'If multiple were specified the first averaging ', &
                    'time will be used'
                WRITE(io,*)
              END DO
              warning_printed = .TRUE.
            END IF
          ELSE
            averaged_var_block(mask_element) = block_number
          END IF
        END IF
      END IF

      IF (is == 1) THEN
        io_block%dumpmask(mask_element) = mask
      ELSE
        subset_list(subset)%dumpmask(block_number,mask_element) = mask
      END IF
    END DO

    DEALLOCATE(subsets)

  END FUNCTION io_block_handle_element



  FUNCTION io_block_check() RESULT(errcode)

    INTEGER :: errcode, io, iu

    ! Just assume that anything not included except for the compulsory
    ! elements is not wanted
    errcode = c_err_none

    ! Averaging info not compulsory unless averaged variable selected
    IF (any_average .AND. .NOT. got_dt_average) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Required output block element "dt_average"', &
              ' absent. Please create this entry in the input deck'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    ! Can't check the io_block if it hasn't been allocated.
    IF (n_io_blocks == 0) RETURN

    IF (io_block%dt_average > io_block%dt_snapshot) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Averaging time is longer than dt_snapshot, will set', &
              ' averaging time equal'
          WRITE(io,*) 'to dt_snapshot.'
        END DO
      END IF
      io_block%dt_average = io_block%dt_snapshot
    END IF

  END FUNCTION io_block_check



  SUBROUTINE init_io_block(io_block)

    TYPE(io_block_type) :: io_block
    INTEGER :: i

    io_block%name = ''
    io_block%dt_snapshot = -1.0_num
    io_block%time_prev = 0.0_num
    io_block%time_first = 0.0_num
    io_block%dt_average = -1.0_num
    io_block%dt_min_average = -1.0_num
    io_block%average_time = -1.0_num
    io_block%average_time_start = -1.0_num
    io_block%nstep_snapshot = -1
    io_block%nstep_prev = 0
    io_block%nstep_first = 0
    io_block%nstep_average = -1
    io_block%restart = .FALSE.
    io_block%dump = .FALSE.
    io_block%any_average = .FALSE.
    io_block%dump_first = .TRUE.
    io_block%dump_last = .TRUE.
    io_block%dump_source_code = .FALSE.
    io_block%dump_input_decks = .FALSE.
    io_block%dump_first_after_restart = .FALSE.
    io_block%dumpmask = c_io_none
    io_block%time_start = -1.0_num
    io_block%time_stop  = HUGE(1.0_num)
    io_block%nstep_start = -1
    io_block%nstep_stop  = HUGE(1)
    io_block%dump_cycle  = HUGE(1)
    io_block%dump_cycle_first_index = 0
    io_block%prefix_index = 1
    io_block%rolling_restart = .FALSE.
    io_block%disabled = .FALSE.
    io_block%use_offset_grid = .FALSE.
    io_block%walltime_interval = -1.0_num
    io_block%walltime_prev = 0.0_num
    io_block%walltime_start = -1.0_num
    io_block%walltime_stop  = HUGE(1.0_num)
    NULLIFY(io_block%dump_at_nsteps)
    NULLIFY(io_block%dump_at_times)
    NULLIFY(io_block%dump_at_walltimes)
    DO i = 1, num_vars_to_dump
      io_block%averaged_data(i)%dump_single = .FALSE.
    END DO

  END SUBROUTINE init_io_block



  SUBROUTINE copy_io_block(io_block, io_block_copy)

    TYPE(io_block_type) :: io_block, io_block_copy
    INTEGER :: i

    io_block_copy = io_block
    NULLIFY(io_block%dump_at_nsteps)
    NULLIFY(io_block%dump_at_times)
    NULLIFY(io_block%dump_at_walltimes)
    DO i = 1, num_vars_to_dump
      io_block_copy%averaged_data(i) = io_block%averaged_data(i)
    END DO

  END SUBROUTINE copy_io_block



  SUBROUTINE set_restart_dumpmasks

    ! Set the dumpmask for variables required to restart

    ! Particles
    io_block%dumpmask(c_dump_part_grid) = &
        IOR(io_block%dumpmask(c_dump_part_grid), c_io_restartable)
    io_block%dumpmask(c_dump_part_species) = &
        IOR(io_block%dumpmask(c_dump_part_species), c_io_restartable)
    io_block%dumpmask(c_dump_part_weight) = &
        IOR(io_block%dumpmask(c_dump_part_weight), c_io_restartable)
    io_block%dumpmask(c_dump_part_px) = &
        IOR(io_block%dumpmask(c_dump_part_px), c_io_restartable)
    io_block%dumpmask(c_dump_part_py) = &
        IOR(io_block%dumpmask(c_dump_part_py), c_io_restartable)
    io_block%dumpmask(c_dump_part_pz) = &
        IOR(io_block%dumpmask(c_dump_part_pz), c_io_restartable)
#ifdef PHOTONS
    io_block%dumpmask(c_dump_part_opdepth) = &
        IOR(io_block%dumpmask(c_dump_part_opdepth), c_io_restartable)
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    io_block%dumpmask(c_dump_part_qed_energy) = &
        IOR(io_block%dumpmask(c_dump_part_qed_energy), c_io_restartable)
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    io_block%dumpmask(c_dump_part_opdepth_tri) = &
        IOR(io_block%dumpmask(c_dump_part_opdepth_tri), c_io_restartable)
#endif
#ifdef BREMSSTRAHLUNG
    io_block%dumpmask(c_dump_part_opdepth_brem) = &
        IOR(io_block%dumpmask(c_dump_part_opdepth_brem), c_io_restartable)
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    io_block%dumpmask(c_dump_part_id) = &
        IOR(io_block%dumpmask(c_dump_part_id), c_io_restartable)
#endif
#ifdef WORK_DONE_INTEGRATED
    io_block%dumpmask(c_dump_part_work_x_total) = &
         IOR(io_block%dumpmask(c_dump_part_work_x_total), c_io_restartable)
    io_block%dumpmask(c_dump_part_work_y_total) = &
         IOR(io_block%dumpmask(c_dump_part_work_y_total), c_io_restartable)
    io_block%dumpmask(c_dump_part_work_z_total) = &
         IOR(io_block%dumpmask(c_dump_part_work_z_total), c_io_restartable)
#endif   
    ! Persistent IDs
    io_block%dumpmask(c_dump_persistent_ids) = &
        IOR(io_block%dumpmask(c_dump_persistent_ids), c_io_restartable)
    ! Fields
    io_block%dumpmask(c_dump_grid) = &
        IOR(io_block%dumpmask(c_dump_grid), c_io_restartable)
    io_block%dumpmask(c_dump_ex) = &
        IOR(io_block%dumpmask(c_dump_ex), c_io_restartable)
    io_block%dumpmask(c_dump_ey) = &
        IOR(io_block%dumpmask(c_dump_ey), c_io_restartable)
    io_block%dumpmask(c_dump_ez) = &
        IOR(io_block%dumpmask(c_dump_ez), c_io_restartable)
    io_block%dumpmask(c_dump_bx) = &
        IOR(io_block%dumpmask(c_dump_bx), c_io_restartable)
    io_block%dumpmask(c_dump_by) = &
        IOR(io_block%dumpmask(c_dump_by), c_io_restartable)
    io_block%dumpmask(c_dump_bz) = &
        IOR(io_block%dumpmask(c_dump_bz), c_io_restartable)
    io_block%dumpmask(c_dump_jx) = &
        IOR(io_block%dumpmask(c_dump_jx), c_io_restartable)
    io_block%dumpmask(c_dump_jy) = &
        IOR(io_block%dumpmask(c_dump_jy), c_io_restartable)
    io_block%dumpmask(c_dump_jz) = &
        IOR(io_block%dumpmask(c_dump_jz), c_io_restartable)
    ! CPML boundaries
    io_block%dumpmask(c_dump_cpml_psi_eyx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_eyx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_ezx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_ezx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_byx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_byx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bzx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bzx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_exy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_exy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_ezy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_ezy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bxy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bxy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bzy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bzy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_exz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_exz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_eyz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_eyz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bxz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bxz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_byz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_byz), c_io_restartable)

  END SUBROUTINE set_restart_dumpmasks

END MODULE deck_io_block
