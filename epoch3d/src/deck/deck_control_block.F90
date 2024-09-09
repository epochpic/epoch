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

MODULE deck_control_block

  USE strings_advanced
  USE fields
  USE timer
  USE sdf

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: control_deck_initialise, control_deck_finalise
  PUBLIC :: control_block_start, control_block_end
  PUBLIC :: control_block_handle_element, control_block_check

  LOGICAL :: got_time, got_grid(2*c_ndims)
  LOGICAL :: got_optimal_layout, got_nproc

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state /= c_ds_first) RETURN

    use_exact_restart = .TRUE.
    use_exact_restart_set = .FALSE.
    allow_cpu_reduce = .TRUE.
    check_walltime = .FALSE.
    simplify_deck = .TRUE.
    print_deck_constants = .FALSE.
    allow_missing_restart = .FALSE.
    print_eta_string = .TRUE.
    use_current_correction = .FALSE.
    use_particle_count_update = .FALSE.
    use_accurate_n_zeros = .FALSE.
    reset_walltime = .FALSE.
    balance_first = .TRUE.
    ic_from_restart = .FALSE.
    neutral_background = .TRUE.
    use_particle_migration = .FALSE.
    use_pre_balance = .TRUE.
    use_optimal_layout = .TRUE.
    got_optimal_layout = .FALSE.
    got_nproc = .FALSE.
    restart_number = 0
    check_stop_frequency = 10
    stop_at_walltime = -1.0_num
    restart_filename = ''
    n_zeros_control = -1
    dlb_maximum_interval = 500
    dlb_force_interval = 2000
    dlb_threshold = -1.0_num
    nx_global = -1
    ny_global = -1
    nz_global = -1
    particle_push_start_time = 0.0_num
    particle_migration_interval = 1
    maxwell_solver = c_maxwell_solver_yee
    got_grid(:) = .FALSE.
    got_time = .FALSE.
    physics_table_location = 'src/physics_packages/TABLES'

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

    CHARACTER(LEN=22) :: filename_fmt, str
    INTEGER :: io, iu
    LOGICAL, SAVE :: warn = .TRUE.
    LOGICAL :: exists

    IF (n_zeros_control > 0) THEN
      IF (n_zeros_control < n_zeros) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(n_zeros, str)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'n_zeros was less than ', TRIM(str), &
                        ' and has been ignored'
          END DO
        END IF
        n_zeros_control = -1
      ELSE
        n_zeros = n_zeros_control
      END IF
    END IF

    IF (ic_from_restart) THEN
      IF (TRIM(restart_filename) == '') THEN
        WRITE(filename_fmt, '(''(i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
            n_zeros, n_zeros
        WRITE(restart_filename, filename_fmt) restart_number
      END IF
      full_restart_filename = TRIM(filesystem) &
          // TRIM(data_dir) // '/' // TRIM(restart_filename)

      CALL check_valid_restart
    END IF

    IF (maxwell_solver == c_maxwell_solver_lehe &
        .OR. maxwell_solver == c_maxwell_solver_lehe_x &
        .OR. maxwell_solver == c_maxwell_solver_lehe_y &
        .OR. maxwell_solver == c_maxwell_solver_lehe_z) THEN
      fng = 2
      IF (maxwell_solver == c_maxwell_solver_lehe) THEN
        maxwell_solver = c_maxwell_solver_lehe_x
        IF (rank == 0 .AND. deck_state == c_ds_first) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Using Lehe solver optimised for the x-direction'
            WRITE(io,*)
          END DO
        END IF
      END IF
    END IF

    IF (.NOT.ic_from_restart) use_exact_restart = .FALSE.

    IF (rank == 0) THEN
      IF (nx_global < 1) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The mandatory parameter "nx" has not been specified.'
        END DO
        CALL abort_code(c_err_missing_elements)
      END IF
      IF (ny_global < 1) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The mandatory parameter "ny" has not been specified.'
        END DO
        CALL abort_code(c_err_missing_elements)
      END IF
      IF (nz_global < 1) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The mandatory parameter "nz" has not been specified.'
        END DO
        CALL abort_code(c_err_missing_elements)
      END IF
    END IF

    IF (got_nproc .AND. got_optimal_layout) THEN
      IF (warn .AND. rank == 0) THEN
        warn = .FALSE.
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Both "use_optimal_layout" and "nprocx/y/z" specified ', &
                      'in the input deck.'
          WRITE(io,*) 'The specified processor layout will be ignored.'
          WRITE(io,*)
        END DO
      END IF
    ELSE IF (got_nproc .AND. .NOT.got_optimal_layout) THEN
      use_optimal_layout = .FALSE.
    END IF

    IF (deck_state == c_ds_first) RETURN

    IF (stop_at_walltime >= 0.0_num) THEN
      check_walltime = .TRUE.
      timer_collect = .TRUE.
    END IF

    ! use_balance only if threshold is positive
    IF (dlb_threshold > 0) use_balance = .TRUE.
    IF (dlb_maximum_interval < 1) dlb_maximum_interval = HUGE(1)
    IF (dlb_force_interval < 1) dlb_force_interval = HUGE(1)

    ! Check physics tables are in the correct place
    IF (rank == 0) THEN
      INQUIRE(file=TRIM(physics_table_location) // &
          '/ionisation_energies.table', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Unable to find ionisation_energies.table in ', &
              'directory "' // TRIM(physics_table_location) // '"'
          WRITE(io,*) 'Either tables have been modified, or the path to the ', &
              'physics tables must be set'
          WRITE(io,*) 'Use key: physics_table_location in the control block.'
        END DO
      END IF
    END IF

    IF (use_field_ionisation) need_random_state = .TRUE.

  END SUBROUTINE control_deck_finalise



  SUBROUTINE control_block_start

  END SUBROUTINE control_block_start



  SUBROUTINE control_block_end

  END SUBROUTINE control_block_end



  FUNCTION control_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    CHARACTER(LEN=string_length) :: str_tmp
    CHARACTER(LEN=1) :: c
    INTEGER :: errcode
    INTEGER :: field_order, ierr, io, iu, i
    LOGICAL :: isnum
    INTEGER, DIMENSION(:), POINTER :: stride_temp

    errcode = c_err_none

    IF (deck_state /= c_ds_first) RETURN

    IF (str_cmp(element, 'nx')) THEN
      nx_global = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'ny')) THEN
      ny_global = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'nz')) THEN
      nz_global = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'x_min') &
        .OR. str_cmp(element, 'x_start')) THEN
      x_min = as_real_print(value, element, errcode)
      got_grid(1) = .TRUE.

    ELSE IF (str_cmp(element, 'x_max') &
        .OR. str_cmp(element, 'x_end')) THEN
      x_max = as_real_print(value, element, errcode)
      got_grid(2) = .TRUE.

    ELSE IF (str_cmp(element, 'y_min') &
        .OR. str_cmp(element, 'y_start')) THEN
      y_min = as_real_print(value, element, errcode)
      got_grid(3) = .TRUE.

    ELSE IF (str_cmp(element, 'y_max') &
        .OR. str_cmp(element, 'y_end')) THEN
      y_max = as_real_print(value, element, errcode)
      got_grid(4) = .TRUE.

    ELSE IF (str_cmp(element, 'z_min') &
        .OR. str_cmp(element, 'z_start')) THEN
      z_min = as_real_print(value, element, errcode)
      got_grid(5) = .TRUE.

    ELSE IF (str_cmp(element, 'z_max') &
        .OR. str_cmp(element, 'z_end')) THEN
      z_max = as_real_print(value, element, errcode)
      got_grid(6) = .TRUE.

    ELSE IF (str_cmp(element, 'nprocx')) THEN
      nprocx = as_integer_print(value, element, errcode)
      got_nproc = .TRUE.

    ELSE IF (str_cmp(element, 'nprocy')) THEN
      nprocy = as_integer_print(value, element, errcode)
      got_nproc = .TRUE.

    ELSE IF (str_cmp(element, 'nprocz')) THEN
      nprocz = as_integer_print(value, element, errcode)
      got_nproc = .TRUE.

    ELSE IF (str_cmp(element, 'npart') &
        .OR. str_cmp(element, 'nparticles')) THEN
      npart_global = as_long_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'nsteps')) THEN
      nsteps = as_integer_print(value, element, errcode)
      got_time = .TRUE.

    ELSE IF (str_cmp(element, 't_end')) THEN
      t_end = as_real_print(value, element, errcode)
      got_time = .TRUE.

    ELSE IF (str_cmp(element, 'dt_multiplier')) THEN
      dt_multiplier = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dlb_threshold') &
        .OR. str_cmp(element, 'balance_threshold')) THEN
      dlb_threshold = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'icfile')) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)

    ELSE IF (str_cmp(element, 'restart_snapshot')) THEN
      isnum = .TRUE.
      str_tmp = TRIM(ADJUSTL(value))
      DO i = 1,LEN_TRIM(str_tmp)
        c = str_tmp(i:i)
        IF (c < '0' .OR. c > '9') THEN
          isnum = .FALSE.
          EXIT
        END IF
      END DO
      IF (isnum) THEN
        restart_number = as_integer_print(value, element, errcode)
      ELSE
        restart_filename = TRIM(str_tmp)
      END IF
      ic_from_restart = .TRUE.

    ELSE IF (str_cmp(element, 'neutral_background')) THEN
      neutral_background = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'field_order')) THEN
      field_order = as_integer_print(value, element, errcode)
      IF (field_order /= 2 .AND. field_order /= 4 &
          .AND. field_order /= 6) THEN
        errcode = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      END IF

    ELSE IF (str_cmp(element, 'stdout_frequency')) THEN
      stdout_frequency = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_random_seed')) THEN
      use_random_seed = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'smooth_currents')) THEN
      smooth_currents = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'field_ionisation') &
        .OR. str_cmp(element, 'field_ionization') &
        .OR. str_cmp(element, 'use_field_ionise')) THEN
      use_field_ionisation = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'physics_table_location')) THEN
      physics_table_location = TRIM(ADJUSTL(value))

    ELSE IF (str_cmp(element, 'use_multiphoton') &
        .OR. str_cmp(element, 'multiphoton')) THEN
      use_multiphoton = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_bsi') &
        .OR. str_cmp(element, 'bsi')) THEN
      use_bsi = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_radiative_recombination')) THEN
      use_radiative_recombination = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_dielectronic_recombination')) THEN
      use_dielectronic_recombination = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_three_body_recombination')) THEN
      use_three_body_recombination = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'particle_tstart')) THEN
      particle_push_start_time = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_migration') &
        .OR. str_cmp(element, 'migrate_particles')) THEN
      use_particle_migration = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'migration_interval')) THEN
      particle_migration_interval = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_exact_restart')) THEN
      use_exact_restart = as_logical_print(value, element, errcode)
      use_exact_restart_set = use_exact_restart

    ELSE IF (str_cmp(element, 'allow_cpu_reduce')) THEN
      allow_cpu_reduce = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'check_stop_file_frequency') &
        .OR. str_cmp(element, 'check_stop_frequency')) THEN
      check_stop_frequency = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'stop_at_walltime')) THEN
      stop_at_walltime = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'stop_at_walltime_file')) THEN
      IF (rank == 0) THEN
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(value))
        IF (ierr == 0) THEN
          READ(lu,*,iostat=ierr) stop_at_walltime
          CLOSE(lu)
        END IF
      END IF
      CALL MPI_BCAST(stop_at_walltime, 1, mpireal, 0, comm, errcode)

    ELSE IF (str_cmp(element, 'simplify_deck')) THEN
      simplify_deck = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'print_constants')) THEN
      print_deck_constants = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'allow_missing_restart')) THEN
      allow_missing_restart = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'print_eta_string')) THEN
      print_eta_string = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'n_zeros')) THEN
      n_zeros_control = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_current_correction')) THEN
      use_current_correction = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'maxwell_solver')) THEN
      maxwell_solver = as_integer_print(value, element, errcode)
      IF (maxwell_solver /= c_maxwell_solver_yee &
          .AND. maxwell_solver /= c_maxwell_solver_lehe &
          .AND. maxwell_solver /= c_maxwell_solver_lehe_x &
          .AND. maxwell_solver /= c_maxwell_solver_lehe_y &
          .AND. maxwell_solver /= c_maxwell_solver_lehe_z &
          .AND. maxwell_solver /= c_maxwell_solver_cowan &
          .AND. maxwell_solver /= c_maxwell_solver_pukhov &
          .AND. maxwell_solver /= c_maxwell_solver_custom) THEN
        errcode = c_err_bad_value
      END IF

    ELSE IF (str_cmp(element, 'use_particle_count_update')) THEN
      use_particle_count_update = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_accurate_n_zeros')) THEN
      use_accurate_n_zeros = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'reset_walltime')) THEN
      reset_walltime = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dlb_maximum_interval') &
        .OR. str_cmp(element, 'balance_maximum_interval')) THEN
      dlb_maximum_interval = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'dlb_force_interval') &
        .OR. str_cmp(element, 'balance_force_interval')) THEN
      dlb_force_interval = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'balance_first')) THEN
      balance_first = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_pre_balance') &
        .OR. str_cmp(element, 'pre_balance')) THEN
      use_pre_balance = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'use_optimal_layout') &
        .OR. str_cmp(element, 'optimal_layout')) THEN
      use_optimal_layout = as_logical_print(value, element, errcode)
      got_optimal_layout = use_optimal_layout

    ELSE IF (str_cmp(element, 'smooth_iterations')) THEN
      smooth_its = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'recombine_n_step')) THEN
      recombine_n_step = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'smooth_compensation')) THEN
      IF (as_logical_print(value, element, errcode)) smooth_comp_its = 1

    ELSE IF (str_cmp(element, 'smooth_strides')) THEN
      IF (str_cmp(value, 'auto')) THEN
        ALLOCATE(smooth_strides(4), SOURCE=[1,2,3,4])
        RETURN
      END IF
      stride_temp => NULL()
      CALL get_allocated_array(value, stride_temp, errcode)
      ALLOCATE(smooth_strides(SIZE(stride_temp)), SOURCE=stride_temp)
      DEALLOCATE(stride_temp)
      sng = MAXVAL(smooth_strides)

    ELSE IF (str_cmp(element, 'use_more_setup_memory')) THEN
      use_more_setup_memory = as_logical_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'deck_warnings_fatal')) THEN
      all_deck_errcodes_fatal = as_logical_print(value, element, errcode)

    ELSE
      errcode = c_err_unknown_element

    END IF

  END FUNCTION control_block_handle_element



  FUNCTION control_block_check() RESULT(errcode)

    INTEGER :: errcode, idx, io, iu
    CHARACTER(LEN=5), DIMENSION(6), PARAMETER :: &
        grid_name = ['x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']

    errcode = c_err_none

    IF (.NOT. got_time) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Required control block element "t_end" or "nsteps"', &
              ' absent. Please create this entry in the input deck'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    DO idx = 1, 2 * c_ndims
      IF (.NOT. got_grid(idx)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required control block element "' // grid_name(idx) &
                // '" absent. Please create this entry in the input deck'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
    END DO

    IF (.NOT. neutral_background) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The option "neutral_background=F" is not supported', &
              ' in this version of EPOCH.'
        END DO
      END IF
      errcode = c_err_terminate
    END IF

    IF (maxwell_solver /= c_maxwell_solver_yee &
        .AND. field_order /= 2) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'For "field_order" > 2 only "maxwell_solver = yee"', &
              ' is supported in this version of EPOCH.'
        END DO
      END IF
      errcode = c_err_terminate
    END IF

  END FUNCTION control_block_check



  SUBROUTINE check_valid_restart

    CHARACTER(LEN=c_id_length) :: code_name
    INTEGER :: step, code_io_version, string_len
    REAL(num) :: time
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle
    LOGICAL :: valid = .TRUE.

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read, &
                  handle_errors=.FALSE.)

    IF (sdf_handle%error_code == 0) THEN
      CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
          string_len, restart_flag)
      CALL sdf_close(sdf_handle)
    ELSE
      restart_flag = .FALSE.
    END IF

    IF (.NOT. restart_flag) THEN
      valid = .FALSE.
      IF (.NOT. allow_missing_restart) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SDF file ', TRIM(full_restart_filename), &
              ' is not a restart dump. Unable to continue.'
        END IF
        CALL abort_code(c_err_io_error)
        STOP
      END IF
    END IF

    IF (valid .AND. .NOT.str_cmp(code_name, 'Epoch3d')) THEN
      valid = .FALSE.
      IF (.NOT. allow_missing_restart) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SDF restart file was not generated by Epoch3d. Unable to ', &
              'continue.'
        END IF
        CALL abort_code(c_err_io_error)
        STOP
      END IF
    END IF

    IF (.NOT. valid) THEN
      ic_from_restart = .FALSE.
      IF (rank == 0) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'No valid restart dump found. Using initial conditions instead.'
      END IF
    END IF

  END SUBROUTINE check_valid_restart

END MODULE deck_control_block
