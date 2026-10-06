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

MODULE deck_laser_block

  USE strings_advanced
  USE laser
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: laser_deck_initialise, laser_deck_finalise
  PUBLIC :: laser_block_start, laser_block_end
  PUBLIC :: laser_block_handle_element, laser_block_check

  TYPE(laser_block), POINTER :: working_laser
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  SUBROUTINE laser_deck_initialise

    n_lasers(:) = 0

  END SUBROUTINE laser_deck_initialise



  SUBROUTINE laser_deck_finalise

  END SUBROUTINE laser_deck_finalise



  SUBROUTINE laser_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.
    working_laser%use_phase_function = .TRUE.
    working_laser%use_profile_function = .TRUE.
    working_laser%use_omega_function = .FALSE.

    ! Explicitly initialise custom-profile flags and filename for each new
    ! laser block. Blank profile_data_file triggers the default filename
    ! ('temporal_spatial_profile.dat' or 'spatial_profile.dat') at load time.
    working_laser%use_custom_profile = .FALSE.
    working_laser%use_spatiotemporal = .TRUE.
    working_laser%profile_data_file = ' '

    ! Phase-from-file defaults: disabled, with a blank filename that triggers
    ! the default 'phase_profile.dat' at load time.
    working_laser%use_phase_from_file = .FALSE.
    working_laser%phase_data_file = ' '

  END SUBROUTINE laser_block_start



  SUBROUTINE laser_block_end

    IF (deck_state == c_ds_first) RETURN

    CALL attach_laser(working_laser)
    boundary_set = .FALSE.

  END SUBROUTINE laser_block_end



  FUNCTION laser_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy
    INTEGER :: io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary') .OR. str_cmp(element, 'direction')) THEN
      IF (rank == 0 .AND. str_cmp(element, 'direction')) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "direction" in the block "laser" is deprecated.'
          WRITE(io,*) 'Please use the element name "boundary" instead.'
        END DO
      END IF
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      CALL init_laser(boundary, working_laser)
      RETURN
    END IF

    IF (.NOT. boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Cannot set laser properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'amp')) THEN
      working_laser%amp = as_real_print(value, element, errcode)
      RETURN
    END IF

    ! SI (W/m^2)
    IF (str_cmp(element, 'irradiance') .OR. str_cmp(element, 'intensity')) THEN
      working_laser%amp = SQRT(as_real_print(value, element, errcode) &
          / (c*epsilon0/2.0_num))
      RETURN
    END IF

    IF (str_cmp(element, 'irradiance_w_cm2') &
        .OR. str_cmp(element, 'intensity_w_cm2')) THEN
      working_laser%amp = SQRT(as_real_print(value, element, errcode) &
          / (c*epsilon0/2.0_num)) * 100_num
      RETURN
    END IF

    IF (str_cmp(element, 'omega') .OR. str_cmp(element, 'freq')) THEN
      IF (rank == 0 .AND. str_cmp(element, 'freq')) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "freq" in the block "laser" is deprecated.'
          WRITE(io,*) 'Please use the element name "omega" instead.'
        END DO
      END IF
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_omega
      CALL laser_update_omega(working_laser)
      IF (working_laser%omega_function%is_time_varying) THEN
        working_laser%use_omega_function = .TRUE.
      ELSE
        CALL deallocate_stack(working_laser%omega_function)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'frequency')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_freq
      CALL laser_update_omega(working_laser)
      IF (working_laser%omega_function%is_time_varying) THEN
        working_laser%use_omega_function = .TRUE.
      ELSE
        CALL deallocate_stack(working_laser%omega_function)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'lambda')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_lambda
      CALL laser_update_omega(working_laser)
      IF (working_laser%omega_function%is_time_varying) THEN
        working_laser%use_omega_function = .TRUE.
      ELSE
        CALL deallocate_stack(working_laser%omega_function)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'profile')) THEN
      CALL initialise_stack(working_laser%profile_function)
      CALL tokenize(value, working_laser%profile_function, errcode)
      working_laser%profile = 0.0_num
      CALL laser_update_profile(working_laser)
      IF (working_laser%profile_function%is_time_varying) THEN
        working_laser%use_profile_function = .TRUE.
      ELSE
        CALL deallocate_stack(working_laser%profile_function)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'phase')) THEN
      CALL initialise_stack(working_laser%phase_function)
      CALL tokenize(value, working_laser%phase_function, errcode)
      working_laser%phase = 0.0_num
      CALL laser_update_phase(working_laser)
      IF (working_laser%phase_function%is_time_varying) THEN
        working_laser%use_phase_function = .TRUE.
      ELSE
        CALL deallocate_stack(working_laser%phase_function)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_laser%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_laser%t_end = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_profile')) THEN
      working_laser%use_time_function = .TRUE.
      CALL initialise_stack(working_laser%time_function)
      CALL tokenize(value, working_laser%time_function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'pol_angle') &
        .OR. str_cmp(element, 'polarisation_angle')) THEN
      working_laser%pol_angle = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'pol') &
        .OR. str_cmp(element, 'polarisation')) THEN
      ! Convert from degrees to radians
      working_laser%pol_angle = &
          pi * as_real_print(value, element, errcode) / 180.0_num
      RETURN
    END IF

    IF (str_cmp(element, 'id')) THEN
      working_laser%id = as_integer_print(value, element, errcode)
      RETURN
    END IF

    ! Custom laser profile: read the amplitude (and optionally phase) from
    ! an external binary file instead of a deck expression. See
    ! src/shared_data.F90 for the laser_block fields these set.
    IF (str_cmp(element, 'use_custom_profile')) THEN
      working_laser%use_custom_profile = as_logical_print(value, element, &
          errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_spatiotemporal_profile')) THEN
      working_laser%use_spatiotemporal = as_logical_print(value, element, &
          errcode)
      RETURN
    END IF

    ! Parse the custom profile data filename. The value can be:
    !   - A plain filename (e.g. 'TS01.dat'), resolved relative to data_dir
    !   - An absolute path (e.g. '/home/user/profiles/TS01.dat'), used as-is
    ! If omitted, the default filenames are used for backward compatibility.
    IF (str_cmp(element, 'profile_data_file')) THEN
      working_laser%profile_data_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    ! Enable reading the spatiotemporal phase from file. When true, EPOCH
    ! ignores any 'phase = ...' deck expression and instead interpolates the
    ! phase from phase_data_file at every time step (handled in
    ! laser_update_phase / custom_laser_phase).
    IF (str_cmp(element, 'use_phase_from_file')) THEN
      working_laser%use_phase_from_file = as_logical_print(value, &
          element, errcode)
      RETURN
    END IF

    ! Parse the phase data filename (same resolution rules as profile_data_file:
    ! plain filenames resolve relative to data_dir, absolute paths used as-is).
    ! If omitted, the default 'phase_profile.dat' is used.
    IF (str_cmp(element, 'phase_data_file')) THEN
      working_laser%phase_data_file = TRIM(ADJUSTL(value))
      RETURN
    END IF

    ! Shape and bounds of the spatiotemporal profile/phase binary files
    ! (use_spatiotemporal_profile = T only). These files carry no embedded
    ! header (per EPOCH's documented binary-file convention), so the grid
    ! must be declared here. The temporal extent reuses t_start/t_end rather
    ! than a separate pair of elements, since the laser is only ever active
    ! within that window anyway.
    IF (str_cmp(element, 'n_t_points') .OR. str_cmp(element, 'n_t')) THEN
      working_laser%n_t_points = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'n_transverse_points') &
        .OR. str_cmp(element, 'n_y')) THEN
      working_laser%n_transverse_points = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile_transverse_min') &
        .OR. str_cmp(element, 'y_min')) THEN
      working_laser%profile_transverse_min = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile_transverse_max') &
        .OR. str_cmp(element, 'y_max')) THEN
      working_laser%profile_transverse_max = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION laser_block_handle_element



  FUNCTION laser_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(laser_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none

    error = 0
    current => lasers
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    END DO

    IF (IAND(error, 1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "amp" or "irradiance" for every laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

  END FUNCTION laser_block_check

END MODULE deck_laser_block
