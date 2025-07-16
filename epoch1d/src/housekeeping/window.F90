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

MODULE window

  USE boundary
  USE partlist
  USE evaluator

  IMPLICIT NONE

  REAL(num) :: density, temperature(3), drift(3)
  REAL(num), SAVE :: window_shift_fraction

CONTAINS

  SUBROUTINE initialise_window

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    window_started = .FALSE.
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

    IF (use_window_stack) CALL deallocate_stack(window_v_x_stack)

  END SUBROUTINE deallocate_window



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow, ix
    REAL(num) :: xb_min

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, window_shift_cells
      CALL insert_particles

      ! Shift the box around
      x_grid_min = x_global(1) + dx
      xb_min = xb_global(1) + dx
      x_min = xb_min + dx * cpml_thickness

      ! Setup global grid
      DO ix = 1-ng, nx_global + ng
        x_global(ix) = x_grid_min + (ix - 1) * dx
        xb_global(ix) = xb_min + (ix - 1) * dx
      END DO
      x_grid_max = x_global(nx_global)
      x_max = xb_global(nx_global+1) - dx * cpml_thickness

      CALL setup_grid_x

      CALL remove_particles

      ! Shift fields around
      CALL shift_fields
    END DO

  END SUBROUTINE shift_window



  SUBROUTINE shift_fields

    CALL shift_field(ex, ng)
    CALL shift_field(ey, ng)
    CALL shift_field(ez, ng)

    CALL shift_field(bx, ng)
    CALL shift_field(by, ng)
    CALL shift_field(bz, ng)

    CALL shift_field(jx, jng)
    CALL shift_field(jy, jng)
    CALL shift_field(jz, jng)

    IF (cpml_boundaries) THEN
      CALL shift_field(cpml_psi_eyx, ng)
      CALL shift_field(cpml_psi_ezx, ng)
      CALL shift_field(cpml_psi_byx, ng)
      CALL shift_field(cpml_psi_bzx, ng)
    END IF

    IF (x_max_boundary) THEN
      ! Fix incoming field cell.
      ex(nx)   = ex_x_max
      ex(nx+1) = ex_x_max
      ey(nx+1) = ey_x_max
      ez(nx+1) = ez_x_max
      ex(nx-1) = 0.5_num * (ex(nx-2) + ex(nx))
      ey(nx)   = 0.5_num * (ey(nx-1) + ey(nx+1))
      ez(nx)   = 0.5_num * (ez(nx-1) + ez(nx+1))
      bx(nx+1) = bx_x_max
      by(nx)   = by_x_max
      bz(nx)   = bz_x_max
      bx(nx)   = 0.5_num * (bx(nx-1) + bx(nx+1))
      by(nx-1) = 0.5_num * (by(nx-2) + by(nx))
      bz(nx-1) = 0.5_num * (bz(nx-2) + bz(nx))

      IF (cpml_boundaries) THEN
        cpml_psi_eyx(nx:nx+1) = cpml_psi_eyx(nx)
        cpml_psi_ezx(nx:nx+1) = cpml_psi_ezx(nx)
        cpml_psi_byx(nx:nx+1) = cpml_psi_byx(nx)
        cpml_psi_bzx(nx:nx+1) = cpml_psi_bzx(nx)
      END IF
    END IF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER :: i

    ! Shift field to the left by one cell
    DO i = 1-ng, nx+ng-1
      field(i) = field(i+1)
    END DO

    CALL field_bc(field, ng)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i, errcode
    INTEGER(i8) :: ipart, npart_per_cell, n_frac
    REAL(num) :: temp_local, drift_local, npart_frac
    REAL(num) :: x0, dmin, dmax, wdata
    TYPE(parameter_pack) :: parameters
    TYPE(particle_species), POINTER :: species
    REAL(num), DIMENSION(c_ndirs, 2) :: ranges

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    errcode = c_err_none

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      IF (species%background_species) CYCLE

      CALL create_empty_partlist(append_list)
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      npart_frac = species%npart_per_cell - npart_per_cell

      dmin = species%initial_conditions%density_min
      dmax = species%initial_conditions%density_max

      parameters%pack_ix = nx
      density = evaluate_with_parameters( &
          species%density_function, parameters, errcode)
      IF (density > dmax) density = dmax
      IF (density < dmin) density = 0.0_num

      IF (density < dmin) CYCLE

      DO i = 1, 3
        temperature(i) = evaluate_with_parameters( &
            species%temperature_function(i), parameters, errcode)
        drift(i) = evaluate_with_parameters( &
            species%drift_function(i), parameters, errcode)
      END DO

      ! Place extra particle based on probability
      n_frac = 0
      IF (npart_frac > 0.0_num) THEN
        IF (random() < npart_frac) n_frac = 1
      END IF

      wdata = dx / (npart_per_cell + n_frac)

      x0 = x_grid_max + 0.5_num * dx
      DO ipart = 1, npart_per_cell + n_frac
        ! Place extra particle based on probability
        IF (ipart == 0) THEN
          IF (npart_frac < random()) CYCLE
        END IF
        CALL create_particle(current)
#ifdef PER_PARTICLE_CHARGE_MASS
        current%mass = species%mass
        current%charge = species%charge
#endif
        current%part_pos = x0 + random() * dx

        IF (species%ic_df_type == c_ic_df_thermal) THEN
          DO i = 1, c_ndirs
            temp_local = temperature(i)
            drift_local = drift(i)
            current%part_p(i) = momentum_from_temperature(species%mass, &
                temp_local, drift_local)
          END DO
        ELSE IF (species%ic_df_type == c_ic_df_relativistic_thermal) THEN
          current%part_p = momentum_from_temperature_relativistic(&
              species%mass, temperature, species%fractional_tail_cutoff, drift)
        ELSE IF (species%ic_df_type == c_ic_df_arbitrary) THEN
          parameters%use_grid_position = .FALSE.
          parameters%pack_pos = current%part_pos
          errcode = c_err_none
          CALL evaluate_with_parameters_to_array(species%dist_fn_range(1), &
              parameters, 2, ranges(1,:), errcode)
          CALL evaluate_with_parameters_to_array(species%dist_fn_range(2), &
              parameters, 2, ranges(2,:), errcode)
          CALL evaluate_with_parameters_to_array(species%dist_fn_range(3), &
              parameters, 2, ranges(3,:), errcode)

          CALL sample_from_deck_expression(current, &
              species%dist_fn, parameters, ranges, species%mass, drift)
        END IF

        current%weight = density * wdata
#ifdef PARTICLE_DEBUG
        current%processor = rank
        current%processor_at_t0 = rank
#endif
        CALL add_particle_to_partlist(append_list, current)
      END DO

      CALL append_partlist(species%attached_list, append_list)
    END DO

  END SUBROUTINE insert_particles



  SUBROUTINE remove_particles

    TYPE(particle), POINTER :: current, next
    INTEGER :: ispecies

    ! Only processors on the left need do anything
    IF (x_min_boundary) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          next => current%next
          IF (current%part_pos < x_min) THEN
            CALL remove_particle_from_partlist(&
                species_list(ispecies)%attached_list, current)
            CALL destroy_particle(current)
          END IF
          current => next
        END DO
      END DO
    END IF

  END SUBROUTINE remove_particles
#endif



  SUBROUTINE moving_window

#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells, errcode = 0
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time >= window_start_time .AND. time < window_stop_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        CALL setup_boundaries
        IF (.NOT.ic_from_restart) window_shift_fraction = 0.0_num
        window_started = .TRUE.
      END IF
    END IF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
      IF (time >= window_stop_time) RETURN
      IF (use_window_stack) window_v_x = evaluate(window_v_x_stack, errcode)
      IF (window_v_x <= 0.0_num) RETURN
      window_shift_fraction = window_shift_fraction + dt * window_v_x / dx
      window_shift_cells = FLOOR(window_shift_fraction)
      ! Allow for posibility of having jumped two cells at once
      IF (window_shift_cells > 0) THEN
        window_shift_real = REAL(window_shift_cells, num)
        window_offset = window_offset + window_shift_real * dx
        CALL shift_window(window_shift_cells)
        CALL setup_bc_lists
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real
      END IF
    END IF
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE moving_window

END MODULE window
