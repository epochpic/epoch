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

  REAL(num), ALLOCATABLE :: density(:,:), temperature(:,:,:), drift(:,:,:)
  REAL(num), SAVE :: window_shift_fraction

CONTAINS

  SUBROUTINE initialise_window

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    ALLOCATE(density(0:ny+1, 0:nz+1))
    ALLOCATE(temperature(0:ny+1, 0:nz+1, 1:3))
    ALLOCATE(drift(0:ny+1, 0:nz+1, 1:3))
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
    IF (ALLOCATED(density)) DEALLOCATE(density)
    IF (ALLOCATED(temperature)) DEALLOCATE(temperature)
    IF (ALLOCATED(drift)) DEALLOCATE(drift)

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

    INTEGER :: j, k

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

      CALL shift_field(cpml_psi_exy, ng)
      CALL shift_field(cpml_psi_ezy, ng)
      CALL shift_field(cpml_psi_bxy, ng)
      CALL shift_field(cpml_psi_bzy, ng)

      CALL shift_field(cpml_psi_exz, ng)
      CALL shift_field(cpml_psi_eyz, ng)
      CALL shift_field(cpml_psi_bxz, ng)
      CALL shift_field(cpml_psi_byz, ng)
    END IF

    IF (x_max_boundary) THEN
      DO k = 1-ng, nz+ng
        DO j = 1-ng, ny+ng
          ! Fix incoming field cell.
          ex(nx,j,k)   = ex_x_max(j,k)
          ex(nx+1,j,k) = ex_x_max(j,k)
          ey(nx+1,j,k) = ey_x_max(j,k)
          ez(nx+1,j,k) = ez_x_max(j,k)
          ex(nx-1,j,k) = 0.5_num * (ex(nx-2,j,k) + ex(nx,j,k))
          ey(nx,j,k)   = 0.5_num * (ey(nx-1,j,k) + ey(nx+1,j,k))
          ez(nx,j,k)   = 0.5_num * (ez(nx-1,j,k) + ez(nx+1,j,k))
          bx(nx+1,j,k) = bx_x_max(j,k)
          by(nx,j,k)   = by_x_max(j,k)
          bz(nx,j,k)   = bz_x_max(j,k)
          bx(nx,j,k)   = 0.5_num * (bx(nx-1,j,k) + bx(nx+1,j,k))
          by(nx-1,j,k) = 0.5_num * (by(nx-2,j,k) + by(nx,j,k))
          bz(nx-1,j,k) = 0.5_num * (bz(nx-2,j,k) + bz(nx,j,k))
        END DO
      END DO

      IF (cpml_boundaries) THEN
        DO k = 1-ng, nz+ng
          DO j = 1-ng, ny+ng
            cpml_psi_eyx(nx:nx+1,j,k) = cpml_psi_eyx(nx,j,k)
            cpml_psi_ezx(nx:nx+1,j,k) = cpml_psi_ezx(nx,j,k)
            cpml_psi_byx(nx:nx+1,j,k) = cpml_psi_byx(nx,j,k)
            cpml_psi_bzx(nx:nx+1,j,k) = cpml_psi_bzx(nx,j,k)

            cpml_psi_exy(nx:nx+1,j,k) = cpml_psi_exy(nx,j,k)
            cpml_psi_ezy(nx:nx+1,j,k) = cpml_psi_ezy(nx,j,k)
            cpml_psi_bxy(nx:nx+1,j,k) = cpml_psi_bxy(nx,j,k)
            cpml_psi_bzy(nx:nx+1,j,k) = cpml_psi_bzy(nx,j,k)

            cpml_psi_exz(nx:nx+1,j,k) = cpml_psi_exz(nx,j,k)
            cpml_psi_eyz(nx:nx+1,j,k) = cpml_psi_eyz(nx,j,k)
            cpml_psi_bxz(nx:nx+1,j,k) = cpml_psi_bxz(nx,j,k)
            cpml_psi_byz(nx:nx+1,j,k) = cpml_psi_byz(nx,j,k)
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, j, k

    ! Shift field to the left by one cell
    DO k = 1-ng, nz+ng
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng-1
      field(i,j,k) = field(i+1,j,k)
    END DO
    END DO
    END DO

    CALL field_bc(field, ng)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i, iy, iz, isuby, isubz, errcode
    INTEGER(i8) :: ipart, npart_per_cell, n_frac
    REAL(num) :: cell_frac_y, cy2
    REAL(num) :: cell_frac_z, cz2
    REAL(num), DIMENSION(-1:1) :: gy, gz
    REAL(num), DIMENSION(c_ndirs) :: temp_local, drift_local
    REAL(num) :: npart_frac
    REAL(num) :: weight_local, x0, dmin, dmax, wdata
    TYPE(parameter_pack) :: parameters
    TYPE(particle_species), POINTER :: species
    REAL(num), DIMENSION(c_ndirs, 2) :: ranges

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    IF (nproc > 1) THEN
      IF (SIZE(density,1) /= ny+2 .OR. SIZE(density,2) /= nz+2) THEN
        DEALLOCATE(density, temperature, drift)
        ALLOCATE(density(0:ny+1, 0:nz+1))
        ALLOCATE(temperature(0:ny+1, 0:nz+1, 1:3))
        ALLOCATE(drift(0:ny+1, 0:nz+1, 1:3))
      END IF
    END IF

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
      DO i = 1, 3
        DO iz = 0, nz+1
          parameters%pack_iz = iz
          DO iy = 0, ny+1
            parameters%pack_iy = iy
            temperature(iy,iz,i) = evaluate_with_parameters( &
                species%temperature_function(i), parameters, errcode)
            drift(iy,iz,i) = evaluate_with_parameters( &
                species%drift_function(i), parameters, errcode)
          END DO
        END DO
      END DO
      DO iz = 0, nz+1
        parameters%pack_iz = iz
        DO iy = 0, ny+1
          parameters%pack_iy = iy
          density(iy,iz) = evaluate_with_parameters(species%density_function, &
              parameters, errcode)
          IF (density(iy,iz) > dmax) density(iy,iz) = dmax
          IF (density(iy,iz) < dmin) density(iy,iz) = 0.0_num
        END DO
      END DO

      x0 = x_grid_max + 0.5_num * dx
      DO iz = 1, nz
        DO iy = 1, ny
          IF (density(iy,iz) < dmin) CYCLE

          ! Place extra particle based on probability
          n_frac = 0
          IF (npart_frac > 0.0_num) THEN
            IF (random() < npart_frac) n_frac = 1
          END IF

          wdata = dx * dy * dz / (npart_per_cell + n_frac)

          DO ipart = 1, npart_per_cell + n_frac
            CALL create_particle(current)
#ifdef PER_PARTICLE_CHARGE_MASS
            current%mass = species%mass
            current%charge = species%charge
#endif
            cell_frac_y = 0.5_num - random()
            cell_frac_z = 0.5_num - random()
            current%part_pos(1) = x0 + random() * dx
            current%part_pos(2) = y(iy) - cell_frac_y * dy
            current%part_pos(3) = z(iz) - cell_frac_z * dz

            ! Always use the triangle particle weighting for simplicity
            cy2 = cell_frac_y**2
            gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
            gy( 0) = 0.75_num - cy2
            gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

            cz2 = cell_frac_z**2
            gz(-1) = 0.5_num * (0.25_num + cz2 + cell_frac_z)
            gz( 0) = 0.75_num - cz2
            gz( 1) = 0.5_num * (0.25_num + cz2 - cell_frac_z)

            temp_local = 0.0_num
            drift_local = 0.0_num
            DO i = 1, c_ndirs
              DO isubz = -1, 1
                DO isuby = -1, 1
                  temp_local(i) = temp_local(i) + gy(isuby) * gz(isubz) &
                      * temperature(iy+isuby, iz+isubz, i)
                  drift_local(i) = drift_local(i) + gy(isuby) * gz(isubz) &
                      * drift(iy+isuby, iz+isubz, i)
                END DO
              END DO
            END DO

            IF (species%ic_df_type == c_ic_df_thermal) THEN
              DO i = 1, c_ndirs
                current%part_p(i) = momentum_from_temperature(species%mass, &
                    temp_local(i), drift_local(i))
              END DO
            ELSE IF (species%ic_df_type == c_ic_df_relativistic_thermal) THEN
              current%part_p = momentum_from_temperature_relativistic(&
                  species%mass, temp_local, species%fractional_tail_cutoff, &
                  drift_local)
            ELSE IF (species%ic_df_type == c_ic_df_arbitrary) THEN
              parameters%use_grid_position = .FALSE.
              parameters%pack_pos = current%part_pos
              parameters%pack_iy = iy
              parameters%pack_iz = iz
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

            weight_local = 0.0_num
            DO isubz = -1, 1
              DO isuby = -1, 1
                weight_local = weight_local + gy(isuby) * gz(isubz) &
                    * density(iy+isuby, iz+isubz)
              END DO
            END DO

            current%weight = weight_local * wdata
#ifdef PARTICLE_DEBUG
            current%processor = rank
            current%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(append_list, current)
          END DO
        END DO
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
          IF (current%part_pos(1) < x_min) THEN
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
        bc_field(c_bd_y_min) = bc_y_min_after_move
        bc_field(c_bd_y_max) = bc_y_max_after_move
        bc_field(c_bd_z_min) = bc_z_min_after_move
        bc_field(c_bd_z_max) = bc_z_max_after_move
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
