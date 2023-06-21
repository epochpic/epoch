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

MODULE injectors

  USE partlist
  USE particle_temperature
  USE evaluator
  USE random_generator
  USE utilities
  USE file_injectors

  IMPLICIT NONE

  REAL(num) :: flow_limit_val = 10.0_num

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = -1.0_num
    injector%species = -1
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
    injector%has_t_end = .FALSE.
    injector%density_min = 0.0_num
    injector%density_max = HUGE(1.0_num)
    injector%use_flux_injector = .TRUE.
    injector_boundary(boundary) = .TRUE.
    NULLIFY(injector%next)

    injector%depth = 1.0_num
    need_random_state = .TRUE.

    ! Additional variables for file injectors
    injector%inject_from_file = .FALSE.
    injector%file_finished = .FALSE.
    injector%px_data_given = .FALSE.
    injector%py_data_given = .FALSE.
    injector%pz_data_given = .FALSE.
    injector%t_data_given = .FALSE.
#ifndef PER_SPECIES_WEIGHT
    injector%w_data_given = .FALSE.
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    injector%id_data_given = .FALSE.
#endif

  END SUBROUTINE init_injector



  SUBROUTINE attach_injector(injector)

    TYPE(injector_block), POINTER :: injector
    TYPE(injector_block), POINTER :: current

    NULLIFY(injector%next)

    IF (ASSOCIATED(injector_list)) THEN
      current => injector_list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => injector
    ELSE
      injector_list => injector
    END IF

  END SUBROUTINE attach_injector



  SUBROUTINE deallocate_injectors

    TYPE(injector_block), POINTER :: current, next
    INTEGER :: i

    current => injector_list
    DO WHILE(ASSOCIATED(current))
      next => current%next
      IF (current%density_function%init) &
          CALL deallocate_stack(current%density_function)
      DO i = 1, 3
        IF (current%temperature_function(i)%init) &
            CALL deallocate_stack(current%temperature_function(i))
        IF (current%drift_function(i)%init) &
            CALL deallocate_stack(current%drift_function(i))
      END DO
      DEALLOCATE(current)
      current => next
    END DO

  END SUBROUTINE deallocate_injectors



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    current => injector_list
    DO WHILE(ASSOCIATED(current))
      IF (.NOT. current%inject_from_file) THEN
        CALL run_single_injector(current)
      ELSE
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector)

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: bdy_pos, cell_size
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, density, vol, p_drift, p_ratio
    REAL(num) :: npart_ideal, itemp, v_inject_s, density_correction, dir_mult
    REAL(num) :: v_inject_dt
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: weight_fac
#endif
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, idir, dir_index, flux_dir, flux_dir_cell
    INTEGER :: direction
    TYPE(parameter_pack) :: parameters
    REAL(num), PARAMETER :: sqrt2 = SQRT(2.0_num)
    REAL(num), PARAMETER :: sqrt2_inv = 1.0_num / sqrt2
    REAL(num), PARAMETER :: sqrt2pi_inv = 1.0_num / SQRT(2.0_num * pi)

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. injector%has_t_end) &
        RETURN

    direction = injector%boundary

    IF (.NOT. is_boundary(direction)) RETURN

    IF (direction == c_bd_x_min) THEN
      bdy_pos = x_min
      parameters%pack_ix = 0
      dir_mult = 1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
    ELSE IF (direction == c_bd_x_max) THEN
      bdy_pos = x_max
      parameters%pack_ix = nx
      dir_mult = -1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
    ELSE
      RETURN
    END IF

    parameters%use_grid_position = .TRUE.

    CALL populate_injector_properties(injector, parameters, density=density)

    IF (density < injector%density_min) RETURN

    CALL populate_injector_properties(injector, parameters, &
        temperature=temperature, drift=drift)

    IF (injector%use_flux_injector) THEN
      flux_dir = dir_index
    ELSE
      flux_dir = -1
    END IF

    vol = dx
    bdy_pos = bdy_pos - 0.5_num * dir_mult * cell_size * png

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
#ifndef PER_SPECIES_WEIGHT
    weight_fac = vol / injector%npart_per_cell
#endif

    ! Assume agressive maximum thermal momentum, all components
    ! like hottest component
    p_therm = SQRT(mass * kb * MAXVAL(temperature))
    p_inject_drift = drift(dir_index)
    flux_dir_cell = flux_dir

    IF (flux_dir_cell /= -1) THEN
      ! Drift adjusted so that +ve is 'inwards' through boundary
      p_drift = p_inject_drift * dir_mult

      ! Average momentum of inflowing part
      ! For large inwards drift, is asymptotic to drift
      ! Otherwise it is a complicated expression
      ! Inwards drift - lhs terms are same sign -> +ve
      IF (p_drift > flow_limit_val * p_therm) THEN
        ! For sufficiently large drifts, net inflow -> p_drift
        gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
        v_inject_s = p_inject_drift / gamma_mass
        density_correction = 1.0_num
        ! Large drift flux Maxwellian can be approximated by a
        ! non-flux Maxwellian
        flux_dir_cell = -1
      ELSE IF (p_drift < -flow_limit_val * p_therm) THEN
        ! Net is outflow - inflow velocity is zero. No particles injected
        RETURN
      ELSE IF (ABS(p_therm) < c_tiny) THEN
        RETURN
      ELSE IF (ABS(p_drift) < p_therm * 1.0e-9_num) THEN
        v_inject_s = 2.0_num * sqrt2pi_inv * p_therm &
            + (1.0_num - 2.0_num * sqrt2 / pi) * p_drift
        gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
        v_inject_s = v_inject_s / gamma_mass
        density_correction = 0.5_num
      ELSE
        p_ratio = sqrt2_inv * p_drift / p_therm

        ! Fraction of the drifting Maxwellian distribution inflowing
        density_correction = 0.5_num * (1.0_num + erf_func(p_ratio))
        IF (density_correction < c_tiny) RETURN

        ! Below is actually MOMENTUM, will correct on next line
        v_inject_s = dir_mult * (p_drift &
            + sqrt2pi_inv * p_therm * EXP(-p_ratio**2) / density_correction)

        gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
        v_inject_s = v_inject_s / gamma_mass
      END IF
    ELSE
      ! User asked for Maxwellian only - no correction to apply
      gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
      v_inject_s = p_inject_drift / gamma_mass
      density_correction = 1.0_num
    END IF

    v_inject = ABS(v_inject_s)
    v_inject_dt = dt * v_inject_s

    npart_ideal = injector%npart_per_cell * v_inject * density_correction &
        * dt / cell_size
    itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
        * (1.0_num - npart_ideal / injector%npart_per_cell))) + npart_ideal
    injector%depth = injector%depth - itemp

    IF (injector%depth >= 0.0_num) RETURN

    parts_this_time = FLOOR(ABS(injector%depth - 1.0_num))
    injector%depth = injector%depth + REAL(parts_this_time, num)

    IF (parts_this_time < 1) RETURN

    CALL create_empty_partlist(plist)

    DO ipart = 1, parts_this_time
      CALL create_particle(new)

      new%part_pos = bdy_pos - random() * v_inject_dt
      parameters%pack_pos = new%part_pos
      parameters%use_grid_position = .FALSE.

#ifdef PER_SPECIES_WEIGHT
      CALL populate_injector_properties(injector, parameters, &
          temperature=temperature, drift=drift)
#else
      CALL populate_injector_properties(injector, parameters, density, &
          temperature, drift)
#endif

      DO idir = 1, 3
        IF (idir == flux_dir_cell) THEN
          ! Drift is signed - dir mult is the direciton we want to get
          new%part_p(idir) = flux_momentum_from_temperature(&
              mass, temperature(idir), drift(idir), dir_mult)
        ELSE
          new%part_p(idir) = momentum_from_temperature(mass, &
              temperature(idir), drift(idir))
        END IF
      END DO
#ifdef PER_PARTICLE_CHARGE_MASS
      new%charge = species_list(injector%species)%charge
      new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
      density = MIN(density, injector%density_max)
      new%weight = weight_fac * density
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      ! For photons, assign additional variable used in photon particle-push
      IF (species_list(injector%species)%species_type == c_species_id_photon) &
          THEN 
        new%particle_energy = SQRT(SUM(new%part_p**2)) * c
      END IF
#endif
      CALL add_particle_to_partlist(plist, new)
    END DO

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

  END SUBROUTINE run_single_injector



  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT), OPTIONAL :: density
    REAL(num), DIMENSION(3), INTENT(OUT), OPTIONAL :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    IF (PRESENT(density)) THEN
      density = 0.0_num
      IF (injector%density_function%init) THEN
        density = MAX(evaluate_with_parameters(injector%density_function, &
            parameters, errcode), 0.0_num)
      END IF
    END IF

    ! Stack can only be time varying if valid. Change if this isn't true
    IF (PRESENT(temperature)) THEN
      temperature(:) = 0.0_num
      DO i = 1, 3
        IF (injector%temperature_function(i)%init) THEN
          temperature(i) = &
              MAX(evaluate_with_parameters(injector%temperature_function(i), &
                  parameters, errcode), 0.0_num)
        END IF
      END DO
    END IF

    IF (PRESENT(drift)) THEN
      drift(:) = 0.0_num
      DO i = 1, 3
        IF (injector%drift_function(i)%init) THEN
          drift(i) = &
              evaluate_with_parameters(injector%drift_function(i), &
                                       parameters, errcode)
        END IF
      END DO
    END IF

    IF (errcode /= c_err_none) CALL abort_code(errcode)

  END SUBROUTINE populate_injector_properties



  SUBROUTINE finish_injector_setup

    TYPE(injector_block), POINTER :: current

    current => injector_list
    DO WHILE(ASSOCIATED(current))
      IF (.NOT. current%inject_from_file) THEN
        CALL finish_single_injector_setup(current)
      END IF
      current => current%next
    END DO

  END SUBROUTINE finish_injector_setup



  SUBROUTINE finish_single_injector_setup(injector)

    TYPE(injector_block), POINTER :: injector
    TYPE(particle_species), POINTER :: species
    INTEGER :: i

    species => species_list(injector%species)
    IF (injector%npart_per_cell < 0.0_num) THEN
      injector%npart_per_cell = species%npart_per_cell
    END IF

    IF (.NOT.injector%density_function%init) THEN
      CALL copy_stack(species%density_function, injector%density_function)
    END IF

    DO i = 1, 3
      IF (.NOT.injector%drift_function(i)%init) THEN
        CALL copy_stack(species%drift_function(i), injector%drift_function(i))
      END IF
      IF (.NOT.injector%temperature_function(i)%init) THEN
        CALL copy_stack(species%temperature_function(i), &
            injector%temperature_function(i))
      END IF
    END DO

  END SUBROUTINE finish_single_injector_setup



  SUBROUTINE create_boundary_injector(ispecies, bnd)

    INTEGER, INTENT(IN) :: ispecies, bnd
    TYPE(injector_block), POINTER :: working_injector

    species_list(ispecies)%bc_particle(bnd) = c_bc_open
    use_injectors = .TRUE.

    ALLOCATE(working_injector)

    CALL init_injector(bnd, working_injector)
    working_injector%species = ispecies

    CALL attach_injector(working_injector)

  END SUBROUTINE create_boundary_injector



  SUBROUTINE setup_injector_depths(boundary, depths, injector_count)

    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(:), INTENT(IN) :: depths
    INTEGER, INTENT(OUT) :: injector_count
    TYPE(injector_block), POINTER :: injector
    INTEGER :: inj, bnd

    inj = 1
    injector => injector_list

    DO WHILE(ASSOCIATED(injector))
      bnd = injector%boundary
      IF (bnd == boundary) THEN
        ! Exclude ghost cells
        injector%depth = depths(inj)
        inj = inj + 1
      END IF
      injector => injector%next
    END DO

    injector_count = inj - 1

  END SUBROUTINE setup_injector_depths

END MODULE injectors
