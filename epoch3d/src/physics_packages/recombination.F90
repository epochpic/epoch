! Copyright (C) 2009-2021 University of Warwick
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
!
!-------------------------------------------------------------------------------
!
! This module contains subroutines used to calculate recombination, provided 
! recombination tables have been provided for a given atomic number and charge
! state.
!
! Here we consider dielectronic, radiative and three-body recombination. The 
! electron temperature in a cell determines the recombination rates, which can 
! be used to obtain a recombination cross-section. Each electron in the cell has
! a chance of triggering recombination. Radiative rates are taken from FLYCHK,
! as are most dielectronic rates (some include contributions from JAC). The
! 3-body recombination rates have been derived following the theory of Kremp's
! "Quantum statistics of non-ideal plasmas", using detailed balance arguments.

MODULE recombination

  USE calc_df
  USE collisions
  USE collision_ionise
  
  IMPLICIT NONE

  TYPE(interpolation_state), SAVE :: last_table_state
  REAL(num) :: alpha_j, el_ne

CONTAINS

  SUBROUTINE setup_recombination_tables

    ! This subroutine extracts di-electric recombination rates from files as a
    ! function of electron temperature, and saves them to the species in 1D 
    ! arrays. Same is done for radiative recombination.
    !
    ! Three-body recombination requires collisional ionisation cross-sections,
    ! so calculate these if collisional ionisation is switched off

    INTEGER :: ispecies, io, iu

    ! Loop through all species which are set to recombine
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%recombine) THEN
        
        use_recombination = .TRUE.
        use_particle_lists = .TRUE.
        species_list(ispecies)%make_secondary_list = .TRUE.
        species_list(species_list(ispecies)%recombine_to_species) &
            %make_secondary_list = .TRUE.

        IF (use_dielectronic_recombination) THEN
          CALL setup_dr_tables(ispecies)
        END IF

        IF (use_radiative_recombination) THEN
          CALL setup_rr_tables(ispecies)
        END IF

      END IF
    END DO

    ! If recombination is present, we must create secondary lists for e-
    IF (use_recombination) THEN
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%electron) &
            species_list(ispecies)%make_secondary_list = .TRUE.
      END DO
    END IF

    ! If three-body recombination is present, but not collisional ionisation, 
    ! load the collisional ionisation cross section scripts anyway
    IF (use_recombination .AND. use_three_body_recombination &
        .AND. .NOT. use_collisional_ionisation) CALL setup_coll_ionise_tables 

  END SUBROUTINE setup_recombination_tables



  SUBROUTINE setup_dr_tables(ispecies)

    ! Read the dielectronic recombination data tables from FLYCHK, converted to
    ! SI units, and populate the recombine_temp_dr and recombine_rate_dr arrays
    ! for the species with ID ispecies

    INTEGER, INTENT(in) :: ispecies
    INTEGER :: entries, z_int, q_int, q_read, i_read
    LOGICAL :: exists
    CHARACTER(LEN=9) :: file_name

    ! Data is in src/physics_packages/TABLES/dr_rates/
    ! Filename format dr_**.dat where "*" is replaced by the atomic number
    z_int = species_list(ispecies)%atomic_no
    WRITE(file_name, '(A3,I2.2,A4)') 'dr_', z_int, '.dat'

    ! Check if the data-file exists for the current ion species
    INQUIRE(FILE=TRIM(physics_table_location)  // "/dr_rates/"&
        // file_name, EXIST=exists)

    ! No data-file, so switch off recombination for this species, and skip
    ! to the next species
    IF (.NOT. exists) THEN
      IF (rank == 0) THEN
        PRINT*, 'No data present for dielectronic recombination from ', &
            TRIM(species_list(ispecies)%name)
      END IF
      species_list(ispecies)%recombine = .FALSE.
      RETURN
    END IF

    ! Attempt to read data file. Typical format (reduced dr_03.dat):
    ! ! n   Temperature [K]
    ! !  36 5.802259e+03 1.160452e+04 1.740678e+04 ...
    ! ! Q   Rate [m3/s]
    ! !   1 7.596000e-59 5.057000e-38 3.589000e-31 ...
    ! !   2 1.000000e-76 2.412000e-49 7.918000e-39 ...
    ! !   3 1.000000e-76 1.000000e-76 1.000000e-76 ...
    ! 
    ! In this example: 
    ! ! Top-left number 36 gives the number of temperature samples
    ! ! Numbers after 36 give temperature sample points in Kelvin
    ! ! Then we have the rates: left column for initial charge of transition
    ! ! Numbers following Q are the rates at each temperature sample point
    OPEN(UNIT = lu, &
        FILE = TRIM(physics_table_location)  // "/dr_rates/"&
        // file_name, &
        STATUS = 'OLD')
    
    ! Read the number of temperature entries, ignoring first header line
    READ(lu,*)
    READ(lu,*) entries
    species_list(ispecies)%recombine_array_size_dr = entries
    CLOSE(lu)

    ALLOCATE(species_list(ispecies)%recombine_temp_dr(entries))
    ALLOCATE(species_list(ispecies)%recombine_rate_dr(entries))

    ! Read the full first row of data
    OPEN(UNIT = lu, &
      FILE = TRIM(physics_table_location)  // "/dr_rates/"&
      // file_name, &
      STATUS = 'OLD')
    READ(lu,*)
    READ(lu,*) entries, species_list(ispecies)%recombine_temp_dr(:)

    ! Skip to the line for the current species charge-state
    q_int = NINT(species_list(ispecies)%charge / q0)
    DO i_read = 1,q_int
      READ(lu,*)
    END DO

    ! Read recombination rate [m³/s]
    READ(lu,*) q_read, species_list(ispecies)%recombine_rate_dr(:)

    CLOSE(lu)

  END SUBROUTINE setup_dr_tables



  SUBROUTINE setup_rr_tables(ispecies)

    ! Read the radiative recombination data tables from FLYCHK, converted to
    ! SI units, and populate the recombine_temp_rr and recombine_rate_rr arrays
    ! for the species with ID ispecies

    INTEGER, INTENT(in) :: ispecies
    INTEGER :: entries, z_int, q_int, q_read, i_read
    LOGICAL :: exists
    CHARACTER(LEN=9) :: file_name

    ! Data is in src/physics_packages/TABLES/dr_rates/
    ! Filename format dr_**.dat where "*" is replaced by the atomic number
    z_int = species_list(ispecies)%atomic_no
    WRITE(file_name, '(A3,I2.2,A4)') 'rr_', z_int, '.dat'

    ! Check if the data-file exists for the current ion species
    INQUIRE(FILE=TRIM(physics_table_location)  // "/rr_rates/"&
        // file_name, EXIST=exists)

    ! No data-file, so switch off recombination for this species, and skip
    ! to the next species
    IF (.NOT. exists) THEN
      IF (rank == 0) THEN
        PRINT*, 'No data present for radiative recombination from ', &
            TRIM(species_list(ispecies)%name)
      END IF
      species_list(ispecies)%recombine = .FALSE.
      RETURN
    END IF

    ! Attempt to read data file. Typical format (reduced dr_03.dat):
    ! ! n   Temperature [K]
    ! !  36 5.802259e+03 1.160452e+04 1.740678e+04 ... 
    ! ! Q   Rate [m3/s]
    ! !   1 3.984000e-19 2.381000e-19 1.724000e-19 ...
    ! !   2 3.318000e-18 2.116000e-18 1.610000e-18 ...
    ! !   3 8.023000e-18 5.281000e-18 4.099000e-18 ...
    ! 
    ! In this example: 
    ! ! Top-left number 36 gives the number of temperature samples
    ! ! Numbers after 36 give temperature sample points in Kelvin
    ! ! Then we have the rates: left column for initial charge of transition
    ! ! Numbers following Q are the rates at each temperature sample point
    OPEN(UNIT = lu, &
        FILE = TRIM(physics_table_location)  // "/rr_rates/"&
        // file_name, &
        STATUS = 'OLD')
    
    ! Read the number of temperature entries, ignoring first header line
    READ(lu,*)
    READ(lu,*) entries
    species_list(ispecies)%recombine_array_size_rr = entries
    CLOSE(lu)

    ALLOCATE(species_list(ispecies)%recombine_temp_rr(entries))
    ALLOCATE(species_list(ispecies)%recombine_rate_rr(entries))

    ! Read the full first row of data
    OPEN(UNIT = lu, &
      FILE = TRIM(physics_table_location)  // "/rr_rates/"&
      // file_name, &
      STATUS = 'OLD')
    READ(lu,*)
    READ(lu,*) entries, species_list(ispecies)%recombine_temp_rr(:)

    ! Skip to the line for the current species charge-state
    q_int = NINT(species_list(ispecies)%charge / q0)
    DO i_read = 1,q_int
      READ(lu,*)
    END DO

    ! Read recombination rate [m³/s]
    READ(lu,*) q_read, species_list(ispecies)%recombine_rate_rr(:)

    CLOSE(lu)

  END SUBROUTINE setup_rr_tables


  SUBROUTINE run_recombination

    ! This subroutine is called by the main PIC loop, and triggers the
    ! recombination calculation. The species_list is cycled through to identify
    ! pairs of valid electron and ion species which can recombine, and a
    ! subroutine is called to sample the recombination. Ions are added to the 
    ! correct species list, and electron macro-particles are removed
    !
    ! This subroutine can perform dielectronic, radiatiave, and 3-body 
    ! recombination

    INTEGER :: ispecies, jspecies
    INTEGER :: electron_species, ion_species, recombined_species
    INTEGER(i8) :: ix, iy, iz
    TYPE(particle_list), POINTER :: p_list
    TYPE(particle_list) :: list_i_recombined
    LOGICAL :: i_is_ion, i_is_electron, calc_recombination

    ! Only run this type of recombination for dielectronic or radiative
    ! recombination
    IF (.NOT. use_dielectronic_recombination .AND. .NOT. &
      use_radiative_recombination) RETURN

    CALL create_empty_partlist(list_i_recombined)

    ! Shuffle ion species lists ahead of calculation
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%ionise .AND. .NOT. &
          species_list(ispecies)%is_shuffled) THEN
        DO iz = 1, nz
        DO iy = 1, ny
        DO ix = 1, nx
          p_list => species_list(ispecies)%secondary_list(ix,iy,iz)
          CALL shuffle_particle_list_random(p_list)
        END DO ! ix
        END DO ! iy
        END DO ! iz
        species_list(ispecies)%is_shuffled = .TRUE.
      END IF
    END DO

    ! Loop over all species and extract the electron-ion pairs
    DO ispecies = 1, n_species

      ! Photons cannot recombine with ions
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE

      ! Variables to track the identity of species ispecies
      i_is_ion = .FALSE.
      i_is_electron = .FALSE.
      IF (species_list(ispecies)%electron) THEN
        ! Current species is an electron species
        electron_species = ispecies
        i_is_electron = .TRUE.
      ELSE IF (species_list(ispecies)%recombine) THEN
        ! Current species is an ion species which can recombine
        ion_species = ispecies
        i_is_ion = .TRUE.
      ELSE
        ! Skip species which are not involved in recombination
        CYCLE
      END IF

      ! Loop over remaining species for recombination partners
      DO jspecies = ispecies, n_species

        ! Determine whether ispecies and jspecies are recombination partners
        calc_recombination = .FALSE.
        IF (species_list(jspecies)%electron .AND. i_is_ion) THEN
          ! jspecies is an electron species which can recombine with ispecies
          electron_species = jspecies
          calc_recombination = .TRUE.
        ELSE IF (species_list(jspecies)%recombine .AND. i_is_electron) THEN
          ! jspecies can recombine with ispecies
          ion_species = jspecies
          calc_recombination = .TRUE.
        ELSE
          ! Skip species which don't form a recombination pair with ispecies
          CYCLE
        END IF

        ! Get corresponding species ID for new list
        recombined_species = species_list(ion_species)%recombine_to_species

        ! Perform recombination for valid electron-ion pairs
        DO iz = 1, nz
        DO iy = 1, ny
        DO ix = 1, nx
        IF (calc_recombination) THEN

            ! Estimate collisional-ionisation rate (for 3-body recombination)
            IF (use_three_body_recombination) THEN
              CALL calculate_alpha_and_ne(ion_species, ix, iy, iz)
            END IF

            ! Apply recombination, moving recombined ions to the new particle 
            ! lists
            CALL recombine_ei(ion_species, electron_species, &
                species_list(electron_species)%secondary_list(ix,iy,iz), &
                species_list(ion_species)%secondary_list(ix,iy,iz), &
                list_i_recombined)

            ! Shift recombined ions to the correct species
            CALL append_partlist( &
                species_list(recombined_species)%secondary_list(ix,iy,iz), &
                list_i_recombined)
          END IF
        END DO ! ix
        END DO ! iy
        END DO ! iz

      END DO ! jspecies
    END DO ! ispecies

  END SUBROUTINE run_recombination



  SUBROUTINE recombine_ei(ion_species, el_species, list_e_in, list_i_in, &
      list_i_recombined)

    INTEGER, INTENT(IN) :: ion_species, el_species
    TYPE(particle_list), INTENT(INOUT) :: list_e_in, list_i_in
    TYPE(particle_list), INTENT(INOUT) :: list_i_recombined
    TYPE(particle), POINTER :: electron, ion
    TYPE(particle), POINTER :: next_ion, next_electron
    INTEGER(i8) :: e_count, ion_count, i_ion, i_el, recombined_count
    LOGICAL, ALLOCATABLE :: e_recombined(:), i_recombined(:)
    LOGICAL :: ion_at_rest
    LOGICAL :: create_recombination, e_collide_again
    REAL(num) :: sum_wi, n_i, inv_cell_volume, ion_mass, el_weight, ion_weight
    REAL(num) :: sum_we, el_temp, el_p_i(3)
    REAL(num) :: el_p2, mean_el_p2, mean_el_pxe
    REAL(num) :: ion_p2, gamma_i, beta_i
    REAL(num) :: recombine_frac, recombine_rate, uncombined_frac, recombine_no
    REAL(num) :: recombine_prob, recombined_ion_e_i

    ! Ignore collisions from empty species
    e_count = list_e_in%count
    ion_count = list_i_in%count
    IF (e_count == 0 .OR. ion_count == 0) RETURN 

    ! Allocate arrays to track which electrons and ions are involved in
    ! recombination events
    ALLOCATE(e_recombined(e_count), i_recombined(ion_count))
    e_recombined = .FALSE.
    i_recombined = .FALSE.

    ! Pre-read weights to support this compiler option
#ifdef PER_SPECIES_WEIGHT
    el_weight = species_list(el_species)%weight
    ion_weight = species_list(ion_species)%weight
#endif

    ! Calculate initial ion number density
    ion => list_i_in%head
    sum_wi = 0
    DO i_ion = 1, ion_count
#ifndef PER_SPECIES_WEIGHT
      ion_weight = ion%weight
#endif
      sum_wi = sum_wi + ion_weight
      ion => ion%next
    END DO
    inv_cell_volume = 1.0_num / (dx * dy * dz)
    n_i = sum_wi * inv_cell_volume

    ! Calculate electron <p^2> and <px E>, which are needed to calculate
    ! temperature in the ion-rest-frame
    electron => list_e_in%head
    sum_we = 0
    mean_el_p2 = 0.0_num
    mean_el_pxe = 0.0_num
    DO i_el = 1, e_count 
#ifndef PER_SPECIES_WEIGHT
      el_weight = electron%weight
#endif
      sum_we = sum_we + electron%weight
      el_p2 = SUM(electron%part_p(:)**2)
      mean_el_p2 = mean_el_p2 + el_p2 * el_weight
      mean_el_pxe = mean_el_pxe + electron%part_p(1) &
          * SQRT(el_p2*c**2 + m0c2**2) * el_weight
      electron => electron%next
    END DO
    mean_el_p2 = mean_el_p2 / sum_we
    mean_el_pxe = mean_el_pxe / sum_we

    ! Each macro-e collision must pair to a macro-ion. To ensure there is always
    ! a macro-ion available, temporarily make the ion list circular
    list_i_in%tail%next => list_i_in%head
    ion => list_i_in%head
    i_ion = 1

    ! Loop over all macro-electrons
    recombined_count = 0
    electron => list_e_in%head
    DO i_el = 1, e_count

      ! Some macro-electrons must collide with many macro-ions to create the
      ! right number of recombined ions, this loop allows for multiple 
      ! collisions
      uncombined_frac = 1.0_num
      e_collide_again = .TRUE.
      DO WHILE (e_collide_again)

#ifndef PER_SPECIES_WEIGHT
        el_weight = electron%weight
        ion_weight = ion%weight
#endif

        ! Check if ion is moving or at rest
        ion_p2 = DOT_PRODUCT(ion%part_p, ion%part_p)
        ion_at_rest = ion_p2 < 1.0e-100_num

        ! If the ion is moving, then perform a Lorentz transform to the ion
        ! rest-frame (rates are evaluted in rest-frame)
        IF (.NOT. ion_at_rest) THEN
#ifdef PER_PARTICLE_CHARGE_MASS
          ion_mass = ion%mass
#else
          ion_mass = species_list(ion_species)%mass
#endif
          gamma_i = SQRT(ion_p2 / (ion_mass * c)**2 + 1.0_num)
          beta_i = SQRT(1.0_num - 1.0_num / gamma_i**2)

          ! Calculate temperature in ion-rest-frame, assuming isotropic e-
          el_temp = (mean_el_p2 &
              * ((gamma_i**2 + 2.0_num) / 3.0_num + (gamma_i * beta_i)**2) &
              + (gamma_i * beta_i * m0 * c)**2 &
              - 2.0_num * beta_i * gamma_i**2 * mean_el_pxe / c) &
              / (3 * m0 * kb)
        ELSE
          ! Ion is already at rest
          el_temp = mean_el_p2 / (3.0_num * m0 * kb)
        END IF

        recombine_rate = 0.0_num
        IF (use_dielectronic_recombination) THEN
          recombine_rate = recombine_rate + &
              find_value_from_table_1d_recombine(el_temp, &
              species_list(ion_species)%recombine_array_size_dr, &
              species_list(ion_species)%recombine_temp_dr, &
              species_list(ion_species)%recombine_rate_dr, last_table_state)
        END IF
        IF (use_radiative_recombination) THEN
          recombine_rate = recombine_rate + &
              find_value_from_table_1d_recombine(el_temp, &
              species_list(ion_species)%recombine_array_size_rr, &
              species_list(ion_species)%recombine_temp_rr, &
              species_list(ion_species)%recombine_rate_rr, last_table_state)
        END IF
        IF (use_three_body_recombination) THEN
          recombine_rate = recombine_rate + rate_3br(ion_species, el_temp, &
              el_ne)
        END IF

        ! Expected number of recombined ions from the incident electron
        recombine_prob = 1.0_num - EXP(-recombine_rate * n_i * dt)
        recombine_no = el_weight * recombine_prob * uncombined_frac

        ! Do we recombine a macro-ion?
        IF (recombine_no > ion_weight) THEN
          ! Weight of macro-ion less than expected number of recombinations
          ! Recombine macro-ion but find a new target for macro-electron
          create_recombination = .TRUE.

          ! Save fraction of recombinations left to make
          recombine_frac = ion_weight / recombine_no
          uncombined_frac = uncombined_frac * (1.0_num - recombine_frac)
        ELSE
          ! Weight of macro-ion >= expected number of recombinations
          ! Sample the probability of recombination of the macro-ion
          create_recombination = (random() <= recombine_no / ion_weight)

          ! Electron has created all the recombinations it will (if any), do not 
          ! re-collide
          e_collide_again = .FALSE.
        END IF

        ! Mark particles as recombined
        IF (create_recombination) THEN
          i_recombined(i_ion) = .TRUE.
          e_recombined(i_el) = .TRUE.
          recombined_count = recombined_count + 1
          n_i = n_i - ion_weight * inv_cell_volume

          ! Conserve momentum of recombined particle
          ion%part_p = ion%part_p + electron%part_p
        END IF

        ! Check if any ions remain for recombination by the current electron
        IF (recombined_count == ion_count) EXIT

        ! Ensure the next ion is a valid target (not previously recombined)
        ion => ion%next
        i_ion = MOD(i_ion, ion_count) + 1
        DO WHILE (i_recombined(i_ion))
          ion => ion%next
          i_ion = MOD(i_ion, ion_count) + 1
        END DO
      END DO

      ! Check if any ions remain for ionisation of the next electrons
      IF (recombined_count == ion_count) EXIT

      electron => electron%next
    END DO

    ! Restore the tail of the ion list
    NULLIFY(list_i_in%tail%next)

    ! Move recombined macro-ions to the recombined-ion particle list
    ion => list_i_in%head
    DO i_ion = 1, ion_count
      next_ion => ion%next
      IF (i_recombined(i_ion)) THEN
        CALL remove_particle_from_partlist(list_i_in, ion)
        CALL add_particle_to_partlist(list_i_recombined, ion)
      END IF
      ion => next_ion
    END DO

    ! Remove recombined macro-electrons from the simulation
    electron => list_e_in%head
    DO i_el = 1, e_count
      next_electron => electron%next
      IF (e_recombined(i_el)) THEN
        CALL remove_particle_from_partlist(list_e_in, electron)
      END IF
      electron => next_electron
    END DO

    DEALLOCATE(e_recombined, i_recombined)

  END SUBROUTINE recombine_ei



  SUBROUTINE calculate_alpha_and_ne(ion_species, ix, iy, iz)

    ! Calculate alpha for the current particle species. This describes the rate 
    ! of collisional ionisation for bound electrons in the outermost shell of
    ! the recombined ion. This implies the recombined electron would go to the 
    ! shell which, under the ground-state transition approximation, would lose
    ! the electron upon ionisation. This is determined through comparison of the
    ! ground-state configurations of the initial and recombined ions.
    !
    ! Subroutine also calculates the total electron number density in the cell

    INTEGER, INTENT(IN) :: ion_species
    INTEGER(i8), INTENT(IN) ::  ix, iy, iz
    INTEGER :: recombined_species, ispecies
    TYPE(particle), POINTER :: current
    REAL(num) :: e_weight_sum, sigma_v_sum, weight, sigma_ci_outer
    REAL(num) :: el_p2, el_e, el_v, el_ke

    recombined_species = species_list(ion_species)%recombine_to_species

    ! Loop over all electrons from all electron species
    e_weight_sum = 0.0_num
    sigma_v_sum = 0.0_num
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%electron) THEN
        current => species_list(ispecies)%secondary_list(ix,iy,iz)%head
        DO WHILE(ASSOCIATED(current))
#ifdef PER_SPECIES_WEIGHT
          weight = species_list(ispecies)%weight
#else
          weight = current%weight
#endif
          e_weight_sum = e_weight_sum + weight

          ! Electron speed
          el_p2 = current%part_p(1)**2 + current%part_p(2)**2 &
              + current%part_p(3)**2
          el_e = SQRT(el_p2*c**2 + m0**2*c**4)
          el_v = SQRT(el_p2) * c**2 / el_e

          ! Calculate cross-section * ve, weighted to macro-weight
          el_ke = el_e - m0*c**2
          sigma_ci_outer = find_value_from_table_1d_coll(el_ke, &
              sample_el_in, &
              species_list(recombined_species)%coll_ion_incident_ke, &
              species_list(recombined_species)%ci_outer_cross_sec, &
              last_table_state)
          sigma_v_sum = sigma_v_sum + weight * el_v * sigma_ci_outer

          current=>current%next
        END DO
      END IF
    END DO

    ! Calculate alpha = <sigma * v>
    alpha_j = sigma_v_sum / e_weight_sum

    ! Save electron number density
    el_ne = e_weight_sum / (dx * dy * dz)

  END SUBROUTINE calculate_alpha_and_ne



  FUNCTION rate_3br(ispecies, e_temp, num_dens_e)

    ! Calculate the three-body recombination rate. From Kremp's "Quantum 
    ! statistics of non-ideal plasmas", we take:
    ! (eq 2.232): dne/dt = (ne*ni*alpha_CI - ne^2*ni*beta_3BR + radiative terms)
    ! (eq 8.47): beta_j = gj * lambda_e^3 * alpha_j * exp(Ieff / kT)
    ! In a PIC code, we assume weakly coupled plasma, so Ieff ~ I
    ! j refers to the nl state capturing the free e- in a groundstate capture
    ! To find j, seek the shell which ionises from the recombined species

    REAL(num) :: rate_3br
    REAL(num), INTENT(IN) :: e_temp, num_dens_e 
    INTEGER, INTENT(IN) :: ispecies
    REAL(num) :: lambda_e, gj, ionise_energy, beta_j
    INTEGER :: ion_l

    ! Thermal de Broglie wavelength - Kremp (2.18)
    lambda_e = SQRT(2.0_num * pi * h_bar**2 / (m0 * kb * e_temp))

    ! Degeneracy of state j
    ion_l = species_list(species_list(ispecies)%recombine_to_species)%l
    gj = 2.0_num * (2.0_num * REAL(ion_l, num) + 1.0_num)

    ! Ionisation energy of state j
    ionise_energy = species_list(species_list(ispecies)%recombine_to_species)&
        %ionisation_energy

    beta_j = gj * lambda_e**3 * alpha_j * exp(ionise_energy/(kb * e_temp))
    rate_3br = num_dens_e * beta_j

  END FUNCTION



  FUNCTION find_value_from_table_1d_recombine(x_in, nx, x, values, state)

    ! For a pair of arrays, x and values, of size nx, this function returns the
    ! interpolated value of "values" corresponding to x_in in the x array. This
    ! uses linear interpolation, unlike in photons.F90 which is logarithmic

    REAL(num) :: find_value_from_table_1d_recombine
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(:), values(:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    REAL(num) :: fx, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.
    LOGICAL :: found

    ! Scan through x to find correct row of table
    i1 = state%ix1
    i2 = state%ix2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = x(i1) - x_in
      xdif2 = x(i2) - x_in
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%ix1+1,nx)
        i2 = MIN(state%ix2+1,nx)
        xdif1 = x(i1) - x_in
        xdif2 = x(i2) - x_in
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%ix1-1,1)
          i2 = MAX(state%ix2-1,1)
          xdif1 = x(i1) - x_in
          xdif2 = x(i2) - x_in
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through x to find correct row of table
      i1 = 1
      i2 = nx
      xdif1 = x(i1) - x_in
      xdif2 = x(i2) - x_in
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = x(im) - x_in
          IF (xdif1 * xdifm < 0) THEN
            i2 = im
          ELSE
            i1 = im
            xdif1 = xdifm
          END IF
          IF (i2 - i1 == 1) EXIT
        END DO
        found = .TRUE.
      END IF
    END IF

    IF (found) THEN
      ! Interpolate in x to find fraction between the cells
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
      state%ix1 = i1
      state%ix2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Electron temperature outside range of known values for ', &
            'recombination data'
        PRINT*,'Using rate at limit of table. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      ! Our x_in value falls outside of the x array - truncate the value
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
      state%ix1 = 1
      state%ix2 = 1
    END IF

    ! Corresponding number from value array, a fraction fx between i1 and i2
    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)
    find_value_from_table_1d_recombine = value_interp
    state%x = x_in
    state%val1d = find_value_from_table_1d_recombine

  END FUNCTION find_value_from_table_1d_recombine

END MODULE recombination
  