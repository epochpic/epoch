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

MODULE bremsstrahlung
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE calc_df
  USE setup
  USE bethe_heitler

  IMPLICIT NONE

CONTAINS

  ! Confirm that the correct particles are present to actually use the
  ! bremsstrahlung routines, and warn the user if they're not. Then call
  ! subroutines to initialise the bremsstrahlung table array, and give all
  ! particles an optical depth for bremsstrahlung emission (if appropriate)

  SUBROUTINE setup_bremsstrahlung_module

    INTEGER :: ispecies, iu, io
    LOGICAL :: found

    IF (produce_bremsstrahlung_photons .AND. .NOT. ic_from_restart) THEN
      ! Look for electron and photon species. This is a secondary sanity check -
      ! full check appears in check_bremsstrahlung_variables
      found = .FALSE.
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%count <= 0) CYCLE
        IF (species_list(ispecies)%species_type == c_species_id_electron &
            .OR. ispecies == bremsstrahlung_photon_species ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO

      ! Print warning if the correct species aren't present
      IF (.NOT. found) THEN
        IF (rank == 0) THEN
          CALL open_status_file
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Electron and photon species are either ', &
                'unspecified or contain no'
            WRITE(io,*) 'particles. Bremsstrahlung routines will do nothing ', &
                'unless electrons are '
            WRITE(io,*) 'injected'
          END DO
          CALL close_status_file
        END IF
      END IF
    END IF

    ! Load the tables for the bremsstrahlung routines
    CALL setup_tables_bremsstrahlung

    IF (.NOT. ic_from_restart) THEN
      DO ispecies = 1, n_species
        CALL initialise_optical_depth(species_list(ispecies))
      END DO
    END IF

  END SUBROUTINE setup_bremsstrahlung_module



  ! Function to ensure all required species are present

  FUNCTION check_bremsstrahlung_variables()

    INTEGER :: check_bremsstrahlung_variables
    INTEGER :: io, iu
    INTEGER :: ispecies
    INTEGER :: first_electron = -1
    INTEGER :: first_positron = -1

    check_bremsstrahlung_variables = c_err_none

    ! No special species required if bremsstrahlung is turned off
    IF (.NOT.use_bremsstrahlung) RETURN

    ! No special species required if we only do radiation reaction
    IF (.NOT.produce_bremsstrahlung_photons) RETURN

    ! Identify if there exists any electron or positron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      END IF

      IF (species_list(ispecies)%species_type == c_species_id_positron &
          .AND. first_positron == -1) THEN
        first_positron = ispecies
      END IF
    END DO

    ! Print warning if there is no electron species
    IF (first_electron < 0 .AND. (first_positron < 0 .AND. positron_brem)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron species specified.'
          WRITE(io,*) 'Specify using "identify:electron".'
          WRITE(io,*) 'Bremsstrahlung routines require at least one ' , &
              'species of electrons.'
          WRITE(io,*) 'Alternatively, use a positron species with', &
              'use_positron_brem.'
        END DO
      END IF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    END IF

    ! Print warning if there is no positron species
    IF (first_positron < 0 .AND. use_bethe_heitler) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No positron species specified.'
          WRITE(io,*) 'Specify using "identify:positron".'
          WRITE(io,*) 'Bethe-Heitler pair production requires at least one ' , &
              'species of positrons.'
        END DO
      END IF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    END IF

#ifdef PHOTONS
    ! photon_species can act as bremsstrahlung_photon_species if no
    ! bremsstrahlung species is defined
    IF (bremsstrahlung_photon_species == -1 &
        .AND. .NOT. photon_species == -1 ) THEN
      bremsstrahlung_photon_species = photon_species
    END IF
#endif

    ! Any electron may act as a Bethe-Heitler species if no electron is
    ! specified
    IF (use_bethe_heitler) THEN
      IF (bethe_heitler_electron_species == -1) THEN
        bethe_heitler_electron_species = first_electron
      END IF

      IF (bethe_heitler_positron_species == -1) THEN
        bethe_heitler_positron_species = first_positron
      END IF
    END IF

    ! Check if there exists a species to populate with bremsstrahlung photons
    IF (bremsstrahlung_photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:brem_photon"'
        END DO
      END IF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    END IF

  END FUNCTION check_bremsstrahlung_variables



  ! Subroutine to generate tables of bremsstrahlung cross sections as a function
  ! of incident electron energy, T. Also computes a table of cumlative
  ! probabilities for each photon energy, k, and each value of T. All tables are
  ! stored in an array (brem_array) as the variable type brem_tables. Each
  ! element of brem_array correpsonds to a different target atomic number.

  SUBROUTINE setup_tables_bremsstrahlung

    INTEGER :: z_temp, iu, io
    INTEGER :: i_species, iz, jz
    INTEGER, ALLOCATABLE :: z_flags(:)
    CHARACTER(LEN=3) :: z_string
    INTEGER :: size_k, size_t
    INTEGER :: i, j
    INTEGER :: buf_index, buf_size
    REAL(num) :: ln_e, ln_e2, ln_e3, ln_e4, ln_e5, ln_e6, ln_e7, pos_fac
    INTEGER, ALLOCATABLE :: int_buf(:)
    REAL(num), ALLOCATABLE :: real_buf(:)

    ! Do all analysis on rank 0, then MPI_BROADCAST to other ranks
    IF (rank == 0) THEN
      ! For each unique atomic number, z, let z_flags(z) = 1
      ALLOCATE(z_flags(100))
      z_flags(:) = 0

      DO i_species = 1, n_species

        z_temp = species_list(i_species)%atomic_no
        IF (z_temp > 0 .AND. z_temp < 101) THEN
          z_flags(z_temp) = 1

        ! We only have tables up to z=100
        ELSE IF (z_temp > 100) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Species with atomic numbers over 100 cannot be ', &
                'modelled using bremsstrahlung libraries, and will be ignored'
          END DO
        END IF
      END DO

      ! Create brem_array to store tables for each atomic number, z. The value
      ! of z_to_index(z) is the brem_array index for tables of atomic number z.
      ! z_values contains the z values present.
      size_brem_array = SUM(z_flags)
      ALLOCATE(brem_array(size_brem_array))
      ALLOCATE(z_values(size_brem_array))
      ALLOCATE(z_to_index(100))
      iz = 1
      DO jz = 1 ,100
        IF (z_flags(jz) == 1) THEN
          z_to_index(jz) = iz
          z_values(iz) = jz
          iz = iz + 1
        ELSE
          z_to_index(jz) = 0
        END IF
      END DO
      DEALLOCATE(z_flags)

      ! Obtain tables for each brem_array value
      DO iz = 1,size_brem_array
        z_temp = z_values(iz)

        ! Open table for given z value
        IF (z_temp < 10) THEN
          WRITE(z_string, '(I1)') z_temp
        ELSE IF (z_temp < 100) THEN
          WRITE(z_string, '(I2)') z_temp
        ELSE IF (z_temp == 100) THEN
          WRITE(z_string, '(I3)') z_temp
        ELSE
          CYCLE
        END IF

        OPEN(unit = lu, file = TRIM(bremsstrahlung_table_location) // '/br' &
            // z_string, status = 'OLD')

        ! Extract data from the table file
        ! File format: (number of k values) (number of e- energy values)
        !              Line of electron eneries [J]
        !              Line of corresponding cross sections [m²]
        !              Table of k values [J]
        !              Table of CDF values for photon emission
        ! In the table, rows correspond to constant electron energy, and columns
        ! at different photon energies, k. Data taken from Geant4 scaled
        ! differential cross section libraries, and processed to include a cut
        ! off k for each electron energy. This ignores a fraction (1e-9) of the
        ! radiated energy produced by the lowest energy photons. This is for
        ! consistency with photons.F90.
        READ(lu,*) size_k, size_t
        ALLOCATE(brem_array(iz)%cross_section(size_t))
        ALLOCATE(brem_array(iz)%cross_section_pos(size_t))
        ALLOCATE(brem_array(iz)%k_table(size_t,size_k))
        ALLOCATE(brem_array(iz)%cdf_table(size_t,size_k))
        ALLOCATE(brem_array(iz)%e_table(size_t))
        brem_array(iz)%size_t = size_t
        brem_array(iz)%size_k = size_k

        ! Read electron energies and cross sections
        READ(lu,*) brem_array(iz)%e_table(1:size_t)
        READ(lu,*) brem_array(iz)%cross_section(1:size_t)

        ! Read table of k values
        DO i = 1, size_t
          READ(lu,*) brem_array(iz)%k_table(i,1:size_k)
        END DO

        ! Read table of cdf values
        DO i = 1, size_t
          READ(lu,*) brem_array(iz)%cdf_table(i,1:size_k)
        END DO

        CLOSE(unit = lu)

        ! Apply scaling factor to electron bremsstrahlung cross sections for
        ! positron bremsstrahlung. This is an analytic approximation to the
        ! ratio of stopping powers between the two particle types:
        ! Kim, L., 1986. Phys. Rev. A, 33(5)
        ! Approximation from the Geant4 Physics Reference Manual
        IF (positron_brem) THEN
          DO i = 1, size_t
            ln_e = LOG(1.0_num + 1.0e6_num * brem_array(iz)%e_table(i) &
                / (REAL(z_temp,num)**2*m0c2))
            ln_e2 = ln_e**2
            ln_e3 = ln_e2*ln_e
            ln_e4 = ln_e3*ln_e
            ln_e5 = ln_e4*ln_e
            ln_e6 = ln_e5*ln_e
            ln_e7 = ln_e6*ln_e
            pos_fac = 1.0_num - EXP(-0.12359_num*ln_e + 0.061274_num*ln_e2 &
                - 0.031516_num*ln_e3 + 7.7446e-3_num*ln_e4 &
                - 1.0595e-3_num*ln_e5 + 7.0568e-5_num*ln_e6 &
                - 1.808e-6_num*ln_e7)
            brem_array(iz)%cross_section_pos(i) = &
                brem_array(iz)%cross_section(i) * pos_fac
          END DO
        ELSE
          brem_array(iz)%cross_section_pos = 0.0_num
        END IF
      END DO

      ! Pack data for broadcast to other ranks
      ! How many brem_tables do we need to load in brem_array?
      CALL MPI_BCAST(size_brem_array, 1, MPI_INTEGER, 0, comm, errcode)

      ! Pack z values and array sizes into int_buf
      buf_size = 3 * size_brem_array + 100
      ALLOCATE(int_buf(buf_size))
      buf_index = 1
      DO i = 1, size_brem_array
        int_buf(buf_index) = z_values(i)
        buf_index = buf_index + 1
      END DO
      DO i = 1, 100
        int_buf(buf_index) = z_to_index(i)
        buf_index = buf_index + 1
      END DO
      DO i = 1, size_brem_array
        int_buf(buf_index) = brem_array(i)%size_k
        buf_index = buf_index + 1
        int_buf(buf_index) = brem_array(i)%size_t
        buf_index = buf_index + 1
      END DO
      CALL MPI_BCAST(int_buf, buf_size, MPI_INTEGER, 0, comm, errcode)

      ! Pack array data into real_buf and broadcast
      buf_size = 0
      DO i = 1, size_brem_array
        buf_size = buf_size &
            + brem_array(i)%size_t * (2 * brem_array(i)%size_k + 3)
      END DO

      ALLOCATE(real_buf(buf_size))
      buf_index = 1
      DO i = 1, size_brem_array
        size_k = brem_array(i)%size_k
        DO j = 1, brem_array(i)%size_t
          real_buf(buf_index:buf_index+size_k-1) = brem_array(i)%cdf_table(j,:)
          buf_index = buf_index + size_k
          real_buf(buf_index:buf_index+size_k-1) = brem_array(i)%k_table(j,:)
          buf_index = buf_index + size_k
          real_buf(buf_index) = brem_array(i)%cross_section(j)
          buf_index = buf_index + 1
          real_buf(buf_index) = brem_array(i)%cross_section_pos(j)
          buf_index = buf_index + 1
          real_buf(buf_index) = brem_array(i)%e_table(j)
          buf_index = buf_index + 1
        END DO
      END DO
      CALL MPI_BCAST(real_buf, buf_size, mpireal, 0, comm, errcode)

    ! All other ranks
    ELSE

      ! First call to get number of z values present in simulation
      CALL MPI_BCAST(size_brem_array, 1, MPI_INTEGER, 0, comm, errcode)

      ALLOCATE(brem_array(size_brem_array))

      ! Second call to get sizes to allocate arrays, and integer z values
      buf_size = 3 * size_brem_array + 100
      ALLOCATE(int_buf(buf_size))
      CALL MPI_BCAST(int_buf, buf_size, MPI_INTEGER, 0, comm, errcode)

      ! Get z arrays
      ALLOCATE(z_values(size_brem_array))
      ALLOCATE(z_to_index(100))
      buf_index = 1
      DO i = 1, size_brem_array
        z_values(i) = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO
      DO i = 1, 100
        z_to_index(i) = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO

      ! Get brem table sizes
      DO i = 1,size_brem_array
        brem_array(i)%size_k = int_buf(buf_index)
        buf_index = buf_index + 1
        brem_array(i)%size_t = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO

      ! Allocate brem tables
      DO i = 1, size_brem_array
        size_t = brem_array(i)%size_t
        size_k = brem_array(i)%size_k
        ALLOCATE(brem_array(i)%cdf_table(size_t,size_k))
        ALLOCATE(brem_array(i)%k_table(size_t,size_k))
        ALLOCATE(brem_array(i)%cross_section(size_t))
        ALLOCATE(brem_array(i)%cross_section_pos(size_t))
        ALLOCATE(brem_array(i)%e_table(size_t))
      END DO

      ! Final MPI call to get data for tables
      buf_size = 0
      DO i = 1,size_brem_array
        buf_size = buf_size &
            + brem_array(i)%size_t * (2 * brem_array(i)%size_k + 3)
      END DO
      ALLOCATE(real_buf(buf_size))
      CALL MPI_BCAST(real_buf, buf_size, mpireal, 0, comm, errcode)

      ! Load real_buf into brem_array
      buf_index = 1
      DO i = 1, size_brem_array
        size_k = brem_array(i)%size_k
        DO j = 1, brem_array(i)%size_t
          brem_array(i)%cdf_table(j,:) = real_buf(buf_index:buf_index+size_k-1)
          buf_index = buf_index + size_k
          brem_array(i)%k_table(j,:) = real_buf(buf_index:buf_index+size_k-1)
          buf_index = buf_index + size_k
          brem_array(i)%cross_section(j) = real_buf(buf_index)
          buf_index = buf_index + 1
          brem_array(i)%cross_section_pos(j) = real_buf(buf_index)
          buf_index = buf_index + 1
          brem_array(i)%e_table(j) = real_buf(buf_index)
          buf_index = buf_index + 1
        END DO
      END DO
    END IF

    ! Initialise Z-dependent Bethe-Heitler functions using this Z list
    IF (use_bethe_heitler) &
        CALL init_bethe_heitler(size_brem_array, z_values, z_to_index)

    DEALLOCATE(int_buf)
    DEALLOCATE(real_buf)

  END SUBROUTINE setup_tables_bremsstrahlung



  ! Loops through all particles in a species and sets a bremsstrahlung optical
  ! depth

  SUBROUTINE initialise_optical_depth(current_species)

    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    ! Randomly sample the optical depth from an exponential distribution
    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth_bremsstrahlung = -LOG(1.0_num - p_tau)
      current => current%next
    END DO

  END SUBROUTINE initialise_optical_depth



  ! Deallocates the shared_data created in setup_tables_bremsstrahlung

  SUBROUTINE deallocate_tables_bremsstrahlung

    INTEGER :: i

    DO i = 1,size_brem_array
      DEALLOCATE(brem_array(i)%k_table)
      DEALLOCATE(brem_array(i)%cdf_table)
      DEALLOCATE(brem_array(i)%cross_section)
      DEALLOCATE(brem_array(i)%cross_section_pos)
      DEALLOCATE(brem_array(i)%e_table)
    END DO

    DEALLOCATE(brem_array)

  END SUBROUTINE deallocate_tables_bremsstrahlung



  ! Updates the optical depth for electrons. This subroutine is responsible for
  ! also calling the function which calculates the optical depth change, and
  ! calling the generate_photon subroutine. This subroutine serves as the main
  ! interface to the bremsstrahlung module for main-loop processes in
  ! epoch3d.F90

  SUBROUTINE bremsstrahlung_update_optical_depth

    INTEGER :: ispecies, iz, z_temp, i, j, k, part_type
    TYPE(particle), POINTER :: current
    REAL(num), ALLOCATABLE :: grid_num_density_electron_temp(:,:,:)
    REAL(num), ALLOCATABLE :: grid_num_density_electron(:,:,:)
    REAL(num), ALLOCATABLE :: grid_temperature_electron_temp(:,:,:)
    REAL(num), ALLOCATABLE :: grid_temperature_electron(:,:,:)
    REAL(num), ALLOCATABLE :: grid_root_temp_over_num(:,:,:)
    REAL(num), ALLOCATABLE :: grid_num_density_ion(:,:,:)
    REAL(num) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: part_x, part_y, part_z, part_e, part_ni, part_v
    REAL(num) :: q_temp, part_root_te_over_ne, plasma_factor

    ! Calculate electron number density and temperature, summed over all
    ! electron species
    IF (use_plasma_screening) THEN
      ALLOCATE(grid_num_density_electron_temp(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(grid_num_density_electron(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(grid_temperature_electron_temp(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(grid_temperature_electron(1-ng:nx+ng, 1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(grid_root_temp_over_num(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      grid_num_density_electron = c_tiny
      grid_temperature_electron = 0.0_num

      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN
          CALL calc_number_density(grid_num_density_electron_temp, ispecies)
          CALL calc_temperature(grid_temperature_electron_temp, ispecies)

          ! Sum the densities, perform a weighted mean to find mean temperature
          grid_num_density_electron = grid_num_density_electron &
              + grid_num_density_electron_temp
          grid_temperature_electron = grid_temperature_electron &
              + grid_temperature_electron_temp*grid_num_density_electron_temp
        END IF
      END DO
      grid_temperature_electron = grid_temperature_electron &
          / grid_num_density_electron

      ! Create a grid of sqrt(Te/ne) values
      DO k = 1-ng, nz+ng
        DO j = 1-ng, ny+ng
          DO i = 1-ng, nx+ng
            IF (grid_num_density_electron(i,j,k) < 1.0e-10_num) THEN
              grid_root_temp_over_num(i,j,k) = 0.0_num
            ELSE IF (grid_temperature_electron(i,j,k) < 1.0e-10_num) THEN
              grid_root_temp_over_num(i,j,k) = 0.0_num
            ELSE
              grid_root_temp_over_num(i,j,k) = &
                  SQRT(grid_temperature_electron(i,j,k) &
                  / grid_num_density_electron(i,j,k))
            END IF
          END DO
        END DO
      END DO

      CALL field_bc(grid_root_temp_over_num, ng)

      DEALLOCATE(grid_num_density_electron_temp)
      DEALLOCATE(grid_temperature_electron_temp)
      DEALLOCATE(grid_num_density_electron)
      DEALLOCATE(grid_temperature_electron)
    ELSE
      ! This will neglect thermal contributions to the bremsstrahlung cross
      ! section
      plasma_factor = 1.0_num
    END IF

    ALLOCATE(grid_num_density_ion(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ! Calculate the number density of each ion species
    DO iz = 1, n_species
      ! Identify if the charge is greater than 1
      z_temp = species_list(iz)%atomic_no

      IF (z_temp < 1 .OR. z_temp > 100) CYCLE

      CALL calc_number_density(grid_num_density_ion, iz)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each species
      DO ispecies = 1, n_species
        ! Only update optical_depth_bremsstrahlung for the electron species, or
        ! positrons if we are using positron bremsstrahlung
        part_type = species_list(ispecies)%species_type
        IF (part_type == c_species_id_electron &
            .OR.(positron_brem .AND. part_type == c_species_id_positron)) THEN
          ! Cycle through all particles in this species
          current => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))
            ! Get electron energy
            part_ux = current%part_p(1) / mc0
            part_uy = current%part_p(2) / mc0
            part_uz = current%part_p(3) / mc0
            gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
            part_e = gamma_rel * m0 * c**2
            current%particle_energy = part_e

            ! Don't update the optical depth if the particle hasn't moved
            IF (gamma_rel - 1.0_num < 1.0e-15_num) THEN
              current => current%next
              CYCLE
            END IF

            ! Get electron speed
            part_v = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
                + current%part_p(3)**2) * c**2 / part_e

            ! Get number density at electron
            part_x = current%part_pos(1) - x_grid_min_local
            part_y = current%part_pos(2) - y_grid_min_local
            part_z = current%part_pos(3) - z_grid_min_local

            CALL grid_centred_var_at_particle(part_x, part_y, part_z, part_ni,&
                grid_num_density_ion)

            ! Update the optical depth for the screening option chosen
            IF (use_plasma_screening) THEN
              ! Obtain extra parameters needed for plasma screening model
              CALL grid_centred_var_at_particle(part_x, part_y, part_z, &
                  part_root_te_over_ne, grid_root_temp_over_num)

              q_temp = species_list(iz)%charge/q0
              plasma_factor = get_plasma_factor(q_temp, REAL(z_temp, num), &
                  part_root_te_over_ne)
            END IF

            current%optical_depth_bremsstrahlung = &
                current%optical_depth_bremsstrahlung &
                - delta_optical_depth(z_temp, part_e, part_v, part_ni, &
                plasma_factor, part_type)

            ! If optical depth dropped below zero generate photon and reset
            ! optical depth
            IF (current%optical_depth_bremsstrahlung < 0.0_num) THEN
              CALL generate_photon(current, z_temp, &
                  bremsstrahlung_photon_species)
              current%optical_depth_bremsstrahlung = reset_optical_depth()
            END IF

            current => current%next
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE(grid_num_density_ion)

    IF (use_plasma_screening) THEN
      DEALLOCATE(grid_root_temp_over_num)
    END IF

    IF (use_bethe_heitler) CALL bethe_heitler_update_depth

  END SUBROUTINE bremsstrahlung_update_optical_depth



  ! Calculate the change in optical depth during this timestep, given by:
  ! (ion number density) * (emission cross section) * (distance traversed)
  ! Here, cross sections are determined from Geant4 look-up tables

  FUNCTION delta_optical_depth(z, part_e, part_v, part_ni, plasma_factor, spec)

    INTEGER, INTENT(IN) :: z, spec
    REAL(num), INTENT(IN) :: part_e, part_v, part_ni, plasma_factor
    INTEGER :: brem_index
    REAL(num) :: cross_sec_val
    REAL(num) :: delta_optical_depth

    brem_index = z_to_index(z)

    ! Switch to interpolate cross section value for correct particle type
    IF (spec == c_species_id_electron) THEN
      cross_sec_val = find_value_from_table_1d(part_e, &
          brem_array(brem_index)%size_t, brem_array(brem_index)%e_table, &
          brem_array(brem_index)%cross_section, brem_array(brem_index)%state) &
          * plasma_factor
    ELSE IF (spec == c_species_id_positron) THEN
      cross_sec_val = find_value_from_table_1d(part_e, &
          brem_array(brem_index)%size_t, brem_array(brem_index)%e_table, &
          brem_array(brem_index)%cross_section_pos, &
          brem_array(brem_index)%state) * plasma_factor
    END IF

    delta_optical_depth = part_ni * cross_sec_val * part_v * dt &
        / photon_weight

  END FUNCTION delta_optical_depth



  ! Calculate the enhancement of the cross section due to plasma screening
  ! z: Charge of ion species
  ! a: Atomic number of ion species
  ! part_root_te_over_ne: Number density of background electron species at the
  !                       electron divided by the temperature of this species,
  !                       evaluated at the electron position.

  FUNCTION get_plasma_factor(z, a, part_root_te_over_ne)

    REAL(num), INTENT(IN) :: a, z
    REAL(num), INTENT(IN) :: part_root_te_over_ne
    REAL(num) :: term1, term2, log_a_third
    REAL(num) :: get_plasma_factor
    REAL(num), PARAMETER :: one_third = 1.0_num / 3.0_num

    log_a_third = one_third * LOG(a)
    term1 = log_plasma_screen_const_1 - log_a_third
    term2 = log_plasma_screen_const_2 + log_a_third &
        + LOG(part_root_te_over_ne + c_tiny)
    get_plasma_factor = 1.0_num + ((z / a)**2 * term2 / term1)
    get_plasma_factor = MAX(1.0_num, get_plasma_factor)

  END FUNCTION get_plasma_factor



  ! Draws a new random number for the exponentially distributed optical depths

  FUNCTION reset_optical_depth()

    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  ! Generates a photon at the same position as the electron "electron". Emission
  ! was triggered due to electron interaction with a nucleus of atomic number z,
  ! and the photon is written to the species with ID iphoton if appropriate.
  ! Photon energy is sampled from bremsstrahlung tables, and momentum is either
  ! parallel to the electron momentum (without brem_scatter), or scattered by
  ! angles drawn from the bremsstrahlung differential cross section (with
  ! brem_scatter). Electron recoil is applied if requested

  SUBROUTINE generate_photon(electron, z, iphoton)

    TYPE(particle), POINTER :: electron
    INTEGER, INTENT(IN) :: iphoton, z
    INTEGER :: brem_index
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: rand_temp, photon_energy, photon_p, part_e
    REAL(num) :: scatter_theta, scatter_phi
    TYPE(particle), POINTER :: new_photon
    LOGICAL :: add_photon

    ! Obtain electron direction (magnitude must be > 0 to prevent 1/0 issues)
    part_e = electron%particle_energy
    mag_p = MAX(SQRT(electron%part_p(1)**2 + electron%part_p(2)**2 &
        + electron%part_p(3)**2), c_tiny)
    dir_x = electron%part_p(1) / mag_p
    dir_y = electron%part_p(2) / mag_p
    dir_z = electron%part_p(3) / mag_p

    ! Determine photon energy
    rand_temp = random()
    brem_index = z_to_index(z)
    photon_energy = find_value_from_table_alt(part_e, rand_temp, &
        brem_array(brem_index)%size_t, brem_array(brem_index)%size_k, &
        brem_array(brem_index)%e_table, brem_array(brem_index)%k_table, &
        brem_array(brem_index)%cdf_table, brem_array(brem_index)%state)

    ! Is this photon to be added into the simulation?
    add_photon = (photon_energy > photon_energy_min_bremsstrahlung &
        .AND. produce_bremsstrahlung_photons)

    ! Ensure photon_energy is a number we can handle at our precision
    IF (photon_energy < c_tiny) photon_energy = c_tiny

    ! Temporarily store photon momentum in the "new_photon" particle
    CALL create_particle(new_photon)
    photon_p = photon_energy / c
    new_photon%part_p(1) = dir_x * photon_p
    new_photon%part_p(2) = dir_y * photon_p
    new_photon%part_p(3) = dir_z * photon_p

    ! Rotate photon according to bremsstrahlung emission angular distribution
    ! Only consider trajectory changes from e- with high enough energy to be
    ! added into the simulation
    IF (photon_energy > photon_energy_min_bremsstrahlung &
        .AND. use_brem_scatter) THEN
      scatter_theta = calc_scatter_theta(part_e)
      scatter_phi = 2.0_num * pi * random()
      CALL rotate_p(new_photon, COS(scatter_theta), scatter_phi, photon_p)
    END IF

    ! Calculate electron recoil (subtract weighted photon momentum)
    IF (use_bremsstrahlung_recoil) THEN
      electron%part_p = electron%part_p - photon_weight * new_photon%part_p
      electron%particle_energy = electron%particle_energy &
          - photon_weight * photon_energy
    END IF

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on.
    IF (add_photon) THEN
      ! Create new photon at the electron position
      new_photon%part_pos = electron%part_pos

#ifdef PHOTONS
      new_photon%optical_depth = reset_optical_depth()
#endif
      new_photon%particle_energy = photon_energy
      new_photon%weight = electron%weight * photon_weight

      ! Add photon to its species
      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)

    ! We are not adding this photon to the simulation, so remove the particle
    ELSE
      CALL destroy_particle(new_photon)
    END IF

  END SUBROUTINE generate_photon



  ! Calculates the angular scatter of photons using the Tsai method outlined in
  ! the Geant4 physics reference manual (release 10.6, section 10.2.1). This
  ! model is accurate for electron energies over 500 keV (as discussed in the
  ! "Comparisons between Tsai, 2BS and 2BN generators" section)

  FUNCTION calc_scatter_theta(part_E)

    REAL(num), INTENT(IN) :: part_E
    REAL(num) :: calc_scatter_theta
    REAL(num) :: gamma_theta_max, gamma_theta, gamma_rel
    REAL(num), PARAMETER :: a1 = 0.625_num
    REAL(num), PARAMETER :: a2 = 1.875_num
    REAL(num), PARAMETER :: border = 0.25_num
    REAL(num) :: r1, r2, r3

    gamma_rel = part_E / m0c2
    gamma_theta_max = pi * gamma_rel

    ! Perform the Tsai algorithm
    DO
      ! We require three random numbers per sample attempt
      r1 = random()
      r2 = random()
      r3 = random()

      ! Random sampling returns a u value, which is equivalent to the product of
      ! the particle gamma factor and the deflection theta
      IF (r1 < border) THEN
        gamma_theta = -LOG(r2 * r3) / a1
      ELSE
        gamma_theta = -LOG(r2 * r3) / a2
      END IF

      ! We require a theta value less than pi
      IF (gamma_theta <= gamma_theta_max) EXIT
    END DO

    calc_scatter_theta = gamma_theta / gamma_rel

  END FUNCTION calc_scatter_theta



  ! For a pair of arrays, x and values, of size nx, this function returns the
  ! interpolated value of "values" corresponding to x_in in the x array. This
  ! uses linear interpolation, unlike in photons.F90 which is logarithmic

  FUNCTION find_value_from_table_1d(x_in, nx, x, values, state)

    REAL(num) :: find_value_from_table_1d
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(:), values(:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    REAL(num) :: fx, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.
    LOGICAL :: found

    ! If x_in matches the x_in of the last state, return the last value
    ! Check fractional difference
    IF (ABS(x_in) > 0.0_num) THEN
      IF (ABS((state%x - x_in)/x_in) < 1e-10_num) THEN
        find_value_from_table_1d = state%val1d
        RETURN
      END IF
    END IF

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
        PRINT*,'Argument to "find_value_from_table_1d" in ', &
            'bremsstrahlung.F90 outside the range of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
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
    find_value_from_table_1d = value_interp
    state%x = x_in
    state%val1d = find_value_from_table_1d

  END FUNCTION find_value_from_table_1d



  ! For each element of x, we have a 1D array of y values and a 1D array of P
  ! values, such that the 1D array x has a corresponding 2D array of y and P.
  ! The P values represent the cumulative probability of finding a particular y
  ! value for a given x. This function returns a y value for a given probability
  ! p_value, at a particular x value x_in, and is used in the code to obtain
  ! photon energies (y) for a given incident electron energy (x).

  FUNCTION find_value_from_table_alt(x_in, p_value, nx, ny, x, y, p_table, &
      state)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(:), y(:,:), p_table(:,:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.
    LOGICAL :: found

    ! If x_in and p_value match previous values, then return previous value
    IF (ABS(x_in) > 0.0_num .AND. ABS(p_value) > 0.0_num) THEN
      ! Calculate fractional differences
      IF (ABS((state%x - x_in)/x_in) < 1e-10_num &
          .AND. ABS((state%y - p_value)/p_value) < 1e-10_num) THEN
        find_value_from_table_alt = state%val2d
        RETURN
      END IF
    END IF

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
      ! Interpolate in x
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
      state%ix1 = i1
      state%ix2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
      state%ix1 = 1
      state%ix2 = 1
    END IF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = state%iy1
    i2 = state%iy2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%iy1+1,ny)
        i2 = MIN(state%iy2+1,ny)
        xdif1 = p_table(ix,i1) - p_value
        xdif2 = p_table(ix,i2) - p_value
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%iy1-1,1)
          i2 = MAX(state%iy2-1,1)
          xdif1 = p_table(ix,i1) - p_value
          xdif2 = p_table(ix,i2) - p_value
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through table row to find p_value
      i1 = 1
      i2 = ny
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = p_table(ix,im) - p_value
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
      state%iy1 = i1
      state%iy2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
      state%iy1 = 1
      state%iy2 = 1
    END IF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = state%iy1
    i2 = state%iy2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%iy1+1,ny)
        i2 = MIN(state%iy2+1,ny)
        xdif1 = p_table(ix,i1) - p_value
        xdif2 = p_table(ix,i2) - p_value
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%iy1-1,1)
          i2 = MAX(state%iy2-1,1)
          xdif1 = p_table(ix,i1) - p_value
          xdif2 = p_table(ix,i2) - p_value
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through table row to find p_value
      i1 = 1
      i2 = ny
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = p_table(ix,im) - p_value
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
      state%iy1 = i1
      state%iy2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
      state%iy1 = 1
      state%iy2 = 1
    END IF

    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x
    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_alt = y_interp
    state%x = x_in
    state%y = p_value
    state%val2d = find_value_from_table_alt

  END FUNCTION find_value_from_table_alt



  ! Shuts down the bremsstrahlung module - to keep consistency with photons.F90

  SUBROUTINE shutdown_bremsstrahlung_module

    CALL deallocate_tables_bremsstrahlung

  END SUBROUTINE

#endif
END MODULE bremsstrahlung
