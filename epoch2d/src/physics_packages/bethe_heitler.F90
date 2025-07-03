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
!
!-------------------------------------------------------------------------------
!
! This module contains scripts used to calculate the Bethe-Heitler decay of high
! energy photons into electron-positron pairs. These routines are called from
! bremsstrahlung.F90, and describe the optical depth update, e-/e+ energy
! sampling, and angular deflections of the pair particles.
!
! These algorithms have been adapted from Geant4, and details are provided in
! the Geant4 Physics Reference Manual. The model we are using is valid for
! photon energies between 1.5 MeV and 100 GeV.
!
! To save on particle variable count, the bremsstrahlung optical depth will be
! used for Bethe-Heiter pair production, which was previously unused for photons

MODULE bethe_heitler
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE particles
  USE calc_df
  USE setup

  IMPLICIT NONE

  ! Pre-calculate constants to be used by these functions/subroutines
  REAL(num), PRIVATE, PARAMETER :: ubarn = 1.0e-34_num        ! 1 microbarn [m²]
  REAL(num), PRIVATE, PARAMETER :: gev_100 = 1.0e11_num * q0  ! 100 GeV [J]
  REAL(num), PRIVATE, PARAMETER :: mev_50 = 50.0e6_num * q0   ! 50 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: mev_2 = 2.0e6_num * q0     ! 2 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: kev_1500 = 1.5e6_num * q0  ! 1.5 MeV [J]
  REAL(num), PRIVATE, PARAMETER :: two_mc2 = 2.0_num * m0c2
  REAL(num), PRIVATE, PARAMETER :: extrap_fac = 1.0_num/(kev_1500 - two_mc2)

  ! Pre-calculated functions
  REAL(num), PRIVATE :: cdt
  INTEGER, PRIVATE :: z2ind(100)
  REAL(num), ALLOCATABLE, PRIVATE :: screen_fac(:), fz_low(:), fz_high(:)

CONTAINS

  SUBROUTINE init_bethe_heitler(size_z, z_values, z_to_index)

    ! The Bethe-Heitler algorithm uses some functions which are only Z-dependent
    ! so it makes no sense to recalculate these each step if the Z materials
    ! present don't change over a simulation. Hence, we will pre-calculate some
    ! functions and store them in private arrays
    !
    ! Here, z_values is an array containing all the unique z values present in
    ! the simulation, z2index is an array such that z2index(z) returns the index
    ! in z_values corresponding to z, and size_z is the size of z_values
    !
    ! The functions fz_low and fz_high describe screening and Coulomb
    ! corrections

    INTEGER, INTENT(IN) :: z_values(:)
    INTEGER, INTENT(IN) :: z_to_index(100)
    INTEGER, INTENT(IN) :: size_z
    INTEGER :: iz
    REAL(num) :: z_ion, iz13
    REAL(num) :: alpha_z2, alpha_z4, alpha_z6, coulomb_correction

    cdt = c*dt
    z2ind = z_to_index

    ! Allocate pre-calculated arrays, same indexing as z_values
    ALLOCATE(screen_fac(size_z), fz_low(size_z), fz_high(size_z))

    ! Loop over all possible atomic numbers, calculate terms if this number is
    ! present
    DO iz = 1,size_z

        z_ion = REAL(z_values(iz), num)

        ! Useful pre-factors
        iz13 = z_ion**(-1.0_num/3.0_num)
        screen_fac(iz) = 136.0_num * iz13

        ! Screening terms
        fz_low(iz) = 8.0_num/3.0_num * LOG(z_ion)
        alpha_z2 = (alpha*z_ion)**2
        alpha_z4 = alpha_z2**2
        alpha_z6 = alpha_z4*alpha_z2
        coulomb_correction = alpha_z2 * (1.0_num/(1.0_num + alpha_z2) &
            + 0.20206_num - 0.0369_num*alpha_z2 + 0.0083_num*alpha_z4 &
            - 0.002_num*(alpha_z6))
        fz_high(iz) = fz_low(iz) + 8.0_num * coulomb_correction

    END DO

  END SUBROUTINE init_bethe_heitler



  SUBROUTINE bethe_heitler_update_depth

    ! Cycles through all photon species and updates the Bethe-Heitler optical
    ! depths (labelled bremsstrahlung for photons)

    INTEGER :: iz, z_temp, ispecies
    REAL(num), ALLOCATABLE :: grid_num_density_ion(:,:)
    TYPE(particle), POINTER :: photon, next_photon
    REAL(num) :: part_x, part_y
    REAL(num) :: cross_sec, part_ni, delta_opdep

    ALLOCATE(grid_num_density_ion(1-ng:nx+ng, 1-ng:ny+ng))

    ! Calculate the number density of each ion species
    DO iz = 1, n_species
      ! Identify if the charge is greater than 1
      z_temp = species_list(iz)%atomic_no
      IF (z_temp < 1 .OR. z_temp > 100) CYCLE

      CALL calc_number_density(grid_num_density_ion, iz)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each photon species
      DO ispecies = 1, n_species
        ! Only update optical_depth_bremsstrahlung for the photon species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
          ! Cycle through all photons in this species
          photon => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(photon))
            ! Remember the next photon in the list, as the current photon could
            ! be removed
            next_photon => photon%next

            ! Photon energy is calculated on creation, and this energy never
            ! changes in the current version of the code (11/Feb/2021). Check it
            ! is high enough to pair produce.
            IF (photon%particle_energy < two_mc2) THEN
              photon => next_photon
              CYCLE
            END IF

            ! Calculate the cross section at this photon energy
            cross_sec = calc_bh_cross_sec(z_temp, photon%particle_energy)

            ! Get background number density at photon
            part_x = photon%part_pos(1) - x_grid_min_local
            part_y = photon%part_pos(2) - y_grid_min_local
            CALL grid_centred_var_at_particle(part_x, part_y, part_ni, &
                grid_num_density_ion)

            ! Update the photon optical depth
            delta_opdep = part_ni * cross_sec * cdt
            photon%optical_depth_bremsstrahlung = &
                photon%optical_depth_bremsstrahlung - delta_opdep

            ! If photon optical depth drops below 0, create an e-/e+ pair and
            ! remove the photon from the simulation
            IF (photon%optical_depth_bremsstrahlung < 0.0_num) THEN
              CALL generate_pair(photon, z_temp, &
                  bethe_heitler_electron_species, &
                  bethe_heitler_positron_species, ispecies)
            END IF

            photon => next_photon
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE(grid_num_density_ion)

  END SUBROUTINE bethe_heitler_update_depth



  function calc_bh_cross_sec(z_int, e_gamma)

    ! Calculates the Bethe-Heitler cross section in the presence of an ion
    ! species with atomic number z_int, for a photon of energy e_gamma.
    !
    ! Parametrisation taken from Geant4 source code, which acts as a fit to
    ! cross section data and is estimated to have less than 5% error.
    !
    ! For e_gamma above 100 GeV, the cross section is assumed constant. Below
    ! 1.5 MeV, results are extrapolated down to e_gamma = 2*m0*c².

    INTEGER, INTENT(IN) :: z_int
    REAL(num), INTENT(IN) :: e_gamma
    REAL(num) :: z_ion, e_gamma_calc
    REAL(num) :: x_fit, x2_fit, x3_fit, x4_fit, x5_fit
    REAL(num) :: f1_fit, f2_fit, f3_fit
    LOGICAL :: extrapolate_low
    REAL(num) :: calc_bh_cross_sec

    ! Convert integer atomic number to real
    z_ion = REAL(z_int, num)

    ! High/low energy treatment
    extrapolate_low = .FALSE.
    IF (e_gamma > gev_100) THEN
      e_gamma_calc = gev_100
    ELSE IF (e_gamma < kev_1500) THEN
      e_gamma_calc = kev_1500
      extrapolate_low = .TRUE.
    ELSE
      e_gamma_calc = e_gamma
    END IF

    ! Obtain fit variables
    x_fit = LOG(e_gamma_calc/m0c2)
    x2_fit = x_fit *x_fit
    x3_fit = x2_fit*x_fit
    x4_fit = x3_fit*x_fit
    x5_fit = x4_fit*x_fit

    ! Calculate fit functions
    f1_fit = (878.42_num - 1962.5_num*x_fit + 1294.9_num*x2_fit &
        - 200.28_num*x3_fit + 12.575_num*x4_fit - 0.28333_num*x5_fit)
    f2_fit = (-10.342_num + 17.692_num*x_fit - 8.2381_num*x2_fit &
        + 1.3063_num*x3_fit - 0.090815_num*x4_fit &
        + 2.3586e-3_num*x5_fit)
    f3_fit =  (-452.63_num + 1116.1_num*x_fit - 867.49_num*x2_fit &
        + 217.73_num*x3_fit - 20.467_num*x4_fit + 0.65372_num*x5_fit)

    ! Combine for the total cross section
    calc_bh_cross_sec = (z_ion + 1.0_num) &
        * (f1_fit*z_ion + f2_fit*z_ion**2 + f3_fit) * ubarn

    ! Consider extrapolation if photon energy is below 1.5 MeV
    IF (extrapolate_low) THEN
      calc_bh_cross_sec = calc_bh_cross_sec &
          * ((e_gamma - two_mc2) * extrap_fac)**2
    END IF

  END FUNCTION calc_bh_cross_sec



  SUBROUTINE generate_pair(photon, z_int, ielectron, ipositron, iphoton)

    ! This subroutine is called when "photon" undergoes Bethe-Heitler pair
    ! production. Macro-particles for e- and e+ are created, and the energies
    ! and directions are sampled using Geant4 algorithms. A full description of
    ! these equations may be found in the Geant4 Physics Reference Manual, and
    ! in the references therein. The photon is removed from the simulation here

    TYPE(particle), POINTER :: photon
    INTEGER, INTENT(IN) ::  z_int, ielectron, ipositron, iphoton
    TYPE(particle), POINTER :: new_electron, new_positron
    REAL(num) :: e_frac
    REAL(num) :: energy_1, momentum_1, theta_1, phi_1
    REAL(num) :: energy_2, momentum_2, theta_2, phi_2
    REAl(num) :: norm, dir_x, dir_y, dir_z
    REAL(num) :: p_save, e_save, theta_save, phi_save

    ! Create pair at photon position
    CALL create_particle(new_electron)
    CALL create_particle(new_positron)
    new_electron%part_pos = photon%part_pos
    new_positron%part_pos = photon%part_pos

    ! e- and e+ have the same weights as generating photon
    new_electron%weight = photon%weight
    new_positron%weight = photon%weight

    ! Calculate fractional energy split
    e_frac = energy_split(photon, z_int)

    ! Calculate momentum magnitude going to each particle
    energy_1 = e_frac * photon%particle_energy
    energy_2 = (1.0_num - e_frac) * photon%particle_energy
    momentum_1 = SQRT(energy_1**2 - m0c2**2)/c
    momentum_2 = SQRT(energy_2**2 - m0c2**2)/c

    ! Calculate scatter angles of each particle
    CALL sample_theta(energy_1, energy_2, theta_1, theta_2)
    phi_1 = random()*2.0_num*pi
    phi_2 = phi_1 + pi

    ! Calculate momentum direction of generating photon
    norm  = c / photon%particle_energy
    dir_x = photon%part_p(1) * norm
    dir_y = photon%part_p(2) * norm
    dir_z = photon%part_p(3) * norm

    ! Randomly swap properties of 1 and 2 (could be either e+ or e-)
    IF (random() > 0.5_num) THEN
      p_save = momentum_1
      e_save = energy_1
      theta_save = theta_1
      phi_save = phi_1

      momentum_1 = momentum_2
      energy_1 = energy_2
      theta_1 = theta_2
      phi_1 = phi_2

      momentum_2 = p_save
      energy_2 = e_save
      theta_2 = theta_save
      phi_2 = phi_save
    END IF

    ! Assign momenta to particles
    new_electron%part_p(1:3) = momentum_1 * (/ dir_x, dir_y, dir_z /)
    new_positron%part_p(1:3) = momentum_2 * (/ dir_x, dir_y, dir_z /)

    ! Rotate momenta for each particle
    CALL rotate_p(new_electron, COS(theta_1), phi_1, momentum_1)
    CALL rotate_p(new_positron, COS(theta_2), phi_2, momentum_2)

    ! Set particle energies
    new_electron%particle_energy = energy_1
    new_positron%particle_energy = energy_2

    ! Save pair to particle list
    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Remove photon
    CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
        photon)
    CALL destroy_particle(photon)

  END SUBROUTINE generate_pair



  FUNCTION energy_split(photon, z_int)

    ! This subroutine calculates the energy fraction given to each particle in
    ! the e-/e+ pair

    TYPE(particle), POINTER, INTENT(IN) :: photon
    INTEGER, INTENT(IN) :: z_int
    REAL(num) :: split_0, split_1, split_min
    REAL(num) :: screen_min, screen_max, screen_term
    REAL(num) :: f10, f20, screen_func
    REAL(num) :: split_n1, split_n2, samp_cut
    LOGICAL :: reject
    REAL(num) :: r1, r2, r3
    REAL(num) :: f1_test, f2_test, fz_test
    REAL(num) :: energy_split

    ! Energy split which corresponds to one particle being at rest
    split_0 = m0c2 / photon%particle_energy

    ! Low energy sampling is uniform between split_0 and 2 MeV
    IF (photon%particle_energy < mev_2) THEN
      energy_split = (0.5_num - split_0)*random() + split_0
      RETURN
    END IF

    ! Parameters describing atomic screening
    screen_min = screen_fac(z2ind(z_int)) * (4.0_num * split_0)
    CALL get_screening_funcs(photon, screen_min, z_int, f10, f20, screen_func)

    ! Energy split which prevents screening from making the cross section go -ve
    screen_max = EXP((42.24_num - screen_func)/8.368_num) - 0.952_num
    split_1 = 0.5_num * (1.0_num - SQRT(1.0_num - screen_min/screen_max))

    ! Calculate the minimum energy fraction which may be transferred to one
    ! particle in the pair
    split_min = MAX(split_0, split_1)

    ! Term for splitting treatment in accept/reject algorithm
    split_n1 = (0.5_num - split_min)**2 * f10
    split_n2 = 1.5_num * f20
    samp_cut = split_n1/(split_n1 + split_n2)

    ! Perform accept/reject algorithm
    reject = .TRUE.
    DO WHILE (reject)
      r1 = random()
      r2 = random()
      r3 = random()

      ! Let r1 choose the sampling method
      IF (r1 < samp_cut) THEN
        energy_split = 0.5_num - (0.5_num - split_min) * r2**(1.0_num/3.0_num)
      ELSE
        energy_split = split_min + (0.5_num - split_min) * r2
      END IF

      ! Calculate corresponding screening functions
      screen_term = screen_fac(z2ind(z_int)) * split_0 &
          / (energy_split * (1.0_num - energy_split))
      CALL get_screening_funcs(photon, screen_term, z_int, f1_test, f2_test, &
          fz_test)

      ! Test for accept/reject condition
      IF (r1 < samp_cut) THEN
        IF (f1_test/f10 > r3) reject = .FALSE.
      ELSE
        IF (f2_test/f20 > r3) reject = .FALSE.
      END IF
    END DO

  END FUNCTION energy_split



  SUBROUTINE get_screening_funcs(photon, screen_term, iz, f1, f2, fz)

    ! The Geant4 physics reference manual describes three functions, F1, F2 and
    ! f(Z) which are used to sample the energy split in pair creation. These
    ! are calculated for each sampled split in the accept/reject algorithm, and
    ! so are moved to this subroutine for conciceness.

    TYPE(particle), POINTER, INTENT(IN) :: photon
    REAL(num), INTENT(IN) :: screen_term
    REAL(num) :: phi_1, phi_2
    REAL(num), INTENT(OUT) :: f1, f2, fz
    INTEGER :: iz

    ! Auxiliary phi functions
    IF (screen_term <= 1.0_num) THEN
      phi_1 = 20.867_num - 3.242_num*screen_term + 0.625_num*screen_term**2
      phi_2 = 20.209_num - 1.930_num*screen_term - 0.086_num*screen_term**2
    ELSE
      phi_1 = 21.12_num - 4.184_num*LOG(screen_term + 0.952_num)
      phi_2 = phi_1
    END IF

    IF (photon%particle_energy < mev_50) THEN
      fz = fz_low(z2ind(iz))
    ELSE
      fz = fz_high(z2ind(iz))
    END IF

    ! Auxiliary screening functions
    f1 = 3.0_num * phi_1 - phi_2 - fz
    f2 = 1.5_num * phi_1 - 0.5_num * phi_2 - fz

  END SUBROUTINE get_screening_funcs



  SUBROUTINE sample_theta(energy_1, energy_2, theta_1, theta_2)

    ! Let the polar angle be parallel to the photon momentum direction. This
    ! subroutine samples the angles theta_1 and theta_2 for the scattered
    ! electron/positron pair, using an algorithm described by the Geant4 physics
    ! reference manual

    REAL(num), INTENT(IN) :: energy_1, energy_2
    REAL(num), INTENT(OUT) :: theta_1, theta_2
    REAL(num) :: r1, r2, r3, u_samp

    r1 = random()
    r2 = random()
    r3 = random()

    IF (r1 < 0.25_num) THEN
      u_samp = -1.6_num * LOG(r2*r3)
    ELSE
      u_samp = -8.0_num * LOG(r2*r3) / 15.0_num
    END IF

    theta_1 = u_samp * m0c2 / energy_1
    theta_2 = u_samp * m0c2 / energy_2

  END SUBROUTINE sample_theta



  ! Calculates the value of a grid-centred variable part_var stored in the grid
  ! grid_var, averaged over the particle shape for a particle at position
  ! (part_x, part_y)

  SUBROUTINE grid_centred_var_at_particle(part_x, part_y, part_var, &
      grid_var)

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: grid_var(1-ng:,1-ng:)
    REAL(num), INTENT(OUT) :: part_var
    INTEGER :: cell_x1, cell_y1
    REAL(num) :: cell_x_r, cell_y_r
    REAL(num) :: cell_frac_x, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! The following method is lifted from photons.F90 (field_at_particle), for
    ! the cell-centered fields, taking into account the various particle shapes
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#include "bspline3/part_var.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#include "tophat/part_var.inc"
#else
#include "triangle/gx.inc"
#include "triangle/part_var.inc"
#endif
    part_var = fac * part_var

  END SUBROUTINE grid_centred_var_at_particle

#endif
END MODULE bethe_heitler
