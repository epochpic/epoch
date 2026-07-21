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

MODULE calc_df

  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array, species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN), OPTIONAL :: species

    CALL processor_summation_bcs(data_array, ng, species=species)

  END SUBROUTINE calc_boundary



  SUBROUTINE calc_mass_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      fac = io_list(ispecies)%weight
      wdata = part_m * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_m * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_m * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          wdata = current%particle_energy * part_w
#else
          wdata = 0.0_num
#endif
        END IF

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_cou_log(data_array, current_species, direction)

    ! Method taken from: "Lee, Y. T., & More, R. M. (1984). The Physics of
    ! fluids, 27(5), 1273-1286", Section IIIB. The equations present have been
    ! converted to SI units for the EPOCH implementation.
    !
    ! Populates the array: data_array with the Coulomb logarithm values in each
    ! cell for a species with ID current_species
    !
    ! Note that while Debye length, and the uncertainty_dx is calculated cell by
    ! cell, the other impact parameters involve quantities averaged over the
    ! whole simulation window.

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER :: ispecies, use_species, ix, iy
    TYPE(particle), POINTER :: current
    REAL(num), ALLOCATABLE :: species_charge(:), species_weight(:)
    REAL(num), ALLOCATABLE :: num_dens(:,:), num_dens_all(:,:)
    REAL(num), ALLOCATABLE :: mean_ke(:,:), temperature(:,:)
    REAL(num), ALLOCATABLE :: debye_length(:,:), atomic_spacing(:,:)
    REAL(num), ALLOCATABLE :: classical_b0(:,:), uncertain_dx(:,:)
    REAL(num) :: total_charge, total_mass, n_part
    REAL(num) :: part_q, part_m, part_w, b_min, b_max
    REAL(num) :: current_charge, current_mass, mean_q1q2, p_rel

    ALLOCATE(num_dens(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(mean_ke(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(temperature(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(num_dens_all(1-ng:nx+ng, 1-ng:ny+ng))

    ALLOCATE(debye_length(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(atomic_spacing(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(classical_b0(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(uncertain_dx(1-ng:nx+ng, 1-ng:ny+ng))

    ! Estimate the mean charge and mass of each particle species, along with the
    ! total weight of each
    ALLOCATE(species_charge(n_species))
    ALLOCATE(species_weight(n_species))

    ! If no species ID has been provided, default to the first species
    IF (current_species <= 0) THEN
      use_species = 1
    ELSE
      use_species = current_species
    END IF

    ! Loop over all particle species
    DO ispecies = 1, n_species
      total_charge = 0.0_num
      total_mass = 0.0_num
      n_part = 0.0_num
      current => species_list(ispecies)%attached_list%head

      ! Loop over all particles in the particle species
      DO WHILE(ASSOCIATED(current))
        ! Determine real-particle charge and mass
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q = current%charge
        part_m = current%mass
#else
        part_q = species_list(ispecies)%charge
        part_m = species_list(ispecies)%mass
#endif
        ! Determine macro-particle weight
#ifdef PER_SPECIES_WEIGHT
        part_w = species_list(ispecies)%weight
#else
        part_w = current%weight
#endif
        ! Sum weights, and weighted charges and masses
        n_part = n_part + part_w
        total_charge = total_charge + part_q * part_w
        IF (ispecies == use_species) total_mass = total_mass &
            + part_m * part_w
        current => current%next
      END DO

      ! Save average charge and mass
      species_charge(ispecies) = total_charge / n_part
      IF (ispecies == use_species) current_mass = total_mass / n_part

      ! Save charge for species of interest
      IF (ispecies == use_species) THEN
        current_charge = species_charge(ispecies)
      END IF

      ! Save total species weight
      species_weight(ispecies) = n_part
    END DO

    ! Estimate for the product of charges between the current species, and
    ! all other species
    mean_q1q2 = 0.0_num
    n_part = 0.0_num
    ! Weighted sum of all q2 (save in mean_q1q2 for now)
    DO ispecies = 1, n_species
      mean_q1q2 = mean_q1q2 &
          + ABS(species_charge(ispecies)) * species_weight(ispecies)
      n_part = n_part + species_weight(ispecies)
    END DO
    ! Mean q1*q2
    mean_q1q2 = ABS(species_charge(use_species) * mean_q1q2) / n_part

    ! Loop over all species to deduce parameters for b_max, b_min lengths
    debye_length = 0.0_num
    num_dens_all = 0.0_num
    DO ispecies = 1, n_species
      CALL calc_number_density(num_dens, ispecies)
      CALL calc_temperature(temperature, ispecies)

      ! Sum(T/nq2) - store in Debye length array to save memory
      ! Prevent 1/0 error for empty cells
      DO iy = 1-ng, ny+ng
      DO ix = 1-ng, nx+ng
        IF (num_dens(ix,iy) > 1.0e-15_num) THEN
          debye_length(ix,iy) = debye_length(ix,iy) &
              + temperature(ix,iy) &
              / (num_dens(ix,iy) * (species_charge(ispecies))**2)
        ELSE
          debye_length(ix,iy) = 0.0_num
        END IF
      END DO
      END DO
      ! Calculate total number density
      num_dens_all = num_dens_all + num_dens
    END DO

    ! KE of current species
    CALL calc_ekbar(mean_ke, use_species)

    ! Calculate the Coulomb logarithm in each cell
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng

      ! Debye-Huckel screening length, Lee-More (18) - ignore Fermi correction
      ! Initially, Debye length is the sum of T/nq2 for all species
      IF (debye_length(ix,iy) > 0.0_num) THEN
        debye_length(ix,iy) = SQRT(debye_length(ix,iy) * epsilon0 * kb)
      ELSE
        debye_length(ix,iy) = 0.0_num
      END IF

      ! Interatomic distance, assuming each particle is a sphere
      IF (num_dens_all(ix,iy) > 1.0e-15_num) THEN
        atomic_spacing(ix,iy) = (3.0_num &
            / (4.0_num * pi * num_dens_all(ix,iy)))**(1.0_num/3.0_num)
      ELSE
        atomic_spacing(ix,iy) = 0.0_num
      END IF

      ! Classical distance of closest approach, Lee-More (20), replacing mv2
      ! with 0.5*<KE>. Initially, classical_b0 = <KE>
      IF (mean_ke(ix,iy) > 0.0_num) THEN
        classical_b0(ix,iy) = mean_q1q2 &
            / (8.0_num * pi * epsilon0 * mean_ke(ix,iy))
      ELSE
        classical_b0(ix,iy) = 0.0_num
      END IF

      ! Uncertainty principle limit, Lee-More (21)
      p_rel = SQRT((mean_ke(ix,iy) + current_mass * c**2)**2 &
          - current_mass**2 * c**4) / c
      IF (p_rel > 0.0_num) THEN
        uncertain_dx(ix,iy) = h_planck / (2.0_num * p_rel)
      ELSE
        uncertain_dx(ix,iy) = 0.0_num
      END IF

      ! Coulomb logarithm, Lee-More (17)
      b_min = MAX(uncertain_dx(ix,iy), classical_b0(ix,iy))
      b_max = MAX(debye_length(ix,iy), atomic_spacing(ix,iy))
      IF (b_min > 0.0_num .AND. b_max >= b_min) THEN
        data_array(ix,iy) = 0.5_num * LOG(1.0_num + (b_max/b_min)**2)
      ELSE
        data_array(ix,iy) = 0.0_num
      END IF
    END DO
    END DO

    DEALLOCATE(num_dens, mean_ke, temperature, num_dens_all)
    DEALLOCATE(debye_length, atomic_spacing, classical_b0, uncertain_dx)
    DEALLOCATE(species_charge, species_weight)

  END SUBROUTINE calc_cou_log



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    xfac = c * dy
    yfac = c * dx
    zfac = c * dx * dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          fac = c / current%particle_energy
          part_ux = current%part_p(1) * fac
          part_uy = current%part_p(2) * fac
          part_uz = current%part_p(3) * fac
          gamma_rel = 1.0_num
          wdata = current%particle_energy * part_w
#else
          gamma_rel = 1.0_num
          wdata = 0.0_num
#endif
        END IF

        SELECT CASE(direction)
        CASE(-c_dir_x)
          ! negative flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_x)
          ! positive flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_y)
          ! negative flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_y)
          ! positive flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_z)
          ! negative flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_z)
          ! positive flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER :: ix, iy
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iy = 1, ny
      DO ix = 1, nx
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy  ))
        ez_cc = ez(ix, iy)
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
      END DO
      END DO
    CASE(c_dir_y)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy  ))
        ez_cc = ez(ix, iy)
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
      END DO
      END DO
    CASE(c_dir_z)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy))
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy))
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy))
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy))
        data_array(ix,iy) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
      END DO
      END DO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx, vol
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    vol = dx * dy
    idx = 1.0_num / vol

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum &
          .AND. io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      IF (io_list(ispecies)%background_species) THEN
        data_array(1:nx, 1:ny) = data_array(1:nx, 1:ny) &
            + io_list(ispecies)%background_density(1:nx, 1:ny) * vol
      ELSE
        current => io_list(ispecies)%attached_list%head
        wdata = io_list(ispecies)%weight

        DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
          wdata = current%weight
#endif

#include "particle_to_grid.inc"

          DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          END DO
          END DO

          current => current%next
        END DO
      END IF
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_ppc(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      DO WHILE (ASSOCIATED(current))
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + 1.0_num

        current => current%next
      END DO
    END DO

  END SUBROUTINE calc_ppc



  SUBROUTINE calc_average_weight(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    part_count = 0.0_num

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + wdata
        part_count(cell_x,cell_y) = part_count(cell_x,cell_y) + 1.0_num

        current => current%next
      END DO
    END DO

    data_array = data_array / MAX(part_count, c_tiny)

    DEALLOCATE(part_count)

  END SUBROUTINE calc_average_weight



  SUBROUTINE calc_temperature(sigma, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: dof, wdata
    INTEGER :: dir

#include "particle_head.inc"

    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    IF (PRESENT(direction)) THEN
      dir = direction
      dof = 1.0_num
    ELSE
      dir = -1
      dof = 3.0_num
    END IF

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      part_w = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        SELECT CASE(dir)
          CASE(c_dir_x)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanx(cell_x+ix, cell_y+iy) = &
                  meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_y)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meany(cell_x+ix, cell_y+iy) = &
                  meany(cell_x+ix, cell_y+iy) + gf * part_pmy
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_z)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanz(cell_x+ix, cell_y+iy) = &
                  meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE DEFAULT
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanx(cell_x+ix, cell_y+iy) = &
                  meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
              meany(cell_x+ix, cell_y+iy) = &
                  meany(cell_x+ix, cell_y+iy) + gf * part_pmy
              meanz(cell_x+ix, cell_y+iy) = &
                  meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
        END SELECT
        current => current%next
      END DO

      SELECT CASE(dir)
        CASE(c_dir_x)
          CALL calc_boundary(meanx, ispecies)
        CASE(c_dir_y)
          CALL calc_boundary(meany, ispecies)
        CASE(c_dir_z)
          CALL calc_boundary(meanz, ispecies)
        CASE DEFAULT
          CALL calc_boundary(meanx, ispecies)
          CALL calc_boundary(meany, ispecies)
          CALL calc_boundary(meanz, ispecies)
      END SELECT

      CALL calc_boundary(part_count, ispecies)
    END DO

    SELECT CASE(dir)
      CASE(c_dir_x)
        CALL calc_boundary(meanx)
      CASE(c_dir_y)
        CALL calc_boundary(meany)
      CASE(c_dir_z)
        CALL calc_boundary(meanz)
      CASE DEFAULT
        CALL calc_boundary(meanx)
        CALL calc_boundary(meany)
        CALL calc_boundary(meanz)
    END SELECT
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    ! Restore ghost cell values for means
    SELECT CASE(dir)
      CASE(c_dir_x)
        CALL field_bc(meanx, ng)
      CASE(c_dir_y)
        CALL field_bc(meany, ng)
      CASE(c_dir_z)
        CALL field_bc(meanz, ng)
      CASE DEFAULT
        CALL field_bc(meanx, ng)
        CALL field_bc(meany, ng)
        CALL field_bc(meanz, ng)
    END SELECT

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      part_w = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        SELECT CASE(dir)
          CASE(c_dir_x)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmx - meanx(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_y)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmy - meany(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_z)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmz - meanz(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE DEFAULT
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
                    + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
                    + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
        END SELECT
        current => current%next
      END DO
      CALL calc_boundary(sigma, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / dof

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx, root
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    IF (.NOT. PRESENT(direction)) THEN
      IF (rank == 0) THEN
        PRINT*, 'Error: No direction argument supplied to ', &
            'calc_per_species_current'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        SELECT CASE (direction)
          CASE(c_dir_x)
            wdata = wdata * part_px * root
          CASE(c_dir_y)
            wdata = wdata * part_py * root
          CASE(c_dir_z)
            wdata = wdata * part_pz * root
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    idx = c * idx
    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_average_momentum(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    ! The data to be weighted onto the grid
    REAL(num) :: wdata, weight
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    IF (.NOT. PRESENT(direction)) THEN
      IF (rank == 0) THEN
        PRINT*, 'Error: No direction argument supplied to calc_average_momentum'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))

    part_count = 0.0_num
    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      weight = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifndef PER_SPECIES_WEIGHT
        weight = current%weight
#endif
        wdata = weight * current%part_p(direction)

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * weight
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    data_array = data_array / MAX(part_count, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(part_count)

  END SUBROUTINE calc_average_momentum



  SUBROUTINE calc_total_energy_sum(per_species)

    LOGICAL, INTENT(IN) :: per_species
    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: part_mc, part_w, fac, gamma_rel, gamma_rel_m1, part_energy
    REAL(num), ALLOCATABLE :: sum_out(:), sum_in(:)
    REAL(num), ALLOCATABLE :: species_energy(:)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j, nsum
    TYPE(particle), POINTER :: current

    particle_energy = 0.0_num
    IF (per_species) THEN
      ALLOCATE(species_energy(n_species))
      nsum = 1 + n_species
    ELSE
      nsum = 2
    END IF

    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%zero_current) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = part_mc * part_w * c
      part_energy = 0.0_num

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          part_energy = part_energy + gamma_rel_m1 * fac
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
        ELSE
          part_energy = part_energy + current%particle_energy * part_w
#endif
        END IF

        current => current%next
      END DO
      IF (per_species) species_energy(ispecies) = part_energy
      particle_energy = particle_energy + part_energy
    END DO

    ! EM field energy
    field_energy = 0.0_num
    DO j = 1, ny
    DO i = 1, nx
      field_energy = field_energy + ex(i,j)**2 + ey(i,j)**2 &
          + ez(i,j)**2 + c2 * (bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
    END DO
    END DO
    field_energy = 0.5_num * epsilon0 * field_energy * dx * dy

    ALLOCATE(sum_out(nsum))
    ALLOCATE(sum_in(nsum))
    sum_out(1) = field_energy
    IF (per_species) THEN
      sum_out(2:1+n_species) = species_energy(:)
    ELSE
      sum_out(2) = particle_energy
    END IF

    CALL MPI_REDUCE(sum_out, sum_in, nsum, mpireal, MPI_SUM, 0, comm, errcode)
    total_field_energy = sum_in(1)
    IF (per_species) THEN
      total_particle_energy_species(:) = sum_in(2:1+n_species)
      total_particle_energy = SUM(sum_in(2:1+n_species))
    ELSE
      total_particle_energy = sum_in(2)
    END IF

  END SUBROUTINE calc_total_energy_sum



  SUBROUTINE calc_initial_current

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: part_jx, part_jy, part_jz
    REAL(num) :: fac, idx, root
    REAL(num) :: sum_in(3), sum_out(3)
    INTEGER :: ispecies, ix, iy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng))

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    idx = 1.0_num / dx / dy

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%zero_current) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_q  = species_list(ispecies)%charge
      fac = species_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)

        part_jx = wdata * part_px * root
        part_jy = wdata * part_py * root
        part_jz = wdata * part_pz * root

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          jx(cell_x+ix, cell_y+iy) = &
              jx(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jx
          jy(cell_x+ix, cell_y+iy) = &
              jy(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jy
          jz(cell_x+ix, cell_y+iy) = &
              jz(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jz
        END DO
        END DO

        current => current%next
      END DO
    END DO

    sum_out(1) = SUM(jx)
    sum_out(2) = SUM(jy)
    sum_out(3) = SUM(jz)
    DEALLOCATE(jx, jy, jz)

    CALL MPI_ALLREDUCE(sum_out, sum_in, 3, mpireal, MPI_SUM, comm, errcode)

    fac = c * idx / nx_global / ny_global

    initial_jx = sum_in(1) * fac
    initial_jy = sum_in(2) * fac
    initial_jz = sum_in(3) * fac

  END SUBROUTINE calc_initial_current

END MODULE calc_df
