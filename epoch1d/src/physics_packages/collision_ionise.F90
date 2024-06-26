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
! This module contains subroutines used to calculate collisional ionisation,
! which provides additional functionality for the subroutines present in
! collisions.F90

MODULE collision_ionise

  USE calc_df
  USE collisions
  USE prefetch

  IMPLICIT NONE

  ! Variables for deriving cross sections
  REAL(num), ALLOCATABLE :: binding_energy(:), bound_ke(:), cross_sec_parts(:,:)
  REAL(num), PARAMETER :: rbeb_const = 2.0_num * pi * a0**2 * alpha**4
  REAL(num) :: mbell_a(3,0:1), mbell_b(3,0:1,7), mbell_m, mbell_lambda(0:1)
  INTEGER :: el_num
  INTEGER :: el_n(100), el_l(100)
  INTEGER :: table_n(29), table_l(29), table_count(29)

  ! Interpolation table variables
  INTEGER, PARAMETER :: sample_el_in = 100
  INTEGER, PARAMETER :: sample_el_out = 20
  TYPE(interpolation_state), SAVE :: last_state

  ! Saved variables for Lorentz transforms
  REAL(num) :: rot_y, cos_y, sin_y, rot_z, cos_z, sin_z, ion_p2, gamma_i, beta_i

  ! Super-cycled time-step
  REAL(num) :: dt_ci

CONTAINS

  SUBROUTINE setup_coll_ionise_tables

    ! Called before the start of the EPOCH PIC loop. Cycles through the species
    ! list, and calls scripts which calculate collisional ionisation tables for
    ! each species which can be ionised. These tables include:
    !
    ! coll_ion_incident_ke: Log-scale range of incident e- kinetic energies
    ! coll_ion_cross_sec: Total collisional ionisation cross sections for e- KE
    ! coll_ion_secondary_ke: Log range of ejected e- KE, for each incident e- KE
    ! coll_ion_secondary_cdf: CDF for each incident and ejected KE pair
    ! coll_ion_mean_bind: Weighted mean of shell binding energy at each KE pair

    INTEGER :: ispecies, i_el

    CALL setup_mbell_tables

    ! n and l quantum numbers in the binding energy tables
    table_n = (/1, 2, 2, 2, 3, 3, 3, 4, 3, 3, 4, 4, 5, 4, 4, 5, 5, 6, 4, 4, 5, &
        5, 6, 6, 7, 5, 5, 6, 6/)
    table_l = (/0, 0, 1, 1, 0, 1, 1, 0, 2, 2, 1, 1, 0, 2, 2, 1, 1, 0, 3, 3, 2, &
        2, 1, 1, 0, 3, 3, 2, 2/)

    ! Number of electrons which can fit in each shell
    DO i_el = 1,29
      IF (table_l(i_el) == 0) THEN
        ! S state degeneracy is 2
        table_count(i_el) = 2
      ELSE
        ! Higher shells are split into two (j = l + 1/2, j = l - 1/2)
        table_count(i_el) = 2*table_l(i_el) + 1
      END IF
    END DO

    ! Loop over all ionising species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%ionise) THEN

        ! Number of orbital electrons in current species
        el_num = NINT(species_list(ispecies)%atomic_no &
            - species_list(ispecies)%charge/q0)

        ! Get binding energies and quantum numbers for each electron in current
        ! species
        ALLOCATE(binding_energy(el_num))
        ALLOCATE(bound_ke(el_num))
        CALL get_electron_data_from_file(species_list(ispecies)%atomic_no, &
            NINT(species_list(ispecies)%charge/q0), binding_energy, bound_ke, &
            el_n, el_l)

        ! Temporary array for cross section contributions of each bound e-, at
        ! each incident e- KE (used for mean binding energy calculation)
        ALLOCATE(cross_sec_parts(sample_el_in, el_num))
        cross_sec_parts = 0.0_num

        ! Calculate collisional ionisation cross section in each species.
        ! Different models are used for low and high atomic numbers
        IF (species_list(ispecies)%atomic_no <= 18) THEN
          CALL calculate_cross_sections_mbell(ispecies)
        ELSE
          CALL calculate_cross_sections_rbeb(ispecies)
        END IF

        ! Create tables for sampling secondary electron kinetic energies
        CALL calculate_secondary_e_ke_tables(ispecies)

        ! Calculate mean binding energy for each incident KE, ejected KE pair
        CALL calculate_mean_binding_energy(ispecies)

        DEALLOCATE(binding_energy, bound_ke, cross_sec_parts)

      END IF
    END DO

  END SUBROUTINE setup_coll_ionise_tables



  SUBROUTINE setup_mbell_tables

    ! Set the parameters of the empirical MBELL (Modified Bell) model for
    ! calculating electron impact ionisation cross-sections. Data was taken from
    ! Table 1 in Haque et al, "Electron impact ionization of M-shell atoms".
    ! Physica Scripta 74.3 (2006).

    ! Parameter A, where mbell_a(n,l) gives A for quantum numbers n,l
    ! Units: 1e-13 * eV^2 * cm^2
    mbell_a(1,0) = 0.525_num
    mbell_a(2,0) = 0.530_num
    mbell_a(2,1) = 0.600_num
    mbell_a(3,0) = 0.130_num
    mbell_a(3,1) = 0.388_num

    ! Parameters B1 to B7, where mbell_b(n,l,i) gives Bi for quantum numbers n,l
    ! Units: 1e-13 * eV^2 * cm^2
    mbell_b(1,0,:) = (/-0.510_num,  0.2000_num,  0.0500_num, -0.025_num, &
        -0.100_num,  0.00_num,  0.00_num/)
    mbell_b(2,0,:) = (/-0.410_num,  0.1500_num,  0.1500_num, -0.200_num, &
        -0.150_num,  0.00_num,  0.00_num/)
    mbell_b(2,1,:) = (/-0.400_num, -0.7100_num,  0.6550_num,  0.425_num, &
        -0.750_num,  0.00_num,  0.00_num/)
    mbell_b(3,0,:) = (/ 0.250_num, -1.5000_num,  2.4000_num,  3.220_num, &
        -3.667_num,  0.00_num,  0.00_num/)
    mbell_b(3,1,:) = (/-0.200_num, -0.2356_num,  0.5355_num,  3.150_num, &
        -8.500_num,  5.05_num,  0.37_num/)

    ! Convert from [1e-13 * eV^2 * cm^2] to [J^2 m^2]
    mbell_a = 1.0e-13_num * mbell_a * q0**2 * 1.0e-4_num
    mbell_b = 1.0e-13_num * mbell_b * q0**2 * 1.0e-4_num

    ! Parameter m has the same value for all shells up to 3p
    mbell_m = 3.000_num

    ! Parameter lambda, where mbell_lambda(l) gives lambda for quantum number l
    mbell_lambda(0:1) = (/ 1.270_num, 0.542_num /)

  END SUBROUTINE setup_mbell_tables



  SUBROUTINE calculate_cross_sections_mbell(ispecies)

    ! This subroutine is called to calculate collisional ionisation cross
    ! sections for the spcies with ID "ispecies". Called for species with atomic
    ! numbers up to 18. Table saved to species_list(ispecies)%coll_ion_cross_sec
    ! at electron energies stored in species_list(ispecies)%coll_ion_incident_ke
    ! Method described in Haque et al, "Electron impact ionization of M-shell
    ! atoms" Physica Scripta 74.3 (2006).

    INTEGER, INTENT(IN) :: ispecies
    INTEGER :: isamp
    INTEGER :: n_qm, l_qm
    INTEGER :: i_el, ik, n_u, mbell_q
    REAL(num) :: min_ke_ev, max_ke_ev, el_ke
    REAL(num) :: mbell_u, mbell_j, sig_beli_b, sig_beli, f_ion
    REAL(num) :: one_j, one_2j, one_u, u_j, u_2j, gryzinski
    REAL(num) :: sig_add, sig_mbell

    ! Kinetic energy range of table. In Haque (2006), agreement between MBELL
    ! and experiment is verified for some elements between ~10 eV to >100 MeV
    min_ke_ev = MINVAL(binding_energy) / q0
    max_ke_ev = 1.0e9_num

    ! Electron kinetic energies [J]
    ALLOCATE(species_list(ispecies)%coll_ion_incident_ke(sample_el_in))

    ! Cross sections corresponding to electron energies [m**2]
    ALLOCATE(species_list(ispecies)%coll_ion_cross_sec(sample_el_in))

    ! Calculate electron energy and cross section at each table sample point
    DO isamp = 1, sample_el_in

      ! Logarithmically space energy sample points between limits
      el_ke = min_ke_ev * q0 * &
          (max_ke_ev / min_ke_ev)**(REAL(isamp-1,num)/REAL(sample_el_in-1,num))
      species_list(ispecies)%coll_ion_incident_ke(isamp) = el_ke

      ! Caclulate cross section contribution from each electron, assuming shells
      ! are filled in the order:
      ! 1s,2s,2p,3s,3p
      sig_mbell = 0.0_num
      DO i_el = 1, el_num

        ! Energy fractions
        mbell_u = el_ke / binding_energy(i_el)
        mbell_j = m0 * c**2 / binding_energy(i_el)

        ! Only consider ionisation in shells where binding energy < incident KE
        IF (mbell_u < 1.0_num) CYCLE

        ! Quantum numbers of current electron
        n_qm = el_n(i_el)
        l_qm = el_l(i_el)

        ! BELL cross section term
        sig_beli_b = 0.0_num
        DO ik = 1, 7
          sig_beli_b = sig_beli_b + mbell_b(n_qm, l_qm, ik) &
              * (1.0_num - 1.0_num / mbell_u)**REAL(ik, num)
        END DO
        sig_beli = (mbell_a(n_qm, l_qm) * LOG(mbell_u) + sig_beli_b) &
            / (binding_energy(i_el) * el_ke)

        ! Gryzinski relativistic factor (Gr)
        one_j = 1.0_num + mbell_j
        one_2j = 1.0_num + 2.0_num * mbell_j
        one_u = 1.0_num + mbell_u
        u_j = mbell_u + mbell_j
        u_2j = mbell_u + 2.0_num * mbell_j
        gryzinski = one_2j / u_2j * (u_j / one_j)**2 &
            * (one_u * u_2j * one_j**2 &
            / (mbell_j**2 * one_2j + mbell_u * u_2j * one_j**2))**1.5_num

        ! Count electrons up to and including the current shells which share n
        ! and l (e.g. 2p and 2p* are considered the same shell here)
        n_u = i_el
        DO WHILE(n_u <= el_num-1)
          IF ((n_qm == el_n(n_u+1)) .AND. (l_qm == el_l(n_u+1))) THEN
            ! n_u+1 is an electron on the same shell
            n_u = n_u + 1
          ELSE
            ! n_u+1 refers to an electron on a higher shell
            EXIT
          END IF
        END DO

        ! Ionic correction factor F_ion
        mbell_q = species_list(ispecies)%atomic_no - n_u
        f_ion = 1.0_num + mbell_m * (REAL(mbell_q, num) &
            / (species_list(ispecies)%atomic_no * mbell_u))**mbell_lambda(l_qm)

        ! Combine to MBELL cross section
        sig_add = f_ion * gryzinski * sig_beli
        sig_mbell = sig_mbell + sig_add

        ! Save cross section contribution from this bound electron
        cross_sec_parts(isamp, i_el) = sig_add
      END DO

      ! Set cross section after summing contributions from all electrons
      species_list(ispecies)%coll_ion_cross_sec(isamp) = sig_mbell
    END DO

  END SUBROUTINE calculate_cross_sections_mbell



  SUBROUTINE calculate_cross_sections_rbeb(ispecies)

    ! This subroutine is called to calculate collisional ionisation cross
    ! sections for the spcies with ID "ispecies". Called for species with atomic
    ! numbers over 18. Incident electron energies and corresponding cross
    ! sections saved to species_list(ispecies)%coll_ion_incident_ke and
    ! species_list(ispecies)%coll_ion_cross_sec respectively.
    !
    ! Method described in Kim et al, "Extension of the binary-encounter-dipole
    ! model to relativistic incident electrons" Phys. Rev. A 62(5) (2000).
    !
    ! This subroutine uses cross section (22) with the pre-factor of (37) from
    ! Kim (2000). This model works best for the inner shells of heavy nuclei.
    ! The approximation B = U (binding energy = average orbital kinetic energy)
    ! has been made, due to a lack of U data

    INTEGER, INTENT(IN) :: ispecies
    INTEGER ::   isamp, i_el
    REAL(num) :: min_ke_ev, max_ke_ev, el_ke
    REAL(num) :: rbeb_t, rbeb_tp, rbeb_bp, rbeb_up
    REAL(num) :: beta_t2, beta_b2, beta_u2, beta2
    REAL(num) :: rbeb_pre, rbeb_fac, log_2bp, rbeb_1_middle, rbeb_1
    REAL(num) :: one_2tp, tp2_frac, ln_t_term, rbeb_2, rbeb_3
    REAL(num) :: sig_add, sig_rbeb

    ! Range of table and number of sample points. In Kim (2000), agreement
    ! between RBEB and experiment is verified for some elements between ~10 eV
    ! to ~10 MeV
    min_ke_ev = MINVAL(binding_energy)/q0
    max_ke_ev = 1.0e9_num

    ! Electron kinetic energies [J]
    ALLOCATE(species_list(ispecies)%coll_ion_incident_ke(sample_el_in))

    ! Cross sections corresponding to electron energies [m**2]
    ALLOCATE(species_list(ispecies)%coll_ion_cross_sec(sample_el_in))

    ! Calculate electron energy and cross section at each table sample point
    DO isamp = 1, sample_el_in

      ! Logarithmically space energy sample points between limits
      el_ke = min_ke_ev * q0 * &
          (max_ke_ev / min_ke_ev)**(REAL(isamp-1,num)/REAL(sample_el_in-1,num))
      species_list(ispecies)%coll_ion_incident_ke(isamp) = el_ke

      ! Caclulate cross section contribution from each electron, assuming shells
      ! are filled in the order:
      ! 1s,2s,2p,3s,3p,4s,3d,4s,3d,4p,5s,4d,5p,6s,4f,5d,6p,7s,5f,6d
      sig_rbeb = 0.0_num
      DO i_el = 1, el_num

        ! Only consider shells with binding energy less than the incident
        ! electron kinetic energy
        IF (el_ke < binding_energy(i_el)) CYCLE

        ! t
        rbeb_t = el_ke / binding_energy(i_el)

        ! t', b', u' (with approximation B = U)
        rbeb_tp = el_ke / (m0*c**2)
        rbeb_bp = binding_energy(i_el) / (m0*c**2)
        rbeb_up = bound_ke(i_el) / (m0*c**2)

        ! beta_t**2, beta_b**2, beta_u**2
        beta_t2 = 1.0_num - 1.0_num/(1.0_num + rbeb_tp)**2
        beta_b2 = 1.0_num - 1.0_num/(1.0_num + rbeb_bp)**2
        beta_u2 = 1.0_num - 1.0_num/(1.0_num + rbeb_up)**2

        ! Averaging pre-factor of (37) [Kim (2000)]
        beta2 = beta_t2 + beta_b2 + beta_u2
        rbeb_pre = 0.5_num * (1.0_num + beta2/beta_t2)

        ! Calculate RBEB cross section terms. Factor of N from (22) in Kim
        ! (2000) is not needed as we loop over electrons, not shells
        rbeb_fac = rbeb_const / (beta2 * rbeb_bp)

        log_2bp = LOG(2.0_num * rbeb_bp)
        rbeb_1_middle = LOG(beta_t2 / (1.0_num - beta_t2)) - beta_t2 - log_2bp
        rbeb_1 = 0.5_num * rbeb_1_middle * (1.0_num - 1.0_num / rbeb_t**2)

        one_2tp = 1.0_num + 2.0_num * rbeb_tp
        tp2_frac = 1.0_num / (1.0_num + 0.5_num * rbeb_tp)**2
        ln_t_term = LOG(rbeb_t) / (rbeb_t + 1.0_num) * one_2tp * tp2_frac
        rbeb_2 = 1.0_num - (1.0_num / rbeb_t) - ln_t_term

        rbeb_3 = rbeb_bp**2 * tp2_frac * (rbeb_t - 1.0_num) / 2.0_num

        ! Combine to RBEB cross section
        sig_add = rbeb_pre * rbeb_fac * (rbeb_1 + rbeb_2 + rbeb_3)
        sig_rbeb = sig_rbeb + sig_add

        ! Save cross section contribution from this bound electron
        cross_sec_parts(isamp, i_el) = sig_add
      END DO

      ! Set cross section after summing contributions from all electrons
      species_list(ispecies)%coll_ion_cross_sec(isamp) = sig_rbeb
    END DO

  END SUBROUTINE calculate_cross_sections_rbeb



  SUBROUTINE calculate_secondary_e_ke_tables(ispecies)

    ! The electron-impact-ionisation-cross-sections are calculated at multiple
    ! incident electron kinetic energies, in_ke. For each in_ke, this subroutine
    ! generates two arrays: one containing a list of secondary electron kinetic
    ! energies, out_ke, the other giving the cumulative distribution function
    ! (CDF) of each out_ke value.
    !
    ! out_ke is sampled logarithmically up to out_ke_max = 0.5*(in_ke-B), where
    ! B is the binding energy, and down to MAX(0.01 eV, 1e-6*Wmax). out_ke is
    ! also given for 0
    !
    ! CDF is calculated using the differential cross sections (DCS) in Kim,
    ! "Extension of the binary-encounter-dipole model to relativistic incident
    ! electrons" Phys. Rev. A 62(5) (2000). The RBED DCS (19) is simplified to
    ! RBEQ using (7) (and following discussion), and then an RBEB DCS is derived
    ! by setting Ni/N = 1. Pre-factor (37) has been applied.

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), ALLOCATABLE :: cdf_unnorm(:)
    INTEGER :: i_in, i_out, i_bnd
    REAL(num) :: in_ke, out_ke, cross_sec
    REAL(num) :: out_ke_max, out_ke_min, log_out_max, log_out_min, delta_log
    REAL(num) :: out_max_shell, rbeb_w, rbeb_t, rbeb_tp, rbeb_bp, rbeb_up
    REAL(num) :: beta_t2, beta_b2, beta_u2, beta2, rbeb_pre, rbeb_fac
    REAL(num) :: inv_tw, inv_w1, inv_t, log_betat_b, inv_1tp2
    REAL(num) :: cdf_sum, cdf_0, cdf_1, cdf_2, cdf_3, cdf_4, cdf_5, cdf_6

    ALLOCATE(species_list(ispecies)%coll_ion_secondary_ke(sample_el_in, &
        sample_el_out))
    ALLOCATE(species_list(ispecies)%coll_ion_secondary_cdf(sample_el_in, &
        sample_el_out))
    ALLOCATE(cdf_unnorm(sample_el_out-1))

    ! Initialise values as first row has in_ke = minval(binding_energy)
    species_list(ispecies)%coll_ion_secondary_ke = 0.0_num
    species_list(ispecies)%coll_ion_secondary_cdf = 1.0_num

    ! Loop over sampled incident electron kinetic energies
    DO i_in = 1, sample_el_in

      ! Incident electron kinetic energy at index i_in
      in_ke = species_list(ispecies)%coll_ion_incident_ke(i_in)

      ! Ignore incident electron energies below the minimum binding energy
      IF (in_ke <= MINVAL(binding_energy)) CYCLE

      ! Calculate range of allowed secondary energies. Minimum value is chosen
      ! to allow the CDF to be well-sampled
      out_ke_max = 0.5_num * (in_ke - MINVAL(binding_energy))
      out_ke_min = MAX(0.01_num * q0, 1.0e-6_num * out_ke_max)

      ! Terms for logarithmic sampling of secondary electron energies
      log_out_max = LOG(out_ke_max)
      log_out_min = LOG(out_ke_min)
      delta_log = (log_out_max - log_out_min) / REAL(sample_el_out - 2, num)

      ! First sampled energy is 0, with corresponding CDF 0
      species_list(ispecies)%coll_ion_secondary_ke(i_in, 1) = 0.0_num
      species_list(ispecies)%coll_ion_secondary_cdf(i_in, 1) = 0.0_num

      ! Loop over remaining secondary electron kinetic energies
      DO i_out = 1, sample_el_out-1

        ! Sample secondary kinetic energy
        out_ke = EXP(log_out_min + REAL(i_out - 1, num) * delta_log)
        species_list(ispecies)%coll_ion_secondary_ke(i_in, i_out+1) = out_ke

        ! Loop over bound electrons for CDF calculation (see Kim (2000))
        cdf_sum = 0
        DO i_bnd = 1, el_num

          ! Only consider shells with binding energy less than the incident
          ! electron kinetic energy
          IF (in_ke < binding_energy(i_bnd)) CYCLE

          ! Only integrate out_ke up to the maximum ke allowed by this shell
          out_max_shell = 0.5_num * (in_ke - binding_energy(i_bnd))
          rbeb_w = MIN(out_ke, out_max_shell) / binding_energy(i_bnd)

          ! t
          rbeb_t = in_ke / binding_energy(i_bnd)

          ! t', b', u' (with approximation B = U)
          rbeb_tp = in_ke / (m0*c**2)
          rbeb_bp = binding_energy(i_bnd) / (m0*c**2)
          rbeb_up = bound_ke(i_bnd) / (m0*c**2)

          ! beta_t**2, beta_b**2, beta_u**2
          beta_t2 = 1.0_num - 1.0_num/(1.0_num + rbeb_tp)**2
          beta_b2 = 1.0_num - 1.0_num/(1.0_num + rbeb_bp)**2
          beta_u2 = 1.0_num - 1.0_num/(1.0_num + rbeb_up)**2

          ! Averaging pre-factor of (37) [Kim (2000)]
          beta2 = beta_t2 + beta_b2 + beta_u2
          rbeb_pre = 0.5_num * (1.0_num + beta2/beta_t2)

          ! Calculate RBEB cross section terms. Factor of N is not needed as we
          ! loop over electrons, not shells
          rbeb_fac = rbeb_const / (beta2 * rbeb_bp)

          ! Derived terms for neatness
          inv_tw = 1.0_num / (rbeb_t - rbeb_w)
          inv_w1 = 1.0_num / (rbeb_w + 1.0_num)
          inv_t = 1.0_num / rbeb_t
          log_betat_b = LOG(beta_t2 / (2.0_num * rbeb_bp * (1.0_num - beta_t2)))
          inv_1tp2 = 1.0_num / (1.0_num + 0.5_num * rbeb_tp)**2

          ! Calculate CDF contribution using a modified form of (19)
          ! [Kim (2000)], see header comment for details
          cdf_0 = rbeb_pre * rbeb_fac
          cdf_1 = 0.5_num * (inv_tw**2 - inv_w1**2 - inv_t**2 + 1.0_num)
          cdf_2 = log_betat_b - beta_t2
          cdf_3 = inv_tw - inv_w1 - inv_t + 1.0_num
          cdf_4 = rbeb_bp**2 * rbeb_w * inv_1tp2
          cdf_5 = LOG(rbeb_t * (rbeb_w + 1.0_num) / (rbeb_t - rbeb_w))
          cdf_6 = (1.0_num + 2.0_num*rbeb_tp)/(rbeb_t + 1.0_num) * inv_1tp2
          cdf_sum = cdf_sum + cdf_0*(cdf_1*cdf_2 + cdf_3 + cdf_4 - cdf_5*cdf_6)
        END DO

        ! Save un-normalised CDF
        cdf_unnorm(i_out) = cdf_sum
      END DO

      ! CDF is normalised by the cross section (equivalent to final CDF)
      cross_sec = cdf_unnorm(sample_el_out-1)
      species_list(ispecies)%coll_ion_secondary_cdf(i_in, 2:sample_el_out) = &
          cdf_unnorm / cross_sec
    END DO

    DEALLOCATE(cdf_unnorm)

  END SUBROUTINE calculate_secondary_e_ke_tables



  SUBROUTINE calculate_mean_binding_energy(ispecies)

    ! Populates a 2D array coll_ion_mean_bind, which gives the mean binding
    ! energy seen by an incident electron during an ionisation event. This value
    ! is calculated for each incident and ejected KE pair, and the average is
    ! performed over all shells which are capable ejecting an e- at this KE.
    !
    ! This ensures the sum of (ejected KE) + (binding energy) never exceeds half
    ! the incident KE, which is required for energy conservation.

    INTEGER, INTENT(IN) :: ispecies
    INTEGER :: i_in, i_out, i_bnd
    REAL(num) :: in_ke, out_ke, max_bind, sum_parts, sum_bind

    ALLOCATE(species_list(ispecies)%coll_ion_mean_bind(sample_el_in, &
        sample_el_out))

    DO i_out = 1, sample_el_out
      DO i_in = 1, sample_el_in

        ! Extract incident and ejected KE values for this i_in, i_out pair
        in_ke = species_list(ispecies)%coll_ion_incident_ke(i_in)
        out_ke = species_list(ispecies)%coll_ion_secondary_ke(i_in, i_out)

        ! A shell cannot produce an ejected electron with KE out_ke if it
        ! exceeds this binding energ, as max(out_ke) = 0.5 * (in_ke - bind)
        max_bind = in_ke - 2.0_num * out_ke

        ! Calculate average binding energy, weighted by cross section
        ! contributions of allowed shells
        sum_parts = 0.0_num
        sum_bind = 0.0_num
        DO i_bnd = 1, el_num
          IF (binding_energy(i_bnd) <= max_bind) THEN
            sum_parts = sum_parts + cross_sec_parts(i_in, i_bnd)
            sum_bind = sum_bind &
                + cross_sec_parts(i_in, i_bnd) * binding_energy(i_bnd)
          END IF
        END DO
        IF (sum_parts > 0.0_num) THEN
          ! Calculate mean binding energy
          species_list(ispecies)%coll_ion_mean_bind(i_in, i_out) = &
              sum_bind / sum_parts
        ELSE
          ! For i_in = 1, in_ke = minval(binding_energy), sum_parts = 0
          species_list(ispecies)%coll_ion_mean_bind(i_in, i_out) = &
              MINVAL(binding_energy)
        END IF
      END DO
    END DO

  END SUBROUTINE calculate_mean_binding_energy



  SUBROUTINE get_electron_data_from_file(atomic_no, ion_state, binding_energy, &
        bound_ke, el_n, el_l)

    ! Populates the array binding_energy with energies taken from the relevant
    ! binding energy table. These files list binding energies [eV] for a given
    ! element, with a line for each ion charge state and a column for each
    ! shell. The binding_energy array is in [J]. Also creates an equivalent
    ! array for the mean orbital kinetic energy, bound_ke. The quantum numbers
    ! of the electrons are also recorded in el_n and el_l
    !
    ! Binding energies for electrons in neutral atoms were taken from Table 2 in
    ! Desclaux (1973) "Relativistic Dirac-Fock expectation values for atoms with
    ! Z = 1 to Z = 120".
    !
    ! Binding energies for electrons in ions were calculated using the method in
    ! Carlson (1970) "Calculated ionization potentials for multiply charged
    ! ions". The formula relied on ionisation energies taken from NIST.
    !
    ! Mean orbital kinetic energies are assumed to be the same as the binding
    ! energies. The user must over-write these tables if they have the data to
    ! do so.

    INTEGER, INTENT(IN) :: atomic_no, ion_state
    REAL(num), INTENT(OUT) :: binding_energy(:), bound_ke(:)
    INTEGER, INTENT(OUT) :: el_n(:), el_l(:)
    CHARACTER(LEN=3) :: z_string
    LOGICAL :: exists, occ_no_from_file
    INTEGER :: occ_no(29)
    INTEGER :: i_file, io, iu, read_ion, i_shell, el_remain, i_el, el_no
    REAL(num), ALLOCATABLE :: be_all_shells(:), u_all_shells(:)

    ! Deduce relevant file
    IF (atomic_no < 10) THEN
      WRITE(z_string, '(I1)') atomic_no
    ELSE IF (atomic_no < 100) THEN
      WRITE(z_string, '(I2)') atomic_no
    ELSE IF (atomic_no == 100) THEN
      WRITE(z_string, '(I3)') atomic_no
    END IF

    ! Check if the tables can be seen, issue warning if not
    INQUIRE(FILE=TRIM(physics_table_location) &
        //'/binding_energy/be_'//z_string, EXIST=exists)
    IF (.NOT.exists) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) ''
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to find the file:'
        WRITE(io,*) TRIM(physics_table_location) &
            //'/binding_energy/be_'//z_string
        WRITE(io,*) ''
      END DO
      CALL abort_code(c_err_io_error)
    END IF

    ! Open binding energy and orbital kinetic energy files
    OPEN(UNIT = lu, &
        FILE = TRIM(physics_table_location)//'/binding_energy/be_'//z_string, &
        STATUS = 'OLD')
    OPEN(UNIT = lu+1, &
        FILE = TRIM(physics_table_location)//'/bound_ke/u_'//z_string, &
        STATUS = 'OLD')

    ! Check if a file is present for occupancy numbers
    INQUIRE(FILE=TRIM(physics_table_location)  // "/occupancy_numbers/occ_no_"&
        // z_string, EXIST=exists)

    ! Open occupancy number file if one is present
    IF (exists) THEN
      occ_no_from_file = .TRUE.
      OPEN(UNIT = lu+2, &
          FILE = TRIM(physics_table_location)//'/occupancy_numbers/occ_no_'//&
          z_string, STATUS = 'OLD')
    ELSE
      occ_no_from_file = .FALSE.
    END IF

    ! Keep reading the file until the correct line is reached (ignore header)
    DO i_file = 1, ion_state+1
      READ(lu,*)
      READ(lu+1,*)
      IF (occ_no_from_file) READ(lu+2,*)
    END DO

    ! Ignore the charge column and read the binding energies and mean oribtial
    ! kinetic energies of the 29 energy shells from 1s to 6d*
    ALLOCATE(be_all_shells(29))
    ALLOCATE(u_all_shells(29))
    READ(lu,*) read_ion, be_all_shells(1:29)
    READ(lu+1,*) read_ion, u_all_shells(1:29)
    CLOSE(lu)
    CLOSE(lu+1)

    ! If the occupancy number has been provided, read the file
    IF (occ_no_from_file) THEN 
      READ(lu+2,*) read_ion, occ_no(1:29)
      CLOSE(lu+2)
    END IF

    ! Convert energies to [J]
    be_all_shells = be_all_shells * q0
    u_all_shells = u_all_shells * q0

    ! Loop over the electrons, deduce the shell and save the binding energy,
    ! bound KE, and quantum numbers
    i_shell = 1
    IF (occ_no_from_file) THEN
      el_remain = occ_no(i_shell)
    ELSE
      el_remain = table_count(i_shell)
    END IF
    el_no = atomic_no - ion_state
    DO i_el = 1, el_no

      ! When we fill the current shell, move onto the next one
      IF (el_remain == 0) THEN
        i_shell = i_shell + 1

        ! Ensure that there are electrons in this shell (binding energy = 0
        ! otherwise)
        DO WHILE(be_all_shells(i_shell) < 1.0e-20)
          i_shell = i_shell + 1
        END DO

        ! Number of vacancies in the new shell
        IF (occ_no_from_file) THEN
          el_remain = occ_no(i_shell)
        ELSE
          el_remain = table_count(i_shell)
        END IF
      END IF

      ! Save shell properties to current electron
      binding_energy(i_el) = be_all_shells(i_shell)
      bound_ke(i_el) = u_all_shells(i_shell)
      el_n(i_el) = table_n(i_shell)
      el_l(i_el) = table_l(i_shell)

      ! Remove vacancy from shell
      el_remain = el_remain - 1
    END DO

    DEALLOCATE(be_all_shells)
    DEALLOCATE(u_all_shells)

  END SUBROUTINE get_electron_data_from_file



  SUBROUTINE run_collisional_ionisation

    ! This subroutine is called by the main PIC loop, and triggers the
    ! collisional ionisation calculation. The species_list is cycled through to
    ! identify pairs of valid electron and ion species which can ionise, and a
    ! subroutine is called to sample the ionisation. Ionised ions, ejected
    ! electrons and ionising incident electrons are moved to the correct species
    ! lists.
    !
    ! Note, this subroutine uses the generic label "ion" to refer to both ions
    ! and atoms

    INTEGER :: ispecies, jspecies
    INTEGER :: electron_species, ion_species, ejected_species, ionised_species
    INTEGER(i8) :: ix
    TYPE(particle_list), POINTER :: p_list
    TYPE(particle_list) :: list_e_ionising, list_i_ionised, list_e_ejected
    LOGICAL :: i_is_ion, i_is_electron, run_coll_ionisation

    dt_ci = dt * REAL(ci_n_step, num)

    CALL create_empty_partlist(list_e_ionising)
    CALL create_empty_partlist(list_i_ionised)
    CALL create_empty_partlist(list_e_ejected)

    ! Shuffle ion species lists ahead of calculation
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%ionise .AND. .NOT. &
          species_list(ispecies)%is_shuffled) THEN
        DO ix = 1, nx
          p_list => species_list(ispecies)%secondary_list(ix)
          CALL shuffle_particle_list_random(p_list)
        END DO ! ix
        species_list(ispecies)%is_shuffled = .TRUE.
      END IF
    END DO

    ! Loop over all species and extract the electron-ion pairs
    DO ispecies = 1, n_species

      ! Currently no support for photon-impact ionisation, so just cycle round
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE

      ! Variables to track the identity of species ispecies
      i_is_ion = .FALSE.
      i_is_electron = .FALSE.
      IF (species_list(ispecies)%electron) THEN
        ! Current species is an electron species
        electron_species = ispecies
        i_is_electron = .TRUE.
      ELSE IF (species_list(ispecies)%ionise) THEN
        ! Current species can be ionised
        ion_species = ispecies
        i_is_ion = .TRUE.
      ELSE
        ! Skip species which are not involved in ionisation
        CYCLE
      END IF

      ! Loop over remaining species for ionisation partners
      DO jspecies = ispecies, n_species

        ! Determine whether ispecies and jspecies are ionisation partners
        run_coll_ionisation = .FALSE.
        IF (species_list(jspecies)%electron .AND. i_is_ion) THEN
          ! jspecies is an electron species which can ionise ispecies
          electron_species = jspecies
          run_coll_ionisation = .TRUE.
        ELSE IF (species_list(jspecies)%ionise .AND. i_is_electron) THEN
          ! jspecies can be ionised by ispecies
          ion_species = jspecies
          run_coll_ionisation = .TRUE.
        ELSE
          ! Skip species which don't form an ionisation pair with ispecies
          CYCLE
        END IF

        ! Get corresponding species ID for new lists
        ejected_species = species_list(ion_species)%release_species
        ionised_species = species_list(ion_species)%ionise_to_species

        ! Perform collisional ionisation for valid electron-ion pairs
        DO ix = 1, nx
        IF (run_coll_ionisation) THEN
            ! Apply collisional ionisation, moving ionising electrons and
            ! ionised ions from the secondary lists to new particle lists, and
            ! filling a list with ejected electrons (collision_ionise.F90)
            CALL collision_ionisation_ei(ion_species, electron_species, &
                species_list(electron_species)%secondary_list(ix), &
                species_list(ion_species)%secondary_list(ix), &
                list_i_ionised, list_e_ionising, list_e_ejected)

            ! Put ions and electrons into respective lists
            CALL append_partlist( &
                species_list(electron_species)%secondary_list(ix), &
                list_e_ionising)
            CALL append_partlist( &
                species_list(ionised_species)%secondary_list(ix), &
                list_i_ionised)
            CALL append_partlist( &
                species_list(ejected_species)%secondary_list(ix), &
                list_e_ejected)
          END IF
        END DO ! ix

      END DO ! jspecies
    END DO ! ispecies

  END SUBROUTINE run_collisional_ionisation



  SUBROUTINE collision_ionisation_ei(ion_species, el_species, list_e_in, &
      list_i_in, list_i_ionised, list_e_ionising, list_e_ejected)

    ! This subroutine reads in particle lists for electrons "list_e_in" and ions
    ! "list_i_in", and performs collisional ionisation between these two
    ! lists. The corresponding species IDs ion_species and el_species are also
    ! required.
    !
    ! The probability of ionisation for an e-/ion pair is calculated using the
    ! cross section tables which were written by setup_coll_ionise_tables at the
    ! start of the simulation. Tables also store the mean binding energy seen at
    ! different electron energies, and CDF tables for ejected electron energies.
    ! Collisions are calculated in the ion rest frame. Even though
    ! macro-electrons are paired to macro-ions to get frame corrections,
    ! electron-energy-loss and ion-ionisation are sampled separately to conserve
    ! total energy and expected secondary electron number
    !
    ! Incident electrons which trigger an ionisation event are removed from
    ! list_e_in and added to "list_e_ionising". Their energy is reduced by the
    ! sum: (ion mean binding energy) + (secondary electron kinetic energy)
    !
    ! Macro-electrons are added to the simulation for orbital electrons ejected
    ! in an ionisation event, and are assigned to "list_e_ejected". e- ejection
    ! is set parallel to the incident e- direction
    !
    ! Ions/atoms which are ionised in an ionisation event are removed from
    ! list_i_in, and are added to "list_i_ionised". Ion recoil is neglected

    INTEGER, INTENT(IN) :: ion_species, el_species
    TYPE(particle_list), INTENT(INOUT) :: list_e_in, list_i_in, list_i_ionised
    TYPE(particle_list), INTENT(INOUT) :: list_e_ionising, list_e_ejected
    TYPE(particle), POINTER :: electron, ion, ejected_electron
    TYPE(particle), POINTER :: next_ion, next_electron
    INTEGER(i8) :: e_count, ion_count, i_ion, i_el, ionised_count
    LOGICAL, ALLOCATABLE :: e_ionised(:), i_ionised(:)
    LOGICAL :: first_pass, ion_at_rest
    LOGICAL :: create_secondary, electron_recoil, e_collide_again
    REAL(num) :: sum_wi, n_i, inv_cell_volume, ion_mass, el_weight, ion_weight
    REAL(num) :: el_p2_i, el_p_mag_i, el_e_i, el_ke_i, el_v_i, p_mag_new
    REAL(num) :: el_p_i(3), el_dir(3), sec_ke_i, eiics, unionised_frac
    REAL(num) :: ionise_prob, sec_no, sec_frac, rand_cdf, mean_bind

    ! Ignore collisions from empty species
    e_count = list_e_in%count
    ion_count = list_i_in%count
    IF (e_count == 0 .OR. ion_count == 0) RETURN

    ! Allocate arrays to track which electrons and ions are involved in
    ! ionisation events
    ALLOCATE(e_ionised(e_count), i_ionised(ion_count))
    e_ionised = .FALSE.
    i_ionised = .FALSE.

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
    inv_cell_volume = 1.0_num / dx
    n_i = sum_wi * inv_cell_volume

    ! Each macro-e collision must pair to a macro-ion. To ensure there is always
    ! a macro-ion available, temporarily make the ion list circular
    list_i_in%tail%next => list_i_in%head
    ion => list_i_in%head
    i_ion = 1

    ! Loop over all macro-electrons
    ionised_count = 0
    electron => list_e_in%head
    DO i_el = 1, e_count

      ! Some macro-electrons must collide with many macro-ions to create the
      ! right number of secondary e-, this loop allows for multiple collisions
      unionised_frac = 1.0_num
      e_collide_again = .TRUE.
      first_pass = .TRUE.
      DO WHILE (e_collide_again)

#ifndef PER_SPECIES_WEIGHT
        el_weight = electron%weight
        ion_weight = ion%weight
#endif

        ! Check if ion is moving or at rest
        ion_p2 = DOT_PRODUCT(ion%part_p, ion%part_p)
        ion_at_rest = ion_p2 < 1.0e-100_num

        ! If the ion is moving, then perform a Lorentz transform to the ion
        ! rest-frame (cross sections are evaluted in rest-frame)
        IF (.NOT. ion_at_rest) THEN
#ifdef PER_PARTICLE_CHARGE_MASS
          ion_mass = ion%mass
#else
          ion_mass = species_list(ion_species)%mass
#endif
          CALL transform_to_ion_frame(ion, ion_mass, electron, el_p_i)
        ELSE
          ! Ion is already at rest
          el_p_i = electron%part_p
        END IF

        ! Calculate the electron KE and speed in the ion rest frame
        el_p2_i = el_p_i(1)**2 + el_p_i(2)**2 + el_p_i(3)**2
        el_p_mag_i = SQRT(el_p2_i)
        el_e_i = SQRT(el_p2_i + (m0*c)**2) * c
        el_ke_i = el_e_i - m0*c**2
        el_v_i = el_p_mag_i * c**2 / el_e_i

        ! Is electron energy higher than the lowest binding energy of the ion?
        IF (el_ke_i < species_list(ion_species)%coll_ion_incident_ke(1)) THEN
          ! Ensure the next ion is a valid target (not previously ionised)
          ion => ion%next
          i_ion = MOD(i_ion, ion_count) + 1
          DO WHILE (i_ionised(i_ion))
            ion => ion%next
            i_ion = MOD(i_ion, ion_count) + 1
          END DO

          ! Consider next electron
          EXIT
        END IF

        ! Find electron impact ionisation cross section (eiics)
        eiics = find_value_from_table_1d_coll(el_ke_i, sample_el_in, &
            species_list(ion_species)%coll_ion_incident_ke, &
            species_list(ion_species)%coll_ion_cross_sec, last_state)

        ! Probability of ionisation from a real electron (not macro-electron)
        ionise_prob = 1.0_num - EXP(- n_i * eiics * el_v_i * dt_ci)

        ! Expected number of secondary electrons created by this macro-electron
        sec_no = el_weight * unionised_frac * ionise_prob

        ! Check if we add a macro-electron to represent ejected electrons
        IF (sec_no > ion_weight) THEN
          ! Weight of macro-ion less than expected number of secondary e-
          ! Ionise macro-ion but find a new target for macro-electron
          create_secondary = .TRUE.

          ! Save fraction of secondary electrons left to make
          sec_frac = ion_weight / sec_no
          unionised_frac = unionised_frac * (1.0_num - sec_frac)
        ELSE
          ! Weight of macro-ion >= expected number of secondary e-
          ! Sample the probability of ionisation for the macro-ion
          create_secondary = (random() <= sec_no / ion_weight)

          ! Electron has created all the secondary e- it will, do not re-collide
          e_collide_again = .FALSE.
        END IF

        ! Separately sample energy loss for macro-electron. This separation is
        ! required for conservation of energy, as a given macro-electron may not
        ! have enough energy to create a secondary macro-electron with the ion
        ! weight, even if the real electrons can ionise some of the real ions in
        ! the macro-particles. Only sample this on the first pass
        IF (random() < ionise_prob .AND. first_pass) THEN
          ! Assume emission of secondary macro-electron with weight of incident
          ! electron (arbitrary weight threshold to reduce sampling)
          electron_recoil = .TRUE.
          e_ionised(i_el) = .TRUE.
        ELSE
          ! Ignore recoil of this electron
          electron_recoil = .FALSE.
        END IF

        ! Compute additional emission variables
        IF (create_secondary .OR. electron_recoil) THEN
          ! Sample energy of secondary macro-electron
          rand_cdf = random()
          sec_ke_i = find_value_from_table_2d_coll(el_ke_i, rand_cdf, &
              sample_el_in, sample_el_out, &
              species_list(ion_species)%coll_ion_incident_ke, &
              species_list(ion_species)%coll_ion_secondary_ke, &
              species_list(ion_species)%coll_ion_secondary_cdf, last_state)

          ! Calculate momentum direction of incident electron
          el_dir = el_p_i / el_p_mag_i
        END IF

        ! Sample the mean binding energy seen by an incident e- of this energy,
        ! averaged over shells which can eject an e- with KE sec_ke_i
        IF (electron_recoil) THEN
          mean_bind = find_value_from_table_2d_coll(el_ke_i, sec_ke_i, &
              sample_el_in, sample_el_out, &
              species_list(ion_species)%coll_ion_incident_ke, &
              species_list(ion_species)%coll_ion_mean_bind, &
              species_list(ion_species)%coll_ion_secondary_ke, last_state)

          ! Interpolation issue when ejected KE is very high, mean_bind can
          ! exceed the maximum binding energy. Replace with maximum
          IF (mean_bind > el_ke_i - 2.0_num*sec_ke_i) THEN
            mean_bind = el_ke_i - 2.0_num*sec_ke_i
          END IF
        END IF

        ! Reduce incident e- energy by (mean binding energy + secondary KE)
        ! Ignore recoil of nucleus
        IF (electron_recoil) THEN
          el_e_i = el_e_i - sec_ke_i - mean_bind

          ! Calculate new electron energy in ion rest frame
          p_mag_new = SQRT((el_e_i/c)**2 - (m0*c)**2)
          el_p_i = p_mag_new * el_dir

          ! Return to simulation frame
          IF (ion_at_rest) THEN
            ! No transformation required
            electron%part_p = el_p_i
          ELSE
            CALL transform_to_sim_frame(el_p_i, el_e_i, electron)
          END IF
        END IF

        ! Mark ion as ionised and create secondary electron macro-particle
        IF (create_secondary) THEN
          i_ionised(i_ion) = .TRUE.
          ionised_count = ionised_count + 1
          CALL generate_secondary_electron(sec_ke_i, el_dir, ejected_electron, &
              ion, ion_at_rest)

          ! Assign macro-electron to ejected particle list
          CALL add_particle_to_partlist(list_e_ejected, ejected_electron)
          NULLIFY(ejected_electron)

          ! Reduce background ion number density
          n_i = n_i - ion_weight * inv_cell_volume
        END IF

        ! Check if any ions remain for ionisation by the current electron
        IF (ionised_count == ion_count) EXIT

        ! Ensure the next ion is a valid target (not previously ionised)
        ion => ion%next
        i_ion = MOD(i_ion, ion_count) + 1
        DO WHILE (i_ionised(i_ion))
          ion => ion%next
          i_ion = MOD(i_ion, ion_count) + 1
        END DO

        ! Only remove incident electron energy once
        IF (first_pass) first_pass = .FALSE.
      END DO

      ! Check if any ions remain for ionisation of the next electrons
      IF (ionised_count == ion_count) EXIT

      electron => electron%next
    END DO

    ! Restore the tail of the ion list
    NULLIFY(list_i_in%tail%next)

    ! Move ionised macro-ions to the ionised-ion particle list
    ion => list_i_in%head
    DO i_ion = 1, ion_count
      next_ion => ion%next
      IF (i_ionised(i_ion)) THEN
        CALL remove_particle_from_partlist(list_i_in, ion)
        CALL add_particle_to_partlist(list_i_ionised, ion)
      END IF
      ion => next_ion
    END DO

    ! Move ionising macro-electrons to the ionising-electron particle list
    electron => list_e_in%head
    DO i_el = 1, e_count
      next_electron => electron%next
      IF (e_ionised(i_el)) THEN
        CALL remove_particle_from_partlist(list_e_in, electron)
        CALL add_particle_to_partlist(list_e_ionising, electron)
      END IF
      electron => next_electron
    END DO

    DEALLOCATE(e_ionised, i_ionised)

  END SUBROUTINE collision_ionisation_ei



  SUBROUTINE transform_to_ion_frame(ion, ion_mass, electron, el_p_i)

    ! Transforms the electron momentum to a frame where the ion motion is in the
    ! x-direction, then performs a Lorentz boost to get the electron momentum in
    ! the ion-rest-frame, el_p_i. Reads in the ion and electron macro-particles,
    ! as well as the ion_mass to determine the ion gamma factor

    TYPE(particle), POINTER, INTENT(IN) :: ion, electron
    REAL(num), INTENT(IN) :: ion_mass
    REAL(num), INTENT(OUT) :: el_p_i(3)
    REAL(num) :: el_px, el_py, el_pz
    REAL(num) :: el_px_rot, el_py_rot, el_pz_rot
    REAL(num) :: el_p2, el_e

    ! Angles for rotation which map the ion velocity into the x direction.
    ! This rotates about the y-axis by rot_y first, then about the z-axis
    ! by rot_z. These are stored as module variables
    rot_y = ATAN2(ion%part_p(3), ion%part_p(1))
    cos_y = COS(rot_y)
    sin_y = SIN(rot_y)
    rot_z = ATAN2(ion%part_p(2), &
        ion%part_p(1) * cos_y + ion%part_p(3) * sin_y)
    cos_z = COS(rot_z)
    sin_z = SIN(rot_z)

    ! Rotate e- momentum into ion frame to simplify Lorentz transform
    el_px = electron%part_p(1)
    el_py = electron%part_p(2)
    el_pz = electron%part_p(3)
    el_px_rot =  (el_px*cos_y + el_pz*sin_y) * cos_z + el_py*sin_z
    el_py_rot = -(el_px*cos_y + el_pz*sin_y) * sin_z + el_py*cos_z
    el_pz_rot = el_pz*cos_y - el_px*sin_y

    ! Ion-frame variables. These are also stored as module variables
    gamma_i = SQRT(ion_p2 / (ion_mass * c)**2 + 1.0_num)
    beta_i = SQRT(1.0_num - 1.0_num / gamma_i**2)

    ! Electron energy in simulation frame
    el_p2 = el_px**2 + el_py**2 + el_pz**2
    el_e = SQRT(el_p2 + (m0*c)**2) * c

    ! Lorentz transform electron momentum to ion rest frame
    el_p_i(1) = gamma_i * (el_px_rot - beta_i * el_e / c)
    el_p_i(2) = el_py_rot
    el_p_i(3) = el_pz_rot

  END SUBROUTINE transform_to_ion_frame



  SUBROUTINE transform_to_sim_frame(el_p_i, el_e_i, electron)

    ! Reverses the transformation performed by transform_to_ion_frame. The
    ! ion-based variables were saved, and do not need to be recalculated. This
    ! subroutine reads the electron momentum and energy, el_p_i and el_e_i,
    ! evaluated in the ion-rest-frame, and saves the simulation-frame electron
    ! momentum to the electron particle

    REAL(num), INTENT(IN) :: el_p_i(3), el_e_i
    TYPE(particle), POINTER, INTENT(INOUT) :: electron
    REAL(num) :: el_px_rot, px_cosz_py_sinz

    ! Inverse Lorentz transform electron momentum back into simulation
    ! frame (py and pz are unchanged by boost)
    el_px_rot = gamma_i * (el_p_i(1) + beta_i * el_e_i / c)

    ! Undo frame rotation
    px_cosz_py_sinz = el_px_rot*cos_z - el_p_i(2)*sin_z
    electron%part_p(1) = px_cosz_py_sinz * cos_y - el_p_i(3) * sin_y
    electron%part_p(2) = el_px_rot * sin_z + el_p_i(2) * cos_z
    electron%part_p(3) = px_cosz_py_sinz * sin_y + el_p_i(3) * cos_y

  END SUBROUTINE transform_to_sim_frame



  SUBROUTINE generate_secondary_electron(sec_ke_i, el_dir, ejected_electron, &
        ion, ion_at_rest)

    ! Creates a new electron macro-particle "ejected_electron", with a
    ! pre-sampled kinetic energy in the ion-rest-frame, ke_i. The momentum
    ! direction in this frame is chosen to match the incident electron
    ! direction, el_dir. The corresponding ion macro-particle defines the
    ! remaining properties. The logical ion_at_rest determines whether a Lorentz
    ! transform is required

    REAL(num), INTENT(IN) :: sec_ke_i, el_dir(3)
    TYPE(particle), POINTER, INTENT(OUT) :: ejected_electron
    TYPE(particle), POINTER, INTENT(IN) :: ion
    LOGICAL, INTENT(IN) :: ion_at_rest
    REAL(num) :: sec_e_i, sec_p_mag_i, sec_p_i(3)

    CALL create_particle(ejected_electron)

    ! Assume secondary electron direction is parallel to incident
    sec_e_i = sec_ke_i + m0*c**2
    sec_p_mag_i = SQRT((sec_e_i/c)**2 - (m0*c)**2)
    sec_p_i = el_dir * sec_p_mag_i

    ! Transform ejected electron momentum to simulation frame
    IF (ion_at_rest) THEN
      ejected_electron%part_p = sec_p_i
    ELSE
      CALL transform_to_sim_frame(sec_p_i, sec_e_i, ejected_electron)
    END IF

    ! Ejected electron gets properties from the macro-ion
    ejected_electron%part_pos = ion%part_pos
#ifndef PER_SPECIES_WEIGHT
    ejected_electron%weight = ion%weight
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    ejected_electron%charge = -q0
    ejected_electron%mass = m0
    ion%charge = ion%charge + q0
    ion%mass = ion%mass - m0
#endif
#ifdef PARTICLE_DEBUG
    ejected_electron%processor = rank
    ejected_electron%processor_at_t0 = rank
#endif

  END SUBROUTINE generate_secondary_electron



  FUNCTION find_value_from_table_1d_coll(x_in, nx, x, values, state)

    ! For a pair of arrays, x and values, of size nx, this function returns the
    ! interpolated value of "values" corresponding to x_in in the x array. This
    ! uses linear interpolation, unlike in photons.F90 which is logarithmic

    REAL(num) :: find_value_from_table_1d_coll
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
        PRINT*,'Argument to "find_value_from_table_1d_coll" in ', &
            'collision_ionise.F90 outside the range of the table.'
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
    find_value_from_table_1d_coll = value_interp
    state%x = x_in
    state%val1d = find_value_from_table_1d_coll

  END FUNCTION find_value_from_table_1d_coll



  FUNCTION find_value_from_table_2d_coll(x_in, p_value, nx, ny, x, y, p_table, &
      state)

    ! For each element of x, we have a 1D array of y values and a 1D array of P
    ! values, such that the 1D array x has corresponding 2D arrays y and
    ! p_table. The 2D arrays y and p_table are of equal size (nx,ny). This
    ! function interpolates in x_in first, creating an interpolated 1D array of
    ! y and p_table values. The second interpolation finds p_value in p_table,
    ! and the function returns the corresponding value in the interpolated 1D y
    ! array. Used for CDF and mean binding energy look-ups.

    REAL(num) :: find_value_from_table_2d_coll
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(:), y(:,:), p_table(:,:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, y_interp, xdif1, xdif2, xdifm
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
      ! Interpolate in x
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
      state%ix1 = i1
      state%ix2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_2d_coll" outside the ', &
            'range of the table.'
        PRINT*,'An incident electron kinetic energy exceeds tabulated values'
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

    find_value_from_table_2d_coll = y_interp
    state%x = x_in
    state%y = p_value
    state%val2d = find_value_from_table_2d_coll

  END FUNCTION find_value_from_table_2d_coll

END MODULE collision_ionise
