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

! Relativistic binary collision module
! written by M. G. Ramsay and H. Schmitz
! based on algorithm by Sentoku & Kemp [J Comput Phys, 227, 6846 (2008)]

MODULE collisions

  USE calc_df
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: particle_collisions, setup_collisions
  PUBLIC :: deallocate_collisions, shuffle_particle_list_random

  ABSTRACT INTERFACE

    SUBROUTINE intra_collisions_proto(p_list, mass, charge, weight, &
        dens, log_lambda, user_factor)

      IMPORT particle_list, num

      TYPE(particle_list), INTENT(INOUT) :: p_list
      REAL(num), INTENT(IN) :: mass, charge, weight
      REAL(num), INTENT(IN) :: user_factor
      REAL(num), INTENT(IN) :: dens, log_lambda

    END SUBROUTINE intra_collisions_proto



    SUBROUTINE inter_collisions_proto(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, log_lambda, user_factor )

      IMPORT particle_list, num

      TYPE(particle_list), INTENT(INOUT) :: p_list1
      TYPE(particle_list), INTENT(INOUT) :: p_list2

      REAL(num), INTENT(IN) :: mass1, charge1, weight1
      REAL(num), INTENT(IN) :: mass2, charge2, weight2

      REAL(num), INTENT(IN) :: idens, jdens
      REAL(num), INTENT(IN) :: log_lambda
      REAL(num), INTENT(IN) :: user_factor

    END SUBROUTINE inter_collisions_proto

  END INTERFACE

  PROCEDURE(intra_collisions_proto), POINTER, SAVE :: intra_coll_fn => NULL()
  PROCEDURE(inter_collisions_proto), POINTER, SAVE :: inter_coll_fn => NULL()

  REAL(num) :: dt_coll

  REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
  REAL(num), PARAMETER :: one_m_2eps = 1.0_num - 2.0_num * eps
  REAL(num), PARAMETER :: one_p_2eps = 1.0_num + 2.0_num * eps

  REAL(num), PARAMETER :: e_rest = m0 * c**2
  REAL(num), PARAMETER :: e_rest_ev = e_rest / ev
  REAL(num), PARAMETER :: cc = c**2

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: meanx, meany, meanz, part_count

CONTAINS

  SUBROUTINE particle_collisions

    ! Loops through all pairs of particle species, identifies colliding pairs,
    ! and calls subroutines to apply collisions between them

    INTEGER :: ispecies, jspecies
    INTEGER(i8) :: ix, iy
    TYPE(particle_list), POINTER :: p_list1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: idens, jdens
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jtemp, log_lambda
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: iekbar
    REAL(num) :: user_factor, q1, q2, m1, m2, w1, w2
    LOGICAL :: collide_species
    INTEGER :: collision_type

    ALLOCATE(idens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jdens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jtemp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(log_lambda(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(iekbar(1-ng:nx+ng,1-ng:ny+ng))

    ! Time-step used by collisions routine, for super-cycle speed-up
    dt_coll = dt * REAL(coll_n_step, num)

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE

      ! No support for neutral collisions at present
      IF (ABS(species_list(ispecies)%charge) <= c_tiny) CYCLE

      ! Does ispecies have any permitted collision partners?
      collide_species = .FALSE.
      DO jspecies = ispecies, n_species
        user_factor = coll_pairs(ispecies, jspecies)
        collision_type = coll_pairs_state(ispecies, jspecies)
        IF (collision_type == c_coll_collide) THEN
          collide_species = .TRUE.
          EXIT
        END IF
      END DO

      ! If collisions with ispecies have been switched off, ignore collisions
      IF (.NOT.collide_species) CYCLE

      CALL calc_coll_number_density(idens, ispecies)

      IF (coulomb_log_auto) THEN
        CALL calc_coll_ekbar(iekbar, ispecies)
      END IF

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight

      IF (.NOT. species_list(ispecies)%is_shuffled) THEN
        DO iy = 1, ny
        DO ix = 1, nx
          p_list1 => species_list(ispecies)%secondary_list(ix,iy)
          CALL shuffle_particle_list_random(p_list1)
        END DO ! ix
        END DO ! iy
        species_list(ispecies)%is_shuffled = .TRUE.
      END IF

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type == c_species_id_photon) &
            CYCLE

        ! No support for neutral collisions at present
        IF (ABS(species_list(jspecies)%charge) <= c_tiny) CYCLE

        ! If collisions between ispecies and jspecies are disabled, then cycle
        user_factor = coll_pairs(ispecies, jspecies)
        collision_type = coll_pairs_state(ispecies, jspecies)

        IF (.NOT. collision_type == c_coll_collide) CYCLE

        IF (ispecies /= jspecies) THEN
          CALL calc_coll_number_density(jdens, jspecies)
        END IF

        m2 = species_list(jspecies)%mass
        q2 = species_list(jspecies)%charge
        w2 = species_list(jspecies)%weight

        IF (coulomb_log_auto) THEN
          CALL calc_coll_temperature_ev(jtemp, jspecies)
          IF (ispecies == jspecies) THEN
            log_lambda = calc_coulomb_log(iekbar, jtemp, idens, idens, &
                q1, q1, m1)
          ELSE
            log_lambda = calc_coulomb_log(iekbar, jtemp, idens, jdens, &
                q1, q2, m1)
          END IF
        ELSE
          log_lambda = coulomb_log
        END IF

        DO iy = 1, ny
        DO ix = 1, nx
          IF (ispecies == jspecies) THEN
            CALL intra_coll_fn( &
                species_list(ispecies)%secondary_list(ix,iy), &
                m1, q1, w1, idens(ix,iy), &
                log_lambda(ix,iy), user_factor)
          ELSE
            CALL inter_coll_fn( &
                species_list(ispecies)%secondary_list(ix,iy), &
                species_list(jspecies)%secondary_list(ix,iy), &
                m1, m2, q1, q2, w1, w2, idens(ix,iy), jdens(ix,iy), &
                log_lambda(ix,iy), user_factor)
          END IF
        END DO ! ix
        END DO ! iy
      END DO ! jspecies
    END DO ! ispecies

    DEALLOCATE(idens, jdens, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)
    DEALLOCATE(iekbar)

  END SUBROUTINE particle_collisions



  SUBROUTINE intra_collisions_sk(p_list, mass, charge, weight, &
      dens, log_lambda, user_factor)
    ! Perform collisions between particles of the same species.

    TYPE(particle_list), INTENT(INOUT) :: p_list
    REAL(num), INTENT(IN) :: mass, charge, weight
    REAL(num), INTENT(IN) :: user_factor
    REAL(num), INTENT(IN) :: dens, log_lambda
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor, np
    INTEGER(i8) :: icount, k, pcount
    REAL(num), DIMENSION(3) :: p1, p2 ! Pre-collision momenta
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm ! Normalised momenta
    REAL(num), DIMENSION(3) :: vc ! Velocity of COM frame wrt lab frame
    REAL(num), DIMENSION(3) :: p3, p4
    REAL(num), DIMENSION(3) :: p5, p6
    REAL(num), DIMENSION(3) :: v3, v4
    REAL(num), DIMENSION(3) :: vr, vcr
    REAL(num), DIMENSION(3) :: c1, c2, c3
    REAL(num) :: m1, m2, q1, q2 ! Masses and charges
    REAL(num) :: e1, e2 ! Pre-collision energies
    REAL(num) :: e3, e4
    REAL(num) :: gamma_rel, gamma_rel2, gamma_rel_m1, gamma_rel_r
    REAL(num) :: tvar ! Dummy variable for temporarily storing values
    REAL(num) :: vc_sq, vc_sq_cc, p1_vc, p2_vc, p3_mag
    REAL(num) :: delta, sin_theta, cos_theta, tan_theta_cm, tan_theta_cm2
    REAL(num) :: vrabs, denominator
    REAL(num) :: nu, ran1, ran2
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: e5, e6
    REAL(num) :: w1, w2, wr
#endif

    factor = 0.0_num
    np = 0.0_num

    ! Intra-species collisions
    icount = p_list%count

    ! If there aren't enough particles to collide, then don't bother
    IF (icount <= 1) RETURN

    ! Number of collisions
    pcount = icount / 2 + MOD(icount, 2_i8)

#ifdef PER_SPECIES_WEIGHT
    np = icount * weight
    ! Factor of 2 due to intra species collisions
    ! See Section 4.1 of Nanbu
    factor = user_factor / (pcount * weight * 2.0_num)
#else
    ! temporarily join tail to the head of the list to make it circular
    p_list%tail%next => p_list%head

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
      np = np + current%weight + impact%weight
      factor = factor + MIN(current%weight, impact%weight)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO
    factor = user_factor / factor / 2.0_num
#endif
    ! If possible, use per-species properties
    m1 = mass
    m2 = mass
    q1 = charge
    q2 = charge

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
      m1 = current%mass
      m2 = impact%mass
      q1 = current%charge
      q2 = impact%charge
#endif
      ! Copy all of the necessary particle data into variables with easier to
      ! read names
      p1 = current%part_p
      p2 = impact%part_p
      p1_norm = p1 / mc0
      p2_norm = p2 / mc0

      ! Two stationary particles can't collide, so don't try
      IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
          .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE

      ! Ditto for two particles with the same momentum
      vc = (p1_norm - p2_norm)
      IF (DOT_PRODUCT(vc, vc) < eps) CYCLE

      ! Pre-collision energies
      e1 = c * SQRT(DOT_PRODUCT(p1, p1) + (m1 * c)**2)
      e2 = c * SQRT(DOT_PRODUCT(p2, p2) + (m2 * c)**2)

      ! Velocity of centre-of-momentum (COM) reference frame
      vc = (p1 + p2) * cc / (e1 + e2)
      vc_sq = DOT_PRODUCT(vc, vc)
      vc_sq_cc = vc_sq / cc

      gamma_rel2 = 1.0_num / (1.0_num - vc_sq_cc)
      gamma_rel = SQRT(gamma_rel2)
      gamma_rel_m1 = gamma_rel2 * vc_sq_cc / (gamma_rel + 1.0_num)

      ! Lorentz momentum transform to get into COM frame
      p1_vc = DOT_PRODUCT(p1, vc)
      p2_vc = DOT_PRODUCT(p2, vc)
      tvar = p1_vc * gamma_rel_m1 / (vc_sq + c_tiny)
      p3 = p1 + vc * (tvar - gamma_rel * e1 / cc)
      tvar = p2_vc * gamma_rel_m1 / (vc_sq + c_tiny)
      p4 = p2 + vc * (tvar - gamma_rel * e2 / cc)

      p3_mag = SQRT(DOT_PRODUCT(p3, p3))

      ! Lorentz energy transform
      e3 = gamma_rel * (e1 - p1_vc)
      e4 = gamma_rel * (e2 - p2_vc)
      ! Pre-collision velocities in COM frame
      v3 = p3 * cc / e3
      v4 = p4 * cc / e4

      ! Relative velocity
      tvar = 1.0_num - (DOT_PRODUCT(v3, v4) / cc)
      vr = (v3 - v4) / tvar
      vrabs = SQRT(DOT_PRODUCT(vr, vr))

      ! Collision frequency
      nu = coll_freq(vrabs, log_lambda, m1, m2, q1, q2, dens)

      ! Calculate number of collisions in the timestep
      ! Limit value according to Sentoku & Kemp
      nu = MIN(nu * factor * np * dt_coll, 0.02_num)

      ! New coordinate system to simplify scattering.
      CALL new_coords(vr, c1, c2, c3)

      ! this is to ensure that 0 < ran1 < 1
      ! ran1=0 gives NaN in logarithm
      ! ran1=1 could give positive logarithm due to rounding errors
      ran1 = (1.0_num - 1.0e-10_num) * random() + 0.5e-10_num
      ran2 = 2.0_num * pi * random()

      ! Box Muller method for random from Gaussian distribution,
      ! mean 0, variance of nu
      ! Possible place for speed up by caching the second Box Muller number
      ! and using it later
      ! SQRT(-2.0_num * nu * LOG(ran1)) * COS(ran2)
      delta = SQRT(-2.0_num * nu * LOG(ran1)) * SIN(ran2)
      ran2 = 2.0_num * pi * random()

      ! angle theta in the One Particle at Rest frame
      sin_theta = 2.0_num * delta / (1.0_num + delta**2)
      cos_theta = (1.0_num - delta**2) / (1.0_num + delta**2)

      ! Transform angles from particle j's rest frame to COM frame
      ! Note azimuthal angle (ran2) is invariant under this transformation
      IF (m1 > m2) THEN
        vcr = v3
      ELSE
        vcr = v4
      END IF
      gamma_rel_r = 1.0_num / SQRT(1.0_num - (DOT_PRODUCT(vcr, vcr) / cc))

      denominator = gamma_rel_r * (cos_theta - SQRT(DOT_PRODUCT(vcr, vcr)) &
          / MAX(vrabs, c_tiny))
      IF (ABS(denominator) > SQRT(c_tiny)) THEN
        tan_theta_cm = sin_theta / denominator
        tan_theta_cm2 = tan_theta_cm**2
      ELSE
        tan_theta_cm = c_largest_number
        tan_theta_cm2 = c_largest_number
      END IF

      sin_theta = tan_theta_cm / SQRT(1.0_num + tan_theta_cm2)
      cos_theta = 1.0_num / SQRT(1.0_num + tan_theta_cm2)

      ! Post-collision momenta in COM frame
      p3 = p3_mag * (c1 * cos_theta + c2 * sin_theta * COS(ran2) &
          + c3 * sin_theta * SIN(ran2))
      p4 = -p3

      ! Lorentz momentum transform to get back to lab frame
      tvar = DOT_PRODUCT(p3, vc) * gamma_rel_m1 / vc_sq
      p5 = p3 + vc * (tvar + gamma_rel * e3 / cc)
      tvar = DOT_PRODUCT(p4, vc) * gamma_rel_m1 / vc_sq
      p6 = p4 + vc * (tvar + gamma_rel * e4 / cc)

#ifndef PER_SPECIES_WEIGHT
      w1 = current%weight
      w2 = impact%weight
      wr = w1 / w2

      ! Post collision energies. Only needed for weighted particle correction
      e5 = c * SQRT(DOT_PRODUCT(p5, p5) + (m1 * c)**2)
      e6 = c * SQRT(DOT_PRODUCT(p6, p6) + (m2 * c)**2)

      IF (wr > one_p_2eps) THEN
        CALL weighted_particles_correction(w2 / w1, p1, p5, e1, e5, m1)
      ELSE IF (wr < one_m_2eps) THEN
        CALL weighted_particles_correction(w1 / w2, p2, p6, e2, e6, m2)
      END IF
#endif
      ! Update particle properties
      current%part_p = p5
      impact%part_p = p6

      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO

    ! restore the tail of the list
    NULLIFY(p_list%tail%next)

  END SUBROUTINE intra_collisions_sk



  SUBROUTINE intra_collisions_np(p_list, mass, charge, weight, &
      dens, log_lambda, user_factor)
    ! Perform collisions between particles of the same species.

    TYPE(particle_list), INTENT(INOUT) :: p_list
    REAL(num), INTENT(IN) :: mass, charge, weight
    REAL(num), INTENT(IN) :: user_factor
    REAL(num), INTENT(IN) :: dens, log_lambda
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor
    INTEGER(i8) :: icount, k, pcount
    REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
    REAL(num) :: sinp_cos, sinp_sin, s_prime, s_fac_prime
    REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
    REAL(num) :: p_perp2, p_perp_inv, cell_fac
    REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm
    REAL(num), DIMENSION(3,3) :: mat
    REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq
    REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
    REAL(num) :: m1, m2, q1, q2, dens_23
    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    REAL(num), PARAMETER :: pi_fac = &
                                (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)
    factor = 0.0_num

    ! Intra-species collisions
    icount = p_list%count

    ! If there aren't enough particles to collide, then don't bother
    IF (icount <= 1) RETURN

    ! Number of collisions
    pcount = icount / 2 + MOD(icount, 2_i8)

#ifdef PER_SPECIES_WEIGHT
    ! Factor of 2 due to intra species collisions
    ! See Section 4.1 of Nanbu
    factor = user_factor / (pcount * weight * 2.0_num)
#else
    ! temporarily join tail to the head of the list to make it circular
    p_list%tail%next => p_list%head

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
      factor = factor + MIN(current%weight, impact%weight)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO
    factor = user_factor / factor / 2.0_num
#endif
    ! If possible, use per-species properties
    m1 = mass
    m2 = mass
    q1 = charge
    q2 = charge

    current => p_list%head
    impact => current%next

    ! Per-cell constant factors
    cell_fac = dens**2 * dt_coll * factor * dx * dy
    s_fac = cell_fac * log_lambda / pi4_eps2_c4
    dens_23 = dens**two_thirds
    s_fac_prime = cell_fac * pi_fac / dens_23

    DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
      m1 = current%mass
      m2 = impact%mass
      q1 = current%charge
      q2 = impact%charge
#endif
      p1 = current%part_p / c
      p2 = impact%part_p / c

      p1_norm = p1 / m0
      p2_norm = p2 / m0

      ! Two stationary particles can't collide, so don't try
      IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
          .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE

      ! Ditto for two particles with the same momentum
      vc = (p1_norm - p2_norm)
      IF (DOT_PRODUCT(vc, vc) < eps) CYCLE

      p1_norm = p1 / m1
      gm1 = SQRT(DOT_PRODUCT(p1_norm, p1_norm) + 1.0_num) * m1

      p2_norm = p2 / m2
      gm2 = SQRT(DOT_PRODUCT(p2_norm, p2_norm) + 1.0_num) * m2

      gm = gm1 + gm2

      ! Pre-collision velocities
      v1 = p1 / gm1
      v2 = p2 / gm2

      ! Velocity of centre-of-momentum (COM) reference frame
      vc = (p1 + p2) / gm
      vc_sq = DOT_PRODUCT(vc, vc)

      gamma_rel_inv = SQRT(1.0_num - vc_sq)
      gc = 1.0_num / gamma_rel_inv

      gc_m1_vc = (gc - 1.0_num) / vc_sq

      p3 = p1 + (gc_m1_vc * DOT_PRODUCT(vc, v1) - gc) * gm1 * vc

      v_sq = DOT_PRODUCT(vc, v1)
      gm3 = (1.0_num - v_sq) * gc * gm1
      v_sq = DOT_PRODUCT(vc, v2)
      gm4 = (1.0_num - v_sq) * gc * gm2

      p_mag2 = DOT_PRODUCT(p3, p3)
      p_mag = SQRT(p_mag2)

      fac = (q1 * q2)**2 * s_fac / (gm1 * gm2)
      s12 = fac * gc * p_mag * c / gm * (gm3 * gm4 / p_mag2 + 1.0_num)**2

      ! Cold plasma upper limit for s12
      v_rel = gm * p_mag * c / (gm3 * gm4 * gc)
      s_prime = s_fac_prime * (m1 + m2) * v_rel / MAX(m1, m2)

      s12 = MIN(s12, s_prime)

      ran1 = random()
      ran2 = random() * 2.0_num * pi

      ! Inversion from Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012)
      IF (s12 < 0.1_num) THEN
        cosp = 1.0_num + s12 * LOG(MAX(ran1, 5e-9_num))
      ELSE IF (s12 >= 0.1_num .AND. s12 < 3.0_num) THEN
        a_inv = 0.0056958_num + (0.9560202_num + (-0.508139_num &
             + (0.47913906_num + (-0.12788975_num + 0.02389567_num &
             * s12) * s12) * s12) * s12) * s12
        a = 1.0_num / a_inv
        cosp = a_inv * LOG(EXP(-a) + 2.0_num * ran1 * SINH(a))
      ELSE IF (s12 >= 3.0_num .AND. s12 < 6.0_num) THEN
        a = 3.0_num * EXP(-s12)
        cosp = LOG(EXP(-a) + 2.0_num * ran1 * SINH(a)) / a
      ELSE
        cosp = 2.0_num * ran1 - 1.0_num
      END IF

      ! Branches 2 and 3 can result in rounding errors
      cosp = MAX(MIN(cosp, 1.0_num), -1.0_num)

      sinp = SIN(ACOS(cosp))

      ! Calculate new momenta according to rotation by angle p
      p_perp2 = p3(1)**2 + p3(2)**2
      p_perp = SQRT(p_perp2)
      p_tot = SQRT(p_perp2 + p3(3)**2)
      p_perp_inv = 1.0_num / (p_perp + c_tiny)

      mat(1,1) =  p3(1) * p3(3) * p_perp_inv
      mat(1,2) = -p3(2) * p_tot * p_perp_inv
      mat(1,3) =  p3(1)
      mat(2,1) =  p3(2) * p3(3) * p_perp_inv
      mat(2,2) =  p3(1) * p_tot * p_perp_inv
      mat(2,3) =  p3(2)
      mat(3,1) = -p_perp
      mat(3,2) =  0.0_num
      mat(3,3) =  p3(3)

      sinp_cos = sinp * COS(ran2)
      sinp_sin = sinp * SIN(ran2)

      p3(1) = mat(1,1) * sinp_cos + mat(1,2) * sinp_sin + mat(1,3) * cosp
      p3(2) = mat(2,1) * sinp_cos + mat(2,2) * sinp_sin + mat(2,3) * cosp
      p3(3) = mat(3,1) * sinp_cos + mat(3,2) * sinp_sin + mat(3,3) * cosp

      p4 = -p3

      p5 = (p3 + (gc_m1_vc * DOT_PRODUCT(vc, p3) + gm3 * gc) * vc) * c
      p6 = (p4 + (gc_m1_vc * DOT_PRODUCT(vc, p4) + gm4 * gc) * vc) * c

      ! Update particle properties
      current%part_p = p5
      impact%part_p = p6

      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO

    ! restore the tail of the list
    NULLIFY(p_list%tail%next)

  END SUBROUTINE intra_collisions_np



  SUBROUTINE inter_collisions_sk(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, log_lambda, user_factor )

    TYPE(particle_list), INTENT(INOUT) :: p_list1
    TYPE(particle_list), INTENT(INOUT) :: p_list2

    REAL(num), INTENT(IN) :: mass1, charge1, weight1
    REAL(num), INTENT(IN) :: mass2, charge2, weight2

    REAL(num), INTENT(IN) :: idens, jdens
    REAL(num), INTENT(IN) :: log_lambda
    REAL(num), INTENT(IN) :: user_factor

    TYPE(particle), POINTER :: current, impact

    REAL(num) :: factor, np
    INTEGER(i8) :: icount, jcount, pcount, k

    REAL(num) :: m1, m2, q1, q2, w1, w2
    REAL(num), DIMENSION(3) :: p1, p2 ! Pre-collision momenta
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm ! Normalised momenta
    REAL(num), DIMENSION(3) :: vc ! Velocity of COM frame wrt lab frame
    REAL(num), DIMENSION(3) :: p3, p4
    REAL(num), DIMENSION(3) :: p5, p6
    REAL(num), DIMENSION(3) :: v3, v4
    REAL(num), DIMENSION(3) :: vr, vcr
    REAL(num), DIMENSION(3) :: c1, c2, c3
    REAL(num) :: e1, e2 ! Pre-collision energies
    REAL(num) :: e3, e4, e5, e6
    REAL(num) :: gamma_rel, gamma_rel2, gamma_rel_m1, gamma_rel_r
    REAL(num) :: tvar ! Dummy variable for temporarily storing values
    REAL(num) :: vc_sq, vc_sq_cc, p1_vc, p2_vc, p3_mag
    REAL(num) :: delta, sin_theta, cos_theta, tan_theta_cm, tan_theta_cm2
    REAL(num) :: vrabs, denominator, wr
    REAL(num) :: nu, ran1, ran2
    factor = 0.0_num
    np = 0.0_num

    ! Inter-species collisions
    icount = p_list1%count
    jcount = p_list2%count
    pcount = MAX(icount, jcount)

    IF (icount > 0 .AND. jcount > 0) THEN
      ! temporarily join tail to the head of the lists to make them circular
      p_list1%tail%next => p_list1%head
      p_list2%tail%next => p_list2%head

#ifdef PER_SPECIES_WEIGHT
      np = icount * weight1
      factor = pcount * MIN(weight1, weight2)
#else
      current => p_list1%head
      impact => p_list2%head

      IF (icount >= jcount) THEN
        DO k = 1, icount
          np = np + current%weight
          current => current%next
        END DO
      ELSE
        DO k = 1, jcount
          np = np + impact%weight
          impact => impact%next
        END DO
      END IF

      current => p_list1%head
      impact => p_list2%head

      DO k = 1, pcount
        factor = factor + MIN(current%weight, impact%weight)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO
#endif
      factor = user_factor / factor

      ! If possible, use per-species properties
      m1 = mass1
      m2 = mass2
      q1 = charge1
      q2 = charge2
      w1 = weight1
      w2 = weight2
      wr = w1 / w2

      current => p_list1%head
      impact => p_list2%head
      DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
        m1 = current%mass
        m2 = impact%mass
        q1 = current%charge
        q2 = impact%charge
#endif
#ifndef PER_SPECIES_WEIGHT
        w1 = current%weight
        w2 = impact%weight
        wr = w1 / w2
#endif
        ! Copy all of the necessary particle data into variables with easier to
        ! read names
        p1 = current%part_p
        p2 = impact%part_p
        p1_norm = p1 / mc0
        p2_norm = p2 / mc0

        ! Two stationary particles can't collide, so don't try
        IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
            .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE

        ! Ditto for two particles with the same momentum
        vc = (p1_norm - p2_norm)
        IF (DOT_PRODUCT(vc, vc) < eps) CYCLE

        ! Pre-collision energies
        e1 = c * SQRT(DOT_PRODUCT(p1, p1) + (m1 * c)**2)
        e2 = c * SQRT(DOT_PRODUCT(p2, p2) + (m2 * c)**2)

        ! Velocity of centre-of-momentum (COM) reference frame
        vc = (p1 + p2) * cc / (e1 + e2)
        vc_sq = DOT_PRODUCT(vc, vc)
        vc_sq_cc = vc_sq / cc

        gamma_rel2 = 1.0_num / (1.0_num - vc_sq_cc)
        gamma_rel = SQRT(gamma_rel2)
        gamma_rel_m1 = gamma_rel2 * vc_sq_cc / (gamma_rel + 1.0_num)

        ! Lorentz momentum transform to get into COM frame
        p1_vc = DOT_PRODUCT(p1, vc)
        p2_vc = DOT_PRODUCT(p2, vc)
        tvar = p1_vc * gamma_rel_m1 / (vc_sq + c_tiny)
        p3 = p1 + vc * (tvar - gamma_rel * e1 / cc)
        tvar = p2_vc * gamma_rel_m1 / (vc_sq + c_tiny)
        p4 = p2 + vc * (tvar - gamma_rel * e2 / cc)

        p3_mag = SQRT(DOT_PRODUCT(p3, p3))

        ! Lorentz energy transform
        e3 = gamma_rel * (e1 - p1_vc)
        e4 = gamma_rel * (e2 - p2_vc)
        ! Pre-collision velocities in COM frame
        v3 = p3 * cc / e3
        v4 = p4 * cc / e4

        ! Relative velocity
        tvar = 1.0_num - (DOT_PRODUCT(v3, v4) / cc)
        vr = (v3 - v4) / tvar
        vrabs = SQRT(DOT_PRODUCT(vr, vr))

        ! Collision frequency
        nu = coll_freq(vrabs, log_lambda, m1, m2, q1, q2, MIN(idens, jdens))
        nu = MIN(nu * factor * np * dt_coll, 0.02_num)

        ! NOTE: nu is now the number of collisions per timestep, NOT collision
        ! frequency

        ! New coordinate system to simplify scattering.
        CALL new_coords(vr, c1, c2, c3)

        ! this is to ensure that 0 < ran1 < 1
        ! ran1=0 gives NaN in logarithm
        ! ran1=1 could give positive logarithm due to rounding errors
        ran1 = (1.0_num - 1.0e-10_num) * random() + 0.5e-10_num
        ran2 = 2.0_num * pi * random()

        ! Box Muller method for random from Gaussian distribution,
        ! mean 0, variance of nu
        ! Possible place for speed up by caching the second Box Muller number
        ! and using it later
        ! SQRT(-2.0_num * nu * LOG(ran1)) * COS(ran2)
        delta = SQRT(-2.0_num * nu * LOG(ran1)) * SIN(ran2)
        ran2 = 2.0_num * pi * random()

        ! angle theta in the One Particle at Rest frame
        sin_theta = 2.0_num * delta / (1.0_num + delta**2)
        cos_theta = (1.0_num - delta**2) / (1.0_num + delta**2)

        ! Transform angles from particle j's rest frame to COM frame
        ! Note azimuthal angle (ran2) is invariant under this transformation
        IF (m1 > m2) THEN
          vcr = v3
        ELSE
          vcr = v4
        END IF
        gamma_rel_r = 1.0_num / SQRT(1.0_num - (DOT_PRODUCT(vcr, vcr) / cc))

        denominator = gamma_rel_r * (cos_theta - SQRT(DOT_PRODUCT(vcr, vcr)) &
            / MAX(vrabs, c_tiny))
        IF (ABS(denominator) > SQRT(c_tiny)) THEN
          tan_theta_cm = sin_theta / denominator
          tan_theta_cm2 = tan_theta_cm**2
        ELSE
          tan_theta_cm = c_largest_number
          tan_theta_cm2 = c_largest_number
        END IF

        sin_theta = tan_theta_cm / SQRT(1.0_num + tan_theta_cm2)
        cos_theta = 1.0_num / SQRT(1.0_num + tan_theta_cm2)

        ! Post-collision momenta in COM frame
        p3 = p3_mag * (c1 * cos_theta + c2 * sin_theta * COS(ran2) &
            + c3 * sin_theta * SIN(ran2))
        p4 = -p3

        ! Lorentz momentum transform to get back to lab frame
        tvar = DOT_PRODUCT(p3, vc) * gamma_rel_m1 / vc_sq
        p5 = p3 + vc * (tvar + gamma_rel * e3 / cc)
        tvar = DOT_PRODUCT(p4, vc) * gamma_rel_m1 / vc_sq
        p6 = p4 + vc * (tvar + gamma_rel * e4 / cc)

        e5 = c * SQRT(DOT_PRODUCT(p5, p5) + (m1 * c)**2)
        e6 = c * SQRT(DOT_PRODUCT(p6, p6) + (m2 * c)**2)

        IF (wr > one_p_2eps) THEN
          CALL weighted_particles_correction(w2 / w1, p1, p5, e1, e5, m1)
        ELSE IF (wr < one_m_2eps) THEN
          CALL weighted_particles_correction(w1 / w2, p2, p6, e2, e6, m2)
        END IF

        ! Update particle properties
        current%part_p = p5
        impact%part_p = p6

        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO

      ! restore the tail of the lists
      NULLIFY(p_list1%tail%next)
      NULLIFY(p_list2%tail%next)
    END IF

  END SUBROUTINE inter_collisions_sk



  ! Binary collision scattering operator based jointly on:
  ! Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012), and
  ! K. Nanbu and S. Yonemura, J. Comput. Phys. 145, 639 (1998)
  SUBROUTINE inter_collisions_np(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, log_lambda, user_factor )

    TYPE(particle_list), INTENT(INOUT) :: p_list1
    TYPE(particle_list), INTENT(INOUT) :: p_list2
    REAL(num), INTENT(IN) :: mass1, charge1, weight1
    REAL(num), INTENT(IN) :: mass2, charge2, weight2
    REAL(num), INTENT(IN) :: idens, jdens
    REAL(num), INTENT(IN) :: log_lambda
    REAL(num), INTENT(IN) :: user_factor
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor
    INTEGER(i8) :: icount, jcount, pcount, k
    REAL(num) :: m1, m2, q1, q2, w1, w2
    REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
    REAL(num) :: sinp_cos, sinp_sin, s_prime, s_fac_prime
    REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
    REAL(num) :: p_perp2, p_perp_inv, cell_fac
    REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm
    REAL(num), DIMENSION(3,3) :: mat
    REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq
    REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    REAL(num), PARAMETER :: pi_fac = &
                                (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)

    factor = 0.0_num

    ! Inter-species collisions
    icount = p_list1%count
    jcount = p_list2%count
    pcount = MAX(icount, jcount)

    IF (icount > 0 .AND. jcount > 0) THEN
      ! temporarily join tail to the head of the lists to make them circular
      p_list1%tail%next => p_list1%head
      p_list2%tail%next => p_list2%head

#ifdef PER_SPECIES_WEIGHT
      factor = pcount * MIN(weight1, weight2)
#else
      current => p_list1%head
      impact => p_list2%head

      DO k = 1, pcount
        factor = factor + MIN(current%weight, impact%weight)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO
#endif
      factor = user_factor / factor

      ! If possible, use per-species properties
      m1 = mass1
      m2 = mass2
      q1 = charge1
      q2 = charge2
      w1 = weight1
      w2 = weight2

      current => p_list1%head
      impact => p_list2%head

      ! Per-cell constant factors
      cell_fac = idens * jdens * dt_coll * factor * dx * dy
      s_fac = cell_fac * log_lambda / pi4_eps2_c4
      s_fac_prime = cell_fac * pi_fac

      DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
        m1 = current%mass
        m2 = impact%mass
        q1 = current%charge
        q2 = impact%charge
#endif
#ifndef PER_SPECIES_WEIGHT
        w1 = current%weight
        w2 = impact%weight
#endif

        p1 = current%part_p / c
        p2 = impact%part_p / c

        p1_norm = p1 / m0
        p2_norm = p2 / m0

        ! Two stationary particles can't collide, so don't try
        IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
            .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE

        ! Ditto for two particles with the same momentum
        vc = (p1_norm - p2_norm)
        IF (DOT_PRODUCT(vc, vc) < eps) CYCLE

        p1_norm = p1 / m1
        gm1 = SQRT(DOT_PRODUCT(p1_norm, p1_norm) + 1.0_num) * m1

        p2_norm = p2 / m2
        gm2 = SQRT(DOT_PRODUCT(p2_norm, p2_norm) + 1.0_num) * m2

        gm = gm1 + gm2

        ! Pre-collision velocities
        v1 = p1 / gm1
        v2 = p2 / gm2

        ! Velocity of centre-of-momentum (COM) reference frame
        vc = (p1 + p2) / gm
        vc_sq = DOT_PRODUCT(vc, vc)

        gamma_rel_inv = SQRT(1.0_num - vc_sq)
        gc = 1.0_num / gamma_rel_inv

        gc_m1_vc = (gc - 1.0_num) / vc_sq

        p3 = p1 + (gc_m1_vc * DOT_PRODUCT(vc, v1) - gc) * gm1 * vc

        v_sq = DOT_PRODUCT(vc, v1)
        gm3 = (1.0_num - v_sq) * gc * gm1
        v_sq = DOT_PRODUCT(vc, v2)
        gm4 = (1.0_num - v_sq) * gc * gm2

        p_mag2 = DOT_PRODUCT(p3, p3)
        p_mag = SQRT(p_mag2)

        fac = (q1 * q2)**2 * s_fac / (gm1 * gm2)
        s12 = fac * gc * p_mag * c / gm * (gm3 * gm4 / p_mag2 + 1.0_num)**2

        ! Cold plasma upper limit for s12
        v_rel = gm * p_mag * c / (gm3 * gm4 * gc)
        s_prime = s_fac_prime * (m1 + m2) * v_rel &
            / MAX(m1 * idens**two_thirds, m2 * jdens**two_thirds)

        s12 = MIN(s12, s_prime)

        ran1 = random()
        ran2 = random() * 2.0_num * pi

        ! Inversion from Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012)
        IF (s12 < 0.1_num) THEN
          cosp = 1.0_num + s12 * LOG(MAX(ran1, 5e-9_num))
        ELSE IF (s12 >= 0.1_num .AND. s12 < 3.0_num) THEN
          a_inv = 0.0056958_num + (0.9560202_num + (-0.508139_num &
               + (0.47913906_num + (-0.12788975_num + 0.02389567_num &
               * s12) * s12) * s12) * s12) * s12
          a = 1.0_num / a_inv
          cosp = a_inv * LOG(EXP(-a) + 2.0_num * ran1 * SINH(a))
        ELSE IF (s12 >= 3.0_num .AND. s12 < 6.0_num) THEN
          a = 3.0_num * EXP(-s12)
          cosp = LOG(EXP(-a) + 2.0_num * ran1 * SINH(a)) / a
        ELSE
          cosp = 2.0_num * ran1 - 1.0_num
        END IF

        ! Branches 2 and 3 can result in rounding errors
        cosp = MAX(MIN(cosp, 1.0_num), -1.0_num)

        sinp = SIN(ACOS(cosp))

        ! Calculate new momenta according to rotation by angle p
        p_perp2 = p3(1)**2 + p3(2)**2
        p_perp = SQRT(p_perp2)
        p_tot = SQRT(p_perp2 + p3(3)**2)
        p_perp_inv = 1.0_num / (p_perp + c_tiny)

        mat(1,1) =  p3(1) * p3(3) * p_perp_inv
        mat(1,2) = -p3(2) * p_tot * p_perp_inv
        mat(1,3) =  p3(1)
        mat(2,1) =  p3(2) * p3(3) * p_perp_inv
        mat(2,2) =  p3(1) * p_tot * p_perp_inv
        mat(2,3) =  p3(2)
        mat(3,1) = -p_perp
        mat(3,2) =  0.0_num
        mat(3,3) =  p3(3)

        sinp_cos = sinp * COS(ran2)
        sinp_sin = sinp * SIN(ran2)

        p3(1) = mat(1,1) * sinp_cos + mat(1,2) * sinp_sin + mat(1,3) * cosp
        p3(2) = mat(2,1) * sinp_cos + mat(2,2) * sinp_sin + mat(2,3) * cosp
        p3(3) = mat(3,1) * sinp_cos + mat(3,2) * sinp_sin + mat(3,3) * cosp

        p4 = -p3

        ran1 = random()
        IF (ran1 < w2 / w1) THEN
          p5 = (p3 + (gc_m1_vc * DOT_PRODUCT(vc, p3) + gm3 * gc) * vc) * c
          current%part_p = p5
        END IF
        IF (ran1 < w1 / w2) THEN
          p6 = (p4 + (gc_m1_vc * DOT_PRODUCT(vc, p4) + gm4 * gc) * vc) * c
          ! Update particle properties
          impact%part_p = p6
        END IF

        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO

      ! restore the tail of the lists
      NULLIFY(p_list1%tail%next)
      NULLIFY(p_list2%tail%next)
    END IF

  END SUBROUTINE inter_collisions_np



  PURE FUNCTION coll_freq(vrabs, log_lambda, m1, m2, q1, q2, jdens)

    REAL(num), INTENT(IN) :: vrabs, log_lambda, m1, m2, q1, q2
    REAL(num), INTENT(IN) :: jdens
    REAL(num) :: mu, coll_freq, numerator, denominator
    REAL(num), PARAMETER :: fac = 4.0_num * pi * epsilon0**2

    mu = (m1 * m2) / (m1 + m2)

    IF (vrabs > 0.0_num) THEN
      numerator = (q1 * q2)**2 * jdens * log_lambda
      denominator = fac * mu**2 * vrabs**3
      IF (denominator <= 0.0_num &
          .OR. EXPONENT(numerator) - EXPONENT(denominator) &
          >= c_maxexponent) THEN
        coll_freq = 0.0_num
      ELSE
        coll_freq = numerator / denominator
      END IF
    ELSE
      coll_freq = 0.0_num
    END IF

  END FUNCTION coll_freq



  SUBROUTINE weighted_particles_correction(wtr, p, p_scat, en, en_scat, mass)

    ! This is the correction to the particle according to
    ! Sentoku and Kemp (2008) formulas 21 to 26.
    REAL(num), INTENT(INOUT) :: p_scat(3)
    REAL(num), INTENT(IN) :: p(3)
    REAL(num), INTENT(IN) :: wtr
    REAL(num), INTENT(IN) :: en, en_scat
    REAL(num), INTENT(IN) :: mass

    REAL(num) :: p_after(3)
    REAL(num) :: en_after
    REAL(num) :: gamma_en, gamma_p
    REAL(num) :: delta_p, p_mag, p_trans_mag
    REAL(num) :: c1(3), c2(3), c3(3)
    REAL(num) :: phi

    en_after = (1.0_num - wtr) * en + wtr * en_scat
    p_after  = (1.0_num - wtr) * p  + wtr * p_scat
    p_mag = SQRT(DOT_PRODUCT(p_after, p_after))
    gamma_en = en_after / (mass * cc)
    gamma_p = SQRT(1.0_num + (p_mag / mass / c)**2)

    ! This if-statement is just to take care of possible rounding errors
    ! gamma_p should always be smaller than gamma_en
    IF (gamma_p < gamma_en) THEN
      ! magnitude of the momentum correction
      delta_p = mass * c * SQRT(gamma_en**2 - gamma_p**2)
      p_trans_mag = SQRT(p_after(2)**2 + p_after(3)**2)

      CALL new_coords(p_after, c1, c2, c3)

      phi = 2.0_num * pi * random()

      ! Correcting for the loss in energy by adding a perpendicular
      ! momentum correction
      p_scat = p_after + delta_p * (c2 * COS(phi) + c3 * SIN(phi))
    END IF

  END SUBROUTINE weighted_particles_correction



  SUBROUTINE new_coords(vector, c1, c2, c3)

    ! New orthonormal basis vectors for calculating scattering angles:
    ! c1 = v / |v|
    ! c2 = (v x e1) / |v x e1|
    ! c3 = ((v x e1) x v) / |(v x e1) x v|
    ! where e1 is a unit vector parallel to x-axis. I.e., e1 = (1,0,0)
    ! New x-axis, c1, is parallel to the input vector
    ! New y-axis, c2, is perpendicular to c1 and the x-axis
    ! New z-axis ,c3, is perpendicular to c1 and c2
    REAL(num), DIMENSION(3), INTENT(IN) :: vector
    REAL(num), DIMENSION(3), INTENT(OUT) :: c1, c2, c3
    REAL(num) :: vtrans, vmag

    vmag = SQRT(DOT_PRODUCT(vector, vector))
    vtrans = SQRT(vector(2)**2 + vector(3)**2)

    IF (vtrans > c_tiny) THEN
      c1 = vector / vmag
      c2 = (/ 0.0_num, vector(3), -vector(2) /)
      c2 = c2 / vtrans
      c3 = (/ vtrans**2, &
          -(vector(1) * vector(2)), &
          -(vector(1) * vector(3)) /)
      c3 = c3 / (vmag * vtrans)
    ELSE
      c1 = (/ 1.0_num, 0.0_num, 0.0_num /)
      c2 = (/ 0.0_num, 1.0_num, 0.0_num /)
      c3 = (/ 0.0_num, 0.0_num, 1.0_num /)
    END IF

  END SUBROUTINE new_coords



  SUBROUTINE shuffle_particle_list_random(p_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list
    TYPE(particle), POINTER :: particle1, particle2

    INTEGER :: i, idx, swap_idx
    INTEGER :: p_num

    p_num = INT(p_list%count)

    ! Nothing to be done
    IF (p_num <= 2) RETURN

    ! First make sure that the sorting array is large enough
    ! This should happen infrequently
    IF (p_num > coll_sort_array_size) THEN
      DEALLOCATE(coll_sort_array)

      ! make the sort array somewhat larger to avoid frequent deallocation
      ! and reallocation
      coll_sort_array_size = (11 * p_num) / 10 + 10
      ALLOCATE(coll_sort_array(coll_sort_array_size))
    END IF

    ! Copy all the particle pointers into the array and create random
    ! sort indices
    particle1 => p_list%head
    DO i = 1,p_num
      coll_sort_array(i)%particle => particle1
      particle1 => particle1%next
    END DO

    ! Shuffle particles using Durstenfeld's algorithm
    DO idx = p_num,2,-1
      swap_idx = FLOOR(idx * random()) + 1
      particle1 => coll_sort_array(idx)%particle
      coll_sort_array(idx)%particle => coll_sort_array(swap_idx)%particle
      coll_sort_array(swap_idx)%particle => particle1
    END DO

    ! Finally we have to copy back to the list
    ! Do head first
    particle1 => coll_sort_array(1)%particle
    p_list%head => particle1
    NULLIFY(particle1%prev)

    ! Then do all the particle links between head and tail.
    DO i = 2, p_num
      particle2 => particle1
      particle1 => coll_sort_array(i)%particle

      particle1%prev => particle2
      particle2%next => particle1
    END DO

    ! Finally set the tail (at the end of the loop, particle is pointing to
    ! the tail)
    p_list%tail => particle1
    NULLIFY(particle1%next)

  END SUBROUTINE shuffle_particle_list_random



  PURE FUNCTION calc_coulomb_log(ekbar1, temp2, dens1, dens2, q1, q2, m1)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: ekbar1, temp2
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: dens1, dens2
    REAL(num), INTENT(IN) :: q1, q2, m1
    REAL(num), DIMENSION(1-ng:nx+ng,1-ng:ny+ng) :: calc_coulomb_log
    REAL(num) :: b0, dB, bmin, bmax
    REAL(num) :: local_ekbar1, local_temp2, gamm
    INTEGER :: i, j

    calc_coulomb_log = 0.0_num
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng
      local_ekbar1 = MAX(ekbar1(i,j), 100.0_num * q0)
      local_temp2 = MAX(temp2(i,j), 100.0_num)
      IF (dens1(i,j) <= 1.0_num .OR. dens2(i,j) <= 1.0_num) THEN
        calc_coulomb_log(i,j) = 1.0_num
      ELSE
        bmax = SQRT(epsilon0 * q0 * local_temp2 / (ABS(q2) * q0 * dens2(i,j)))
        b0 = ABS(q1 * q2) / (8.0_num * pi * epsilon0 * local_ekbar1)
        gamm = (local_ekbar1 / (m1 * cc)) + 1.0_num
        dB = 2.0_num * pi * h_bar / (SQRT(gamm**2 - 1.0_num) * m1 * c)
        bmin = MAX(b0, dB)
        calc_coulomb_log(i,j) = MAX(1.0_num, LOG(bmax / bmin))
      END IF
    END DO
    END DO

  END FUNCTION calc_coulomb_log



  SUBROUTINE calc_coll_number_density(data_array, ispecies)

    ! This routine calculates an approximate number density
    ! It assumes each particle only contributes to a single cell
    ! Returns zero if species_type == c_species_id_photon

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: jx, jy, cell_x, cell_y
    TYPE(particle), POINTER :: current

    data_array = 0.0_num

    IF (species_list(ispecies)%species_type == c_species_id_photon) RETURN

    idx = 1.0_num / dx / dy

    wdata = species_list(ispecies)%weight

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)

        data_array(cell_x, cell_y) = &
            data_array(cell_x, cell_y) + wdata

        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    data_array = data_array * idx

  END SUBROUTINE calc_coll_number_density



  SUBROUTINE calc_coll_temperature_ev(sigma, ispecies)

    ! This subroutine calculates the grid-based temperature of a given
    ! particle species.
    ! It is almost identical to the calc_temperature subroutine in calc_df,
    ! except it uses the secondary_list rather than the attached_list.

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    INTEGER :: ix, iy
    INTEGER :: jx, jy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    sqrt_part_m  = SQRT(species_list(ispecies)%mass)
    part_w = species_list(ispecies)%weight

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

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
        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    sqrt_part_m  = SQRT(species_list(ispecies)%mass)
    part_w = species_list(ispecies)%weight

    part_count = 0.0_num
    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

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

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * part_w
          sigma(cell_x+ix, cell_y+iy) = &
              sigma(cell_x+ix, cell_y+iy) + gf &
              * ((part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
              + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
              + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2)
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        END DO
        END DO
        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / q0 / 3.0_num

  END SUBROUTINE calc_coll_temperature_ev



  SUBROUTINE calc_coll_ekbar(data_array, ispecies)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1
    INTEGER :: ix, iy
    INTEGER :: jx, jy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    data_array = 0.0_num
    part_count = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    part_mc = c * species_list(ispecies)%mass
    part_w = species_list(ispecies)%weight
    fac = part_mc * part_w * c

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

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
          gf = gx(ix) * gy(iy)
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * wdata
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf * part_w
        END DO
        END DO

        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    data_array = data_array / MAX(part_count, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_coll_ekbar



  SUBROUTINE setup_collisions

    ALLOCATE(coll_pairs(n_species, n_species))
    ALLOCATE(coll_pairs_state(n_species, n_species))
    coll_pairs = 1.0_num
    coll_pairs_state = c_coll_collide
    coll_sort_array_size = 1
    ALLOCATE(coll_sort_array(coll_sort_array_size))

    IF (use_nanbu) THEN
      intra_coll_fn => intra_collisions_np
      inter_coll_fn => inter_collisions_np
    ELSE
      intra_coll_fn => intra_collisions_sk
      inter_coll_fn => inter_collisions_sk
    END IF

  END SUBROUTINE setup_collisions



  SUBROUTINE deallocate_collisions

    INTEGER :: stat

    DEALLOCATE(coll_pairs, coll_pairs_state, coll_sort_array, STAT=stat)

  END SUBROUTINE deallocate_collisions

END MODULE collisions
