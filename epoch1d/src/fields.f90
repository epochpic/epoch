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

MODULE fields

  USE boundary

  IMPLICIT NONE

  INTEGER :: field_order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx
  REAL(num) :: cnx
  REAL(num) :: alphax = 1.0_num
  REAL(num) :: deltax = 0.0_num

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

    field_order = order
    fng = field_order / 2

    IF (field_order == 2) THEN
      cfl = 1.0_num
    ELSE IF (field_order == 4) THEN
      cfl = 6.0_num / 7.0_num
    ELSE
      cfl = 120.0_num / 149.0_num
    END IF

  END SUBROUTINE set_field_order



  SUBROUTINE set_maxwell_solver

    REAL(num) :: dx_cdt

    IF (maxwell_solver == c_maxwell_solver_custom) THEN
      alphax = 1.0_num - 3.0_num * deltax

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_x) THEN
      ! R. Lehe et al., Phys. Rev. ST Accel. Beams 16, 021301 (2013)
      dx_cdt = dx / (c * dt)
      deltax = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      alphax = 1.0_num - 3.0_num * deltax
    END IF

    IF (rank == 0 .AND. maxwell_solver /= c_maxwell_solver_yee) THEN
      PRINT*, 'Maxwell solver set to the following parameters:'
      PRINT'(A9, 1F14.9)', 'alpha =', alphax
      PRINT'(A9, 1F14.9)', 'delta =', deltax
      PRINT'(A9, 1F14.9)', 'c*dt/dx = ', dt * c / dx
      PRINT*
    END IF

  END SUBROUTINE set_maxwell_solver



  SUBROUTINE update_e_field

    INTEGER :: ix
    REAL(num) :: cpml_x
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        DO ix = 0, nx
          cx1 = cnx / cpml_kappa_ex(ix)

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              - fac * jz(ix)
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 0, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              - fac * jz(ix)
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 0, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x
          cx3 = c3 * cpml_x

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - cx3 * (bz(ix+2) - bz(ix-3)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              + cx3 * (by(ix+2) - by(ix-3)) &
              - fac * jz(ix)
        END DO
      END IF

      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = cnx

        DO ix = 0, nx
          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              - fac * jz(ix)
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx

        DO ix = 0, nx
          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              - fac * jz(ix)
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx
        cx3 = c3 * cnx

        DO ix = 0, nx
          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - cx3 * (bz(ix+2) - bz(ix-3)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              + cx3 * (by(ix+2) - by(ix-3)) &
              - fac * jz(ix)
        END DO
      END IF
    END IF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix
    REAL(num) :: cpml_x
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO ix = 0, nx
            cx1 = hdtx / cpml_kappa_bx(ix)

            by(ix) = by(ix) &
                + cx1 * (ez(ix+1) - ez(ix  ))

            bz(ix) = bz(ix) &
                - cx1 * (ey(ix+1) - ey(ix  ))
          END DO
        ELSE
          DO ix = 0, nx
            cx1 = hdtx / cpml_kappa_bx(ix)

            by(ix) = by(ix) &
                + cx1 * (alphax * (ez(ix+1) - ez(ix  ))  &
                       + deltax * (ez(ix+2) - ez(ix-1)))

            bz(ix) = bz(ix) &
                - cx1 * (alphax * (ey(ix+1) - ey(ix  ))  &
                       + deltax * (ey(ix+2) - ey(ix-1)))
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 0, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1))
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 0, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x
          cx3 = c3 * cpml_x

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1)) &
              + cx3 * (ez(ix+3) - ez(ix-2))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1)) &
              - cx3 * (ey(ix+3) - ey(ix-2))
        END DO
      END IF

      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = hdtx

        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO ix = 0, nx
            by(ix) = by(ix) &
                + cx1 * (ez(ix+1) - ez(ix  ))

            bz(ix) = bz(ix) &
                - cx1 * (ey(ix+1) - ey(ix  ))
          END DO
        ELSE
          DO ix = 0, nx
            by(ix) = by(ix) &
                + cx1 * (alphax * (ez(ix+1) - ez(ix  ))  &
                       + deltax * (ez(ix+2) - ez(ix-1)))

            bz(ix) = bz(ix) &
                - cx1 * (alphax * (ey(ix+1) - ey(ix  ))  &
                       + deltax * (ey(ix+2) - ey(ix-1)))
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx

        DO ix = 0, nx
          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1))
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx
        cx3 = c3 * hdtx

        DO ix = 0, nx
          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1)) &
              + cx3 * (ez(ix+3) - ez(ix-2))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1)) &
              - cx3 * (ey(ix+3) - ey(ix-2))
        END DO
      END IF
    END IF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx

    cnx = hdtx * c**2

    fac = hdt / epsilon0

    ! Update E field to t+dt/2
    CALL update_e_field

    ! Now have E(t+dt/2), do boundary conditions on E
    CALL efield_bcs

    ! Energy lost through boundaries using E(t+dt/2) and B(t)
    IF (dump_field_energy_lost) CALL calc_field_energy_lost

    ! Update B field to t+dt/2 using E(t+dt/2)
    CALL update_b_field

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    hdt  = 0.5_num * dt
    hdtx = hdt / dx

    cnx = hdtx * c**2

    fac = hdt / epsilon0

    CALL update_b_field

    CALL bfield_final_bcs

    IF (dump_field_energy_lost) CALL calc_field_energy_lost

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final



  ! Returns the electromagnetic power leaving the domain through boundary
  FUNCTION boundary_field_power(bd) RESULT(power)

    INTEGER, INTENT(IN) :: bd
    REAL(num) :: power

    SELECT CASE(bd)
      CASE(c_bd_x_min)
        power = -(ey(1) * bz(0)  - ez(1) * by(0))
      CASE(c_bd_x_max)
        power =  ey(nx) * bz(nx) - ez(nx) * by(nx)
      CASE DEFAULT
        power = 0.0_num
    END SELECT

    ! Negative values represent field energy entering the simulation, since we
    ! only want to track energy lost, set these to 0.
    IF (power < 0.0_num) THEN
      power = 0.0_num
    END IF

    power = power / mu0

  END FUNCTION boundary_field_power



  ! Accumulates the net field energy leaving the domain through each
  ! simple_laser or simple_outflow boundary over half a timestep.
  SUBROUTINE calc_field_energy_lost

    INTEGER :: i

    DO i = 1, 2*c_ndims
      IF (.NOT.is_boundary(i)) CYCLE
      IF (bc_field(i) /= c_bc_simple_laser &
          .AND. bc_field(i) /= c_bc_simple_outflow) CYCLE
      field_energy_lost_local(i) = field_energy_lost_local(i) &
          + 0.5_num * dt * boundary_field_power(i)
    END DO

  END SUBROUTINE calc_field_energy_lost

END MODULE fields
