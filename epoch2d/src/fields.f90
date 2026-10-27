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
  REAL(num) :: hdtx, hdty
  REAL(num) :: cnx, cny
  REAL(num) :: alphax = 1.0_num, alphay = 1.0_num
  REAL(num) :: betaxy = 0.0_num, betayx = 0.0_num
  REAL(num) :: deltax = 0.0_num, deltay = 0.0_num

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

    REAL(num) :: delta, dx_cdt

    IF (maxwell_solver == c_maxwell_solver_custom) THEN
      alphax = 1.0_num - 2.0_num * betaxy - 3.0_num * deltax
      alphay = 1.0_num - 2.0_num * betayx - 3.0_num * deltay

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_x) THEN
      ! R. Lehe et al., Phys. Rev. ST Accel. Beams 16, 021301 (2013)
      dx_cdt = dx / (c * dt)
      betaxy = 0.125_num * (dx / dy)**2
      betayx = 0.125_num
      deltax = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      deltay = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy - 3.0_num * deltax
      alphay = 1.0_num - 2.0_num * betayx

    ELSE IF (maxwell_solver == c_maxwell_solver_lehe_y) THEN
      dx_cdt = dy / (c * dt)
      betayx = 0.125_num * (dy / dx)**2
      betaxy = 0.125_num
      deltax = 0.0_num
      deltay = 0.25_num * (1.0_num - dx_cdt**2 * SIN(0.5_num * pi / dx_cdt)**2)
      alphax = 1.0_num - 2.0_num * betaxy
      alphay = 1.0_num - 2.0_num * betayx - 3.0_num * deltay

    ELSE IF (maxwell_solver == c_maxwell_solver_pukhov) THEN
      ! A. Pukhov, Journal of Plasma Physics 61, 425-433 (1999)
      delta = MIN(dx, dy)

      betayx = 0.125_num * (delta / dx)**2
      betaxy = 0.125_num * (delta / dy)**2
      deltax = 0.0_num
      deltay = 0.0_num
      alphax = 1.0_num - 2.0_num * betaxy
      alphay = 1.0_num - 2.0_num * betayx
    END IF

    IF (rank == 0 .AND. maxwell_solver /= c_maxwell_solver_yee) THEN
      PRINT*, 'Maxwell solver set to the following parameters:'
      PRINT'(A9, 2F14.9)', 'alpha =', alphax, alphay
      PRINT'(A9, 1F14.9)', 'betax =', betaxy
      PRINT'(A9, 1F14.9)', 'betay =', betayx
      PRINT'(A9, 2F14.9)', 'delta =', deltax, deltay
      PRINT'(A9, 1F14.9)', 'c*dt/dx = ', dt * c / dx
      PRINT*
    END IF

  END SUBROUTINE set_maxwell_solver



  SUBROUTINE update_e_field

    INTEGER :: ix, iy
    REAL(num) :: cpml_x, cpml_y
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        DO iy = 0, ny
          cy1 = cny / cpml_kappa_ey(iy)
          DO ix = 0, nx
            cx1 = cnx / cpml_kappa_ex(ix)

            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iy = 0, ny
          cpml_y = cny / cpml_kappa_ey(iy)
          cy1 = c1 * cpml_y
          cy2 = c2 * cpml_y
          DO ix = 0, nx
            cpml_x = cnx / cpml_kappa_ex(ix)
            cx1 = c1 * cpml_x
            cx2 = c2 * cpml_x

            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                + cy2 * (bz(ix  , iy+1) - bz(ix  , iy-2)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - cx2 * (bz(ix+1, iy  ) - bz(ix-2, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                + cx2 * (by(ix+1, iy  ) - by(ix-2, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - cy2 * (bx(ix  , iy+1) - bx(ix  , iy-2)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iy = 0, ny
          cpml_y = cny / cpml_kappa_ey(iy)
          cy1 = c1 * cpml_y
          cy2 = c2 * cpml_y
          cy3 = c3 * cpml_y
          DO ix = 0, nx
            cpml_x = cnx / cpml_kappa_ex(ix)
            cx1 = c1 * cpml_x
            cx2 = c2 * cpml_x
            cx3 = c3 * cpml_x

            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                + cy2 * (bz(ix  , iy+1) - bz(ix  , iy-2)) &
                + cy3 * (bz(ix  , iy+2) - bz(ix  , iy-3)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - cx2 * (bz(ix+1, iy  ) - bz(ix-2, iy  )) &
                - cx3 * (bz(ix+2, iy  ) - bz(ix-3, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                + cx2 * (by(ix+1, iy  ) - by(ix-2, iy  )) &
                + cx3 * (by(ix+2, iy  ) - by(ix-3, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - cy2 * (bx(ix  , iy+1) - bx(ix  , iy-2)) &
                - cy3 * (bx(ix  , iy+2) - bx(ix  , iy-3)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      END IF

      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = cnx
        cy1 = cny

        DO iy = 0, ny
          DO ix = 0, nx
            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx
        cy1 = c1 * cny
        cy2 = c2 * cny

        DO iy = 0, ny
          DO ix = 0, nx
            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                + cy2 * (bz(ix  , iy+1) - bz(ix  , iy-2)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - cx2 * (bz(ix+1, iy  ) - bz(ix-2, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                + cx2 * (by(ix+1, iy  ) - by(ix-2, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - cy2 * (bx(ix  , iy+1) - bx(ix  , iy-2)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * cnx
        cx2 = c2 * cnx
        cx3 = c3 * cnx
        cy1 = c1 * cny
        cy2 = c2 * cny
        cy3 = c3 * cny

        DO iy = 0, ny
          DO ix = 0, nx
            ex(ix, iy) = ex(ix, iy) &
                + cy1 * (bz(ix  , iy  ) - bz(ix  , iy-1)) &
                + cy2 * (bz(ix  , iy+1) - bz(ix  , iy-2)) &
                + cy3 * (bz(ix  , iy+2) - bz(ix  , iy-3)) &
                - fac * jx(ix, iy)

            ey(ix, iy) = ey(ix, iy) &
                - cx1 * (bz(ix  , iy  ) - bz(ix-1, iy  )) &
                - cx2 * (bz(ix+1, iy  ) - bz(ix-2, iy  )) &
                - cx3 * (bz(ix+2, iy  ) - bz(ix-3, iy  )) &
                - fac * jy(ix, iy)

            ez(ix, iy) = ez(ix, iy) &
                + cx1 * (by(ix  , iy  ) - by(ix-1, iy  )) &
                + cx2 * (by(ix+1, iy  ) - by(ix-2, iy  )) &
                + cx3 * (by(ix+2, iy  ) - by(ix-3, iy  )) &
                - cy1 * (bx(ix  , iy  ) - bx(ix  , iy-1)) &
                - cy2 * (bx(ix  , iy+1) - bx(ix  , iy-2)) &
                - cy3 * (bx(ix  , iy+2) - bx(ix  , iy-3)) &
                - fac * jz(ix, iy)
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy
    REAL(num) :: cpml_x, cpml_y
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3

    IF (cpml_boundaries) THEN
      IF (field_order == 2) THEN
        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO iy = 0, ny
            cy1 = hdty / cpml_kappa_by(iy)
            DO ix = 0, nx
              cx1 = hdtx / cpml_kappa_bx(ix)

              bx(ix, iy) = bx(ix, iy) &
                  - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  ))

              by(ix, iy) = by(ix, iy) &
                  + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  ))

              bz(ix, iy) = bz(ix, iy) &
                  - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                  + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  ))
            END DO
          END DO
        ELSE
          DO iy = 0, ny
            cy1 = hdty / cpml_kappa_by(iy)
            DO ix = 0, nx
              cx1 = hdtx / cpml_kappa_bx(ix)

              bx(ix, iy) = bx(ix, iy) &
                  - cy1 * (alphay * (ez(ix  , iy+1) - ez(ix  , iy  ))  &
                         + betayx * (ez(ix+1, iy+1) - ez(ix+1, iy  )   &
                                   + ez(ix-1, iy+1) - ez(ix-1, iy  ))  &
                         + deltay * (ez(ix  , iy+2) - ez(ix  , iy-1)))

              by(ix, iy) = by(ix, iy) &
                  + cx1 * (alphax * (ez(ix+1, iy  ) - ez(ix  , iy  ))  &
                         + betaxy * (ez(ix+1, iy+1) - ez(ix  , iy+1)   &
                                   + ez(ix+1, iy-1) - ez(ix  , iy-1))  &
                         + deltax * (ez(ix+2, iy  ) - ez(ix-1, iy  )))

              bz(ix, iy) = bz(ix, iy) &
                  - cx1 * (alphax * (ey(ix+1, iy  ) - ey(ix  , iy  ))  &
                         + betaxy * (ey(ix+1, iy+1) - ey(ix  , iy+1)   &
                                   + ey(ix+1, iy-1) - ey(ix  , iy-1))  &
                         + deltax * (ey(ix+2, iy  ) - ey(ix-1, iy  ))) &
                  + cy1 * (alphay * (ex(ix  , iy+1) - ex(ix  , iy  ))  &
                         + betayx * (ex(ix+1, iy+1) - ex(ix+1, iy  )   &
                                   + ex(ix-1, iy+1) - ex(ix-1, iy  ))  &
                         + deltay * (ex(ix  , iy+2) - ex(ix  , iy-1)))
            END DO
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iy = 0, ny
          cpml_y = hdty / cpml_kappa_by(iy)
          cy1 = c1 * cpml_y
          cy2 = c2 * cpml_y
          DO ix = 0, nx
            cpml_x = hdtx / cpml_kappa_bx(ix)
            cx1 = c1 * cpml_x
            cx2 = c2 * cpml_x

            bx(ix, iy) = bx(ix, iy) &
                - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  )) &
                - cy2 * (ez(ix  , iy+2) - ez(ix  , iy-1))

            by(ix, iy) = by(ix, iy) &
                + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  )) &
                + cx2 * (ez(ix+2, iy  ) - ez(ix-1, iy  ))

            bz(ix, iy) = bz(ix, iy) &
                - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                - cx2 * (ey(ix+2, iy  ) - ey(ix-1, iy  )) &
                + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  )) &
                + cy2 * (ex(ix  , iy+2) - ex(ix  , iy-1))
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iy = 0, ny
          cpml_y = hdty / cpml_kappa_by(iy)
          cy1 = c1 * cpml_y
          cy2 = c2 * cpml_y
          cy3 = c3 * cpml_y
          DO ix = 0, nx
            cpml_x = hdtx / cpml_kappa_bx(ix)
            cx1 = c1 * cpml_x
            cx2 = c2 * cpml_x
            cx3 = c3 * cpml_x

            bx(ix, iy) = bx(ix, iy) &
                - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  )) &
                - cy2 * (ez(ix  , iy+2) - ez(ix  , iy-1)) &
                - cy3 * (ez(ix  , iy+3) - ez(ix  , iy-2))

            by(ix, iy) = by(ix, iy) &
                + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  )) &
                + cx2 * (ez(ix+2, iy  ) - ez(ix-1, iy  )) &
                + cx3 * (ez(ix+3, iy  ) - ez(ix-2, iy  ))

            bz(ix, iy) = bz(ix, iy) &
                - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                - cx2 * (ey(ix+2, iy  ) - ey(ix-1, iy  )) &
                - cx3 * (ey(ix+3, iy  ) - ey(ix-2, iy  )) &
                + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  )) &
                + cy2 * (ex(ix  , iy+2) - ex(ix  , iy-1)) &
                + cy3 * (ex(ix  , iy+3) - ex(ix  , iy-2))
          END DO
        END DO
      END IF

      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        cx1 = hdtx
        cy1 = hdty

        IF (maxwell_solver == c_maxwell_solver_yee) THEN
          DO iy = 0, ny
            DO ix = 0, nx
              bx(ix, iy) = bx(ix, iy) &
                  - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  ))

              by(ix, iy) = by(ix, iy) &
                  + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  ))

              bz(ix, iy) = bz(ix, iy) &
                  - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                  + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  ))
            END DO
          END DO
        ELSE
          DO iy = 0, ny
            DO ix = 0, nx
              bx(ix, iy) = bx(ix, iy) &
                  - cy1 * (alphay * (ez(ix  , iy+1) - ez(ix  , iy  ))  &
                         + betayx * (ez(ix+1, iy+1) - ez(ix+1, iy  )   &
                                   + ez(ix-1, iy+1) - ez(ix-1, iy  ))  &
                         + deltay * (ez(ix  , iy+2) - ez(ix  , iy-1)))

              by(ix, iy) = by(ix, iy) &
                  + cx1 * (alphax * (ez(ix+1, iy  ) - ez(ix  , iy  ))  &
                         + betaxy * (ez(ix+1, iy+1) - ez(ix  , iy+1)   &
                                   + ez(ix+1, iy-1) - ez(ix  , iy-1))  &
                         + deltax * (ez(ix+2, iy  ) - ez(ix-1, iy  )))

              bz(ix, iy) = bz(ix, iy) &
                  - cx1 * (alphax * (ey(ix+1, iy  ) - ey(ix  , iy  ))  &
                         + betaxy * (ey(ix+1, iy+1) - ey(ix  , iy+1)   &
                                   + ey(ix+1, iy-1) - ey(ix  , iy-1))  &
                         + deltax * (ey(ix+2, iy  ) - ey(ix-1, iy  ))) &
                  + cy1 * (alphay * (ex(ix  , iy+1) - ex(ix  , iy  ))  &
                         + betayx * (ex(ix+1, iy+1) - ex(ix+1, iy  )   &
                                   + ex(ix-1, iy+1) - ex(ix-1, iy  ))  &
                         + deltay * (ex(ix  , iy+2) - ex(ix  , iy-1)))
            END DO
          END DO
        END IF
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx
        cy1 = c1 * hdty
        cy2 = c2 * hdty

        DO iy = 0, ny
          DO ix = 0, nx
            bx(ix, iy) = bx(ix, iy) &
                - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  )) &
                - cy2 * (ez(ix  , iy+2) - ez(ix  , iy-1))

            by(ix, iy) = by(ix, iy) &
                + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  )) &
                + cx2 * (ez(ix+2, iy  ) - ez(ix-1, iy  ))

            bz(ix, iy) = bz(ix, iy) &
                - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                - cx2 * (ey(ix+2, iy  ) - ey(ix-1, iy  )) &
                + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  )) &
                + cy2 * (ex(ix  , iy+2) - ex(ix  , iy-1))
          END DO
        END DO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        cx1 = c1 * hdtx
        cx2 = c2 * hdtx
        cx3 = c3 * hdtx
        cy1 = c1 * hdty
        cy2 = c2 * hdty
        cy3 = c3 * hdty

        DO iy = 0, ny
          DO ix = 0, nx
            bx(ix, iy) = bx(ix, iy) &
                - cy1 * (ez(ix  , iy+1) - ez(ix  , iy  )) &
                - cy2 * (ez(ix  , iy+2) - ez(ix  , iy-1)) &
                - cy3 * (ez(ix  , iy+3) - ez(ix  , iy-2))

            by(ix, iy) = by(ix, iy) &
                + cx1 * (ez(ix+1, iy  ) - ez(ix  , iy  )) &
                + cx2 * (ez(ix+2, iy  ) - ez(ix-1, iy  )) &
                + cx3 * (ez(ix+3, iy  ) - ez(ix-2, iy  ))

            bz(ix, iy) = bz(ix, iy) &
                - cx1 * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
                - cx2 * (ey(ix+2, iy  ) - ey(ix-1, iy  )) &
                - cx3 * (ey(ix+3, iy  ) - ey(ix-2, iy  )) &
                + cy1 * (ex(ix  , iy+1) - ex(ix  , iy  )) &
                + cy2 * (ex(ix  , iy+2) - ex(ix  , iy-1)) &
                + cy3 * (ex(ix  , iy+3) - ex(ix  , iy-2))
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy

    cnx = hdtx * c**2
    cny = hdty * c**2

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
    hdty = hdt / dy

    cnx = hdtx * c**2
    cny = hdty * c**2

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
        power = -dy * SUM(ey(1,1:ny)  * bz(0,0:ny)  - ez(1,1:ny)  * by(0,0:ny))
      CASE(c_bd_x_max)
        power =  dy * SUM(ey(nx,1:ny) * bz(nx,0:ny) - ez(nx,1:ny) * by(nx,0:ny))
      CASE(c_bd_y_min)
        power = -dx * SUM(ez(1:nx,1)  * bx(0:nx,0)  - ex(1:nx,1)  * bz(0:nx,0))
      CASE(c_bd_y_max)
        power =  dx * SUM(ez(1:nx,ny) * bx(0:nx,ny) - ex(1:nx,ny) * bz(0:nx,ny))
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
