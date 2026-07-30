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

MODULE laser

  USE custom_laser
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%use_phase_function = .FALSE.
    laser%use_profile_function = .FALSE.
    laser%use_omega_function = .FALSE.
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    laser%current_integral_phase = 0.0_num
    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    laser%profile = 1.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  SUBROUTINE setup_laser_phases(phases)

    REAL(num), DIMENSION(:), INTENT(IN) :: phases
    TYPE(laser_block), POINTER :: laser
    INTEGER :: ilas

    ilas = 1
    laser => lasers
    DO WHILE(ASSOCIATED(laser))
      laser%current_integral_phase = phases(ilas)
      ilas = ilas + 1
      laser => laser%next
    END DO

  END SUBROUTINE setup_laser_phases



  SUBROUTINE deallocate_laser(laser)

    TYPE(laser_block), POINTER :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)
    IF (ASSOCIATED(laser%file_field_matrix)) &
        DEALLOCATE(laser%file_field_matrix)
    IF (ASSOCIATED(laser%file_phase_matrix)) &
        DEALLOCATE(laser%file_phase_matrix)
    IF (laser%use_profile_function) &
        CALL deallocate_stack(laser%profile_function)
    IF (laser%use_phase_function) &
        CALL deallocate_stack(laser%phase_function)
    IF (laser%use_time_function) &
        CALL deallocate_stack(laser%time_function)
    IF (laser%use_omega_function) &
        CALL deallocate_stack(laser%omega_function)
    DEALLOCATE(laser)

  END SUBROUTINE deallocate_laser



  SUBROUTINE deallocate_lasers

    TYPE(laser_block), POINTER :: current, next

    current => lasers
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_lasers



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    TYPE(laser_block), POINTER :: laser
    TYPE(laser_block), POINTER :: current
    INTEGER :: boundary

    boundary = laser%boundary

    n_lasers(boundary) = n_lasers(boundary) + 1

    IF (ASSOCIATED(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => laser
    ELSE
      lasers => laser
    END IF

    ! Load any custom file-based profile/phase data now, during setup, so
    ! that the MPI_BCAST calls involve ALL ranks (deferring the load to the
    ! per-boundary-cell timestepping loop risks an MPI deadlock). The
    ! routine dispatches internally on use_custom_profile,
    ! use_spatiotemporal and use_phase_from_file. allow_defer = .TRUE.
    ! lets it instead postpone a spatiotemporal load until
    ! finalize_custom_laser_domain confirms the domain has settled, when
    ! use_pre_balance or use_balance means it might otherwise move (see
    ! local_slab_window and reslab_custom_laser_files in
    ! custom_laser.f90).
    CALL custom_laser_spatial_setup(laser, allow_defer=.TRUE.)

  END SUBROUTINE attach_laser



  ! This routine populates the constant elements of a parameter pack
  ! from a laser

  SUBROUTINE populate_pack_from_laser(laser, parameters)

    TYPE(laser_block), POINTER :: laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0
    parameters%pack_iy = 0
    parameters%pack_iz = 0

    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
      CASE(c_bd_y_min)
        parameters%pack_iy = 0
      CASE(c_bd_y_max)
        parameters%pack_iy = ny
      CASE(c_bd_z_min)
        parameters%pack_iz = 0
      CASE(c_bd_z_max)
        parameters%pack_iz = nz
    END SELECT

  END SUBROUTINE populate_pack_from_laser



  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    IF (laser%use_time_function) THEN
      laser_time_profile = evaluate_with_parameters(laser%time_function, &
          parameters, err)
      RETURN
    END IF

    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, j, err
    TYPE(parameter_pack) :: parameters

    ! On the static spatial path the file phase is interpolated onto
    ! laser%phase once at setup (custom_laser_spatial_setup); there is
    ! nothing to update per step, and any deck 'phase = ...' expression
    ! must not overwrite it.
    IF (laser%use_phase_from_file .AND. .NOT. laser%use_spatiotemporal) RETURN

    err = 0
    CALL populate_pack_from_laser(laser, parameters)

    ! Mirrors laser_update_profile: if the phase comes from a file,
    ! interpolate it from this laser's phase matrix (custom_laser_phase) at
    ! every grid point; otherwise evaluate the deck 'phase = ...'
    ! expression as before. Cell-centre coordinates are used for the
    ! file path, matching the amplitude profile and the analytical
    ! evaluator, which resolves deck variables at y(i), z(j).
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,ny
            IF (laser%use_phase_from_file) THEN
              laser%phase(i,j) = custom_laser_phase(laser, y(i), z(j))
            ELSE
              parameters%pack_iy = i
              laser%phase(i,j) = evaluate_with_parameters( &
                  laser%phase_function, parameters, err)
            END IF
          END DO
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO j = 0,nz
          parameters%pack_iz = j
          DO i = 0,nx
            IF (laser%use_phase_from_file) THEN
              laser%phase(i,j) = custom_laser_phase(laser, x(i), z(j))
            ELSE
              parameters%pack_ix = i
              laser%phase(i,j) = evaluate_with_parameters( &
                  laser%phase_function, parameters, err)
            END IF
          END DO
        END DO
      CASE(c_bd_z_min, c_bd_z_max)
        DO j = 0,ny
          parameters%pack_iy = j
          DO i = 0,nx
            IF (laser%use_phase_from_file) THEN
              laser%phase(i,j) = custom_laser_phase(laser, x(i), y(j))
            ELSE
              parameters%pack_ix = i
              laser%phase(i,j) = evaluate_with_parameters( &
                  laser%phase_function, parameters, err)
            END IF
          END DO
        END DO
    END SELECT

  END SUBROUTINE laser_update_phase



  SUBROUTINE laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, j, err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO j = 0,nz
          DO i = 0,ny
            IF (laser%use_custom_profile .AND. laser%use_spatiotemporal) THEN
              ! Cell-centre coordinates, consistent with the analytical
              ! evaluator which resolves deck variables at y(i), z(j)
              laser%profile(i,j) = custom_laser_profile(laser, y(i), z(j))
            ELSE IF (laser%use_profile_function &
                .OR. laser%profile_function%init) THEN
              parameters%pack_iy = i
              parameters%pack_iz = j
              laser%profile(i,j) = evaluate_with_parameters( &
                  laser%profile_function, parameters, err)
            END IF
          END DO
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO j = 0,nz
          DO i = 0,nx
            IF (laser%use_custom_profile .AND. laser%use_spatiotemporal) THEN
              laser%profile(i,j) = custom_laser_profile(laser, x(i), z(j))
            ELSE IF (laser%use_profile_function &
                .OR. laser%profile_function%init) THEN
              parameters%pack_ix = i
              parameters%pack_iz = j
              laser%profile(i,j) = evaluate_with_parameters( &
                  laser%profile_function, parameters, err)
            END IF
          END DO
        END DO
      CASE(c_bd_z_min, c_bd_z_max)
        DO j = 0,ny
          DO i = 0,nx
            IF (laser%use_custom_profile .AND. laser%use_spatiotemporal) THEN
              laser%profile(i,j) = custom_laser_profile(laser, x(i), y(j))
            ELSE IF (laser%use_profile_function &
                .OR. laser%profile_function%init) THEN
              parameters%pack_ix = i
              parameters%pack_iy = j
              laser%profile(i,j) = evaluate_with_parameters( &
                  laser%profile_function, parameters, err)
            END IF
          END DO
        END DO
    END SELECT

  END SUBROUTINE laser_update_profile



  SUBROUTINE laser_update_omega(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    laser%omega = &
        evaluate_with_parameters(laser%omega_function, parameters, err)
    IF (laser%omega_func_type == c_of_freq) &
        laser%omega = 2.0_num * pi * laser%omega
    IF (laser%omega_func_type == c_of_lambda) &
        laser%omega = 2.0_num * pi * c / laser%omega

  END SUBROUTINE laser_update_omega



  SUBROUTINE update_laser_omegas

    TYPE(laser_block), POINTER :: current

    current => lasers
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_laser_omegas



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(array(1-ng:ny+ng, 1-ng:nz+ng))
    ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(array(1-ng:nx+ng, 1-ng:nz+ng))
    ELSE IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
      ALLOCATE(array(1-ng:nx+ng, 1-ng:ny+ng))
    END IF

  END SUBROUTINE allocate_with_boundary



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => lasers
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny,0:nz))
    ALLOCATE(source2(0:ny,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos-1, 0:ny, 0:nz) = bx_x_min(0:ny, 0:nz)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_min) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,ny
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bz(laserpos-1, 0:ny, 0:nz) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min(0:ny, 0:nz) + c * bz_x_min(0:ny, 0:nz)) &
        - 2.0_num * ey(laserpos, 0:ny, 0:nz) &
        - lz * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, 0:ny, -1:nz-1)) &
        + dt_eps * jy(laserpos, 0:ny, 0:nz) &
        + diff * bz(laserpos, 0:ny, 0:nz))

    by(laserpos-1, 0:ny, 0:nz) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min(0:ny, 0:nz) - c * by_x_min(0:ny, 0:nz)) &
        + 2.0_num * ez(laserpos, 0:ny, 0:nz) &
        - ly * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, -1:ny-1, 0:nz)) &
        - dt_eps * jz(laserpos, 0:ny, 0:nz) &
        + diff * by(laserpos, 0:ny, 0:nz))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:ny,0:nz))
    ALLOCATE(source2(0:ny,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos+1, 0:ny, 0:nz) = bx_x_max(0:ny, 0:nz)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_x_max) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,ny
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bz(laserpos, 0:ny, 0:nz) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(0:ny, 0:nz) - c * bz_x_max(0:ny, 0:nz)) &
        + 2.0_num * ey(laserpos, 0:ny, 0:nz) &
        + lz * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, 0:ny, -1:nz-1)) &
        - dt_eps * jy(laserpos, 0:ny, 0:nz) &
        + diff * bz(laserpos-1, 0:ny, 0:nz))

    by(laserpos, 0:ny, 0:nz) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(0:ny, 0:nz) + c * by_x_max(0:ny, 0:nz)) &
        - 2.0_num * ez(laserpos, 0:ny, 0:nz) &
        + ly * (bx(laserpos, 0:ny, 0:nz) - bx(laserpos, -1:ny-1, 0:nz)) &
        + dt_eps * jz(laserpos, 0:ny, 0:nz) &
        + diff * by(laserpos-1, 0:ny, 0:nz))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:nz))
    ALLOCATE(source2(0:nx,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    by(0:nx, laserpos-1, 0:nz) = by_y_min(0:nx, 0:nz)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_min) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,nx
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bx(0:nx, laserpos-1, 0:nz) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ez_y_min(0:nx, 0:nz) + c * bx_y_min(0:nx, 0:nz)) &
        - 2.0_num * ez(0:nx, laserpos, 0:nz) &
        - lx * (by(0:nx, laserpos, 0:nz) - by(-1:nx-1, laserpos, 0:nz)) &
        + dt_eps * jz(0:nx, laserpos, 0:nz) &
        + diff * bx(0:nx, laserpos, 0:nz))

    bz(0:nx, laserpos-1, 0:nz) = sum * (-4.0_num * source2 &
        - 2.0_num * (ex_y_min(0:nx, 0:nz) - c * bz_y_min(0:nx, 0:nz)) &
        + 2.0_num * ex(0:nx, laserpos, 0:nz) &
        - lz * (by(0:nx, laserpos, 0:nz) - by(0:nx, laserpos, -1:nz-1)) &
        - dt_eps * jx(0:nx, laserpos, 0:nz) &
        + diff * bz(0:nx, laserpos, 0:nz))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_max

    laserpos = ny
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:nz))
    ALLOCATE(source2(0:nx,0:nz))
    source1 = 0.0_num
    source2 = 0.0_num

    by(0:nx, laserpos+1, 0:nz) = by_y_max(0:nx, 0:nz)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_y_max) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,nz
            DO i = 0,nx
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    bx(0:nx, laserpos, 0:nz) = sum * (-4.0_num * source1 &
        - 2.0_num * (ez_y_max(0:nx, 0:nz) - c * bx_y_max(0:nx, 0:nz)) &
        + 2.0_num * ez(0:nx, laserpos, 0:nz) &
        + lx * (by(0:nx, laserpos, 0:nz) - by(-1:nx-1, laserpos, 0:nz)) &
        - dt_eps * jz(0:nx, laserpos, 0:nz) &
        + diff * bx(0:nx, laserpos-1, 0:nz))

    bz(0:nx, laserpos, 0:nz) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ex_y_max(0:nx, 0:nz) + c * bz_y_max(0:nx, 0:nz)) &
        - 2.0_num * ex(0:nx, laserpos, 0:nz) &
        + lz * (by(0:nx, laserpos, 0:nz) - by(0:nx, laserpos, -1:nz-1)) &
        + dt_eps * jx(0:nx, laserpos, 0:nz) &
        + diff * bz(0:nx, laserpos-1, 0:nz))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_y_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE outflow_bcs_z_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_z_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_z_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:ny))
    ALLOCATE(source2(0:nx,0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bz(0:nx, 0:ny, laserpos-1) = bz_z_min(0:nx, 0:ny)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_z_min) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,ny
            DO i = 0,nx
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    by(0:nx, 0:ny, laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ex_z_min(0:nx, 0:ny) + c * by_z_min(0:nx, 0:ny)) &
        - 2.0_num * ex(0:nx, 0:ny, laserpos) &
        - ly * (bz(0:nx, 0:ny, laserpos) - bz(0:nx, -1:ny-1, laserpos)) &
        + dt_eps * jx(0:nx, 0:ny, laserpos) &
        + diff * by(0:nx, 0:ny, laserpos))

    bx(0:nx, 0:ny, laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ey_z_min(0:nx, 0:ny) - c * bx_z_min(0:nx, 0:ny)) &
        + 2.0_num * ey(0:nx, 0:ny, laserpos) &
        - lx * (bz(0:nx, 0:ny, laserpos) - bz(-1:nx-1, 0:ny, laserpos)) &
        - dt_eps * jy(0:nx, 0:ny, laserpos) &
        + diff * bx(0:nx, 0:ny, laserpos))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_z_min, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_z_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_z_min



  SUBROUTINE outflow_bcs_z_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, lz, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i, j
    TYPE(laser_block), POINTER :: current

    n = c_bd_z_max

    laserpos = nz
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_z_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    lz = dtc2 / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(0:nx,0:ny))
    ALLOCATE(source2(0:nx,0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bz(0:nx, 0:ny, laserpos+1) = bz_z_max(0:nx, 0:ny)

    IF (add_laser(n)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == c_bd_z_max) THEN
          ! evaluate the temporal evolution of the laser
          IF (time >= current%t_start .AND. time <= current%t_end) THEN
            IF (current%use_phase_function .OR. current%use_phase_from_file) &
                CALL laser_update_phase(current)
            IF (current%use_profile_function .OR. &
                (current%use_custom_profile .AND. current%use_spatiotemporal)) &
                CALL laser_update_profile(current)
            t_env = laser_time_profile(current) * current%amp
            DO j = 0,ny
            DO i = 0,nx
              base = t_env * current%profile(i,j) &
                * SIN(current%current_integral_phase + current%phase(i,j))
              source1(i,j) = source1(i,j) + base * COS(current%pol_angle)
              source2(i,j) = source2(i,j) + base * SIN(current%pol_angle)
            END DO
            END DO
          END IF
        END IF
        current => current%next
      END DO
    END IF

    by(0:nx, 0:ny, laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ex_z_max(0:nx, 0:ny) - c * by_z_max(0:nx, 0:ny)) &
        + 2.0_num * ex(0:nx, 0:ny, laserpos) &
        + ly * (bz(0:nx, 0:ny, laserpos) - bz(0:nx, -1:ny-1, laserpos)) &
        - dt_eps * jx(0:nx, 0:ny, laserpos) &
        + diff * by(0:nx, 0:ny, laserpos-1))

    bx(0:nx, 0:ny, laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ey_z_max(0:nx, 0:ny) + c * bx_z_max(0:nx, 0:ny)) &
        - 2.0_num * ey(0:nx, 0:ny, laserpos) &
        + lx * (bz(0:nx, 0:ny, laserpos) - bz(-1:nx-1, 0:ny, laserpos)) &
        + dt_eps * jy(0:nx, 0:ny, laserpos) &
        + diff * bx(0:nx, 0:ny, laserpos-1))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_z_max, lasers=lasers)
      ELSE
        CALL calc_absorption(c_bd_z_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_z_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd, factor, lfactor, laser_inject_sum
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, nn, ibc, icell, jcell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num
    mm = 1
    nn = 1

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        dd = dy * dz
        mm = ny
        nn = nz
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        END IF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1, 1:nz  ) + ey(ibc, 1:ny  , 1:nz  ))
        e2 = 0.5_num  * (ez(ibc  , 1:ny  , 0:nz-1) + ez(ibc, 1:ny  , 1:nz  ))
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1, 1:nz  ) + bz(ibc, 0:ny-1, 1:nz  ) &
                       + bz(ibc-1, 1:ny  , 1:nz  ) + bz(ibc, 1:ny  , 1:nz  ))
        b2 = 0.25_num * (by(ibc-1, 1:ny  , 0:nz-1) + by(ibc, 1:ny  , 0:nz-1) &
                       + by(ibc-1, 1:ny  , 1:nz  ) + by(ibc, 1:ny  , 1:nz  ))

      CASE(c_bd_y_min, c_bd_y_max)
        dd = dx * dz
        mm = nx
        nn = nz
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        END IF

        e1 = 0.5_num  * (ez(1:nx  , ibc  , 0:nz-1) + ez(1:nx  , ibc, 1:nz  ))
        e2 = 0.5_num  * (ex(0:nx-1, ibc  , 1:nz  ) + ex(1:nx  , ibc, 1:nz  ))
        b1 = 0.25_num * (bx(1:nx  , ibc-1, 0:nz-1) + bx(1:nx  , ibc, 0:nz-1) &
                       + bx(1:nx  , ibc-1, 1:nz  ) + bx(1:nx  , ibc, 1:nz  ))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1, 1:nz  ) + bz(0:nx-1, ibc, 1:nz  ) &
                       + bz(1:nx  , ibc-1, 1:nz  ) + bz(1:nx  , ibc, 1:nz  ))

      CASE(c_bd_z_min, c_bd_z_max)
        dd = dx * dy
        mm = nx
        nn = ny
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        ibc = 1
        IF (bd == c_bd_z_max) THEN
          dir = -1.0_num
          ibc = nz
        END IF

        e1 = 0.5_num  * (ex(0:nx-1, 1:ny  , ibc  ) + ex(1:nx  , 1:ny  , ibc))
        e2 = 0.5_num  * (ey(1:nx  , 0:ny-1, ibc  ) + ey(1:nx  , 1:ny  , ibc))
        b1 = 0.25_num * (by(0:nx-1, 1:ny  , ibc-1) + by(0:nx-1, 1:ny  , ibc) &
                       + by(1:nx  , 1:ny  , ibc-1) + by(1:nx  , 1:ny  , ibc))
        b2 = 0.25_num * (bx(1:nx  , 0:ny-1, ibc-1) + bx(1:nx  , 0:ny-1, ibc) &
                       + bx(1:nx  , 1:ny  , ibc-1) + bx(1:nx  , 1:ny  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        ALLOCATE(e1(mm, nn), e2(mm, nn), b1(mm, nn), b2(mm, nn))

        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    factor = dt * dd * dir
    laser_absorb_local = laser_absorb_local &
        + (factor / mu0) * SUM(e1 * b1 - e2 * b2)

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        IF (current%boundary == bd) THEN
          laser_inject_sum = 0.0_num
          DO jcell = 1, nn
          DO icell = 1, mm
            laser_inject_sum = laser_inject_sum &
                + current%profile(icell, jcell)**2
          END DO
          END DO
          t_env = laser_time_profile(current)
          lfactor = 0.5_num * epsilon0 * c * ABS(factor) &
              * (t_env * current%amp)**2
          laser_inject_local = laser_inject_local + lfactor * laser_inject_sum
        END IF
        current => current%next
      END DO
    END IF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser
