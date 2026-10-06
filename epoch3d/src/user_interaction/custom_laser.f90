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

MODULE custom_laser

  USE shared_data

  IMPLICIT NONE

  INTEGER, PARAMETER :: custom_laser_lu = 150

  ! EPOCH's "num" kind is hardcoded to KIND(1.d0) (constants.F90), so the
  ! binary profile/phase files are always 8 bytes per value. Not derived via
  ! STORAGE_SIZE (F2008) since the rest of this codebase targets F2003.
  INTEGER, PARAMETER :: real_bytes = 8

  ! Set once, by finalize_custom_laser_domain, after the one-off startup
  ! load balance (pre_load_balance and the particle-count-triggered
  ! balance_workload in epoch3d.F90) has run and the domain decomposition
  ! is settled for the rest of the run. use_pre_balance itself never
  ! changes, so local_slab_window needs this separate flag to tell "the
  ! startup rebalance has not happened yet" (domain may still move --
  ! unsafe to size a slab) from "it already happened" (safe).
  LOGICAL, SAVE :: startup_balance_done = .FALSE.

CONTAINS

  FUNCTION custom_laser_time_profile(laser)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile

    custom_laser_time_profile = 1.0_num

  END FUNCTION custom_laser_time_profile



  ! Entry point called from attach_laser for every laser at deck-parse
  ! time (allow_defer = .TRUE.), and again from
  ! finalize_custom_laser_domain once the domain has settled
  ! (allow_defer absent, i.e. .FALSE.). Dispatches to the spatiotemporal
  ! or static spatial loader; both read raw binary files (see
  ! load_binary_file). Loading happens at deck-parse time so that the
  ! MPI_BCAST calls run when ALL ranks participate, avoiding the deadlock
  ! that occurs if loading is deferred to the per-boundary-cell
  ! timestepping loop -- EXCEPT that a spatiotemporal load is itself
  ! deferred past deck-parse time whenever use_pre_balance or
  ! use_balance might still move the domain (see local_slab_window),
  ! since loading the (potentially huge) file immediately would force
  ! it to be stored in full on every rank, and loading it twice would
  ! be wasteful.
  SUBROUTINE custom_laser_spatial_setup(laser, allow_defer)

    TYPE(laser_block), INTENT(INOUT) :: laser
    LOGICAL, INTENT(IN), OPTIONAL :: allow_defer
    CHARACTER(LEN=c_max_path_length) :: filename
    LOGICAL :: defer_ok

    IF (.NOT. laser%use_custom_profile) RETURN

    defer_ok = .FALSE.
    IF (PRESENT(allow_defer)) defer_ok = allow_defer

    ! Under use_balance the domain keeps moving, but
    ! reslab_custom_laser_files re-enters this routine after every
    ! redistribution, so deferring the first load is still worthwhile.
    IF (defer_ok .AND. laser%use_spatiotemporal &
        .AND. (use_pre_balance .OR. use_balance)) RETURN

    IF (laser%use_spatiotemporal) THEN
      IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
        filename = laser%profile_data_file
      ELSE
        filename = 'temporal_spatial_profile.dat'
      END IF
      CALL load_temporal_spatial_profile(laser, filename)

      IF (laser%use_phase_from_file) THEN
        IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
          filename = laser%phase_data_file
        ELSE
          filename = 'phase_profile.dat'
        END IF
        CALL load_phase_profile(laser, filename)
      END IF

      RETURN
    END IF

    CALL load_spatial_fields(laser)

  END SUBROUTINE custom_laser_spatial_setup



  ! Abort with a clear error unless the deck declared a valid uniform grid
  ! for this laser's binary profile/phase files. Required because the files
  ! carry no embedded shape header (per EPOCH's documented binary-file
  ! convention). need_time selects the spatiotemporal variant, which
  ! additionally requires n_t_points and t_start < t_end.
  SUBROUTINE check_file_grid_declared(laser, need_time)

    TYPE(laser_block), INTENT(IN) :: laser
    LOGICAL, INTENT(IN) :: need_time
    LOGICAL :: ok
    INTEGER :: mpi_err

    ok = laser%n_tr1_points >= 2 .AND. laser%n_tr2_points >= 2 &
        .AND. laser%profile_tr1_max > laser%profile_tr1_min &
        .AND. laser%profile_tr2_max > laser%profile_tr2_min
    IF (need_time) THEN
      ok = ok .AND. laser%n_t_points >= 2 .AND. laser%t_end > laser%t_start
    END IF
    IF (ok) RETURN

    IF (rank == 0) THEN
      IF (need_time) THEN
        PRINT *, 'ERROR: use_spatiotemporal_profile = T requires ' &
            // 'n_transverse1_points, n_transverse2_points and ' &
            // 'n_t_points (each >= 2), both pairs of transverse bounds ' &
            // '(max > min) and t_start < t_end in the laser block.'
      ELSE
        PRINT *, 'ERROR: use_custom_profile = T requires ' &
            // 'n_transverse1_points and n_transverse2_points ' &
            // '(each >= 2) and both pairs of transverse bounds ' &
            // '(max > min) in the laser block.'
      END IF
    END IF
    CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)

  END SUBROUTINE check_file_grid_declared



  ! Read a headerless raw binary array (access='stream', n_values of
  ! EPOCH's 8-byte REAL(num)) on rank 0 and broadcast it to all ranks.
  ! Aborts with a clear error if the file is missing or its size does not
  ! match the deck-declared shape -- the only available sanity check, since
  ! the file embeds no shape metadata. Absolute paths are used directly;
  ! relative paths are resolved from data_dir.
  SUBROUTINE load_binary_file(filename, array, n_values)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n_values
    REAL(num), DIMENSION(n_values), INTENT(OUT) :: array
    INTEGER :: io_err, mpi_err
    INTEGER(KIND=8) :: expected_bytes, actual_bytes
    LOGICAL :: file_exists
    CHARACTER(LEN=c_max_path_length) :: full_filename

    IF (filename(1:1) == '/') THEN
      full_filename = TRIM(filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(filename)
    END IF

    IF (rank == 0) THEN
      expected_bytes = INT(n_values, 8) * INT(real_bytes, 8)

      INQUIRE(FILE=TRIM(full_filename), EXIST=file_exists, SIZE=actual_bytes)
      IF (.NOT. file_exists) THEN
        PRINT *, 'ERROR: Could not find ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
      IF (actual_bytes /= expected_bytes) THEN
        PRINT *, 'ERROR: ', TRIM(full_filename), ' is ', actual_bytes, &
            ' bytes; expected ', expected_bytes, &
            ' (product of the deck-declared point counts * 8 bytes)'
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
          ACCESS='STREAM', FORM='UNFORMATTED', ACTION='READ', IOSTAT=io_err)
      IF (io_err /= 0) THEN
        PRINT *, 'ERROR: Could not open ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      READ(custom_laser_lu) array
      CLOSE(custom_laser_lu)
    END IF

    CALL MPI_BCAST(array, n_values, mpireal, 0, mpi_comm_world, mpi_err)

  END SUBROUTINE load_binary_file



  ! Determine which window of the spatiotemporal file grid this rank
  ! needs to keep in memory. The sampler is only ever called for this
  ! rank's own boundary cells, so it is enough to store a slab covering
  ! the local transverse coordinate range plus a one-cell margin, rather
  ! than the full plane on every rank. Ranks that do not own the laser's
  ! boundary face never sample it at all and store an empty slab. Fall
  ! back to the full plane for a moving window shifting x on a y/z-
  ! boundary laser, or until the domain settles for the first time under
  ! use_pre_balance and/or use_balance (startup_balance_done, set by
  ! finalize_custom_laser_domain). After that, use_balance stays
  ! windowed: reslab_custom_laser_files re-derives this window and
  ! reloads the slab after every subsequent redistribution.
  SUBROUTINE local_slab_window(laser, i1_lo, i1_hi, i2_lo, i2_hi)

    TYPE(laser_block), INTENT(IN) :: laser
    INTEGER, INTENT(OUT) :: i1_lo, i1_hi, i2_lo, i2_hi
    REAL(num) :: p1_lo, p1_hi, p2_lo, p2_hi, d1, d2
    LOGICAL :: x_is_transverse
    INTEGER :: idx

    IF (.NOT. is_boundary(laser%boundary)) THEN
      ! Empty slab: this rank never samples this laser
      i1_lo = 1
      i1_hi = 0
      i2_lo = 1
      i2_hi = 0
      RETURN
    END IF

    x_is_transverse = laser%boundary /= c_bd_x_min &
        .AND. laser%boundary /= c_bd_x_max
    IF (((use_balance .OR. use_pre_balance) .AND. .NOT. startup_balance_done) &
        .OR. (move_window .AND. x_is_transverse)) THEN
      i1_lo = 1
      i1_hi = laser%n_tr1_points
      i2_lo = 1
      i2_hi = laser%n_tr2_points
      RETURN
    END IF

    ! Local transverse patch actually sampled (cell centres 0:n)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        p1_lo = y(0)
        p1_hi = y(ny)
        p2_lo = z(0)
        p2_hi = z(nz)
      CASE(c_bd_y_min, c_bd_y_max)
        p1_lo = x(0)
        p1_hi = x(nx)
        p2_lo = z(0)
        p2_hi = z(nz)
      CASE(c_bd_z_min, c_bd_z_max)
        p1_lo = x(0)
        p1_hi = x(nx)
        p2_lo = y(0)
        p2_hi = y(ny)
    END SELECT

    d1 = (laser%profile_tr1_max - laser%profile_tr1_min) &
        / REAL(laser%n_tr1_points - 1, num)
    d2 = (laser%profile_tr2_max - laser%profile_tr2_min) &
        / REAL(laser%n_tr2_points - 1, num)

    ! Same index formula and clamping as sample_file_matrix; the sampler
    ! touches indices [idx, idx+1], widened here by one cell each side
    idx = MAX(1, MIN(INT((p1_lo - laser%profile_tr1_min) / d1) + 1, &
        laser%n_tr1_points - 1))
    i1_lo = MAX(1, idx - 1)
    idx = MAX(1, MIN(INT((p1_hi - laser%profile_tr1_min) / d1) + 1, &
        laser%n_tr1_points - 1))
    i1_hi = MIN(laser%n_tr1_points, idx + 2)

    idx = MAX(1, MIN(INT((p2_lo - laser%profile_tr2_min) / d2) + 1, &
        laser%n_tr2_points - 1))
    i2_lo = MAX(1, idx - 1)
    idx = MAX(1, MIN(INT((p2_hi - laser%profile_tr2_min) / d2) + 1, &
        laser%n_tr2_points - 1))
    i2_hi = MIN(laser%n_tr2_points, idx + 2)

  END SUBROUTINE local_slab_window



  ! Called once from epoch3d.F90, after the one-off startup load balance
  ! (pre_load_balance and the particle-count-triggered balance_workload)
  ! has run and the domain settles for the first time (use_balance may
  ! move it again later -- see reslab_custom_laser_files).
  ! custom_laser_spatial_setup defers a spatiotemporal laser's file load
  ! past deck-parse time whenever use_pre_balance or use_balance might
  ! still move the domain, rather than loading the (potentially huge)
  ! file immediately and forcing local_slab_window to store it in full
  ! on every rank. This performs those deferred loads now, against the
  ! settled domain, so the true per-rank window can be used instead. A
  ! no-op for every other laser: no custom profile, the static spatial
  ! path (never deferred), or already loaded eagerly at deck-parse time
  ! (use_pre_balance = F and use_balance = F).
  SUBROUTINE finalize_custom_laser_domain

    TYPE(laser_block), POINTER :: laser

    startup_balance_done = .TRUE.

    laser => lasers
    DO WHILE (ASSOCIATED(laser))
      IF (laser%use_custom_profile .AND. laser%use_spatiotemporal &
          .AND. .NOT. laser%profile_loaded) THEN
        CALL custom_laser_spatial_setup(laser)
      END IF
      laser => laser%next
    END DO

  END SUBROUTINE finalize_custom_laser_domain



  ! Called from balance_workload (balance.F90) immediately after
  ! redistribute_domain, on every rank, whenever use_balance has just
  ! moved the domain. Re-derives each spatiotemporal custom-laser's
  ! per-rank window and reloads its file against the new domain, so the
  ! windowed slab stays valid instead of falling back to the full plane
  ! for the rest of the run. A no-op for every laser without a per-rank
  ! slab already loaded.
  SUBROUTINE reslab_custom_laser_files

    TYPE(laser_block), POINTER :: laser

    laser => lasers
    DO WHILE (ASSOCIATED(laser))
      IF (laser%use_custom_profile .AND. laser%use_spatiotemporal &
          .AND. laser%profile_loaded) THEN
        DEALLOCATE(laser%file_field_matrix)
        laser%profile_loaded = .FALSE.
        IF (laser%phase_loaded) THEN
          DEALLOCATE(laser%file_phase_matrix)
          laser%phase_loaded = .FALSE.
        END IF
        CALL custom_laser_spatial_setup(laser)
      END IF
      laser => laser%next
    END DO

  END SUBROUTINE reslab_custom_laser_files



  ! Load one spatiotemporal binary file into a per-rank slab, stored on
  ! laser%file_phase_matrix (load_phase = .TRUE.) or
  ! laser%file_field_matrix (.FALSE.). The slab is allocated with its
  ! global index bounds, which the pointer keeps, so the sampler needs no
  ! index translation. Rank 0 reads the file one time-slice at a time and
  ! broadcasts it; each rank copies out only its local window (see
  ! local_slab_window), so no rank ever has to hold the full 3D array --
  ! just its slab plus one transient slice. File existence/size checks as
  ! in load_binary_file.
  SUBROUTINE load_spatiotemporal_file(laser, filename, load_phase)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: load_phase
    REAL(num), DIMENSION(:,:,:), POINTER :: matrix
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: slice
    INTEGER :: i1_lo, i1_hi, i2_lo, i2_hi, it, n1, n2
    INTEGER :: io_err, mpi_err
    INTEGER(KIND=8) :: expected_bytes, actual_bytes
    LOGICAL :: file_exists
    CHARACTER(LEN=c_max_path_length) :: full_filename

    n1 = laser%n_tr1_points
    n2 = laser%n_tr2_points

    IF (filename(1:1) == '/') THEN
      full_filename = TRIM(filename)
    ELSE
      full_filename = TRIM(data_dir) // '/' // TRIM(filename)
    END IF

    IF (rank == 0) THEN
      expected_bytes = INT(n1, 8) * INT(n2, 8) &
          * INT(laser%n_t_points, 8) * INT(real_bytes, 8)

      INQUIRE(FILE=TRIM(full_filename), EXIST=file_exists, SIZE=actual_bytes)
      IF (.NOT. file_exists) THEN
        PRINT *, 'ERROR: Could not find ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
      IF (actual_bytes /= expected_bytes) THEN
        PRINT *, 'ERROR: ', TRIM(full_filename), ' is ', actual_bytes, &
            ' bytes; expected ', expected_bytes, &
            ' (product of the deck-declared point counts * 8 bytes)'
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF

      OPEN(UNIT=custom_laser_lu, FILE=TRIM(full_filename), STATUS='OLD', &
          ACCESS='STREAM', FORM='UNFORMATTED', ACTION='READ', IOSTAT=io_err)
      IF (io_err /= 0) THEN
        PRINT *, 'ERROR: Could not open ', TRIM(full_filename)
        CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
      END IF
    END IF

    CALL local_slab_window(laser, i1_lo, i1_hi, i2_lo, i2_hi)
    ALLOCATE(matrix(i1_lo:i1_hi, i2_lo:i2_hi, laser%n_t_points))

    ALLOCATE(slice(n1, n2))
    DO it = 1, laser%n_t_points
      IF (rank == 0) READ(custom_laser_lu) slice
      CALL MPI_BCAST(slice, n1 * n2, mpireal, 0, mpi_comm_world, mpi_err)
      matrix(:, :, it) = slice(i1_lo:i1_hi, i2_lo:i2_hi)
    END DO
    DEALLOCATE(slice)

    IF (rank == 0) CLOSE(custom_laser_lu)

    IF (load_phase) THEN
      laser%file_phase_matrix => matrix
    ELSE
      laser%file_field_matrix => matrix
    END IF

  END SUBROUTINE load_spatiotemporal_file



  ! Load a 3D spatiotemporal amplitude profile from a raw binary file into
  ! the given laser block: Fortran column-major order with tr1
  ! fastest-varying and time slowest, i.e. n_tr1 * n_tr2 * n_t values
  ! written as a single array. Stored as a per-rank slab (see
  ! load_spatiotemporal_file). Only loads once per laser (guarded by
  ! laser%profile_loaded).
  SUBROUTINE load_temporal_spatial_profile(laser, profile_filename)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: profile_filename

    IF (laser%profile_loaded) RETURN

    CALL check_file_grid_declared(laser, .TRUE.)

    CALL load_spatiotemporal_file(laser, profile_filename, .FALSE.)

    laser%profile_loaded = .TRUE.

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 3D Spatiotemporal Profile Loaded Successfully! <<<'
      PRINT *, '    Grid Size: ', laser%n_tr1_points, ' x ', &
          laser%n_tr2_points, ' (Spatial) x ', laser%n_t_points, &
          ' (Temporal)'
    END IF

  END SUBROUTINE load_temporal_spatial_profile



  ! Load a 3D spatiotemporal phase profile from a raw binary file into the
  ! given laser block. Identical file convention and per-rank slab storage
  ! to the amplitude profile (see load_temporal_spatial_profile); shares
  ! the same deck-declared grid. Phase values are read as-is -- the
  ! Python-side converter writes them already in EPOCH's sign/offset
  ! convention (phase = -phi + pi/2).
  SUBROUTINE load_phase_profile(laser, phase_filename)

    TYPE(laser_block), INTENT(INOUT) :: laser
    CHARACTER(LEN=*), INTENT(IN) :: phase_filename

    IF (laser%phase_loaded) RETURN

    CALL check_file_grid_declared(laser, .TRUE.)

    CALL load_spatiotemporal_file(laser, phase_filename, .TRUE.)

    laser%phase_loaded = .TRUE.

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 3D Spatiotemporal Phase Profile Loaded ' &
          // 'Successfully! <<<'
      PRINT *, '    Grid Size: ', laser%n_tr1_points, ' x ', &
          laser%n_tr2_points, ' (Spatial) x ', laser%n_t_points, &
          ' (Temporal)'
    END IF

  END SUBROUTINE load_phase_profile



  ! Static spatial path: load a 2D amplitude plane (and optionally a phase
  ! plane) from raw binary files and bilinearly interpolate them onto
  ! laser%profile / laser%phase once at setup. The file grid is declared in
  ! the deck exactly as on the spatiotemporal path, minus the temporal
  ! elements; the files themselves are headerless (n_tr1 * n_tr2 values,
  ! tr1 fastest-varying).
  SUBROUTINE load_spatial_fields(laser)

    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: plane
    REAL(num), ALLOCATABLE, DIMENSION(:) :: c1, c2
    CHARACTER(LEN=c_max_path_length) :: filename
    INTEGER :: i, n1, n2
    REAL(num) :: d1, d2

    CALL check_file_grid_declared(laser, .FALSE.)

    n1 = laser%n_tr1_points
    n2 = laser%n_tr2_points
    ALLOCATE(plane(n1, n2), c1(n1), c2(n2))

    ! Reconstruct the uniform file-grid axes from the deck declaration
    d1 = (laser%profile_tr1_max - laser%profile_tr1_min) / REAL(n1-1, num)
    d2 = (laser%profile_tr2_max - laser%profile_tr2_min) / REAL(n2-1, num)
    DO i = 1, n1
      c1(i) = laser%profile_tr1_min + REAL(i-1, num) * d1
    END DO
    DO i = 1, n2
      c2(i) = laser%profile_tr2_min + REAL(i-1, num) * d2
    END DO

    IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
      filename = laser%profile_data_file
    ELSE
      filename = 'spatial_profile.dat'
    END IF
    CALL load_binary_file(filename, plane, n1 * n2)
    CALL interp_plane_to_boundary(laser, c1, c2, plane, laser%profile)

    IF (rank == 0) THEN
      PRINT *, '>>> Custom 2D Spatial Profile Loaded Successfully! <<<'
      PRINT *, '    Grid Size: ', n1, ' x ', n2
    END IF

    IF (laser%use_phase_from_file) THEN
      IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
        filename = laser%phase_data_file
      ELSE
        filename = 'phase_profile.dat'
      END IF
      CALL load_binary_file(filename, plane, n1 * n2)
      CALL interp_plane_to_boundary(laser, c1, c2, plane, laser%phase)

      IF (rank == 0) THEN
        PRINT *, '>>> Custom 2D Spatial Phase Profile Loaded ' &
            // 'Successfully! <<<'
        PRINT *, '    Grid Size: ', n1, ' x ', n2
      END IF
    END IF

    DEALLOCATE(plane, c1, c2)

  END SUBROUTINE load_spatial_fields



  ! Bilinearly interpolate a file-grid plane onto the local section of the
  ! boundary this laser is attached to, writing into 'dest' (laser%profile
  ! or laser%phase, passed as a pointer so its original bounds are kept).
  SUBROUTINE interp_plane_to_boundary(laser, c1, c2, plane, dest)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), DIMENSION(:), INTENT(IN) :: c1, c2
    REAL(num), DIMENSION(:,:), INTENT(IN) :: plane
    REAL(num), DIMENSION(:,:), POINTER :: dest
    INTEGER :: i, j, n1, n2

    n1 = SIZE(c1)
    n2 = SIZE(c2)

    ! The dest array index convention matches allocate_with_boundary:
    !   x_min/x_max -> dest(0:ny, 0:nz), coord1=y, coord2=z
    !   y_min/y_max -> dest(0:nx, 0:nz), coord1=x, coord2=z
    !   z_min/z_max -> dest(0:nx, 0:ny), coord1=x, coord2=y
    ! Cell-centre coordinates are used, consistent with the analytical
    ! evaluator which resolves deck variables at y(i), z(j).
    SELECT CASE(laser%boundary)

      CASE(c_bd_x_min, c_bd_x_max)
        DO j = 0, nz
          DO i = 0, ny
            dest(i, j) = interp2d(y(i), z(j), c1, c2, plane, n1, n2)
          END DO
        END DO

      CASE(c_bd_y_min, c_bd_y_max)
        DO j = 0, nz
          DO i = 0, nx
            dest(i, j) = interp2d(x(i), z(j), c1, c2, plane, n1, n2)
          END DO
        END DO

      CASE(c_bd_z_min, c_bd_z_max)
        DO j = 0, ny
          DO i = 0, nx
            dest(i, j) = interp2d(x(i), y(j), c1, c2, plane, n1, n2)
          END DO
        END DO

    END SELECT

  END SUBROUTINE interp_plane_to_boundary



  ! Bilinear interpolation on a 2D regular grid.
  ! Returns the interpolated value at (p1, p2). Clamps to boundary values
  ! for points outside the data range.
  REAL(num) FUNCTION interp2d(p1, p2, c1, c2, vals, n1, n2)

    REAL(num), INTENT(IN) :: p1, p2
    INTEGER, INTENT(IN) :: n1, n2
    REAL(num), DIMENSION(n1), INTENT(IN) :: c1
    REAL(num), DIMENSION(n2), INTENT(IN) :: c2
    REAL(num), DIMENSION(n1, n2), INTENT(IN) :: vals

    INTEGER :: i1, i2
    REAL(num) :: u, v, q11, q21, q12, q22

    ! Direct index calculation (O(1)) — valid for uniform grids
    i1 = INT((p1 - c1(1)) / (c1(2) - c1(1))) + 1
    i2 = INT((p2 - c2(1)) / (c2(2) - c2(1))) + 1

    ! Clamp to valid interpolation range [1, n-1]
    i1 = MAX(1, MIN(i1, n1 - 1))
    i2 = MAX(1, MIN(i2, n2 - 1))

    u = (p1 - c1(i1)) / (c1(i1+1) - c1(i1))
    v = (p2 - c2(i2)) / (c2(i2+1) - c2(i2))

    ! Clamp fractional positions for points outside the data range
    u = MAX(0.0_num, MIN(u, 1.0_num))
    v = MAX(0.0_num, MIN(v, 1.0_num))

    q11 = vals(i1,   i2)
    q21 = vals(i1+1, i2)
    q12 = vals(i1,   i2+1)
    q22 = vals(i1+1, i2+1)

    interp2d = (1.0_num - u) * (1.0_num - v) * q11 &
             + u * (1.0_num - v) * q21             &
             + (1.0_num - u) * v * q12             &
             + u * v * q22

  END FUNCTION interp2d



  ! Trilinear sample of a laser's spatiotemporal file matrix at transverse
  ! position (pos1, pos2) and the current simulation time. Returns zero
  ! outside the deck-declared grid: the amplitude envelope is zero there
  ! too, so a zero phase is immaterial. The uniform grid is reconstructed
  ! from the deck-declared bounds/counts and t_start/t_end -- there is no
  ! stored coordinate array to search. The matrix is passed as a pointer
  ! so that the per-rank slab's global index bounds (see
  ! load_spatiotemporal_file) are preserved.
  REAL(num) FUNCTION sample_file_matrix(laser, matrix, pos1, pos2)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), DIMENSION(:,:,:), POINTER :: matrix
    REAL(num), INTENT(IN) :: pos1, pos2
    INTEGER :: i1, i2, it, mpi_err
    REAL(num) :: d1, d2, dtf, p10, p20, t0, u, v, w
    REAL(num) :: f00, f10, f01, f11

    sample_file_matrix = 0.0_num

    ! Ranks that do not own this laser's boundary face store an empty
    ! slab and never legitimately need a sample
    IF (SIZE(matrix) == 0) RETURN

    ! --- 1. Boundary & Guard Checks ---
    IF (pos1 < laser%profile_tr1_min &
        .OR. pos1 > laser%profile_tr1_max) RETURN
    IF (pos2 < laser%profile_tr2_min &
        .OR. pos2 > laser%profile_tr2_max) RETURN
    IF (time < laser%t_start .OR. time > laser%t_end) RETURN

    ! --- 2. Locate the Bounding Cell Box ---
    d1 = (laser%profile_tr1_max - laser%profile_tr1_min) &
        / REAL(laser%n_tr1_points - 1, num)
    d2 = (laser%profile_tr2_max - laser%profile_tr2_min) &
        / REAL(laser%n_tr2_points - 1, num)
    dtf = (laser%t_end - laser%t_start) / REAL(laser%n_t_points - 1, num)

    i1 = INT((pos1 - laser%profile_tr1_min) / d1) + 1
    i2 = INT((pos2 - laser%profile_tr2_min) / d2) + 1
    it = INT((time - laser%t_start) / dtf) + 1

    ! Clamp to valid interpolation range [1, n-1]
    i1 = MAX(1, MIN(i1, laser%n_tr1_points - 1))
    i2 = MAX(1, MIN(i2, laser%n_tr2_points - 1))
    it = MAX(1, MIN(it, laser%n_t_points - 1))

    ! The slab covers this rank's patch by construction
    ! (local_slab_window); an index outside it means the domain was
    ! reconfigured in a way that routine did not anticipate -- fail
    ! loudly rather than silently inject wrong fields
    IF (i1 < LBOUND(matrix, 1) .OR. i1 + 1 > UBOUND(matrix, 1) &
        .OR. i2 < LBOUND(matrix, 2) .OR. i2 + 1 > UBOUND(matrix, 2)) THEN
      PRINT *, 'ERROR: custom laser profile sampled outside the stored ', &
          'per-rank slab on rank ', rank
      CALL MPI_ABORT(mpi_comm_world, 1, mpi_err)
    END IF

    ! --- 3. Trilinear Interpolation Math ---
    p10 = laser%profile_tr1_min + REAL(i1 - 1, num) * d1
    p20 = laser%profile_tr2_min + REAL(i2 - 1, num) * d2
    t0 = laser%t_start + REAL(it - 1, num) * dtf
    u = (pos1 - p10) / d1
    v = (pos2 - p20) / d2
    w = (time - t0) / dtf

    ! Interpolate along tr1 on the four (tr2, t) corner lines...
    f00 = (1.0_num - u) * matrix(i1, i2,   it  ) + u * matrix(i1+1, i2,   it)
    f10 = (1.0_num - u) * matrix(i1, i2+1, it  ) + u * matrix(i1+1, i2+1, it)
    f01 = (1.0_num - u) * matrix(i1, i2,   it+1) &
        + u * matrix(i1+1, i2,   it+1)
    f11 = (1.0_num - u) * matrix(i1, i2+1, it+1) &
        + u * matrix(i1+1, i2+1, it+1)

    ! ...then bilinearly across tr2 and time
    sample_file_matrix = (1.0_num - w) * ((1.0_num - v) * f00 + v * f10) &
        + w * ((1.0_num - v) * f01 + v * f11)

  END FUNCTION sample_file_matrix



  ! Trilinear interpolation of this laser's spatiotemporal amplitude
  ! profile at transverse position (pos1, pos2) and the current simulation
  ! time. pos1/pos2 are the two in-plane boundary coordinates (e.g. y and z
  ! for an x_min/x_max laser). The data is normally preloaded by
  ! custom_laser_spatial_setup; the load here is a fallback guard only.
  REAL(num) FUNCTION custom_laser_profile(laser, pos1, pos2)

    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos1, pos2
    CHARACTER(LEN=c_max_path_length) :: fname

    IF (.NOT. laser%profile_loaded) THEN
      IF (LEN_TRIM(laser%profile_data_file) > 0) THEN
        fname = laser%profile_data_file
      ELSE
        fname = 'temporal_spatial_profile.dat'
      END IF
      CALL load_temporal_spatial_profile(laser, fname)
    END IF

    custom_laser_profile = sample_file_matrix(laser, &
        laser%file_field_matrix, pos1, pos2)

  END FUNCTION custom_laser_profile



  ! Trilinear interpolation of this laser's spatiotemporal phase profile at
  ! transverse position (pos1, pos2) and the current simulation time.
  ! Mirrors custom_laser_profile exactly, but reads from
  ! laser%file_phase_matrix (loaded by load_phase_profile).
  REAL(num) FUNCTION custom_laser_phase(laser, pos1, pos2)

    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num), INTENT(IN) :: pos1, pos2
    CHARACTER(LEN=c_max_path_length) :: fname

    IF (.NOT. laser%phase_loaded) THEN
      IF (LEN_TRIM(laser%phase_data_file) > 0) THEN
        fname = laser%phase_data_file
      ELSE
        fname = 'phase_profile.dat'
      END IF
      CALL load_phase_profile(laser, fname)
    END IF

    custom_laser_phase = sample_file_matrix(laser, &
        laser%file_phase_matrix, pos1, pos2)

  END FUNCTION custom_laser_phase

END MODULE custom_laser
