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

MODULE helper

  USE balance
  USE boundary
  USE partlist
  USE simple_io
  USE deltaf_loader

  IMPLICIT NONE

  REAL(num), POINTER :: species_density(:)
  REAL(num), POINTER :: species_temp(:,:)
  REAL(num), POINTER :: species_drift(:,:)

CONTAINS

  SUBROUTINE set_thermal_bcs(ispecies)

    INTEGER, INTENT(IN) :: ispecies
    TYPE(particle_species), POINTER :: species

    ! Set temperature at boundary for thermal bcs.

    species => species_list(ispecies)

    IF (species%bc_particle(c_bd_x_min) == c_bc_thermal) THEN
      species%ext_temp_x_min(:) = species_temp(1,:)
    END IF
    IF (species%bc_particle(c_bd_x_max) == c_bc_thermal) THEN
      species%ext_temp_x_max(:) = species_temp(nx,:)
    END IF

  END SUBROUTINE set_thermal_bcs



  SUBROUTINE set_thermal_bcs_all

    INTEGER :: ispecies

    ! Set temperature at boundary for thermal bcs.

    DO ispecies = 1, n_species
      CALL setup_ic_density(ispecies)
      CALL setup_ic_temp(ispecies)
      CALL setup_ic_drift(ispecies)
      CALL set_thermal_bcs(ispecies)
    END DO

  END SUBROUTINE set_thermal_bcs_all



  SUBROUTINE setup_background_species

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      IF (.NOT.species%background_species) CYCLE

      CALL setup_ic_density(ispecies)

      ALLOCATE(species%background_density(1-ng:nx+ng))
      species%background_density  = species_density
    END DO

  END SUBROUTINE setup_background_species



  SUBROUTINE auto_load

    INTEGER :: ispecies, n
    TYPE(particle_species), POINTER :: species
    TYPE(particle), POINTER :: current
    INTEGER :: i0, i1, iu, io
    TYPE(initial_condition_block), POINTER :: ic

    IF (pre_loading .AND. n_species > 0) THEN
      i0 = 1 - ng
      IF (use_field_ionisation) i0 = -ng
      i1 = 1 - i0

      ALLOCATE(npart_per_cell_array(i0:nx+i1))
      npart_per_cell_array = 0
    ELSE IF (n_species > 0) THEN
      IF (rank == 0) WRITE(*,*) 'Attempting to load particles'
    END IF

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      ic => species%initial_conditions

      CALL setup_ic_density(ispecies)

      IF (species%background_species) THEN
        ALLOCATE(species%background_density(1-ng:nx+ng))
        species%background_density  = species_density
        CYCLE
      END IF

#ifdef PER_SPECIES_WEIGHT
      CALL non_uniform_load_particles(species_density, species, &
          ic%density_min, ic%density_max)
#else
      CALL setup_particle_density(species_density, species, &
          ic%density_min, ic%density_max)
#endif
      IF (pre_loading) CYCLE

      CALL setup_ic_temp(ispecies)
      CALL setup_ic_drift(ispecies)
      CALL set_thermal_bcs(ispecies)

      IF (species%ic_df_type == c_ic_df_thermal) THEN
        DO n = 1, 3
          CALL setup_particle_temperature(species_temp(:,n), n, species, &
              species_drift(:,n))
        END DO
        CALL deltaf_load(ispecies, species_temp, species_drift)
      ELSE IF (species%ic_df_type == c_ic_df_relativistic_thermal) THEN
        CALL setup_particle_temperature_relativistic(species_temp, species, &
            species_drift)
      ELSE IF (species%ic_df_type == c_ic_df_arbitrary) THEN
        CALL setup_particle_dist_fn(species, species_drift)
      END IF

#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      ! For photons, assign additional variable used in photon particle-push
      IF (species_list(ispecies)%species_type == c_species_id_photon) THEN 
        current => species%attached_list%head 
        DO WHILE (ASSOCIATED(current))
          current%particle_energy = SQRT(SUM(current%part_p**2)) * c
          current => current%next
        END DO
      END IF
#endif
    END DO

    IF (pre_loading) RETURN

    IF (rank == 0) THEN
      DO ispecies = 1, n_species
        species => species_list(ispecies)
        IF (species%background_species) CYCLE
        IF (species%count < 0) THEN
          DO iu = 1, nio_units
            io = ios_units(iu)
            WRITE(io,*) 'No particles specified for species ', &
                '"' // TRIM(species%name) // '"'
          END DO
          species%count = 0
        END IF
      END DO
    END IF

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies
    TYPE(initial_condition_block), POINTER :: ic

    DO ispecies = 1, n_species
      ic => species_list(ispecies)%initial_conditions
      NULLIFY(ic%density)
      NULLIFY(ic%temp)
      NULLIFY(ic%drift)

      ic%density_min = EPSILON(1.0_num)
      ic%density_max = HUGE(1.0_num)
      ic%density_back = 0.0_num
      ic%temp_back = 0.0_num
      ic%drift_back = 0.0_num
    END DO

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies
    TYPE(initial_condition_block), POINTER :: ic

    DO ispecies = 1, n_species
      ic => species_list(ispecies)%initial_conditions
      IF (ASSOCIATED(ic%density)) DEALLOCATE(ic%density)
      IF (ASSOCIATED(ic%temp)) DEALLOCATE(ic%temp)
      IF (ASSOCIATED(ic%drift)) DEALLOCATE(ic%drift)
      IF (ALLOCATED(global_species_density)) DEALLOCATE(global_species_density)
      IF (ALLOCATED(global_species_temp)) DEALLOCATE(global_species_temp)
      IF (ALLOCATED(global_species_drift)) DEALLOCATE(global_species_drift)
    END DO

  END SUBROUTINE deallocate_ic



#ifdef PER_SPECIES_WEIGHT
  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(i8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER :: iu, io

    partlist => species%attached_list

    num_valid_cells_local = 0
    density_total = 0.0_num

    DO ix = 1-ng, nx+ng
      IF (density(ix) > density_max) density(ix) = density_max
    END DO ! ix

    DO ix = 1, nx
      IF (density(ix) >= density_min) THEN
        num_valid_cells_local = num_valid_cells_local + 1
        density_total = density_total + density(ix)
      END IF
    END DO ! ix

    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell_average = FLOOR(species%npart_per_cell, num)
    ELSE
      npart_per_cell_average = REAL(species%count, num) &
          / REAL(num_valid_cells_global, num)
    END IF

    IF (npart_per_cell_average <= 0) RETURN

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global / REAL(num_valid_cells_global, num)

    IF (pre_loading) THEN
      DO ix = 1, nx
        npart_per_cell = NINT(density(ix) / density_average &
            * npart_per_cell_average)
        npart_per_cell_array(ix) = &
            npart_per_cell_array(ix) + INT(npart_per_cell)
      END DO ! ix

      RETURN
    END IF

    npart_this_proc_new = 0
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    END DO ! ix

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current => partlist%head
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix) / density_average &
          * npart_per_cell_average)

      ipart = 0
      DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
        ! Even if particles have per particle charge and mass, assume
        ! that initially they all have the same charge and mass (user
        ! can easily over_ride)
        current%charge = species%charge
        current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
        ! Store the number of particles per cell to allow calculation
        ! of phase space volume later
        current%pvol = npart_per_cell
#endif
        current%part_pos = x(ix) + (random() - 0.5_num) * dx

        ipart = ipart + 1
        current => current%next
      END DO
    END DO ! ix

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      CALL destroy_particle(current)
      current => next
    END DO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species
    species%weight = density_total_global * dx / npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      DO iu = 1, nio_units
        io = ios_units(iu)
        WRITE(io,*) 'Loaded ', TRIM(ADJUSTL(string)), &
            ' particles of species ', '"' // TRIM(species%name) // '"'
      END DO
    END IF

    CALL setup_bc_lists
    CALL particle_bcs

  END SUBROUTINE non_uniform_load_particles
#endif



#ifndef PER_SPECIES_WEIGHT
  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(1-ng:), INTENT(IN) :: load_list
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart, npart_per_cell, num_int, num_total, idx
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_this_species, num_new_particles, npart_left
    INTEGER(i8), ALLOCATABLE :: num_valid_cells_all(:), num_idx(:)
    REAL(num) :: valid_cell_frac, num_real, f0, f1
    REAL(num), ALLOCATABLE :: num_frac(:)
    INTEGER(i8) :: cell_x
    INTEGER(i8) :: i, ipos
    INTEGER :: ix
    INTEGER :: ix_min, ix_max
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep
    INTEGER :: iu, io

    npart_this_species = species%count
    IF (npart_this_species <= 0) RETURN

    ix_min = 1
    ix_max = nx

    IF (species%fill_ghosts) THEN
      IF (x_min_boundary) THEN
        IF (injector_boundary(c_bd_x_min) &
            .OR. species%bc_particle(c_bd_x_min) == c_bc_thermal) THEN
          ix_min = ix_min - png
        END IF
      END IF
      IF (x_max_boundary) THEN
        IF (injector_boundary(c_bd_x_max) &
            .OR. species%bc_particle(c_bd_x_max) == c_bc_thermal) THEN
          ix_max = ix_max + png
        END IF
      END IF
    END IF

    num_valid_cells_local = 0
    DO ix = ix_min, ix_max
      IF (load_list(ix)) num_valid_cells_local = num_valid_cells_local + 1
    END DO ! ix

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      num_new_particles = &
          FLOOR(species%npart_per_cell * num_valid_cells_local, KIND=i8)
    ELSE
      ALLOCATE(num_valid_cells_all(nproc), num_idx(nproc), num_frac(nproc))

      ! Calculate global number of particles per cell
      CALL MPI_ALLGATHER(num_valid_cells_local, 1, MPI_INTEGER8, &
          num_valid_cells_all, 1, MPI_INTEGER8, comm, errcode)

      num_valid_cells_global = 0
      DO i = 1,nproc
        num_valid_cells_global = num_valid_cells_global + num_valid_cells_all(i)
      END DO

      IF (num_valid_cells_global == 0) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
              // 'where particles may'
          WRITE(*,*) 'validly be placed for species "' // TRIM(species%name) &
              // '". ', 'Code will now terminate.'
          CALL abort_code(c_err_bad_setup)
        END IF
      END IF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = FLOOR(num_real, KIND=i8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = FLOOR(num_real, KIND=i8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i - 1
        num_total = num_total + num_int
      END DO
      num_total = npart_this_species - num_total

      IF (num_total > 0) THEN
        ! Sort the list of fractions into decreasing order using bubble sort
        sweep = .TRUE.
        DO WHILE(sweep)
          sweep = .FALSE.
          f0 = num_frac(1)
          DO i = 2,nproc
            f1 = num_frac(i)
            IF (f1 > f0) THEN
              num_frac(i-1) = f1
              num_frac(i) = f0
              f1 = f0
              idx = num_idx(i-1)
              num_idx(i-1) = num_idx(i)
              num_idx(i) = idx
              sweep = .TRUE.
            END IF
            f0 = f1
          END DO
        END DO

        ! Accumulate fractional particles until they have all been accounted
        ! for. If any of them have been assigned to the current processor,
        ! add them and exit the loop.

        DO i = 1,nproc
          IF (num_idx(i) == rank) THEN
            num_new_particles = num_new_particles + 1
            EXIT
          END IF
          num_total = num_total - 1
          IF (num_total <= 0) EXIT
        END DO
      END IF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
    END IF

    IF (pre_loading) THEN
      IF (npart_per_cell <= 0) RETURN

      DO ix = ix_min, ix_max
        IF (.NOT. load_list(ix)) CYCLE

        npart_per_cell_array(ix) = &
            npart_per_cell_array(ix) + INT(npart_per_cell)
      END DO ! ix

      RETURN
    END IF

    partlist => species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current => partlist%head
    IF (npart_per_cell > 0) THEN

      DO ix = ix_min, ix_max
        IF (.NOT. load_list(ix)) CYCLE

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
          ! Store the number of particles per cell to allow calculation
          ! of phase space volume later
          current%pvol = npart_per_cell
#endif
          current%part_pos = x(ix) + (random() - 0.5_num) * dx

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        END DO
      END DO ! ix

    END IF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left > 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO ix = ix_min, ix_max
        IF (load_list(ix)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - ix_min
        END IF
      END DO ! ix

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_x = ipos + ix_min

#ifdef PER_PARTICLE_CHARGE_MASS
        ! Even if particles have per particle charge and mass, assume
        ! that initially they all have the same charge and mass (user
        ! can easily over_ride)
        current%charge = species%charge
        current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
        ! Store the number of particles per cell to allow calculation
        ! of phase space volume later
        current%pvol = npart_per_cell
#endif
        current%part_pos = x(cell_x) + (random() - 0.5_num) * dx

        current => current%next
      END DO

      DEALLOCATE(valid_cell_list)
    END IF

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      CALL destroy_particle(current)
      current => next
    END DO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      DO iu = 1, nio_units
        io = ios_units(iu)
        WRITE(io,*) 'Loaded ', TRIM(ADJUSTL(string)), &
            ' particles of species ', '"' // TRIM(species%name) // '"'
      END DO
    END IF

    CALL setup_bc_lists
    CALL particle_bcs

  END SUBROUTINE load_particles
#endif



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart
    INTEGER, DIMENSION(:), ALLOCATABLE :: npart_in_cell
    REAL(num) :: wdata, x0, x1
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, i, isubx
    REAL(num), DIMENSION(:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:), ALLOCATABLE :: density_map
#include "particle_head.inc"

    ALLOCATE(density(1-ng:nx+ng))
    ALLOCATE(density_map(1-ng:nx+ng))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density, ng)

    DO ix = 1-ng, nx+ng
      IF (density(ix) > density_max) density(ix) = density_max
      IF (density(ix) >= density_min) THEN
        density_map(ix) = .TRUE.
      ELSE
        density(ix) = 0.0_num
      END IF
    END DO ! ix

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(npart_in_cell(1-ng:nx+ng))
    npart_in_cell = 0

    partlist => species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, 'Bad Particle'

#include "particle_to_grid.inc"

      ! Calculate density at the particle position
      wdata = 0.0_num
      DO isubx = sf_min, sf_max
        i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i)) i = cell_x + 1 - isubx
#else
        IF (.NOT. density_map(i)) THEN
          i = cell_x + isubx / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
          IF (.NOT. density_map(i)) i = cell_x - isubx / 2
#endif
        END IF
#endif
        wdata = wdata + gx(isubx) * density(i)
      END DO ! isubx

      current%weight = wdata
#ifdef PARTICLE_SHAPE_TOPHAT
      ! For a TOPHAT shape function, cell_x may not be the cell
      ! containing the particle position
      IF (gx(1) > gx(0)) cell_x = cell_x + 1
#endif
      npart_in_cell(cell_x) = npart_in_cell(cell_x) + 1

      current => current%next
      ipart = ipart + 1
    END DO
    DEALLOCATE(density_map)
    DEALLOCATE(density)

    wdata = dx

    partlist => species%attached_list
    ! Second loop renormalises particle weights
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx + 1.5_num)

      current%weight = current%weight * wdata / npart_in_cell(cell_x)

      current => current%next
      ipart = ipart + 1
    END DO

    DEALLOCATE(npart_in_cell)

    ! If you are filling ghost cells to meet an injector
    ! Then you have overfilled by half a cell but need those particles
    ! To calculate weights correctly. Now delete those particles that
    ! Overlap with the injection region
    IF (species%fill_ghosts .AND. use_injectors) THEN
      x0 = x_min
      IF (injector_boundary(c_bd_x_min)) x0 = x0 - 0.5_num * dx * png
      x1 = x_max
      IF (injector_boundary(c_bd_x_max)) x1 = x1 + 0.5_num * dx * png

      current => partlist%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
        IF (current%part_pos < x0 .OR. current%part_pos >= x1) THEN
          CALL remove_particle_from_partlist(partlist, current)
          CALL destroy_particle(current)
        END IF
        current => next
      END DO
    END IF

  END SUBROUTINE setup_particle_density
#endif



  FUNCTION sample_dist_function(axis, dist_fn)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(n_points))

    cdf(1) = dist_fn(1)
    DO ipoint = 2, n_points
      cdf(ipoint) = cdf(ipoint-1) + dist_fn(ipoint)
    END DO

    cdf = cdf / cdf(n_points)

    position = random()
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start + endpoint) / 2

    DO current = 1, n_points-1
      IF (cdf(current) <= position .AND. cdf(current+1) >= position) THEN
        d_cdf = cdf(current+1) - cdf(current)
        sample_dist_function = (axis(current) * (position - cdf(current)) &
            + axis(current+1) * (cdf(current+1) - position)) / d_cdf
        EXIT
      END IF
    END DO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function



  SUBROUTINE custom_particle_load

    LOGICAL :: file_inconsistencies
    INTEGER :: current_loader_num
    INTEGER :: part_count, read_count, iu, io
    CHARACTER(LEN=string_length) :: stra
    REAL(num), DIMENSION(:), POINTER :: xbuf
    REAL(num), DIMENSION(:), POINTER :: pxbuf, pybuf, pzbuf
#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS) || defined(BREMSSTRAHLUNG)
    REAL(num), DIMENSION(:), POINTER :: wbuf
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    INTEGER(KIND=i4), DIMENSION(:), POINTER :: idbuf4
    INTEGER(KIND=i8), DIMENSION(:), POINTER :: idbuf8
    INTEGER :: i, id_offset
    INTEGER, DIMENSION(:), POINTER :: part_counts
#endif
    TYPE(particle_species), POINTER :: species
    TYPE(custom_particle_loader), POINTER :: curr_loader
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: new_particle

    ! Only needed if there are custom loaders to act on
    IF (n_custom_loaders < 1) RETURN

    ! For every custom loader
    DO current_loader_num = 1, n_custom_loaders
      file_inconsistencies = .FALSE.

      curr_loader => custom_loaders_list(current_loader_num)

      ! Grab associated particle lists
      species => species_list(curr_loader%species_id)
      partlist => species%attached_list

      ! Just to be sure
      CALL destroy_partlist(partlist)
      CALL create_empty_partlist(partlist)

      ! MPI read files
      part_count = load_1d_real_array(curr_loader%x_data, xbuf, &
          curr_loader%x_data_offset, errcode)

#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS) || defined(BREMSSTRAHLUNG)
      read_count = load_1d_real_array(curr_loader%w_data, wbuf, &
          curr_loader%w_data_offset, errcode)
      IF (part_count /= read_count) file_inconsistencies = .TRUE.
#endif
      IF (curr_loader%px_data_given) THEN
        read_count = load_1d_real_array(curr_loader%px_data, pxbuf, &
            curr_loader%px_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      END IF

      IF (curr_loader%py_data_given) THEN
        read_count = load_1d_real_array(curr_loader%py_data, pybuf, &
            curr_loader%py_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      END IF

      IF (curr_loader%pz_data_given) THEN
        read_count = load_1d_real_array(curr_loader%pz_data, pzbuf, &
            curr_loader%pz_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      END IF

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      IF (curr_loader%id_data_given) THEN
        IF (curr_loader%id_data_4byte) THEN
          read_count = load_1d_integer4_array(curr_loader%id_data, idbuf4, &
              curr_loader%id_data_offset, errcode)
        ELSE
          read_count = load_1d_integer8_array(curr_loader%id_data, idbuf8, &
              curr_loader%id_data_offset, errcode)
        END IF
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      END IF
#endif

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, file_inconsistencies, 1, MPI_LOGICAL, &
          MPI_LOR, comm, errcode)

      IF (file_inconsistencies) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Error while loading particles_from_file for species ', &
              TRIM(species%name)
        END IF
        CALL abort_code(c_err_bad_setup)
      END IF

! This is needed to get the IDs assigned properly
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      IF (.NOT.curr_loader%id_data_given) THEN
        ALLOCATE(part_counts(0:nproc-1))
        CALL MPI_ALLGATHER(part_count, 1, MPI_INTEGER4, part_counts, 1, &
            MPI_INTEGER4, comm, errcode)
        id_offset = 0
        DO i = 0, rank
          id_offset = id_offset + part_counts(i)
        END DO
      END IF
#endif

      DO read_count = 1, part_count
        CALL create_particle(new_particle)
        CALL add_particle_to_partlist(partlist, new_particle)

        ! Insert data to particle
        new_particle%part_pos = xbuf(read_count)
#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS) || defined(BREMSSTRAHLUNG)
        new_particle%weight = wbuf(read_count)
#endif
        IF (curr_loader%px_data_given) THEN
          new_particle%part_p(1) = pxbuf(read_count)
        END IF
        IF (curr_loader%py_data_given) THEN
          new_particle%part_p(2) = pybuf(read_count)
        END IF
        IF (curr_loader%pz_data_given) THEN
          new_particle%part_p(3) = pzbuf(read_count)
        END IF
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
        IF (curr_loader%id_data_given) THEN
#if defined(PARTICLE_ID4)
          IF (curr_loader%id_data_4byte) THEN
            new_particle%id = idbuf4(read_count)
          ELSE
            new_particle%id = INT(idbuf8(read_count), i4)
          END IF
#else
          IF (curr_loader%id_data_4byte) THEN
            new_particle%id = INT(idbuf4(read_count), i8)
          ELSE
            new_particle%id = idbuf8(read_count)
          END IF
#endif
        ELSE
#if defined(PARTICLE_ID4)
          new_particle%id = INT(id_offset + read_count, i4)
#else
          new_particle%id = INT(id_offset + read_count, i8)
#endif
        END IF
#endif
        ! Just being careful
        NULLIFY(new_particle)
      END DO

      ! Need to keep totals accurate
      CALL MPI_ALLREDUCE(partlist%count, species%count, 1, MPI_INTEGER8, &
          MPI_SUM, comm, errcode)

      IF (rank == 0) THEN
        CALL integer_as_string(species%count, stra)
        DO iu = 1, nio_units
          io = ios_units(iu)
          WRITE(io,*) 'Inserted ', TRIM(stra), &
              ' custom particles of species "', TRIM(species%name), '"'
        END DO
      END IF
    END DO

    DEALLOCATE(custom_loaders_list)

    CALL distribute_particles

  END SUBROUTINE custom_particle_load



  SUBROUTINE setup_ic_density(ispecies)

    INTEGER, INTENT(IN) :: ispecies
    TYPE(particle_species), POINTER :: species
    TYPE(initial_condition_block), POINTER :: ic
    TYPE(parameter_pack) :: parameters
    INTEGER :: ix

    species => species_list(ispecies)
    ic => species%initial_conditions

    IF (ASSOCIATED(ic%density)) THEN
      species_density => ic%density
      RETURN
    END IF

    IF (use_more_setup_memory) THEN
      ALLOCATE(ic%density(1-ng:nx+ng))
      species_density => ic%density
    ELSE
      IF (.NOT. ALLOCATED(global_species_density)) THEN
        ALLOCATE(global_species_density(1-ng:nx+ng))
      END IF
      species_density => global_species_density
    END IF

    DO ix = 1-ng, nx+ng
      parameters%pack_ix = ix
      species_density(ix) = &
          evaluate_with_parameters(species%density_function, &
              parameters, errcode)
    END DO

  END SUBROUTINE setup_ic_density



  SUBROUTINE setup_ic_temp(ispecies)

    INTEGER, INTENT(IN) :: ispecies
    TYPE(particle_species), POINTER :: species
    TYPE(initial_condition_block), POINTER :: ic
    TYPE(parameter_pack) :: parameters
    INTEGER :: ix, n

    species => species_list(ispecies)
    ic => species%initial_conditions

    IF (ASSOCIATED(ic%temp)) THEN
      species_temp => ic%temp
      RETURN
    END IF

    IF (use_more_setup_memory) THEN
      ALLOCATE(ic%temp(1-ng:nx+ng,3))
      species_temp => ic%temp
    ELSE
      IF (.NOT. ALLOCATED(global_species_temp)) THEN
        ALLOCATE(global_species_temp(1-ng:nx+ng,3))
      END IF
      species_temp => global_species_temp
    END IF

    DO n = 1, 3
      DO ix = 1-ng, nx+ng
        parameters%pack_ix = ix
        species_temp(ix,n) = &
            evaluate_with_parameters(species%temperature_function(n), &
                parameters, errcode)
      END DO
    END DO

  END SUBROUTINE setup_ic_temp



  SUBROUTINE setup_ic_drift(ispecies)

    INTEGER, INTENT(IN) :: ispecies
    TYPE(particle_species), POINTER :: species
    TYPE(initial_condition_block), POINTER :: ic
    TYPE(parameter_pack) :: parameters
    INTEGER :: ix, n

    species => species_list(ispecies)
    ic => species%initial_conditions

    IF (ASSOCIATED(ic%drift)) THEN
      species_drift => ic%drift
      RETURN
    END IF

    IF (use_more_setup_memory) THEN
      ALLOCATE(ic%drift(1-ng:nx+ng,3))
      species_drift => ic%drift
    ELSE
      IF (.NOT. ALLOCATED(global_species_drift)) THEN
        ALLOCATE(global_species_drift(1-ng:nx+ng,3))
      END IF
      species_drift => global_species_drift
    END IF

    DO n = 1, 3
      DO ix = 1-ng, nx+ng
        parameters%pack_ix = ix
        species_drift(ix,n) = &
            evaluate_with_parameters(species%drift_function(n), &
                parameters, errcode)
      END DO
    END DO

  END SUBROUTINE setup_ic_drift

END MODULE helper
