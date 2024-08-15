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

MODULE deck_bremsstrahlung_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bremsstrahlung_deck_initialise, bremsstrahlung_deck_finalise
  PUBLIC :: bremsstrahlung_block_start, bremsstrahlung_block_end
  PUBLIC :: bremsstrahlung_block_handle_element, bremsstrahlung_block_check

CONTAINS

  SUBROUTINE bremsstrahlung_deck_initialise

#ifdef BREMSSTRAHLUNG
    IF (deck_state /= c_ds_first) RETURN
    bremsstrahlung_table_location = 'src/physics_packages/TABLES/br'
    use_bremsstrahlung_recoil = .TRUE.
    use_bremsstrahlung = .FALSE.
    bremsstrahlung_photon_species = -1
#ifndef PHOTONS
    photon_species = -1
#endif
    photon_energy_min_bremsstrahlung = EPSILON(1.0_num)
    bremsstrahlung_start_time = 0.0_num
    photon_weight = 1.0_num
    produce_bremsstrahlung_photons = .FALSE.
    bremsstrahlung_photon_dynamics = .FALSE.
    use_plasma_screening = .FALSE.
#endif

  END SUBROUTINE bremsstrahlung_deck_initialise



  SUBROUTINE bremsstrahlung_deck_finalise

    INTEGER :: io, iu
#ifdef BREMSSTRAHLUNG
    LOGICAL :: exists

    IF (deck_state == c_ds_first) RETURN

    IF (rank == 0 .AND. use_bremsstrahlung) THEN
      INQUIRE(file=TRIM(bremsstrahlung_table_location) // '/br1', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find bremsstrahlung tables in the ', &
              'directory "' // TRIM(bremsstrahlung_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
    END IF

    IF (use_bremsstrahlung) need_random_state = .TRUE.
#else
    IF (use_bremsstrahlung) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_bremsstrahlung=T" in the ', &
              '"bremsstrahlung" block.'
          WRITE(io,*) 'Please recompile with the -DBREMSSTRAHLUNG ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE bremsstrahlung_deck_finalise



  SUBROUTINE bremsstrahlung_block_start

#ifdef BREMSSTRAHLUNG
    ! Enable by default if any bremsstrahlung block is defined
    use_bremsstrahlung = .TRUE.
#endif

  END SUBROUTINE bremsstrahlung_block_start



  SUBROUTINE bremsstrahlung_block_end

  END SUBROUTINE bremsstrahlung_block_end



  FUNCTION bremsstrahlung_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_bremsstrahlung') &
        .OR. str_cmp(element, 'bremsstrahlung') &
        .OR. str_cmp(element, 'enable')) THEN
      use_bremsstrahlung = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef BREMSSTRAHLUNG
    IF (str_cmp(element, 'bremsstrahlung_start_time') &
        .OR. str_cmp(element, 'start_time')) THEN
      bremsstrahlung_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_bremsstrahlung_photons') &
        .OR. str_cmp(element, 'produce_photons')) THEN
      produce_bremsstrahlung_photons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_bremsstrahlung_recoil') &
        .OR. str_cmp(element, 'use_radiation_reaction')) THEN
      use_bremsstrahlung_recoil = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy') &
        .OR. str_cmp(element, 'photon_energy_min_bremsstrahlung')) THEN
      photon_energy_min_bremsstrahlung = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_weight') &
        .OR. str_cmp(element, 'photon_weight_multiplier')) THEN
      photon_weight = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bremsstrahlung_table_location') &
        .OR. str_cmp(element, 'table_location')) THEN
      bremsstrahlung_table_location = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'photon_dynamics') &
        .OR. str_cmp(element, 'bremsstrahlung_photon_dynamics')) THEN
      bremsstrahlung_photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_plasma_screening')) THEN
      use_plasma_screening = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_brem_scatter') &
        .OR. str_cmp(element, 'use_bremsstrahlung_scatter')) THEN
      use_brem_scatter = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_bethe_heitler')) THEN
      use_bethe_heitler = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betheheitler_upscaling')) THEN
      betheheitler_upscaling = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_positron_brem') .OR. &
        str_cmp(element, 'positron_brem')) THEN
      positron_brem = as_logical_print(value, element, errcode)
      RETURN
    END IF
    
    errcode = c_err_unknown_element
#endif

  END FUNCTION bremsstrahlung_block_handle_element



  FUNCTION bremsstrahlung_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef BREMSSTRAHLUNG
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef BREMSSTRAHLUNG
    IF (photon_weight <= 0.0_num) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set the photon_weight to less than or ', &
              'equal to zero. To prevent bremsstrahlung photons from being ', &
              'emitted, set use_bremsstrahlung = F.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    IF (photon_weight > 1.0_num) THEN
      photon_weight = 1.0_num
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'You cannot set photon_weight > 1.0. This variable ', &
              'has been truncated to 1.0.'
        END DO
      END IF
    END IF
#endif

  END FUNCTION bremsstrahlung_block_check

END MODULE deck_bremsstrahlung_block
