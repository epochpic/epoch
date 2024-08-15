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

MODULE deck_qed_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qed_deck_initialise, qed_deck_finalise
  PUBLIC :: qed_block_start, qed_block_end
  PUBLIC :: qed_block_handle_element, qed_block_check

CONTAINS

  SUBROUTINE qed_deck_initialise

#ifdef PHOTONS
    IF (deck_state == c_ds_first) THEN
      qed_table_location = 'src/physics_packages/TABLES'
      use_radiation_reaction = .TRUE.
      use_qed = .FALSE.
      photon_species = -1
      trident_electron_species = -1
      breit_wheeler_electron_species = -1
      trident_positron_species = -1
      breit_wheeler_positron_species = -1
      photon_energy_min = EPSILON(1.0_num)
      qed_start_time = 0.0_num
      produce_pairs = .FALSE.
      use_radiation_reaction = .TRUE.
      produce_photons = .FALSE.
      photon_dynamics = .FALSE.
    END IF
#endif

  END SUBROUTINE qed_deck_initialise



  SUBROUTINE qed_deck_finalise

    INTEGER :: io, iu
#ifdef PHOTONS
    LOGICAL :: exists

    IF (deck_state == c_ds_first) RETURN

    IF (rank == 0 .AND. use_qed) THEN
      INQUIRE(file=TRIM(qed_table_location)//'/hsokolov.table', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find QED tables in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
    END IF

    IF (use_qed) need_random_state = .TRUE.
#else
    IF (use_qed) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_qed=T" in the "qed" block.'
          WRITE(io,*) 'Please recompile with the -DPHOTONS preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE qed_deck_finalise



  SUBROUTINE qed_block_start

  END SUBROUTINE qed_block_start



  SUBROUTINE qed_block_end

  END SUBROUTINE qed_block_end



  FUNCTION qed_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_qed') .OR. str_cmp(element, 'qed')) THEN
      use_qed = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef PHOTONS
    IF (str_cmp(element, 'qed_start_time')) THEN
      qed_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_photons')) THEN
      produce_photons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_radiation_reaction')) THEN
      use_radiation_reaction = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy')) THEN
      photon_energy_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_pairs')) THEN
      produce_pairs = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'pair_upscaling')) THEN
      pair_upscaling = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'qed_table_location')) THEN
      qed_table_location = TRIM(ADJUSTL(value))
      RETURN
    END IF

    IF (str_cmp(element, 'photon_dynamics')) THEN
      photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element
#endif

  END FUNCTION qed_block_handle_element



  FUNCTION qed_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef PHOTONS
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef PHOTONS
    IF (produce_pairs .AND. .NOT. photon_dynamics) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set photon_dynamics=F when ', &
            'produce_pairs=T. Without ', 'photon motion, pair ', &
            'creation will be incorrect.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF
#endif

  END FUNCTION qed_block_check

END MODULE deck_qed_block
