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

MODULE deck_collision_block

  USE strings_advanced
  USE collisions

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: collision_deck_initialise, collision_deck_finalise
  PUBLIC :: collision_block_start, collision_block_end
  PUBLIC :: collision_block_handle_element, collision_block_check

  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: coll_pairs_touched
  LOGICAL :: got_nanbu

CONTAINS

  SUBROUTINE collision_deck_initialise

    IF (deck_state == c_ds_first) THEN
      use_collisions = .FALSE.
      use_collisional_ionisation = .FALSE.
      got_nanbu = .FALSE.
    ELSE
      ALLOCATE(coll_pairs_touched(1:n_species, 1:n_species))
      coll_pairs_touched = .FALSE.
      CALL setup_collisions
    END IF

  END SUBROUTINE collision_deck_initialise



  SUBROUTINE collision_deck_finalise

    INTEGER :: i, j
    LOGICAL, SAVE :: first = .TRUE.

    IF (deck_state == c_ds_first) RETURN
    DEALLOCATE(coll_pairs_touched)

    IF (use_collisions) THEN

      ! Switch off collisions if there are no colliding speices
      use_collisions = .FALSE.
      DO j = 1, n_species
        DO i = 1, n_species
          IF (coll_pairs_state(i,j) == c_coll_collide) THEN
            use_collisions = .TRUE.
            EXIT
          END IF
        END DO
      END DO
      use_particle_lists = use_particle_lists .OR. use_collisions
      need_random_state = .TRUE.

      ! Mark particle species with collision types
      DO j = 1, n_species
        DO i = 1, n_species
          IF (coll_pairs_state(i,j) == c_coll_collide) THEN
            species_list(i)%coll_pairwise = .TRUE.
            species_list(j)%coll_pairwise = .TRUE.
            species_list(i)%make_secondary_list = .TRUE.
            species_list(j)%make_secondary_list = .TRUE.
          ELSE IF (coll_pairs_state(i,j) == c_coll_background_1st) THEN
            species_list(i)%coll_background = .TRUE.
            species_list(j)%coll_fast = .TRUE.
          ELSE IF (coll_pairs_state(i,j) == c_coll_background_2nd) THEN
            species_list(i)%coll_fast = .TRUE.
            species_list(j)%coll_background = .TRUE.
          END IF
        END DO
      END DO

      ! Identify species to make secondary lists for in collisional ionisation
      IF (use_collisional_ionisation) THEN
        ! Tag all release species as electrons (required for 
        ! unique_electron_species)
        DO i = 1, n_species
          IF (species_list(i)%ionise) THEN 
            species_list(species_list(i)%release_species)%electron = .TRUE. 
          END IF 
        END DO 

        ! Make secondary lists for ionising particles, and electrons
        DO i = 1, n_species
          IF (species_list(i)%ionise .OR. species_list(i)%electron) THEN
            species_list(i)%make_secondary_list = .TRUE.
          END IF
        END DO

        ! Secondary lists also needed for the final ionisation state, which does 
        ! not ionise, as ionised ions are written to this list
        DO i = 1, n_species 
          IF (species_list(i)%ionise) THEN 
            species_list( &
                species_list(i)%ionise_to_species)%make_secondary_list = .TRUE.
          END IF
        END DO
      END IF

      IF (first) THEN
        first = .FALSE.
        IF (rank == 0 .AND. .NOT.got_nanbu) THEN
          IF (use_nanbu) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'The collision routine now uses the Nanbu-Perez scheme ', &
                'by default, rather than'
            PRINT*, 'Sentoku-Kemp. This method is faster and does not ', &
                'appear to suffer from some'
            PRINT*, 'unusual behaviour exhibited by Sentoku-Kemp under ', &
                'certain conditions.'
            PRINT*, 'To revert to Sentoku-Kemp, specify "use_nanbu = F" in ', &
                'the collisions block.'
            PRINT*, 'To remove this warning message, specify "use_nanbu = T".'
            PRINT*
          END IF
        END IF
      END IF

      coll_n_step = MAX(coll_n_step, 1)
      ci_n_step = MAX(ci_n_step, 1)

    END IF

    use_particle_lists = use_particle_lists .OR. use_collisional_ionisation

  END SUBROUTINE collision_deck_finalise



  SUBROUTINE collision_block_start

  END SUBROUTINE collision_block_start



  SUBROUTINE collision_block_end

  END SUBROUTINE collision_block_end



  FUNCTION collision_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (element == blank .OR. value == blank) RETURN

    ! Performed on second parse to ensure that species are set up first.
    IF (str_cmp(element, 'collide')) THEN
      IF (deck_state /= c_ds_first) THEN
        CALL set_collision_matrix(TRIM(ADJUSTL(value)), errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'use_collisions')) THEN
      use_collisions = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_nanbu')) THEN
      use_nanbu = as_logical_print(value, element, errcode)
      got_nanbu = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'use_cold_correction')) THEN
      use_cold_correction = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'rel_cutoff')) THEN
      rel_cutoff = as_real_print(value, element, errcode)
      use_rel_cutoff = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'coll_n_step') .OR. &
        str_cmp(element, 'n_coll_steps')) THEN
      coll_n_step = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'ci_n_step') .OR. &
        str_cmp(element, 'n_ci_steps')) THEN
      ci_n_step = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'back_update_dt')) THEN
      back_update_dt = as_real_print(value, element, errcode)
      coll_subcycle_back = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'coulomb_log')) THEN
      IF (str_cmp(value, 'auto')) THEN
        coulomb_log_auto = .TRUE.
      ELSE
        coulomb_log_auto = .FALSE.
        coulomb_log = as_real_print(value, element, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'collisional_ionisation') &
        .OR. str_cmp(element, 'collisional_ionization')) THEN
      use_collisional_ionisation = as_logical_print(value, element, errcode)
      IF (use_collisional_ionisation) THEN
#ifndef PER_SPECIES_WEIGHT
        use_collisions = .TRUE.
#else
        IF (as_logical_print(value, element, errcode)) THEN
          errcode = c_err_pp_options_wrong
          extended_error_string = '-DPER_SPECIES_WEIGHT'
        END IF
        use_collisional_ionisation = .FALSE.
#endif
      END IF
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION collision_block_handle_element



  FUNCTION collision_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION collision_block_check



! The following code is all about reading the coll_pairs from the input deck

  SUBROUTINE get_token(str_in, str_out, token_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    CHARACTER(*), INTENT(OUT) :: token_out
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, chr, pos
    CHARACTER(1) :: c

    str_len = LEN(str_in)
    pos = str_len

    DO chr = 1, str_len
      c = str_in(chr:chr)
      IF (c == ' ')  THEN
        pos = chr
        EXIT
      END IF
    END DO

    IF (pos < str_len) THEN
      str_out = TRIM(ADJUSTL(str_in(pos+1:str_len)))
    ELSE
      str_out = ''
    END IF

    token_out = TRIM(str_in(1:pos))

  END SUBROUTINE get_token



  SUBROUTINE set_collision_matrix(str_in, errcode)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: errcode
    CHARACTER(LEN=string_length) :: tstr1, tstr2
    CHARACTER(LEN=string_length) :: species1, species2
    REAL(num) :: collfreq
    INTEGER :: io, iu, sp1, sp2, collstate

    IF (deck_state /= c_ds_last) RETURN

    IF (str_cmp(TRIM(str_in), 'all')) THEN
      coll_pairs = 1.0_num
      coll_pairs_state = c_coll_collide
      RETURN
    END IF

    IF (str_cmp(TRIM(str_in), 'none')) THEN
      coll_pairs = -1.0_num
      coll_pairs_state = c_coll_ignore
      RETURN
    END IF

    CALL get_token(str_in, tstr1, species1, errcode)
    IF (errcode /= 0) RETURN

    sp1 = as_integer(species1, errcode)
    IF (errcode /= 0) RETURN

    CALL get_token(tstr1, tstr2, species2, errcode)
    IF (errcode /= 0) RETURN

    sp2 = as_integer(species2, errcode)
    IF (errcode /= 0) RETURN

    ! Collfreq is the collision frequency, which can be artificially modified by
    ! the user. collstate describes the kind of collisions used
    collfreq = 1.0_num
    collstate = c_coll_collide
    IF (str_cmp(TRIM(tstr2), 'on') .OR. str_cmp(TRIM(tstr2), '')) THEN
      collfreq = 1.0_num
      collstate = c_coll_collide
    ELSE IF (str_cmp(TRIM(tstr2), 'off')) THEN
      collfreq = -1.0_num
      collstate = c_coll_ignore
    ELSE IF (str_cmp(TRIM(tstr2), 'background')) THEN
      collfreq = -1.0_num
      collstate = c_coll_background_2nd
    ELSE
      collfreq = as_real(tstr2, errcode)
      IF (errcode /= 0) RETURN
      IF (collfreq < 0.0_num) collstate = c_coll_ignore
    END IF

    IF (coll_pairs_touched(sp1, sp2) .AND. rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'The collide parameter for ' // TRIM(species1) // ' <-> ' &
            // TRIM(species2)
        WRITE(io,*) 'has been set multiple times!'
        WRITE(io,*) 'Collisions will only be carried out once per species pair.'
        WRITE(io,*) 'Later specifications will always override earlier ones.'
        WRITE(io,*)
      END DO
    END IF

    coll_pairs(sp1, sp2) = collfreq
    coll_pairs(sp2, sp1) = collfreq
    coll_pairs_touched(sp1, sp2) = .TRUE.
    coll_pairs_touched(sp2, sp1) = .TRUE.

    ! Save type of collision
    coll_pairs_state(sp1, sp2) = collstate
    IF (collstate == c_coll_background_2nd) THEN
      coll_pairs_state(sp2, sp1) = c_coll_background_1st
    ELSE
      coll_pairs_state(sp2, sp1) = collstate
    END IF

  END SUBROUTINE set_collision_matrix

END MODULE deck_collision_block
