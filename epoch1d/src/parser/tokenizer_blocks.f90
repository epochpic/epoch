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

MODULE tokenizer_blocks

  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: c_max_new_func = 256
  INTEGER :: n_new_func = 0
  TYPE(string_type), DIMENSION(c_max_new_func) :: new_func_name
  INTEGER, DIMENSION(c_max_new_func) :: new_func_code

  INTEGER, PARAMETER :: c_max_new_const = 256
  INTEGER :: n_new_constant = 0
  TYPE(string_type), DIMENSION(c_max_new_const) :: new_constant_name
  INTEGER, DIMENSION(c_max_new_const) :: new_constant_code
  INTEGER :: last_block_type, tokenize_stagger

CONTAINS

  ! Functions to register new functions and constants
  FUNCTION register_function(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: register_function

    IF (n_new_func == c_max_new_func) THEN
      register_function = -1
      RETURN
    END IF
    n_new_func = n_new_func+1
    new_func_name(n_new_func)%value = name
    new_func_code(n_new_func) = c_func_custom_lowbound+n_new_func
    register_function = c_func_custom_lowbound+n_new_func

  END FUNCTION register_function



  FUNCTION register_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: register_constant

    IF (n_new_constant == c_max_new_const) THEN
      register_constant = -1
      RETURN
    END IF
    n_new_constant = n_new_constant+1
    new_constant_name(n_new_constant)%value = name
    new_constant_code(n_new_constant) = c_const_custom_lowbound+n_new_constant
    register_constant = c_const_custom_lowbound+n_new_constant

  END FUNCTION register_constant



  FUNCTION as_species(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_species
    INTEGER :: i

    as_species = c_prc_not_this_type
    IF (.NOT. ALLOCATED(species_list)) RETURN

    DO i = 1, n_species
      IF (str_cmp(name, species_list(i)%name)) THEN
        as_species = i
        RETURN
      END IF
    END DO

  END FUNCTION as_species



  FUNCTION as_subset(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_subset
    INTEGER :: i

    as_subset = c_prc_not_this_type
    IF (.NOT. ALLOCATED(subset_list)) RETURN

    DO i = 1, n_subsets
      IF (str_cmp(name, subset_list(i)%name)) THEN
        as_subset = i
        RETURN
      END IF
    END DO

  END FUNCTION as_subset



  FUNCTION as_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_constant
    INTEGER :: i

    as_constant = c_prc_not_this_type

    ! Numeric constants
    IF (str_cmp(name, 'pi')) THEN
      as_constant = c_const_pi

    ELSE IF (str_cmp(name, 'kb')) THEN
      as_constant = c_const_kb

    ELSE IF (str_cmp(name, 'me')) THEN
      as_constant = c_const_me

    ELSE IF (str_cmp(name, 'qe')) THEN
      as_constant = c_const_qe

    ELSE IF (str_cmp(lowercase(name), 't') &
        .OR. str_cmp(lowercase(name), 'true')) THEN
      as_constant = c_const_true

    ELSE IF (str_cmp(lowercase(name), 'f') &
        .OR. str_cmp(lowercase(name), 'false')) THEN
      as_constant = c_const_false

    ELSE IF (str_cmp(name, 'c')) THEN
      as_constant = c_const_c

    ELSE IF (str_cmp(name, 'eps0') &
        .OR. str_cmp(name, 'epsilon0') &
        .OR. str_cmp(name, 'epsilonnought')) THEN
      as_constant = c_const_eps0

    ELSE IF (str_cmp(name, 'mu0') &
        .OR. str_cmp(name, 'munought')) THEN
      as_constant = c_const_mu0

    ELSE IF (str_cmp(name, 'ev')) THEN
      as_constant = c_const_ev

    ELSE IF (str_cmp(name, 'kev')) THEN
      as_constant = c_const_kev

    ELSE IF (str_cmp(name, 'mev')) THEN
      as_constant = c_const_mev

    ELSE IF (str_cmp(name, 'milli')) THEN
      as_constant = c_const_milli

    ELSE IF (str_cmp(name, 'micro') &
        .OR. str_cmp(name, 'micron') &
        .OR. str_cmp(name, 'cm3') &
        .OR. str_cmp(name, 'cc')) THEN
      as_constant = c_const_micro

    ELSE IF (str_cmp(name, 'nano')) THEN
      as_constant = c_const_nano

    ELSE IF (str_cmp(name, 'pico')) THEN
      as_constant = c_const_pico

    ELSE IF (str_cmp(name, 'femto')) THEN
      as_constant = c_const_femto

    ELSE IF (str_cmp(name, 'atto')) THEN
      as_constant = c_const_atto

    ELSE IF (str_cmp(name, 'time')) THEN
      as_constant = c_const_time

    ELSE IF (str_cmp(name, 'x')) THEN
      IF (stagger(c_dir_x,tokenize_stagger)) THEN
        as_constant = c_const_xb
      ELSE
        as_constant = c_const_x
      END IF

    ELSE IF (str_cmp(name, 'xb')) THEN
      as_constant = c_const_xb

    ELSE IF (str_cmp(name, 'ix')) THEN
      as_constant = c_const_ix

    ELSE IF (str_cmp(name, 'nx')) THEN
      as_constant = c_const_nx

    ELSE IF (str_cmp(name, 'dx')) THEN
      as_constant = c_const_dx

    ELSE IF (str_cmp(name, 'x_min') .OR. str_cmp(name, 'x_start')) THEN
      as_constant = c_const_x_min

    ELSE IF (str_cmp(name, 'x_max') .OR. str_cmp(name, 'x_end')) THEN
      as_constant = c_const_x_max

    ELSE IF (str_cmp(name, 'lengthx') .OR. str_cmp(name, 'length_x')) THEN
      as_constant = c_const_lx

    ELSE IF (str_cmp(name, 'nprocx') .OR. str_cmp(name, 'nproc_x')) THEN
      as_constant = c_const_nprocx

    ELSE IF (str_cmp(name, 'nsteps')) THEN
      as_constant = c_const_nsteps

    ELSE IF (str_cmp(name, 't_end')) THEN
      as_constant = c_const_t_end

    ELSE IF (str_cmp(name, 'ndims')) THEN
      as_constant = c_const_ndims

    ! Dumpmask constants
    ELSE IF (str_cmp(name, 'never')) THEN
      as_constant = c_const_io_never

    ELSE IF (str_cmp(name, 'always')) THEN
      as_constant = c_const_io_always

    ELSE IF (str_cmp(name, 'full')) THEN
      as_constant = c_const_io_full

    ELSE IF (str_cmp(name, 'restartable')) THEN
      as_constant = c_const_io_restartable

    ELSE IF (str_cmp(name, 'restart')) THEN
      as_constant = c_const_io_restartable

    ELSE IF (str_cmp(name, 'average')) THEN
      as_constant = c_const_io_average

    ELSE IF (str_cmp(name, 'snapshot')) THEN
      as_constant = c_const_io_snapshot

    ELSE IF (str_cmp(name, 'species')) THEN
      as_constant = c_const_io_species

    ELSE IF (str_cmp(name, 'no_sum')) THEN
      as_constant = c_const_io_no_sum

    ELSE IF (str_cmp(name, 'single')) THEN
      as_constant = c_const_io_dump_single

    ELSE IF (str_cmp(name, 'average_single')) THEN
      as_constant = c_const_io_average_single

    ! Distribution function constants
    ELSE IF (str_cmp(name, 'dir_x')) THEN
      as_constant = c_const_dir_x

    ELSE IF (str_cmp(name, 'dir_px')) THEN
      as_constant = c_const_dir_px

    ELSE IF (str_cmp(name, 'dir_py')) THEN
      as_constant = c_const_dir_py

    ELSE IF (str_cmp(name, 'dir_pz')) THEN
      as_constant = c_const_dir_pz

    ELSE IF (str_cmp(name, 'dir_en') &
        .OR. str_cmp(name, 'dir_energy')) THEN
      as_constant = c_const_dir_en

    ELSE IF (str_cmp(name, 'dir_gamma_m1') &
        .OR. str_cmp(name, 'dir_gamma_minus_one')) THEN
      as_constant = c_const_dir_gamma_m1

    ELSE IF (str_cmp(name, 'dir_xy_angle')) THEN
      as_constant = c_const_dir_xy_angle

    ELSE IF (str_cmp(name, 'dir_yz_angle')) THEN
      as_constant = c_const_dir_yz_angle

    ELSE IF (str_cmp(name, 'dir_zx_angle')) THEN
      as_constant = c_const_dir_zx_angle

    ELSE IF (str_cmp(name, 'dir_mod_p')) THEN
      as_constant = c_const_dir_mod_p

    ELSE IF (str_cmp(name, 'px')) THEN
      as_constant = c_const_px

    ELSE IF (str_cmp(name, 'py')) THEN
      as_constant = c_const_py

    ELSE IF (str_cmp(name, 'pz')) THEN
      as_constant = c_const_pz

    ELSE IF (str_cmp(name, 'yee')) THEN
      as_constant = c_const_maxwell_solver_yee

    ELSE IF (str_cmp(name, 'cowan')) THEN
      as_constant = c_const_maxwell_solver_cowan

    ELSE IF (str_cmp(name, 'pukhov')) THEN
      as_constant = c_const_maxwell_solver_pukhov

    ELSE IF (str_cmp(name, 'lehe_x')) THEN
      as_constant = c_const_maxwell_solver_lehe_x

    ELSE IF (str_cmp(name, 'lehe_y')) THEN
      as_constant = c_const_maxwell_solver_lehe_y

    ELSE IF (str_cmp(name, 'lehe_z')) THEN
      as_constant = c_const_maxwell_solver_lehe_z

    ELSE IF (str_cmp(name, 'lehe')) THEN
      as_constant = c_const_maxwell_solver_lehe

    ELSE IF (str_cmp(name, 'custom')) THEN
      as_constant = c_const_maxwell_solver_custom

    END IF

    ! User submitted constant using 'Register'
    DO i = 1, n_new_constant
      IF (str_cmp(TRIM(name), TRIM(new_constant_name(i)%value))) &
          as_constant = new_constant_code(i)
    END DO

  END FUNCTION as_constant



  FUNCTION as_default_constant(name) RESULT(as_constant)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_constant
    INTEGER :: io, iu
    LOGICAL, SAVE :: warn = .TRUE.

    as_constant = c_prc_not_this_type

    IF (str_cmp(name, 'y')) THEN
      IF (stagger(c_dir_y,tokenize_stagger)) THEN
        as_constant = c_const_yb
      ELSE
        as_constant = c_const_y
      END IF

    ELSE IF (str_cmp(name, 'z')) THEN
      IF (stagger(c_dir_z,tokenize_stagger)) THEN
        as_constant = c_const_zb
      ELSE
        as_constant = c_const_z
      END IF

    ELSE IF (str_cmp(name, 'yb')) THEN
      as_constant = c_const_yb

    ELSE IF (str_cmp(name, 'zb')) THEN
      as_constant = c_const_zb

    ELSE IF (str_cmp(name, 'iy')) THEN
      as_constant = c_const_iy

    ELSE IF (str_cmp(name, 'iz')) THEN
      as_constant = c_const_iz

    ELSE IF (str_cmp(name, 'ny')) THEN
      as_constant = c_const_ny

    ELSE IF (str_cmp(name, 'nz')) THEN
      as_constant = c_const_nz

    ELSE IF (str_cmp(name, 'dy')) THEN
      as_constant = c_const_dy

    ELSE IF (str_cmp(name, 'dz')) THEN
      as_constant = c_const_dz

    ELSE IF (str_cmp(name, 'y_min') .OR. str_cmp(name, 'y_start')) THEN
      as_constant = c_const_y_min

    ELSE IF (str_cmp(name, 'z_min') .OR. str_cmp(name, 'z_start')) THEN
      as_constant = c_const_z_min

    ELSE IF (str_cmp(name, 'y_max') .OR. str_cmp(name, 'y_end')) THEN
      as_constant = c_const_y_max

    ELSE IF (str_cmp(name, 'z_max') .OR. str_cmp(name, 'z_end')) THEN
      as_constant = c_const_z_max

    ELSE IF (str_cmp(name, 'lengthy') .OR. str_cmp(name, 'length_y')) THEN
      as_constant = c_const_ly

    ELSE IF (str_cmp(name, 'lengthz') .OR. str_cmp(name, 'length_z')) THEN
      as_constant = c_const_lz

    ELSE IF (str_cmp(name, 'r_xy')) THEN
      as_constant = c_const_r_xy

    ELSE IF (str_cmp(name, 'r_yz')) THEN
      as_constant = c_const_r_yz

    ELSE IF (str_cmp(name, 'r_xz')) THEN
      as_constant = c_const_r_xz

    ELSE IF (str_cmp(name, 'r_xyz')) THEN
      as_constant = c_const_r_xyz

    ELSE IF (str_cmp(name, 'nprocy') .OR. str_cmp(name, 'nproc_y')) THEN
      as_constant = c_const_nprocy

    ELSE IF (str_cmp(name, 'nprocz') .OR. str_cmp(name, 'nproc_z')) THEN
      as_constant = c_const_nprocz

    ELSE IF (str_cmp(name, 'dir_y')) THEN
      as_constant = c_const_dir_y

    ELSE IF (str_cmp(name, 'dir_z')) THEN
      as_constant = c_const_dir_z

    END IF

    IF (warn .AND. as_constant /= c_prc_not_this_type) THEN
      warn = .FALSE.
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'A default value (set y and z to 0) was used for the ', &
              'constant ', '"' // TRIM(name) // '"'
          WRITE(io,*)
        END DO
      END IF
    END IF

  END FUNCTION as_default_constant



  FUNCTION as_deck_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_deck_constant
    INTEGER :: i

    as_deck_constant = 0

    DO i = 1, n_deck_constants
      IF (str_cmp(TRIM(name), TRIM(deck_constant_list(i)%name))) THEN
        as_deck_constant = i
        RETURN
      END IF
    END DO

  END FUNCTION as_deck_constant



  FUNCTION as_function(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_function
    INTEGER :: i

    as_function = c_prc_not_this_type

    IF (str_cmp(name, '-')) THEN
      as_function = c_func_neg

    ELSE IF (str_cmp(name, 'abs')) THEN
      as_function = c_func_abs

    ELSE IF (str_cmp(name, 'floor')) THEN
      as_function = c_func_floor

    ELSE IF (str_cmp(name, 'ceil')) THEN
      as_function = c_func_ceil

    ELSE IF (str_cmp(name, 'nint')) THEN
      as_function = c_func_nint

    ELSE IF (str_cmp(name, 'sqrt')) THEN
      as_function = c_func_sqrt

    ELSE IF (str_cmp(name, 'sin')) THEN
      as_function = c_func_sine

    ELSE IF (str_cmp(name, 'cos')) THEN
      as_function = c_func_cosine

    ELSE IF (str_cmp(name, 'tan')) THEN
      as_function = c_func_tan

    ELSE IF (str_cmp(name, 'asin')) THEN
      as_function = c_func_arcsine

    ELSE IF (str_cmp(name, 'acos')) THEN
      as_function = c_func_arccosine

    ELSE IF (str_cmp(name, 'atan')) THEN
      as_function = c_func_arctan

    ELSE IF (str_cmp(name, 'atan2')) THEN
      as_function = c_func_arctan2

    ELSE IF (str_cmp(name, 'sinh')) THEN
      as_function = c_func_sinh

    ELSE IF (str_cmp(name, 'cosh')) THEN
      as_function = c_func_cosh

    ELSE IF (str_cmp(name, 'tanh')) THEN
      as_function = c_func_tanh

    ELSE IF (str_cmp(name, 'exp')) THEN
      as_function = c_func_exp

    ELSE IF (str_cmp(name, 'loge')) THEN
      as_function = c_func_loge

    ELSE IF (str_cmp(name, 'log10')) THEN
      as_function = c_func_log10

    ELSE IF (str_cmp(name, 'log_base')) THEN
      as_function = c_func_log_base

    ELSE IF (str_cmp(name, 'min')) THEN
      as_function = c_func_min

    ELSE IF (str_cmp(name, 'max')) THEN
      as_function = c_func_max

    ELSE IF (str_cmp(name, 'gauss')) THEN
      as_function = c_func_gauss

    ELSE IF (str_cmp(name, 'semigauss')) THEN
      as_function = c_func_semigauss

    ELSE IF (str_cmp(name, 'supergauss')) THEN
      as_function = c_func_supergauss

    ELSE IF (str_cmp(name, 'interpolate')) THEN
      as_function = c_func_interpolate

    ELSE IF (str_cmp(name, 'if')) THEN
      as_function = c_func_if

    ELSE IF (str_cmp(name, 'density') &
        .OR. str_cmp(name, 'rho') &
        .OR. str_cmp(name, 'number_density')) THEN
      as_function = c_func_rho

    ELSE IF (str_cmp(name, 'temp_x') &
        .OR. str_cmp(name, 'temp_x_k') &
        .OR. str_cmp(name, 'temperature_x') &
        .OR. str_cmp(name, 'temperature_x_k')) THEN
      as_function = c_func_tempx

    ELSE IF (str_cmp(name, 'temp_y') &
        .OR. str_cmp(name, 'temp_y_k') &
        .OR. str_cmp(name, 'temperature_y') &
        .OR. str_cmp(name, 'temperature_y_k')) THEN
      as_function = c_func_tempy

    ELSE IF (str_cmp(name, 'temp_z') &
        .OR. str_cmp(name, 'temp_z_k') &
        .OR. str_cmp(name, 'temperature_z') &
        .OR. str_cmp(name, 'temperature_z_k')) THEN
      as_function = c_func_tempz

    ELSE IF (str_cmp(name, 'temp_x_ev') &
        .OR. str_cmp(name, 'temperature_x_ev')) THEN
      as_function = c_func_tempx_ev

    ELSE IF (str_cmp(name, 'temp_y_ev') &
        .OR. str_cmp(name, 'temperature_y_ev')) THEN
      as_function = c_func_tempy_ev

    ELSE IF (str_cmp(name, 'temp_z_ev') &
        .OR. str_cmp(name, 'temperature_z_ev')) THEN
      as_function = c_func_tempz_ev

    ELSE IF (str_cmp(name, 'ex')) THEN
      as_function = c_func_ex

    ELSE IF (str_cmp(name, 'ey')) THEN
      as_function = c_func_ey

    ELSE IF (str_cmp(name, 'ez')) THEN
      as_function = c_func_ez

    ELSE IF (str_cmp(name, 'bx')) THEN
      as_function = c_func_bx

    ELSE IF (str_cmp(name, 'by')) THEN
      as_function = c_func_by

    ELSE IF (str_cmp(name, 'bz')) THEN
      as_function = c_func_bz

    ELSE IF (str_cmp(name, 'critical')) THEN
      as_function = c_func_crit

    END IF

    DO i = 1, n_new_func
      IF (str_cmp(TRIM(name), TRIM(new_func_name(i)%value))) THEN
        as_function = new_func_code(i)
      END IF
    END DO

  END FUNCTION as_function



  FUNCTION as_operator(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_operator

    as_operator = c_prc_not_this_type

    IF (str_cmp(name, '+')) THEN
      IF (last_block_type == c_pt_variable &
          .OR. last_block_type == c_pt_constant &
          .OR. last_block_type == c_pt_default_constant &
          .OR. last_block_type == c_pt_deck_constant &
          .OR. last_block_type == c_pt_species &
          .OR. last_block_type == c_pt_subset) THEN
        as_operator = c_opcode_plus
      ELSE
        as_operator = c_opcode_unary_plus
      END IF

    ELSE IF (str_cmp(name, '-'))  THEN
      IF (last_block_type == c_pt_variable &
          .OR. last_block_type == c_pt_constant &
          .OR. last_block_type == c_pt_default_constant &
          .OR. last_block_type == c_pt_deck_constant &
          .OR. last_block_type == c_pt_species &
          .OR. last_block_type == c_pt_subset) THEN
        as_operator = c_opcode_minus
      ELSE
        as_operator = c_opcode_unary_minus
      END IF

    ELSE IF (str_cmp(name, '*')) THEN
      as_operator = c_opcode_times

    ELSE IF (str_cmp(name, '/')) THEN
      as_operator = c_opcode_divide

    ELSE IF (str_cmp(name, '^')) THEN
      as_operator = c_opcode_power

    ELSE IF (str_cmp(name, 'e')) THEN
      as_operator = c_opcode_expo

    ELSE IF (str_cmp(name, 'lt')) THEN
      as_operator = c_opcode_lt

    ELSE IF (str_cmp(name, 'gt')) THEN
      as_operator = c_opcode_gt

    ELSE IF (str_cmp(name, 'eq')) THEN
      as_operator = c_opcode_eq

    ELSE IF (str_cmp(name, 'and')) THEN
      as_operator = c_opcode_and

    ELSE IF (str_cmp(name, 'or')) THEN
      as_operator = c_opcode_or

    END IF

  END FUNCTION as_operator



  SUBROUTINE check_deprecated(name)

    CHARACTER(LEN=*), INTENT(IN) :: name

    IF (str_cmp(name, 'ln')) THEN
      WRITE(*,*) '"' // TRIM(name) // '" is deprecated.'
      WRITE(*,*) 'Use "loge" instead.'
      RETURN
    END IF

  END SUBROUTINE check_deprecated

END MODULE tokenizer_blocks
