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

MODULE constants

  USE sdf, ONLY : c_stagger_edge_x, c_stagger_edge_y, c_stagger_edge_z, &
      c_stagger_face_x, c_stagger_face_y, c_stagger_face_z, &
      c_stagger_cell_centre, c_stagger_vertex, c_id_length

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(1.d0)
  INTEGER, PARAMETER :: dbl = KIND(1.d0)
  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18
  REAL(num), PARAMETER :: c_tiny = TINY(1.0_num)
  REAL(num), PARAMETER :: c_largest_number = HUGE(1.0_num)
  REAL(num), PARAMETER :: c_maxexponent = MAXEXPONENT(1.0_num)
  REAL(num), PARAMETER :: c_log2 = 0.69314718055994530941723212145817657_num
  REAL(num), PARAMETER :: c_largest_exp = c_maxexponent * c_log2
  REAL(num), PARAMETER :: &
      c_smallest_exp = (MINEXPONENT(1.0_num) - 1.0_num) * c_log2

  INTEGER, PARAMETER :: c_ndims = 2
  INTEGER, PARAMETER :: c_ndirs = 3 ! Number of directions for velocity

  ! File unit for deck parser diagnostic output
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: du = 40
  INTEGER, PARAMETER :: lu = 41
  INTEGER, PARAMETER :: duc = 42
  INTEGER, PARAMETER :: stat_unit = 20
#ifdef NO_IO
  INTEGER, PARAMETER :: io_units(1) = (/ stdout /)
  INTEGER, PARAMETER :: ios_units(1) = (/ stdout /)
#else
  INTEGER, PARAMETER :: io_units(2) = (/ stdout,du /)
  INTEGER, PARAMETER :: ios_units(2) = (/ stdout,stat_unit /)
#endif
  INTEGER, PARAMETER :: nio_units = SIZE(io_units)

  ! File units for the injectors which read particle info from files
  ! Each injector can open up to 8 files (in 2D)
  ! Multiple file injectors can be specified (up to custom_injector_limit)
  ! File units start at custom_base_unit, then 8 units for injector 1, the next
  ! 8 for injector 2, etc
  ! File units 10000000 to 17999999 are reserved for this routine
  INTEGER, PARAMETER :: custom_injector_limit = 1.0e6
  INTEGER, PARAMETER :: custom_base_unit = 1.0e7
  INTEGER, PARAMETER :: unit_px = 0
  INTEGER, PARAMETER :: unit_py = 1
  INTEGER, PARAMETER :: unit_pz = 2
  INTEGER, PARAMETER :: unit_t  = 3
  INTEGER, PARAMETER :: unit_w  = 4
  INTEGER, PARAMETER :: unit_id = 5
  INTEGER, PARAMETER :: unit_x = 6
  INTEGER, PARAMETER :: unit_y = 7
  INTEGER, PARAMETER :: custom_var_num = 8

  ! Boundary type codes
  INTEGER, PARAMETER :: c_bc_null = -1
  INTEGER, PARAMETER :: c_bc_periodic = 1
  INTEGER, PARAMETER :: c_bc_other = 2
  INTEGER, PARAMETER :: c_bc_simple_laser = 3
  INTEGER, PARAMETER :: c_bc_simple_outflow = 4
  INTEGER, PARAMETER :: c_bc_open = 5
  INTEGER, PARAMETER :: c_bc_dump = 6
  INTEGER, PARAMETER :: c_bc_zero_gradient = 7
  INTEGER, PARAMETER :: c_bc_clamp = 8
  INTEGER, PARAMETER :: c_bc_reflect = 9
  INTEGER, PARAMETER :: c_bc_conduct = 10
  INTEGER, PARAMETER :: c_bc_thermal = 11
  INTEGER, PARAMETER :: c_bc_cpml_laser = 12
  INTEGER, PARAMETER :: c_bc_cpml_outflow = 13
  INTEGER, PARAMETER :: c_bc_mixed = 14
  INTEGER, PARAMETER :: c_bc_heat_bath = 15

  ! Collision types
  INTEGER, PARAMETER :: c_coll_ignore = 0
  INTEGER, PARAMETER :: c_coll_collide = 1
  INTEGER, PARAMETER :: c_coll_background_1st = 2
  INTEGER, PARAMETER :: c_coll_background_2nd = 3

  ! Boundary location codes
  INTEGER, PARAMETER :: c_bd_x_min = 1
  INTEGER, PARAMETER :: c_bd_x_max = 2
  INTEGER, PARAMETER :: c_bd_y_min = 3
  INTEGER, PARAMETER :: c_bd_y_max = 4
  INTEGER, PARAMETER :: c_bd_z_min = 5
  INTEGER, PARAMETER :: c_bd_z_max = 6

  ! Frequency function type codes
  INTEGER, PARAMETER :: c_of_omega = 1
  INTEGER, PARAMETER :: c_of_freq = 2
  INTEGER, PARAMETER :: c_of_lambda = 3

  ! Error codes
  INTEGER, PARAMETER :: c_err_none = 0
  INTEGER, PARAMETER :: c_err_unknown_block = 2**0
  INTEGER, PARAMETER :: c_err_unknown_element = 2**1
  INTEGER, PARAMETER :: c_err_preset_element = 2**2
  INTEGER, PARAMETER :: c_err_preset_element_use_later = 2**3
  INTEGER, PARAMETER :: c_err_bad_value = 2**4
  INTEGER, PARAMETER :: c_err_missing_elements = 2**5
  INTEGER, PARAMETER :: c_err_terminate = 2**6
  INTEGER, PARAMETER :: c_err_required_element_not_set = 2**7
  INTEGER, PARAMETER :: c_err_pp_options_missing = 2**8
  INTEGER, PARAMETER :: c_err_bad_array_length = 2**9
  INTEGER, PARAMETER :: c_err_other = 2**10
  INTEGER, PARAMETER :: c_err_warn_bad_value = 2**11
  INTEGER, PARAMETER :: c_err_generic_warning = 2**12
  INTEGER, PARAMETER :: c_err_generic_error = 2**13
  INTEGER, PARAMETER :: c_err_pp_options_wrong = 2**14
  INTEGER, PARAMETER :: c_err_io_error = 2**15
  INTEGER, PARAMETER :: c_err_bad_setup = 2**16
  INTEGER, PARAMETER :: c_err_window = 2**17
  INTEGER, PARAMETER :: c_err_max = 17

  CHARACTER(LEN=*), PARAMETER :: c_err_char(0:c_err_max) = (/ &
      'unknown_block           ', &
      'unknown_element         ', &
      'preset_element          ', &
      'preset_element_use_later', &
      'bad_value               ', &
      'missing_elements        ', &
      'terminate               ', &
      'required_element_not_set', &
      'pp_options_missing      ', &
      'bad_array_length        ', &
      'other                   ', &
      'warn_bad_value          ', &
      'generic_warning         ', &
      'generic_error           ', &
      'pp_options_wrong        ', &
      'io_error                ', &
      'bad_setup               ', &
      'window                  '/)

  INTEGER, PARAMETER :: c_ds_first = 1
  INTEGER, PARAMETER :: c_ds_last = 2

  INTEGER, PARAMETER :: c_true_value  = 123454321
  INTEGER, PARAMETER :: c_false_value = -123454321

  ! IO codes
  INTEGER, PARAMETER :: c_io_none = 0
  INTEGER, PARAMETER :: c_io_always = 2**0
  INTEGER, PARAMETER :: c_io_full = 2**1
  INTEGER, PARAMETER :: c_io_restartable = 2**2
  INTEGER, PARAMETER :: c_io_species = 2**3
  INTEGER, PARAMETER :: c_io_no_sum = 2**4
  INTEGER, PARAMETER :: c_io_averaged = 2**5
  INTEGER, PARAMETER :: c_io_snapshot = 2**6
  INTEGER, PARAMETER :: c_io_field = 2**7
  INTEGER, PARAMETER :: c_io_dump_single = 2**8
  INTEGER, PARAMETER :: c_io_average_single = 2**9
  INTEGER, PARAMETER :: c_io_never = 2**10

  ! Maxwell Solvers
  INTEGER, PARAMETER :: c_maxwell_solver_custom = -1
  INTEGER, PARAMETER :: c_maxwell_solver_yee = 0
  INTEGER, PARAMETER :: c_maxwell_solver_lehe = 1
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_x = 2
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_y = 3
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_z = 4
  INTEGER, PARAMETER :: c_maxwell_solver_cowan = 5
  INTEGER, PARAMETER :: c_maxwell_solver_pukhov = 6

  ! Particle distribution type codes
  INTEGER, PARAMETER :: c_ic_df_thermal = 1
  INTEGER, PARAMETER :: c_ic_df_relativistic_thermal = 2
  INTEGER, PARAMETER :: c_ic_df_arbitrary = 3

  ! domain codes
  INTEGER, PARAMETER :: c_do_full = 0
  INTEGER, PARAMETER :: c_do_decomposed = 1

  ! Taken from http://physics.nist.gov/cuu/Constants (05/07/2012)
  REAL(num), PARAMETER :: pi = 3.141592653589793238462643383279503_num
  REAL(num), PARAMETER :: q0 = 1.602176565e-19_num ! C (+/- 3.5e-27)
  REAL(num), PARAMETER :: m0 = 9.10938291e-31_num ! kg (+/- 4e-38)
  REAL(num), PARAMETER :: c  = 2.99792458e8_num   ! m/s^2 (exact)
  REAL(num), PARAMETER :: kb = 1.3806488e-23_num  ! J/K (+/- 1.3e-29)
  REAL(num), PARAMETER :: mu0 = 4.e-7_num * pi ! N/A^2 (exact)
  ! epsilon0 = 1.0_num / mu0 / c**2 ! F/m (exact)
  REAL(num), PARAMETER :: epsilon0 = 8.854187817620389850536563031710750e-12_num
  REAL(num), PARAMETER :: h_planck = 6.62606957e-34_num ! J s (+/- 2.9e-41)
  REAL(num), PARAMETER :: ev = q0 ! J
  ! Derived physical parameters used in ionisation
  ! h_bar = h_planck / 2.0_num / pi
  REAL(num), PARAMETER :: h_bar = 1.054571725336289397963133257349698e-34_num
  ! Bohr radius
  ! a0 = 4.0_num * pi * epsilon0 * (h_bar / q0)**2.0_num / m0
  REAL(num), PARAMETER :: a0 = 5.291772101121111395947216558438464e-11_num
  ! Hartree energy = 2 * Rydberg energy
  ! hartree = (h_bar / a0)**2.0_num / m0
  REAL(num), PARAMETER :: hartree = 4.359744350823120007758594450644308e-18_num
  ! Fine structure constant, alpha = h_bar / a0 / m0 / c
  REAL(num), PARAMETER :: alpha = 7.2973525755230202568508027295271584628e-3_num
  ! atomic_time = h_bar / hartree
  REAL(num), PARAMETER :: &
      atomic_time = 2.418884320905619591809404261549867e-17_num
  ! atomic_electric_field = hartree / q0 / a0
  REAL(num), PARAMETER :: &
      atomic_electric_field = 5.142206538736485312185213306837419e11_num
  ! m0 * c
  REAL(num), PARAMETER :: mc0 = 2.73092429345209278e-22_num
  ! m0 * c**2
  REAL(num), PARAMETER :: m0c2 = mc0 * c

  ! Constants used in pair production
#ifdef PHOTONS
  ! b_s = mc0**2 / (h_bar * q0)
  REAL(num), PARAMETER :: b_s = 4.414005028109566589829741352306303e9_num
  ! e_s = b_s * c
  REAL(num), PARAMETER :: e_s = 1.323285417001326061279735961512150e18_num
  ! alpha_f = q0**2 / (2.0_num * epsilon0 * h_planck * c)
  REAL(num), PARAMETER :: alpha_f = 7.297352575523020256850802729527158e-3_num
  ! tau_c = h_bar / (m0 * c**2)
  REAL(num), PARAMETER :: tau_c = 1.288088667367242662108649212042082e-21_num
#endif

  ! Constants used for bremsstrahlung with plasma screening
#ifdef BREMSSTRAHLUNG
  REAL(num), PARAMETER :: e_radius = 0.25_num / pi / epsilon0 / m0 * (q0 / c)**2
  REAL(num), PARAMETER :: log_plasma_screen_const_1 = LOG(1.4_num / alpha)
  REAL(num), PARAMETER :: log_plasma_screen_const_2 = &
      LOG(SQRT(epsilon0 * kb) / q0 * m0 * c * alpha / 1.4_num / h_bar)
#endif

  ! define special particle IDs
  INTEGER, PARAMETER :: c_species_id_generic = 0
  INTEGER, PARAMETER :: c_species_id_photon = 1
  INTEGER, PARAMETER :: c_species_id_electron = 2
  INTEGER, PARAMETER :: c_species_id_positron = 3
  INTEGER, PARAMETER :: c_species_id_proton = 4

  ! direction parameters
  INTEGER, PARAMETER :: c_dir_x = 1
  INTEGER, PARAMETER :: c_dir_y = 2
  INTEGER, PARAMETER :: c_dir_z = 3
  INTEGER, PARAMETER :: c_dir_px = c_ndims + 1
  INTEGER, PARAMETER :: c_dir_py = c_ndims + 2
  INTEGER, PARAMETER :: c_dir_pz = c_ndims + 3
  INTEGER, PARAMETER :: c_dir_en = c_ndims + 4
  INTEGER, PARAMETER :: c_dir_gamma_m1 = c_ndims + 5
  INTEGER, PARAMETER :: c_dir_xy_angle = c_ndims + 6
  INTEGER, PARAMETER :: c_dir_yz_angle = c_ndims + 7
  INTEGER, PARAMETER :: c_dir_zx_angle = c_ndims + 8
  INTEGER, PARAMETER :: c_dir_mod_p = c_ndims + 9

  ! constants defining the maximum number of dimensions and directions
  ! in a distribution function
  INTEGER, PARAMETER :: c_df_maxdirs = c_dir_mod_p
  INTEGER, PARAMETER :: c_df_maxdims = 3

  ! define flags
  INTEGER(i8), PARAMETER :: c_def_particle_debug = 2**0
  INTEGER(i8), PARAMETER :: c_def_field_debug = 2**1
  INTEGER(i8), PARAMETER :: c_def_particle_shape_bspline3 = 2**2
  INTEGER(i8), PARAMETER :: c_def_split_part_after_push = 2**3
  INTEGER(i8), PARAMETER :: c_def_per_particle_weight = 2**4
  INTEGER(i8), PARAMETER :: c_def_particle_count_update = 2**5
  INTEGER(i8), PARAMETER :: c_def_zero_current_particles = 2**6
  INTEGER(i8), PARAMETER :: c_def_particle_probes = 2**7
  INTEGER(i8), PARAMETER :: c_def_per_particle_chargemass = 2**8
  INTEGER(i8), PARAMETER :: c_def_particle_ionise = 2**9
  INTEGER(i8), PARAMETER :: c_def_high_order_smoothing = 2**10
  INTEGER(i8), PARAMETER :: c_def_particle_shape_tophat = 2**11
  INTEGER(i8), PARAMETER :: c_def_parser_debug = 2**12
  INTEGER(i8), PARAMETER :: c_def_particle_id4 = 2**13
  INTEGER(i8), PARAMETER :: c_def_particle_id  = 2**14
  INTEGER(i8), PARAMETER :: c_def_photons = 2**15
  INTEGER(i8), PARAMETER :: c_def_trident_photons = 2**16
  INTEGER(i8), PARAMETER :: c_def_prefetch = 2**17
  INTEGER(i8), PARAMETER :: c_def_mpi_debug = 2**18
  INTEGER(i8), PARAMETER :: c_def_parser_checking = 2**19
  INTEGER(i8), PARAMETER :: c_def_deltaf_method = 2**20
  INTEGER(i8), PARAMETER :: c_def_deltaf_debug = 2**21
  INTEGER(i8), PARAMETER :: c_def_work_done_integrated = 2**22
  INTEGER(i8), PARAMETER :: c_def_hc_push = 2**23
  INTEGER(i8), PARAMETER :: c_def_use_isatty = 2**24
  INTEGER(i8), PARAMETER :: c_def_use_mpi3 = 2**25
  INTEGER(i8), PARAMETER :: c_def_bremsstrahlung = 2**26
  INTEGER(i8), PARAMETER :: c_def_probe_time = 2**27
  INTEGER(i8), PARAMETER :: c_def_transition_rates = 2**28

  ! Stagger types
  INTEGER, PARAMETER :: c_stagger_ex = c_stagger_face_x
  INTEGER, PARAMETER :: c_stagger_ey = c_stagger_face_y
  INTEGER, PARAMETER :: c_stagger_ez = c_stagger_face_z
  INTEGER, PARAMETER :: c_stagger_bx = c_stagger_edge_x
  INTEGER, PARAMETER :: c_stagger_by = c_stagger_edge_y
  INTEGER, PARAMETER :: c_stagger_bz = c_stagger_edge_z
  INTEGER, PARAMETER :: c_stagger_centre = c_stagger_cell_centre
  INTEGER, PARAMETER :: c_stagger_max = c_stagger_vertex
  INTEGER, PARAMETER :: c_stagger_jx = c_stagger_ex
  INTEGER, PARAMETER :: c_stagger_jy = c_stagger_ey
  INTEGER, PARAMETER :: c_stagger_jz = c_stagger_ez

  ! Length of a standard string
  INTEGER, PARAMETER :: string_length = 256

#ifdef PARTICLE_ID4
  INTEGER, PARAMETER :: idkind = i4
#else
  INTEGER, PARAMETER :: idkind = i8
#endif

  !----------------------------------------------------------------------------
  ! Parser data
  !----------------------------------------------------------------------------

  INTEGER, PARAMETER :: c_char_numeric = 1
  INTEGER, PARAMETER :: c_char_alpha = 2
  INTEGER, PARAMETER :: c_char_delimiter = 3
  INTEGER, PARAMETER :: c_char_space = 4
  INTEGER, PARAMETER :: c_char_opcode = 5
  INTEGER, PARAMETER :: c_char_unknown = 1024

  INTEGER, PARAMETER :: c_prc_not_this_type = 0

  ! block type constants
  INTEGER, PARAMETER :: c_pt_variable = 1
  INTEGER, PARAMETER :: c_pt_constant = 2
  INTEGER, PARAMETER :: c_pt_operator = 3
  INTEGER, PARAMETER :: c_pt_function = 4
  INTEGER, PARAMETER :: c_pt_parenthesis = 5
  INTEGER, PARAMETER :: c_pt_separator = 6
  INTEGER, PARAMETER :: c_pt_character = 7
  INTEGER, PARAMETER :: c_pt_deck_constant = 8
  INTEGER, PARAMETER :: c_pt_species = 9
  INTEGER, PARAMETER :: c_pt_subset = 10
  INTEGER, PARAMETER :: c_pt_default_constant = 11
  INTEGER, PARAMETER :: c_pt_max = 11
  INTEGER, PARAMETER :: c_pt_bad = 1024
  INTEGER, PARAMETER :: c_pt_null = 1025

  CHARACTER(LEN=*), PARAMETER :: c_pt_char(c_pt_max) = (/ &
      'variable        ', &
      'constant        ', &
      'operator        ', &
      'function        ', &
      'parenthesis     ', &
      'separator       ', &
      'character       ', &
      'deck_constant   ', &
      'species         ', &
      'subset          ', &
      'default_constant'/)

  ! Opcode constants
  INTEGER, PARAMETER :: c_opcode_plus = 1
  INTEGER, PARAMETER :: c_opcode_minus = 2
  INTEGER, PARAMETER :: c_opcode_times = 3
  INTEGER, PARAMETER :: c_opcode_divide = 4
  INTEGER, PARAMETER :: c_opcode_power = 5
  INTEGER, PARAMETER :: c_opcode_expo = 6
  INTEGER, PARAMETER :: c_opcode_lt = 7
  INTEGER, PARAMETER :: c_opcode_gt = 8
  INTEGER, PARAMETER :: c_opcode_eq = 9
  INTEGER, PARAMETER :: c_opcode_and = 10
  INTEGER, PARAMETER :: c_opcode_or = 11
  INTEGER, PARAMETER :: c_opcode_unary_minus = 12
  INTEGER, PARAMETER :: c_opcode_unary_plus = 13
  INTEGER, PARAMETER :: c_num_ops = 13

  ! Associativity constants
  INTEGER, PARAMETER :: c_assoc_a = 1
  INTEGER, PARAMETER :: c_assoc_la = 2
  INTEGER, PARAMETER :: c_assoc_ra = 3

  INTEGER, DIMENSION(c_num_ops), PARAMETER :: &
      opcode_precedence = (/2, 2, 3, 3, 4, 4, 1, 1, 1, 0, 0, 4, 4/)
  INTEGER, DIMENSION(c_num_ops), PARAMETER :: &
      opcode_assoc = (/c_assoc_a, c_assoc_la, c_assoc_a, c_assoc_la, &
          c_assoc_ra, c_assoc_ra, c_assoc_la, c_assoc_la, c_assoc_la, &
          c_assoc_la, c_assoc_la, c_assoc_ra, c_assoc_ra/)

  INTEGER, PARAMETER :: c_paren_left_bracket = 1
  INTEGER, PARAMETER :: c_paren_right_bracket = 2

  ! Actual constants
  INTEGER, PARAMETER :: c_const_pi = 1
  INTEGER, PARAMETER :: c_const_kb = 2
  INTEGER, PARAMETER :: c_const_me = 3
  INTEGER, PARAMETER :: c_const_qe = 4
  INTEGER, PARAMETER :: c_const_eps0 = 5
  INTEGER, PARAMETER :: c_const_mu0 = 6
  INTEGER, PARAMETER :: c_const_c = 7
  INTEGER, PARAMETER :: c_const_ev = 8
  INTEGER, PARAMETER :: c_const_kev = 9
  INTEGER, PARAMETER :: c_const_mev = 10
  INTEGER, PARAMETER :: c_const_milli = 11
  INTEGER, PARAMETER :: c_const_micro = 12
  INTEGER, PARAMETER :: c_const_nano = 13
  INTEGER, PARAMETER :: c_const_pico = 14
  INTEGER, PARAMETER :: c_const_femto = 15
  INTEGER, PARAMETER :: c_const_atto = 16
  INTEGER, PARAMETER :: c_const_true = 17
  INTEGER, PARAMETER :: c_const_false = 18

  ! Constants refering to grid properties
  INTEGER, PARAMETER :: c_const_xb = 22
  INTEGER, PARAMETER :: c_const_yb = 23
  INTEGER, PARAMETER :: c_const_zb = 24
  INTEGER, PARAMETER :: c_const_x = 25
  INTEGER, PARAMETER :: c_const_y = 26
  INTEGER, PARAMETER :: c_const_z = 27
  INTEGER, PARAMETER :: c_const_lx = 28
  INTEGER, PARAMETER :: c_const_ly = 29
  INTEGER, PARAMETER :: c_const_lz = 30
  INTEGER, PARAMETER :: c_const_dx = 31
  INTEGER, PARAMETER :: c_const_dy = 32
  INTEGER, PARAMETER :: c_const_dz = 33
  INTEGER, PARAMETER :: c_const_x_min = 34
  INTEGER, PARAMETER :: c_const_y_min = 35
  INTEGER, PARAMETER :: c_const_z_min = 36
  INTEGER, PARAMETER :: c_const_x_max = 37
  INTEGER, PARAMETER :: c_const_y_max = 38
  INTEGER, PARAMETER :: c_const_z_max = 39
  INTEGER, PARAMETER :: c_const_ix = 40
  INTEGER, PARAMETER :: c_const_iy = 41
  INTEGER, PARAMETER :: c_const_iz = 42
  INTEGER, PARAMETER :: c_const_nx = 43
  INTEGER, PARAMETER :: c_const_ny = 44
  INTEGER, PARAMETER :: c_const_nz = 45
  INTEGER, PARAMETER :: c_const_time = 46
  INTEGER, PARAMETER :: c_const_r_xy = 47
  INTEGER, PARAMETER :: c_const_r_yz = 48
  INTEGER, PARAMETER :: c_const_r_xz = 49
  INTEGER, PARAMETER :: c_const_nprocx = 50
  INTEGER, PARAMETER :: c_const_nprocy = 51
  INTEGER, PARAMETER :: c_const_nprocz = 52
  INTEGER, PARAMETER :: c_const_nsteps = 53
  INTEGER, PARAMETER :: c_const_t_end = 54
  INTEGER, PARAMETER :: c_const_ndims = 55
  INTEGER, PARAMETER :: c_const_r_xyz = 56

  INTEGER, PARAMETER :: c_const_io_never = 60
  INTEGER, PARAMETER :: c_const_io_always = 61
  INTEGER, PARAMETER :: c_const_io_full = 62
  INTEGER, PARAMETER :: c_const_io_restartable = 63
  INTEGER, PARAMETER :: c_const_io_species = 64
  INTEGER, PARAMETER :: c_const_io_no_sum = 65
  INTEGER, PARAMETER :: c_const_io_average = 66
  INTEGER, PARAMETER :: c_const_io_snapshot = 67
  INTEGER, PARAMETER :: c_const_io_dump_single = 68
  INTEGER, PARAMETER :: c_const_io_average_single = 69

  INTEGER, PARAMETER :: c_const_dir_x = 80
  INTEGER, PARAMETER :: c_const_dir_y = 81
  INTEGER, PARAMETER :: c_const_dir_z = 82
  INTEGER, PARAMETER :: c_const_dir_px = 83
  INTEGER, PARAMETER :: c_const_dir_py = 84
  INTEGER, PARAMETER :: c_const_dir_pz = 85
  INTEGER, PARAMETER :: c_const_dir_en = 86
  INTEGER, PARAMETER :: c_const_dir_gamma_m1 = 87
  INTEGER, PARAMETER :: c_const_dir_xy_angle = 88
  INTEGER, PARAMETER :: c_const_dir_yz_angle = 89
  INTEGER, PARAMETER :: c_const_dir_zx_angle = 90
  INTEGER, PARAMETER :: c_const_dir_mod_p = 91

  INTEGER, PARAMETER :: c_const_maxwell_solver_yee = 100
  INTEGER, PARAMETER :: c_const_maxwell_solver_lehe = 101
  INTEGER, PARAMETER :: c_const_maxwell_solver_lehe_x = 102
  INTEGER, PARAMETER :: c_const_maxwell_solver_lehe_y = 103
  INTEGER, PARAMETER :: c_const_maxwell_solver_lehe_z = 104
  INTEGER, PARAMETER :: c_const_maxwell_solver_cowan = 105
  INTEGER, PARAMETER :: c_const_maxwell_solver_pukhov = 106
  INTEGER, PARAMETER :: c_const_maxwell_solver_custom = 107

  INTEGER, PARAMETER :: c_const_px = 108
  INTEGER, PARAMETER :: c_const_py = 109
  INTEGER, PARAMETER :: c_const_pz = 110

  ! Custom constants
  INTEGER, PARAMETER :: c_const_deck_lowbound = 4096
  INTEGER, PARAMETER :: c_const_custom_lowbound = 8192

  INTEGER, PARAMETER :: c_func_sine = 1
  INTEGER, PARAMETER :: c_func_cosine = 2
  INTEGER, PARAMETER :: c_func_tan = 3
  INTEGER, PARAMETER :: c_func_exp = 4
  INTEGER, PARAMETER :: c_func_arcsine = 10
  INTEGER, PARAMETER :: c_func_arccosine = 11
  INTEGER, PARAMETER :: c_func_arctan = 12
  INTEGER, PARAMETER :: c_func_neg = 13
  INTEGER, PARAMETER :: c_func_if = 14
  INTEGER, PARAMETER :: c_func_floor = 15
  INTEGER, PARAMETER :: c_func_ceil = 16
  INTEGER, PARAMETER :: c_func_nint = 17
  INTEGER, PARAMETER :: c_func_rho = 18
  INTEGER, PARAMETER :: c_func_tempx = 19
  INTEGER, PARAMETER :: c_func_tempy = 20
  INTEGER, PARAMETER :: c_func_tempz = 21
  INTEGER, PARAMETER :: c_func_interpolate = 22
  INTEGER, PARAMETER :: c_func_tanh = 23
  INTEGER, PARAMETER :: c_func_sinh = 24
  INTEGER, PARAMETER :: c_func_cosh = 25
  INTEGER, PARAMETER :: c_func_ex = 26
  INTEGER, PARAMETER :: c_func_ey = 27
  INTEGER, PARAMETER :: c_func_ez = 28
  INTEGER, PARAMETER :: c_func_bx = 29
  INTEGER, PARAMETER :: c_func_by = 30
  INTEGER, PARAMETER :: c_func_bz = 31
  INTEGER, PARAMETER :: c_func_sqrt = 32
  INTEGER, PARAMETER :: c_func_gauss = 33
  INTEGER, PARAMETER :: c_func_semigauss = 34
  INTEGER, PARAMETER :: c_func_crit = 35
  INTEGER, PARAMETER :: c_func_abs = 36
  INTEGER, PARAMETER :: c_func_loge = 37
  INTEGER, PARAMETER :: c_func_log10 = 38
  INTEGER, PARAMETER :: c_func_log_base = 39
  INTEGER, PARAMETER :: c_func_supergauss = 40
  INTEGER, PARAMETER :: c_func_tempx_ev = 41
  INTEGER, PARAMETER :: c_func_tempy_ev = 42
  INTEGER, PARAMETER :: c_func_tempz_ev = 43
  INTEGER, PARAMETER :: c_func_driftx = 44
  INTEGER, PARAMETER :: c_func_drifty = 45
  INTEGER, PARAMETER :: c_func_driftz = 46
  INTEGER, PARAMETER :: c_func_arctan2 = 47
  INTEGER, PARAMETER :: c_func_min = 48
  INTEGER, PARAMETER :: c_func_max = 49

  INTEGER, PARAMETER :: c_func_custom_lowbound = 4096

  ! The order for the spline interpolation used as a particle representation.
  ! png is the number of ghost cells needed by the particles
#ifdef PARTICLE_SHAPE_BSPLINE3
  INTEGER, PARAMETER :: sf_min = -2
  INTEGER, PARAMETER :: sf_max =  2
  INTEGER, PARAMETER :: png =  4
#elif  PARTICLE_SHAPE_TOPHAT
  INTEGER, PARAMETER :: sf_min =  0
  INTEGER, PARAMETER :: sf_max =  1
  INTEGER, PARAMETER :: png =  2
#else
  INTEGER, PARAMETER :: sf_min = -1
  INTEGER, PARAMETER :: sf_max =  1
  INTEGER, PARAMETER :: png =  3
#endif
  ! ng is the number of ghost cells allocated in the arrays
  ! jng is the number of ghost cells needed by the current arrays
  ! In shared_data:
  ! fng is the number of ghost cells needed by the field solver
  ! sng is the number of ghost cells needed by the current smoother
  INTEGER, PARAMETER :: ng = png + 2
  INTEGER, PARAMETER :: jng = MAX(ng,png)
  ! ncell_min is the number of cells needed by the domain before subcyling
  ! of communications is required
  INTEGER, PARAMETER :: ncell_min = (png + 1) / 2 + 1

  INTEGER, PARAMETER :: c_max_zeros = 9
  INTEGER, PARAMETER :: c_dump_part_grid         = 1
  INTEGER, PARAMETER :: c_dump_grid              = 2
  INTEGER, PARAMETER :: c_dump_part_species      = 3
  INTEGER, PARAMETER :: c_dump_part_weight       = 4
  INTEGER, PARAMETER :: c_dump_part_px           = 5
  INTEGER, PARAMETER :: c_dump_part_py           = 6
  INTEGER, PARAMETER :: c_dump_part_pz           = 7
  INTEGER, PARAMETER :: c_dump_part_vx           = 8
  INTEGER, PARAMETER :: c_dump_part_vy           = 9
  INTEGER, PARAMETER :: c_dump_part_vz           = 10
  INTEGER, PARAMETER :: c_dump_part_charge       = 11
  INTEGER, PARAMETER :: c_dump_part_mass         = 12
  INTEGER, PARAMETER :: c_dump_part_id           = 13
  INTEGER, PARAMETER :: c_dump_ex                = 14
  INTEGER, PARAMETER :: c_dump_ey                = 15
  INTEGER, PARAMETER :: c_dump_ez                = 16
  INTEGER, PARAMETER :: c_dump_bx                = 17
  INTEGER, PARAMETER :: c_dump_by                = 18
  INTEGER, PARAMETER :: c_dump_bz                = 19
  INTEGER, PARAMETER :: c_dump_jx                = 20
  INTEGER, PARAMETER :: c_dump_jy                = 21
  INTEGER, PARAMETER :: c_dump_jz                = 22
  INTEGER, PARAMETER :: c_dump_ekbar             = 23
  INTEGER, PARAMETER :: c_dump_mass_density      = 24
  INTEGER, PARAMETER :: c_dump_charge_density    = 25
  INTEGER, PARAMETER :: c_dump_number_density    = 26
  INTEGER, PARAMETER :: c_dump_temperature       = 27
  INTEGER, PARAMETER :: c_dump_dist_fns          = 28
  INTEGER, PARAMETER :: c_dump_probes            = 29
  INTEGER, PARAMETER :: c_dump_ejected_particles = 30
  INTEGER, PARAMETER :: c_dump_ekflux            = 31
  INTEGER, PARAMETER :: c_dump_poynt_flux        = 32
  INTEGER, PARAMETER :: c_dump_cpml_psi_eyx      = 33
  INTEGER, PARAMETER :: c_dump_cpml_psi_ezx      = 34
  INTEGER, PARAMETER :: c_dump_cpml_psi_byx      = 35
  INTEGER, PARAMETER :: c_dump_cpml_psi_bzx      = 36
  INTEGER, PARAMETER :: c_dump_cpml_psi_exy      = 37
  INTEGER, PARAMETER :: c_dump_cpml_psi_ezy      = 38
  INTEGER, PARAMETER :: c_dump_cpml_psi_bxy      = 39
  INTEGER, PARAMETER :: c_dump_cpml_psi_bzy      = 40
  INTEGER, PARAMETER :: c_dump_cpml_psi_exz      = 41
  INTEGER, PARAMETER :: c_dump_cpml_psi_eyz      = 42
  INTEGER, PARAMETER :: c_dump_cpml_psi_bxz      = 43
  INTEGER, PARAMETER :: c_dump_cpml_psi_byz      = 44
  INTEGER, PARAMETER :: c_dump_absorption        = 45
  INTEGER, PARAMETER :: c_dump_part_ek           = 46
  INTEGER, PARAMETER :: c_dump_part_opdepth      = 47
  INTEGER, PARAMETER :: c_dump_part_qed_energy   = 48
  INTEGER, PARAMETER :: c_dump_part_opdepth_tri  = 49
  INTEGER, PARAMETER :: c_dump_total_energy      = 50
  INTEGER, PARAMETER :: c_dump_total_energy_sum  = 51
  INTEGER, PARAMETER :: c_dump_part_rel_mass     = 52
  INTEGER, PARAMETER :: c_dump_part_gamma        = 53
  INTEGER, PARAMETER :: c_dump_part_proc         = 54
  INTEGER, PARAMETER :: c_dump_part_proc0        = 55
  INTEGER, PARAMETER :: c_dump_ppc               = 56
  INTEGER, PARAMETER :: c_dump_average_weight    = 57
  INTEGER, PARAMETER :: c_dump_persistent_ids    = 58
  INTEGER, PARAMETER :: c_dump_temperature_x     = 59
  INTEGER, PARAMETER :: c_dump_temperature_y     = 60
  INTEGER, PARAMETER :: c_dump_temperature_z     = 61
  INTEGER, PARAMETER :: c_dump_average_px        = 62
  INTEGER, PARAMETER :: c_dump_average_py        = 63
  INTEGER, PARAMETER :: c_dump_average_pz        = 64
  INTEGER, PARAMETER :: c_dump_part_work_x       = 65
  INTEGER, PARAMETER :: c_dump_part_work_y       = 66
  INTEGER, PARAMETER :: c_dump_part_work_z       = 67
  INTEGER, PARAMETER :: c_dump_part_work_x_total = 68
  INTEGER, PARAMETER :: c_dump_part_work_y_total = 69
  INTEGER, PARAMETER :: c_dump_part_work_z_total = 70
  INTEGER, PARAMETER :: c_dump_part_opdepth_brem = 71
  INTEGER, PARAMETER :: c_dump_probe_time        = 72
  INTEGER, PARAMETER :: c_dump_cou_log           = 73
  INTEGER, PARAMETER :: c_dump_part_rate_fi      = 74
  INTEGER, PARAMETER :: c_dump_part_rate_ci      = 75
  INTEGER, PARAMETER :: c_dump_part_rate_dr      = 76
  INTEGER, PARAMETER :: c_dump_part_rate_rr      = 77
  INTEGER, PARAMETER :: c_dump_part_rate_3br     = 78
  INTEGER, PARAMETER :: num_vars_to_dump         = 78

  INTEGER, PARAMETER :: c_subset_random     = 1
  INTEGER, PARAMETER :: c_subset_gamma_min  = 2
  INTEGER, PARAMETER :: c_subset_gamma_max  = 3
  INTEGER, PARAMETER :: c_subset_x_min      = 4
  INTEGER, PARAMETER :: c_subset_x_max      = 5
  INTEGER, PARAMETER :: c_subset_y_min      = 6
  INTEGER, PARAMETER :: c_subset_y_max      = 7
  INTEGER, PARAMETER :: c_subset_z_min      = 8
  INTEGER, PARAMETER :: c_subset_z_max      = 9
  INTEGER, PARAMETER :: c_subset_px_min     = 10
  INTEGER, PARAMETER :: c_subset_px_max     = 11
  INTEGER, PARAMETER :: c_subset_py_min     = 12
  INTEGER, PARAMETER :: c_subset_py_max     = 13
  INTEGER, PARAMETER :: c_subset_pz_min     = 14
  INTEGER, PARAMETER :: c_subset_pz_max     = 15
  INTEGER, PARAMETER :: c_subset_weight_min = 16
  INTEGER, PARAMETER :: c_subset_weight_max = 17
  INTEGER, PARAMETER :: c_subset_charge_min = 18
  INTEGER, PARAMETER :: c_subset_charge_max = 19
  INTEGER, PARAMETER :: c_subset_mass_min   = 20
  INTEGER, PARAMETER :: c_subset_mass_max   = 21
  INTEGER, PARAMETER :: c_subset_id_min     = 22
  INTEGER, PARAMETER :: c_subset_id_max     = 23
  INTEGER, PARAMETER :: c_subset_max        = 23

  INTEGER, PARAMETER :: c_max_string_length = 64
  INTEGER, PARAMETER :: c_max_prefix = 16
  ! Maximum path length on Linux machines
  INTEGER, PARAMETER :: c_max_path_length = 4096 + c_max_prefix

END MODULE constants
