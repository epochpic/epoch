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

MODULE random_generator

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: random, random_init, get_random_state, set_random_state
  PUBLIC :: random_box_muller, random_poisson, random_state_type, random_flush_cache

  INTEGER, PARAMETER :: init_x = 123456789
  INTEGER, PARAMETER :: init_y = 362436069
  INTEGER, PARAMETER :: init_z = 521288629
  INTEGER, PARAMETER :: init_w = 916191069

  TYPE :: random_state_type
    INTEGER :: x, y, z, w
    LOGICAL :: box_muller_cached
    DOUBLE PRECISION :: cached_random_value
  END TYPE random_state_type

  TYPE(random_state_type), TARGET, SAVE :: global_random

CONTAINS

  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators,
  !     period 597273182964842497>2^59
  ! Overall period>2^123;  Default seeds x,y,z,w.
  ! Set your own seeds with random_init(ix,iy,iz,iw).
  FUNCTION random(state)

    DOUBLE PRECISION :: random
    TYPE(random_state_type), INTENT(INOUT), OPTIONAL, TARGET :: state
    INTEGER :: kiss, a1, b1, a2, b2, a3, b3
    TYPE(random_state_type), POINTER :: current_state

    IF (PRESENT(state)) THEN
      current_state => state
    ELSE
      current_state => global_random
    END IF

    current_state%x = 69069 * current_state%x + 1327217885

    a1 = current_state%y
    b1 = 13
    a2 = IEOR(a1, ISHFT(a1, b1))
    b2 = -17
    a3 = IEOR(a2, ISHFT(a2, b2))
    b3 = 5
    current_state%y  = IEOR(a3, ISHFT(a3, b3))

    current_state%z = 18000 * IAND(current_state%z, 65535) &
        + ISHFT(current_state%z, - 16)
    current_state%w = 30903 * IAND(current_state%w, 65535) &
        + ISHFT(current_state%w, - 16)

    kiss = current_state%x + current_state%y + ISHFT(current_state%z, 16) &
        + current_state%w

    random = (DBLE(kiss) + 2147483648.d0) / 4294967296.d0

  END FUNCTION random



  SUBROUTINE random_init(seed, state)

    INTEGER, INTENT(IN) :: seed
    TYPE(random_state_type), INTENT(INOUT), OPTIONAL, TARGET :: state
    TYPE(random_state_type), POINTER :: current_state
    INTEGER :: i
    DOUBLE PRECISION :: dummy

    IF (PRESENT(state)) THEN
      current_state => state
    ELSE
      current_state => global_random
    END IF

    current_state%x = init_x + seed
    current_state%y = init_y + seed
    current_state%z = init_z + seed
    current_state%w = init_w + seed
    current_state%box_muller_cached = .FALSE.
    current_state%cached_random_value = 0.0d0

    ! 'Warm-up' the generator by cycling through a few times
    DO i = 1, 1000
      dummy = random(state)
    END DO

  END SUBROUTINE random_init



  FUNCTION random_box_muller(stdev, mu, state)

    DOUBLE PRECISION, INTENT(IN) :: stdev
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: mu
    TYPE(random_state_type), INTENT(INOUT), TARGET, OPTIONAL :: state
    DOUBLE PRECISION :: random_box_muller

    DOUBLE PRECISION :: rand1, rand2, w, mu_val
    DOUBLE PRECISION, PARAMETER :: c_tiny = TINY(1.0D0)
    DOUBLE PRECISION, SAVE :: cached_random_value
    LOGICAL, SAVE :: cached = .FALSE.

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers

    IF (PRESENT(mu)) THEN
      mu_val = mu
    ELSE
      mu_val = 0.0D0
    END IF

    IF (PRESENT(state)) THEN
      cached = state%box_muller_cached
      cached_random_value = state%cached_random_value
    ELSE
      cached = global_random%box_muller_cached
      cached_random_value = global_random%cached_random_value
    END IF

    IF (cached) THEN
      cached = .FALSE.
      random_box_muller = cached_random_value * stdev + mu_val
    ELSE
      cached = .TRUE.

      DO
        rand1 = random(state)
        rand2 = random(state)

        rand1 = 2.0D0 * rand1 - 1.0D0
        rand2 = 2.0D0 * rand2 - 1.0D0

        w = rand1**2 + rand2**2

        IF (w > c_tiny .AND. w < 1.0D0) EXIT
      END DO

      w = SQRT((-2.0D0 * LOG(w)) / w)

      random_box_muller = rand1 * w * stdev + mu_val
      cached_random_value = rand2 * w
    END IF

    IF (PRESENT(state)) THEN
      state%box_muller_cached = cached
      state%cached_random_value = cached_random_value
    ELSE
      global_random%box_muller_cached = cached
      global_random%cached_random_value = cached_random_value
    END IF

  END FUNCTION random_box_muller
  
  
  
  FUNCTION factorial(n)
  ! Quick recursive factorial function
  ! C. Arran 2023
  
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION :: factorial
    INTEGER :: i
    
    DO i = 2,INT(n)
      factorial = factorial * i
    END DO
    
  END FUNCTION factorial



  FUNCTION random_poisson(lambda)
  ! This samples a random number from a Poisson distribution with mean lambda
  ! Based on A.C. Atkinson Appl. Statist. 28, 29 (1979)
  ! C. Arran 2023
  
    DOUBLE PRECISION, INTENT(IN) :: lambda
    DOUBLE PRECISION :: random_poisson
    INTEGER :: n
    DOUBLE PRECISION :: target, total
    DOUBLE PRECISION :: alpha, beta, k, x, u

    IF (lambda.LT.20) THEN
    ! At low lambda < 20 it uses the multiplication method of successive independent Poisson processes
      total = 1
      n = -1
      target = EXP(-lambda)
      DO WHILE (total > target)
        total = total * (1.0 - random())
        n = n + 1
      END DO           
      random_poisson = n 

    ELSE
    ! At high lambda >= 20 it uses a rejection sampling algorithm enveloped by the logistic distribution
    ! The shape of the logistic envelope varies with lambda - estimate this with alpha, beta, k
      beta = 4*ATAN(1.d0) / SQRT(3*lambda) ! pi = 4*ATAN(1.d0)
      alpha = beta*lambda
      k = LOG(0.767-3.36/lambda) - LOG(beta)
      
      total = 1
      target = 0
      DO WHILE (total > target)
        x = 0
        DO WHILE (x < -0.5d0)
          u = random()
          x = (alpha - LOG((1-u)/u)) / beta
        END DO
        n = INT(x+0.5d0)
        u = random()
        total = alpha - beta*x + LOG(u/(1+EXP(alpha-beta*x))**2) - k
        target = - lambda + n*LOG(lambda) - LOG(factorial(n))
      END DO
      random_poisson = n

    END IF

  END FUNCTION random_poisson



  SUBROUTINE get_random_state(state)

    INTEGER, INTENT(OUT) :: state(:)

    state(1) = global_random%x
    state(2) = global_random%y
    state(3) = global_random%z
    state(4) = global_random%w

  END SUBROUTINE get_random_state



  SUBROUTINE set_random_state(state)

    INTEGER, INTENT(IN) :: state(:)

    global_random%x = state(1)
    global_random%y = state(2)
    global_random%z = state(3)
    global_random%w = state(4)

  END SUBROUTINE set_random_state



  SUBROUTINE random_flush_cache(state)

    TYPE(random_state_type), INTENT(INOUT), OPTIONAL :: state

    IF (PRESENT(state)) THEN
      state%box_muller_cached = .FALSE.
      state%cached_random_value = 0.0d0
    ELSE
      global_random%box_muller_cached = .FALSE.
      global_random%cached_random_value = 0.0d0
    END IF

  END SUBROUTINE random_flush_cache

END MODULE random_generator
