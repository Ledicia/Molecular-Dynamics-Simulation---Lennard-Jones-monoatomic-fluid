!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! GEOMETRY + PERIODIC BOUNDARY CONDITIONS (PBC)
! - Geometry utilities for a cubic simulation box with periodic boundary conditions.
! - PBC mimic an infinite system by tiling the simulation box in all directions.
!
! Key ideas:
! - Wrapping: keep coordinates inside [0, L) for clean trajectories and numerical stability.
!   - With PBC, the simulation box is replicated infinitely in all directions.
!     A particle leaving the box re-enters from the opposite face.
! - Minimum-image convention: use the nearest periodic image for pair distances
!   - For each particle pair (i,j), the closest periodic image of j relative to i is used.
!     This is implemented by shifting dx,dy,dz by +/-L so that each component lies in [-L/2, L/2].
!
! Notes:
! - Wrapping positions is not strictly required for correct forces (if minimum-image is used),
!   but it prevents coordinates from drifting to large values over long runs.
! - This module does not track unwrapped trajectories or image counters (needed for MSD).
!******************************************************************************************************
module geometry_pbc
  use define_precision, only: dp_kind, int_kind
  implicit none

contains

  ! ------------------------------------------------------------
  ! Wrap particle positions into the simulation box [0, L)
  !  - Maps any coordinate x to x = x - L * floor(x/L)
  !  - Guarantees x in [0, L) for any real x (including negative values)
  !
  ! 1D examples (L = 10):
  !  - x =  9.8  ->  9.8
  !  - x = 10.2  ->  0.2
  !  - x = -0.4  ->  9.6
  !  - x = 20.7  ->  0.7
  ! ------------------------------------------------------------
  subroutine wrap_positions(rx, ry, rz, box_length)
    real(kind=dp_kind), intent(inout) :: rx(:), ry(:), rz(:)    ! Particle coordinates to be wrapped
    real(kind=dp_kind), intent(in)    :: box_length             ! Box length L

    integer(kind=int_kind) :: i                                 ! Particle index
    integer(kind=int_kind) :: n                                 ! Number of particles (array size)
    real(kind=dp_kind)     :: invL                              ! 1/L (used to avoid repeated division)

    if (box_length <= 0.d0) stop 'wrap_positions(): box_length must be > 0.'
    if (size(ry) /= size(rx) .or. size(rz) /= size(rx)) stop 'wrap_positions(): inconsistent array sizes.'

    n    = size(rx)                                              ! Number of particles in the arrays
    invL = 1.d0 / box_length                                     ! Precompute inverse box length

    do i = 1, n
      rx(i) = rx(i) - box_length * floor(rx(i) * invL)           ! Wrap x into [0, L)
      ry(i) = ry(i) - box_length * floor(ry(i) * invL)           ! Wrap y into [0, L)
      rz(i) = rz(i) - box_length * floor(rz(i) * invL)           ! Wrap z into [0, L)
    end do

  end subroutine wrap_positions


  ! ------------------------------------------------------------
  ! Minimum-image convention for a displacement component
  !  - For a displacement dx = x_i - x_j between two particles, periodicity implies dx is equivalent to:
  !      dx, dx ± L, dx ± 2L, ...
  !  - Choose the nearest image so that dx lies in [-L/2, L/2]:
  !      dx_mic = dx - L * dnint(dx/L)
  !
  ! Example (L = 10):
  !  - x_i = 9.8,  x_j = 0.3 where x_j is the wrapped position, one possible unwrapped position would be xu_j = 10.5
  !    Raw displacement: dx = x_i - x_j = 9.5 (looks far)
  !    Apply MIC: dx_mic = dx - L*nint(dx/L) = 9.5 - 10*nint(0.95) = 9.5 - 10*1 = -0.5 (nearest image)
  !    Interpretation:
  !     - Particle j is actually 0.5 units to the right of particle i across the boundary.
  !
  ! Notes:
  ! - dnint returns the nearest integer as a real(dp_kind).
  ! - Passing inv_box_length avoids repeated division inside pair loops.
  ! ------------------------------------------------------------
  pure function minimum_image(dx, box_length, inv_box_length) result(dx_mic)
    real(kind=dp_kind), intent(in) :: dx                         ! Raw displacement component
    real(kind=dp_kind), intent(in) :: box_length                 ! Box length L
    real(kind=dp_kind), intent(in) :: inv_box_length             ! 1/L
    real(kind=dp_kind)             :: dx_mic                     ! Displacement mapped to [-L/2, L/2]

    dx_mic = dx - box_length * dnint(dx * inv_box_length)        ! Apply minimum-image convention

  end function minimum_image

end module geometry_pbc
