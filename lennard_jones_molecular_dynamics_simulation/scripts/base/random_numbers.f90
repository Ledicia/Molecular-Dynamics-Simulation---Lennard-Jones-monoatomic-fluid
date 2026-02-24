!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! RANDOM NUMBER MODULE
! - Provide a simple, reproducible uniform random number generator in (0,1).
! - This generator is deterministic: using the same initial seed produces
!   exactly the same sequence of random numbers.
!
! Typical usage:
! --------------
!   integer(kind=int_kind) :: seed
!   real(kind=dp_kind)     :: r
!
!   seed = -12345        ! initialize (negative seed resets generator)
!   r = random_uniform(seed)
!
! Notes:
! ------
! - If seed <= 0, the generator is (re)initialized.
! - After initialization, seed is internally updated and should NOT be reset.
!******************************************************************************************************

module random_numbers
  use define_precision, only: dp_kind, int_kind
  implicit none
  private

  public :: random_uniform

contains

  ! -----------------------------------------------------------------------------------
  ! Uniform random number generator in (0,1).
  ! - Based on a classic subtractive method.
  ! - Uses an internal table of 55 numbers.
  !
  ! Input / Output:
  !  - seed (inout):
  !     - If seed <= 0, the generator is initialized.
  !     - On output, seed is set to a positive value and should be reused.
  !  - r (out): real number uniformly distributed in (0,1).
  !
  ! IMPORTANT:
  ! - This function is NOT thread-safe (uses SAVE variables).
  ! - Suitable for serial MD simulations.
  ! -----------------------------------------------------------------------------------
  function random_uniform(seed) result(r)

    integer(kind=int_kind), intent(inout) :: seed
    real(kind=dp_kind)                   :: r

    ! Constants for the generator
    real(kind=dp_kind), parameter :: mbig  = 4.0d6
    real(kind=dp_kind), parameter :: mseed = 1618033.d0
    real(kind=dp_kind), parameter :: mz    = 0.d0
    real(kind=dp_kind), parameter :: fac   = 1.d0 / mbig

    ! Internal state (saved between calls)
    integer(kind=int_kind) :: i, k, ii
    integer(kind=int_kind) :: inext, inextp
    real(kind=dp_kind)     :: mj, mk
    real(kind=dp_kind)     :: ma(55)
    integer(kind=int_kind) :: iff

    save ma, inext, inextp, iff
    data iff /0/

    ! ------------------------------------------------------------
    ! Initialization block (runs if seed <= 0 or first call)
    ! ------------------------------------------------------------
    if (seed <= 0 .or. iff == 0) then

      iff = 1
      mj  = abs(mseed - abs(seed))
      mj  = mod(mj, mbig)
      ma(55) = mj
      mk  = 1.d0

      do i = 1, 54
        ii = mod(21*i, 55)
        ma(ii) = mk
        mk = mj - mk
        if (mk < mz) mk = mk + mbig
        mj = ma(ii)
      end do

      do k = 1, 4
        do i = 1, 55
          ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
          if (ma(i) < mz) ma(i) = ma(i) + mbig
        end do
      end do

      inext  = 0
      inextp = 31
      seed   = 1

    end if

    ! ------------------------------------------------------------
    ! Generate next random number
    ! ------------------------------------------------------------
    inext = inext + 1
    if (inext == 56) inext = 1

    inextp = inextp + 1
    if (inextp == 56) inextp = 1

    mj = ma(inext) - ma(inextp)
    if (mj < mz) mj = mj + mbig
    ma(inext) = mj

    r = mj * fac

  end function random_uniform

end module random_numbers
    


