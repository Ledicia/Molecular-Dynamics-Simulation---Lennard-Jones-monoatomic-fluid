!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! STATS_MATH MODULE (pure numerical utilities for scalar statistics)
!
! Goal:
! - Provide a small collection of STATISTICS / NUMERICAL ROUTINES with NO internal state.
! - These routines do NOT store samples, do NOT allocate internal series, and do NOT depend on MD globals.
! - They are meant to be used by higher-level accumulator modules:
!     * md_means                (time averages)
!     * md_correlations         (time-series storage + autocorrelation)
!     * md_correlation_means    (block-averaging of correlation curves)
!
! What this module provides:
! - Robust standard deviation from <x^2> and <x>
! - Autocorrelation of one scalar series A(k) for k = 1...N: C(L) = (1/(N-L)) * sum_{k=1..N-L} A(k) * A(k+L)
! - Normalized autocorrelation: C_norm(L) = C(L) / C(0) (with protection if C(0) ~ 0)
!
! Design notes:
! - Everything is PRIVATE by default to avoid name conflicts.
! - Routines are written to be:
!     * deterministic
!     * independent of global state
!     * safe against common indexing/size errors
! - Whenever possible, routines are PURE to enable compiler optimizations and easier testing.
!******************************************************************************************************
module stats_math

  use define_precision, only: dp_kind, int_kind
  implicit none

  private

  ! Public API
  public :: stats_std_from_moments
  public :: autocorr_scalar, autocorr_scalar_centered
  public :: normalize_corr

contains

  !====================================================================================================
  ! Robust standard deviation from the first two moments:
  !
  ! Inputs:
  !   mean_x2 : <x^2>
  !   mean_x  : <x>
  !
  ! Output:
  !   std = sqrt( <x^2> - <x>^2 )
  !
  ! Notes:
  ! - Protects against tiny negative values due to round-off:
  !     <x^2> - <x>^2 can become slightly negative (~ -1e-16).
  !====================================================================================================
  pure function stats_std_from_moments(mean_x2, mean_x) result(std)
    real(kind=dp_kind), intent(in) :: mean_x2, mean_x
    real(kind=dp_kind) :: std

    std = dsqrt(max(0.d0, mean_x2 - mean_x*mean_x))
  end function stats_std_from_moments


  !====================================================================================================
  ! Autocorrelation of a scalar series using dot_product slices (O(N*lag_max)).
  !
  ! Series:
  !   A(k), k = 1...N <=> A: scalar quantity sampled during the MD simulation (E.g. E_pot(t), E_kin(t), E_tot(t))
  !
  ! Raw autocorrelation:
  !   C(L) = (1/(N-L)) * sum_{k=1...N-L} A(k) * A(k+L)
  !
  ! Inputs:
  !   series_in    : array containing at least n_samples_in values
  !   n_samples_in : number of valid samples N to use (N >= 1)
  !   lag_max      : maximum lag L (0..lag_max), must satisfy lag_max < n_samples_in
  !
  ! Output:
  !   corr_out(0:) : array where corr_out(L) is C(L) for L = 0..lag_max
  !
  ! Notes:
  ! - This routine does NOT normalize; call normalize_corr() for C(L)/C(0).
  ! - Caller must provide corr_out with ubound >= lag_max.
  !====================================================================================================
  subroutine autocorr_scalar(series_in, n_samples_in, lag_max, corr_out)
    integer(kind=int_kind), intent(in)  :: n_samples_in, lag_max
    real(kind=dp_kind),     intent(in)  :: series_in(:)
    real(kind=dp_kind),     intent(out) :: corr_out(0:)

    integer(kind=int_kind) :: lag
    integer(kind=int_kind) :: n_valid

    if (n_samples_in <= 0_int_kind) stop 'autocorr_scalar(): n_samples must be > 0.'
    if (lag_max < 0_int_kind) stop 'autocorr_scalar(): lag_max must be >= 0.'
    if (lag_max >= n_samples_in) stop 'autocorr_scalar(): lag_max must be < n_samples.'
    if (size(series_in) < n_samples_in) stop 'autocorr_scalar(): series size < n_samples.'
    if (ubound(corr_out,1) < lag_max) stop 'autocorr_scalar(): corr_out upper bound < lag_max.'

    do lag = 0_int_kind, lag_max
      n_valid = n_samples_in - lag
      corr_out(lag) = dot_product(series_in(1:n_valid), series_in(1+lag : lag+n_valid)) / dble(n_valid)
    end do
  end subroutine autocorr_scalar
  
  
  !====================================================================================================
  ! Centered autocorrelation (autocovariance) of a scalar series using dot_product slices (O(N*lag_max)).
  !
  ! Series: A(k), k = 1...N <=> A: scalar quantity sampled during the MD simulation (E.g. E_pot(t), E_kin(t), E_tot(t), T(t), P(t))
  !
  ! Centered autocorrelation (autocovariance):
  !   Let <A> = (1/N) * sum_{k=1...N} A(k)  =>  C_c(L) = (1/(N-L)) * sum_{k=1...N-L} [ (A(k) - <A>) * (A(k+L) - <A>) ]
  !
  ! Inputs:
  !   series_in    : array containing at least n_samples_in values
  !   n_samples_in : number of valid samples N to use (N >= 1)
  !   lag_max      : maximum lag L (0..lag_max), must satisfy lag_max < n_samples_in
  !
  ! Output:
  !   corr_out(0:) : array where corr_out(L) is C_c(L) for L = 0..lag_max
  !
  ! Notes:
  ! - This routine does NOT normalize; call normalize_corr() for C_c(L)/C_c(0).
  ! - Caller must provide corr_out with ubound >= lag_max.
  ! - The series is centered using the MEAN computed from the same N samples.
  !====================================================================================================
  subroutine autocorr_scalar_centered(series_in, n_samples_in, lag_max, corr_out)
    integer(kind=int_kind), intent(in)  :: n_samples_in, lag_max
    real(kind=dp_kind),     intent(in)  :: series_in(:)
    real(kind=dp_kind),     intent(out) :: corr_out(0:)

    integer(kind=int_kind) :: lag
    integer(kind=int_kind) :: n_valid
    real(kind=dp_kind)     :: mean_a

    if (n_samples_in <= 0_int_kind) stop 'autocorr_scalar_centered(): n_samples must be > 0.'
    if (lag_max < 0_int_kind) stop 'autocorr_scalar_centered(): lag_max must be >= 0.'
    if (lag_max >= n_samples_in) stop 'autocorr_scalar_centered(): lag_max must be < n_samples.'
    if (size(series_in) < n_samples_in) stop 'autocorr_scalar_centered(): series size < n_samples.'
    if (ubound(corr_out,1) < lag_max) stop 'autocorr_scalar_centered(): corr_out upper bound < lag_max.'

    mean_a = sum(series_in(1:n_samples_in)) / dble(n_samples_in)

    do lag = 0_int_kind, lag_max
      n_valid = n_samples_in - lag
      corr_out(lag) = dot_product( series_in(1:n_valid) - mean_a, &
                                   series_in(1+lag : lag+n_valid) - mean_a ) / dble(n_valid)
    end do
  end subroutine autocorr_scalar_centered


  !====================================================================================================
  ! Normalize a correlation curve by its value at lag = 0:
  !
  !   C_norm(L) = C(L) / C(0)
  !
  ! Inputs:
  !   lag_max    : maximum lag index used (0..lag_max)
  !   corr_in(0:): raw correlation curve C(L)
  !
  ! Output:
  !   corr_out(0:): normalized correlation curve C_norm(L)
  !
  ! Notes:
  ! - If |C(0)| is extremely small, normalization is ill-defined -> return zeros.
  ! - Caller must provide corr_out with ubound >= lag_max.
  !====================================================================================================
  subroutine normalize_corr(lag_max, corr_in, corr_out)
    integer(kind=int_kind), intent(in)  :: lag_max
    real(kind=dp_kind),     intent(in)  :: corr_in(0:)
    real(kind=dp_kind),     intent(out) :: corr_out(0:)

    real(kind=dp_kind) :: c0
    integer(kind=int_kind) :: lag

    if (lag_max < 0_int_kind) stop 'normalize_corr(): lag_max must be >= 0.'
    if (ubound(corr_in,1)  < lag_max) stop 'normalize_corr(): corr_in upper bound < lag_max.'
    if (ubound(corr_out,1) < lag_max) stop 'normalize_corr(): corr_out upper bound < lag_max.'

    c0 = corr_in(0)

    if (abs(c0) <= 1.d-14) then
      corr_out(0:lag_max) = 0.d0
      return
    end if

    do lag = 0_int_kind, lag_max
      corr_out(lag) = corr_in(lag) / c0
    end do
  end subroutine normalize_corr

end module stats_math
