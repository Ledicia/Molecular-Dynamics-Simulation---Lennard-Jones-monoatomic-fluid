!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! MD_CORRELATIONS MODULE (time-series storage + autocorrelations for scalar observables + Correlation means)
!
! Goal:
! - Store scalar time series sampled during MD (one value per sampling time).
! - Compute raw autocorrelation C(tau) and normalized autocorrelation C(tau)/C(0).
! - Compute correlation means over blocks of time series
!
! 1) Time-Series Storage (Observables stored (reduced units)):
!   - E_pot(t), E_kin(t), E_tot(t) = E_pot + E_kin
!   - T(t) from kinetic energy (if provided): T = 2 E_kin / (3 N)
!   - P(t) (optional, if virial/pressure is available and provided by the caller)
!
! 2) Autocorrelation of a scalar series A(k), k = 1...N where A is any observable
!   - Raw autocorrelation at lag L (discrete): C(L) = < A(k) * A(k+L) >  for k = 1..(N-L)
!   - Normalized autocorrelation: C_norm(L) = C(L) / C(0)
!
! 3) Correlation means:
!   - Given a LONG sampled time series A(1...N) of one (or several) scalar observables,
!     compute a SMOOTHER correlation estimate by:
!     1) Splitting the full series into B contiguous blocks (approximately independent).
!         A_b(tau): of observable A at time index tau INSIDE block b

!     2) Computing ONE autocorrelation curve per block:
!          C_b(tau) = (1/(N_b - tau)) * sum_{k=1..N_b-tau} A_b(k) * A_b(k+tau)

!     3) Averaging the curves over blocks:
!          <C(tau)>      = (1/B) * sum_b C_b(tau)
!          <C_norm(tau)> = (1/B) * sum_b [C_b(tau)/C_b(0)]
!          C_b(tau) = (1/(N_b - tau)) * sum_{k=1..N_b-tau} A_b(k) * A_b(k+tau)
!
! Design notes:
! - This module stores NO global state: everything lives inside type(md_corr_accum).
! - Pattern:
!     1) md_corr_init()
!     2) md_corr_add_sample()
!     3) md_corr_compute()
!     4) md_corr_get_*()
!     5) md_corr_free()
!     6) md_corr_means()
!
!
! Extended:
! - This module also provides block-averaged correlation curves (correlation means)
!   computed FROM the same internally stored time series (no duplicated memory).
!******************************************************************************************************
module md_correlations

  use define_precision, only: dp_kind, int_kind
  use md_types,         only: sim_params
  use stats_math,       only: autocorr_scalar, normalize_corr, autocorr_scalar_centered
  implicit none

  private

  ! Public API (correlations)
  public :: md_corr_accum
  public :: md_corr_init, md_corr_free
  public :: md_corr_add_sample
  public :: md_corr_compute
  public :: md_corr_get_info
  public :: md_corr_get_by_name, md_corr_get_all

  ! Public API (correlation means / block averaging)
  public :: md_corr_cm_init, md_corr_cm_free
  public :: md_corr_cm_compute
  public :: md_corr_cm_get_by_name, md_corr_cm_get_all
  public :: md_corr_cm_num_blocks, md_corr_cm_max_lag

  ! ==================================================================================
  ! Accumulator container (NO module-level state)
  ! ==================================================================================
  type :: md_corr_accum

    ! ==================================================================================
    ! Correlations
    ! ==================================================================================

    ! -----------------------
    ! Run configuration snapshot (copied at init)
    ! -----------------------
    integer(kind=int_kind) :: n      = 0_int_kind     ! N (not strictly needed for correlations, but kept for consistency)
    real(kind=dp_kind)     :: volume = 0.d0           ! V (idem)

    ! -----------------------
    ! Storage / control
    ! -----------------------
    logical :: enabled = .false.                              ! Correlations enabled?
    integer(kind=int_kind) :: n_alloc = 0_int_kind            ! Capacity requested by the user (max number of samples)
    integer(kind=int_kind) :: n_samples_used  = 0_int_kind    ! Samples pushed so far
    integer(kind=int_kind) :: max_lag = 0_int_kind            ! Maximum lag requested (0...max_lag)

    ! -----------------------
    ! Which observables are being stored
    ! -----------------------
    logical :: store_temp  = .false.
    logical :: store_press = .false.

    ! -----------------------
    ! Time series
    ! -----------------------
    real(kind=dp_kind), allocatable :: epot_series(:)   ! E_pot(k), k=1...n_samples_used
    real(kind=dp_kind), allocatable :: ekin_series(:)   ! E_kin(k), k=1...n_samples_used
    real(kind=dp_kind), allocatable :: etot_series(:)   ! E_tot(k), k=1...n_samples_used
    real(kind=dp_kind), allocatable :: temp_series(:)   ! T(k),     k=1...n_samples_used
    real(kind=dp_kind), allocatable :: press_series(:)  ! P(k),     k=1...n_samples_used

    ! -----------------------
    ! Autocorrelation buffers (raw and normalized)
    !  - corr_epot(0:max_lag): E_pot autocorrelation curve (C_epot(0), C_epot(1), ..., C_epot(max_lag)), where E_pot(k), k = 1..N_samples_used:
    !    C_epot(L) = (1 / (N_samples_used - L)) * sum_{k=1}^{N_samples_used-L} [ E_pot(k) * E_pot(k+L) ],   L = 0...max_lag
    ! 
    ! Important:
    ! - n_samples: defines the allocated size of the series arrays.
    ! - n_samples_used: true length of the time series
    !
    ! - Example for N_samples_used and N_samples: 
    !   total_steps     = 10000         |  n_samples = total_steps / output_interval = 100
    !   warmup_steps    = 2000       => |
    !   output_interval = 100           |  n_samples_used = (total_steps - warmup_steps) / output_interval = 80
    ! -----------------------
    real(kind=dp_kind), allocatable :: corr_epot(:), corrn_epot(:)
    real(kind=dp_kind), allocatable :: corr_ekin(:), corrn_ekin(:)
    real(kind=dp_kind), allocatable :: corr_etot(:), corrn_etot(:)
    real(kind=dp_kind), allocatable :: corr_temp(:), corrn_temp(:)     ! optional
    real(kind=dp_kind), allocatable :: corr_press(:), corrn_press(:)   ! optional

    logical :: computed = .false.                      ! correlations already computed?

    ! ==================================================================================
    ! Correlation means (block averaging) state
    ! ==================================================================================

    ! -----------------------
    ! Control flags
    ! -----------------------
    logical :: cm_enabled  = .false.
    logical :: cm_computed = .false.

    ! -----------------------
    ! Block settings
    ! -----------------------
    integer(kind=int_kind) :: cm_n_blocks  = 0_int_kind
    integer(kind=int_kind) :: cm_block_len = 0_int_kind

    ! -----------------------
    ! Accumulators (sum over blocks) so as to compute averages later
    ! Arrays of sums of correlations for each block (b = 0, ..., B)
    !   (C_0(0) + C_1(0) + ... + C_B(0), C_0(1) + C_1(1) + ... + C_B(1), ..., C_0(n_valid) + C_1(n_valid) + ... + C_B(n_valid))
    !
    ! C_b(tau): Correlation for time series (epot, ekin, etot, temp, press) in block b and lag tau
    ! -----------------------
    real(kind=dp_kind), allocatable :: cm_sum_corr_epot(:),  cm_sum_corrn_epot(:)
    real(kind=dp_kind), allocatable :: cm_sum_corr_ekin(:),  cm_sum_corrn_ekin(:)
    real(kind=dp_kind), allocatable :: cm_sum_corr_etot(:),  cm_sum_corrn_etot(:)
    real(kind=dp_kind), allocatable :: cm_sum_corr_temp(:),  cm_sum_corrn_temp(:)
    real(kind=dp_kind), allocatable :: cm_sum_corr_press(:), cm_sum_corrn_press(:)

    ! -----------------------
    ! Final averaged curves (mean over blocks): Array of correlation means (1/B sum_{b=1}^{B} c_b(0), 1/B sum_{b=1}^{B} c_b(1), ..., 1/B sum_{b=1}^{B} c_b(n_valid))
    ! -----------------------
    real(kind=dp_kind), allocatable :: cm_mean_corr_epot(:),  cm_mean_corrn_epot(:)
    real(kind=dp_kind), allocatable :: cm_mean_corr_ekin(:),  cm_mean_corrn_ekin(:)
    real(kind=dp_kind), allocatable :: cm_mean_corr_etot(:),  cm_mean_corrn_etot(:)

    real(kind=dp_kind), allocatable :: cm_mean_corr_temp(:),  cm_mean_corrn_temp(:)
    real(kind=dp_kind), allocatable :: cm_mean_corr_press(:), cm_mean_corrn_press(:)

  end type md_corr_accum

contains

  !====================================================================================================
  ! Allocate buffers for scalar series and their autocorrelations.
  !
  ! Inputs:
  !   params     : used only to store N and V for consistency / logging.
  !   n_samples  : maximum number of samples that will be stored (capacity).
  !   lag_max    : maximum lag to compute (0..lag_max), must satisfy lag_max < n_samples.
  !   want_temp  : if .true., store T series and compute its correlation (caller must pass temp in add_sample()).
  !   want_press : if .true., store P series and compute its correlation (caller must pass press in add_sample()).
  !
  ! Notes:
  ! - Always stores Epot/Ekin/Etot correlations because they are the core scalars.
  ! - Temperature and pressure are optional to avoid overhead if you never use them.
  !====================================================================================================
  subroutine md_corr_init(acc, params, n_samples, lag_max, want_temp, want_press)

    type(md_corr_accum), intent(inout) :: acc
    type(sim_params),    intent(in)    :: params
    integer(kind=int_kind), intent(in) :: n_samples, lag_max
    logical, intent(in) :: want_temp, want_press

    if (n_samples <= 0_int_kind) stop 'md_corr_init(): n_samples must be > 0.'
    if (lag_max < 0_int_kind) stop 'md_corr_init(): lag_max must be >= 0.'
    if (lag_max >= n_samples) stop 'md_corr_init(): lag_max must be < n_samples.'

    ! -----------------------
    ! Time series allocation
    ! -----------------------
    call md_corr_free(acc)  ! Ensure clean slate (deallocates and resets), then re-init below

    acc%enabled         = .true.
    acc%computed        = .false.
    acc%n_alloc         = n_samples
    acc%n_samples_used  = 0_int_kind
    acc%max_lag         = lag_max
    acc%n               = params%n
    acc%volume          = params%volume
    acc%store_temp      = want_temp
    acc%store_press     = want_press

    allocate(acc%epot_series(1:acc%n_alloc)); acc%epot_series = 0.d0
    allocate(acc%ekin_series(1:acc%n_alloc)); acc%ekin_series = 0.d0
    allocate(acc%etot_series(1:acc%n_alloc)); acc%etot_series = 0.d0   

    if (acc%store_temp) then
      allocate(acc%temp_series(1:acc%n_alloc)); acc%temp_series = 0.d0
    end if

    if (acc%store_press) then
      allocate(acc%press_series(1:acc%n_alloc)); acc%press_series = 0.d0
    end if

    ! -----------------------
    ! Correlation buffers allocation => corr_epot, corrn_epot, etc. are vecotrs/arrays (0...max_lag)
    ! -----------------------
    allocate(acc%corr_epot(0:acc%max_lag), acc%corrn_epot(0:acc%max_lag))
    allocate(acc%corr_ekin(0:acc%max_lag), acc%corrn_ekin(0:acc%max_lag))
    allocate(acc%corr_etot(0:acc%max_lag), acc%corrn_etot(0:acc%max_lag))

    acc%corr_epot = 0.d0; acc%corrn_epot = 0.d0
    acc%corr_ekin = 0.d0; acc%corrn_ekin = 0.d0
    acc%corr_etot = 0.d0; acc%corrn_etot = 0.d0

    if (acc%store_temp) then
      allocate(acc%corr_temp(0:acc%max_lag), acc%corrn_temp(0:acc%max_lag))
      acc%corr_temp = 0.d0; acc%corrn_temp = 0.d0
    end if

    if (acc%store_press) then
      allocate(acc%corr_press(0:acc%max_lag), acc%corrn_press(0:acc%max_lag))
      acc%corr_press = 0.d0; acc%corrn_press = 0.d0
    end if

    ! Correlation means state is independent: it stays disabled until md_corr_cm_init()
    acc%cm_enabled  = .false.
    acc%cm_computed = .false.
    acc%cm_n_blocks  = 0_int_kind
    acc%cm_block_len = 0_int_kind

  end subroutine md_corr_init


  !====================================================================================================
  ! Add one sample to the stored scalar series. Correlations are not computed here yet.
  !
  ! Inputs:
  !   epot, ekin : energies at sampling time k
  ! Optional:
  !   temp       : temperature at sampling time k (required if store_temp=.true.)
  !   press      : pressure at sampling time k (required if store_press=.true.)
  !
  ! Notes:
  ! - Call this once per sampling point (e.g., every output_interval steps).
  ! - The spacing between samples defines the physical time associated with one "lag".
  !====================================================================================================
  subroutine md_corr_add_sample(acc, epot, ekin, temp, press)

    type(md_corr_accum), intent(inout) :: acc
    real(kind=dp_kind),  intent(in)    :: epot, ekin
    real(kind=dp_kind),  intent(in), optional :: temp, press

    real(kind=dp_kind) :: etot

    if (.not. acc%enabled) stop 'md_corr_add_sample(): correlations not enabled (call md_corr_init()).'
    if (acc%n_samples_used >= acc%n_alloc) stop 'md_corr_add_sample(): series buffer is full.'

    if (acc%store_temp) then
      if (.not. present(temp)) stop 'md_corr_add_sample(): temp expected but not provided.'
    end if
    if (acc%store_press) then
      if (.not. present(press)) stop 'md_corr_add_sample(): press expected but not provided.'
    end if

    etot = epot + ekin

    acc%n_samples_used = acc%n_samples_used + 1_int_kind

    acc%epot_series(acc%n_samples_used) = epot
    acc%ekin_series(acc%n_samples_used) = ekin
    acc%etot_series(acc%n_samples_used) = etot

    if (acc%store_temp) acc%temp_series(acc%n_samples_used) = temp
    if (acc%store_press) acc%press_series(acc%n_samples_used) = press

    acc%computed = .false.
    acc%cm_computed = .false.

  end subroutine md_corr_add_sample

  
  !====================================================================================================
  ! Compute autocorrelations for all stored observables.
  !
  ! Outputs (stored in accumulator arrays):
  !   corr_*(0:max_lag)  : autocorrelation curve C(L)
  !   corrn_*(0:max_lag) : normalized autocorrelation C(L)/C(0)
  !
  ! Options:
  !   centered = .false. (default): raw autocorrelation             C(L) = (1/(N-L)) * sum_{k=1...N-L} A(k) * A(k+L)
  !   centered = .true.: centered autocorrelation (autocovariance)  C_c(L) = (1/(N-L)) * sum_{k=1...N-L} [ (A(k) - <A>) * (A(k+L) - <A>) ]     
  !
  ! Requirements:
  ! - Need at least 2 samples.
  ! - Need max_lag < n_samples_used so there are valid pairs for every lag.
  !====================================================================================================
  subroutine md_corr_compute(acc, centered)

    type(md_corr_accum), intent(inout)        :: acc
    logical,             intent(in), optional :: centered
    logical                                   :: do_center
    
    do_center = .false.
    if (present(centered)) do_center = centered

    if (.not. acc%enabled) stop 'md_corr_compute(): correlations not enabled.'
    if (acc%n_samples_used <= 1_int_kind) stop 'md_corr_compute(): not enough samples (need at least 2).'
    if (acc%max_lag >= acc%n_samples_used) stop 'md_corr_compute(): max_lag too large for stored samples.'

    ! Energies (always)
    ! NOTE: autocorr_scalar* effectively works on acc%*_series(1:n_samples_used); passing the full array acc%*_series is safe because n_samples_in limits the range used.
    if (do_center) then

      call autocorr_scalar_centered(acc%epot_series,  acc%n_samples_used, acc%max_lag, acc%corr_epot)
      call normalize_corr(acc%max_lag, acc%corr_epot, acc%corrn_epot)

      call autocorr_scalar_centered(acc%ekin_series,  acc%n_samples_used, acc%max_lag, acc%corr_ekin)
      call normalize_corr(acc%max_lag, acc%corr_ekin, acc%corrn_ekin)

      call autocorr_scalar_centered(acc%etot_series,  acc%n_samples_used, acc%max_lag, acc%corr_etot)
      call normalize_corr(acc%max_lag, acc%corr_etot, acc%corrn_etot)

      if (acc%store_temp) then
        call autocorr_scalar_centered(acc%temp_series, acc%n_samples_used, acc%max_lag, acc%corr_temp)
        call normalize_corr(acc%max_lag, acc%corr_temp, acc%corrn_temp)
      end if

      if (acc%store_press) then
        call autocorr_scalar_centered(acc%press_series, acc%n_samples_used, acc%max_lag, acc%corr_press)
        call normalize_corr(acc%max_lag, acc%corr_press, acc%corrn_press)
      end if

    else

      call autocorr_scalar(acc%epot_series,  acc%n_samples_used, acc%max_lag, acc%corr_epot)
      call normalize_corr(acc%max_lag, acc%corr_epot, acc%corrn_epot)

      call autocorr_scalar(acc%ekin_series,  acc%n_samples_used, acc%max_lag, acc%corr_ekin)
      call normalize_corr(acc%max_lag, acc%corr_ekin, acc%corrn_ekin)

      call autocorr_scalar(acc%etot_series,  acc%n_samples_used, acc%max_lag, acc%corr_etot)
      call normalize_corr(acc%max_lag, acc%corr_etot, acc%corrn_etot)

      if (acc%store_temp) then
        call autocorr_scalar(acc%temp_series, acc%n_samples_used, acc%max_lag, acc%corr_temp)
        call normalize_corr(acc%max_lag, acc%corr_temp, acc%corrn_temp)
      end if

      if (acc%store_press) then
        call autocorr_scalar(acc%press_series, acc%n_samples_used, acc%max_lag, acc%corr_press)
        call normalize_corr(acc%max_lag, acc%corr_press, acc%corrn_press)
      end if

    end if

    acc%computed = .true.

  end subroutine md_corr_compute


  !====================================================================================================
  ! Return metadata (useful for allocation and I/O).
  !====================================================================================================
  subroutine md_corr_get_info(acc, n_samples_used_out, lag_max_out, has_temp, has_press)

    type(md_corr_accum), intent(in) :: acc
    integer(kind=int_kind), intent(out) :: n_samples_used_out, lag_max_out
    logical, intent(out) :: has_temp, has_press

    if (.not. acc%enabled) stop 'md_corr_get_info(): correlations not enabled.'

    n_samples_used_out  = acc%n_samples_used
    lag_max_out         = acc%max_lag
    has_temp            = acc%store_temp
    has_press           = acc%store_press

  end subroutine md_corr_get_info


  !====================================================================================================
  ! Get ONE observable correlation curve (raw + normalized) by NAME.
  !  - Names accepted: 'epot', 'ekin', 'etot', 'temp', 'temperature', 'press', 'pressure'
  !
  ! Outputs:
  !   corr_out(0:max_lag), corrn_out(0:max_lag)
  !
  ! Notes:
  ! - Caller must provide arrays corr_out, corrn_out with upper bound >= max_lag.
  ! - md_corr_compute() must have been called first.
  !====================================================================================================
  subroutine md_corr_get_by_name(acc, name, corr_out, corrn_out)

    type(md_corr_accum), intent(in) :: acc
    character(len=*),    intent(in) :: name
    real(kind=dp_kind),  intent(out) :: corr_out(0:)
    real(kind=dp_kind),  intent(out) :: corrn_out(0:)

    if (.not. acc%enabled)  stop 'md_corr_get_by_name(): correlations not enabled.'
    if (.not. acc%computed) stop 'md_corr_get_by_name(): correlations not computed (call md_corr_compute()).'
    if (ubound(corr_out,1)  < acc%max_lag) stop 'md_corr_get_by_name(): corr_out too small.'
    if (ubound(corrn_out,1) < acc%max_lag) stop 'md_corr_get_by_name(): corrn_out too small.'

    select case (name)
    case ('epot')
      corr_out(0:acc%max_lag)  = acc%corr_epot(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%corrn_epot(0:acc%max_lag)
    case ('ekin')
      corr_out(0:acc%max_lag)  = acc%corr_ekin(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%corrn_ekin(0:acc%max_lag)
    case ('etot')
      corr_out(0:acc%max_lag)  = acc%corr_etot(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%corrn_etot(0:acc%max_lag)
    case ('temp','temperature')
      if (.not. acc%store_temp) stop 'md_corr_get_by_name(): temperature correlation not enabled at init.'
      corr_out(0:acc%max_lag)  = acc%corr_temp(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%corrn_temp(0:acc%max_lag)
    case ('press','pressure')
      if (.not. acc%store_press) stop 'md_corr_get_by_name(): pressure correlation not enabled at init.'
      corr_out(0:acc%max_lag)  = acc%corr_press(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%corrn_press(0:acc%max_lag)
    case default
      stop 'md_corr_get_by_name(): unknown name. Use epot/ekin/etot/temperature/pressure.'
    end select

  end subroutine md_corr_get_by_name


  !====================================================================================================
  ! Get ALL correlations at once.
  !
  ! Outputs (caller must provide arrays with upper bound >= max_lag):
  !   corr_epot_out, corrn_epot_out, corr_ekin_out, corrn_ekin_out, corr_etot_out, corrn_etot_out
  !
  ! Outputs Optional (only if enabled at init):
  !   corr_temp_out, corrn_temp_out, corr_press_out, corrn_press_out
  !
  ! Notes:
  ! - md_corr_compute() must have been called before.
  ! - This routine does NOT compute anything, it only copies stored results.
  !====================================================================================================
  subroutine md_corr_get_all(acc, &
                          corr_epot_out,  corrn_epot_out, &
                          corr_ekin_out,  corrn_ekin_out, &
                          corr_etot_out,  corrn_etot_out, &
                          corr_temp_out,  corrn_temp_out, &
                          corr_press_out, corrn_press_out)

    type(md_corr_accum), intent(in) :: acc

    ! Mandatory outputs
    real(kind=dp_kind), intent(out) :: corr_epot_out(0:),  corrn_epot_out(0:)
    real(kind=dp_kind), intent(out) :: corr_ekin_out(0:),  corrn_ekin_out(0:)
    real(kind=dp_kind), intent(out) :: corr_etot_out(0:),  corrn_etot_out(0:)

    ! Optional outputs
    real(kind=dp_kind), intent(out), optional :: corr_temp_out(0:),  corrn_temp_out(0:)
    real(kind=dp_kind), intent(out), optional :: corr_press_out(0:), corrn_press_out(0:)

    ! Common safety checks (fast fail before doing anything)
    if (.not. acc%enabled)  stop 'md_corr_get_all(): correlations not enabled.'
    if (.not. acc%computed) stop 'md_corr_get_all(): correlations not computed (call md_corr_compute()).'

    ! Mandatory correlations
    call md_corr_get_by_name(acc, 'epot', corr_epot_out, corrn_epot_out)
    call md_corr_get_by_name(acc, 'ekin', corr_ekin_out, corrn_ekin_out)
    call md_corr_get_by_name(acc, 'etot', corr_etot_out, corrn_etot_out)

    ! Optional: temperature (must come as a pair)
    if (present(corr_temp_out) .or. present(corrn_temp_out)) then
      if (.not. present(corr_temp_out) .or. .not. present(corrn_temp_out)) stop 'md_corr_get_all(): both corr_temp_out and corrn_temp_out must be present.'
      call md_corr_get_by_name(acc, 'temperature', corr_temp_out, corrn_temp_out)
    end if

    ! Optional: pressure (must come as a pair)
    if (present(corr_press_out) .or. present(corrn_press_out)) then
      if (.not. present(corr_press_out) .or. .not. present(corrn_press_out)) stop 'md_corr_get_all(): both corr_press_out and corrn_press_out must be present.'
      call md_corr_get_by_name(acc, 'pressure', corr_press_out, corrn_press_out)
    end if

  end subroutine md_corr_get_all


  !====================================================================================================
  ! Free memory and reset accumulator.
  !====================================================================================================
  subroutine md_corr_free(acc)

    type(md_corr_accum), intent(inout) :: acc

    if (allocated(acc%epot_series)) deallocate(acc%epot_series)
    if (allocated(acc%ekin_series)) deallocate(acc%ekin_series)
    if (allocated(acc%etot_series)) deallocate(acc%etot_series)
    if (allocated(acc%temp_series)) deallocate(acc%temp_series)
    if (allocated(acc%press_series)) deallocate(acc%press_series)

    if (allocated(acc%corr_epot))  deallocate(acc%corr_epot)
    if (allocated(acc%corrn_epot)) deallocate(acc%corrn_epot)
    if (allocated(acc%corr_ekin))  deallocate(acc%corr_ekin)
    if (allocated(acc%corrn_ekin)) deallocate(acc%corrn_ekin)
    if (allocated(acc%corr_etot))  deallocate(acc%corr_etot)
    if (allocated(acc%corrn_etot)) deallocate(acc%corrn_etot)

    if (allocated(acc%corr_temp))  deallocate(acc%corr_temp)
    if (allocated(acc%corrn_temp)) deallocate(acc%corrn_temp)
    if (allocated(acc%corr_press))  deallocate(acc%corr_press)
    if (allocated(acc%corrn_press)) deallocate(acc%corrn_press)

    call md_corr_cm_free(acc)

    acc%enabled         = .false.
    acc%computed        = .false.
    acc%n_alloc         = 0_int_kind
    acc%n_samples_used  = 0_int_kind
    acc%max_lag         = 0_int_kind
    acc%store_temp      = .false.
    acc%store_press     = .false.

    acc%n      = 0_int_kind
    acc%volume = 0.d0

  end subroutine md_corr_free


  !****************************************************************************************************
  ! CORRELATION MEANS (block averaging) API
  !
  ! Goal:
  ! - Given the internally stored time series, compute a smoother correlation estimate by:
  !     1) Splitting the full series into B contiguous blocks (approximately independent).
  !     2) Computing ONE autocorrelation curve per block:
  !          C_b(tau) = (1/(N_b - tau)) * sum_{k=1..N_b-tau} A_b(k) * A_b(k+tau)
  !     3) Averaging the curves over blocks:
  !          <C(tau)>      = (1/B) * sum_b C_b(tau)
  !          <C_norm(tau)> = (1/B) * sum_b [C_b(tau)/C_b(0)]
  !
  ! Observables supported (same pattern as md_correlations):
  ! - Mandatory: Epot, Ekin, Etot
  ! - Optional : Temperature T, Pressure P
  !
  ! Design notes:
  ! - Correlation means use the stored time series as full-series correlations and divides them into different blocks
  !
  ! Pattern:
  !   md_corr_cm_init   -> md_corr_cm_compute -> md_corr_cm_get (individual or all) -> md_corr_cm_free
  !****************************************************************************************************


  !====================================================================================================
  ! Initialize / reset correlation-means buffers.
  !
  ! Inputs:
  !  - n_blocks_in : number of blocks B (>= 1)
  !
  ! Notes:
  ! - Uses acc%max_lag as the lag size for the averaged curves (same as full-series correlations).
  ! - block_len is defined later in md_corr_cm_compute() because it depends on acc%n_samples_used.
  !====================================================================================================
  subroutine md_corr_cm_init(acc, n_blocks_in)

    type(md_corr_accum),      intent(inout) :: acc
    integer(kind=int_kind),   intent(in)    :: n_blocks_in

    if (.not. acc%enabled) stop 'md_corr_cm_init(): correlations not enabled (call md_corr_init()).'
    if (n_blocks_in <= 0_int_kind) stop 'md_corr_cm_init(): n_blocks_in must be >= 1.'
    if (acc%max_lag < 0_int_kind) stop 'md_corr_cm_init(): invalid max_lag.'

    call md_corr_cm_free(acc)

    acc%cm_enabled  = .true.
    acc%cm_computed = .false.
    acc%cm_n_blocks = n_blocks_in
    acc%cm_block_len = 0_int_kind

    ! -----------------------
    ! Allocate sums (per-lag): 
    !   Array of sums (C_0(0) + C_1(0) + ... + C_B(0), C_0(1) + C_1(1) + ... + C_B(1), ..., C_0(n_valid) + C_1(n_valid) + ... + C_B(n_valid))
    !   C_b(tau) is the correlation for (epot, ekin, etot, temp, press) for block b and lag tau
    ! -----------------------
    allocate(acc%cm_sum_corr_epot(0:acc%max_lag),  acc%cm_sum_corrn_epot(0:acc%max_lag))
    allocate(acc%cm_sum_corr_ekin(0:acc%max_lag),  acc%cm_sum_corrn_ekin(0:acc%max_lag))
    allocate(acc%cm_sum_corr_etot(0:acc%max_lag),  acc%cm_sum_corrn_etot(0:acc%max_lag))

    acc%cm_sum_corr_epot  = 0.d0; acc%cm_sum_corrn_epot = 0.d0
    acc%cm_sum_corr_ekin  = 0.d0; acc%cm_sum_corrn_ekin = 0.d0
    acc%cm_sum_corr_etot  = 0.d0; acc%cm_sum_corrn_etot = 0.d0

    if (acc%store_temp) then
      allocate(acc%cm_sum_corr_temp(0:acc%max_lag), acc%cm_sum_corrn_temp(0:acc%max_lag))
      acc%cm_sum_corr_temp  = 0.d0
      acc%cm_sum_corrn_temp = 0.d0
    end if

    if (acc%store_press) then
      allocate(acc%cm_sum_corr_press(0:acc%max_lag), acc%cm_sum_corrn_press(0:acc%max_lag))
      acc%cm_sum_corr_press  = 0.d0
      acc%cm_sum_corrn_press = 0.d0
    end if

    ! -----------------------
    ! Allocate output means (per-lag): 
    !  Array of correlation means (1/B sum_{b=1}^{B} c_b(0), 1/B sum_{b=1}^{B} c_b(1), ..., 1/B sum_{b=1}^{B} c_b(n_valid))
    ! -----------------------
    allocate(acc%cm_mean_corr_epot(0:acc%max_lag),  acc%cm_mean_corrn_epot(0:acc%max_lag))
    allocate(acc%cm_mean_corr_ekin(0:acc%max_lag),  acc%cm_mean_corrn_ekin(0:acc%max_lag))
    allocate(acc%cm_mean_corr_etot(0:acc%max_lag),  acc%cm_mean_corrn_etot(0:acc%max_lag))

    acc%cm_mean_corr_epot  = 0.d0; acc%cm_mean_corrn_epot = 0.d0
    acc%cm_mean_corr_ekin  = 0.d0; acc%cm_mean_corrn_ekin = 0.d0
    acc%cm_mean_corr_etot  = 0.d0; acc%cm_mean_corrn_etot = 0.d0

    if (acc%store_temp) then
      allocate(acc%cm_mean_corr_temp(0:acc%max_lag), acc%cm_mean_corrn_temp(0:acc%max_lag))
      acc%cm_mean_corr_temp  = 0.d0
      acc%cm_mean_corrn_temp = 0.d0
    end if

    if (acc%store_press) then
      allocate(acc%cm_mean_corr_press(0:acc%max_lag), acc%cm_mean_corrn_press(0:acc%max_lag))
      acc%cm_mean_corr_press  = 0.d0
      acc%cm_mean_corrn_press = 0.d0
    end if

  end subroutine md_corr_cm_init


  !====================================================================================================
  ! Compute block-averaged correlations from the internally stored series.
  !
  ! NOTE:
  ! - Blocks are built as contiguous segments of equal size: block_len = floor(N_samples_used / n_blocks)
  !   Remaining samples (if N_samples_used is not divisible by n_blocks) are ignored at the end.
  !
  ! - A block is only valid if block_len >= max_lag + 1 (needs at least one valid pair at max lag).
  !
  ! Correlations inside one block:
  !  - For each block b, we compute an autocorrelation function C_b(tau) = < A_b(k) * A_b(k+tau) >_k
  !    where the average is over all valid k inside THAT block only.
  !  - This gives ONE estimate of the correlation function per block.
  !
  ! Options:
  !   centered = .false. (default): raw autocorrelation
  !   centered = .true.           : centered autocorrelation (autocovariance)
  !====================================================================================================
  subroutine md_corr_cm_compute(acc, centered)

    type(md_corr_accum), intent(inout)        :: acc
    logical,             intent(in), optional :: centered

    integer(kind=int_kind) :: b                ! Block counter
    integer(kind=int_kind) :: b_start, b_end   ! Indexes defining where a block starts and ends
    integer(kind=int_kind) :: n_block_samples  ! Num samples per block
    real(kind=dp_kind)     :: inv_nb           ! 1 / n_blocks

    logical :: do_center

    ! Work arrays for ONE block (raw + normalized)
    real(kind=dp_kind), allocatable :: corr_block(:), corrn_block(:)

    if (.not. acc%cm_enabled) stop 'md_corr_cm_compute(): correlation means not enabled (call md_corr_cm_init()).'
    if (acc%n_samples_used <= 1_int_kind) stop 'md_corr_cm_compute(): need at least 2 samples.'
    if (acc%cm_n_blocks <= 0_int_kind) stop 'md_corr_cm_compute(): invalid number of blocks.'
    if (acc%max_lag >= acc%n_samples_used) stop 'md_corr_cm_compute(): max_lag too large for stored samples.'

    do_center = .false.
    if (present(centered)) do_center = centered

    ! Build equal-size blocks
    acc%cm_block_len = acc%n_samples_used / acc%cm_n_blocks
    if (acc%cm_block_len <= 0_int_kind) stop 'md_corr_cm_compute(): block_len <= 0 (too many blocks).'
    if (acc%max_lag >= acc%cm_block_len) stop 'md_corr_cm_compute(): max_lag must be < block_len.'

    ! Reset sums (so compute can be called again safely)
    acc%cm_sum_corr_epot   = 0.d0; acc%cm_sum_corrn_epot  = 0.d0
    acc%cm_sum_corr_ekin   = 0.d0; acc%cm_sum_corrn_ekin  = 0.d0
    acc%cm_sum_corr_etot   = 0.d0; acc%cm_sum_corrn_etot  = 0.d0
    if (acc%store_temp)  then; acc%cm_sum_corr_temp  = 0.d0; acc%cm_sum_corrn_temp  = 0.d0; end if
    if (acc%store_press) then; acc%cm_sum_corr_press = 0.d0; acc%cm_sum_corrn_press = 0.d0; end if

    allocate(corr_block(0:acc%max_lag), corrn_block(0:acc%max_lag))

    ! -----------------------
    ! Loop over blocks and accumulate one correlation curve per block
    ! -----------------------
    do b = 1_int_kind, acc%cm_n_blocks

      b_start = (b-1_int_kind)*acc%cm_block_len + 1_int_kind
      b_end   = b_start + acc%cm_block_len - 1_int_kind
      n_block_samples = acc%cm_block_len

      ! Epot
      if (do_center) then
        call autocorr_scalar_centered(acc%epot_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      else
        call autocorr_scalar(acc%epot_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      end if
      call normalize_corr(acc%max_lag, corr_block, corrn_block)
      ! For each block b, compute epot correlation: 
      ! b = 1 (0, 0, ..., 0) + (C_0(0), C_0(1), ..., C_0(n_valid)) ! b = 2 (C_0(0) + C_1(0), C_0(1) + C_1(1), ..., C_0(n_valid) + C_1(n_valid)) 
      ! b = 3 (C_0(0) + C_1(0) + C_2(0), C_0(1) + C_1(1) + C_2(1), ..., C_0(n_valid) + C_1(n_valid) + C_3(n_valid)) 
      ! 
      ! At the end: (C_0(0) + C_1(0) + ... + C_B(0), C_0(1) + C_1(1) + ... + C_B(1), ..., C_0(n_valid) + C_1(n_valid) + ... + C_B(n_valid)) where B is total number of Blocks
      acc%cm_sum_corr_epot(0:acc%max_lag)   = acc%cm_sum_corr_epot(0:acc%max_lag)   + corr_block(0:acc%max_lag)
      acc%cm_sum_corrn_epot(0:acc%max_lag)  = acc%cm_sum_corrn_epot(0:acc%max_lag)  + corrn_block(0:acc%max_lag)

      ! Ekin
      if (do_center) then
        call autocorr_scalar_centered(acc%ekin_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      else
        call autocorr_scalar(acc%ekin_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      end if
      call normalize_corr(acc%max_lag, corr_block, corrn_block)
      acc%cm_sum_corr_ekin(0:acc%max_lag)   = acc%cm_sum_corr_ekin(0:acc%max_lag)   + corr_block(0:acc%max_lag)
      acc%cm_sum_corrn_ekin(0:acc%max_lag)  = acc%cm_sum_corrn_ekin(0:acc%max_lag)  + corrn_block(0:acc%max_lag)

      ! Etot
      if (do_center) then
        call autocorr_scalar_centered(acc%etot_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      else
        call autocorr_scalar(acc%etot_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
      end if
      call normalize_corr(acc%max_lag, corr_block, corrn_block)
      acc%cm_sum_corr_etot(0:acc%max_lag)   = acc%cm_sum_corr_etot(0:acc%max_lag)   + corr_block(0:acc%max_lag)
      acc%cm_sum_corrn_etot(0:acc%max_lag)  = acc%cm_sum_corrn_etot(0:acc%max_lag)  + corrn_block(0:acc%max_lag)

      ! Optional: Temperature
      if (acc%store_temp) then
        if (do_center) then
          call autocorr_scalar_centered(acc%temp_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
        else
          call autocorr_scalar(acc%temp_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
        end if
        call normalize_corr(acc%max_lag, corr_block, corrn_block)
        acc%cm_sum_corr_temp(0:acc%max_lag)   = acc%cm_sum_corr_temp(0:acc%max_lag)   + corr_block(0:acc%max_lag)
        acc%cm_sum_corrn_temp(0:acc%max_lag)  = acc%cm_sum_corrn_temp(0:acc%max_lag)  + corrn_block(0:acc%max_lag)
      end if

      ! Optional: Pressure
      if (acc%store_press) then
        if (do_center) then
          call autocorr_scalar_centered(acc%press_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
        else
          call autocorr_scalar(acc%press_series(b_start:b_end), n_block_samples, acc%max_lag, corr_block)
        end if
        call normalize_corr(acc%max_lag, corr_block, corrn_block)
        acc%cm_sum_corr_press(0:acc%max_lag)  = acc%cm_sum_corr_press(0:acc%max_lag)  + corr_block(0:acc%max_lag)
        acc%cm_sum_corrn_press(0:acc%max_lag) = acc%cm_sum_corrn_press(0:acc%max_lag) + corrn_block(0:acc%max_lag)
      end if

    end do

    deallocate(corr_block, corrn_block)

    ! -----------------------
    ! Final averaging step (divide sums by B)
    ! -----------------------
    inv_nb = 1.d0 / dble(acc%cm_n_blocks)

    acc%cm_mean_corr_epot(0:acc%max_lag)   = acc%cm_sum_corr_epot(0:acc%max_lag)   * inv_nb
    acc%cm_mean_corrn_epot(0:acc%max_lag)  = acc%cm_sum_corrn_epot(0:acc%max_lag)  * inv_nb

    acc%cm_mean_corr_ekin(0:acc%max_lag)   = acc%cm_sum_corr_ekin(0:acc%max_lag)   * inv_nb
    acc%cm_mean_corrn_ekin(0:acc%max_lag)  = acc%cm_sum_corrn_ekin(0:acc%max_lag)  * inv_nb

    acc%cm_mean_corr_etot(0:acc%max_lag)   = acc%cm_sum_corr_etot(0:acc%max_lag)   * inv_nb
    acc%cm_mean_corrn_etot(0:acc%max_lag)  = acc%cm_sum_corrn_etot(0:acc%max_lag)  * inv_nb

    if (acc%store_temp) then
      acc%cm_mean_corr_temp(0:acc%max_lag)   = acc%cm_sum_corr_temp(0:acc%max_lag)   * inv_nb
      acc%cm_mean_corrn_temp(0:acc%max_lag)  = acc%cm_sum_corrn_temp(0:acc%max_lag)  * inv_nb
    end if

    if (acc%store_press) then
      acc%cm_mean_corr_press(0:acc%max_lag)  = acc%cm_sum_corr_press(0:acc%max_lag)  * inv_nb
      acc%cm_mean_corrn_press(0:acc%max_lag) = acc%cm_sum_corrn_press(0:acc%max_lag) * inv_nb
    end if

    acc%cm_computed = .true.

  end subroutine md_corr_cm_compute

  
  
  !====================================================================================================
  ! Get ONE observable averaged correlation curve (by name)
  !  - Names accepted: 'epot', 'ekin', 'etot', 'temp', 'temperature', 'press', 'pressure'
  !
  ! Outputs:
  !   corr_out(0:max_lag), corrn_out(0:max_lag)
  !====================================================================================================
  subroutine md_corr_cm_get_by_name(acc, name, corr_out, corrn_out)

    type(md_corr_accum), intent(in) :: acc
    character(len=*),    intent(in) :: name
    real(kind=dp_kind),  intent(out) :: corr_out(0:)
    real(kind=dp_kind),  intent(out) :: corrn_out(0:)

    if (.not. acc%cm_enabled)  stop 'md_corr_cm_get_by_name(): not enabled.'
    if (.not. acc%cm_computed) stop 'md_corr_cm_get_by_name(): not computed (call md_corr_cm_compute()).'
    if (ubound(corr_out,1)  < acc%max_lag) stop 'md_corr_cm_get_by_name(): corr_out too small.'
    if (ubound(corrn_out,1) < acc%max_lag) stop 'md_corr_cm_get_by_name(): corrn_out too small.'

    select case (name)
    case ('epot')
      corr_out(0:acc%max_lag)  = acc%cm_mean_corr_epot(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%cm_mean_corrn_epot(0:acc%max_lag)
    case ('ekin')
      corr_out(0:acc%max_lag)  = acc%cm_mean_corr_ekin(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%cm_mean_corrn_ekin(0:acc%max_lag)
    case ('etot')
      corr_out(0:acc%max_lag)  = acc%cm_mean_corr_etot(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%cm_mean_corrn_etot(0:acc%max_lag)
    case ('temp', 'temperature')
      if (.not. acc%store_temp) stop 'md_corr_cm_get_by_name(): temperature means not enabled.'
      corr_out(0:acc%max_lag)  = acc%cm_mean_corr_temp(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%cm_mean_corrn_temp(0:acc%max_lag)
    case ('press', 'pressure')
      if (.not. acc%store_press) stop 'md_corr_cm_get_by_name(): pressure means not enabled.'
      corr_out(0:acc%max_lag)  = acc%cm_mean_corr_press(0:acc%max_lag)
      corrn_out(0:acc%max_lag) = acc%cm_mean_corrn_press(0:acc%max_lag)
    case default
      stop 'md_corr_cm_get_by_name(): unknown observable name.'
    end select

  end subroutine md_corr_cm_get_by_name


  !====================================================================================================
  ! Get ALL averaged curves at once (mandatory + optional)
  !
  ! Mandatory outputs:
  !   epot, ekin, etot (raw + normalized)
  !
  ! Optional outputs:
  !   temp, press (raw + normalized)  -> must be passed in pairs   temp, press (raw + normalized)  -> must be passed in pairs
  !
  ! Notes:
  ! - This routine does NOT compute anything; it only copies stored results.
  !====================================================================================================
  subroutine md_corr_cm_get_all(acc, &
                             corr_epot_out, corrn_epot_out, &
                             corr_ekin_out, corrn_ekin_out, &
                             corr_etot_out, corrn_etot_out, &
                             corr_temp_out, corrn_temp_out, &
                             corr_press_out, corrn_press_out)

    type(md_corr_accum), intent(in) :: acc

    real(kind=dp_kind), intent(out) :: corr_epot_out(0:), corrn_epot_out(0:)
    real(kind=dp_kind), intent(out) :: corr_ekin_out(0:), corrn_ekin_out(0:)
    real(kind=dp_kind), intent(out) :: corr_etot_out(0:), corrn_etot_out(0:)

    real(kind=dp_kind), intent(out), optional :: corr_temp_out(0:),  corrn_temp_out(0:)
    real(kind=dp_kind), intent(out), optional :: corr_press_out(0:), corrn_press_out(0:)

    if (.not. acc%cm_enabled)  stop 'md_corr_cm_get_all(): not enabled.'
    if (.not. acc%cm_computed) stop 'md_corr_cm_get_all(): not computed.'

    if (ubound(corr_epot_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corr_epot_out too small.'
    if (ubound(corrn_epot_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corrn_epot_out too small.'
    if (ubound(corr_ekin_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corr_ekin_out too small.'
    if (ubound(corrn_ekin_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corrn_ekin_out too small.'
    if (ubound(corr_etot_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corr_etot_out too small.'
    if (ubound(corrn_etot_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corrn_etot_out too small.'

    ! Mandatory curves
    call md_corr_cm_get_by_name(acc, 'epot', corr_epot_out, corrn_epot_out)
    call md_corr_cm_get_by_name(acc, 'ekin', corr_ekin_out, corrn_ekin_out)
    call md_corr_cm_get_by_name(acc, 'etot', corr_etot_out, corrn_etot_out)

    ! Optional: temperature (must come as a pair)
    if (present(corr_temp_out) .or. present(corrn_temp_out)) then
    if (.not. present(corr_temp_out) .or. .not. present(corrn_temp_out)) &
      stop 'md_corr_cm_get_all(): both corr_temp_out and corrn_temp_out must be present.'
    if (ubound(corr_temp_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corr_temp_out too small.'
    if (ubound(corrn_temp_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corrn_temp_out too small.'
    call md_corr_cm_get_by_name(acc, 'temperature', corr_temp_out, corrn_temp_out)
    end if

    ! Optional: pressure (must come as a pair)
    if (present(corr_press_out) .or. present(corrn_press_out)) then
    if (.not. present(corr_press_out) .or. .not. present(corrn_press_out)) &
      stop 'md_corr_cm_get_all(): both corr_press_out and corrn_press_out must be present.'
    if (ubound(corr_press_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corr_press_out too small.'
    if (ubound(corrn_press_out,1) < acc%max_lag) stop 'md_corr_cm_get_all(): corrn_press_out too small.'
    call md_corr_cm_get_by_name(acc, 'pressure', corr_press_out, corrn_press_out)
    end if

  end subroutine md_corr_cm_get_all


  !====================================================================================================
  ! Free / reset correlation means buffers
  !====================================================================================================
  subroutine md_corr_cm_free(acc)

    type(md_corr_accum), intent(inout) :: acc

    if (allocated(acc%cm_sum_corr_epot))   deallocate(acc%cm_sum_corr_epot)
    if (allocated(acc%cm_sum_corrn_epot))  deallocate(acc%cm_sum_corrn_epot)
    if (allocated(acc%cm_sum_corr_ekin))   deallocate(acc%cm_sum_corr_ekin)
    if (allocated(acc%cm_sum_corrn_ekin))  deallocate(acc%cm_sum_corrn_ekin)
    if (allocated(acc%cm_sum_corr_etot))   deallocate(acc%cm_sum_corr_etot)
    if (allocated(acc%cm_sum_corrn_etot))  deallocate(acc%cm_sum_corrn_etot)

    if (allocated(acc%cm_sum_corr_temp))   deallocate(acc%cm_sum_corr_temp)
    if (allocated(acc%cm_sum_corrn_temp))  deallocate(acc%cm_sum_corrn_temp)
    if (allocated(acc%cm_sum_corr_press))  deallocate(acc%cm_sum_corr_press)
    if (allocated(acc%cm_sum_corrn_press)) deallocate(acc%cm_sum_corrn_press)

    if (allocated(acc%cm_mean_corr_epot))  deallocate(acc%cm_mean_corr_epot)
    if (allocated(acc%cm_mean_corrn_epot)) deallocate(acc%cm_mean_corrn_epot)
    if (allocated(acc%cm_mean_corr_ekin))  deallocate(acc%cm_mean_corr_ekin)
    if (allocated(acc%cm_mean_corrn_ekin)) deallocate(acc%cm_mean_corrn_ekin)
    if (allocated(acc%cm_mean_corr_etot))  deallocate(acc%cm_mean_corr_etot)
    if (allocated(acc%cm_mean_corrn_etot)) deallocate(acc%cm_mean_corrn_etot)

    if (allocated(acc%cm_mean_corr_temp))  deallocate(acc%cm_mean_corr_temp)
    if (allocated(acc%cm_mean_corrn_temp)) deallocate(acc%cm_mean_corrn_temp)
    if (allocated(acc%cm_mean_corr_press))  deallocate(acc%cm_mean_corr_press)
    if (allocated(acc%cm_mean_corrn_press)) deallocate(acc%cm_mean_corrn_press)

    acc%cm_enabled   = .false.
    acc%cm_computed  = .false.
    acc%cm_n_blocks  = 0_int_kind
    acc%cm_block_len = 0_int_kind

  end subroutine md_corr_cm_free


  !----------------------------------------------------------------------------------------------------
  ! Return number of blocks used.
  !----------------------------------------------------------------------------------------------------
  function md_corr_cm_num_blocks(acc) result(nb)
    type(md_corr_accum), intent(in) :: acc
    integer(kind=int_kind) :: nb
    nb = acc%cm_n_blocks
  end function md_corr_cm_num_blocks


  !----------------------------------------------------------------------------------------------------
  ! Return maximum lag.
  !----------------------------------------------------------------------------------------------------
  function md_corr_cm_max_lag(acc) result(ml)
    type(md_corr_accum), intent(in) :: acc
    integer(kind=int_kind) :: ml
    ml = acc%max_lag
  end function md_corr_cm_max_lag

end module md_correlations
