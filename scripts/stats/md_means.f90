!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! MD_MEANS MODULE (time averages and fluctuations for scalar observables sampled during MD)
!
! Goal:
! - Accumulate (streaming) time averages (means) and fluctuations (std dev) of scalar observables
!   needed by md_simulation.f90.
! - Avoid storing full time series (only sums and counters are kept).
! - Keep md_simulation clean: all bookkeeping for means stays here.
!
! Quantities handled (reduced units, kB = 1, m = 1):
! - Energies: E_pot, E_kin, E_tot = E_pot + E_kin
! - Temperature from kinetic energy:  T = 2 E_kin / (3 N)
! - Optional pressure if virial is available: P = rho*T + W/(3V), where rho = N/V and W is the virial
! - Other quantities needed for thermodynamic coefficients calculations
!
! Mean and mean of squares:
! - Mean: <x> = (1/n) sum_k x_k
! - Fluctuations (std dev): std = sqrt( <x^2> - <x>^2 )
!
! Design notes:
! - This module stores NO global state: everything lives inside type(md_means_accum).
! - Pattern:
!     1) md_means_init()
!     2) md_means_add_sample()  [called once per sampling time]
!     3) md_means_get()         [retrieve means + std devs]
!
! Naming legend (THIS MODULE ONLY):
!   Basics:
!     U     : potential energy (E_pot)
!     K     : kinetic energy (E_kin)
!     E     : total energy (E_tot = U + K)
!     T     : temperature from equipartition (T = 2K/(3N))
!     P     : pressure (if virial available)
!     Kinv  : inverse kinetic energy, 1/K
!
!   Potential-derivative observables (provided by LJ module):
!     dU    : d_epot  = sum_{i<j} r_ij dU/dr_ij
!     ddU   : dd_epot = sum_{i<j} r_ij^2 d^2U/dr_ij^2
!     dU2   : (dU)^2  (square of dU; NOT ddU)
!
!   Accumulator naming convention:
!     sum_X     : sum over samples of X
!     sum_sq_X  : sum over samples of X^2
!     *_Kinv    : product with (1/K), i.e. X * (1/K)
!******************************************************************************************************
module md_means

  use define_precision, only: dp_kind, int_kind
  use md_types,         only: sim_params
  use stats_math,       only: stats_std_from_moments
  implicit none

  private

  ! Public API
  public :: md_means_accum
  public :: md_means_init, md_means_reset
  public :: md_means_add_sample
  public :: md_means_get

  ! ==================================================================================
  ! Accumulator container (NO module-level state)
  ! ==================================================================================
  type :: md_means_accum

    ! ==================================================================================
    ! 1) Run configuration snapshot
    ! ==================================================================================
    integer(kind=int_kind) :: n      = 0_int_kind   ! N: number of particles (needed for T and rho)
    real(kind=dp_kind)     :: volume = 0.d0         ! V: box volume (needed for rho and W/(3V))

    ! ==================================================================================
    ! 2) MEANS: sampling counters and accumulators
    ! Sums are over samples
    ! ==================================================================================
    integer(kind=int_kind) :: n_samples = 0_int_kind  ! Number of scalar samples accumulated so far for mean calculation

    ! Energies and temperature:
    real(kind=dp_kind) :: sum_U      = 0.d0   ! sum U
    real(kind=dp_kind) :: sum_sq_U   = 0.d0   ! sum [U^2]
    real(kind=dp_kind) :: sum_K      = 0.d0   ! sum K
    real(kind=dp_kind) :: sum_sq_K   = 0.d0   ! sum [K^2]
    real(kind=dp_kind) :: sum_E      = 0.d0   ! sum E = sum(U + K)
    real(kind=dp_kind) :: sum_sq_E   = 0.d0   ! sum [E^2]
    real(kind=dp_kind) :: sum_T      = 0.d0   ! sum T, with T = 2 K / (3 N)
    real(kind=dp_kind) :: sum_sq_T   = 0.d0   ! sum [T^2]

    ! Optional pressure accumulation (requires virial W):
    logical :: virial_present = .false.              ! If .true., pressure is accumulated and virial must be provided
    real(kind=dp_kind) :: sum_P    = 0.d0            ! sum P
    real(kind=dp_kind) :: sum_sq_P = 0.d0            ! sum [P^2]

    ! Scalar observables for thermodynamic coefficients calculations
    real(kind=dp_kind) :: sum_Kinv     = 0.d0  ! sum (1/K)
    real(kind=dp_kind) :: sum_sq_Kinv  = 0.d0  ! sum [(1/K)^2]

    real(kind=dp_kind) :: sum_dU     = 0.d0    ! sum dU
    real(kind=dp_kind) :: sum_sq_dU  = 0.d0    ! sum [dU^2]
    real(kind=dp_kind) :: sum_ddU    = 0.d0    ! sum ddU
    real(kind=dp_kind) :: sum_sq_ddU = 0.d0    ! sum [ddU^2]

    real(kind=dp_kind) :: sum_dU_Kinv     = 0.d0  ! sum (dU * 1/K)
    real(kind=dp_kind) :: sum_sq_dU_Kinv  = 0.d0  ! sum [(dU * 1/K)^2]
    real(kind=dp_kind) :: sum_ddU_Kinv    = 0.d0  ! sum (ddU * 1/K)
    real(kind=dp_kind) :: sum_sq_ddU_Kinv = 0.d0  ! sum [(ddU * 1/K)^2]

    ! Fluctuation term needed in K_S: < dU^2 * (1/K) >
    real(kind=dp_kind) :: sum_dU2_Kinv     = 0.d0  ! sum ((dU^2) * 1/K)
    real(kind=dp_kind) :: sum_sq_dU2_Kinv  = 0.d0  ! sum [((dU^2) * 1/K)^2]

  end type md_means_accum

contains

  !====================================================================================================
  ! Initialize / reset MEANS (and store params locally)
  !
  ! Inputs:
  !  - params: provides N and V needed to compute derived quantities (T, rho, P)
  !  - use_virial:
  !      .true.  -> pressure will be accumulated, and virial W must be provided
  !      .false. -> pressure is ignored (less overhead)
  !====================================================================================================
  subroutine md_means_init(acc, params, use_virial)
    type(md_means_accum), intent(inout) :: acc         ! Internal type
    type(sim_params),     intent(in)    :: params
    logical,              intent(in)    :: use_virial

    ! Store run constants locally so this module does not depend on globals.
    acc%n      = params%n
    acc%volume = params%volume

    ! Sanity checks
    if (acc%n <= 0_int_kind) stop 'md_means_init(): params%n must be > 0.'
    if (acc%volume <= 0.d0)  stop 'md_means_init(): params%volume must be > 0.'

    acc%virial_present = use_virial

    call md_means_reset(acc)

  end subroutine md_means_init


  !====================================================================================================
  ! Reset counters and accumulators (keeps stored n, volume, virial_present)
  !====================================================================================================
  subroutine md_means_reset(acc)
    type(md_means_accum), intent(inout) :: acc

    acc%n_samples = 0_int_kind

    acc%sum_U    = 0.d0; acc%sum_sq_U    = 0.d0
    acc%sum_K    = 0.d0; acc%sum_sq_K    = 0.d0
    acc%sum_E    = 0.d0; acc%sum_sq_E    = 0.d0
    acc%sum_T    = 0.d0; acc%sum_sq_T    = 0.d0

    acc%sum_P    = 0.d0; acc%sum_sq_P    = 0.d0

    acc%sum_Kinv     = 0.d0; acc%sum_sq_Kinv     = 0.d0
    acc%sum_dU       = 0.d0; acc%sum_sq_dU       = 0.d0
    acc%sum_ddU      = 0.d0; acc%sum_sq_ddU      = 0.d0
    acc%sum_dU_Kinv  = 0.d0; acc%sum_sq_dU_Kinv  = 0.d0
    acc%sum_ddU_Kinv = 0.d0; acc%sum_sq_ddU_Kinv = 0.d0

    ! NEW (but now in the new naming scheme)
    acc%sum_dU2_Kinv    = 0.d0
    acc%sum_sq_dU2_Kinv = 0.d0

  end subroutine md_means_reset


  !====================================================================================================
  ! Add one instantaneous sample to MEANS. No averages are computed here yet, although the sums needed
  ! for computing averages later are what is stored
  !
  ! Inputs:
  !  - epot (U), ekin (K) : instantaneous energies
  !  - virial (optional): instantaneous virial W, only required if virial_present = .true.
  !
  ! Optional outputs:
  !  - temp_inst  : instantaneous temperature computed here: T = 2 E_kin / (3 N)
  !  - press_inst : instantaneous pressure computed here (only if virial_present)
  !
  ! Optional inputs:
  !  - d_epot (dU), dd_epot (ddU) : provided by LJ module in md_simulation.f90 (used in coefficient formulas)
  !
  ! NOTE: epot, ekin, temp_inst, virial, press_inst, d_epot, dd_epot names are kept as in the rest of the program
  !====================================================================================================
  subroutine md_means_add_sample(acc, epot, ekin, temp_inst, virial, press_inst, d_epot, dd_epot)

    type(md_means_accum), intent(inout) :: acc
    real(kind=dp_kind),   intent(in)    :: epot, ekin

    real(kind=dp_kind), intent(out), optional :: temp_inst, press_inst
    real(kind=dp_kind), intent(in),  optional :: virial
    real(kind=dp_kind), intent(in),  optional :: d_epot, dd_epot

    real(kind=dp_kind) :: E        ! E for this sample
    real(kind=dp_kind) :: T        ! T from kinetic energy
    real(kind=dp_kind) :: P        ! P if virial is enabled
    real(kind=dp_kind) :: npd      ! N as double precision
    real(kind=dp_kind) :: rho      ! rho = N/V

    real(kind=dp_kind) :: Kinv     ! 1/K for this sample
    real(kind=dp_kind) :: dU, ddU
    real(kind=dp_kind) :: dU2_Kinv ! (dU^2) * (1/K) for this sample
    real(kind=dp_kind) :: U, K

    if (acc%n <= 0_int_kind) stop 'md_means_add_sample(): accumulator not initialized (n <= 0).'
    if (acc%volume <= 0.d0)  stop 'md_means_add_sample(): invalid volume (accumulator not initialized?).'

    npd = dble(acc%n)              ! Number of particles in double precision for real valued formulas
    rho = npd / acc%volume         ! Density

    U = epot
    K = ekin
    E = U + K
    T = 2.d0 * K / (3.d0 * npd)       ! T = 2*K/(3 N)

    ! Pressure (only if enabled): P = rho*T + W/(3V)
    P = 0.d0
    if (acc%virial_present) then
      if (.not. present(virial)) stop 'md_means_add_sample(): virial expected but not provided.'
      P = rho*T + virial/(3.d0*acc%volume)
    end if

    ! Update sample count
    acc%n_samples = acc%n_samples + 1_int_kind

    ! Accumulate sums and sums of squares
    acc%sum_U    = acc%sum_U    + U; acc%sum_sq_U = acc%sum_sq_U + U*U
    acc%sum_K    = acc%sum_K    + K; acc%sum_sq_K = acc%sum_sq_K + K*K
    acc%sum_E    = acc%sum_E    + E; acc%sum_sq_E = acc%sum_sq_E + E*E
    acc%sum_T    = acc%sum_T    + T; acc%sum_sq_T = acc%sum_sq_T + T*T

    if (acc%virial_present) then
      acc%sum_P    = acc%sum_P    + P
      acc%sum_sq_P = acc%sum_sq_P + P*P
    end if

    if (K <= 0.d0) stop 'md_means_add_sample(): ekin must be > 0 to accumulate 1/ekin terms.'
    Kinv = 1.d0 / K

    acc%sum_Kinv    = acc%sum_Kinv    + Kinv
    acc%sum_sq_Kinv = acc%sum_sq_Kinv + Kinv*Kinv

    if (present(d_epot)) then
      dU = d_epot
      acc%sum_dU    = acc%sum_dU    + dU
      acc%sum_sq_dU = acc%sum_sq_dU + dU*dU

      acc%sum_dU_Kinv    = acc%sum_dU_Kinv    + dU*Kinv
      acc%sum_sq_dU_Kinv = acc%sum_sq_dU_Kinv + (dU*Kinv)*(dU*Kinv)
      
      dU2_Kinv = (dU*dU) * Kinv
      acc%sum_dU2_Kinv    = acc%sum_dU2_Kinv    + dU2_Kinv
      acc%sum_sq_dU2_Kinv = acc%sum_sq_dU2_Kinv + dU2_Kinv*dU2_Kinv
    end if

    if (present(dd_epot)) then
      ddU = dd_epot
      acc%sum_ddU    = acc%sum_ddU    + ddU
      acc%sum_sq_ddU = acc%sum_sq_ddU + ddU*ddU

      acc%sum_ddU_Kinv    = acc%sum_ddU_Kinv    + ddU*Kinv
      acc%sum_sq_ddU_Kinv = acc%sum_sq_ddU_Kinv + (ddU*Kinv)*(ddU*Kinv)
    end if

    if (present(temp_inst))  temp_inst  = T
    if (present(press_inst)) press_inst = P

  end subroutine md_means_add_sample


  !====================================================================================================
  ! Return means and standard deviations FOR ALL variables accumulated.
  !
  ! Mean(x) = 1/nsamples sum_{i=1}^{nsamples}(x)
  ! Std = sqrt( <x^2> - <x>^2 )
  !====================================================================================================
  subroutine md_means_get(acc, &
                         U_mean, K_mean, E_mean, T_mean, P_mean, &
                         Kinv_mean, dU_mean, ddU_mean, dU_Kinv_mean, dU2_Kinv_mean, ddU_Kinv_mean, &
                         U_std,  K_std,  E_std,  T_std,  P_std, &
                         Kinv_std, dU_std, ddU_std, dU_Kinv_std, dU2_Kinv_std, ddU_Kinv_std)

    type(md_means_accum), intent(in) :: acc

    !* Returns *!
    ! Means
    real(kind=dp_kind), intent(out) :: U_mean, K_mean, E_mean, T_mean, P_mean
    real(kind=dp_kind), intent(out) :: Kinv_mean, dU_mean, ddU_mean
    real(kind=dp_kind), intent(out) :: dU_Kinv_mean, dU2_Kinv_mean, ddU_Kinv_mean

    ! Standard deviations
    real(kind=dp_kind), intent(out) :: U_std, K_std, E_std, T_std, P_std
    real(kind=dp_kind), intent(out) :: Kinv_std, dU_std, ddU_std
    real(kind=dp_kind), intent(out) :: dU_Kinv_std, dU2_Kinv_std, ddU_Kinv_std

    ! Internal variables
    real(kind=dp_kind) :: inv_ns
    real(kind=dp_kind) :: U_m2, K_m2, E_m2, T_m2
    real(kind=dp_kind) :: P_m2
    real(kind=dp_kind) :: Kinv_m2, dU_m2, ddU_m2
    real(kind=dp_kind) :: dU_Kinv_m2, ddU_Kinv_m2
    real(kind=dp_kind) :: dU2_Kinv_m2

    if (acc%n_samples <= 0_int_kind) stop 'md_means_get(): no samples accumulated.'
    inv_ns = 1.d0 / dble(acc%n_samples)

    ! Means
    U_mean = acc%sum_U * inv_ns
    K_mean = acc%sum_K * inv_ns
    E_mean = acc%sum_E * inv_ns
    T_mean = acc%sum_T * inv_ns
    if (acc%virial_present) then
      P_mean = acc%sum_P * inv_ns
    else
      P_mean = 0.d0
    end if

    Kinv_mean    = acc%sum_Kinv * inv_ns
    dU_mean      = acc%sum_dU   * inv_ns
    ddU_mean     = acc%sum_ddU  * inv_ns
    dU_Kinv_mean = acc%sum_dU_Kinv * inv_ns
    ddU_Kinv_mean = acc%sum_ddU_Kinv * inv_ns
    dU2_Kinv_mean = acc%sum_dU2_Kinv * inv_ns

    ! <x^2> (moments)
    U_m2 = acc%sum_sq_U * inv_ns
    K_m2 = acc%sum_sq_K * inv_ns
    E_m2 = acc%sum_sq_E * inv_ns
    T_m2 = acc%sum_sq_T * inv_ns
    if (acc%virial_present) then
      P_m2 = acc%sum_sq_P * inv_ns
    else
      P_m2 = 0.d0
    end if

    Kinv_m2 = acc%sum_sq_Kinv * inv_ns
    dU_m2   = acc%sum_sq_dU   * inv_ns
    ddU_m2  = acc%sum_sq_ddU  * inv_ns

    dU_Kinv_m2  = acc%sum_sq_dU_Kinv  * inv_ns
    ddU_Kinv_m2 = acc%sum_sq_ddU_Kinv * inv_ns
    dU2_Kinv_m2 = acc%sum_sq_dU2_Kinv * inv_ns

    ! Standard deviations (delegated to stats_math)
    U_std = stats_std_from_moments(U_m2, U_mean)
    K_std = stats_std_from_moments(K_m2, K_mean)
    E_std = stats_std_from_moments(E_m2, E_mean)
    T_std = stats_std_from_moments(T_m2, T_mean)
    P_std = stats_std_from_moments(P_m2, P_mean)

    Kinv_std = stats_std_from_moments(Kinv_m2, Kinv_mean)
    dU_std   = stats_std_from_moments(dU_m2,   dU_mean)
    ddU_std  = stats_std_from_moments(ddU_m2,  ddU_mean)

    dU_Kinv_std  = stats_std_from_moments(dU_Kinv_m2,  dU_Kinv_mean)
    ddU_Kinv_std = stats_std_from_moments(ddU_Kinv_m2, ddU_Kinv_mean)
    dU2_Kinv_std = stats_std_from_moments(dU2_Kinv_m2, dU2_Kinv_mean)

  end subroutine md_means_get

end module md_means
