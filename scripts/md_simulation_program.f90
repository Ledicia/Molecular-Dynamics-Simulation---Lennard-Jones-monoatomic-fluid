!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! Molecular Dynamics (MD) SIMULATION
!
! Goal:
! - Read an already prepared initial configuration (positions + velocities) from outputs/rv_init.dat
! - Propagate the system in time using the velocity-Verlet algorithm
! - Sample energies every output_interval steps and store a time series
! - Compute time averages and an autocorrelation of total energy
! - Compute thermodynamic coefficients
!
! Files:
! - Input : inputs/input_simulation_parameters.txt
! - Input : outputs/rv_init.dat                         => unformatted: record1 rx,ry,rz; record2 vx,vy,vz
! - Output: outputs/one_run/instantaneous_energies.dat
! - Output: outputs/one_run/rva.dat                     => unformatted snapshots at sampling times
! - Output: outputs/one_run/corr_*.dat                  => Save correlations of epot, ekin, etot, temp, press
! - Output: outputs/one_run/corrmean_*.dat              => Save block-averaged correlation curves of epot, ekin, etot, temp, press
! - Output: outputs/one_run/md_final_results.txt        => summary
!
! Notes:
! - Everything is stored in:
!     type(sim_params) :: params   => n = N, L, V, dt, rc, etc.
!     type(sim_state)  :: state    => rx,ry,rz,vx,vy,vz,ax,ay,az
! - Lennard-Jones potential (reduced units): U(r) = 4 * (r^-12 - r^-6)
! - Velocity-Verlet (mass m = 1):
!     r(t+dt)   = r(t) + v(t)*dt + (1/2)*a(t)*dt^2
!     v(t+dt/2) = v(t) + (1/2)*a(t)*dt
!     a(t+dt)   = F(r(t+dt))  (from LJ forces)
!     v(t+dt)   = v(t+dt/2) + (1/2)*a(t+dt)*dt
!
! - Trajectory post-processing:
!   * RDF g(r): needs MANY snapshots of positions r(t). Wrapped positions are OK (MIC distances).
!   * VACF: needs velocities v(t).
!   * MSD/DCM: needs UNWRAPPED positions ru(t), i.e., positions continuous across PBC.
!
!   To build unwrapped coordinates we integrate the MIC displacement each MD step:
!     dx      = x_new_wrapped - x_old_wrapped
!     dx_mic  = dx - L * nint(dx/L)
!     x_unwrapped += dx_mic
!
!   This assumes dt is small enough that a particle does not move more than L/2 per step.
!******************************************************************************************************
program md_simulation

  use define_precision,    only: dp_kind, int_kind
  use md_types,            only: sim_params, sim_state, init_state
  use read_input_files,    only: read_simulation_parameters
  use thermodynamic_coefs, only: thermodynamic_coefs_results, thermodynamic_coefs_compute

  ! Scalar statistics (means, correlations and correlation means)
  ! - ALL scalar statistics needed by MD post-processing are accumulated
  !   inside dedicated MD modules to avoid duplicated bookkeeping in the main program.
  !
  ! Pattern used for the 3 stats paths:
  !   1) INIT      -> allocate/reset internal accumulators
  !   2) ADD       -> push one sampled value (only inside sampling if)
  !   3) COMPUTE   -> finalize derived arrays (only for correlations / corr_means)
  !   4) GET       -> retrieve results to local arrays/variables
  !   5) FREE      -> release internal memory (optional, but recommended)
  use md_means,            only: md_means_accum, md_means_init, md_means_add_sample, md_means_get
  use md_correlations,     only: md_corr_accum,  md_corr_init,  md_corr_add_sample,  md_corr_compute, md_corr_get_all, md_corr_free, &
                                 md_corr_cm_init, md_corr_cm_compute, md_corr_cm_get_all, md_corr_cm_free, md_corr_cm_max_lag

  use lj_potential_energy, only: compute_lj_potential_energy
  use verlet,              only: verlet_step   ! velocity-Verlet integrator

  implicit none

  ! -----------------------
  ! Simulation containers
  ! -----------------------
  type(sim_params) :: params        ! Run configuration (N, L, V, dt, rc, ...)
  type(sim_state)  :: state         ! Dynamic state (r, v, a arrays)
  
  ! -----------------------
  ! Control parameters read from (some parameters are kept for compatibility with input format, but are not used here)
  !------------------------
  ! Output_steps usually =! 1 => Usually we do not store one value per MD step since:
  !   - Consecutive MD steps are strongly correlated and would add little statistical information
  !   - It would also increase memory usage and noise in correlation functions
  !   - Instead, sample the system every output_interval steps: step = output_interval, 2*output_interval, ...
  ! Therefore:
  !   n_samples         = total_steps / output_interval                 <= total possible sampling instants 
  !   n_samples_warmup  = warmup_steps / output_interval                <= warmup sampling instants (discarded)
  !   n_samples_md      = (total_steps - warmup_steps)/output_interval  <= production samples
  ! -----------------------
  integer(kind=int_kind) :: total_steps
  integer(kind=int_kind) :: warmup_steps
  integer(kind=int_kind) :: output_interval
  real(kind=dp_kind)     :: rc_over_L
  real(kind=dp_kind)     :: target_total_energy

  ! -----------------------
  ! MD loop counters
  ! -----------------------
  integer(kind=int_kind) :: step    ! Time step
  integer(kind=int_kind) :: i_part  ! For the particle loop

  ! -----------------------
  ! Instantaneous observables
  ! -----------------------
  real(kind=dp_kind) :: time
  real(kind=dp_kind) :: epot, ekin, etot
  real(kind=dp_kind) :: temp_inst          ! Temperature
  real(kind=dp_kind) :: virial_inst        ! Virial  W = sum_{i<j} r_ij · F_ij
  real(kind=dp_kind) :: press_inst         ! Pressure

  ! LJ radial derivative sums (returned by LJ module):
  !   d_epot  = sum_{i<j} [ r_ij * dU/dr_ij ]
  !   dd_epot = sum_{i<j} [ r_ij^2 * d^2U/dr_ij^2 ]   (see lj_potential_energy.f90)
  real(kind=dp_kind) :: d_epot, dd_epot

  ! -----------------------
  ! Sampling counters
  ! -----------------------
  integer(kind=int_kind) :: num_samples   ! Number of stored samples (one sample every output_interval steps)

  ! -----------------------
  ! MEANS outputs (with std)  returned by md_means
  ! -----------------------
  real(kind=dp_kind) :: epot_mean, ekin_mean, etot_mean, temp_mean, press_mean
  real(kind=dp_kind) :: epot_std,  ekin_std,  etot_std,  temp_std,  press_std

  ! Means needed by thermodynamic coefficient formulas returned by md_means:
  real(kind=dp_kind) :: ekinInv_mean                               ! < (1/ekin) >
  real(kind=dp_kind) :: d_epot_mean, dd_epot_mean                  ! < d_epot > & < dd_epot >
  real(kind=dp_kind) :: d_epot_ekinInv_mean, dd_epot_ekinInv_mean  ! < d_epot_ekinInv > & < dd_epot_ekinInv >
  real(kind=dp_kind) :: d_epot2_ekinInv_mean                       ! < d_epot^2 / ekin >

  ! Optional std of extra coefficient inputs (useful for diagnostics / uncertainty analysis) returned by md_means
  real(kind=dp_kind) :: ekinInv_std
  real(kind=dp_kind) :: d_epot_std, dd_epot_std
  real(kind=dp_kind) :: d_epot_ekinInv_std, dd_epot_ekinInv_std
  real(kind=dp_kind) :: d_epot2_ekinInv_std
  
  ! -----------------------
  ! THERMODYNAMIC COEFICIENTS outputs returned by thermodynamic_coefs.f90
  ! -----------------------
  type(thermodynamic_coefs_results) :: thermodynamic_coefs  ! Stores thermodynamic coefcientes

  ! -----------------------
  ! Correlation settings (multi-observable autocorrelations)
  ! -----------------------
  integer(kind=int_kind) :: n_samples_md, corr_max_lag     ! n_samples_md is the number of samples stored during simulation and also correspond to the number of samples available for correlations
  integer(kind=int_kind) :: lag_max_out
  logical :: corr_enabled                                  ! Flag to know if correlations were actually initialized

  ! Variables for writing full-series correlations and normalized correlations returned by md_correlations
  !    corr_vals_*  & corrn_vals_*: Autocorrelation arrays for sampled observables at lags L = 0,1,2,...,max_lag
  !
  ! E.g. For epot: (C(0), C(1), ..., C(max_lag))
  !      1) Correlations : C(L) = (1/n_samples_md-L)*sum_{k=1,...,n_samples_md}(E_pot(k)*E_pot(k+L)) if is not centralized
  !      2) Covariance   : C(L) = (1/n_samples_md-L)*sum_{k=1,...,n_samples_md}[(E_pot(k) - <E_pot(k)>) * (E_pot(k+L) - <E_pot(k+L)>) 
  !
  ! NOTE: The exact estimator (with or without mean subtraction) is implemented inside md_correlations.
  real(kind=dp_kind), allocatable :: corr_vals_epot(:),  corrn_vals_epot(:)
  real(kind=dp_kind), allocatable :: corr_vals_ekin(:),  corrn_vals_ekin(:)
  real(kind=dp_kind), allocatable :: corr_vals_etot(:),  corrn_vals_etot(:)
  real(kind=dp_kind), allocatable :: corr_vals_temp(:),  corrn_vals_temp(:)
  real(kind=dp_kind), allocatable :: corr_vals_press(:), corrn_vals_press(:)

  ! -----------------------
  ! Correlation means settings (block averaging)
  !    corr_mean_* & corrn_mean_*: Means of the correlations
  !
  ! The whole time series A(1), ..., A(n_samples_md) is divided in B blocks o lenght n_block_samples:
  !  1) Split time series into B blocks    : A_1(1), ..., A_1(n_block_samples), ..., A_B((B-1)*n_block_samples+1), ..., A_B(n_samples_md)  => A_b(1), ..., A_b(n_block_samples),  b = 1...B
  !  2) Compute correlations in each block : (C_1(0), ..., C_1(lag)), ..., (C_B(0), ..., C_B(lag))                                         => C_b(lag) ,  lag = 0...lag_max; b = 1...B
  !  3) Block-averaged correlation         : ((1/B)*sum_{b=1...lag}(C_b(0)), ..., (1/B)*sum_{b=1...lag}(C_b(lag)))                         => <C(lag)>_blocks = (1/B) * sum_{b=1..B} C_b(lag)
  ! -----------------------
  integer(kind=int_kind) :: n_blocks
  real(kind=dp_kind), allocatable :: corr_mean_epot(:),  corrn_mean_epot(:)
  real(kind=dp_kind), allocatable :: corr_mean_ekin(:),  corrn_mean_ekin(:)
  real(kind=dp_kind), allocatable :: corr_mean_etot(:),  corrn_mean_etot(:)
  real(kind=dp_kind), allocatable :: corr_mean_temp(:),  corrn_mean_temp(:)
  real(kind=dp_kind), allocatable :: corr_mean_press(:), corrn_mean_press(:)

  ! -----------------------
  ! Stats accumulators (stateful objects)
  ! -----------------------
  type(md_means_accum) :: means_acc
  type(md_corr_accum)  :: corr_acc

  ! -----------------------
  ! I/O units
  !  - iu_out:  Unit number for writing the instantaneous energy time series (outputs/one_run/instantaneous_energies.dat)
  !  - iu_sum:  Unit number for appending the final summary of averages and thermodynamic coefficients (outputs/one_run/md_final_results.txt)
  !  - iu_rva:  Unit number for writing the trajectory file (outputs/one_run/rva.dat) containing wrapped positions, unwrapped positions, velocities and accelerations
  !  - ios: I/O status flag returned by OPEN statements (=0: successful, !=0: error)
  ! -----------------------
  integer :: iu_out, iu_sum, iu_rva, ios

  ! ===================================================================================
  ! TRAJECTORY POST-PROCESSING (RDF / VACF / MSD)
  ! ===================================================================================
  real(kind=dp_kind), allocatable :: rux(:), ruy(:), ruz(:)               ! Unwrapped coordinates (required for MSD/DCM)
  real(kind=dp_kind), allocatable :: rx_prev(:), ry_prev(:), rz_prev(:)   ! Previous wrapped positions (to compute step displacement)
  real(kind=dp_kind)              :: dx, dy, dz                           ! Step displacement and MIC-corrected displacement
  integer(kind=int_kind)          :: n_snapshots_expected                 ! Expected number of snapshots written to rva.dat

  ! ===================================================================================
  ! 1) Read simulation parameters from file 
  !    - Updates variables inside sim_params (defined in md_types.f90)
  ! ===================================================================================
  call read_simulation_parameters('inputs/input_simulation_parameters.txt', &
                                  params, total_steps, output_interval, warmup_steps, &
                                  rc_over_L, target_total_energy)

  ! ===================================================================================
  ! 2) Allocate and initialize the simulation state in memory
  !    - Save memory for sim_state variables (defined in md_types.f90) and initialize its values to 0
  ! ===================================================================================
  call init_state(params, state)

  ! ===================================================================================
  ! 3) Read initial configuration generated using initial_config.f90: outputs/rv_init.dat
  ! ===================================================================================
  call read_rv_init('outputs/rv_init.dat', params, state)

  ! ===================================================================================
  ! 4) Allocate unwrapped tracking arrays and initialize ru(t=0) = r(t=0)
  ! ===================================================================================
  allocate(rux(params%n), ruy(params%n), ruz(params%n))
  allocate(rx_prev(params%n), ry_prev(params%n), rz_prev(params%n))

  rux(1:params%n) = state%rx(1:params%n)
  ruy(1:params%n) = state%ry(1:params%n)
  ruz(1:params%n) = state%rz(1:params%n)

  ! ===================================================================================
  ! 5) Initial forces and energies at t = 0
  ! ===================================================================================
  call compute_lj_potential_energy(params, state, epot, d_epot, dd_epot)

  ekin = 0.5d0 * sum(state%vx(1:params%n)*state%vx(1:params%n) + &
                     state%vy(1:params%n)*state%vy(1:params%n) + &
                     state%vz(1:params%n)*state%vz(1:params%n))

  etot = epot + ekin
  time = 0.d0

  ! ===================================================================================
  ! 6) Open rva.dat ONCE to save unformatted trajectory snapshots: outputs/one_run/rva.dat
  ! ===================================================================================
  iu_rva = 31
  open(iu_rva, file='outputs/one_run/rva.dat', form='unformatted', status='replace', action='write', iostat=ios)
  if (ios /= 0) stop 'md_simulation: cannot open outputs/one_run/rva.dat'

  ! Number of snapshots written by the sampling condition used inside the MD loop:
  !  - if (step > warmup_steps) .and. (mod(step, output_interval) == 0)
  n_snapshots_expected = (total_steps / output_interval) - (warmup_steps / output_interval)
  if (n_snapshots_expected < 0_int_kind) n_snapshots_expected = 0_int_kind

  write(iu_rva) params%n, params%box_length, params%dt, output_interval, n_snapshots_expected

  ! ===================================================================================
  ! 7) STATS INIT: MEANS (md_means.f90) & CORRELATIONS (md_correations.f90)
  !    - n_samples_md = (total_steps - warmup_steps) / output_interval: 
  !        - REAL number of samples stored after warmup
  !        - Sample only when: (step > warmup_steps) and mod(step, output_interval)==0
  !    - n_samples = total_steps/output_interval:
  !        - Number of samples comsidering also the warmup steps, which are discarded because the system
  !          is not yet in a stationary (equilibrated) regime.
  !        - Warmup steps are used to generate initial config rv_init.dat and MD simulation steps are (total_steps - warmup_steps)
  !
  !    - Use corr_enabled to decide if we can push samples and compute correlations safely
  ! ===================================================================================  
  ! MEANS !
  call md_means_init(means_acc, params, use_virial=.true.)
  num_samples   = 0_int_kind
  corr_enabled  = .false.
  
  ! CORRELATIONS !
  n_samples_md = (total_steps / output_interval) - (warmup_steps / output_interval)
  if (n_samples_md < 0_int_kind) n_samples_md = 0_int_kind                          ! Safety check in case warmup_steps > total_steps (shouldnt happen)

  if (n_samples_md >= 2_int_kind) then
    corr_max_lag = min(1000_int_kind, n_samples_md - 1_int_kind) ! lag <= (n_samples_md - 1) and also max lag should not be greater than 1000
    corr_max_lag = min(corr_max_lag, n_samples_md / 2_int_kind)  ! Avoid having lags larger than ~ half the series

    call md_corr_init(corr_acc, params, n_samples_md, corr_max_lag, want_temp=.true., want_press=.true.)
    corr_enabled = .true.
  else
    corr_max_lag = 0_int_kind
  end if

  ! Energies time series (text) written here (md_simulation) !
  iu_out = 30
  open(iu_out, file='outputs/one_run/instantaneous_energies.dat', status='replace', action='write', iostat=ios)
  if (ios /= 0) stop 'md_simulation: cannot open outputs/one_run/instantaneous_energies.dat'
  write(iu_out,'(a)') '# time   epot   ekin   etot   T   P'


  ! ===================================================================================
  ! 8) MD loop to advance system in time
  ! ===================================================================================
  do step = 1, total_steps

    ! Save wrapped positions BEFORE the integrator step (needed to compute the MIC displacement and update ru(t))
    rx_prev(1:params%n) = state%rx(1:params%n)
    ry_prev(1:params%n) = state%ry(1:params%n)
    rz_prev(1:params%n) = state%rz(1:params%n)

    ! ------------------------------------------------------------
    ! One velocity-Verlet microscopic step:
    !  - Updates the microscopic state IN PLACE:
    !      * state%rx, state%ry, state%rz  -> positions at t + dt
    !      * state%vx, state%vy, state%vz  -> velocities at t + dt
    !      * state%ax, state%ay, state%az  -> accelerations at t + dt
    !
    !  - Recomputes interparticle forces from the new positions and returns instantaneous observables corresponding to the updated state:
    !      * epot    : potential energy at t + dt
    !      * ekin    : kinetic energy at t + dt
    !      * d_epot  : configurational virial term (r * dU/dr), used for pressure and equation of state
    !      * dd_epot : second-order configurational derivative, used in response functions (e.g. bulk modulus)
    !
    ! Note: The previous values of epot, ekin, d_epot and dd_epot are overwritten at each call.
    ! ------------------------------------------------------------
    call verlet_step(params, state, epot, ekin, d_epot, dd_epot)

    ! ------------------------------------------------------------
    ! Update unwrapped coordinates for each particle in the system(required for MSD/DCM).
    !  - r(t)  : wrapped position in [0, L), used for forces and energies
    !  - ru(t) : unwrapped position, continuous in time, used for MSD/diffusion
    !
    ! For each particle:
    !   dx = x_new_wrapped - x_old_wrapped
    !   dx_mic = dx - L * nint(dx/L)   (MIC applied to displacement)
    !   x_unwrapped += dx_mic
    !
    ! Example (1D, L = 10):
    !   x_old = 9.9,  x_new = 0.1  -> dx = -9.8
    !   dx_mic = -9.8 - 10*nint(-0.98) = -9.8 - 10*(-1) = +0.2
    !   xu_new = xu_old + 0.2  -> continuous trajectory across the boundary.
    ! ------------------------------------------------------------
    do i_part = 1, params%n

      dx = state%rx(i_part) - rx_prev(i_part)                          ! dx = x_new_wrapped - x_old_wrapped
      dy = state%ry(i_part) - ry_prev(i_part)
      dz = state%rz(i_part) - rz_prev(i_part)

      dx = dx - params%box_length * dnint(dx * params%inv_box_length)  ! dx_mic = dx - L * nint(dx/L) (MIC: Minimum Image Convention <=> see geometry_pbc.f90 for a better explanation)
      dy = dy - params%box_length * dnint(dy * params%inv_box_length)
      dz = dz - params%box_length * dnint(dz * params%inv_box_length)

      rux(i_part) = rux(i_part) + dx                                   ! xu_new = xu_old + dx_mic
      ruy(i_part) = ruy(i_part) + dy
      ruz(i_part) = ruz(i_part) + dz

    end do

    etot = epot + ekin
    time = time + params%dt

    ! ------------------------------------------------------------
    ! Sampling: warmup is NOT sampled since is just for equilibration
    ! ------------------------------------------------------------
    if ( (step > warmup_steps) .and. (mod(step, output_interval) == 0) ) then

      num_samples = num_samples + 1_int_kind

      ! Virial: W = sum r·F = sum r *(-dU/dr) = - d_epot
      virial_inst = -d_epot


      ! MEANS ADD:
      !  - temp_inst and press_inst are returned (computed internally from ekin and virial)
      call md_means_add_sample(means_acc, epot, ekin, temp_inst=temp_inst, press_inst=press_inst, &
                               virial=virial_inst, d_epot=d_epot, dd_epot=dd_epot)

      write(iu_out,'(1pe13.6,5(2x,1pe13.6))') time, epot, ekin, etot, temp_inst, press_inst

      ! CORRELATIONS ADD:
      ! - Only if correlations were initialized (corr_enabled=.true.)
      ! - We push ONE value per sampling time (same cadence as the means and output series)
      if (corr_enabled) then
        call md_corr_add_sample(corr_acc, epot, ekin, temp=temp_inst, press=press_inst)
      end if

      ! Trajectory snapshot to outputs/one_run/rva.dat
      write(iu_rva) state%rx(1:params%n),  state%ry(1:params%n),  state%rz(1:params%n)
      write(iu_rva) rux(1:params%n),       ruy(1:params%n),       ruz(1:params%n)
      write(iu_rva) state%vx(1:params%n),  state%vy(1:params%n),  state%vz(1:params%n)
      write(iu_rva) state%ax(1:params%n),  state%ay(1:params%n),  state%az(1:params%n)

    end if

  end do

  close(iu_out)
  close(iu_rva)

  ! ===================================================================================
  ! 9) STATS GET (MEANS)
  ! ===================================================================================
  if (num_samples <= 0_int_kind) stop 'md_simulation: no samples were taken (check warmup_steps/output_interval).'

  call md_means_get(means_acc, &
     epot_mean, ekin_mean, etot_mean, temp_mean, press_mean, &
     ekinInv_mean, d_epot_mean, dd_epot_mean, d_epot_ekinInv_mean, d_epot2_ekinInv_mean, dd_epot_ekinInv_mean, &
     epot_std,  ekin_std,  etot_std,  temp_std,  press_std, &
     ekinInv_std, d_epot_std, dd_epot_std, d_epot_ekinInv_std, d_epot2_ekinInv_std, dd_epot_ekinInv_std )
        

  ! ===================================================================================
  ! 10) Thermodynamic coefficients (same algebra as original project)
  ! ===================================================================================
  call thermodynamic_coefs_compute(params, ekin_mean, press_mean, ekinInv_mean, &
                                 d_epot_mean, dd_epot_mean, d_epot_ekinInv_mean, d_epot2_ekinInv_mean, &
                                 thermodynamic_coefs)
                                 

  ! ===================================================================================
  ! 11) CORRELATIONS (full-series autocorrelations)
  ! ===================================================================================
  if (corr_enabled) then

    ! IMPORTANT:
    ! - Here we can choose raw autocorrelation or centered autocorrelation:
    !     centered=.false. -> C(L) = <A(k) A(k+L)>
    !     centered=.true.  -> C(L) = <(A(k)-<A>)(A(k+L)-<A>)>  (autocovariance)
    !
    ! Centered is often preferred for fluctuations, but raw is OK if you want to match older outputs.
    call md_corr_compute(corr_acc, centered=.true.)

    allocate(corr_vals_epot(0:corr_max_lag),  corrn_vals_epot(0:corr_max_lag))
    allocate(corr_vals_ekin(0:corr_max_lag),  corrn_vals_ekin(0:corr_max_lag))
    allocate(corr_vals_etot(0:corr_max_lag),  corrn_vals_etot(0:corr_max_lag))
    allocate(corr_vals_temp(0:corr_max_lag),  corrn_vals_temp(0:corr_max_lag))
    allocate(corr_vals_press(0:corr_max_lag), corrn_vals_press(0:corr_max_lag))

    call md_corr_get_all(corr_acc, &
                         corr_vals_epot,  corrn_vals_epot, &
                         corr_vals_ekin,  corrn_vals_ekin, &
                         corr_vals_etot,  corrn_vals_etot, &
                         corr_vals_temp,  corrn_vals_temp, &
                         corr_vals_press, corrn_vals_press)

    lag_max_out = corr_max_lag

    ! -----------------------
    ! Write full-series correlations (one file per observable)
    ! -----------------------
    call write_corr_file('outputs/one_run/corr_epot.dat',  lag_max_out, corr_vals_epot,  corrn_vals_epot)
    call write_corr_file('outputs/one_run/corr_ekin.dat',  lag_max_out, corr_vals_ekin,  corrn_vals_ekin)
    call write_corr_file('outputs/one_run/corr_etot.dat',  lag_max_out, corr_vals_etot,  corrn_vals_etot)
    call write_corr_file('outputs/one_run/corr_temp.dat',  lag_max_out, corr_vals_temp,  corrn_vals_temp)
    call write_corr_file('outputs/one_run/corr_press.dat', lag_max_out, corr_vals_press, corrn_vals_press)

    deallocate(corr_vals_epot,  corrn_vals_epot)
    deallocate(corr_vals_ekin,  corrn_vals_ekin)
    deallocate(corr_vals_etot,  corrn_vals_etot)
    deallocate(corr_vals_temp,  corrn_vals_temp)
    deallocate(corr_vals_press, corrn_vals_press)

    ! ===================================================================================
    ! 11b) CORR_MEANS (block-averaged correlation curves)
    !  - Computed from the same stored series inside md_correlations.
    ! ===================================================================================

    ! Conservative default:
    ! - up to 5 blocks
    ! - each block must have at least (corr_max_lag + 1) samples
    !
    ! NOTE:
    ! - num_samples is the REAL number of samples pushed into corr_acc (after warmup).
    ! - block_len = floor(num_samples / n_blocks) inside md_corr_cm_compute()
    ! - Requirement: block_len >= corr_max_lag + 1  <=>  (num_samples / n_blocks) >= corr_max_lag + 1
    !
    ! We choose:
    !   n_blocks = min(5, floor(num_samples / (corr_max_lag+1)))
    n_blocks = num_samples / (corr_max_lag + 1_int_kind)
    n_blocks = min(5_int_kind, n_blocks)

    if (n_blocks >= 1_int_kind) then

      call md_corr_cm_init(corr_acc, n_blocks)

      ! Same centered/raw choice as above (must be consistent if you want comparable shapes)
      !
      ! NOTE ABOUT INTERFACE:
      ! - If your md_corr_cm_compute() does NOT accept centered, change this call to:
      !     call md_corr_cm_compute(corr_acc)
      ! - If it DOES accept it (as you planned), keep it as is.
      call md_corr_cm_compute(corr_acc, centered=.true.)

      ! Note:
      ! - The maximum lag for correlation means is stored internally (it equals the correlations lag)
      allocate(corr_mean_epot(0:corr_max_lag),  corrn_mean_epot(0:corr_max_lag))
      allocate(corr_mean_ekin(0:corr_max_lag),  corrn_mean_ekin(0:corr_max_lag))
      allocate(corr_mean_etot(0:corr_max_lag),  corrn_mean_etot(0:corr_max_lag))
      allocate(corr_mean_temp(0:corr_max_lag),  corrn_mean_temp(0:corr_max_lag))
      allocate(corr_mean_press(0:corr_max_lag), corrn_mean_press(0:corr_max_lag))

      call md_corr_cm_get_all(corr_acc, &
                              corr_mean_epot,  corrn_mean_epot, &
                              corr_mean_ekin,  corrn_mean_ekin, &
                              corr_mean_etot,  corrn_mean_etot, &
                              corr_mean_temp,  corrn_mean_temp, &
                              corr_mean_press, corrn_mean_press)

      call write_corrmean_file('outputs/one_run/corrmean_epot.dat',  corr_max_lag, corr_mean_epot,  corrn_mean_epot)
      call write_corrmean_file('outputs/one_run/corrmean_ekin.dat',  corr_max_lag, corr_mean_ekin,  corrn_mean_ekin)
      call write_corrmean_file('outputs/one_run/corrmean_etot.dat',  corr_max_lag, corr_mean_etot,  corrn_mean_etot)
      call write_corrmean_file('outputs/one_run/corrmean_temp.dat',  corr_max_lag, corr_mean_temp,  corrn_mean_temp)
      call write_corrmean_file('outputs/one_run/corrmean_press.dat', corr_max_lag, corr_mean_press, corrn_mean_press)

      deallocate(corr_mean_epot,  corrn_mean_epot)
      deallocate(corr_mean_ekin,  corrn_mean_ekin)
      deallocate(corr_mean_etot,  corrn_mean_etot)
      deallocate(corr_mean_temp,  corrn_mean_temp)
      deallocate(corr_mean_press, corrn_mean_press)

      ! NOTE:
      ! - Not strictly necessary because md_corr_free() calls md_corr_cm_free() internally,
      !   but it is OK to keep it if you want explicit separation.
      call md_corr_cm_free(corr_acc)

    end if

    call md_corr_free(corr_acc)

  end if

  ! ===================================================================================
  ! 12) Write final summary (append)
  ! ===================================================================================
  iu_sum = 10
  open(iu_sum, file='outputs/one_run/md_final_results.txt', access='append', action='write', iostat=ios)
  if (ios /= 0) stop 'md_simulation: cannot open outputs/one_run/md_final_results.txt'

  write(iu_sum,'(a)') '************** MD PRODUCTION RESULTS **************'
  write(iu_sum,'(a,1x,i8)')         'num_particles:', params%n
  write(iu_sum,'(a,1x,i8)')         'num_cells:',     params%num_cells
  write(iu_sum,'(a,1x,1pe19.12)')   'box_length:',    params%box_length
  write(iu_sum,'(a,1x,1pe19.12)')   'volume:',        params%volume
  write(iu_sum,'(a,1x,1pe19.12)')   'density:',       dble(params%n)/params%volume
  write(iu_sum,'(a,1x,1pe19.12)')   'time_step:',     params%dt
  write(iu_sum,'(a,1x,i8)')         'output_interval:', output_interval
  write(iu_sum,'(a,1x,i10)')        'total_steps:',   total_steps
  write(iu_sum,'(a,1x,i10)')        'warmup_steps:',  warmup_steps
  write(iu_sum,'(a)')               '-------------------- Averages --------------------'
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)')   '<Epot>:', epot_mean, 'std:', epot_std
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)')   '<Ekin>:', ekin_mean, 'std:', ekin_std
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)')   '<Etot>:', etot_mean, 'std:', etot_std
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)')   '<T>   :', temp_mean, 'std:', temp_std
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)')   '<P>   :', press_mean,'std:', press_std
  write(iu_sum,'(a)')               '-------------- Thermodynamic coefficients --------------'
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Temperature:', thermodynamic_coefs%temperature, 'Pressure:', thermodynamic_coefs%pressure
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Ca_v:', thermodynamic_coefs%Ca_v, 'Ce_v:', thermodynamic_coefs%Ce_v
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Ca_p:', thermodynamic_coefs%Ca_p, 'Ce_p:', thermodynamic_coefs%Ce_p
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'kappa_S:',  thermodynamic_coefs%K_S_inv, 'kappa_T:',  thermodynamic_coefs%K_T_inv, 'Gamma:', thermodynamic_coefs%gamma
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Alpha_E1:', thermodynamic_coefs%alpha_E1, 'Alpha_E2:', thermodynamic_coefs%alpha_E2
  write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Alpha_S:', thermodynamic_coefs%alpha_S, 'Alpha_P:', thermodynamic_coefs%alpha_P
  write(iu_sum,'(a)')               '--------------------------------------------------------'
  write(iu_sum,*)
  close(iu_sum)

  deallocate(rux, ruy, ruz)
  deallocate(rx_prev, ry_prev, rz_prev)

contains

  ! ------------------------------------------------------------
  ! Read initial configuration snapshot (binary, unformatted)
  !  - File outputs/rv_init.dat contains:
  !     Record 1: rx, ry, rz  (particle positions)
  !     Record 2: vx, vy, vz  (particle velocities)
  ! ------------------------------------------------------------
  subroutine read_rv_init(filename, p, s)
    character(len=*), intent(in) :: filename
    type(sim_params), intent(in) :: p
    type(sim_state),  intent(inout) :: s

    integer :: iu, ios

    iu = 21
    open(iu, file=filename, form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'read_rv_init(): cannot open rv_init file.'

    read(iu) s%rx(1:p%n), s%ry(1:p%n), s%rz(1:p%n)
    read(iu) s%vx(1:p%n), s%vy(1:p%n), s%vz(1:p%n)

    close(iu)
  end subroutine read_rv_init


  ! ------------------------------------------------------------
  ! INTERNAL: write one correlation file (raw + normalized).
  ! ------------------------------------------------------------
  subroutine write_corr_file(filename, lag_max, corr_in, corrn_in)
    character(len=*), intent(in) :: filename
    integer(kind=int_kind), intent(in) :: lag_max
    real(kind=dp_kind), intent(in) :: corr_in(0:), corrn_in(0:)

    integer :: iu, ios
    integer(kind=int_kind) :: lag

    iu = 50
    open(iu, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'write_corr_file(): cannot open output file.'

    write(iu,'(a)') '# lag   C(lag)   C_norm(lag)'
    do lag = 0, lag_max
      write(iu,'(i8,2(2x,1pe19.12))') lag, corr_in(lag), corrn_in(lag)
    end do
    close(iu)
  end subroutine write_corr_file


  ! ------------------------------------------------------------
  ! INTERNAL: write one block-averaged correlation file.
  ! ------------------------------------------------------------
  subroutine write_corrmean_file(filename, lag_max, corr_in, corrn_in)
    character(len=*), intent(in) :: filename
    integer(kind=int_kind), intent(in) :: lag_max
    real(kind=dp_kind), intent(in) :: corr_in(0:), corrn_in(0:)

    integer :: iu, ios
    integer(kind=int_kind) :: lag

    iu = 51
    open(iu, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'write_corrmean_file(): cannot open output file.'

    write(iu,'(a)') '# lag   <C(lag)>_blocks   <C_norm(lag)>_blocks'
    do lag = 0, lag_max
      write(iu,'(i8,2(2x,1pe19.12))') lag, corr_in(lag), corrn_in(lag)
    end do
    close(iu)
  end subroutine write_corrmean_file

end program md_simulation