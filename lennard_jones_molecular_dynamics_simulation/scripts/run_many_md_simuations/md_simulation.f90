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
! - Input : outputs/rv_init.dat => unformatted: record1 rx,ry,rz; record2 vx,vy,vz
! - Output: outputs/instantaneous_energies.dat
! - Output: outputs/rva.dat => unformatted snapshots at sampling times
! - Output: outputs/corr_*.dat => Save correlations of epot, ekin, etot, temp, press
! - Output: outputs/corrmean_*.dat => Save block-averaged correlation curves of epot, ekin, etot, temp, press
! - Output: md_final_results.txt => summary
!
! Notes:
! - Everything is stored in:
!   type(sim_params) :: params => n = N, L, V, dt, rc, etc.
!   type(sim_state)  :: state  => rx,ry,rz,vx,vy,vz,ax,ay,az
!
! - Lennard-Jones potential (reduced units): U(r) = 4 * (r^-12 - r^-6)
!
! - Velocity-Verlet (mass m = 1):
!   r(t+dt)   = r(t) + v(t)*dt + (1/2)*a(t)*dt^2
!   v(t+dt/2) = v(t) + (1/2)*a(t)*dt
!   a(t+dt)   = F(r(t+dt)) (from LJ forces)
!   v(t+dt)   = v(t+dt/2) + (1/2)*a(t+dt)*dt
!
! - Trajectory post-processing:
!   * RDF g(r): needs MANY snapshots of positions r(t). Wrapped positions are OK (MIC distances).
!   * VACF   : needs velocities v(t).
!   * MSD/DCM: needs UNWRAPPED positions ru(t), i.e., positions continuous across PBC.
!
!   To build unwrapped coordinates we integrate the MIC displacement each MD step:
!     dx      = x_new_wrapped - x_old_wrapped
!     dx_mic  = dx - L * nint(dx/L)
!     x_unwrapped += dx_mic
!
!   This assumes dt is small enough that a particle does not move more than L/2 per step.
!******************************************************************************************************

module md_simulation

    use define_precision, only: dp_kind, int_kind
    use md_types,         only: sim_params, sim_state, init_state
    use read_input_files, only: read_simulation_parameters

    ! Scalar statistics (means, correlations and correlation means)
    ! - ALL scalar statistics needed by MD post-processing are accumulated inside dedicated MD modules
    !   to avoid duplicated bookkeeping in the main program.
    !
    ! Pattern used for the 3 stats paths:
    ! 1) INIT    -> allocate/reset internal accumulators
    ! 2) ADD     -> push one sampled value (only inside sampling if)
    ! 3) COMPUTE -> finalize derived arrays (only for correlations / corr_means)
    ! 4) GET     -> retrieve results to local arrays/variables
    ! 5) FREE    -> release internal memory (optional, but recommended)

    use md_means, only: md_means_accum, md_means_init, md_means_add_sample, md_means_get
    use md_correlations, only: md_corr_accum, md_corr_init, md_corr_add_sample, md_corr_compute, &
                               md_corr_get_all, md_corr_free, &
                               md_corr_cm_init, md_corr_cm_compute, md_corr_cm_get_all, md_corr_cm_free, &
                               md_corr_cm_max_lag
    use lj_potential_energy, only: compute_lj_potential_energy
    use verlet,              only: verlet_step

    implicit none
    private

    public :: run_md_simulation

contains

    ! ===================================================================================
    ! ierr   : To keep track of things that should not happen
    ! out_dir: output directory for THIS run (created beforehand by simulation_results)
    ! ===================================================================================
    subroutine run_md_simulation(out_dir, ierr)

        implicit none

        character(len=*),   intent(in)  :: out_dir
        integer(kind=int_kind), intent(out) :: ierr

        ! -----------------------
        ! Simulation containers
        ! -----------------------
        type(sim_params) :: params   ! Run configuration (N, L, V, dt, rc, ...)
        type(sim_state)  :: state    ! Dynamic state (r, v, a arrays)

        ! -----------------------
        ! Control parameters read from input
        ! -----------------------
        integer(kind=int_kind) :: total_steps
        integer(kind=int_kind) :: warmup_steps
        integer(kind=int_kind) :: output_interval
        real(kind=dp_kind)     :: rc_over_L
        real(kind=dp_kind)     :: target_total_energy

        ! -----------------------
        ! MD loop counters
        ! -----------------------
        integer(kind=int_kind) :: step
        integer(kind=int_kind) :: i_part

        ! -----------------------
        ! Instantaneous observables
        ! -----------------------
        real(kind=dp_kind) :: time
        real(kind=dp_kind) :: epot, ekin, etot
        real(kind=dp_kind) :: temp_inst
        real(kind=dp_kind) :: virial_inst
        real(kind=dp_kind) :: press_inst

        ! LJ radial derivative sums (returned by LJ module):
        ! d_epot  = sum_{i<j} [ r_ij * dU/dr_ij ]
        ! dd_epot = sum_{i<j} [ r_ij^2 * d^2U/dr_ij^2 ]
        real(kind=dp_kind) :: d_epot, dd_epot

        ! -----------------------
        ! Sampling counters
        ! -----------------------
        integer(kind=int_kind) :: num_samples

        ! -----------------------
        ! MEANS outputs (with std) returned by md_means
        ! -----------------------
        real(kind=dp_kind) :: epot_mean,  ekin_mean,  etot_mean,  temp_mean,  press_mean
        real(kind=dp_kind) :: epot_std,   ekin_std,   etot_std,   temp_std,   press_std

        ! Means needed by thermodynamic coefficient formulas returned by md_means
        real(kind=dp_kind) :: ekin_inv_mean
        real(kind=dp_kind) :: dpot_mean, ddpot_mean
        real(kind=dp_kind) :: dpot_ekin_inv_mean, ddpot_ekin_inv_mean

        ! Optional std of extra coefficient inputs (diagnostics)
        real(kind=dp_kind) :: ekin_inv_std
        real(kind=dp_kind) :: dpot_std, ddpot_std
        real(kind=dp_kind) :: dpot_ekin_inv_std, ddpot_ekin_inv_std

        ! -----------------------
        ! Thermodynamic coefficients
        ! -----------------------
        real(kind=dp_kind) :: npd
        real(kind=dp_kind) :: degrees_of_freedom
        real(kind=dp_kind) :: aux1, aux2
        real(kind=dp_kind) :: temperature
        real(kind=dp_kind) :: pressure
        real(kind=dp_kind) :: Ca_v, Ce_v
        real(kind=dp_kind) :: gamma
        real(kind=dp_kind) :: k_s, k_s1, k_s2
        real(kind=dp_kind) :: alpha_E1, alpha_E2

        ! -----------------------
        ! Correlation settings
        ! -----------------------
        integer(kind=int_kind) :: n_samples_md, corr_max_lag
        integer(kind=int_kind) :: lag_max_out
        logical                :: corr_enabled

        real(kind=dp_kind), allocatable :: corr_vals_epot(:),  corrn_vals_epot(:)
        real(kind=dp_kind), allocatable :: corr_vals_ekin(:),  corrn_vals_ekin(:)
        real(kind=dp_kind), allocatable :: corr_vals_etot(:),  corrn_vals_etot(:)
        real(kind=dp_kind), allocatable :: corr_vals_temp(:),  corrn_vals_temp(:)
        real(kind=dp_kind), allocatable :: corr_vals_press(:), corrn_vals_press(:)

        ! -----------------------
        ! Correlation means settings (block averaging)
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
        ! -----------------------
        integer :: iu_out, iu_sum, iu_rva, ios
        integer :: iu_means

        ! ===================================================================================
        ! TRAJECTORY POST-PROCESSING (RDF / VACF / MSD)
        ! ===================================================================================
        real(kind=dp_kind), allocatable :: rux(:), ruy(:), ruz(:)
        real(kind=dp_kind), allocatable :: rx_prev(:), ry_prev(:), rz_prev(:)
        real(kind=dp_kind) :: dx, dy, dz
        integer(kind=int_kind) :: n_snapshots_expected

        ! Track which resources were successfully created (cleanup-safe)
        logical :: out_open, rva_open

        ! -----------------------
        ! Init
        ! -----------------------
        out_open     = .false.
        rva_open     = .false.
        ierr         = 0_int_kind
        corr_enabled = .false.

        ! ---------------------------------------------------------------------------
        ! Check out_dir exists (creation should be done by simulation_results)
        ! ---------------------------------------------------------------------------
        if (.not. dir_exists(out_dir)) then
            ierr = 60_int_kind
            return
        end if

        ! ===================================================================================
        ! 1) Read simulation parameters from file
        ! ===================================================================================
        call read_simulation_parameters( &
            'inputs/input_simulation_parameters.txt', &
            params, total_steps, output_interval, warmup_steps, &
            rc_over_L, target_total_energy )

        ! ===================================================================================
        ! 2) Allocate and initialize the simulation state
        ! ===================================================================================
        call init_state(params, state)

        ! ===================================================================================
        ! 3) Read initial configuration: outputs/rv_init.dat
        ! ===================================================================================
        call read_rv_init('outputs/rv_init.dat', params, state, ierr)
        if (ierr /= 0_int_kind) then
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

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

        ekin = 0.5d0 * sum( &
            state%vx(1:params%n)*state%vx(1:params%n) + &
            state%vy(1:params%n)*state%vy(1:params%n) + &
            state%vz(1:params%n)*state%vz(1:params%n) )

        etot = epot + ekin
        time = 0.d0

        ! ===================================================================================
        ! 6) Open rva.dat once (unformatted trajectory snapshots): out_dir/rva.dat
        ! ===================================================================================
        iu_rva = 31
        open(iu_rva, file=join_path(out_dir, 'rva.dat'), form='unformatted', status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if
        rva_open = .true.

        n_snapshots_expected = (total_steps / output_interval) - (warmup_steps / output_interval)
        if (n_snapshots_expected < 0_int_kind) n_snapshots_expected = 0_int_kind

        write(iu_rva) params%n, params%box_length, params%dt, output_interval, n_snapshots_expected

        ! ===================================================================================
        ! 7) STATS INIT: MEANS & CORRELATIONS
        ! ===================================================================================
        call md_means_init(means_acc, params, use_virial=.true.)
        num_samples   = 0_int_kind
        corr_enabled  = .false.

        n_samples_md = (total_steps / output_interval) - (warmup_steps / output_interval)
        if (n_samples_md < 0_int_kind) n_samples_md = 0_int_kind

        if (n_samples_md >= 2_int_kind) then
            corr_max_lag = min(1000_int_kind, n_samples_md - 1_int_kind)
            corr_max_lag = min(corr_max_lag, n_samples_md / 2_int_kind)

            call md_corr_init(corr_acc, params, n_samples_md, corr_max_lag, want_temp=.true., want_press=.true.)
            corr_enabled = .true.
        else
            corr_max_lag = 0_int_kind
        end if

        ! Energies time series (text)
        iu_out = 30
        open(iu_out, file=join_path(out_dir, 'instantaneous_energies.dat'), status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 2_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if
        out_open = .true.

        write(iu_out, '(a)') '# time epot ekin etot T P'

        ! ===================================================================================
        ! 8) Precompute constants for coefficient algebra
        ! ===================================================================================
        npd = dble(params%n)
        degrees_of_freedom = 3.d0*npd - 3.d0

        if (abs(degrees_of_freedom) < 1.d-14) then
            ierr = 20_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        aux1 = 1.d0 - 2.d0 / degrees_of_freedom
        aux2 = degrees_of_freedom/2.d0 - 1.d0

        ! ===================================================================================
        ! 9) MD loop
        ! ===================================================================================
        do step = 1, total_steps

            rx_prev(1:params%n) = state%rx(1:params%n)
            ry_prev(1:params%n) = state%ry(1:params%n)
            rz_prev(1:params%n) = state%rz(1:params%n)

            call verlet_step(params, state, epot, ekin, d_epot, dd_epot)

            do i_part = 1, params%n
                dx = state%rx(i_part) - rx_prev(i_part)
                dy = state%ry(i_part) - ry_prev(i_part)
                dz = state%rz(i_part) - rz_prev(i_part)

                dx = dx - params%box_length * dnint(dx * params%inv_box_length)
                dy = dy - params%box_length * dnint(dy * params%inv_box_length)
                dz = dz - params%box_length * dnint(dz * params%inv_box_length)

                rux(i_part) = rux(i_part) + dx
                ruy(i_part) = ruy(i_part) + dy
                ruz(i_part) = ruz(i_part) + dz
            end do

            etot = epot + ekin
            time = time + params%dt

            if ( (step > warmup_steps) .and. (mod(step, output_interval) == 0) ) then

                num_samples = num_samples + 1_int_kind
                virial_inst = -d_epot

                call md_means_add_sample( &
                    means_acc, epot, ekin, &
                    temp_inst=temp_inst, press_inst=press_inst, &
                    virial=virial_inst, d_epot=d_epot, dd_epot=dd_epot )

                write(iu_out,'(1pe13.6,5(2x,1pe13.6))') time, epot, ekin, etot, temp_inst, press_inst

                if (corr_enabled) then
                    call md_corr_add_sample(corr_acc, epot, ekin, temp=temp_inst, press=press_inst)
                end if

                write(iu_rva) state%rx(1:params%n), state%ry(1:params%n), state%rz(1:params%n)
                write(iu_rva) rux(1:params%n),      ruy(1:params%n),      ruz(1:params%n)
                write(iu_rva) state%vx(1:params%n), state%vy(1:params%n), state%vz(1:params%n)
                write(iu_rva) state%ax(1:params%n), state%ay(1:params%n), state%az(1:params%n)

            end if

        end do

        close(iu_out); out_open = .false.
        close(iu_rva); rva_open = .false.

        ! ===================================================================================
        ! 10) STATS GET (MEANS)
        ! ===================================================================================
        if (num_samples <= 0_int_kind) then
            ierr = 3_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        call md_means_get( &
            means_acc, &
            epot_mean, ekin_mean, etot_mean, temp_mean, press_mean, &
            epot_std,  ekin_std,  etot_std,  temp_std,  press_std, &
            ekin_inv_mean, ekin_inv_std, &
            dpot_mean, dpot_std, ddpot_mean, ddpot_std, &
            dpot_ekin_inv_mean, dpot_ekin_inv_std, &
            ddpot_ekin_inv_mean, ddpot_ekin_inv_std )

        ! ===================================================================================
        ! 11) Thermodynamic coefficients
        ! ===================================================================================
        temperature = 2.d0 * ekin_mean / (3.d0 * npd)
        pressure    = press_mean

        if (abs(1.d0 - aux1 * ekin_mean * ekin_inv_mean) < 1.d-14) then
            ierr = 30_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        Ca_v = 1.d0 / (1.d0 - aux1 * ekin_mean * ekin_inv_mean)
        Ce_v = Ca_v / npd

        gamma = 1.d0/Ce_v + params%volume * aux2 * (dpot_mean * ekin_inv_mean - dpot_ekin_inv_mean)

        k_s1 = (npd * temperature * (1.d0 + 2.d0 * gamma - 1.d0 / Ce_v) / params%volume) + (params%volume * ddpot_mean)
        k_s2 = k_s1 - (params%volume * aux2 * ( &
               ddpot_ekin_inv_mean - 2.d0 * dpot_mean * dpot_ekin_inv_mean + &
               dpot_mean * dpot_mean * ekin_inv_mean))

        if (abs(k_s2) < 1.d-14) then
            ierr = 31_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        k_s = 1.d0 / k_s2

        if (abs(pressure * params%volume / Ca_v - gamma * temperature) < 1.d-14) then
            ierr = 32_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        alpha_E1 = 1.d0 / (pressure * params%volume / Ca_v - gamma * temperature)

        if (abs(params%volume * (aux1 * ekin_mean * dpot_ekin_inv_mean - dpot_mean)) < 1.d-14) then
            ierr = 33_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        alpha_E2 = 1.d0 / (params%volume * (aux1 * ekin_mean * dpot_ekin_inv_mean - dpot_mean))

        ! ===================================================================================
        ! 12) CORRELATIONS (full-series autocorrelations)
        ! ===================================================================================
        if (corr_enabled) then

            call md_corr_compute(corr_acc, centered=.true.)

            allocate(corr_vals_epot(0:corr_max_lag),  corrn_vals_epot(0:corr_max_lag))
            allocate(corr_vals_ekin(0:corr_max_lag),  corrn_vals_ekin(0:corr_max_lag))
            allocate(corr_vals_etot(0:corr_max_lag),  corrn_vals_etot(0:corr_max_lag))
            allocate(corr_vals_temp(0:corr_max_lag),  corrn_vals_temp(0:corr_max_lag))
            allocate(corr_vals_press(0:corr_max_lag), corrn_vals_press(0:corr_max_lag))

            call md_corr_get_all( &
                corr_acc, &
                corr_vals_epot,  corrn_vals_epot, &
                corr_vals_ekin,  corrn_vals_ekin, &
                corr_vals_etot,  corrn_vals_etot, &
                corr_vals_temp,  corrn_vals_temp, &
                corr_vals_press, corrn_vals_press )

            lag_max_out = corr_max_lag

            call write_corr_file(join_path(out_dir, 'corr_epot.dat'),  lag_max_out, corr_vals_epot,  corrn_vals_epot,  ierr)
            if (ierr /= 0_int_kind) return

            call write_corr_file(join_path(out_dir, 'corr_ekin.dat'),  lag_max_out, corr_vals_ekin,  corrn_vals_ekin,  ierr)
            if (ierr /= 0_int_kind) return

            call write_corr_file(join_path(out_dir, 'corr_etot.dat'),  lag_max_out, corr_vals_etot,  corrn_vals_etot,  ierr)
            if (ierr /= 0_int_kind) return

            call write_corr_file(join_path(out_dir, 'corr_temp.dat'),  lag_max_out, corr_vals_temp,  corrn_vals_temp,  ierr)
            if (ierr /= 0_int_kind) return

            call write_corr_file(join_path(out_dir, 'corr_press.dat'), lag_max_out, corr_vals_press, corrn_vals_press, ierr)
            if (ierr /= 0_int_kind) return

            deallocate(corr_vals_epot,  corrn_vals_epot)
            deallocate(corr_vals_ekin,  corrn_vals_ekin)
            deallocate(corr_vals_etot,  corrn_vals_etot)
            deallocate(corr_vals_temp,  corrn_vals_temp)
            deallocate(corr_vals_press, corrn_vals_press)

            ! ===================================================================================
            ! 12b) CORR_MEANS (block-averaged correlation curves)
            ! ===================================================================================
            n_blocks = num_samples / (corr_max_lag + 1_int_kind)
            n_blocks = min(5_int_kind, n_blocks)

            if (n_blocks >= 1_int_kind) then

                call md_corr_cm_init(corr_acc, n_blocks)
                call md_corr_cm_compute(corr_acc, centered=.true.)

                allocate(corr_mean_epot(0:corr_max_lag),  corrn_mean_epot(0:corr_max_lag))
                allocate(corr_mean_ekin(0:corr_max_lag),  corrn_mean_ekin(0:corr_max_lag))
                allocate(corr_mean_etot(0:corr_max_lag),  corrn_mean_etot(0:corr_max_lag))
                allocate(corr_mean_temp(0:corr_max_lag),  corrn_mean_temp(0:corr_max_lag))
                allocate(corr_mean_press(0:corr_max_lag), corrn_mean_press(0:corr_max_lag))

                call md_corr_cm_get_all( &
                    corr_acc, &
                    corr_mean_epot,  corrn_mean_epot, &
                    corr_mean_ekin,  corrn_mean_ekin, &
                    corr_mean_etot,  corrn_mean_etot, &
                    corr_mean_temp,  corrn_mean_temp, &
                    corr_mean_press, corrn_mean_press )

                call write_corrmean_file(join_path(out_dir, 'corrmean_epot.dat'),  corr_max_lag, corr_mean_epot,  corrn_mean_epot,  ierr)
                if (ierr /= 0_int_kind) return

                call write_corrmean_file(join_path(out_dir, 'corrmean_ekin.dat'),  corr_max_lag, corr_mean_ekin,  corrn_mean_ekin,  ierr)
                if (ierr /= 0_int_kind) return

                call write_corrmean_file(join_path(out_dir, 'corrmean_etot.dat'),  corr_max_lag, corr_mean_etot,  corrn_mean_etot,  ierr)
                if (ierr /= 0_int_kind) return

                call write_corrmean_file(join_path(out_dir, 'corrmean_temp.dat'),  corr_max_lag, corr_mean_temp,  corrn_mean_temp,  ierr)
                if (ierr /= 0_int_kind) return

                call write_corrmean_file(join_path(out_dir, 'corrmean_press.dat'), corr_max_lag, corr_mean_press, corrn_mean_press, ierr)
                if (ierr /= 0_int_kind) return

                deallocate(corr_mean_epot,  corrn_mean_epot)
                deallocate(corr_mean_ekin,  corrn_mean_ekin)
                deallocate(corr_mean_etot,  corrn_mean_etot)
                deallocate(corr_mean_temp,  corrn_mean_temp)
                deallocate(corr_mean_press, corrn_mean_press)

                call md_corr_cm_free(corr_acc)

            end if

            call md_corr_free(corr_acc)

        end if

        ! ===================================================================================
        ! 13) Write final summary (per-run): out_dir/md_final_results.txt
        ! ===================================================================================
        iu_sum = 10
        open(iu_sum, file=join_path(out_dir, 'md_final_results.txt'), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 4_int_kind
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        write(iu_sum,'(a)')            '************** MD PRODUCTION RESULTS **************'
        write(iu_sum,'(a,1x,i8)')      'num_particles:', params%n
        write(iu_sum,'(a,1x,i8)')      'num_cells:',     params%num_cells
        write(iu_sum,'(a,1x,1pe19.12)')'box_length:',    params%box_length
        write(iu_sum,'(a,1x,1pe19.12)')'volume:',        params%volume
        write(iu_sum,'(a,1x,1pe19.12)')'density:',       dble(params%n)/params%volume
        write(iu_sum,'(a,1x,1pe19.12)')'time_step:',     params%dt
        write(iu_sum,'(a,1x,i8)')      'output_interval:', output_interval
        write(iu_sum,'(a,1x,i10)')     'total_steps:',     total_steps
        write(iu_sum,'(a,1x,i10)')     'warmup_steps:',    warmup_steps

        write(iu_sum,'(a)')            '-------------------- Averages --------------------'
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Epot>:', epot_mean,  'std:', epot_std
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Ekin>:', ekin_mean,  'std:', ekin_std
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Etot>:', etot_mean,  'std:', etot_std
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<T>  :',  temp_mean,  'std:', temp_std
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<P>  :',  press_mean, 'std:', press_std

        write(iu_sum,'(a)')            '-------------- Thermodynamic coefficients --------------'
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Temperature:', temperature, 'Pressure:', pressure
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Ca_v:',       Ca_v,        'Ce_v:',    Ce_v
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'K_s:',        k_s,         'Gamma:',   gamma
        write(iu_sum,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Alpha_E1:',   alpha_E1,    'Alpha_E2:', alpha_E2
        write(iu_sum,'(a)')            '--------------------------------------------------------'

        close(iu_sum)

        ! ===================================================================================
        ! 14) Write means file (machine-friendly, per-run): out_dir/means.dat
        ! ===================================================================================
        call write_means_file( &
            join_path(out_dir, 'means.dat'), params, total_steps, warmup_steps, output_interval, rc_over_L, target_total_energy, &
            n_samples_md, num_samples, corr_max_lag, &
            epot_mean, epot_std, ekin_mean, ekin_std, etot_mean, etot_std, temp_mean, temp_std, press_mean, press_std, &
            ekin_inv_mean, ekin_inv_std, dpot_mean, dpot_std, ddpot_mean, ddpot_std, &
            dpot_ekin_inv_mean, dpot_ekin_inv_std, ddpot_ekin_inv_mean, ddpot_ekin_inv_std, &
            temperature, pressure, Ca_v, Ce_v, gamma, k_s, alpha_E1, alpha_E2, ierr )

        if (ierr /= 0_int_kind) then
            call cleanup_on_error(out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
                                  rux, ruy, ruz, rx_prev, ry_prev, rz_prev)
            return
        end if

        deallocate(rux, ruy, ruz)
        deallocate(rx_prev, ry_prev, rz_prev)

    end subroutine run_md_simulation

    ! ------------------------------------------------------------
    ! Read initial configuration snapshot (binary, unformatted)
    ! - File outputs/rv_init.dat contains:
    !   Record 1: rx, ry, rz (particle positions)
    !   Record 2: vx, vy, vz (particle velocities)
    ! ------------------------------------------------------------
    subroutine read_rv_init(filename, p, s, ierr)

        implicit none

        character(len=*), intent(in)    :: filename
        type(sim_params), intent(in)    :: p
        type(sim_state),  intent(inout) :: s
        integer(kind=int_kind), intent(out) :: ierr

        integer :: iu, ios

        ierr = 0_int_kind

        iu = 21
        open(iu, file=filename, form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 10_int_kind
            return
        end if

        read(iu, iostat=ios) s%rx(1:p%n), s%ry(1:p%n), s%rz(1:p%n)
        if (ios /= 0) then
            ierr = 11_int_kind
            close(iu)
            return
        end if

        read(iu, iostat=ios) s%vx(1:p%n), s%vy(1:p%n), s%vz(1:p%n)
        if (ios /= 0) then
            ierr = 12_int_kind
            close(iu)
            return
        end if

        close(iu)

    end subroutine read_rv_init

    ! ------------------------------------------------------------
    ! INTERNAL: write one correlation file (raw + normalized)
    ! ------------------------------------------------------------
    subroutine write_corr_file(filename, lag_max, corr_in, corrn_in, ierr)

        implicit none

        character(len=*), intent(in) :: filename
        integer(kind=int_kind), intent(in) :: lag_max
        real(kind=dp_kind), intent(in) :: corr_in(0:), corrn_in(0:)
        integer(kind=int_kind), intent(out) :: ierr

        integer :: iu, ios
        integer(kind=int_kind) :: lag

        ierr = 0_int_kind

        iu = 50
        open(iu, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 40_int_kind
            return
        end if

        write(iu,'(a)') '# lag C(lag) C_norm(lag)'
        do lag = 0, lag_max
            write(iu,'(i8,2(2x,1pe19.12))') lag, corr_in(lag), corrn_in(lag)
        end do

        close(iu)

    end subroutine write_corr_file

    ! ------------------------------------------------------------
    ! INTERNAL: write one block-averaged correlation file
    ! ------------------------------------------------------------
    subroutine write_corrmean_file(filename, lag_max, corr_in, corrn_in, ierr)

        implicit none

        character(len=*), intent(in) :: filename
        integer(kind=int_kind), intent(in) :: lag_max
        real(kind=dp_kind), intent(in) :: corr_in(0:), corrn_in(0:)
        integer(kind=int_kind), intent(out) :: ierr

        integer :: iu, ios
        integer(kind=int_kind) :: lag

        ierr = 0_int_kind

        iu = 51
        open(iu, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 41_int_kind
            return
        end if

        write(iu,'(a)') '# lag <C(lag)>_blocks <C_norm(lag)>_blocks'
        do lag = 0, lag_max
            write(iu,'(i8,2(2x,1pe19.12))') lag, corr_in(lag), corrn_in(lag)
        end do

        close(iu)

    end subroutine write_corrmean_file

    ! ------------------------------------------------------------
    ! INTERNAL: write means + std + run parameters (per-run file)
    ! ------------------------------------------------------------
    subroutine write_means_file( &
        filename, params, total_steps, warmup_steps, output_interval, rc_over_L, target_total_energy, &
        n_samples_md, num_samples, corr_max_lag, &
        epot_mean, epot_std, ekin_mean, ekin_std, etot_mean, etot_std, temp_mean, temp_std, press_mean, press_std, &
        ekin_inv_mean, ekin_inv_std, dpot_mean, dpot_std, ddpot_mean, ddpot_std, &
        dpot_ekin_inv_mean, dpot_ekin_inv_std, ddpot_ekin_inv_mean, ddpot_ekin_inv_std, &
        temperature, pressure, Ca_v, Ce_v, gamma, k_s, alpha_E1, alpha_E2, ierr )

        implicit none

        character(len=*), intent(in) :: filename
        type(sim_params), intent(in) :: params

        integer(kind=int_kind), intent(in) :: total_steps, warmup_steps, output_interval
        real(kind=dp_kind),     intent(in) :: rc_over_L, target_total_energy
        integer(kind=int_kind), intent(in) :: n_samples_md, num_samples, corr_max_lag

        real(kind=dp_kind), intent(in) :: epot_mean, epot_std, ekin_mean, ekin_std, etot_mean, etot_std
        real(kind=dp_kind), intent(in) :: temp_mean, temp_std, press_mean, press_std

        real(kind=dp_kind), intent(in) :: ekin_inv_mean, ekin_inv_std
        real(kind=dp_kind), intent(in) :: dpot_mean, dpot_std, ddpot_mean, ddpot_std
        real(kind=dp_kind), intent(in) :: dpot_ekin_inv_mean, dpot_ekin_inv_std
        real(kind=dp_kind), intent(in) :: ddpot_ekin_inv_mean, ddpot_ekin_inv_std

        real(kind=dp_kind), intent(in) :: temperature, pressure, Ca_v, Ce_v, gamma, k_s, alpha_E1, alpha_E2
        integer(kind=int_kind), intent(out) :: ierr

        integer :: iu, ios

        ierr = 0_int_kind

        iu = 70
        open(iu, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 70_int_kind
            return
        end if

        write(iu,'(a)') '# -------------------- Run parameters --------------------'
        write(iu,'(a,1x,i12)')      'num_particles:', params%n
        write(iu,'(a,1x,i12)')      'num_cells:',     params%num_cells
        write(iu,'(a,1x,1pe19.12)') 'box_length:',    params%box_length
        write(iu,'(a,1x,1pe19.12)') 'volume:',        params%volume
        write(iu,'(a,1x,1pe19.12)') 'density:',       dble(params%n)/params%volume
        write(iu,'(a,1x,1pe19.12)') 'time_step:',     params%dt
        write(iu,'(a,1x,1pe19.12)') 'rc_over_L:',     rc_over_L
        write(iu,'(a,1x,1pe19.12)') 'target_total_energy:', target_total_energy
        write(iu,'(a,1x,i12)')      'output_interval:', output_interval
        write(iu,'(a,1x,i12)')      'total_steps:',    total_steps
        write(iu,'(a,1x,i12)')      'warmup_steps:',   warmup_steps
        write(iu,'(a,1x,i12)')      'n_samples_md_expected:', n_samples_md
        write(iu,'(a,1x,i12)')      'num_samples_taken:',     num_samples
        write(iu,'(a,1x,i12)')      'corr_max_lag:',          corr_max_lag

        write(iu,'(a)') '# -------------------- Means and std --------------------'
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Epot>:', epot_mean,  'std:', epot_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Ekin>:', ekin_mean,  'std:', ekin_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<Etot>:', etot_mean,  'std:', etot_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<T>  :',  temp_mean,  'std:', temp_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<P>  :',  press_mean, 'std:', press_std

        write(iu,'(a)') '# -------------------- Extra means and std (diagnostics) --------------------'
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<1/Ekin>:',            ekin_inv_mean,        'std:', ekin_inv_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<d_epot>:',            dpot_mean,            'std:', dpot_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<dd_epot>:',           ddpot_mean,           'std:', ddpot_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<d_epot*(1/Ekin)>:',   dpot_ekin_inv_mean,   'std:', dpot_ekin_inv_std
        write(iu,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') '<dd_epot*(1/Ekin)>:',  ddpot_ekin_inv_mean,  'std:', ddpot_ekin_inv_std

        write(iu,'(a)') '# -------------------- Thermodynamic coefficients --------------------'
        write(iu,'(a,1x,1pe19.12)') 'Temperature:', temperature
        write(iu,'(a,1x,1pe19.12)') 'Pressure:',    pressure
        write(iu,'(a,1x,1pe19.12)') 'Ca_v:',        Ca_v
        write(iu,'(a,1x,1pe19.12)') 'Ce_v:',        Ce_v
        write(iu,'(a,1x,1pe19.12)') 'Gamma:',       gamma
        write(iu,'(a,1x,1pe19.12)') 'K_s:',         k_s
        write(iu,'(a,1x,1pe19.12)') 'Alpha_E1:',    alpha_E1
        write(iu,'(a,1x,1pe19.12)') 'Alpha_E2:',    alpha_E2

        close(iu)

    end subroutine write_means_file

    ! ------------------------------------------------------------
    ! INTERNAL: cleanup helper for early returns
    ! ------------------------------------------------------------
    subroutine cleanup_on_error( &
        out_open, rva_open, iu_out, iu_rva, corr_enabled, corr_acc, means_acc, &
        rux, ruy, ruz, rx_prev, ry_prev, rz_prev )

        implicit none

        logical, intent(inout) :: out_open, rva_open
        integer, intent(in)    :: iu_out, iu_rva

        logical, intent(in) :: corr_enabled
        type(md_corr_accum),  intent(inout) :: corr_acc
        type(md_means_accum), intent(inout) :: means_acc

        real(kind=dp_kind), allocatable, intent(inout) :: rux(:), ruy(:), ruz(:)
        real(kind=dp_kind), allocatable, intent(inout) :: rx_prev(:), ry_prev(:), rz_prev(:)

        if (out_open) then
            close(iu_out)
            out_open = .false.
        end if

        if (rva_open) then
            close(iu_rva)
            rva_open = .false.
        end if

        if (allocated(rux))     deallocate(rux)
        if (allocated(ruy))     deallocate(ruy)
        if (allocated(ruz))     deallocate(ruz)
        if (allocated(rx_prev)) deallocate(rx_prev)
        if (allocated(ry_prev)) deallocate(ry_prev)
        if (allocated(rz_prev)) deallocate(rz_prev)

        if (corr_enabled) then
            call md_corr_free(corr_acc)
        end if

    end subroutine cleanup_on_error

    ! ------------------------------------------------------------
    ! INTERNAL: join paths safely (supports 'dir/file' and 'dir\file')
    ! ------------------------------------------------------------
    pure function join_path(dir, file) result(path)

        implicit none

        character(len=*), intent(in) :: dir, file
        character(len=:), allocatable :: path

        integer :: n
        character(len=1) :: lastc

        n = len_trim(dir)
        if (n <= 0) then
            path = trim(file)
            return
        end if

        lastc = dir(n:n)
        if (lastc == '/' .or. lastc == '\') then
            path = trim(dir)//trim(file)
        else
            path = trim(dir)//'/'//trim(file)
        end if

    end function join_path

    ! ------------------------------------------------------------
    ! INTERNAL: directory exists check using INQUIRE(file=...)
    ! ------------------------------------------------------------
    logical function dir_exists(dir)

        implicit none

        character(len=*), intent(in) :: dir
        logical :: ex

        inquire(file=trim(dir), exist=ex)
        dir_exists = ex

    end function dir_exists

end module md_simulation
