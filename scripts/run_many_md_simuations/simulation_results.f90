!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! SIMULATION_RESULTS (post-processing of time-averaged data)
!
! - Analyzes data and computes stats between independent simulations.
! - Reads "./outputs/means.txt" containing one line per independent run (or block).
! - Each line stores time-averaged magnitudes produced by the MD code:
!
!     <ekin>  <1/ekin>  <d_epot>  <dd_epot>  <d_epot*(1/ekin)>  <dd_epot*(1/ekin)>
!
!   where:
!     <ekin>             : mean kinetic energy
!     <1/ekin>           : mean of the inverse kinetic energy
!     <d_epot>           : mean of the LJ derivative term used in the pressure formula
!     <dd_epot>          : mean of the LJ second-derivative term used in K_s
!     <d_epot*(1/ekin)>  : mean of d_epot / ekin
!     <dd_epot*(1/ekin)> : mean of dd_epot / ekin
!
! - For each run, computes thermodynamic coefficients using the same algebra as in
!   molecular_dynamics_simulation.f90:
!
!   Temperature: T = 2 <ekin> / (3 N)
!
!   Pressure   : P = (N T)/V - <d_epot>
!
!   Heat capacities (microcanonical estimator):
!     f    = 3N - 3
!     aux1 = 1 - 2/f
!     C_A,v = 1 / [ 1 - aux1 <ekin><1/ekin> ]
!     C_E,v = C_A,v / N
!
!   Grüneisen parameter:
!     aux2 = f/2 - 1
!     gamma = 1/C_E,v + V*aux2 [ <d_epot><1/ekin> - <d_epot/ekin> ]
!
!   Isentropic compressibility-like coefficient:
!     K_s = 1 / k_s2
!
!   Thermal expansion estimators:
!     alpha_E1 = 1 / ( P V / C_A,v - gamma T )
!     alpha_E2 = 1 / ( V [ aux1 <ekin><d_epot/ekin> - <d_epot> ] )
!
! - After computing run-by-run values, computes:
!   mean across runs
!   standard deviation across runs (no error propagation)
!
! - Writes report to "./outputs/simulation_results.txt".
!
! IMPORTANT:
! - The normalization/meaning of d_epot and dd_epot must match the MD code that produced means.txt.
! - This program assumes num_particles, box_length, volume, density, etc. are already initialized
!   (e.g., by reading the same input file and calling update_derived()).
!******************************************************************************************************

program simulation_results

    use define_precision,  only: dp_kind, int_kind
    use common_variables,  only: num_particles, box_length, volume, density

    implicit none

    ! -----------------------
    ! File paths
    ! -----------------------
    character(len=*), parameter :: input_file  = './outputs/means.txt'
    character(len=*), parameter :: output_file = './outputs/simulation_results.txt'

    ! -----------------------
    ! Run-level data read from means.txt (one value per run)
    ! -----------------------
    integer(kind=int_kind) :: n_runs, i

    real(kind=dp_kind), allocatable :: ekin_mean(:)
    real(kind=dp_kind), allocatable :: ekin_inv_mean(:)
    real(kind=dp_kind), allocatable :: d_epot_mean(:)
    real(kind=dp_kind), allocatable :: dd_epot_mean(:)
    real(kind=dp_kind), allocatable :: d_epot_ekin_inv_mean(:)
    real(kind=dp_kind), allocatable :: dd_epot_ekin_inv_mean(:)

    ! -----------------------
    ! Computed magnitudes per run
    ! -----------------------
    real(kind=dp_kind), allocatable :: temperature(:)
    real(kind=dp_kind), allocatable :: pressure(:)
    real(kind=dp_kind), allocatable :: Ca_v(:), Ce_v(:)
    real(kind=dp_kind), allocatable :: gamma(:)
    real(kind=dp_kind), allocatable :: K_s(:)
    real(kind=dp_kind), allocatable :: alpha_E1(:), alpha_E2(:)

    ! -----------------------
    ! Constants for formulas (common to all runs)
    ! -----------------------
    real(kind=dp_kind) :: npd
    real(kind=dp_kind) :: dof, aux1, aux2
    real(kind=dp_kind) :: ks1, ks2

    ! ===================================================================================
    ! 0) Minimal guard: formulas need N and V
    ! ===================================================================================
    if (num_particles <= 0 .or. volume <= 0.d0) then
        stop 'simulation_results(): num_particles/volume are not initialized. Read inputs and call update_derived() before.'
    end if

    npd = dble(num_particles)
    dof = 3.d0*npd - 3.d0

    if (dof <= 0.d0) then
        stop 'simulation_results(): invalid degrees of freedom.'
    end if

    aux1 = 1.d0 - 2.d0 / dof
    aux2 = dof/2.d0 - 1.d0

    ! ===================================================================================
    ! 1) Count independent runs (data lines)
    ! ===================================================================================
    n_runs = count_data_lines(input_file)
    if (n_runs <= 0_int_kind) then
        stop 'simulation_results(): means.txt contains no data lines.'
    end if

    ! ===================================================================================
    ! 2) Allocate arrays (one entry per run)
    ! ===================================================================================
    allocate( &
        ekin_mean(n_runs), ekin_inv_mean(n_runs), d_epot_mean(n_runs), dd_epot_mean(n_runs), &
        d_epot_ekin_inv_mean(n_runs), dd_epot_ekin_inv_mean(n_runs) )

    allocate( &
        temperature(n_runs), pressure(n_runs), Ca_v(n_runs), Ce_v(n_runs), gamma(n_runs), &
        K_s(n_runs), alpha_E1(n_runs), alpha_E2(n_runs) )

    ! ===================================================================================
    ! 3) Read means.txt into arrays
    ! Expected columns per run:
    !   <ekin> <1/ekin> <d_epot> <dd_epot> <d_epot*(1/ekin)> <dd_epot*(1/ekin)>
    ! ===================================================================================
    call read_means_file( &
        input_file, n_runs, &
        ekin_mean, ekin_inv_mean, d_epot_mean, dd_epot_mean, &
        d_epot_ekin_inv_mean, dd_epot_ekin_inv_mean )

    ! ===================================================================================
    ! 4) Compute thermodynamic coefficients run by run
    ! ===================================================================================
    do i = 1, n_runs

        if (ekin_mean(i) <= 0.d0)     stop 'simulation_results(): <ekin> must be > 0.'
        if (ekin_inv_mean(i) <= 0.d0) stop 'simulation_results(): <1/ekin> must be > 0.'

        ! Temperature: T = 2 <ekin> / (3 N)
        temperature(i) = 2.d0 * ekin_mean(i) / (3.d0 * npd)

        ! Pressure: P = (N T)/V - <d_epot>
        pressure(i) = npd * temperature(i) / volume - d_epot_mean(i)

        ! Heat capacities:
        ! C_A,v = 1 / [ 1 - aux1 <ekin><1/ekin> ]
        ! C_E,v = C_A,v / N
        if (abs(1.d0 - aux1 * ekin_mean(i) * ekin_inv_mean(i)) < 1.d-14) then
            stop 'simulation_results(): Ca_v denominator too small.'
        end if

        Ca_v(i) = 1.d0 / (1.d0 - aux1 * ekin_mean(i) * ekin_inv_mean(i))
        Ce_v(i) = Ca_v(i) / npd

        if (abs(Ce_v(i)) < 1.d-14) then
            stop 'simulation_results(): Ce_v too small (division by zero risk).'
        end if

        ! Grüneisen parameter:
        ! gamma = 1/C_E,v + V*aux2 [ <d_epot><1/ekin> - <d_epot/ekin> ]
        gamma(i) = 1.d0 / Ce_v(i) + volume * aux2 * &
                   (d_epot_mean(i) * ekin_inv_mean(i) - d_epot_ekin_inv_mean(i))

        ! Isentropic compressibility-like coefficient:
        ! K_s = 1 / ks2
        !
        ! ks1 = (N T / V) * (1 + 2 gamma - 1/C_E,v) + V <dd_epot>
        ! ks2 = ks1 - V*aux2 [ <dd_epot/ekin> - 2 <d_epot><d_epot/ekin> + <d_epot>^2 <1/ekin> ]
        ks1 = (npd * temperature(i) * (1.d0 + 2.d0*gamma(i) - 1.d0/Ce_v(i)) / volume) + &
              (volume * dd_epot_mean(i))

        ks2 = ks1 - volume * aux2 * &
              (dd_epot_ekin_inv_mean(i) - 2.d0 * d_epot_mean(i) * d_epot_ekin_inv_mean(i) + &
               d_epot_mean(i) * d_epot_mean(i) * ekin_inv_mean(i))

        if (abs(ks2) < 1.d-14) then
            stop 'simulation_results(): ks2 too small (division by zero risk).'
        end if

        K_s(i) = 1.d0 / ks2

        ! Thermal expansion coefficient estimators:
        ! alpha_E1 = 1 / ( P V / C_A,v - gamma T )
        ! alpha_E2 = 1 / ( V [ aux1 <ekin><d_epot/ekin> - <d_epot> ] )
        if (abs(pressure(i)*volume/Ca_v(i) - gamma(i)*temperature(i)) < 1.d-14) then
            stop 'simulation_results(): alpha_E1 denominator too small.'
        end if

        alpha_E1(i) = 1.d0 / (pressure(i)*volume/Ca_v(i) - gamma(i)*temperature(i))

        if (abs(volume * (aux1*ekin_mean(i)*d_epot_ekin_inv_mean(i) - d_epot_mean(i))) < 1.d-14) then
            stop 'simulation_results(): alpha_E2 denominator too small.'
        end if

        alpha_E2(i) = 1.d0 / (volume * (aux1*ekin_mean(i)*d_epot_ekin_inv_mean(i) - d_epot_mean(i)))

    end do

    ! ===================================================================================
    ! 5) Write report: per-run values + global mean/std across runs
    ! ===================================================================================
    call write_results(output_file, n_runs, temperature, pressure, Ca_v, Ce_v, gamma, K_s, alpha_E1, alpha_E2)

    write(*,'(a)') 'simulation_results: done.'

contains

    ! ------------------------------------------------------------
    ! Count non-comment, non-empty lines in the input file.
    ! - Lines starting with '#' are ignored.
    ! - Used to know how many runs are present.
    ! ------------------------------------------------------------
    integer(kind=int_kind) function count_data_lines(fname)

        implicit none

        character(len=*), intent(in) :: fname
        integer :: unit, ios
        character(len=512) :: line

        count_data_lines = 0_int_kind
        unit = 99

        open(unit, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            stop 'count_data_lines(): cannot open means.txt.'
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            count_data_lines = count_data_lines + 1_int_kind
        end do

        close(unit)

    end function count_data_lines

    ! ------------------------------------------------------------
    ! Read means.txt into arrays.
    !
    ! Expected columns per data line (run):
    !   a1 a2 a3 a4 a5 a6
    !
    ! mapped as:
    !   a1 -> <ekin>
    !   a2 -> <1/ekin>
    !   a3 -> <d_epot>
    !   a4 -> <dd_epot>
    !   a5 -> <d_epot*(1/ekin)>
    !   a6 -> <dd_epot*(1/ekin)>
    ! ------------------------------------------------------------
    subroutine read_means_file(fname, n, ek, ek_inv, dp, ddp, dp_einv, ddp_einv)

        implicit none

        character(len=*), intent(in) :: fname
        integer(kind=int_kind), intent(in) :: n
        real(kind=dp_kind), intent(out) :: ek(:), ek_inv(:), dp(:), ddp(:), dp_einv(:), ddp_einv(:)

        integer(kind=int_kind) :: idx
        integer :: unit, ios
        character(len=512) :: line
        real(kind=dp_kind) :: a1, a2, a3, a4, a5, a6

        unit = 98
        idx  = 0_int_kind

        open(unit, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            stop 'read_means_file(): cannot open means.txt.'
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            read(line, *, iostat=ios) a1, a2, a3, a4, a5, a6
            if (ios /= 0) then
                stop 'read_means_file(): bad data line format in means.txt.'
            end if

            idx = idx + 1_int_kind
            if (idx > n) exit

            ek(idx)       = a1
            ek_inv(idx)   = a2
            dp(idx)       = a3
            ddp(idx)      = a4
            dp_einv(idx)  = a5
            ddp_einv(idx) = a6
        end do

        close(unit)

        if (idx /= n) then
            stop 'read_means_file(): number of read lines does not match expected count.'
        end if

    end subroutine read_means_file

    ! ------------------------------------------------------------
    ! Compute mean and standard deviation for an array.
    !
    ! Mean:
    !   <x> = (1/n) sum_i x_i
    !
    ! Standard deviation (population form, no Bessel correction):
    !   std = sqrt( <x^2> - <x>^2 )
    ! ------------------------------------------------------------
    subroutine mean_and_std(x, mean_x, std_x)

        implicit none

        real(kind=dp_kind), intent(in)  :: x(:)
        real(kind=dp_kind), intent(out) :: mean_x, std_x

        integer(kind=int_kind) :: n
        real(kind=dp_kind) :: m2

        n = size(x)
        if (n <= 0) then
            stop 'mean_and_std(): empty array.'
        end if

        mean_x = sum(x) / dble(n)
        m2     = sum(x*x) / dble(n)

        std_x = dsqrt(max(0.d0, m2 - mean_x*mean_x))

    end subroutine mean_and_std

    ! ------------------------------------------------------------
    ! Write per-run values and final mean/std across runs.
    ! - Output includes:
    !   * constants (N, L, V, rho)
    !   * per-run coefficients
    !   * global mean and std for each coefficient
    ! ------------------------------------------------------------
    subroutine write_results(fname, n, T, P, Cav, Cev, g, Ks, a1, a2)

        implicit none

        character(len=*), intent(in) :: fname
        integer(kind=int_kind), intent(in) :: n
        real(kind=dp_kind), intent(in) :: T(:), P(:), Cav(:), Cev(:), g(:), Ks(:), a1(:), a2(:)

        integer :: unit
        integer(kind=int_kind) :: i
        real(kind=dp_kind) :: m, s

        unit = 77
        open(unit, file=fname, status='replace', action='write')

        write(unit,'(a)') '************** SIMULATION RESULTS **************'
        write(unit,'(a,1x,i8)')       'num_particles:', num_particles
        write(unit,'(a,1x,1pe19.12)') 'box_length:',    box_length
        write(unit,'(a,1x,1pe19.12)') 'volume:',        volume
        write(unit,'(a,1x,1pe19.12)') 'density:',       density
        write(unit,'(a)') '------------------------------------------------'

        write(unit,'(a)') 'Per-run values (reduced units):'
        write(unit,'(a)') '# run T P Ca_v Ce_v gamma K_s alpha_E1 alpha_E2'

        do i = 1, n
            write(unit,'(i5,1x,8(1x,1pe11.4))') i, T(i), P(i), Cav(i), Cev(i), g(i), Ks(i), a1(i), a2(i)
        end do

        write(unit,'(a)') '------------------------------------------------'
        write(unit,'(a)') 'Final mean and statistical uncertainty (std across runs):'

        call mean_and_std(T,   m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'T:',        m, 'std:', s
        call mean_and_std(P,   m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'P:',        m, 'std:', s
        call mean_and_std(Cav, m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Ca_v:',     m, 'std:', s
        call mean_and_std(Cev, m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'Ce_v:',     m, 'std:', s
        call mean_and_std(g,   m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'gamma:',    m, 'std:', s
        call mean_and_std(Ks,  m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'K_s:',      m, 'std:', s
        call mean_and_std(a1,  m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'alpha_E1:', m, 'std:', s
        call mean_and_std(a2,  m, s); write(unit,'(a,1x,1pe19.12,2x,a,1x,1pe19.12)') 'alpha_E2:', m, 'std:', s

        write(unit,'(a)') '************************************************'
        close(unit)

    end subroutine write_results

end program simulation_results
