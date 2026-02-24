!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! MD_TYPES MODULE
! - Centralizes the main "containers" used by the Molecular Dynamics (MD) code:
!     (1) sim_params   : run configuration (box, dt, rc, sizes, switches)
!     (2) sim_state    : dynamic state (positions, velocities, accelerations)
!     (3) inst_obs     : instantaneous observables (epot, ekin, virial, ...)
!     (4) accum_means  : running sums for time-averages (optional for now)
!******************************************************************************************************

module md_types
  use define_precision
  implicit none
  private

  public :: sim_params, sim_state, inst_obs, accum_means
  public :: init_params, compute_derived_params
  public :: allocate_state, deallocate_state, zero_state, init_state

  real(kind=dp_kind), parameter, public :: pi = 3.1415926535897932384626433832795d0

  ! ============================================================
  ! 1) Parameters (mostly constant during a run)
  ! ============================================================
  type :: sim_params
    ! Sizes
    integer(kind=int_kind) :: n          = 0          ! number of particles (N)
    integer(kind=int_kind) :: num_cells  = 0          ! FCC cells per edge

    ! Box
    real(kind=dp_kind) :: box_length     = 0.d0       ! L
    real(kind=dp_kind) :: inv_box_length = 0.d0       ! 1/L
    real(kind=dp_kind) :: volume         = 0.d0       ! V = L^3
    real(kind=dp_kind) :: density        = 0.d0       ! rho = N/V

    ! Integration => needed for velocity verlet calculations in verlet.f90
    real(kind=dp_kind) :: dt               = 0.d0     ! time step
    real(kind=dp_kind) :: dt_half          = 0.d0     ! 0.5*dt
    real(kind=dp_kind) :: dt_square_half   = 0.d0     ! 0.5*dt^2

    ! Lennard Jones Potential (LJ) cutoff
    real(kind=dp_kind) :: rc  = 0.d0                  ! Cutoff Radius
    real(kind=dp_kind) :: rc_square = 0.d0            ! Cutoff Radius squared

    ! Legacy scaling
    real(kind=dp_kind) :: length_scale     = 0.d0     ! pl
    real(kind=dp_kind) :: inv_length_scale = 0.d0     ! 1/pl
  end type sim_params


  ! ============================================================
  ! 2) Dynamic state (changes each time step)
  ! ============================================================
  type :: sim_state
    real(kind=dp_kind), allocatable :: rx(:), ry(:), rz(:)   ! Positions    
    real(kind=dp_kind), allocatable :: vx(:), vy(:), vz(:)   ! Velocities    
    real(kind=dp_kind), allocatable :: ax(:), ay(:), az(:)   ! Accelerations (forces / mass)
  end type sim_state


  ! ============================================================
  ! 3) Instantaneous observables per time step (frame)
  !    - Instantaneous observables correspond to a single MD configuration (frame).
  !    - They are recomputed and overwritten at every time step and are not accumulated over time.
  ! ============================================================
  type :: inst_obs
    real(kind=dp_kind) :: epot   = 0.d0   ! Potential Energy
    real(kind=dp_kind) :: ekin   = 0.d0   ! Kinetic Energy
    real(kind=dp_kind) :: etot   = 0.d0   ! Total Energy = epot + ekin
    real(kind=dp_kind) :: virial = 0.d0   ! Virial (W): sum over interacting particle pairs (sum _{i<j} r_ij*F_ij). It accounts for the contribution of interparticle forces to the pressure.
    real(kind=dp_kind) :: temp   = 0.d0   ! Temperature
    real(kind=dp_kind) :: press  = 0.d0   ! Pressure
  end type inst_obs


  ! ============================================================
  ! 4) Accumulators for means -> initialize
  !    - Variables are used to compute statistical averages over many time steps of a single MD run.
  !    - The averages are obtained by accumulating the sum of each quantity, and the sum of its square,
  !      which allows estimation of mean values and (variance, standard deviation) without storing all samples.
  !
  ! The accumulators are reset at the beginning of each run and updated every sampling step.
  ! ============================================================
  type :: accum_means
    integer(kind=int_kind) :: nsamples = 0   ! Number of sampled configurations used in the averages

    real(kind=dp_kind) :: sum_epot  = 0.d0   ! Sum of potential energy samples
    real(kind=dp_kind) :: sum_epot2 = 0.d0   ! Sum of squared potential energy samples

    real(kind=dp_kind) :: sum_ekin  = 0.d0   ! Sum of kinetic energy samples
    real(kind=dp_kind) :: sum_ekin2 = 0.d0   ! Sum of squared kinetic energy samples

    real(kind=dp_kind) :: sum_vir   = 0.d0   ! Sum of virial samples
    real(kind=dp_kind) :: sum_vir2  = 0.d0   ! Sum of squared virial samples
  end type accum_means
 

contains

  ! ------------------------------------------------------------
  ! Initialize params object: Given each variable its correspondant value defined in inputs/simulation_params.txt
  ! ------------------------------------------------------------
  subroutine init_params(p, n, box_length, dt, rc, num_cells)
    type(sim_params), intent(inout) :: p   ! NOTE: p%box_length access box_length saved in sim_params which has been given the name p
    
    integer(kind=int_kind), intent(in) :: n                     ! Number of particles
    real(kind=dp_kind), intent(in) :: box_length, dt, rc        ! L, rho, time step, cutoff radius
    integer(kind=int_kind), intent(in), optional :: num_cells   ! Number of cells: number of FCC unit cells per spatial direction (x, y, z).
                                                                ! A unit cell is the smallest repeating building block that generates the entire lattice.

    p%n          = n
    p%box_length = box_length
    p%dt         = dt 
    p%rc         = rc
    if (present(num_cells)) p%num_cells = num_cells

    call compute_derived_params(p)
  end subroutine init_params
  
  ! ------------------------------------------------------------
  ! Compute/update quantities derived from the input parameters that:
  !  - are used repeatedly during the simulation
  !  - depend only on the basic input parameters
  !  - should NOT be recomputed inside time–integration loops
  !
  ! Examples:
  !  - Inverse box length and volume are derived from box_length L => inverse_box_length = 1/L, Volume = L^3
  !  - rc_squared is derived from rc => rc_squared = rc^2
  ! ------------------------------------------------------------
  subroutine compute_derived_params(p)
    type(sim_params), intent(inout) :: p

    ! Box-related derived quantities
    if (p%box_length > 0.d0) then
      p%inv_box_length = 1.d0 / p%box_length     ! 1 / L
      p%volume         = p%box_length**3         ! Simulation box volume
      
      ! Density is derived from N and V (policy: L is authoritative)
      if (p%n > 0) p%density = p%n / p%volume
    else
      stop 'compute_derived_params(): box_length must be > 0.'
    end if

    ! Lennard-Jones cutoff
    if (p%rc > 0.d0) then
      p%rc_square = p%rc * p%rc   ! rc^2
    else
      stop 'compute_derived_params(): rc (cutoff_radius) must be > 0.'
    end if
    if (p%rc >= 0.5d0 * p%box_length) then
      stop 'compute_derived_params(): rc (cutoff_radius) must be < L/2 (minimum image convention).'
    end if

    ! Time–step combinations: Used in Velocity Verlet integration
    if (p%dt > 0.d0) then
      p%dt_half        = 0.5d0 * p%dt            ! dt / 2
      p%dt_square_half = p%dt_half * p%dt        ! dt^2 / 2
    else
      stop 'compute_derived_params(): dt must be > 0.'
    end if

    ! Legacy scaling
    if (p%length_scale > 0.d0) then
      p%inv_length_scale = 1.d0 / p%length_scale
    end if

  end subroutine compute_derived_params

  ! ------------------------------------------------------------
  ! Initialize simulation state object: 
  !  - Given the number of particles defined in sim_params%n, allocate and initlize to 0 all the values of the allocatable arrays
  ! ------------------------------------------------------------
  subroutine init_state(p, s)
    type(sim_params), intent(in)    :: p
    type(sim_state),  intent(inout) :: s
    
    call allocate_state(p, s)
    call zero_state(s)
    
  end subroutine init_state


  ! ------------------------------------------------------------
  ! Allocate state arrays given params%n (number of particles)
  ! ------------------------------------------------------------
  subroutine allocate_state(p, s)
    type(sim_params), intent(in)    :: p
    type(sim_state),  intent(inout) :: s

    if (p%n <= 0) stop 'allocate_state(): params%n must be > 0.'

    ! Deallocate memory previously allocated for allocatable arrays.
    call deallocate_state(s)
    
    ! Allocate memory for allocatable arrays once the system size (N = number of particles) is known.
    allocate(s%rx(1:p%n), s%ry(1:p%n), s%rz(1:p%n), &
             s%vx(1:p%n), s%vy(1:p%n), s%vz(1:p%n), &
             s%ax(1:p%n), s%ay(1:p%n), s%az(1:p%n))
  end subroutine allocate_state


  ! ------------------------------------------------------------
  ! Deallocate state arrays:
  !  - Deallocate memory previously allocated for allocatable arrays.
  ! ------------------------------------------------------------
  subroutine deallocate_state(s)
    type(sim_state), intent(inout) :: s
    
    if (allocated(s%rx)) then
      deallocate(s%rx, s%ry, s%rz, s%vx, s%vy, s%vz, s%ax, s%ay, s%az)
    end if
  end subroutine deallocate_state


  ! ------------------------------------------------------------
  ! Zero all state arrays:
  !  - Initialize all state arrays to zero.
  !  - This avoids using uninitialized memory after allocation.
  ! ------------------------------------------------------------
  subroutine zero_state(s)
    type(sim_state), intent(inout) :: s

    if (allocated(s%rx)) then
      s%rx = 0.d0; s%ry = 0.d0; s%rz = 0.d0
      s%vx = 0.d0; s%vy = 0.d0; s%vz = 0.d0
      s%ax = 0.d0; s%ay = 0.d0; s%az = 0.d0
    end if
  end subroutine zero_state

end module md_types
