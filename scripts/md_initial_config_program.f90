!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! INITIAL CONFIGURATION PROGRAM
! - Reads simulation parameters from inputs/input_simulation_parameters.txt (block-style file)
! - Builds an initial FCC configuration inside a cubic box
! - Assigns initial velocities and relaxes the configuration (warmup)
! - Writes outputs/rv_init.dat (binary snapshot with initial positions r and velocities v)
!
! Notes:
! - Particle data are stored in:
!     type(sim_state) :: state   => rx,ry,rz,vx,vy,vz,ax,ay,az
! - Simulation parameters are stored in:
!     type(sim_params) :: params => N, L, dt, rc, etc.
!******************************************************************************************************
program md_initial_config_program

  use define_precision,      only: dp_kind, int_kind
  use md_types,              only: sim_params, sim_state, init_state
  use read_input_files,      only: read_simulation_parameters
  use random_numbers,        only: random_uniform
  use lj_potential_energy,   only: compute_lj_potential_energy
  use verlet,                only: verlet_step  ! Velocity-Verlet time integrator (advances the system by one time step integrating the movement equations of the particles)

  implicit none

  ! -----------------------
  ! Simulation containers
  ! -----------------------
  type(sim_params) :: params   ! Run configuration (N, L, dt, rc, ...)
  type(sim_state)  :: state    ! Dynamic state (r, v, a arrays)

  ! -----------------------
  ! Input parameters (read from file)
  ! -----------------------
  integer(kind=int_kind) :: total_steps         ! Total number of MD integration steps (not used here but needed since we are reading input file)
  integer(kind=int_kind) :: warmup_steps        ! Number of equilibration steps (warmup)
  integer(kind=int_kind) :: output_interval     ! Output frequency (not used here but needed since we are reading input file))
  real(kind=dp_kind)     :: rc_over_L           ! Cutoff ratio rc/L
  real(kind=dp_kind)     :: target_total_energy ! Target total energy for velocity rescaling

  ! -----------------------
  ! Local variables
  ! -----------------------
  integer(kind=int_kind) :: step              ! MD time-step counter (used in warmup loop)
  integer(kind=int_kind) :: seed              ! Seed for the random number generator (initial velocities)
  real(kind=dp_kind)     :: lattice_spacing   ! FCC lattice spacing: a = L / num_cells

  real(kind=dp_kind) :: epot, ekin, etot
  real(kind=dp_kind) :: d_epot, dd_epot       ! First and Second derivative of Lennard Jones (LJ) epot


  ! ------------------------------------------------------------
  ! 1) Read parameters (block-format file)
  !    This routine from read_input_files.f90 builds "params" and returns control parameters.
  ! ------------------------------------------------------------
  call read_simulation_parameters('inputs/input_simulation_parameters.txt', &
                                  params, total_steps, output_interval, warmup_steps, &
                                  rc_over_L, target_total_energy)

  ! ------------------------------------------------------------
  ! 2) Allocate and initialize the simulation state in memory
  !    - init_state() from md_types.f90 performs: allocate_state() + zero_state()
  ! ------------------------------------------------------------
  call init_state(params, state)

  ! ------------------------------------------------------------
  ! 3) Build FCC lattice positions in [0, L)
  !    FCC basis in a cubic cell:
  !      (0,0,0), (0,1/2,1/2), (1/2,0,1/2), (1/2,1/2,0)
  ! ------------------------------------------------------------
  lattice_spacing = params%box_length / dble(params%num_cells)
  call build_fcc_lattice(params, state, lattice_spacing)

  ! ------------------------------------------------------------
  ! 4) Initial velocities
  !  - Assign random velocities (uniform in [-0.5, 0.5]) to each component.
  !  - Then remove center-of-mass velocity to enforce zero total momentum:
  !      sum_i v_i = 0
  ! ------------------------------------------------------------
  seed = -12345_int_kind
  call assign_random_velocities(params, state, seed)
  call remove_center_of_mass_velocity(params, state)

  ! ------------------------------------------------------------
  ! 5) Initial forces and energies at t = 0
  !  - Compute LJ forces/accelerations and potential energy
  !  - Compute kinetic energy and total energy
  ! ------------------------------------------------------------
  call compute_lj_potential_energy(params, state, epot, d_epot, dd_epot)
  ekin = 0.5d0 * sum(state%vx*state%vx + state%vy*state%vy + state%vz*state%vz)
  etot = epot + ekin

  ! ------------------------------------------------------------
  ! 6) Rescale velocities (and therefore E_kin) to match target total energy
  !  - Impose E_tot_target = E_pot + E_kin_new
  !  - Therefore: E_kin_new = E_tot_target - E_pot
  !  - Rescale uniformly: v_i = alpha * v_i,  alpha = sqrt(E_kin_new / E_kin_old)
  ! ------------------------------------------------------------
  call rescale_velocities_to_target_energy(params, state, target_total_energy, epot)

  ! Recompute energies after rescaling
  call compute_lj_potential_energy(params, state, epot, d_epot, dd_epot)
  ekin = 0.5d0 * sum(state%vx*state%vx + state%vy*state%vy + state%vz*state%vz)
  etot = epot + ekin

  ! ------------------------------------------------------------
  ! 7) Warmup / equilibration (velocity-Verlet)
  !  - The system is evolved for "warmup_steps" time steps.
  !  - No averages are accumulated here (equilibration only).
  ! ------------------------------------------------------------
  do step = 1, warmup_steps
    call verlet_step(params, state, epot, ekin, d_epot, dd_epot)

  end do

  ! ------------------------------------------------------------
  ! 8) Write snapshot for the production MD code
  ! ------------------------------------------------------------
  call write_rv_init('outputs/rv_init.dat', params, state)

contains

  ! ------------------------------------------------------------
  ! Build an FCC lattice inside the simulation box [0, L).
  ! - The box is divided into p%num_cells unit cells per direction.
  ! - Each FCC unit cell contains 4 particles located at:
  !     FCC basis: (0,0,0), (0,1/2,1/2), (1/2,0,1/2), (1/2,1/2,0)
  ! - Particle positions are stored in the state arrays (rx, ry, rz).
  ! ------------------------------------------------------------
  subroutine build_fcc_lattice(p, s, a)
    type(sim_params), intent(in)    :: p    ! Simulation parameters (contains num_cells and box length)
    type(sim_state),  intent(inout) :: s    ! Simulation state (particle positions will be written here)
    real(kind=dp_kind), intent(in)  :: a    ! Lattice spacing (unit cell length): a = L / num_cells

    integer(kind=int_kind) :: ix, iy, iz    ! Unit cell indices along x, y, z directions
    integer(kind=int_kind) :: idx           ! Global particle index
    real(kind=dp_kind)     :: x0, y0, z0    ! Origin of the current FCC unit cell

    idx = 0                                 ! Initialize particle counter

    ! Loop over all FCC unit cells in the box
    do ix = 0, p%num_cells - 1
      do iy = 0, p%num_cells - 1
        do iz = 0, p%num_cells - 1

          ! Coordinates of the lower corner of the current unit cell
          x0 = dble(ix) * a
          y0 = dble(iy) * a
          z0 = dble(iz) * a

          ! ----------------------------------------------------
          ! Place the 4 particles of the FCC basis in this cell
          ! ----------------------------------------------------
          ! P1rticle 1 (0, 0, 0)
          idx = idx + 1
          s%rx(idx) = x0
          s%ry(idx) = y0
          s%rz(idx) = z0

          ! Particle 2 (0, 1/2, 1/2)
          idx = idx + 1
          s%rx(idx) = x0
          s%ry(idx) = y0 + 0.5d0*a
          s%rz(idx) = z0 + 0.5d0*a

          ! Particle 3 (1/2, 0, 1/2)
          idx = idx + 1
          s%rx(idx) = x0 + 0.5d0*a
          s%ry(idx) = y0
          s%rz(idx) = z0 + 0.5d0*a

          ! Particle 4 (1/2, 1/2, 0)
          idx = idx + 1
          s%rx(idx) = x0 + 0.5d0*a
          s%ry(idx) = y0 + 0.5d0*a
          s%rz(idx) = z0

        end do
      end do
    end do

    ! Sanity check: total number of particles must match N = 4*num_cells^3
    if (idx /= p%n) stop 'build_fcc_lattice(): unexpected particle count.'

  end subroutine build_fcc_lattice


  ! ------------------------------------------------------------
  ! INITIAL VELOCITIES
  !  - Assign random velocities (uniform in [-0.5, 0.5]) to each component.
  !  - This provides a random initial kinetic energy and avoids any preferred
  !    direction in velocity space.
  ! ------------------------------------------------------------
  subroutine assign_random_velocities(p, s, seed)
    type(sim_params), intent(in)    :: p          ! Simulation parameters (contains number of particles n = N)
    type(sim_state),  intent(inout) :: s          ! Simulation state (velocity arrays will be modified)
    integer(kind=int_kind), intent(inout) :: seed ! Seed for the random number generator (updated internally)

    integer(kind=int_kind) :: i                   ! Particle index

    do i = 1, p%n                                 ! Loop over all particles
      s%vx(i) = random_uniform(seed) - 0.5d0      ! Random x-velocity in the interval (-0.5, 0.5)
      s%vy(i) = random_uniform(seed) - 0.5d0      ! Random y-velocity in the interval (-0.5, 0.5)
      s%vz(i) = random_uniform(seed) - 0.5d0      ! Random z-velocity in the interval (-0.5, 0.5)
    end do

  end subroutine assign_random_velocities


  ! ------------------------------------------------------------
  ! REMOVE CENTER-OF-MASS VELOCITY
  !  - Enforce zero total linear momentum: sum_i v_i = 0
  !    - Center-of-mass velocity:  v_cm = (1/N) sum_i v_i
  !    - Shift velocities:         v_i  = v_i - v_cm
  ! ------------------------------------------------------------
  subroutine remove_center_of_mass_velocity(p, s)
    type(sim_params), intent(in)    :: p      ! Simulation parameters (contains number of particles n = N)
    type(sim_state),  intent(inout) :: s      ! Simulation state (velocity arrays will be modified)

    real(kind=dp_kind) :: vx_cm               ! Center-of-mass velocity component along x
    real(kind=dp_kind) :: vy_cm               ! Center-of-mass velocity component along y
    real(kind=dp_kind) :: vz_cm               ! Center-of-mass velocity component along z

    vx_cm = sum(s%vx(1:p%n)) / dble(p%n)      ! Compute x-component of center-of-mass velocity
    vy_cm = sum(s%vy(1:p%n)) / dble(p%n)      ! Compute y-component of center-of-mass velocity
    vz_cm = sum(s%vz(1:p%n)) / dble(p%n)      ! Compute z-component of center-of-mass velocity

    s%vx(1:p%n) = s%vx(1:p%n) - vx_cm         ! Subtract center-of-mass velocity from all x-velocities
    s%vy(1:p%n) = s%vy(1:p%n) - vy_cm         ! Subtract center-of-mass velocity from all y-velocities
    s%vz(1:p%n) = s%vz(1:p%n) - vz_cm         ! Subtract center-of-mass velocity from all z-velocities

  end subroutine remove_center_of_mass_velocity


  ! -------------------------------------------------------------------------
  ! VELOCITY RESCALING TO MATCH TARGET TOTAL ENERGY
  ! - Impose a desired total energy:
  !     E_tot_target = E_pot + E_kin_new  =>  E_kin_new = E_tot_target - E_pot
  ! - Rescale uniformly:
  !     v_i = alpha * v_i,  alpha = sqrt(E_kin_new / E_kin_old)
  ! -------------------------------------------------------------------------
  subroutine rescale_velocities_to_target_energy(p, s, target_energy, epot)
    type(sim_params), intent(in)    :: p              ! Simulation parameters (contains number of particles N)
    type(sim_state),  intent(inout) :: s              ! Simulation state (velocity arrays will be rescaled)
    real(kind=dp_kind), intent(in)  :: target_energy  ! Desired total energy E_tot_target
    real(kind=dp_kind), intent(in)  :: epot           ! Current potential energy E_pot

    real(kind=dp_kind) :: ekin_old                    ! Kinetic energy before rescaling
    real(kind=dp_kind) :: ekin_new                    ! Target kinetic energy after rescaling
    real(kind=dp_kind) :: scale                       ! Velocity rescaling factor

    ekin_old = 0.5d0 * sum(s%vx*s%vx + s%vy*s%vy + s%vz*s%vz)  ! Current kinetic energy: (1/2) sum_i |v_i|^2
    ekin_new = target_energy - epot                            ! Required kinetic energy to reach target total energy

    if (ekin_new <= 0.d0) stop 'rescale_velocities_to_target_energy(): target energy too low (zero or negative kinetic).'
    if (ekin_old <= 0.d0) stop 'rescale_velocities_to_target_energy(): ekin_old <= 0 (cannot rescale).'
    
    scale = sqrt(ekin_new / ekin_old)             ! Compute uniform velocity scaling factor

    s%vx(1:p%n) = s%vx(1:p%n) * scale             ! Rescale x-components of all particle velocities
    s%vy(1:p%n) = s%vy(1:p%n) * scale             ! Rescale y-components of all particle velocities
    s%vz(1:p%n) = s%vz(1:p%n) * scale             ! Rescale z-components of all particle velocities

  end subroutine rescale_velocities_to_target_energy


  ! ------------------------------------------------------------
  ! WRITE INITIAL CONFIGURATION SNAPSHOT (binary, unformatted)
  !  - File outputs/rv_init.dat contains:
  !     Record 1: rx, ry, rz  (particle positions)
  !     Record 2: vx, vy, vz  (particle velocities)
  ! ------------------------------------------------------------
  subroutine write_rv_init(filename, p, s)
    character(len=*), intent(in) :: filename   ! Name of the output file to store the initial configuration
    type(sim_params), intent(in) :: p          ! Simulation parameters (contains number of particles N)
    type(sim_state),  intent(in) :: s          ! Simulation state (positions and velocities to be written)
    integer :: iu, ios                         ! File unit number and I/O status flag

    iu = 2001                                                                                       ! Assign a file unit number for output
    open(unit=iu, file=filename, form='unformatted', status='replace', action='write', iostat=ios)  ! Open binary (unformatted) file for writing
    if (ios /= 0) stop 'write_rv_init(): cannot open output file.'                                  ! Abort if file cannot be opened

    write(iu) s%rx(1:p%n), s%ry(1:p%n), s%rz(1:p%n)  ! Write particle positions (x, y, z) as first binary record
    write(iu) s%vx(1:p%n), s%vy(1:p%n), s%vz(1:p%n)  ! Write particle velocities (x, y, z) as second binary record

    close(iu)    ! Close the output file
    
  end subroutine write_rv_init

end program md_initial_config_program
