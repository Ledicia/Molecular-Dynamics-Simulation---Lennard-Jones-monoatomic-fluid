!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! VERLET MODULE (velocity-Verlet integrator)
! - Advances the system by one time step (frame) using the velocity-Verlet scheme.
!
! Algorithm (velocity-Verlet): updates r, v, and a from t to t+dt
!   1) r(t+dt)   = r(t) + v(t)*dt + (1/2)*a(t)*dt^2+ O(dt^3)  => Taylor expansion of position eq. r(t+dt) = r(t) + integral_{t}^{t+dt}(v(t')dt')
!   2) v(t+dt/2) = v(t) + (1/2)*a(t)*dt                       => Half of the Taylor expansion of velocity eq.
!   3) a(t+dt)   = F(r(t+dt)) / m                             => Computed from LJ potential
!   4) v(t+dt)   = v(t+dt/2) + (1/2)*a(t+dt)*dt               => Other half of the Taylor expansion of velocity eq.
!
!   NOTE: Velocity update from a second-order Taylor expansion: v(t+dt) = v(t) + a(t)*dt + (1/2)*a_dot(t)*dt^2
!     - Expanding the acceleration: a(t+dt) = a(t) + a_dot(t)*dt + O(dt^2)  =>  a_dot(t)*dt ≈ a(t+dt) - a(t)
!     - Substituting into the velocity expansion: v(t+dt) ≈ v(t) + a(t)*dt + (1/2)[a(t+dt) - a(t)]*dt
!     - Final form: v(t+dt) = v(t) + (1/2)[a(t) + a(t+dt)]*dt
!
! Notes:
! - Particle data are stored in:
!     type(sim_state)  :: state  => rx,ry,rz,vx,vy,vz,ax,ay,az
! - Time-step combinations are stored in:
!     type(sim_params) :: params => dt, dt_half, dt_square_half
! - Positions are wrapped into [0, L) for clean trajectories.
!******************************************************************************************************
module verlet
  use define_precision,    only: dp_kind, int_kind
  use md_types,            only: sim_params, sim_state
  use lj_potential_energy, only: compute_lj_potential_energy
  use geometry_pbc,        only: wrap_positions
  implicit none

contains

  ! ------------------------------------------------------------
  ! Advance the system by one velocity-Verlet step
  !  - Updates positions and velocities in-place
  !  - Recomputes accelerations and potential energy after moving particles
  !  - Returns instantaneous energies (epot, ekin)
  ! ------------------------------------------------------------
  subroutine verlet_step(params, state, epot, ekin, d_epot, dd_epot)
    type(sim_params), intent(in)    :: params   ! Simulation parameters (N, L, dt, dt_half, dt_square_half, rc, ...)
    type(sim_state),  intent(inout) :: state    ! Simulation state (r, v, a arrays updated in-place)

    real(kind=dp_kind), intent(out) :: epot     ! Potential energy after force evaluation at t+dt
    real(kind=dp_kind), intent(out) :: ekin     ! Kinetic energy after full velocity update
    real(kind=dp_kind), intent(out) :: d_epot   ! First derivative of the epot
    real(kind=dp_kind), intent(out) :: dd_epot  ! Second derivative of the epot

    ! Minimal safety checks
    if (params%n <= 0_int_kind) stop 'verlet_step(): params%n must be > 0.'
    if (.not. allocated(state%rx)) stop 'verlet_step(): state arrays are not allocated.'

    ! ------------------------------------------------------------
    ! 1) Position update using current accelerations a(t) (Taylor expansion)
    !    r(t+dt) = r(t) + v(t)*dt + (1/2)*a(t)*dt^2
    ! ------------------------------------------------------------
    state%rx(1:params%n) = state%rx(1:params%n) + state%vx(1:params%n)*params%dt + state%ax(1:params%n)*params%dt_square_half
    state%ry(1:params%n) = state%ry(1:params%n) + state%vy(1:params%n)*params%dt + state%ay(1:params%n)*params%dt_square_half
    state%rz(1:params%n) = state%rz(1:params%n) + state%vz(1:params%n)*params%dt + state%az(1:params%n)*params%dt_square_half

    ! Keep coordinates inside [0, L) (see geometry_pbc.f90 for more information)
    call wrap_positions(state%rx, state%ry, state%rz, params%box_length)

    ! ------------------------------------------------------------
    ! 2) Half-step velocity update: v(t+dt/2)
    !    v(t+dt) is unknown , but it is possible to compute v(t+dt/2) using Taylor approximation: v(t+dt/2) = v(t) + (1/2)*a(t)*dt
    !
    ! Forces at t+dt are still unknown, so velocities can only be advanced by dt/2 using the current acceleration a(t). 
    ! The second half-step will be applied after computing a(t+dt).
    ! ------------------------------------------------------------
    state%vx(1:params%n) = state%vx(1:params%n) + state%ax(1:params%n)*params%dt_half
    state%vy(1:params%n) = state%vy(1:params%n) + state%ay(1:params%n)*params%dt_half
    state%vz(1:params%n) = state%vz(1:params%n) + state%az(1:params%n)*params%dt_half

    ! ------------------------------------------------------------
    ! 3) Compute new accelerations a(t+dt) from updated positions
    !    Evaluate the potential energy and interparticle forces, and update accelerations a = F/m.
    ! ------------------------------------------------------------
    call compute_lj_potential_energy(params, state, epot, d_epot, dd_epot)

    ! ------------------------------------------------------------
    ! 4) Finish velocity update: v(t+dt)
    !    v(t+dt) = v(t+dt/2) + (1/2)*a(t+dt)*dt
    ! ------------------------------------------------------------
    state%vx(1:params%n) = state%vx(1:params%n) + state%ax(1:params%n)*params%dt_half
    state%vy(1:params%n) = state%vy(1:params%n) + state%ay(1:params%n)*params%dt_half
    state%vz(1:params%n) = state%vz(1:params%n) + state%az(1:params%n)*params%dt_half

    ! ------------------------------------------------------------
    ! Kinetic energy after full velocity update (reduced units, m = 1)
    ! ------------------------------------------------------------
    ekin = 0.5d0 * ( sum(state%vx(1:params%n)*state%vx(1:params%n)) + &
                    sum(state%vy(1:params%n)*state%vy(1:params%n)) + &
                    sum(state%vz(1:params%n)*state%vz(1:params%n)) )

  end subroutine verlet_step

end module verlet
