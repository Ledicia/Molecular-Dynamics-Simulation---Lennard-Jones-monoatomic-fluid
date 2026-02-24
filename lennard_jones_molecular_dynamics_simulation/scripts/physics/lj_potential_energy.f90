!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! LENNARD-JONES POTENTIAL + FORCES MODULE
! - Computes Lennard-Jones (LJ) potential energy and forces from particle positions.
!
! Outputs:
! - epot    : total LJ potential energy of the system
! - d_epot  : raw first-order radial sum  sum_{i<j} [ r_ij * dU/dr_ij ] (used for virial/pressure)
! - dd_epot : raw second-order radial sum used in response functions (e.g. bulk modulus)
! - state%ax, state%ay, state%az : accelerations (forces in reduced units, mass = 1)
!
! Model:
! - LJ pair potential (reduced units): U(r) = 4 * (r^-12 - r^-6)
! - Total potential energy: epot = sum_{i<j} U(r_ij)
!
! Periodic boundary conditions (PBC) + minimum-image convention (MIC):
! - Displacements are mapped to the nearest periodic image:
!     dx_mic = dx - L * nint(dx/L)   =>   dx_mic in [-L/2, L/2]
! - This ensures distances are computed consistently in a periodic cubic box.
!
! Cutoff:
! - Only pairs with r^2 < rc^2 contribute to energy/forces (truncation).
! - Optional tail corrections can be added assuming g(r)=1 for r > rc.
!******************************************************************************************************
module lj_potential_energy
  use define_precision, only: dp_kind, int_kind
  use md_types,         only: sim_params, sim_state, pi
  use geometry_pbc,     only: minimum_image
  implicit none

  ! ------------------------------------------------------------
  ! Switch for tail corrections (mean-field approximation beyond rc)
  ! ------------------------------------------------------------
  logical, parameter :: use_tail_corrections = .true.


contains

  ! ------------------------------------------------------------
  ! Compute LJ potential energy and forces under PBC (MIC)
  !  - Resets accelerations to zero and accumulates pair forces
  !  - Returns epot and legacy derivative sums (d_epot, dd_epot)
  ! ------------------------------------------------------------
  subroutine compute_lj_potential_energy(params, state, epot, d_epot, dd_epot)
    type(sim_params), intent(in)    :: params   ! Simulation parameters (N, L, 1/L, V, rc, rc^2, ...)
    type(sim_state),  intent(inout) :: state    ! Simulation state (positions in, accelerations out)

    real(kind=dp_kind), intent(out) :: epot     ! Total potential energy
    real(kind=dp_kind), intent(out) :: d_epot   ! Normalized sum of the first-order radial derivative of the epot
    real(kind=dp_kind), intent(out) :: dd_epot  ! Normalized sum of the second-order radial derivative of the epot (response functions)

    ! -----------------------
    ! Local variables
    ! -----------------------
    integer(kind=int_kind) :: i, j              ! Particle indices

    real(kind=dp_kind) :: xi, yi, zi            ! Position of particle i
    real(kind=dp_kind) :: dx, dy, dz            ! Displacement components (MIC)
    real(kind=dp_kind) :: rij2                  ! Squared distance r_ij^2

    real(kind=dp_kind) :: inv_r2, inv_r6, inv_r12
    real(kind=dp_kind) :: dU_r                  ! Stores (-2*r^-12 + r^-6) (prefactor restored later)
    real(kind=dp_kind) :: fx, fy, fz            ! Force components (prefactor restored later)
    
    ! Tail correction terms for r_ij > rc
    real(kind=dp_kind) :: tail_factor
    real(kind=dp_kind) :: tail_corr_epot
    real(kind=dp_kind) :: tail_corr_d_epot
    real(kind=dp_kind) :: tail_corr_dd_epot
    real(kind=dp_kind) :: npd                   ! Convert number of particles N from integer to double precision for real-valued formulas

    ! ------------------------------------------------------------
    ! Minimal safety checks
    ! ------------------------------------------------------------
    if (params%n <= 0_int_kind) stop 'compute_lj_potential_energy(): params%n must be > 0.'
    if (params%box_length <= 0.d0) stop 'compute_lj_potential_energy(): params%box_length must be > 0.'
    if (params%volume <= 0.d0) stop 'compute_lj_potential_energy(): params%volume must be > 0.'
    if (params%rc <= 0.d0) stop 'compute_lj_potential_energy(): params%rc must be > 0.'
    if (params%rc_square <= 0.d0) stop 'compute_lj_potential_energy(): params%rc_square must be > 0.'
    if (.not. allocated(state%rx)) stop 'compute_lj_potential_energy(): state arrays are not allocated.'

    ! ------------------------------------------------------------
    ! Initialize accumulators
    ! ------------------------------------------------------------
    epot    = 0.d0
    d_epot  = 0.d0
    dd_epot = 0.d0

    state%ax(1:params%n) = 0.d0
    state%ay(1:params%n) = 0.d0
    state%az(1:params%n) = 0.d0


    ! -------------------------------------------------------------------------
    ! Pair interaction loop (i < j) (O(N^2)): Lennard-Jones potential + forces with PBC
    ! - Apply MIC to dx,dy,dz
    ! - If r^2 < rc^2, accumulate LJ energy and forces
    !
    ! Total potential energy: U = sum_{i<j} U_LJ(r_ij)
    !   - Lennard-Jones pair potential (reduced units): U_LJ(r) = 4 [ r^{-12} - r^{-6} ]
    !   - r_ij = | r_i - r_j |   (after applying minimum-image convention)
    !
    ! Implementation detail:
    ! - Accumulate energy without factor 4 and forces without factor 24.
    !   Prefactors are restored at the end.
    ! -------------------------------------------------------------------------
    do i = 1, params%n - 1

      ! Position of particle i
      xi = state%rx(i)
      yi = state%ry(i)
      zi = state%rz(i)

      do j = i + 1, params%n

        ! Raw displacement vector between particles i and j: r_ij = r_i - r_j
        dx = xi - state%rx(j)
        dy = yi - state%ry(j)
        dz = zi - state%rz(j)

        ! Minimum-image convention (MIC): dx = dx - L * nint(dx / L)
        dx = minimum_image(dx, params%box_length, params%inv_box_length)
        dy = minimum_image(dy, params%box_length, params%inv_box_length)
        dz = minimum_image(dz, params%box_length, params%inv_box_length)

        ! Squared interparticle distance: r_ij^2 = dx^2 + dy^2 + dz^2
        rij2 = dx*dx + dy*dy + dz*dz

        ! Cutoff condition: Only interactions with r_ij < r_c are included
        if (rij2 < params%rc_square) then

          ! Inverse powers of distance: r^{-2}, r^{-6}, r^{-12}
          inv_r2  = 1.d0 / rij2
          inv_r6  = inv_r2 * inv_r2 * inv_r2
          inv_r12 = inv_r6 * inv_r6

          ! Lennard-Jones potential contribution (without prefactor 4): U_LJ(r) / 4 = r^{-12} - r^{-6}
          epot = epot + (inv_r12 - inv_r6)

          ! Radial derivative term (without prefactor 24): dU_r = (dU/dr) * r / 24 = ([ -2 r^{-13} + r^{-7} ]) * r =  -2 * r^{-12} + r^{-6}
          dU_r = (-2.d0 * inv_r12 + inv_r6)

          ! Force vector (without prefactor 24): F_ij/24 = [- (dU/dr) * r_hat] / 24 = - (dU_r) * (r_vec / r^2)
          !    dU_r = [(dU/dr) * r] / 24
          !    r_hat = r/|r| => unitary direction of r vector
          fx = -dU_r * dx * inv_r2
          fy = -dU_r * dy * inv_r2
          fz = -dU_r * dz * inv_r2

          ! Newton's 3rd law (action-reaction): F_i += F_ij & F_j -= F_ij => Accelerations are updated since m = 1.
          state%ax(i) = state%ax(i) + fx
          state%ay(i) = state%ay(i) + fy
          state%az(i) = state%az(i) + fz

          state%ax(j) = state%ax(j) - fx
          state%ay(j) = state%ay(j) - fy
          state%az(j) = state%az(j) - fz

          ! Accumulators for thermodynamic derivatives:
          ! - Some macroscopic thermodynamic quantities require radial derivatives of the pair potential (r * dU/dr).
          !
          ! Accumulate auxiliary quantities d_epot and dd_epot:
          !  1) d_epot = sum_{i<j}[r*dU/dr]
          !       - Used to compute fundamental configurational observables in molecular dynamics. E.g.:
          !         - Virial:     W = sum_{i<j} ( r_ij · F_ij ) = - sum_{i<j} ( r_ij * dU/dr_ij )
          !         - Pressure:   P = (N k_B T)/V + W / (3V)
          !         - Equation of state P(ρ, T), through the configurational contribution.
          !         - Thermodynamic averages related to mechanical stability.
          !  2) dd_epot = sum_{i<j} [ r_ij^2 * d^2U/dr_ij^2 ]
          !       - Used to compute thermodynamic response functions. E.g.
          !         - Bulk modulus (isothermal/adiabatic compressibility)
          !           K_conf = (1 / 9V^2) * [ sum r_ij^2 * d^2U/dr_ij^2  - 2 sum r_ij * dU/dr_ij ] = (1 / 9V^2) * [ dd_epot - 2*d_epot ]
          !
          ! As with d_epot, constant prefactors are omitted inside the pair loop and restored afterwards for efficiency.
          d_epot  = d_epot  + dU_r
          dd_epot = dd_epot + (26.d0*inv_r12 - 7.d0*inv_r6)

        end if

      end do
    end do

    ! ------------------------------------------------------------
    ! Restore prefactors
    ! ------------------------------------------------------------
    epot       =  4.d0 * epot
    state%ax   = 24.d0 * state%ax
    state%ay   = 24.d0 * state%ay
    state%az   = 24.d0 * state%az
    d_epot     = 24.d0 * d_epot
    dd_epot    = 24.d0 * dd_epot
    
    ! ------------------------------------------------------------
    ! Tail corrections (optional)
    ! - Mean-field correction for the truncated Lennard-Jones potential
    ! - Accounts for interactions beyond the cutoff radius rc
    ! - Assumes g(r) = 1 for r > rc (homogeneous fluid approximation)
    ! - Corrections are added to:
    !     * epot    : potential energy
    !     * d_epot  : virial-related sum
    !     * dd_epot : higher-order thermodynamic derivative
    ! ------------------------------------------------------------
    if (use_tail_corrections) then
      npd = dble(params%n)    ! Number of particles as double precision (needed for real-valued formulas)

      tail_factor = 8.d0 * pi * (npd*npd) / (params%volume * params%rc**3)   ! Common prefactor appearing in all tail corrections: proportional to particle density squared and inverse cutoff volume
      
      tail_corr_epot = tail_factor * ( (1.d0 / (3.d0*params%rc**6)) - 1.d0 ) / 3.d0    ! Integral of LJ potential from rc to infinity assuming g(r)=1
      tail_corr_d_epot = 2.d0 * tail_factor * ( -2.d0 / (3.d0*params%rc**6) + 1.d0 )   ! Corresponds to the missing contribution to the virial (first radial derivative term)
      tail_corr_dd_epot = 2.d0 * tail_factor * ( 26.d0 / (3.d0*params%rc**6) - 7.d0 )  ! Related to second radial derivative of the potential (used in response functions such as bulk modulus)

    else
      ! No tail corrections: explicitly set all contributions to zero
      tail_corr_epot    = 0.d0
      tail_corr_d_epot  = 0.d0
      tail_corr_dd_epot = 0.d0
    end if

    epot    = epot    + tail_corr_epot
    d_epot  = d_epot  + tail_corr_d_epot
    dd_epot = dd_epot + tail_corr_dd_epot

  end subroutine compute_lj_potential_energy

end module lj_potential_energy
