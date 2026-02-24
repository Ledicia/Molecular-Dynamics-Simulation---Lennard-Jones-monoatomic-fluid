!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! Thermodynamic coefficients computed from time-averaged quantities accumulated during a
! Molecular Dynamics (MD) microcanonical (NVE) production run.
!
! Inputs (time averages over MD samples, computed in md_means.f90):
!   - ekin_mean              = <Ekin>
!   - press_mean             = <P>
!   - ekinInv_mean           = <1/Ekin>
!   - d_epot_mean            = <d_epot>             where d_epot  = sum_{i<j} r_ij dU/dr_ij
!   - dd_epot_mean           = <dd_epot>            where dd_epot = sum_{i<j} r_ij^2 d^2U/dr_ij^2
!   - d_epot_ekinInv_mean    = <d_epot * 1/Ekin>
!   - d_epot2_ekinInv_mean   = <d_epot^2 * 1/Ekin>
!
! Outputs (reduced units, kB=1, m=1):
!   - temperature            : T = 2 <Ekin> / (3 N)
!   - pressure               : P = <P> (passed through for convenience)
!
!   - Ca_v, Ce_v             : heat capacity at constant volume (total / per particle)
!   - Ca_p, Ce_p             : heat capacity at constant pressure (total / per particle)
!
!   - gamma                  : Grüneisen parameter
!
!   - K_S                    : adiabatic (isentropic) bulk modulus
!   - K_S_inv (kappa_S)      : isentropic compressibility  (K_S_inv = 1/K_S)
!   - K_T                    : isothermal bulk modulus
!   - K_T_inv (kappa_T)      : isothermal compressibility  (K_T_inv = 1/K_T)
!
!   - alpha_E1, alpha_E2     : thermal expansion at constant energy (two estimators that should agree)
!   - alpha_S                : thermal expansion at constant entropy  (adiabatic)
!   - alpha_p                : thermal expansion at constant pressure (isobaric)
!
! Notes:
! - Degrees of freedom: f = 3N - 3 (removes center-of-mass momentum conservation).
! - Notation:
!     K_*     : bulk modulus
!     K_*_inv : compressibility (inverse bulk modulus), i.e. kappa_* = K_*_inv = 1/K_*
!******************************************************************************************************

module thermodynamic_coefs

  use define_precision, only: dp_kind, int_kind
  use md_types,         only: sim_params          ! Parameters of the simulation defined in inputs/input_simulation_parameters.txt

  implicit none
  private

  public :: thermodynamic_coefs_results
  public :: thermodynamic_coefs_compute

  type :: thermodynamic_coefs_results
    real(kind=dp_kind) :: npd                 = 0.0_dp_kind  ! Number of particles in double precision for real valued formulas
    real(kind=dp_kind) :: degrees_of_freedom  = 0.0_dp_kind  ! f = 3*npd − 3 (remove center-of-mass momentum conservation)
    real(kind=dp_kind) :: aux1                = 0.0_dp_kind  ! Aux factor to simplify algebra                                   => aux1 = 1 − 2/f
    real(kind=dp_kind) :: aux2                = 0.0_dp_kind  ! Aux factor to simplify algebra                                   => aux2 = f/2 − 1

    real(kind=dp_kind) :: temperature         = 0.0_dp_kind  ! T (reduced units kB=1, m=1)                                      => T = 2 <Ekin> / (3 npd)
    real(kind=dp_kind) :: pressure            = 0.0_dp_kind  ! Mean pressure                                                    => P = <P>  (input already includes virial)

    real(kind=dp_kind) :: Ca_v                = 0.0_dp_kind  ! Total  microcanonical / adiabatic Cv at constant volume          => Ca_v = 1 / [ 1 − (1 − 2/f) <Ekin><1/Ekin> ]
    real(kind=dp_kind) :: Ce_v                = 0.0_dp_kind  ! Per-particle Cv at constant volume                               => Ce_v = Ca_v / npd

    real(kind=dp_kind) :: gamma               = 0.0_dp_kind  ! Grüneisen parameter                                              => gamma = 1/Ce_v + V*(f/2-1)*( <d><1/Ekin> - <d*(1/Ekin)> )

    real(kind=dp_kind) :: Ca_p                = 0.0_dp_kind  ! Total Cp at constant pressure                                    => Ca_p = Ca_v * (K_S / K_T)
    real(kind=dp_kind) :: Ce_p                = 0.0_dp_kind  ! Per-particle Cp at constant pressure                             => Ce_p = Ca_p / npd

    real(kind=dp_kind) :: K_T_inv             = 0.0_dp_kind  ! Isothermal compressibility                                       => kappa_T = 1/K_T
    real(kind=dp_kind) :: K_S_inv             = 0.0_dp_kind  ! Isentropic compressibility                                       => kappa_S = 1/K_S
    real(kind=dp_kind) :: K_T                 = 0.0_dp_kind  ! Isothermal bulk modulus                                          => K_T = K_S - (T*Ca_v*gamma^2)/V
    real(kind=dp_kind) :: K_S                 = 0.0_dp_kind  ! Isentropic bulk modulus                                          => See formula below (computed from means)
    real(kind=dp_kind) :: K_S_aux             = 0.0_dp_kind  ! Intermediate quantity to build K_S                               => K_S_aux = (N T /V)*(1+2γ−1/Ce_v) + V<dd>

    ! Thermal expansion coefficient (two estimators, should agree within statistical error):
    real(kind=dp_kind) :: alpha_E1            = 0.0_dp_kind  ! alpha_E (estimator 1)                                            => alpha_E1 = 1 / ( P V / Ca_v − gamma T )
    real(kind=dp_kind) :: alpha_E2            = 0.0_dp_kind  ! alpha_E (estimator 2)                                            => alpha_E2 = 1 / ( V*( (1−2/f)<Ekin><d*(1/Ekin)> − <d> ) )

    real(kind=dp_kind) :: alpha_S             = 0.0_dp_kind  ! alpha_S = (1/V)(dV/dT)_S                                         => alpha_S = -1/(gamma*T)
    real(kind=dp_kind) :: alpha_P             = 0.0_dp_kind  ! alpha_P = (1/V)(dV/dT)_P                                         => alpha_P = (Ca_v*gamma)/(V*K_T) = (Ca_v*gamma/V)*K_T_inv
  end type thermodynamic_coefs_results


contains

  subroutine thermodynamic_coefs_compute(params, ekin_mean, press_mean, ekinInv_mean, d_epot_mean, dd_epot_mean, d_epot_ekinInv_mean, d_epot2_ekinInv_mean, out)

    type(sim_params),   intent(in)  :: params                                      ! Parameters related to MD configuration (N, L, V, dt, rc, ...)
    
    ! Quantities needed for the calculation of the thermodynamic coefficients (time averages):
    real(kind=dp_kind), intent(in)  :: ekin_mean, press_mean                       ! <Ekin>, <P>
    real(kind=dp_kind), intent(in)  :: ekinInv_mean                                ! <1/Ekin>
    real(kind=dp_kind), intent(in)  :: d_epot_mean, dd_epot_mean                   ! <d_epot>, <dd_epot>
    real(kind=dp_kind), intent(in)  :: d_epot_ekinInv_mean, d_epot2_ekinInv_mean   ! <d_epot * (1/Ekin)>, <(d_epot^2) * (1/Ekin)>

    type(thermodynamic_coefs_results), intent(out) :: out

    real(kind=dp_kind) :: denom
    
    ! -----------------------
    ! Basic constants
    ! -----------------------
    out%npd                = dble(params%n)                      ! Number of particles for real valued formulas
    out%degrees_of_freedom = 3.0_dp_kind*out%npd - 3.0_dp_kind   ! f = 3N - 3

    if (out%degrees_of_freedom <= 0.0_dp_kind) stop 'thermodynamic_compute(): degrees_of_freedom <= 0 (check N).'

    out%aux1 = 1.0_dp_kind - 2.0_dp_kind/out%degrees_of_freedom  ! aux1 = 1 - 2/f
    out%aux2 = out%degrees_of_freedom/2.0_dp_kind - 1.0_dp_kind  ! aux2 = f/2 - 1

    ! ===================================
    ! Thermodynamic coefficient estimators
    ! ===================================

    ! Temperature: T = 2<Ekin>/(3N)
    !out%temperature = 2.0_dp_kind * ekin_mean / (3.0_dp_kind * out%npd)
    out%temperature = 2.0_dp_kind * ekin_mean / out%degrees_of_freedom

    ! Pressure (already computed in md_means
    out%pressure = press_mean

    ! -----------------------
    ! Heat capacities at constant volume:
    !   Ca_v = 1 / [ 1 - (1 - 2/f) <Ekin><1/Ekin> ] = 1 / [ 1 - aux1 <Ekin><1/Ekin> ]
    !   Ce_v = Ca_v / N
    ! -----------------------
    denom = 1.0_dp_kind - out%aux1 * ekin_mean * ekinInv_mean
    if (abs(denom) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): Ca_v denominator ~ 0 (numerical instability).'
    out%Ca_v = 1.0_dp_kind / denom

    out%Ce_v = out%Ca_v / out%npd
    if (abs(out%Ce_v) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): Ce_v ~ 0 (check inputs).'

    ! Grüneisen gamma:
    !   gamma = 1/Ce_v + (f/2 - 1)/3 * ( <d_epot><1/Ekin> - <d_epot*(1/Ekin)> )
    out%gamma = 1.0_dp_kind/out%Ce_v + (out%aux2/3.0_dp_kind) * ( d_epot_mean*ekinInv_mean - d_epot_ekinInv_mean )

    ! -----------------------
    ! Compressibility / bulk modulus estimators:
    ! -----------------------
    ! Isentropic bulk modulus K_S and isentropic compressibility kappa_S = 1/K_S
    !
    ! NOTE:
    !   d  = d_epot  = sum_{i<j} [ r_ij * dU/dr_ij ]
    !   dd = dd_epot = sum_{i<j} [ r_ij^2 * d^2U/dr_ij^2 ]
    !   f  = 3N - 3               (degrees of freedom, removes COM momentum)
    !   aux2 = f/2 - 1
    !   C_V = Ce_v                (heat capacity at constant volume per particle)
    !   gamma                     (Grüneisen parameter, computed above)
    !
    ! Born-like (configurational) contribution:
    !   K_conf,Born = ( <dd_epot> - 2<d_epot> ) / (9 V)
    !
    ! Aux part (ideal + coupling + Born):
    !   K_S_aux = (N T / V) * ( 1 + 2*gamma - 1/c_V )  +  ( <dd_epot> - 2<d_epot> ) / (9 V)
    !
    ! Microcanonical fluctuation correction (coupling to kinetic energy):
    !   K_S_fluct = - (f/2 - 1) / (9 V^2) * < d_epot^2 / K > - 2 <d_epot> * < d_epot / K > + <d_epot>^2 * < 1 / K > ]
    !
    !   K_S = K_S_aux + K_S_fluct
    ! -----------------------
    out%K_S_aux = ((out%npd * out%temperature * ( 1.0_dp_kind + 2.0_dp_kind*out%gamma - 1.0_dp_kind/out%Ce_v ) ) / params%volume) + (dd_epot_mean - 2.0_dp_kind*d_epot_mean) / (9.0_dp_kind * params%volume)
    out%K_S = out%K_S_aux - (out%aux2 * (d_epot2_ekinInv_mean - 2.0_dp_kind*d_epot_mean*d_epot_ekinInv_mean + (d_epot_mean*d_epot_mean)*ekinInv_mean)) / (9.0_dp_kind * params%volume * params%volume)

    if (abs(out%K_S) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): K_S ~ 0 (cannot invert).'
    out%K_S_inv = 1.0_dp_kind / out%K_S   ! kappa_S = 1/K_S

    ! Isothermal bulk modulus K_T and isothermal compressibility kappa_T:
    !   K_T = K_S - (T * Ca_v * gamma^2) / V
    !   kappa_T = 1/K_T
    out%K_T = out%K_S - (out%temperature * out%Ca_v * (out%gamma*out%gamma)) / params%volume
    if (abs(out%K_T) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): K_T ~ 0 (cannot invert / compute Cp, alpha_P).'
    out%K_T_inv = 1.0_dp_kind / out%K_T   ! kappa_T = 1/K_T

    ! -----------------------
    ! Heat capacities at constant pressure:
    !   Ca_p = Ca_v * (K_S / K_T)
    !   Ce_p = Ca_p / N
    ! -----------------------
    out%Ca_p = out%Ca_v * (out%K_S / out%K_T)
    out%Ce_p = out%Ca_p / out%npd

    ! -----------------------
    ! Expansion estimators:
    ! -----------------------
    ! Thermal expansion at fixed energy: alpha_E1 = 1 / ( P*V/Ca_v - gamma*T )
    denom = (out%pressure * params%volume / out%Ca_v) - (out%gamma * out%temperature)
    if (abs(denom) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): alpha_E1 denominator ~ 0.'
    out%alpha_E1 = 1.0_dp_kind / denom

    ! Thermal expansion at fixed energy: alpha_E2 = 1 / [ (1/3) * ( (1 - 2/f) <Ekin><d_epot*(1/Ekin)> - <d_epot> ) ]
    denom = (1.0_dp_kind/3.0_dp_kind) * ( out%aux1*ekin_mean*d_epot_ekinInv_mean - d_epot_mean )
    if (abs(denom) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): alpha_E2 denominator ~ 0.'
    out%alpha_E2 = 1.0_dp_kind / denom

    ! Thermal expansion at fixed entropy: alpha_S = -1/(gamma*T)
    denom = out%gamma * out%temperature
    if (abs(denom) < 1.0e-14_dp_kind) stop 'thermodynamic_compute(): gamma*T ~ 0 (alpha_S undefined).'
    out%alpha_S = -1.0_dp_kind / denom

    ! Thermal expansion at fixed pressure: alpha_P = (Ca_v * gamma) / (V * K_T) = (Ca_v*gamma/V)*K_T_inv
    out%alpha_P = (out%Ca_v * out%gamma) / params%volume * out%K_T_inv

  end subroutine thermodynamic_coefs_compute

end module thermodynamic_coefs
