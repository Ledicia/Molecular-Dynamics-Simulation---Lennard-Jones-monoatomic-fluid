!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! PRECISION DEFINITION
! - Defines integer and floating-point kinds used across the project.
!******************************************************************************************************

module define_precision

    implicit none

    ! Integer kind for counters and indices.
    integer, parameter :: int_kind = selected_int_kind(9)

    ! Real kind for double precision computations.
    integer, parameter :: dp_kind  = selected_real_kind(15, 307)

end module define_precision
