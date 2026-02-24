!******************************************************************************************************
! Ledicia Díaz Lago
! Physical Materials Simulation 2018/2019
!
! READ SIMULATION INPUT (block-style text file)
! - Reads the human-friendly input file "inputs/input_simulation_parameters.txt" written in blocks,
!   skipping comments and headers, and extracting ONLY the numeric lines.
!
! Expected numeric blocks (in this order):
!   - Block 1 numeric line:  k  total_steps  output_interval  warmup_steps
!   - Block 2 numeric line:  dt  L  rc_over_L (rc/L)
!   - Block 3 numeric line:  target_total_energy
!
! This routine:
!   - Builds a sim_params object (md_types.f90) from the parsed values.
!   - Assumes an FCC lattice with n = N = 4*k^3 particles.
!   - Converts rc_over_L to an absolute cutoff radius: rc = rc_over_L * L.
!******************************************************************************************************
module read_input_files
  use define_precision, only: dp_kind, int_kind
  use md_types,         only: sim_params, init_params
  implicit none
  
  
contains

  subroutine read_simulation_parameters(filename, params, total_steps, output_interval,  &
                                      warmup_steps, rc_over_L, target_total_energy)
    implicit none
    
    character(len=*), intent(in)  :: filename   ! Input file name containing the simulation parameters
    type(sim_params), intent(out) :: params     ! Structure holding all simulation parameters (box, dt, rc, N, ...)

    ! Global MD control parameters
    integer(kind=int_kind), intent(out) :: total_steps        ! Total number of MD integration steps
    integer(kind=int_kind), intent(out) :: output_interval    ! Frequency of output writing
    integer(kind=int_kind), intent(out) :: warmup_steps       ! Number of equilibration steps

    ! Additional physical/control parameters
    real(kind=dp_kind), intent(out) :: rc_over_L              ! Cutoff radius in units of box length (rc / L)
    real(kind=dp_kind), intent(out) :: target_total_energy    ! Target total energy (e.g. for validation)
 
    ! File handling variables
    integer :: iu   ! Input file unit number
    integer :: ios  ! I/O status flag (used to detect read/open errors)

    ! Line buffer for reading the input file line by line
    character(len=512) :: line

    ! Temporary variables used while parsing the input file blocks
    integer(kind=int_kind) :: k_read           ! Number of FCC cells per edge
    real(kind=dp_kind)     :: dt_read          ! Time step read from input
    real(kind=dp_kind)     :: L_read           ! Simulation box length read from input
    real(kind=dp_kind)     :: rc_over_L_read   ! Cutoff radius scaled by box length
    integer(kind=int_kind) :: n_read           ! Total number of particles (derived from k_read)

    ! Flags indicating whether each input block has been successfully read
    logical :: got_block1  ! Block 1: lattice size and MD step counters
    logical :: got_block2  ! Block 2: dt, box length, cutoff
    logical :: got_block3  ! Block 3: target energy

    got_block1 = .false.
    got_block2 = .false.
    got_block3 = .false.

    ! Initialization of output variables to avoid using undefined values if input reading fails
    total_steps         = 0_int_kind
    output_interval     = 0_int_kind
    warmup_steps        = 0_int_kind
    rc_over_L           = 0.d0
    target_total_energy = 0.d0

    ! Initialization of simulation parameters structure
    params%n          = 0_int_kind     ! Number of particles
    params%num_cells  = 0_int_kind     ! FCC cells per edge
    params%box_length = 0.d0           ! Simulation box length
    params%dt         = 0.d0           ! Time step
    params%rc         = 0.d0           ! Cutoff radius

    ! ------------------------------------------------------------
    ! Open input file for reading
    ! ------------------------------------------------------------
    iu = 1020
    open(unit=iu, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'read_simulation_parameters(): cannot open input file.'

    do
      read(iu,'(A)', iostat=ios) line
      if (ios /= 0) exit               ! EOF
      if (len_trim(line) == 0) cycle   ! Skip empty lines      
      if (line(1:1) == '#') cycle      ! Skip comment lines (comments start at column 1 with '#')

      ! -------------------------
      ! Block 1: 4 integer values
      ! -------------------------
      if (.not. got_block1) then
        ! Initialization of output variables to avoid using undefined values if input reading fails and for sanity checks
        k_read          = 0_int_kind
        total_steps     = 0_int_kind
        output_interval = 0_int_kind
        warmup_steps    = 0_int_kind

        read(line, *, iostat=ios) k_read, total_steps, output_interval, warmup_steps  ! Read values in line and assign them to variables in list
        if (ios == 0) then
          if (k_read <= 0_int_kind)          stop 'read_simulation_parameters(): k must be > 0.'
          if (total_steps <= 0_int_kind)     stop 'read_simulation_parameters(): total_steps must be > 0.'
          if (output_interval <= 0_int_kind) stop 'read_simulation_parameters(): output_interval must be > 0.'
          if (warmup_steps < 0_int_kind)     stop 'read_simulation_parameters(): warmup_steps must be >= 0.'

          got_block1 = .true.
          cycle
        else
          ios = 0
          cycle
        end if
      end if

      ! -------------------------
      ! Block 2: 3 real values
      ! -------------------------
      if (.not. got_block2) then        
        ! Initialization of output variables to avoid using undefined values if input reading fails and for sanity checks
        dt_read         = 0.d0
        L_read          = 0.d0
        rc_over_L_read  = 0.d0

        read(line, *, iostat=ios) dt_read, L_read, rc_over_L_read  ! Read values in line and assign them to variables in list
        if (ios == 0) then
          if (dt_read <= 0.d0)         stop 'read_simulation_parameters(): dt must be > 0.'
          if (L_read  <= 0.d0)         stop 'read_simulation_parameters(): L must be > 0.'
          if (rc_over_L_read <= 0.d0)  stop 'read_simulation_parameters(): rc_over_L must be > 0.'
          if (rc_over_L_read >  0.5d0) stop 'read_simulation_parameters(): rc_over_L must be <= 0.5 (minimum image).'

          got_block2  = .true.
          rc_over_L   = rc_over_L_read
          cycle
        else
          ios = 0
          cycle
        end if
      end if

      ! -------------------------
      ! Block 3: 1 real value
      ! -------------------------
      if (.not. got_block3) then  
        ! Initialization of output variables to avoid using undefined values if input reading fails and for sanity checks
        target_total_energy = 0.d0
        read(line, *, iostat=ios) target_total_energy  ! Read values in line and assign them to variables in list
        if (ios == 0) then
          got_block3 = .true.
          exit
        else
          ios = 0
          cycle
        end if
      end if

    end do

    close(iu)

    if (.not. got_block1) stop 'read_simulation_parameters(): missing Block 1 numeric line.'
    if (.not. got_block2) stop 'read_simulation_parameters(): missing Block 2 numeric line.'
    if (.not. got_block3) stop 'read_simulation_parameters(): missing Block 3 numeric line.'

    ! FCC lattice particle count: n = N = 4*k^3
    n_read = 4_int_kind * k_read * k_read * k_read

    ! Initialize parameters (function defined in md_types.f90) and convert rc/L to absolute cutoff radius rc = (rc_over_L) * L
    call init_params(params, n_read, L_read, dt_read, rc_over_L * L_read, num_cells=k_read)

  end subroutine read_simulation_parameters

end module read_input_files
