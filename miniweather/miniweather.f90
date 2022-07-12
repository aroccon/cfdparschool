!! miniWeather
!! Author: Matt Norman <normanmr@ornl.gov>  , Oak Ridge National Laboratory
!! This code simulates dry, stratified, compressible, non-hydrostatic fluid flows
!! For documentation, please see the attached documentation in the "documentation" folder

program miniweather
  implicit none
  !Declare the precision for the real constants (at least 15 digits of accuracy = double precision)
  integer , parameter :: rp = selected_real_kind(15)
  !Define some physical constants to use throughout the simulation
  real(rp), parameter :: pi        = 3.14159265358979323846264338327_rp   !Pi
  real(rp), parameter :: grav      = 9.8_rp                               !Gravitational acceleration (m / s^2)
  real(rp), parameter :: cp        = 1004._rp                             !Specific heat of dry air at constant pressure
  real(rp), parameter :: cv        = 717._rp                              !Specific heat of dry air at constant volume
  real(rp), parameter :: rd        = 287._rp                              !Dry air constant for equation of state (P=rho*rd*T)
  real(rp), parameter :: p0        = 1.e5_rp                              !Standard pressure at the surface in Pascals
  real(rp), parameter :: C0        = 27.5629410929725921310572974482_rp   !Constant to translate potential temperature into pressure (P=C0*(rho*theta)**gamma)
  real(rp), parameter :: gamma     = 1.40027894002789400278940027894_rp   !gamma=cp/Rd
  !Define domain and stability-related constants
  real(rp), parameter :: xlen      = 2.e4_rp    !Length of the domain in the x-direction (meters)
  real(rp), parameter :: zlen      = 1.e4_rp    !Length of the domain in the z-direction (meters)
  real(rp), parameter :: hv_beta   = 0.05_rp    !How strong to diffuse the solution: hv_beta \in [0:1]
  real(rp), parameter :: cfl       = 1.50_rp    !"Courant, Friedrichs, Lewy" number (for numerical stability)
  real(rp), parameter :: max_speed = 450        !Assumed maximum wave speed during the simulation (speed of sound + speed of wind) (meter / sec)
  integer , parameter :: hs        = 2          !"Halo" size: number of cells beyond the MPI tasks's domain needed for a full "stencil" of information for reconstruction
  integer , parameter :: sten_size = 4          !Size of the stencil used for interpolation
  !Parameters for indexing and flags
  integer , parameter :: NUM_VARS = 4           !Number of fluid state variables
  integer , parameter :: ID_DENS  = 1           !index for density ("rho")
  integer , parameter :: ID_UMOM  = 2           !index for momentum in the x-direction ("rho * u")
  integer , parameter :: ID_WMOM  = 3           !index for momentum in the z-direction ("rho * w")
  integer , parameter :: ID_RHOT  = 4           !index for density * potential temperature ("rho * theta")
  integer , parameter :: DIR_X = 1              !Integer constant to express that this operation is in the x-direction
  integer , parameter :: DIR_Z = 2              !Integer constant to express that this operation is in the z-direction
  integer , parameter :: DATA_SPEC_COLLISION       = 1
  integer , parameter :: DATA_SPEC_THERMAL         = 2
  integer , parameter :: DATA_SPEC_GRAVITY_WAVES   = 3
  integer , parameter :: DATA_SPEC_DENSITY_CURRENT = 5
  integer , parameter :: DATA_SPEC_INJECTION       = 6
  !Gauss-Legendre quadrature points and weights on the domain [0:1]
  integer , parameter :: nqpoints = 3
  real(rp), parameter :: qpoints (nqpoints) = (/ 0.1127E0_rp , 0.5000E0_rp , 0.8872E0_rp /)
  real(rp), parameter :: qweights(nqpoints) = (/ 0.2777E0_rp , 0.4444E0_rp , 0.2777E0_rp /)
  !! BEGIN USER-CONFIGURABLE PARAMETERS
  !The x-direction length is twice as long as the z-direction length
  integer , parameter :: nx = 800              !Number of total grid cells in the x dimension
  integer , parameter :: nz = 400              !Number of total grid cells in the z dimension
  real(rp), parameter :: sim_time = 10       !How many seconds to run the simulation
  real(rp), parameter :: output_freq = 10    !How frequently to output data to file (in seconds)
  integer , parameter :: data_spec_int = 5  !How to initialize the data (from 1 to 6)
  !! END USER-CONFIGURABLE PARAMETERS

  !! Variables that are initialized but remain static over the course of the simulation
  real(rp) :: dt                                  !Model time step (seconds)
  real(rp) :: dx, dz                              !Grid space length in x- and z-dimension (meters)
  real(rp), allocatable :: hy_dens_cell      (:)  !hydrostatic density (vert cell avgs).   Dimensions: (1-hs:nz+hs)
  real(rp), allocatable :: hy_dens_theta_cell(:)  !hydrostatic rho*t (vert cell avgs).     Dimensions: (1-hs:nz+hs)
  real(rp), allocatable :: hy_dens_int       (:)  !hydrostatic density (vert cell interf). Dimensions: (1:nz+1)
  real(rp), allocatable :: hy_dens_theta_int (:)  !hydrostatic rho*t (vert cell interf).   Dimensions: (1:nz+1)
  real(rp), allocatable :: hy_pressure_int   (:)  !hydrostatic press (vert cell interf).   Dimensions: (1:nz+1)
  !! Variables that are dynamic over the course of the simulation
  real(rp) :: etime                             !Elapsed model time
  real(rp) :: output_counter                    !Helps determine when it's time to do output
  !Runtime variable arrays
  real(rp), allocatable :: state    (:,:,:)     !Fluid state.                            Dimensions: (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
  real(rp), allocatable :: state_tmp(:,:,:)     !Fluid state.                            Dimensions: (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
  real(rp), allocatable :: flux     (:,:,:)     !Cell interface fluxes.                  Dimensions: (nx+1,nz+1,NUM_VARS)
  real(rp), allocatable :: tend     (:,:,:)     !Fluid state tendencies.                 Dimensions: (nx,nz,NUM_VARS)
  integer(8) :: t1, t2, rate, outputn=0         !For CPU Timings
  integer :: filenum=0                          !For counting output
  real(rp) :: mass0, te0, mass ,te                  !Initial domain totals for mass and total energy


  !! Main program
  !Initialize MPI, allocate arrays, initialize the grid and the data
  call init(dt)
  !Initial reductions for mass, kinetic energy, and total energy
  call reductions(mass0,te0)
  !Output the initial state
  call output(state,nx,nz,filenum)
  filenum=filenum+1

  !! Temporal loop
  call system_clock(t1)
  do while (etime < sim_time)
    call perform_timestep(state,state_tmp,flux,tend,dt)
    !Update the elapsed time and output counter
    etime = etime + dt
    output_counter = output_counter + dt
    !If it's time for output, reset the t, and do output
    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
      call output(state,nx,nz,filenum)
      filenum=filenum+1
    endif
  enddo
  call system_clock(t2,rate)
  write(*,*) "CPU Time: ",dble(t2-t1)/dble(rate)
  !Final reductions for mass, kinetic energy, and total energy
  call reductions(mass,te)
  deallocate(state,state_tmp,flux,tend,hy_dens_cell,hy_dens_theta_cell,hy_dens_int,hy_dens_theta_int,hy_pressure_int)

contains


  !Performs a single dimensionally split time step using a simple low-storate three-stage Runge-Kutta time integrator
  !The dimensional splitting is a second-order-accurate alternating Strang splitting in which the
  !order of directions is alternated each time step.
  !The Runge-Kutta method used here is defined as follows:
  ! q*     = q[n] + dt/3 * rhs(q[n])
  ! q**    = q[n] + dt/2 * rhs(q*  )
  ! q[n+1] = q[n] + dt/1 * rhs(q** )
  subroutine perform_timestep(state,state_tmp,flux,tend,dt)
    implicit none
    real(rp), intent(inout) :: state    (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: state_tmp(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux     (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend     (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    logical, save :: direction_switch = .true.
    if (direction_switch) then
      !x-direction first
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_X , flux , tend )
      !z-direction second
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_Z , flux , tend )
    else
      !z-direction second
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_Z , flux , tend )
      !x-direction first
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_X , flux , tend )
    endif
    direction_switch = .not. direction_switch
  end subroutine perform_timestep


  !Perform a single semi-discretized step in time with the form:
  !state_out = state_init + dt * rhs(state_forcing)
  !Meaning the step starts from state_init, computes the rhs using state_forcing, and stores the result in state_out
  subroutine semi_discrete_step( state_init , state_forcing , state_out , dt , dir , flux , tend )
    implicit none
    real(rp), intent(in   ) :: state_init   (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(inout) :: state_forcing(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: state_out    (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux         (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend         (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    integer , intent(in   ) :: dir
    integer :: i,k,ll
    real(rp) :: x, z, wpert, dist, x0, z0, xrad, zrad, amp

    if     (dir == DIR_X) then
      !Set the halo values for this MPI task's fluid state in the x-direction
      call set_halo_values_x(state_forcing)
      !Compute the time tendencies for the fluid state in the x-direction
      call compute_tendencies_x(state_forcing,flux,tend,dt)
    elseif (dir == DIR_Z) then
      !Set the halo values for this MPI task's fluid state in the z-direction
      call set_halo_values_z(state_forcing)
      !Compute the time tendencies for the fluid state in the z-direction
      call compute_tendencies_z(state_forcing,flux,tend,dt)
    endif

    !Apply the tendencies to the fluid state
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          if (data_spec_int == DATA_SPEC_GRAVITY_WAVES) then
            x = (i-0.5_rp)*dx
            z = (k-0.5_rp)*dz
            ! The following requires "acc routine" in OpenACC and "declare target" in OpenMP offload
            ! Neither of these are particularly well supported by compilers, so I'm manually inlining
            ! wpert = sample_ellipse_cosine( x,z , 0.01_rp , xlen/8,1000._rp, 500._rp,500._rp )
            x0 = xlen/8
            z0 = 1000
            xrad = 500
            zrad = 500
            amp = 0.01_rp
            !Compute distance from bubble center
            dist = sqrt( ((x-x0)/xrad)**2 + ((z-z0)/zrad)**2 ) * pi / 2._rp
            !If the distance from bubble center is less than the radius, create a cos**2 profile
            if (dist <= pi / 2._rp) then
              wpert = amp * cos(dist)**2
            else
              wpert = 0._rp
            endif
            tend(i,k,ID_WMOM) = tend(i,k,ID_WMOM) + wpert*hy_dens_cell(k)
          endif
          state_out(i,k,ll) = state_init(i,k,ll) + dt * tend(i,k,ll)
        enddo
      enddo
    enddo
  end subroutine semi_discrete_step














!Compute the time tendencies of the fluid state using forcing in the x-direction
subroutine compute_tendencies_x(state,flux,tend,dt)
    implicit none
    real(rp), intent(in   ) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dx / (16*dt)
    !Compute fluxes in the x-direction for each cell
    do k = 1 , nz
      do i = 1 , nx+1
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i-hs-1+s,k,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state (for artificial viscosity)
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo

        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_cell(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_cell(k) ) / r
        p = C0*(r*t)**gamma

        !Compute the flux vector
        flux(i,k,ID_DENS) = r*u     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*u*u+p - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*u*w   - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*u*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo
    !Use the fluxes to compute tendencies for each cell
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i+1,k,ll) - flux(i,k,ll) ) / dx
        enddo
      enddo
    enddo
  end subroutine compute_tendencies_x









  
!Compute the time tendencies of the fluid state using forcing in the z-direction
subroutine compute_tendencies_z(state,flux,tend,dt)
    implicit none
    real(rp), intent(in   ) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dz / (16*dt)
    !Compute fluxes in the x-direction for each cell
    do k = 1 , nz+1
      do i = 1 , nx
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i,k-hs-1+s,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo

        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_int(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_int(k) ) / r
        p = C0*(r*t)**gamma - hy_pressure_int(k)
        !Enforce vertical boundary condition and exact mass conservation
        if (k == 1 .or. k == nz+1) then
          w                = 0
          d3_vals(ID_DENS) = 0
        endif

        !Compute the flux vector with hyperviscosity
        flux(i,k,ID_DENS) = r*w     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*w*u   - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*w*w+p - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*w*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo

    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i,k+1,ll) - flux(i,k,ll) ) / dz
          if (ll == ID_WMOM) then
            tend(i,k,ID_WMOM) = tend(i,k,ID_WMOM) - state(i,k,ID_DENS)*grav
          endif
        enddo
      enddo
    enddo
  end subroutine


  subroutine set_halo_values_x(state)
    implicit none
    real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    integer :: k, ll
    real(rp) :: z
    if (data_spec_int == DATA_SPEC_INJECTION) then
      do k = 1 , nz
        z = (-1 + k-0.5_rp)*dz
        if (abs(z-3*zlen/4) <= zlen/16) then
          state(-1:0,k,ID_UMOM) = (state(-1:0,k,ID_DENS)+hy_dens_cell(k)) * 50._rp
          state(-1:0,k,ID_RHOT) = (state(-1:0,k,ID_DENS)+hy_dens_cell(k)) * 298._rp - hy_dens_theta_cell(k)
        endif
      enddo
    endif
  end subroutine set_halo_values_x



  subroutine set_halo_values_z(state)
    implicit none
    real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    integer :: i, ll
    real(rp), parameter :: mnt_width = xlen/8
    real(rp) :: x, xloc, mnt_deriv
    do ll = 1 , NUM_VARS
      do i = 1-hs,nx+hs
        if (ll == ID_WMOM) then
          state(i,-1  ,ll) = 0
          state(i,0   ,ll) = 0
          state(i,nz+1,ll) = 0
          state(i,nz+2,ll) = 0
        else if (ll == ID_UMOM) then
          state(i,-1  ,ll) = state(i,1 ,ll) / hy_dens_cell( 1) * hy_dens_cell(-1  )
          state(i,0   ,ll) = state(i,1 ,ll) / hy_dens_cell( 1) * hy_dens_cell( 0  )
          state(i,nz+1,ll) = state(i,nz,ll) / hy_dens_cell(nz) * hy_dens_cell(nz+1)
          state(i,nz+2,ll) = state(i,nz,ll) / hy_dens_cell(nz) * hy_dens_cell(nz+2)
        else
          state(i,-1  ,ll) = state(i,1 ,ll)
          state(i,0   ,ll) = state(i,1 ,ll)
          state(i,nz+1,ll) = state(i,nz,ll)
          state(i,nz+2,ll) = state(i,nz,ll)
        endif
      enddo
    enddo
  end subroutine set_halo_values_z


!Initialize some grid settings, allocate data, and initialize the full grid and data
subroutine init(dt)
  implicit none
  real(rp), intent(  out) :: dt
  integer :: i, k, ii, kk, ll, ierr
  real(rp) :: x, z, r, u, w, t, hr, ht
  !Set the cell grid size
  dx = xlen / nx
  dz = zlen / nz
  !Allocate the model data
  allocate(state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS))
  allocate(state_tmp(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS))
  allocate(flux              (nx+1,nz+1,NUM_VARS))
  allocate(tend              (nx,nz,NUM_VARS))
  allocate(hy_dens_cell      (1-hs:nz+hs))
  allocate(hy_dens_theta_cell(1-hs:nz+hs))
  allocate(hy_dens_int       (nz+1))
  allocate(hy_dens_theta_int (nz+1))
  allocate(hy_pressure_int   (nz+1))
  !Define the maximum stable time step based on an assumed maximum wind speed
  dt = min(dx,dz) / max_speed * cfl
  !Set initial elapsed model time and output_counter to zero
  etime = 0
  output_counter = 0
  print*,'nx, nz:', nx, nz
  print*,'dx,dz: ',dx,dz
  print*,'dt: ',dt

  !! Initialize the cell-averaged fluid state via Gauss-Legendre quadrature
  state = 0
  do k = 1-hs , nz+hs
    do i = 1-hs , nx+hs
      !Initialize the state to zero
      do kk = 1 , nqpoints
        do ii = 1 , nqpoints
          !Compute the x,z location within the global domain based on cell and quadrature index
          x = (i-0.5_rp) * dx + (qpoints(ii)-0.5_rp)*dx
          z = (k-0.5_rp) * dz + (qpoints(kk)-0.5_rp)*dz
          !Set the fluid state based on the user's specification
          if (data_spec_int == DATA_SPEC_COLLISION) call collision(x,z,r,u,w,t,hr,ht)
          if (data_spec_int == DATA_SPEC_THERMAL) call thermal(x,z,r,u,w,t,hr,ht)
          if (data_spec_int == DATA_SPEC_GRAVITY_WAVES) call gravity_waves(x,z,r,u,w,t,hr,ht)
          if (data_spec_int == DATA_SPEC_DENSITY_CURRENT) call density_current(x,z,r,u,w,t,hr,ht)
          if (data_spec_int == DATA_SPEC_INJECTION) call injection(x,z,r,u,w,t,hr,ht)
          !Store into the fluid state array
          state(i,k,ID_DENS) = state(i,k,ID_DENS) + r* qweights(ii)*qweights(kk)
          state(i,k,ID_UMOM) = state(i,k,ID_UMOM) + (r+hr)*u*qweights(ii)*qweights(kk)
          state(i,k,ID_WMOM) = state(i,k,ID_WMOM) + (r+hr)*w*qweights(ii)*qweights(kk)
          state(i,k,ID_RHOT) = state(i,k,ID_RHOT) + ( (r+hr)*(t+ht)-hr*ht)*qweights(ii)*qweights(kk)
        enddo
      enddo
      do ll = 1 , NUM_VARS
          state_tmp(i,k,ll) = state(i,k,ll)
        enddo
      enddo
    enddo
    !Compute the hydrostatic background state over vertical cell averages
    hy_dens_cell = 0
    hy_dens_theta_cell = 0
    do k = 1-hs , nz+hs
      do kk = 1 , nqpoints
        z = (k-0.5_rp) * dz + (qpoints(kk)-0.5_rp)*dz
        !Set the fluid state based on the user's specification
        if (data_spec_int == DATA_SPEC_COLLISION      ) call collision      (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA_SPEC_THERMAL        ) call thermal        (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA_SPEC_GRAVITY_WAVES  ) call gravity_waves  (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA_SPEC_DENSITY_CURRENT) call density_current(0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA_SPEC_INJECTION      ) call injection      (0._rp,z,r,u,w,t,hr,ht)
        hy_dens_cell(k)       = hy_dens_cell(k)       + hr    * qweights(kk)
        hy_dens_theta_cell(k) = hy_dens_theta_cell(k) + hr*ht * qweights(kk)
      enddo
    enddo
    !Compute the hydrostatic background state at vertical cell interfaces
    do k = 1 , nz+1
      z = (-1 + k-1) * dz
      if (data_spec_int == DATA_SPEC_COLLISION) call collision(0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA_SPEC_THERMAL) call thermal(0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA_SPEC_GRAVITY_WAVES) call gravity_waves(0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA_SPEC_DENSITY_CURRENT) call density_current(0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA_SPEC_INJECTION) call injection(0._rp,z,r,u,w,t,hr,ht)
      hy_dens_int(k) = hr
      hy_dens_theta_int(k) = hr*ht
      hy_pressure_int(k) = C0*(hr*ht)**gamma
    enddo
end subroutine init




  !This test case is initially balanced but injects fast, cold air from the left boundary near the model top
  subroutine injection(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
  end subroutine injection


  subroutine density_current(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z,-20._rp ,xlen/2,5000._rp,4000._rp,2000._rp)
  end subroutine density_current


  subroutine gravity_waves(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_bvfreq(z,0.02_rp,hr,ht)
    r = 0
    t = 0
    u = 15
    w = 0
  end subroutine gravity_waves


  !Rising thermal
  subroutine thermal(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z, 3._rp ,xlen/2,2000._rp,2000._rp,2000._rp)
  end subroutine thermal


  !Colliding thermals
  subroutine collision(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z, 20._rp,xlen/2,2000._rp,2000._rp,2000._rp)
    t = t + sample_ellipse_cosine(x,z,-20._rp,xlen/2,8000._rp,2000._rp,2000._rp)
  end subroutine collision


  !Establish hydrstatic balance using constant potential temperature (thermally neutral atmosphere)
  subroutine hydro_const_theta(z,r,t)
    implicit none
    real(rp), intent(in   ) :: z  !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, t  !Density and potential temperature at this point
    real(rp), parameter :: theta0 = 300._rp  !Background potential temperature
    real(rp), parameter :: exner0 = 1._rp    !Surface-level Exner pressure
    real(rp) :: p,exner,rt
    !Establish hydrostatic balance first using Exner pressure
    t = theta0                                  !Potential Temperature at z
    exner = exner0 - grav * z / (cp * theta0)   !Exner pressure at z
    p = p0 * exner**(cp/rd)                     !Pressure at z
    rt = (p / c0)**(1._rp / gamma)              !rho*theta at z
    r = rt / t                                  !Density at z
  end subroutine hydro_const_theta


  !Establish hydrstatic balance using constant Brunt-Vaisala frequency
  subroutine hydro_const_bvfreq( z , bv_freq0 , r , t )
    implicit none
    real(rp), intent(in   ) :: z , bv_freq0
    real(rp), intent(  out) :: r , t
    real(rp) :: theta0 = 300, exner0 = 1
    real(rp) :: p, exner, rt
    t = theta0 * exp( bv_freq0**2 / grav * z )                                  !Pot temp at z
    exner = exner0 - grav**2 / (cp * bv_freq0**2) * (t - theta0) / (t * theta0) !Exner pressure at z
    p = p0 * exner**(cp/rd)                                                     !Pressure at z
    rt = (p / c0)**(1._rp / gamma)                                              !rho*theta at z
    r = rt / t                                                                  !Density at z
  end subroutine hydro_const_bvfreq


  !Sample from an ellipse of a specified center, radius, and amplitude at a specified location
  function sample_ellipse_cosine( x , z , amp , x0 , z0 , xrad , zrad )   result(val)
    implicit none
    real(rp), intent(in) :: x , z         !Location to sample
    real(rp), intent(in) :: amp           !Amplitude of the bubble
    real(rp), intent(in) :: x0 , z0       !Center of the bubble
    real(rp), intent(in) :: xrad , zrad   !Radius of the bubble
    real(rp)             :: val           !The output sampled value
    real(rp) :: dist
    !Compute distance from bubble center
    dist = sqrt( ((x-x0)/xrad)**2 + ((z-z0)/zrad)**2 ) * pi / 2._rp
    !If the distance from bubble center is less than the radius, create a cos**2 profile
    if (dist <= pi / 2._rp) then
      val = amp * cos(dist)**2
    else
      val = 0._rp
    endif
  end function sample_ellipse_cosine












subroutine output(state,nx,nz,filenum)
  implicit none
  real(rp), intent(in) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
  character(len=40) :: namefile
  character(len=10) :: lf
  integer :: nx, ny=1, nz, nfields=4, filenum
  integer :: i,j,k
  character(len=200) :: filename
  character(len=80) :: buffer
  character(len=8) :: str1,str2,str3
  character(len=16) :: str4
  real(rp) :: dens(nx,nz),uwnd(nx,nz),wwnd(nx,nz),theta(nx,nz)
  !compute variables
  do k=1,nz
    do i=1,nx
      dens (i,k) = state(i,k,ID_DENS)
      uwnd (i,k) = state(i,k,ID_UMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
      wwnd (i,k) = state(i,k,ID_WMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
      theta(i,k) = (state(i,k,ID_RHOT)+hy_dens_theta_cell(k))/(hy_dens_cell(k)+state(i,k,ID_DENS))-hy_dens_theta_cell(k)/hy_dens_cell(k)
    enddo
  enddo
  ! end of line character
  lf=achar(10)
  write(filename,'(a,i8.8,a)') './vtk/pview_',filenum,'.vtk'
  open(666,file=trim(filename),status='replace',form='unformatted',access='stream',convert='big_endian')
  ! start writing vtk file
  buffer='# vtk DataFile Version 3.0'//lf
  write(666) trim(buffer)
  buffer='Advection'//lf
  write(666) trim(buffer)
  buffer='BINARY'//lf
  write(666) trim(buffer)
  buffer='DATASET RECTILINEAR_GRID'//lf
  write(666) trim(buffer)
  !write grid
  write(str1(1:8),'(i8)') nx
  write(str2(1:8),'(i8)') ny
  write(str3(1:8),'(i8)') nz
  buffer='DIMENSIONS '//str1//str2//str3//lf
  write(666) trim(buffer)
  buffer='X_COORDINATES '//str1//'  float'//lf
  write(666) trim(buffer)
  do i=1,nx
    write(666) real(i)
  enddo
  buffer='Y_COORDINATES '//str2//'  float'//lf ;
  write(666) trim(buffer)
  write(666) real(0.d0)
  buffer='Z_COORDINATES '//str3//'  float'//lf ;
  write(666) trim(buffer)
  do k=1,nz
    write(666) real(k)
  enddo
  ! write content (data format)
  write(str4(1:16),'(i16)') nx*ny*nz
  buffer='POINT_DATA '//str4//lf
  write(666) trim(buffer)
  write(str1(1:8),'(i8)') nfields
  buffer='FIELD FieldData '//str1//lf
  write(666) trim(buffer)
  ! write state 1 (scalar)
  write(str4(1:16),'(i16)') nx*nz
  buffer='dens 1'//str4//' float'//lf
  write(666) trim(buffer)
  do k=1,nz
    do i=1,nx
      write(666) real(dens(i,k))
    enddo
  enddo
  ! write state 2 (scalar)
  write(str4(1:16),'(i16)') nx*nz
  buffer='uwnd 1'//str4//' float'//lf
  write(666) trim(buffer)
  do k=1,nz
    do i=1,nx
      write(666) real(uwnd(i,k))
    enddo
  enddo
  ! write state 3 (scalar)
  write(str4(1:16),'(i16)') nx*nz
  buffer='wwnd 1'//str4//' float'//lf
  write(666) trim(buffer)
  do k=1,nz
    do i=1,nx
      write(666) real(wwnd(i,k))
    enddo
  enddo
  ! write state 4 (scalar)
  write(str4(1:16),'(i16)') nx*nz
  buffer='theta 1'//str4//' float'//lf
  write(666) trim(buffer)
  do k=1,nz
    do i=1,nx
      write(666) real(theta(i,k))
    enddo
  enddo

buffer=lf
write(666) trim(buffer)
close(666,status='keep')
end subroutine output











!Compute reduced quantities for error checking 
subroutine reductions( mass , te )
  implicit none
  real(rp), intent(out) :: mass, te
  integer :: i, k, ierr
  real(rp) :: r,u,w,th,p,t,ke,ie
  real(rp) :: glob(2)
  mass = 0
  te   = 0
  do k = 1 , nz
    do i = 1 , nx
      r  =   state(i,k,ID_DENS) + hy_dens_cell(k)             ! Density
      u  =   state(i,k,ID_UMOM) / r                           ! U-wind
      w  =   state(i,k,ID_WMOM) / r                           ! W-wind
      th = ( state(i,k,ID_RHOT) + hy_dens_theta_cell(k) ) / r ! Potential Temperature (theta)
      p  = C0*(r*th)**gamma      ! Pressure
      t  = th / (p0/p)**(rd/cp)  ! Temperature
      ke = r*(u*u+w*w)           ! Kinetic Energy
      ie = r*cv*t                ! Internal Energy
      mass = mass + r            *dx*dz ! Accumulate domain mass
      te   = te   + (ke + r*cv*t)*dx*dz ! Accumulate domain total energy
    enddo
  enddo
  mass = glob(1)
  te   = glob(2)
end subroutine


end program
