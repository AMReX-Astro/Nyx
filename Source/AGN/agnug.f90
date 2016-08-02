! Program to set up gas on a uniform grid and apply an AGN accretion and
! feedback model to it.
!
! Paul Ricker (2/2013)

program agnug

!==============================================================================

use agn_models

implicit none

integer                             :: N, i, j, k, m, mmax, accmodel, fbmodel
integer                             :: fnum
logical                             :: just_wrote
real                                :: gamma, t, dt, tmax, rhoinit, Pinit, tout
type(agn_particle)                  :: agn
type(agn_mesh_data)                 :: mesh
real                                :: h, L, Mdot, Eint, tnext
real, allocatable, dimension(:,:,:) :: rho, P, vx, vy, vz, E, cs, rhoold
real, allocatable, dimension(:,:,:) :: rhodot, rhoEdot, pxdot, pydot, pzdot

!==============================================================================

! Get input from user

write(*,*) 'agnug.'
write(*,*)
write(*,*) 'Enter size of grid in CGS units: '
read(*,*)  L
write(*,*) 'Enter number of zones per side of grid: '
read(*,*)  N
write(*,*) 'Enter initial density and pressure in CGS units: '
read(*,*)  rhoinit, Pinit
write(*,*) 'Enter adiabatic index: '
read(*,*)  gamma
write(*,*) 'Enter black hole mass in CGS units: '
read(*,*)  agn%m
write(*,*) 'Enter end time in seconds and maximum number of steps: '
read(*,*)  tmax, mmax
write(*,*) 'Enter choice of accretion and feedback models: '
read(*,*)  accmodel, fbmodel
write(*,*) 'Enter output time interval in seconds: '
read(*,*)  tout

!------------------------------------------------------------------------------

! Set up initial conditions

allocate(rho(N,N,N))
allocate(rhoold(N,N,N))
allocate(P(N,N,N))
allocate(vx(N,N,N))
allocate(vy(N,N,N))
allocate(vz(N,N,N))
allocate(E(N,N,N))
allocate(cs(N,N,N))

allocate(rhodot(N,N,N))
allocate(pxdot(N,N,N))
allocate(pydot(N,N,N))
allocate(pzdot(N,N,N))
allocate(rhoEdot(N,N,N))

allocate(mesh%x(N))
allocate(mesh%y(N))
allocate(mesh%z(N))
allocate(mesh%dx(N))
allocate(mesh%dy(N))
allocate(mesh%dz(N))

agn%mlast = agn%m ! black hole mass at time of last outburst

agn%x = 0.5 * L   ! position of black hole
agn%y = 0.5 * L
agn%z = 0.5 * L

agn%vx = 0.       ! velocity of black hole
agn%vy = 0.
agn%vz = 0.

agn%nx = 0.       ! direction cosines for jet angle
agn%ny = 0.
agn%nz = 1.

h = L / N

mesh%Nx   = N
mesh%xmin = 0.
mesh%xmax = L
mesh%x    = (/ ((i-0.5)*h, i = 1, N) /)
mesh%dx   = (/ (h, i = 1, N) /)

mesh%Ny   = N
mesh%ymin = 0.
mesh%ymax = L
mesh%y    = (/ ((i-0.5)*h, i = 1, N) /)
mesh%dy   = (/ (h, i = 1, N) /)

mesh%Nz   = N
mesh%zmin = 0.
mesh%zmax = L
mesh%z    = (/ ((i-0.5)*h, i = 1, N) /)
mesh%dz   = (/ (h, i = 1, N) /)

rho = rhoinit
P   = Pinit
vx  = 0.
vy  = 0.
vz  = 0.
E   = P/((gamma-1)*rho) + 0.5*(vx**2 + vy**2 + vz**2)

t   = 0.

call output(0, mesh, rho, P, vx, vy, vz, E)
just_wrote = .true.

fnum  = 1
tnext = tout

!------------------------------------------------------------------------------

! Timestep loop

write(*,*)
write(*,'(A1,5X,A1,13X,A2,12X,A3,11X,A4)') 'n', 't', 'dt', 'Mbh', 'Mdot'

do m = 1, mmax

!   Compute timestep

  cs = sqrt(gamma*P/rho)
  dt = 0.5 * h / maxval(cs + sqrt(vx**2 + vy**2 + vz**2))
  if (t+dt > tnext) dt = (tnext-t) * 1.000001

!   Determine accretion rate onto BH

  call agn_compute_acc_rate(agn, Mdot, t, mesh, rho, P, vx, vy, vz, E, gamma, &
                            accmodel)

!   Write log message

  write(*,'(I5,4(1X,ES13.6))') m, t, dt, agn%m, Mdot

!   Apply feedback to the gas

  call agn_compute_feedback(agn, Mdot, t, dt, mesh, rho, P, vx, vy, vz, E, &
                            rhodot, pxdot, pydot, pzdot, rhoEdot, fbmodel)

  rhoold = rho

  rho = max(rho + rhodot*dt, agn_small_rho)
  vx  = (rhoold*vx + pxdot*dt) / rho
  vy  = (rhoold*vy + pydot*dt) / rho
  vz  = (rhoold*vz + pzdot*dt) / rho
  E   = max( (rhoold*E + rhoEdot*dt) / rho, agn_small_e )

!   Update pressure by differencing total and kinetic energy; if we were
!   tracking temperature we would also call EOS here

  do k = 1, N
    do j = 1, N
      do i = 1, N
        Eint = E(i,j,k) - 0.5*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2)
        if (Eint > agn_small_e) then
          P(i,j,k) = (gamma-1) * rho(i,j,k) * Eint
        else
          P(i,j,k) = (gamma-1) * rho(i,j,k) * agn_small_e
          E(i,j,k) = agn_small_e + 0.5*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2)
        endif
      enddo
    enddo
  enddo

!   Grow the BH

  agn%m = agn%m + Mdot*dt

!   Update counters and output if necessary

  t = t + dt
  just_wrote = .false.

  if (t >= tnext) then
    call output(fnum, mesh, rho, P, vx, vy, vz, E)
    fnum = fnum + 1
    tnext = tnext + tout
    just_wrote = .true.
  endif

  if (t >= tmax) exit

enddo

if (.not. just_wrote) call output(fnum, mesh, rho, P, vx, vy, vz, E)

!==============================================================================

end program agnug




! output: A simple output routine that writes bricks of bytes.

subroutine output(n, m, rho, P, vx, vy, vz, E)

use agn_mesh

implicit none

integer, intent(in)             :: n
type(agn_mesh_data), intent(in) :: m
real, dimension(:,:,:)          :: rho, P, vx, vy, vz, E

character(len=8)                :: fname
integer                         :: nbytes

!==============================================================================

nbytes = m%Nx * m%Ny * m%Nz * 8

write(fname(1:8),'(A4,I4.4)') 'dens', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) rho
close(1)

write(fname(1:8),'(A4,I4.4)') 'pres', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) P
close(1)

write(fname(1:8),'(A4,I4.4)') 'velx', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) vx
close(1)

write(fname(1:8),'(A4,I4.4)') 'vely', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) vy
close(1)

write(fname(1:8),'(A4,I4.4)') 'velz', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) vz
close(1)

write(fname(1:8),'(A4,I4.4)') 'ener', n
open(1, file=fname, form='unformatted', access='direct', recl=nbytes)
write(1,rec=1) E
close(1)

!==============================================================================

return
end subroutine output
