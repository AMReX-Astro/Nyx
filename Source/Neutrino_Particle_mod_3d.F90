module neutrino_particle_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

  private

  public :: neutrino_deposit_relativistic_cic, neutrino_deposit_particle_dx_relativistic_cic

contains

  subroutine neutrino_deposit_relativistic_cic(particles, ns, np, nc, rho, lo, hi, plo, dx, csq) &
       bind(c,name='neutrino_deposit_relativistic_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    real(amrex_real)              :: csq

    integer          :: i, j, k, n, comp
    real(amrex_real) :: wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) :: lx, ly, lz
    real(amrex_real) :: inv_dx(3)
    real(amrex_real) :: vel_x, vel_y, vel_z, vsq, gamma

    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       vel_x = particles(5, n)
       vel_y = particles(6, n)
       vel_z = particles(7, n)
       vsq = vel_x*vel_x + vel_y*vel_y + vel_z*vel_z
       gamma = 1.d0 / dsqrt ( 1.d0 - vsq/csq)

       rho(i-1, j-1, k-1, 1) = rho(i-1, j-1, k-1, 1) + wx_lo*wy_lo*wz_lo*particles(4, n)*gamma
       rho(i-1, j-1, k  , 1) = rho(i-1, j-1, k  , 1) + wx_lo*wy_lo*wz_hi*particles(4, n)*gamma
       rho(i-1, j,   k-1, 1) = rho(i-1, j,   k-1, 1) + wx_lo*wy_hi*wz_lo*particles(4, n)*gamma
       rho(i-1, j,   k  , 1) = rho(i-1, j,   k,   1) + wx_lo*wy_hi*wz_hi*particles(4, n)*gamma
       rho(i,   j-1, k-1, 1) = rho(i,   j-1, k-1, 1) + wx_hi*wy_lo*wz_lo*particles(4, n)*gamma
       rho(i,   j-1, k  , 1) = rho(i,   j-1, k  , 1) + wx_hi*wy_lo*wz_hi*particles(4, n)*gamma
       rho(i,   j,   k-1, 1) = rho(i,   j,   k-1, 1) + wx_hi*wy_hi*wz_lo*particles(4, n)*gamma
       rho(i,   j,   k  , 1) = rho(i,   j,   k  , 1) + wx_hi*wy_hi*wz_hi*particles(4, n)*gamma

       rho(i-1, j-1, k-1, 2) = rho(i-1, j-1, k-1, comp) + wx_lo*wy_lo*wz_lo*particles(4, n)*vel_x
       rho(i-1, j-1, k  , 2) = rho(i-1, j-1, k  , comp) + wx_lo*wy_lo*wz_hi*particles(4, n)*vel_x
       rho(i-1, j,   k-1, 2) = rho(i-1, j,   k-1, comp) + wx_lo*wy_hi*wz_lo*particles(4, n)*vel_x
       rho(i-1, j,   k  , 2) = rho(i-1, j,   k,   comp) + wx_lo*wy_hi*wz_hi*particles(4, n)*vel_x
       rho(i,   j-1, k-1, 2) = rho(i,   j-1, k-1, comp) + wx_hi*wy_lo*wz_lo*particles(4, n)*vel_x
       rho(i,   j-1, k  , 2) = rho(i,   j-1, k  , comp) + wx_hi*wy_lo*wz_hi*particles(4, n)*vel_x
       rho(i,   j,   k-1, 2) = rho(i,   j,   k-1, comp) + wx_hi*wy_hi*wz_lo*particles(4, n)*vel_x
       rho(i,   j,   k  , 2) = rho(i,   j,   k  , comp) + wx_hi*wy_hi*wz_hi*particles(4, n)*vel_x

       rho(i-1, j-1, k-1, 3) = rho(i-1, j-1, k-1, comp) + wx_lo*wy_lo*wz_lo*particles(4, n)*vel_y
       rho(i-1, j-1, k  , 3) = rho(i-1, j-1, k  , comp) + wx_lo*wy_lo*wz_hi*particles(4, n)*vel_y
       rho(i-1, j,   k-1, 3) = rho(i-1, j,   k-1, comp) + wx_lo*wy_hi*wz_lo*particles(4, n)*vel_y
       rho(i-1, j,   k  , 3) = rho(i-1, j,   k,   comp) + wx_lo*wy_hi*wz_hi*particles(4, n)*vel_y
       rho(i,   j-1, k-1, 3) = rho(i,   j-1, k-1, comp) + wx_hi*wy_lo*wz_lo*particles(4, n)*vel_y
       rho(i,   j-1, k  , 3) = rho(i,   j-1, k  , comp) + wx_hi*wy_lo*wz_hi*particles(4, n)*vel_y
       rho(i,   j,   k-1, 3) = rho(i,   j,   k-1, comp) + wx_hi*wy_hi*wz_lo*particles(4, n)*vel_y
       rho(i,   j,   k  , 3) = rho(i,   j,   k  , comp) + wx_hi*wy_hi*wz_hi*particles(4, n)*vel_y

       rho(i-1, j-1, k-1, 4) = rho(i-1, j-1, k-1, comp) + wx_lo*wy_lo*wz_lo*particles(4, n)*vel_z
       rho(i-1, j-1, k  , 4) = rho(i-1, j-1, k  , comp) + wx_lo*wy_lo*wz_hi*particles(4, n)*vel_z
       rho(i-1, j,   k-1, 4) = rho(i-1, j,   k-1, comp) + wx_lo*wy_hi*wz_lo*particles(4, n)*vel_z
       rho(i-1, j,   k  , 4) = rho(i-1, j,   k,   comp) + wx_lo*wy_hi*wz_hi*particles(4, n)*vel_z
       rho(i,   j-1, k-1, 4) = rho(i,   j-1, k-1, comp) + wx_hi*wy_lo*wz_lo*particles(4, n)*vel_z
       rho(i,   j-1, k  , 4) = rho(i,   j-1, k  , comp) + wx_hi*wy_lo*wz_hi*particles(4, n)*vel_z
       rho(i,   j,   k-1, 4) = rho(i,   j,   k-1, comp) + wx_hi*wy_hi*wz_lo*particles(4, n)*vel_z
       rho(i,   j,   k  , 4) = rho(i,   j,   k  , comp) + wx_hi*wy_hi*wz_hi*particles(4, n)*vel_z

    end do

  end subroutine neutrino_deposit_relativistic_cic

  subroutine neutrino_deposit_particle_dx_relativistic_cic(particles, ns, np, nc, & 
                                                           rho, lo, hi, plo, dx, dx_particle) &
       bind(c,name='neutrino_deposit_particle_dx_relativistic_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    real(amrex_real)              :: dx_particle(3)

    integer i, j, k, n, comp
    real(amrex_real) lx, ly, lz, hx, hy, hz
    integer lo_x, lo_y, lo_z, hi_x, hi_y, hi_z
    real(amrex_real) wx, wy, wz
    real(amrex_real) inv_dx(3)
    real (amrex_real) factor, weight

    factor = (dx(1)/dx_particle(1))*(dx(2)/dx_particle(2))*(dx(3)/dx_particle(3))
    inv_dx = 1.0d0/dx

    do n = 1, np

       lx = (particles(1, n) - plo(1) - 0.5d0*dx_particle(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2) - 0.5d0*dx_particle(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3) - 0.5d0*dx_particle(3))*inv_dx(3)

       hx = (particles(1, n) - plo(1) + 0.5d0*dx_particle(1))*inv_dx(1)
       hy = (particles(2, n) - plo(2) + 0.5d0*dx_particle(2))*inv_dx(2)
       hz = (particles(3, n) - plo(3) + 0.5d0*dx_particle(3))*inv_dx(3)

       lo_x = floor(lx)
       lo_y = floor(ly)
       lo_z = floor(lz)

       hi_x = floor(hx)
       hi_y = floor(hy)
       hi_z = floor(hz)

       do i = lo_x, hi_x
          if (i < lo(1) .or. i > hi(1)) then
             cycle
          end if
          wx = min(hx - i, 1.d0) - max(lx - i, 0.d0)
          do j = lo_y, hi_y
             if (j < lo(2) .or. j > hi(2)) then
                cycle
             end if
             wy = min(hy - j, 1.d0) - max(ly - j, 0.d0)
             do k = lo_z, hi_z
                if (k < lo(3) .or. k > hi(3)) then
                   cycle
                end if
                wz = min(hz - k, 1.d0) - max(lz - k, 0.d0)

                weight = wx*wy*wz*factor

                rho(i, j, k, 1) = rho(i, j, k, 1) + weight*particles(4, n)

                do comp = 2, nc
                   rho(i, j, k, comp) = rho(i, j, k, comp) + weight*particles(4, n)*particles(3+comp, n) 
                end do
             end do
          end do
       end do
    end do

  end subroutine neutrino_deposit_particle_dx_relativistic_cic

end module neutrino_particle_module
