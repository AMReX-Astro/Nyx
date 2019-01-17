! This module initializes the coefficients for the turbulent forcing

module turbinit_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_paralleldescriptor_module
  use turbforce_module

  implicit none

contains

  subroutine turbforce_init(prob_lo,prob_hi)

    use amrex_fort_module  , only: amrex_random
    use amrex_constants_module, only: TWO, ONE, HALF, ZERO, M_PI

    real(rt), intent(in) :: prob_lo(:), prob_hi(:)
    
    integer :: kx,ky,kz
    integer :: xstep,ystep,zstep
    integer :: mode_count

    real(rt) :: Lx, Ly, Lz
    real(rt) :: kxd,kyd,kzd
    real(rt) :: kappa,kappaMax,Lmin
    real(rt) :: twicePi
    real(rt) :: thetaTmp,cosThetaTmp,sinThetaTmp
    real(rt) :: phiTmp,cosPhiTmp,sinPhiTmp
    real(rt) :: px,py,pz,mp2,Ekh
    real(rt) :: freqMin,freqMax,freqDiff
    real(rt) :: rn

    if (amrex_pd_ioprocessor()) &
       write (*,*) "Initialising random number generator..."

    twicePi = TWO*M_PI

!   if (blrandseed.gt.0) then
!      call init_genrand(20908)
!      call blutilinitrand(blrandseed)
!      rn = genrand_real1()
!      call blutilinitrand(blrandseed)
!      if (amrex_pd_ioprocessor()) then
!         write (*,*) "blrandseed = ",blrandseed
!         write (*,*) "first random number = ",rn
!      endif
!   else
!      call init_genrand(111397)
!      rn = genrand_real1()
!   endif

    rn = amrex_random()

    if (amrex_pd_ioprocessor()) then
       print *,"first random number = ",rn
    endif

    Lx = prob_hi(1)-prob_lo(1)
    Ly = prob_hi(2)-prob_lo(2)
    Lz = prob_hi(3)-prob_lo(3)

    if (amrex_pd_ioprocessor()) then
       write(*,*) "Lx = ",Lx
       write(*,*) "Ly = ",Ly
       write(*,*) "Lz = ",Lz
    endif

         Lmin = min(Lx,Ly,Lz)
         kappaMax = dble(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5+Lx/Lmin)
         nymodes = nmodes*int(0.5+Ly/Lmin)
         nzmodes = nmodes*int(0.5+Lz/Lmin)
         if (amrex_pd_ioprocessor()) then
            write(*,*) "Lmin = ",Lmin
            write(*,*) "kappaMax = ",kappaMax
            write(*,*) "nxmodes = ",nxmodes
            write(*,*) "nymodes = ",nymodes
            write(*,*) "nzmodes = ",nzmodes
         endif

         if (forcing_time_scale_min.eq.ZERO) then
            forcing_time_scale_min = HALF
         endif
         if (forcing_time_scale_max.eq.ZERO) then
            forcing_time_scale_max = ONE
         endif

         freqMin = one/forcing_time_scale_max
         freqMax = one/forcing_time_scale_min
         freqDiff= freqMax-freqMin

         if (amrex_pd_ioprocessor()) then
            write(*,*) "forcing_time_scale_min = ",forcing_time_scale_min
            write(*,*) "forcing_time_scale_max = ",forcing_time_scale_max
            write(*,*) "freqMin = ",freqMin
            write(*,*) "freqMax = ",freqMax
            write(*,*) "freqDiff = ",freqDiff
         endif

         mode_count = 0

         xstep = int(Lx/Lmin+0.5)
         ystep = int(Ly/Lmin+0.5)
         zstep = int(Lz/Lmin+0.5)
         if (amrex_pd_ioprocessor()) then
            write (*,*) "Mode step ",xstep, ystep, zstep
         endif

         do kz = mode_start*zstep, nzmodes, zstep
            kzd = dble(kz)

            do ky = mode_start*ystep, nymodes, ystep
               kyd = dble(ky)

               do kx = mode_start*xstep, nxmodes, xstep
                  kxd = dble(kx)

                  kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

                  if (kappa.le.kappaMax) then
                     rn = amrex_random()
                     FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     rn = amrex_random()
                     FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     rn = amrex_random()
                     FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!     Translation angles, theta=0..2Pi and phi=0..Pi
                     rn = amrex_random()
                     TAT(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     TAP(kx,ky,kz) = rn*M_PI
!     Phases
                     rn = amrex_random()
                     FPXX(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPYX(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPZX(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPXY(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPYY(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPZY(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPXZ(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPYZ(kx,ky,kz) = rn*twicePi
                     rn = amrex_random()
                     FPZZ(kx,ky,kz) = rn*twicePi
!     Amplitudes (alpha)
                     rn = amrex_random()
                     thetaTmp      = rn*twicePi
                     cosThetaTmp   = cos(thetaTmp)
                     sinThetaTmp   = sin(thetaTmp)
                     rn = amrex_random()
                     phiTmp        = rn*M_PI
                     cosPhiTmp     = cos(phiTmp)
                     sinPhiTmp     = sin(phiTmp)

                     px = cosThetaTmp * sinPhiTmp
                     py = sinThetaTmp * sinPhiTmp
                     pz =               cosPhiTmp

                     mp2           = px*px + py*py + pz*pz
                     if (mp2 .lt. 0.000001) then
                        write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
!     Count modes that contribute
                        mode_count = mode_count + 1
!     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
                        endif
                        if (force_scale.ge.zero) then
                           FAX(kx,ky,kz) = force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2
                        else
                           FAX(kx,ky,kz) = px * Ekh / mp2
                           FAY(kx,ky,kz) = py * Ekh / mp2
                           FAZ(kx,ky,kz) = pz * Ekh / mp2
                        endif

!                       if (amrex_pd_ioprocessor()) then
!                          write (*,*) "Mode"
!                          write (*,*) "kappa = ",kx,ky,kz,kappa
!                          write (*,*) "Amplitudes"
!                          write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
!                          write (*,*) "Frequencies"
!                          write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
!                       endif
                     endif
                  endif
               enddo
            enddo
         enddo

         call flush(6)

  end subroutine turbforce_init

end module turbinit_module
