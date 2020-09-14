
      subroutine derforcex(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the x-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use turbforce_module
      use amrex_constants_module, only : TWO, HALF, ZERO, M_PI
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer :: i,j,k
      integer :: kx,ky,kz
      integer :: xstep,ystep,zstep
      real(rt) :: x,  y,  z
      real(rt) :: Lx, Ly, Lz, freqx, freqy, freqz
      real(rt) :: cosx,cosy,cosz,sinx,siny,sinz
      real(rt) :: HLx,HLy,HLz
      real(rt) :: kxd,kyd,kzd
      real(rt) :: kappa,kappaMax,Lmin,xT
      real(rt) :: f1,twicePi
      
      twicePi = TWO*M_PI

      Lx = prob_hi(1)-prob_lo(1)
      Ly = prob_hi(2)-prob_lo(2)
      Lz = prob_hi(3)-prob_lo(3)
      
      Lmin = min(Lx,Ly,Lz)
      kappaMax = dble(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(HALF+Lx/Lmin)
      nymodes = nmodes*int(HALF+Ly/Lmin)
      nzmodes = nmodes*int(HALF+Lz/Lmin)
      
      xstep = int(Lx/Lmin+HALF)
      ystep = int(Ly/Lmin+HALF)
      zstep = int(Lz/Lmin+HALF)
      
      HLx = Lx
      HLy = Ly
      HLz = Lz
      
      do k = lo(3),hi(3)
         z = (dble(k) + HALF) * delta(3)
         
         do j = lo(2),hi(2)
            y =  (dble(j) + HALF) * delta(2)

            do i = lo(1),hi(1)
               x = (dble(i) + HALF) * delta(1)
               
               f1 = ZERO
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dble(kz)
                  freqz = twicePi*kzd*HLz
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dble(ky)
                     freqy=twicePi*kyd/HLy
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dble(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        freqx = twicePi*kxd/HLx
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz))
                           f1 = f1 + xT * ( FAZ(kx,ky,kz)*freqy*sin(freqx*x+FPZX(kx,ky,kz)) * cos(freqy*y+FPZY(kx,ky,kz)) * &
                                                                sin(freqz*z+FPZZ(kx,ky,kz)) &
                                           -FAY(kx,ky,kz)*freqz*sin(freqx*x+FPYX(kx,ky,kz)) * sin(freqy*y+FPYY(kx,ky,kz)) * &
                                                                cos(freqz*z+FPYZ(kx,ky,kz)) )
                        endif
                     enddo
                  enddo
               enddo
               
               force(i,j,k,1) = dat(i,j,k,1)*f1

            enddo
         enddo
      enddo
      
      end subroutine derforcex

!-----------------------------------------------------------------------

      subroutine derforcey(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the y-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use turbforce_module
      use amrex_constants_module, only : TWO, HALF, ZERO, M_PI
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer :: i,j,k
      integer :: kx,ky,kz
      integer :: xstep,ystep,zstep
      real(rt) :: x,  y,  z
      real(rt) :: Lx, Ly, Lz, freqx, freqy, freqz
      real(rt) :: cosx,cosy,cosz,sinx,siny,sinz
      real(rt) :: HLx,HLy,HLz
      real(rt) :: kxd,kyd,kzd
      real(rt) :: kappa,kappaMax,Lmin,xT
      real(rt) :: f2,twicePi
      
      twicePi = TWO*M_PI

      Lx = prob_hi(1)-prob_lo(1)
      Ly = prob_hi(2)-prob_lo(2)
      Lz = prob_hi(3)-prob_lo(3)
      
      Lmin = min(Lx,Ly,Lz)
      kappaMax = dble(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(HALF+Lx/Lmin)
      nymodes = nmodes*int(HALF+Ly/Lmin)
      nzmodes = nmodes*int(HALF+Lz/Lmin)
      
      xstep = int(Lx/Lmin+HALF)
      ystep = int(Ly/Lmin+HALF)
      zstep = int(Lz/Lmin+HALF)
      
      HLx = Lx
      HLy = Ly
      HLz = Lz
      
      do k = lo(3),hi(3)
         z = (dble(k) + HALF) * delta(3)
         
         do j = lo(2),hi(2)
            y =  (dble(j) + HALF) * delta(2)

            do i = lo(1),hi(1)
               x = (dble(i) + HALF) * delta(1)
               
               f2 = ZERO
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dble(kz)
                  freqz = twicePi*kzd*HLz
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dble(ky)
                     freqy=twicePi*kyd/HLy
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dble(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        freqx = twicePi*kxd/HLx
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz))
                           f2 = f2 + xT * ( FAX(kx,ky,kz)*freqz*sin(freqx*x+FPXX(kx,ky,kz)) * sin(freqy*y+FPXY(kx,ky,kz)) * &
                                                                cos(freqz*z+FPXZ(kx,ky,kz)) &
                                           -FAZ(kx,ky,kz)*freqx*cos(freqx*x+FPZX(kx,ky,kz)) * sin(freqy*y+FPZY(kx,ky,kz)) * &
                                                                sin(freqz*z+FPZZ(kx,ky,kz)) )
                        endif
                     enddo
                  enddo
               enddo
               
               force(i,j,k,1) = dat(i,j,k,1)*f2

            enddo
         enddo
      enddo
      
      end subroutine derforcey


!-----------------------------------------------------------------------

      subroutine derforcez(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the z-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use turbforce_module
      use amrex_constants_module, only : TWO, HALF, ZERO, M_PI
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer :: i,j,k
      integer :: kx,ky,kz
      integer :: xstep,ystep,zstep
      real(rt) :: x,  y,  z
      real(rt) :: Lx, Ly, Lz, freqx, freqy, freqz
      real(rt) :: cosx,cosy,cosz,sinx,siny,sinz
      real(rt) :: HLx,HLy,HLz
      real(rt) :: kxd,kyd,kzd
      real(rt) :: kappa,kappaMax,Lmin,xT
      real(rt) :: f3,twicePi
      
      twicePi = TWO*M_PI

      Lx = prob_hi(1)-prob_lo(1)
      Ly = prob_hi(2)-prob_lo(2)
      Lz = prob_hi(3)-prob_lo(3)
      
      Lmin = min(Lx,Ly,Lz)
      kappaMax = dble(nmodes)/Lmin + 1.0d-8
      nxmodes = nmodes*int(HALF+Lx/Lmin)
      nymodes = nmodes*int(HALF+Ly/Lmin)
      nzmodes = nmodes*int(HALF+Lz/Lmin)
      
      xstep = int(Lx/Lmin+HALF)
      ystep = int(Ly/Lmin+HALF)
      zstep = int(Lz/Lmin+HALF)
      
      HLx = Lx
      HLy = Ly
      HLz = Lz

      do k = lo(3),hi(3)
         z = (dble(k) + HALF) * delta(3)
         
         do j = lo(2),hi(2)
            y =  (dble(j) + HALF) * delta(2)

            do i = lo(1),hi(1)
               x = (dble(i) + HALF) * delta(1)
               
               f3 = ZERO
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dble(kz)
                  freqz = twicePi*kzd*HLz
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dble(ky)
                     freqy=twicePi*kyd/HLy
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dble(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        freqx = twicePi*kxd/HLx
                        if (kappa.le.kappaMax) then
                           xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz))
                           f3 = f3 + xT * ( FAY(kx,ky,kz)*freqx*cos(freqx*x+FPYX(kx,ky,kz)) * sin(freqy*y+FPYY(kx,ky,kz)) * &
                                                                sin(freqz*z+FPYZ(kx,ky,kz)) &
                                           -FAX(kx,ky,kz)*freqy*sin(freqx*x+FPXX(kx,ky,kz)) * cos(freqy*y+FPXY(kx,ky,kz)) * &
                                                                sin(freqz*z+FPXZ(kx,ky,kz)) ) 
                        endif
                     enddo
                  enddo
               enddo
               
               force(i,j,k,1) = dat(i,j,k,1)*f3

            enddo
         enddo
      enddo
      
      end subroutine derforcez

