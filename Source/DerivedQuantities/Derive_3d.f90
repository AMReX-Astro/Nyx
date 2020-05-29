
      subroutine dersoundspeed(c,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UEINT
      use  eos_params_module
      implicit none

      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      real(rt) c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp_c)
      real(rt) u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      real(rt) dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      real(rt) :: e
      integer          :: i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom, (rho E), (rho e))
      ! 

      !
      ! Compute soundspeed from the EOS.
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               e      = u(i,j,k,UEINT) / u(i,j,k,URHO)

               if (e .gt. 0.d0) then
                  call nyx_eos_soundspeed(c(i,j,k,1), u(i,j,k,URHO), e)
               else
                  c(i,j,k,1) = 0.d0
               end if

            enddo
         enddo
      enddo

      end subroutine dersoundspeed

!-----------------------------------------------------------------------

      subroutine dermachnumber(mach,mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
      use  eos_params_module
      implicit none

      integer          :: mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      real(rt) :: mach(mach_l1:mach_h1,mach_l2:mach_h2,mach_l3:mach_h3,ncomp_mach)
      real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      real(rt) :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      real(rt) :: rhoInv,ux,uy,uz,e,c
      integer          :: i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom, (rho E), (rho e))
      ! 

      ! 
      ! Compute Mach number of the flow.
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoInv = 1.d0 / u(i,j,k,URHO)
               ux     = u(i,j,k,UMX)*rhoInv
               uy     = u(i,j,k,UMY)*rhoInv
               uz     = u(i,j,k,UMZ)*rhoInv
               e      = u(i,j,k,UEINT)*rhoInv

               if (e .gt. 0.d0) then
                  call nyx_eos_soundspeed(c, u(i,j,k,URHO), e)
                  mach(i,j,k,1) = sqrt(ux**2 + uy**2 + uz**2) / c
               else
                  mach(i,j,k,1) = 0.d0
               end if

            enddo
         enddo
      enddo

      end subroutine dermachnumber

!-----------------------------------------------------------------------

      subroutine derentropy(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s, &
                               u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi, &
                               domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)
      !
      ! Compute entropy from the EOS.
      !
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UEINT
      use  eos_params_module
      implicit none

      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      real(rt) s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ncomp_s)
      real(rt) u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      real(rt) dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      real(rt) :: e, rhoInv
      integer i,j,k

      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom, (rho E), (rho e), Temp, Ne)
      ! 
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rhoInv = 1.d0/u(i,j,k,URHO)
               e  = u(i,j,k,UEINT)*rhoInv

!              if (e .gt. 0.d0) then
!                 call nyx_eos_S_given_Re(s(i,j,k,1), u(i,j,k,URHO), e, &
!                                         u(i,j,k,7), u(i,j,k,8), &
!                                         comoving_a = 1.d0)
!              else
!                 s(i,j,k,1) = 0.d0
!              end if
            enddo
         enddo
      enddo

      end subroutine derentropy

!-----------------------------------------------------------------------

      subroutine derlogden(logden,ld_l1,ld_l2,ld_l3,ld_h1,ld_h2,ld_h3,nd, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer           ld_l1, ld_l2, ld_l3, ld_h1, ld_h2, ld_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) logden( ld_l1: ld_h1, ld_l2: ld_h2, ld_l3: ld_h3,nd)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
 
      integer    i,j,k
      ! 
      ! Here dat contains (Density)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               logden(i,j,k,1) = dlog10(dat(i,j,k,1))
            end do
         end do
      end do
 
      end subroutine derlogden

!-----------------------------------------------------------------------

      subroutine dermagvort(vort,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,nv, & 
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate vorticity
      !     
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer            v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) vort(  v_l1:  v_h1,  v_l2:  v_h2,  v_l3:  v_h3,nv)
      real(rt)  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

      integer          :: i,j,k
      real(rt) :: uy,uz,vx,vz,wx,wy,v1,v2,v3
      real(rt) :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)

      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom)
      ! 

      !
      ! Convert momentum to velocity.
      !
      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               ldat(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,1)
               ldat(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
               ldat(i,j,k,4) = dat(i,j,k,4) / dat(i,j,k,1)
            end do
         end do
      end do
      !
      ! Calculate vorticity.
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = (dat(i,j+1,k,2) - dat(i,j-1,k,2)) / delta(2)
               uz = (dat(i,j,k+1,2) - dat(i,j,k-1,2)) / delta(3)
               vx = (dat(i+1,j,k,3) - dat(i-1,j,k,3)) / delta(1)
               vz = (dat(i,j,k+1,3) - dat(i,j,k-1,3)) / delta(3)
               wx = (dat(i+1,j,k,4) - dat(i-1,j,k,4)) / delta(1)
               wy = (dat(i,j+1,k,4) - dat(i,j-1,k,4)) / delta(2)
               v1 = 0.5d0 * abs(wy - vz)
               v2 = 0.5d0 * abs(uz - wx)
               v3 = 0.5d0 * abs(vx - uy)
               vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)
            end do
         end do
      end do

      end subroutine dermagvort

!-----------------------------------------------------------------------

      subroutine derdivu(divu,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will divergence of velocity.
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) divu(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3,nd)
      real(rt)  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer          :: i,j,k
      real(rt) :: ulo,uhi,vlo,vhi,wlo,whi
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uhi = dat(i+1,j,k,2) / dat(i+1,j,k,1)
               ulo = dat(i-1,j,k,2) / dat(i-1,j,k,1)
               vhi = dat(i,j+1,k,3) / dat(i,j+1,k,1)
               vlo = dat(i,j-1,k,3) / dat(i,j-1,k,1)
               whi = dat(i,j,k+1,4) / dat(i,j,k+1,1)
               wlo = dat(i,j,k-1,4) / dat(i,j,k-1,1)
               divu(i,j,k,1) = 0.5d0 * ( (uhi-ulo) / delta(1) + &
                                         (vhi-vlo) / delta(2) + &
                                         (whi-wlo) / delta(3) )
            end do
         end do
      end do

      end subroutine derdivu

!-----------------------------------------------------------------------

      subroutine dermomt(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine computes Mom + Mom*Sdens/Density
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) vel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k

      ! 
      ! Here dat contains (Density, Single Component of Momentum, Sdens)
      ! 

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = dat(i,j,k,2) + dat(i,j,k,2)*dat(i,j,k,3)/dat(i,j,k,1)
            end do
         end do
      end do

      end subroutine dermomt
