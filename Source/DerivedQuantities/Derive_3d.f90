
!-----------------------------------------------------------------------

      subroutine ca_derstate(state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv, &
                             dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
      ! The outgoing "state" vector contains (rho,T,X_1)
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_3d.f90 :: ca_derstate')
      end if
      !
      ! Density
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,1) = dat(i,j,k,1)
            end do
         end do
      end do
      !
      ! Temperature
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,2) = dat(i,j,k,2)
            end do
         end do
      end do
      !
      ! (rho X)_1 --> X_1
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_derstate

!-----------------------------------------------------------------------

      subroutine ca_dervel(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! Derive velocity from momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k
      ! 
      ! Here dat contains (Density, Single Component of Momentum)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_dervel

!-----------------------------------------------------------------------

      subroutine ca_dermagvel(magvel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! Derive magnitude of velocity.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision magvel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               magvel(i,j,k,1) = sqrt( (dat(i,j,k,2) / dat(i,j,k,1))**2 + &
                                       (dat(i,j,k,3) / dat(i,j,k,1))**2 + &
                                       (dat(i,j,k,4) / dat(i,j,k,1))**2 )
            end do
         end do
      end do

      end subroutine ca_dermagvel

!-----------------------------------------------------------------------

      subroutine ca_dermaggrav(maggrav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng, &
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! Derive magnitude of the gravity vector.
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision maggrav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3,ng)
      double precision     dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (grav_x, grav_y, grav_z)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               maggrav(i,j,k,1) = sqrt( dat(i,j,k,1)**2  + &
                                        dat(i,j,k,2)**2  + &
                                        dat(i,j,k,3)**2 )
            end do
         end do
      end do

      end subroutine ca_dermaggrav

!-----------------------------------------------------------------------

      subroutine ca_dermagmom(magmom,mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision magmom(mom_l1:mom_h1,mom_l2:mom_h2,mom_l3:mom_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (Xmom, Ymom, Zmom)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
            end do
         end do
      end do

      end subroutine ca_dermagmom

!-----------------------------------------------------------------------

      subroutine ca_derpres(p,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)
      !
      ! Compute pressure from (rho e)
      !
      use meth_params_module, only : UEINT, gamma_minus_1
      use  eos_params_module

      implicit none

      integer p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision p(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp_p)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      integer          :: i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom, (rho E), (rho e))
      ! 
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               !
               ! Protect against negative internal energy.
               !
               if (u(i,j,k,UEINT) .le. 0.d0) then
                  print *,'   '
                  print *,'>>> Error: deriving pressure at ',i,j,k
                  print *,'>>> but rho*eint is negative: ', u(i,j,k,UEINT)
                  print *,'    '
                  call bl_error("Error:: Derive_3d.f90 :: derpres")
               else
                  p(i,j,k,1) = gamma_minus_1 * u(i,j,k,UEINT)
               end if

            enddo
         enddo
      enddo

      end subroutine ca_derpres

!-----------------------------------------------------------------------

      subroutine ca_dereint1(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)
      !
      ! Compute internal energy from (rho E).
      !
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN 

      implicit none

      integer e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy,uz
      integer          :: i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom, (rho E), (rho e))
      ! 
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rhoInv = 1.d0/u(i,j,k,URHO)
               ux = u(i,j,k,UMX)*rhoInv
               uy = u(i,j,k,UMY)*rhoInv
               uz = u(i,j,k,UMZ)*rhoInv
               e(i,j,k,1) = u(i,j,k,UEDEN)*rhoInv-0.5d0*(ux**2+uy**2+uz**2)
            enddo
         enddo
      enddo

      end subroutine ca_dereint1

!-----------------------------------------------------------------------

      subroutine ca_dereint2(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UEINT

      implicit none

      integer e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      integer :: i,j,k
      !
      ! Compute internal energy from (rho e).
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               e(i,j,k,1) = u(i,j,k,UEINT) / u(i,j,k,URHO)
            enddo
         enddo
      enddo

      end subroutine ca_dereint2

!-----------------------------------------------------------------------

      subroutine ca_dersoundspeed(c,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use eos_module
      use meth_params_module, only : URHO, UEINT
      use  eos_params_module
      implicit none

      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp_c)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: e
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

      end subroutine ca_dersoundspeed

!-----------------------------------------------------------------------

      subroutine ca_dermachnumber(mach,mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
      use  eos_params_module
      implicit none

      integer          :: mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      double precision :: mach(mach_l1:mach_h1,mach_l2:mach_h2,mach_l3:mach_h3,ncomp_mach)
      double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy,uz,e,c
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

      end subroutine ca_dermachnumber

!-----------------------------------------------------------------------

      subroutine ca_derentropy(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s, &
                               u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi, &
                               domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)
      !
      ! Compute entropy from the EOS.
      !
      use eos_module
      use meth_params_module, only : URHO, UEINT
      use  eos_params_module
      implicit none

      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ncomp_s)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: e, rhoInv
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

      end subroutine ca_derentropy

!-----------------------------------------------------------------------

      subroutine ca_derspec(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive X_i from (rho X)_i
      !
      implicit none

      integer          lo(3), hi(3)
      integer          spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k
      ! 
      ! Here dat contains (Density, (rho X)_i)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_derspec

!-----------------------------------------------------------------------

      subroutine ca_derlogden(logden,ld_l1,ld_l2,ld_l3,ld_h1,ld_h2,ld_h3,nd, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      implicit none

      integer          lo(3), hi(3)
      integer           ld_l1, ld_l2, ld_l3, ld_h1, ld_h2, ld_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision logden( ld_l1: ld_h1, ld_l2: ld_h2, ld_l3: ld_h3,nd)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
 
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
 
      end subroutine ca_derlogden

!-----------------------------------------------------------------------

      subroutine ca_dermagvort(vort,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,nv, & 
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate vorticity
      !     
      implicit none

      integer          lo(3), hi(3)
      integer            v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vort(  v_l1:  v_h1,  v_l2:  v_h2,  v_l3:  v_h3,nv)
      double precision  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

      integer          :: i,j,k
      double precision :: uy,uz,vx,vz,wx,wy,v1,v2,v3
      double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)

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

      end subroutine ca_dermagvort

!-----------------------------------------------------------------------

      subroutine ca_derdivu(divu,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will divergence of velocity.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision divu(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3,nd)
      double precision  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer          :: i,j,k
      double precision :: ulo,uhi,vlo,vhi,wlo,whi
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

      end subroutine ca_derdivu

!-----------------------------------------------------------------------

      subroutine ca_derkineng(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2)
      !
      implicit none

      integer          lo(3), hi(3)
      integer          ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (Density, Xmom, Ymom, Zmom)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               kineng(i,j,k,1) = 0.5d0 / dat(i,j,k,1) * ( dat(i,j,k,2)**2 + &
                                                          dat(i,j,k,3)**2 + &
                                                          dat(i,j,k,4)**2 )
            end do
         end do
      end do

      end subroutine ca_derkineng

!-----------------------------------------------------------------------

      subroutine ca_dernull(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine is used by particle_count.  Yes it does nothing.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      end subroutine ca_dernull

!-----------------------------------------------------------------------

      subroutine ca_dermomt(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine computes Mom + Mom*Sdens/Density
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
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

      end subroutine ca_dermomt
