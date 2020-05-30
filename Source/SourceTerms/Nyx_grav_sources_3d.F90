
! :::
! ::: ------------------------------------------------------------------
! :::

      !===========================================================================
      ! This is called from the C++ so the threading happens here...
      !===========================================================================
      AMREX_CUDA_FORT_DEVICE subroutine fort_correct_gsrc(lo,hi, &
                              gold,gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3, &
                              gnew,gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3, &
                              uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                              unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                              a_old,a_new,dt) &
                              bind(C, name="fort_correct_gsrc")

#ifndef AMREX_USE_CUDA
      use amrex_error_module
#endif
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      real(rt)   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      real(rt)   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      real(rt)  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      real(rt)  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      real(rt), intent(in), value ::  a_old,a_new,dt

      integer i,j,k
      real(rt) SrU_old, SrV_old, SrW_old
      real(rt) SrU_new, SrV_new, SrW_new
      real(rt) SrUcorr, SrVcorr, SrWcorr, SrEcorr
      real(rt) rhoo, Upo, Vpo, Wpo
      real(rt) rhon, Upn, Vpn, Wpn

      real(rt) a_half, a_newsq, rhooinv, rhoninv, a_new_inv
      real(rt) old_ke, old_rhoeint
      real(rt) new_ke, new_rhoeint

      a_half    = 0.5d0 * (a_old + a_new)
      a_newsq   = a_new*a_new
      a_new_inv = 1.0d0 / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rhoo    = uold(i,j,k,URHO)
               rhooinv = 1.0d0 / uold(i,j,k,URHO)
               Upo     = uold(i,j,k,UMX) * rhooinv
               Vpo     = uold(i,j,k,UMY) * rhooinv
               Wpo     = uold(i,j,k,UMZ) * rhooinv

               ! Define old source terms
               SrU_old = rhoo * gold(i,j,k,1)
               SrV_old = rhoo * gold(i,j,k,2)
               SrW_old = rhoo * gold(i,j,k,3)

               rhon    = unew(i,j,k,URHO)
               rhoninv = 1.0d0 / unew(i,j,k,URHO)
               Upn     = unew(i,j,k,UMX) * rhoninv
               Vpn     = unew(i,j,k,UMY) * rhoninv
               Wpn     = unew(i,j,k,UMZ) * rhoninv

               ! Define new source terms
               SrU_new = rhon * gnew(i,j,k,1)
               SrV_new = rhon * gnew(i,j,k,2)
               SrW_new = rhon * gnew(i,j,k,3)

               ! Define corrections to source terms
               SrUcorr = 0.5d0*(SrU_new - SrU_old)
               SrVcorr = 0.5d0*(SrV_new - SrV_old)
               SrWcorr = 0.5d0*(SrW_new - SrW_old)

               ! This does work (in 1-d)
               if (grav_source_type .eq. 1) then
                   SrEcorr =  0.5d0 * ( (SrU_new * Upn - SrU_old * Upo) + &
                                        (SrV_new * Vpn - SrV_old * Vpo) + &
                                        (SrW_new * Wpn - SrW_old * Wpo) )
               end if

               ! Correct state with correction terms
               unew(i,j,k,UMX)   = unew(i,j,k,UMX)   + SrUcorr*dt * a_new_inv
               unew(i,j,k,UMY)   = unew(i,j,k,UMY)   + SrVcorr*dt * a_new_inv
               unew(i,j,k,UMZ)   = unew(i,j,k,UMZ)   + SrWcorr*dt * a_new_inv

               if (grav_source_type .eq. 1) then
                   unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr*dt * (a_half / a_newsq)
               else if (grav_source_type .eq. 3) then
                   new_ke = 0.5d0 * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                     unew(i,j,k,URHO) 
                   unew(i,j,k,UEDEN) = old_rhoeint + new_ke
#ifndef AMREX_USE_CUDA
               else 
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
#endif
               end if

            enddo
         enddo
      enddo

      end subroutine fort_correct_gsrc
