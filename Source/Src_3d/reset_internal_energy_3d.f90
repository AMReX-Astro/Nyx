      subroutine reset_internal_e(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3, &
                                  d,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3,lo,hi, &
                                  print_fortran_warnings,&
                                  comoving_a,sum_energy_added,sum_energy_total)

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, small_temp, &
                                     NE_COMP
      use  eos_params_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: print_fortran_warnings
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
      integer          :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
      real(rt) :: d(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)

      real(rt), intent(in   ) :: comoving_a
      real(rt), intent(inout) :: sum_energy_added
      real(rt), intent(inout) :: sum_energy_total

      ! Local variables
      integer          :: i,j,k
      real(rt) :: Up, Vp, Wp, ke, rho_eint, eint_new
      real(rt) :: dummy_pres, rhoInv

      ! Reset internal energy if necessary
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

           rhoInv = 1.0d0 / u(i,j,k,URHO)
           Up     = u(i,j,k,UMX) * rhoInv
           Vp     = u(i,j,k,UMY) * rhoInv
           Wp     = u(i,j,k,UMZ) * rhoInv
           ke     = 0.5d0 * u(i,j,k,URHO) * (Up*Up + Vp*Vp + Wp*Wp)

           rho_eint = u(i,j,k,UEDEN) - ke

           ! Reset (e from e) if it's greater than 0.01% of big E.
           if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,j,k,UEDEN) .gt. 1.d-6) then

               u(i,j,k,UEINT) = rho_eint

           ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
           else if (u(i,j,k,UEINT) .gt. 0.d0) then

              ! Keep track of how much energy we are adding to (rho E)
              sum_energy_added = sum_energy_added + (u(i,j,k,UEINT) + ke - u(i,j,k,UEDEN))

              u(i,j,k,UEDEN) = u(i,j,k,UEINT) + ke

           ! If not resetting and little e is negative ...
           else if (u(i,j,k,UEINT) .le. 0.d0) then

              call nyx_eos_given_RT(eint_new, dummy_pres, u(i,j,k,URHO), small_temp, &
                                    d(i,j,k,NE_COMP),comoving_a)

              if (print_fortran_warnings .gt. 0) then
                 print *,'   '
                 print *,'>>> Warning: Nyx_3d::reset_internal_energy  ',i,j,k
                 print *,'>>> ... Resetting neg. e from EOS using small_temp: ',small_temp,&
                         ' from ',u(i,j,k,UEINT)/u(i,j,k,URHO),' to ', eint_new
                 call flush(6)
              end if

              u(i,j,k,UEINT) = u(i,j,k,URHO) *  eint_new

              ! Keep track of how much energy we are adding to (rho E)
              sum_energy_added = sum_energy_added + (u(i,j,k,UEINT) + ke - u(i,j,k,UEDEN))

              u(i,j,k,UEDEN) = u(i,j,k,UEINT) + ke

           end if

           sum_energy_total = sum_energy_total + u(i,j,k,UEDEN)

      enddo
      enddo
      enddo

      end subroutine reset_internal_e
