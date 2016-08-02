
      ! *************************************************************************************

      subroutine compute_temp(lo,hi, &
                              state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                              diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                              comoving_a, print_fortran_warnings)

      use eos_module
      use atomic_rates_module, only: this_z, interp_to_this_z
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, &
                                     TEMP_COMP, NE_COMP, small_temp, heat_cool_type
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      double precision, intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(in   ) :: comoving_a

      integer          :: i,j,k
      double precision :: rhoInv,eint
      double precision :: ke,dummy_pres
      double precision :: z
      integer          :: pt_index(3)

      z = 1.d0/comoving_a - 1.d0

      if (heat_cool_type.gt.0) then
          if (z .ne. this_z) &
             call interp_to_this_z(z)
      end if

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,URHO) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,URHO)
                  print *,'    '
                  call bl_error("Error:: compute_temp_3d.f90 :: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoInv = 1.d0 / state(i,j,k,URHO)

               if (state(i,j,k,UEINT) > 0.d0) then

                   eint = state(i,j,k,UEINT) * rhoInv

                   pt_index(1) = i
                   pt_index(2) = j
                   pt_index(3) = k

                   call nyx_eos_T_given_Re(diag_eos(i,j,k,TEMP_COMP), diag_eos(i,j,k,NE_COMP), &
                                           state(i,j,k,URHO), eint, comoving_a)

               else
                  if (print_fortran_warnings .gt. 0) then
                     print *,'   '
                     print *,'>>> Warning: (rho e) is negative in compute_temp: ',i,j,k
                  end if
                   ! Set temp to small_temp and compute corresponding internal energy
                   call nyx_eos_given_RT(eint, dummy_pres, state(i,j,k,URHO), small_temp, &
                                         diag_eos(i,j,k,NE_COMP), comoving_a)

                   ke = 0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) * rhoInv

                   diag_eos(i,j,k,TEMP_COMP) = small_temp
                   state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
                   state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ke

               end if

            enddo
         enddo
      enddo

      end subroutine compute_temp

      subroutine compute_rho_temp(lo,hi,dx, &
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                  diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  rho_ave,rho_T_sum, &
                                  T_sum,T_meanrho_sum,rho_sum,vol_sum,vol_mn_sum)

      use meth_params_module, only : NVAR, URHO, TEMP_COMP

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      double precision, intent(in   ) :: dx(3)
      double precision, intent(in   ) :: rho_ave
      double precision, intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(in   ) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(inout) :: rho_T_sum, rho_sum, T_sum, T_meanrho_sum
      double precision, intent(inout) :: vol_sum, vol_mn_sum

      integer          :: i,j,k
      double precision :: rho_hi, rho_lo, vol

      vol = dx(1)*dx(2)*dx(3)
      rho_hi = 1.1d0*rho_ave
      rho_lo = 0.9d0*rho_ave
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                   T_sum =     T_sum + vol*diag_eos(i,j,k,TEMP_COMP)
               rho_T_sum = rho_T_sum + state(i,j,k,URHO)*diag_eos(i,j,k,TEMP_COMP)
                 rho_sum =   rho_sum + state(i,j,k,URHO)
                 if ( (state(i,j,k,URHO) .lt. rho_hi) .and. &
                      (state(i,j,k,URHO) .gt. rho_lo) .and. &
                      (diag_eos(i,j,k,TEMP_COMP) .le. 1.0e5) ) then
                         T_meanrho_sum = T_meanrho_sum + vol*dlog10(diag_eos(i,j,k,TEMP_COMP))
                         vol_mn_sum = vol_mn_sum + vol
                 endif
                 vol_sum = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine compute_rho_temp

      subroutine compute_max_temp_loc(lo,hi, &
                                      state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                      diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                      max_temp, den_maxt, imax, jmax, kmax)

      use meth_params_module, only : TEMP_COMP, NVAR, URHO

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      double precision, intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(in   ) :: max_temp
      double precision, intent(  out) :: den_maxt
      integer         , intent(inout) :: imax,jmax,kmax

      integer                         :: i,j,k
      double precision                :: one_minus_eps

      one_minus_eps = 1.d0 - 1.d-12

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (diag_eos(i,j,k,TEMP_COMP) .ge. one_minus_eps*max_temp) then
                  imax = i
                  jmax = j
                  kmax = k
                  den_maxt = state(i,j,k,URHO)
               end if
            enddo
         enddo
      enddo

      end subroutine compute_max_temp_loc
