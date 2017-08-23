
      ! *************************************************************************************

      subroutine fort_compute_temp(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp")

      use amrex_fort_module, only : rt => amrex_real
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
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k
      real(rt) :: rhoInv,eint
      real(rt) :: ke,dummy_pres
      real(rt) :: z

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

      end subroutine fort_compute_temp


      subroutine fort_compute_temp_vec(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp_vec")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use atomic_rates_module, only: this_z, interp_to_this_z
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, &
                                     TEMP_COMP, NE_COMP, small_temp, heat_cool_type
      use  eos_params_module
      use misc_params, only: simd_width

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k
      real(rt) :: rhoInv,eint
      real(rt) :: ke,dummy_pres
      real(rt) :: z

      ! Data structures holding temporary state data so we can call the EOS and
      ! HeatCool routines using SIMD. One data structure holds all cells with positive
      ! energies (which call nyx_eos_T_given_Re()), and another holds cells with
      ! negative energies, which call nyx_eos_given_RT()).
      type eint_neg_type
        integer, allocatable :: i(:), j(:), k(:)
        double precision, allocatable :: temp(:), ne(:), urho(:), eint(:)
      end type eint_neg_type
      type eint_pos_type
        integer, allocatable :: i(:), j(:), k(:)
        double precision, allocatable :: temp(:), ne(:), urho(:), eint(:)
      end type eint_pos_type

      type (eint_neg_type) :: eint_neg_list
      type (eint_pos_type) :: eint_pos_list

      integer :: loop_size, eint_neg_counter, eint_pos_counter

      loop_size = (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1) * (hi(3)-lo(3)+1)

      eint_neg_counter = 1
      eint_pos_counter = 1

      allocate(eint_neg_list%i(loop_size))
      allocate(eint_neg_list%j(loop_size))
      allocate(eint_neg_list%k(loop_size))
      allocate(eint_neg_list%temp(loop_size))
      allocate(eint_neg_list%ne(loop_size))
      allocate(eint_neg_list%urho(loop_size))
      allocate(eint_neg_list%eint(loop_size))

      allocate(eint_pos_list%i(loop_size))
      allocate(eint_pos_list%j(loop_size))
      allocate(eint_pos_list%k(loop_size))
      allocate(eint_pos_list%temp(loop_size))
      allocate(eint_pos_list%ne(loop_size))
      allocate(eint_pos_list%urho(loop_size))
      allocate(eint_pos_list%eint(loop_size))

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

      ! Cells which call nyx_eos_given_RT() and nyx_eos_T_given_Re to get equation of
      ! state data now need to call it with vector arguments, since both the integrator
      ! (RKF45) and the RHS functions (in particular, iterate_ne() and ion_n()) are now
      ! vectorized. So we do this in 2 steps:
      ! 
      ! 1.) Build 2 lists of cells, one which contains cells with
      ! positive internal energy, and another with cells with negative energy.
      ! 
      ! 2.) Loop through each list separately and call the two EOS functions with
      ! vector arguments.

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rhoInv = 1.d0 / state(i,j,k,URHO)
               if (state(i,j,k,UEINT) > 0.d0) then
                 eint = state(i,j,k,UEINT) * rhoInv
                 eint_pos_list%i(eint_pos_counter) = i
                 eint_pos_list%j(eint_pos_counter) = j
                 eint_pos_list%k(eint_pos_counter) = k
                 eint_pos_list%temp(eint_pos_counter) = diag_eos(i,j,k,TEMP_COMP)
                 eint_pos_list%ne(eint_pos_counter) = diag_eos(i,j,k,NE_COMP)
                 eint_pos_list%urho(eint_pos_counter) = state(i,j,k,URHO)
                 eint_pos_list%eint(eint_pos_counter) = eint
                 eint_pos_counter = eint_pos_counter + 1
               else
                 eint_neg_list%i(eint_neg_counter) = i
                 eint_neg_list%j(eint_neg_counter) = j
                 eint_neg_list%k(eint_neg_counter) = k
                 eint_neg_list%temp(eint_neg_counter) = diag_eos(i,j,k,TEMP_COMP)
                 eint_neg_list%ne(eint_neg_counter) = diag_eos(i,j,k,NE_COMP)
                 eint_neg_list%urho(eint_neg_counter) = state(i,j,k,URHO)
                 eint_neg_list%eint(eint_neg_counter) = eint
                 eint_neg_counter = eint_neg_counter + 1
               end if
            end do
         end do
      end do

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoInv = 1.d0 / state(i,j,k,URHO)

               if (state(i,j,k,UEINT) > 0.d0) then

                   eint = state(i,j,k,UEINT) * rhoInv

               else
                  if (print_fortran_warnings .gt. 0) then
                     print *,'   '
                     print *,'>>> Warning: (rho e) is negative in compute_temp: ',i,j,k
                  end if
                   ! Set temp to small_temp and compute corresponding internal energy

                   ke = 0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) * rhoInv

                   diag_eos(i,j,k,TEMP_COMP) = small_temp
                   state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
                   state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ke

               end if

            enddo
         enddo
      enddo

      deallocate(eint_neg_list%i)
      deallocate(eint_neg_list%j)
      deallocate(eint_neg_list%k)
      deallocate(eint_neg_list%temp)
      deallocate(eint_neg_list%ne)
      deallocate(eint_neg_list%urho)
      deallocate(eint_neg_list%eint)

      deallocate(eint_pos_list%i)
      deallocate(eint_pos_list%j)
      deallocate(eint_pos_list%k)
      deallocate(eint_pos_list%temp)
      deallocate(eint_pos_list%ne)
      deallocate(eint_pos_list%urho)
      deallocate(eint_pos_list%eint)

      end subroutine fort_compute_temp_vec

      subroutine fort_compute_rho_temp(lo,hi,dx, &
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                  diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  rho_ave,rho_T_sum, &
                                  T_sum,T_meanrho_sum,rho_sum,vol_sum,vol_mn_sum) &
      bind(C, name = "fort_compute_rho_temp")

      use meth_params_module, only : NVAR, URHO, TEMP_COMP

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho_ave
      real(rt), intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(in   ) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      real(rt), intent(inout) :: rho_T_sum, rho_sum, T_sum, T_meanrho_sum
      real(rt), intent(inout) :: vol_sum, vol_mn_sum

      integer          :: i,j,k
      real(rt) :: rho_hi, rho_lo, vol

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

      end subroutine fort_compute_rho_temp

      subroutine fort_compute_max_temp_loc(lo,hi, &
                                           state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                           diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                           max_temp, den_maxt, imax, jmax, kmax) &
      bind(C, name = "fort_compute_max_temp_loc")

      use meth_params_module, only : TEMP_COMP, NVAR, URHO

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      real(rt), intent(in   ) :: max_temp
      real(rt), intent(  out) :: den_maxt
      integer         , intent(inout) :: imax,jmax,kmax

      integer                         :: i,j,k
      real(rt)                :: one_minus_eps

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

      end subroutine fort_compute_max_temp_loc
