subroutine integrate_state_force(lo, hi, &
                                 state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                 diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                 dx, time, a, half_dt)
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   state_* : double arrays
!       The state vars
!   diag_eos_* : double arrays
!       Temp and Ne
!   src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   half_dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_fort_module, only : rt => amrex_real

    use forcing_spect_module
    use eos_params_module
    use eos_module, only: nyx_eos_given_RT, nyx_eos_T_given_Re
    use probdata_module, only: prob_lo, prob_hi, alpha, rho0, temp0
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   TEMP_COMP, NE_COMP, small_pres, small_temp, gamma_minus_1
    use bl_constants_module, only : TWO, ONE, HALF, ZERO, M_PI, M_SQRT_2
    use fundamental_constants_module
 
    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    real(rt), intent(in)    :: dx(3), time, a, half_dt

    integer :: i, j, k
    integer :: m, mi, mj, mk, n
    integer :: alloc
    real(rt) :: rho, rho_e_orig, rho_K_res, T_orig, ne
    real(rt) :: delta, delta_re, eint0, press, small_eint

    integer :: num_phases(3)
    real(rt) :: delta_phase(3), phase_lo(3)
    real(rt) :: accel(3), buf(num_modes) 
    real(rt) :: phasefct_init_even(num_modes), phasefct_init_odd(num_modes)
    real(rt) :: phasefct_mult_even(num_modes,3), phasefct_mult_odd(num_modes,3)
    real(rt), allocatable :: phasefct_even_x(:), phasefct_even_y(:), phasefct_even_z(:)
    real(rt), allocatable :: phasefct_odd_x(:), phasefct_odd_y(:), phasefct_odd_z(:)

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    print *, "integrate_state_force lo = ", lo
    print *, "integrate_state_force hi = ", hi

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                ! Original values
                rho        = state(i,j,k,URHO)
                rho_e_orig = state(i,j,k,UEINT)
		rho_K_res  = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
                T_orig     = diag_eos(i,j,k,TEMP_COMP)
                ne         = diag_eos(i,j,k,  NE_COMP)

                if (rho_e_orig .lt. 0.d0) then
                    print *,'(rho e) entering strang integration negative ',i,j,k, rho_e_orig
                    call bl_abort('bad rho e in strang')
                end if

                ! Compute temperature increment and ensure that new temperature is positive
		delta = half_dt * alpha * (temp0 - T_orig) / a
		diag_eos(i,j,k,TEMP_COMP) = max(T_orig + delta, small_temp)

                ! Call EOS to get internal energy for constant equilibrium temperature
                call nyx_eos_given_RT(eint0, press, rho, temp0, ne, a)
		delta_re = half_dt * alpha * (rho*eint0 - rho_e_orig) / a

                ! Call EOS to get the internal energy floor
                call nyx_eos_given_RT(small_eint, press, rho, small_temp, ne, a)
		
                ! Update cell quantities
 		state(i,j,k,UEINT) = max(rho_e_orig + delta_re, rho*small_eint)
                state(i,j,k,UEDEN) = state(i,j,k,UEINT) + rho_K_res

		!if ((i.eq.16).and.(j.eq.16)) then
                !   print *, "temp: ", k, ne, temp0, T_orig, diag_eos(i,j,k,TEMP_COMP), delta
                !   print *, "rhoe: ", k, rho, rho*eint0, rho_e_orig, state(i,j,k,UEINT), delta_re
                !endif
            end do ! i
        end do ! j
    end do ! k

    !print *, " --- integrate_state_force --- "
    delta_phase(:) = TWO*M_PI * dx(:) / (prob_hi(:) - prob_lo(:)) ! phase increment per cell
    !print *, "dx = ", dx
    !print *, "prob_hi = ", prob_hi
    !print *, "prob_lo = ", prob_lo
    !print *, "delta_phase = ", delta_phase
    phase_lo(:) = (dble(lo(:)) + HALF) * delta_phase(:)              ! phase of low corner
    !print *, "lo = ", lo
    !print *, "hi = ", hi
    !print *, "phase_lo = ", phase_lo

    ! compute initial phase factors and multiplying factors
    do m = 1, num_modes
       i = wavevectors(1,m)
       j = wavevectors(2,m)
       k = wavevectors(3,m)

       phasefct_init_even(m) = &
          (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
           sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * cos(k*phase_lo(3)) - &
          (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
           sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * sin(k*phase_lo(3))

       phasefct_init_odd(m) = &
          (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
           sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * sin(k*phase_lo(3)) + &
          (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
           sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * cos(k*phase_lo(3))

       !print *, m, i*phase_lo(1), j*phase_lo(2), k*phase_lo(3), phasefct_init_even(m), phasefct_init_odd(m)

       phasefct_mult_even(m,1) = cos(i*delta_phase(1));
       phasefct_mult_odd (m,1) = sin(i*delta_phase(1));

       phasefct_mult_even(m,2) = cos(j*delta_phase(2));
       phasefct_mult_odd (m,2) = sin(j*delta_phase(2));

       phasefct_mult_even(m,3) = cos(k*delta_phase(3));
       phasefct_mult_odd (m,3) = sin(k*delta_phase(3));
       !print *, m, i*delta_phase(1), phasefct_mult_even(m,1), phasefct_mult_even(m,1)
       !print *, m, j*delta_phase(2), phasefct_mult_even(m,2), phasefct_mult_even(m,2)
       !print *, m, k*delta_phase(3), phasefct_mult_even(m,3), phasefct_mult_even(m,3)
    end do

    num_phases(:) = (hi(:)-lo(:)+1)*num_modes
    print *, "integrate_state_force num_phases = ", num_phases

    allocate(phasefct_even_x(num_phases(1)), phasefct_even_y(num_phases(2)), phasefct_even_z(num_phases(3)), &
             phasefct_odd_x(num_phases(1)),  phasefct_odd_y(num_phases(2)),  phasefct_odd_z(num_phases(3)), &
             STAT=alloc)

    if (alloc > 0) call bl_abort('failed to allocate arrays for phase factors')      
 
    ! initialize phase factors for each coordinate axis
    do m = 1, num_modes
       phasefct_even_x(m) = ONE 
       phasefct_odd_x(m)  = ZERO
       phasefct_even_y(m) = ONE 
       phasefct_odd_y(m)  = ZERO 
       phasefct_even_z(m) = phasefct_init_even(m)  
       phasefct_odd_z(m)  = phasefct_init_odd(m)
    end do

    do i = lo(1)+1,hi(1)
       mi = (i-lo(1))*num_modes + 1
       do m = 1, num_modes
            buf(m) = phasefct_even_x(mi-num_modes);
            phasefct_even_x(mi) = phasefct_mult_even(m,1) * phasefct_even_x(mi-num_modes) - &
                                  phasefct_mult_odd (m,1) * phasefct_odd_x(mi-num_modes)
            phasefct_odd_x(mi)  = phasefct_mult_even(m,1) * phasefct_odd_x(mi-num_modes) + &
                                  phasefct_mult_odd (m,1) * buf(m)
            mi = mi + 1
       end do
    end do         

    do j = lo(2)+1,hi(2)
       mj = (j-lo(2))*num_modes + 1
       do m = 1, num_modes
            buf(m) = phasefct_even_y(mj-num_modes);
            phasefct_even_y(mj) = phasefct_mult_even(m,2) * phasefct_even_y(mj-num_modes) - &
                                  phasefct_mult_odd (m,2) * phasefct_odd_y(mj-num_modes)
            phasefct_odd_y(mj)  = phasefct_mult_even(m,2) * phasefct_odd_y(mj-num_modes) + &
                                  phasefct_mult_odd (m,2) * buf(m)
            mj = mj + 1
       end do
    end do         

    do k = lo(3)+1, hi(3)
       mk = (k-lo(3))*num_modes + 1
       do m = 1, num_modes
            buf(m) = phasefct_even_z(mk-num_modes);
            phasefct_even_z(mk) = phasefct_mult_even(m,3) * phasefct_even_z(mk-num_modes) - &
                                  phasefct_mult_odd (m,3) * phasefct_odd_z(mk-num_modes)
            phasefct_odd_z(mk)  = phasefct_mult_even(m,3) * phasefct_odd_z(mk-num_modes) + &
                                  phasefct_mult_odd (m,3) * buf(m)
            mk = mk + 1
       end do
    end do

    ! apply forcing in physical space
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             accel(:) = (/ ZERO, ZERO, ZERO /)

             ! compute components of acceleration via inverse FT
             do n = 1, 3
                mi = (i-lo(1))*num_modes + 1 ! offset in x-direction
                mj = (j-lo(2))*num_modes + 1 ! offset in y-direction
                mk = (k-lo(3))*num_modes + 1 ! offset in z-direction
  
                do m = 1, num_modes
                   ! sum up even modes
                   accel(n) = accel(n) + &
                              ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                                phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_even_z(mk) - &
                               (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                                phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_odd_z(mk))  * modes_even(m,n)
                   ! sum up odd modes
                   accel(n) = accel(n) - &
                              ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                                phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_odd_z(mk)  + &
                               (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                                phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_even_z(mk)) * modes_odd(m,n)
                   mi = mi + 1
                   mj = mj + 1
                   mk = mk + 1
                end do
             end do

             accel(:) = M_SQRT_2 * accel(:)

             ! add forcing to state		
             state(i,j,k,UMX) = state(i,j,k,UMX) + half_dt * state(i,j,k,URHO)*accel(1) / a
             state(i,j,k,UMY) = state(i,j,k,UMY) + half_dt * state(i,j,k,URHO)*accel(2) / a
             state(i,j,k,UMZ) = state(i,j,k,UMZ) + half_dt * state(i,j,k,URHO)*accel(3) / a
          end do
       end do
    end do

    ! update total energy
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5d0*(state(i,j,k,UMX)*state(i,j,k,UMX) + &
                                                              state(i,j,k,UMY)*state(i,j,k,UMY) + &
                                                              state(i,j,k,UMZ)*state(i,j,k,UMZ))/state(i,j,k,URHO)
          end do
       end do
    end do

    deallocate(phasefct_even_x, phasefct_even_y, phasefct_even_z, &
               phasefct_odd_x,  phasefct_odd_y,  phasefct_odd_z)

  end subroutine integrate_state_force

