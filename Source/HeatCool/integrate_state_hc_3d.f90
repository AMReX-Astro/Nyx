subroutine integrate_state_hc(lo, hi, &
                              state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                              diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                              a, half_dt, min_iter, max_iter)
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
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, small_pres, gamma_minus_1
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use atomic_rates_module, only: tabulate_rates
    use heating_cooling_module, only: hc_rates

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer, parameter :: NITERS = 20
    real(rt), parameter :: xacc = 1.0d-3

    integer :: i, j, k, n, iter, nsteps, cnt
    real(rt) :: z, rho, T, ne
    real(rt) :: T_orig, rho_e_orig, ne_orig, e_int_old, De_int
    real(rt) :: T_first, ne_first, src_first
    real(rt) :: src_old, src_new, delta_re, delta_t, rho_e, e_int, prev_soln
    real(rt) :: b_fac
    logical          :: do_diag, prnt_cell, done_iter
    logical          :: went_negative, went_negative_at_first

    z = 1.d0/a - 1.d0
    do_diag   = .false.
    prnt_cell = .false.

    b_fac = 0.0d0
    max_iter = 0
    min_iter = NITERS+1

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                ! Original values
                rho        = state(i,j,k,URHO)
                rho_e_orig = state(i,j,k,UEINT)
                T_orig     = diag_eos(i,j,k,TEMP_COMP)
                ne_orig    = diag_eos(i,j,k,  NE_COMP)

                if (rho_e_orig .lt. 0.d0) then
                    print *,'(rho e) entering strang integration negative ',i,j,k, rho_e_orig
                    call bl_abort('bad rho e in strang')
                end if

                e_int = rho_e_orig/rho
                call hc_rates(z, rho, e_int, T, ne, src_new, prnt_cell)
                T_first   = T
                ne_first  = ne
                src_first = src_new

                went_negative_at_first = .false.

                prev_soln = HUGE(prev_soln)
                do iter = 1, NITERS  ! max allowed iterations

                  nsteps  = 2**(iter-1)
                  delta_t = half_dt / nsteps
                  rho_e   = rho_e_orig
                  e_int   = rho_e/rho

                  delta_re = 0.0d0
                  src_old  = 0.0d0

                  do n = 1, nsteps

                    done_iter = .false.
                    e_int_old = e_int
                    e_int     = rho_e/rho

                    if (n.eq.1) then
                        T      = T_first
                       ne      = ne_first
                       src_new = src_first
                    else
                       call hc_rates(z, rho, e_int, T, ne, src_new, prnt_cell)
                    end if

                    if ( (rho_e+delta_t*src_new/a) .gt. 0.0d0) then 

                       went_negative          = .false.     

                       src_new = delta_t * src_new / a
                       rho_e   = rho_e + src_new

                       if (src_old*src_new .lt. 0.0d0) then ! src=0 in between
                          if (rho_e .le. 0.0d0) then 
                             De_int = e_int/2.0d0
                             e_int  = e_int/2.0d0
                          else
                             De_int = abs(e_int_old - e_int)/2.0d0
                             e_int  = e_int + sign(De_int, src_new)
                          endif
                          cnt = 0
                          do
                             cnt = cnt + 1
                             call hc_rates(z, rho, e_int, T, ne, src_new, prnt_cell)
                             if (abs(delta_t*src_new/a)/rho .lt. xacc) EXIT
                             if (cnt .gt. 40) then 
                                print*, 'BISECTION problem in cell:',i,j,k,iter,n
                                call flush(6)
                                call bl_error("Problem in bisection in integrate_hc_3d.f90")
                             endif
                             De_int = De_int / 2.0d0
                             e_int = e_int + sign(De_int, src_new)
                          enddo

                          rho_e     = e_int * rho
                          delta_re  = rho_e - rho_e_orig
                          done_iter = .true.
                          b_fac     = b_fac+1 ! just for diagnostics
                          EXIT
                       endif

                       delta_re = delta_re + src_new ! Cumulative update
                       src_old  = src_new

                    ! Here we just leave rho_e alone and proceed to the next iter
                    else   ! (rho_e + src_new) <= 0
                       went_negative = .true.     
                       went_negative_at_first = .true.     
                       ! print *,'WENT NEGATIVE n, nsteps, iter ',n, nsteps, iter
                       ! print *,' at cell ',i,j,k
                       ! print *,' rho_e_orig       ',rho_e_orig 
                       ! print *,' src              ',src_new
                       ! print *,' dt/a             ',delta_t/a
                       ! print *,' src*dt/a         ',src_new*delta_t/a
                       ! print *,'  '
                       ! call flush(6)
                       exit  ! Exit the n loop to go to higher nsteps
                    end if

                  enddo ! n loop

                  if (.not. went_negative) then
                     if (abs(1.0d0-rho_e/prev_soln) .lt. xacc .or. done_iter) EXIT
                  end if

                  if (iter .ge. NITERS-2) then 
                     print*, 'INTEGRATE_HC ITERATIONS:', i,j,k, iter, rho_e_orig, rho_e, (1.0d0-rho_e/prev_soln)
                     call flush(6)
                  endif

                  if (iter .eq. NITERS) then 
                     print*, 'MAXITER too small!', i,j,k, rho, T_orig
                     call bl_abort('too small MAXITER')
                  endif

                  if (.not. went_negative) &
                      prev_soln = rho_e
                enddo ! iter loop

                ! Update cell quantities
                state(i,j,k,UEINT) = state(i,j,k,UEINT) + delta_re
                state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + delta_re
                diag_eos(i,j,k,TEMP_COMP) = T
                diag_eos(i,j,k,  NE_COMP) = ne

                if (state(i,j,k,UEINT) .lt. 0.d0) then
                    print *,'(rho e) exiting strang integration negative ',i,j,k, rho, rho_e_orig/rho, state(i,j,k,UEINT)/rho
                    call bl_abort('negative rho e exiting strang')
                end if

                if (state(i,j,k,UEINT) .lt. small_pres/gamma_minus_1) then
                   print *,'!!!!! Pressure and (rho e) are too small coming out of integrate !!!!!'
                   print *,'!!!!! (i,j,k) !!!!! ' ,i,j,k
                   print *,'!!!!! pressure      ',state(i,j,k,UEINT) * gamma_minus_1
                   call flush(6)
                end if

                if (max_iter .le. iter) max_iter = iter
                if (min_iter .ge. iter) min_iter = iter

                ! if (went_negative_at_first) print *,'MAX NSTEPS OF NEG AT ',i,j,k,nsteps

            end do ! i
        end do ! j
    end do ! k

    if (do_diag) then 
      print*, 'HC_ITERATIONS: ', z, max_iter, min_iter, b_fac/((hi(3)-lo(3))*(hi(2)-lo(2))*(hi(1)-lo(1)))
      call flush(6)
    endif

end subroutine integrate_state_hc

