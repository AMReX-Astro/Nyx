subroutine integrate_state_force(lo, hi, &
                                 state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                 diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                 a, half_dt)
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
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   TEMP_COMP, NE_COMP, small_pres, gamma_minus_1
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use atomic_rates_module, only: tabulate_rates, interp_to_this_z
 
    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    double precision, intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    double precision, intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    double precision, intent(in)    :: a, half_dt

    integer, parameter :: NITERS = 20
    double precision, parameter :: xacc = 1.0d-3

    integer :: i, j, k, n, iter, nsteps, cnt
    double precision :: z, rho, T, ne
    double precision :: T_orig, rho_e_orig, ne_orig, e_int_old, De_int
    double precision :: T_first, ne_first, src_first
    double precision :: src_old, src_new, delta_re, delta_t, rho_e, e_int, prev_soln
    double precision :: b_fac
    logical          :: do_diag, prnt_cell, done_iter
    logical          :: went_negative, went_negative_at_first

    z = 1.d0/a - 1.d0
    do_diag   = .false.
    prnt_cell = .false.

    ! Interpolate from the table to this redshift
    call interp_to_this_z(z)

    b_fac = 0.0d0
 
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

		src_new = 0.0d0

                e_int = rho_e_orig/rho
                if ( (rho_e+delta_t*src_new/a) .gt. 0.0d0) then 
                    src_new = delta_t * src_new / a
                    e_int   = e_int + src_new
                end if
                rho_e     = e_int * rho
                delta_re  = rho_e - rho_e_orig

		T = T_orig
		ne = ne_orig

                ! Update cell quantities
                state(i,j,k,UEINT) = state(i,j,k,UEINT) + delta_re
                state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + delta_re
                diag_eos(i,j,k,TEMP_COMP) = T
                diag_eos(i,j,k,  NE_COMP) = ne
            end do ! i
        end do ! j
    end do ! k

end subroutine integrate_state_force

