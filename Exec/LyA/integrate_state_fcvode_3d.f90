subroutine integrate_state_fcvode(lo, hi, &
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
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   TEMP_COMP, NE_COMP, gamma_minus_1
    use bl_constants_module, only: M_PI
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use atomic_rates_module, only: tabulate_rates, interp_to_this_z, YHELIUM
    use vode_aux_module    , only: z_vode, i_vode, j_vode, k_vode
    use cvode_interface
    use fnvector_serial
    use fcvode_extras
    use, intrinsic :: iso_c_binding

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    real(rt) :: z, rho
    real(rt) :: T_orig, ne_orig, e_orig
    real(rt) :: T_out , ne_out , e_out, mu, mean_rhob
    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: tstart     ! initial time
    real(c_double) :: atol, rtol
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: CVmem         ! CVODE memory
    integer(c_long), parameter :: neq = 1
    real(c_double), pointer :: yvec(:)

    allocate(yvec(neq))

    z = 1.d0/a - 1.d0

    z_vode = z
    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Interpolate from the table to this redshift
    call interp_to_this_z(z)

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    sunvec_y = N_VMake_Serial(NEQ, yvec)
    if (.not. c_associated(sunvec_y)) then
        call amrex_abort('integrate_state_fcvode: sunvec = NULL')
    end if

    CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
    if (.not. c_associated(CVmem)) then
        call amrex_abort('integrate_state_fcvode: CVmem = NULL')
    end if

    tstart = 0.0
    ! CVodeMalloc allocates variables and initialize the solver. We can initialize the solver with junk because once we enter the
    ! (i,j,k) loop we will immediately call fcvreinit which reuses the same memory allocated from CVodeMalloc but sets up new
    ! initial conditions.
    ierr = FCVodeInit(CVmem, c_funloc(RhsFn), tstart, sunvec_y)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVodeInit() failed')
    end if

    ! Set dummy tolerances. These will be overwritten as soon as we enter the loop and reinitialize the solver.
    rtol = 1.0d-5
    atol = 1.0d-10
    ierr = FCVodeSStolerances(CVmem, rtol, atol)
    if (ierr /= 0) then
      call amrex_abort('integrate_state_fcvode: FCVodeSStolerances() failed')
    end if

    ierr = FCVDense(CVmem, neq)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVDense() failed')
    end if

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)

                ! Original values
                rho     = state(i,j,k,URHO)
                e_orig  = state(i,j,k,UEINT) / rho
                T_orig  = diag_eos(i,j,k,TEMP_COMP)
                ne_orig = diag_eos(i,j,k,  NE_COMP)

                if (e_orig .lt. 0.d0) then
                    print *,'negative e entering strang integration ',z, i,j,k, rho/mean_rhob, e_orig
                    call bl_abort('bad e in strang')
                end if

                i_vode = i
                j_vode = j
                k_vode = k

                call fcvode_wrapper(half_dt,rho,T_orig,ne_orig,e_orig,neq,CVmem,sunvec_y,yvec, &
                                              T_out ,ne_out ,e_out)

                if (e_out .lt. 0.d0) then
                    print *,'negative e exiting strang integration ',z, i,j,k, rho/mean_rhob, e_out
                    T_out  = 10.0
                    ne_out = 0.0
                    mu     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)
                    e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                    call flush(6)
!                    call bl_abort('bad e out of strang')
                end if

                ! Update (rho e) and (rho E)
                state(i,j,k,UEINT) = state(i,j,k,UEINT) + rho * (e_out-e_orig)
                state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + rho * (e_out-e_orig)

                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_T_given_Re(T_out, ne_out, rho, e_out, a)
                diag_eos(i,j,k,TEMP_COMP) = T_out
                diag_eos(i,j,k,  NE_COMP) = ne_out

            end do ! i
        end do ! j
    end do ! k

    call N_VDestroy_Serial(sunvec_y)
    call FCVodeFree(cvmem)

    deallocate(yvec)

end subroutine integrate_state_fcvode
