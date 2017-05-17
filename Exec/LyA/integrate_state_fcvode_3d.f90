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

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    double precision :: rout(10), rpar
    integer :: i, j, k
    real(rt) :: z, rho
    real(rt) :: T_orig, ne_orig, e_orig
    real(rt) :: T_out , ne_out , e_out, mu, mean_rhob
    integer*8 :: NEQ = 1, ipar, iout(25)
    integer :: meth, itmeth, iatol
    integer :: ier
    double precision   :: time
    double precision :: y(1)
    double precision :: atol(1), rtol(1)

    z = 1.d0/a - 1.d0

    z_vode = z
    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Interpolate from the table to this redshift
    call interp_to_this_z(z)

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    call fnvinits(1, NEQ, ier)

    ! fcvmalloc calls CVodeMalloc to allocate variables and initialize the solver. We can initialize the solver with junk because
    ! once we enter the (i,j,k) loop we will immediately call fcvreinit which reuses the same memory allocated from fcvmalloc but
     ! sets up new initial conditions.
    meth = 1   ! 1 = Adams (non-stiff); 2 = BDF (stiff)
    itmeth = 1 ! 1 = functional iteration; 2 = Newton iteration
    iatol = 1
    rpar = 42.0d0
    ipar = 42
    time = 0.0

    call fcvmalloc(time, y, meth, itmeth, iatol, rtol, atol, iout, rout, ipar, rpar, ier)
    call fcvdense(NEQ, ier)

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

                call fcvode_wrapper(half_dt,rho,T_orig,ne_orig,e_orig, &
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

    call fcvfree

end subroutine integrate_state_fcvode

subroutine fcvode_wrapper(dt, rho_in, T_in, ne_in, e_in, T_out, ne_out, e_out)

    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    real(rt), intent(  out) ::         T_out,ne_out,e_out

    ! Set the number of independent variables -- this should be just "e"
    integer*8 :: NEQ = 1, ipar, iout(25)
  
    ! Allocate storage for the input state
    double precision :: y(1)

    double precision :: atol(1), rtol(1)
    double precision   :: time
    
    integer :: ier

    EXTERNAL jac, f_rhs

    double precision :: t_soln

    double precision :: rout(10), rpar
    integer :: meth, itmeth, iatol

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" in time. 
    y(1) = e_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    rtol(1) = 1.d-4

    call fcvreinit(time, y, 1, rtol(1), atol(1), ier)

    call fcvode(dt, t_soln, y, 1, ier)

    e_out  = y(1)
    T_out  = T_vode
    ne_out = ne_vode

end subroutine fcvode_wrapper
