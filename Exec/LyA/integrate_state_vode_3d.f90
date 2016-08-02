subroutine integrate_state_vode(lo, hi, &
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
    double precision, intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    double precision, intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    double precision, intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    double precision :: z, rho
    double precision :: T_orig, ne_orig, e_orig
    double precision :: T_out , ne_out , e_out, mu, mean_rhob

    z = 1.d0/a - 1.d0

    z_vode = z
    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Interpolate from the table to this redshift
    call interp_to_this_z(z)

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    !$OMP PARALLEL DO PRIVATE(i,j,k,rho,e_orig,T_orig,ne_orig,T_out,ne_out,e_out)
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

                call vode_wrapper(half_dt,rho,T_orig,ne_orig,e_orig, &
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
    !$OMP END PARALLEL DO

end subroutine integrate_state_vode

subroutine vode_wrapper(dt, rho_in, T_in, ne_in, e_in, T_out, ne_out, e_out)

    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode

    implicit none

    double precision, intent(in   ) :: dt
    double precision, intent(in   ) :: rho_in, T_in, ne_in, e_in
    double precision, intent(  out) ::         T_out,ne_out,e_out

    ! Set the number of independent variables -- this should be just "e"
    integer, parameter :: NEQ = 1
  
    ! Allocate storage for the input state
    double precision :: y(NEQ)

    ! Our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

    ! Tolerance parameters:
    !
    !  itol specifies whether to use an single absolute tolerance for
    !  all variables (1), or to pass an array of absolute tolerances, one
    !  for each variable with a scalar relative tol (2), a scalar absolute
    !  and array of relative tolerances (3), or arrays for both (4)
    !  
    !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
    !  be > 0.  
    !
    ! We will use arrays for both the absolute and relative tolerances, 
    ! since we want to be easier on the temperature than the species

    integer, parameter :: ITOL = 1
    double precision :: atol(NEQ), rtol(NEQ)
    
    ! We want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1
  
    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate

    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ

    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    double precision   :: rwork(LRW)
    double precision   :: time
    ! double precision   :: dt4
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    double precision :: rpar
    integer          :: ipar

    EXTERNAL jac, f_rhs
    
    logical, save :: firstCall = .true.

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in

    ! We want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = 0.d0
    iwork(:) = 0
    
    ! Set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 2000
    
    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" in time. 
    y(1) = e_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    rtol(1) = 1.d-4

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    e_out  = y(1)
    T_out  = T_vode
    ne_out = ne_vode

    if (istate < 0) then
       print *, 'istate = ', istate, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call bl_error("ERROR in vode_wrapper: integration failed")
    endif

!      print *,'Calling vode with 1/4 the time step'
!      dt4 = 0.25d0  * dt
!      y(1) = e_in

!      do n = 1,4
!         call dvode(f_rhs, NEQ, y, time, dt4, ITOL, rtol, atol, ITASK, &
!                    istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
!                    rpar, ipar)
!         if (istate < 0) then
!            print *, 'doing subiteration ',n
!            print *, 'istate = ', istate, 'at (i,j,k) ',i,j,k
!            call bl_error("ERROR in vode_wrapper: sub-integration failed")
!         end if

!      end do
!   endif

end subroutine vode_wrapper
