subroutine integrate_state(lo, hi, &
                           state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                           diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                           dx, time, a, half_dt, min_iter, max_iter) &
                           bind(C, name="integrate_state")

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
!   diag_eos* : double arrays
!       Temp and Ne
!   src_* : double arrays
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
    use meth_params_module, only : NVAR
 
    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in   ) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    real(rt), intent(in   ) :: dx(3), time, a, half_dt
    integer         , intent(inout) :: min_iter, max_iter

    min_iter = 1
    max_iter = 1

    call integrate_state_force(lo, hi, state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                       diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                               dx, time, a, half_dt)

end subroutine integrate_state


! unused VODE stubs if we are not doing heating/cooling
module vode_aux_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt) :: z_vode
end module vode_aux_module
