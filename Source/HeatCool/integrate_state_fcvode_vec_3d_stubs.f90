subroutine integrate_state_fcvode_vec(lo, hi, &
                                  state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
!
    use amrex_error_module, only : amrex_abort
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, gamma_minus_1
    use bl_constants_module, only: M_PI
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: z_vode, i_vode, j_vode, k_vode, firstcall

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    real(rt) :: z, rho
    real(rt) :: T_orig, ne_orig, e_orig
    real(rt) :: T_out , ne_out , e_out, mu, mean_rhob

    call amrex_abort("Cannot call fcvode without compiling with USE_CVODE=TRUE")

end subroutine integrate_state_fcvode_vec
