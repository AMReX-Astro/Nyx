#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <f_rhs.H>
#include <f_rhs_struct.H>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP   */

#include <nvector/nvector_serial.h>
#ifdef _OPENMP
#include <nvector/nvector_openmp.h>
#endif
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif
//#define MAKE_MANAGED 1
using namespace amrex;

RhsData* f_rhs_data;
/*
amrex::Vector<void*> ptr_lst;
//static amrex::Arena* Managed_Arena;

void* sunalloc(size_t mem_size)
{
  amrex::MultiFab::updateMemUsage ("Sunalloc", mem_size, nullptr);
  amrex::MultiFab::updateMemUsage ("All", mem_size, nullptr);
  void * ptr = (void*) The_Arena()->alloc(mem_size);
  ptr_lst.push_back(ptr);
  return ptr;
}

void sunfree(void* ptr)
{
  size_t mem_size = dynamic_cast<CArena*>(The_Arena())->sizeOf(ptr);
  ptr_lst.erase(std::remove_if(ptr_lst.begin(), ptr_lst.end(), [ptr](void* x) { return x == ptr; }));
  The_Arena()->free(ptr);
  amrex::MultiFab::updateMemUsage ("Sunalloc", -mem_size, nullptr);
  amrex::MultiFab::updateMemUsage ("All", -mem_size, nullptr);
}
*/
/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
/*
subroutine integrate_state_with_source(lo, hi, &
                                state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                state_n ,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                hydro_src, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                reset_src,srcr_l1,srcr_l2,srcr_l3,srcr_h1,srcr_h2,srcr_h3, &
                                I_R, ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3, &
                                a, delta_time, min_iter, max_iter) &
                                bind(C, name="integrate_state_with_source")
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
!   hydro_src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   reset_src_* : doubles arrays
!       The source terms based on the reset correction
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   delta_time : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_error_module, only : amrex_abort
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, gamma_minus_1
    use amrex_constants_module, only: M_PI, ONE, HALF
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: JH_vode, JHe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   T_zhi, T_zheii, inhomogeneous_on

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: sn_l1, sn_l2, sn_l3, sn_h1, sn_h2, sn_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    integer         , intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer         , intent(in) :: srcr_l1, srcr_l2, srcr_l3, srcr_h1, srcr_h2, srcr_h3
    integer         , intent(in) :: ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3
    real(rt), intent(in   ) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) ::  state_n(sn_l1:sn_h1, sn_l2:sn_h2,sn_l3:sn_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in   ) :: hydro_src(src_l1:src_h1, src_l2:src_h2,src_l3:src_h3, NVAR)
    real(rt), intent(in   ) :: reset_src(srcr_l1:srcr_h1, srcr_l2:srcr_h2,srcr_l3:srcr_h3, 1)
    real(rt), intent(inout) :: I_R(ir_l1:ir_h1, ir_l2:ir_h2,ir_l3:ir_h3)
    real(rt), intent(in)    :: a, delta_time
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    real(rt) :: asq,aendsq,ahalf,ahalf_inv,delta_rho,delta_e,delta_rhoe
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
    real(rt) :: rho_orig, T_orig, ne_orig, e_orig
    real(rt) :: rho_out, T_out, ne_out, e_out
    real(rt) :: rho_src, rhoe_src, e_src
    real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)

    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, delta_time)
    z_end = 1.0d0/a_end - 1.0d0
 
    asq = a*a
    aendsq = a_end*a_end
    ahalf     = HALF * (a + a_end)
    ahalf_inv  = ONE / ahalf

    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)

                ! Original values
                rho_orig  = state(i,j,k,URHO)
                e_orig    = state(i,j,k,UEINT) / rho_orig
                T_orig    = diag_eos(i,j,k,TEMP_COMP)
                ne_orig   = diag_eos(i,j,k,  NE_COMP)

                if (inhomogeneous_on) then
                   H_reion_z = diag_eos(i,j,k,ZHI_COMP)
                   if (z .gt. H_reion_z) then
                      JH_vode = 0
                   else
                      JH_vode = 1
                   endif
                endif

                if (e_orig .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e entering strang integration ', z, i,j,k, rho_orig/mean_rhob, e_orig
                    call amrex_abort('bad e in strang')
                    !$OMP END CRITICAL
                end if

                rho_src  = hydro_src(i,j,k,URHO)
                rhoe_src = hydro_src(i,j,k,UEINT)

                !                             
                ! This term satisfies the equation anewsq * (rho_new e_new) = 
                !                                  aoldsq * (rho_old e_old) + dt * H_{rho e} + anewsq * (reset_src)
                !  where e_new = e_old + dt * e_src                           
                e_src = ( ((asq*state(i,j,k,UEINT) + delta_time * rhoe_src ) / aendsq )/ &
                        (      state(i,j,k,URHO ) + delta_time * rho_src  ) - e_orig) / delta_time
                e_src = ( ((asq*state(i,j,k,UEINT) + delta_time * rhoe_src ) / aendsq + reset_src(i,j,k,1))/ &
                        (      state(i,j,k,URHO ) + delta_time * rho_src  ) - e_orig) / delta_time

                i_vode = i
                j_vode = j
                k_vode = k


!                call vode_wrapper_with_source_single(delta_time,rho_orig,T_orig,ne_orig,e_orig,rho_src,e_src, &
!                                                         rho_out ,T_out ,ne_out ,e_out)
                
                call vode_wrapper_with_source(delta_time,rho_orig,T_orig,ne_orig,e_orig,rho_src,e_src, &
                                                         rho_out ,T_out ,ne_out ,e_out)
                !                             
                ! I_R satisfies the equation anewsq * (rho_out  e_out ) = 
                !                            aoldsq * (rho_orig e_orig) + dt * a_half * I_R + dt * H_{rho e} + anewsq * reset_src
                I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*rhoe_src) ) ) / &
                              (delta_time * ahalf) - aendsq * reset_src(i,j,k,1) / (delta_time * ahalf)
                
                if ((state_n(i,j,k,UEINT) + delta_time * ahalf * I_R(i,j,k) / aendsq) / state_n(i,j,k,URHO) .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e exiting single step integration ', z, i,j,k, rho_orig/mean_rhob, e_out, &
                             state_n(i,j,k,UEINT)/state_n(i,j,k,URHO), (state_n(i,j,k,UEINT) + delta_time * ahalf * I_R(i,j,k) / aendsq) / state_n(i,j,k,URHO)/state_n(i,j,k,URHO)
                    !call flush(6)
                    !$OMP END CRITICAL

                    !!!!! FIXME !!!!!!
!                    call vode_wrapper_with_source(delta_time,rho_orig,T_orig,ne_orig,e_orig,rho_src,e_src, &
!                         rho_out ,T_out ,ne_out ,e_out)

                    T_out  = 10.0
                    ne_out = 0.0
                    mu     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)
                    e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                    
!                    I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*rhoe_src) ) ) / &
!                         (delta_time * ahalf) - aendsq * reset_src(i,j,k,1) / (delta_time * ahalf)
                    !call amrex_abort('bad e out of strang')
                end if


                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho_out, e_out, a, species)

                ! Instanteneous heating from reionization:
                T_H = 0.0d0
                if (inhomogeneous_on .or. flash_h) then
                   if ((H_reion_z  .lt. z) .and. (H_reion_z  .ge. z_end)) T_H  = (1.0d0 - species(2))*max((T_zhi-T_out), 0.0d0)
                endif

                T_He = 0.0d0
                if (flash_he) then
                   if ((He_reion_z .lt. z) .and. (He_reion_z .ge. z_end)) T_He = (1.0d0 - species(5))*max((T_zheii-T_out), 0.0d0)
                endif

                !!!!!!!!! FIXME !!!!!!!!!!!!!!!
                if ((T_H .gt. 0.0d0) .or. (T_He .gt. 0.0d0)) then
                   T_out = T_out + T_H + T_He                            ! For simplicity, we assume
                   ne_out = 1.0d0 + YHELIUM                              !    completely ionized medium at
                   if (T_He .gt. 0.0d0) ne_out = ne_out + YHELIUM        !    this point.  It's a very minor
                   mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)   !    detail compared to the overall approximation.
                   e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                   I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*rhoe_src) ) ) / &
                         (delta_time * ahalf) - aendsq * reset_src(i,j,k,1) / (delta_time * ahalf)
                   call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho_out, e_out, a, species)
                endif

                ! Update (rho e) and (rho E)
                ! Note that we add to state_n because those already have hydro_source in them
                state_n(i,j,k,UEINT) = state_n(i,j,k,UEINT) + delta_time * ahalf * I_R(i,j,k) / aendsq
                state_n(i,j,k,UEDEN) = state_n(i,j,k,UEDEN) + delta_time * ahalf * I_R(i,j,k) / aendsq
                 
                ! Update T and ne
                diag_eos(i,j,k,TEMP_COMP) = T_out
                diag_eos(i,j,k,  NE_COMP) = ne_out

            end do ! i
        end do ! j
    end do ! k

end subroutine integrate_state_with_source

subroutine vode_wrapper_with_source_single(dt, rho_in, T_in, ne_in, e_in, rho_src, e_src, rho_out, T_out, ne_out, e_out)

    use amrex_error_module, only : amrex_abort
    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode, NR_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in,  T_in, ne_in, e_in
    real(rt), intent(  out) :: rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    ! Set the number of independent variables -- this should be just "e"
    integer, parameter :: NEQ = 2
  
    ! Allocate storage for the input state
    real(rt) :: y(NEQ)

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

    integer, parameter :: ITOL = 2
    real(rt) :: atol(NEQ), rtol(NEQ)
    
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
    real(rt)   :: rwork(LRW)
    real(rt)   :: time
    ! real(rt)   :: dt4
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    real(rt) :: rpar
    integer          :: ipar

    EXTERNAL jac, f_rhs, f_rhs_split
    
    logical, save :: firstCall = .true.

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    NR_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! We want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = 0.d0
    iwork(:) = 0
    
    ! Set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 2000

    ! Set the initial H0 value
    rwork(5) = dt

    ! Set the maximum order value
    iwork(5) = 1

    ! Set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 1
    
    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    y(1) = e_in
    y(2) = rho_in


    !!!!!!!!!!!! Consider loosening tolerances !!!!!!!!!!!!!!!!!
    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    rtol(1) = 1.d-4

    atol(2) = 1.d-4 * rho_in
    rtol(2) = 1.d-4

    ! call the integration routine
    call dvode(f_rhs_split, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    e_out  = y(1)
    rho_out = y(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    T_out  = T_vode
    ne_out = ne_vode

    if (istate < 0) then
       print *, 'istate = ', istate, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call vode_wrapper_with_source(dt,rho_in,T_in,ne_in,e_in,rho_src,e_src, &
            rho_out ,T_out ,ne_out ,e_out)
!       call amrex_error("ERROR in vode_wrapper: integration failed")
    endif

  end subroutine vode_wrapper_with_source_single

  subroutine vode_wrapper_with_source(dt, rho_in, T_in, ne_in, e_in, rho_src, e_src, rho_out, T_out, ne_out, e_out)

    use amrex_error_module, only : amrex_error
    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode, NR_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in,  T_in, ne_in, e_in
    real(rt), intent(  out) :: rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    ! Set the number of independent variables -- this should be just "e"
    integer, parameter :: NEQ = 2
  
    ! Allocate storage for the input state
    real(rt) :: y(NEQ)

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

    integer, parameter :: ITOL = 2
    real(rt) :: atol(NEQ), rtol(NEQ)
    
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
    real(rt)   :: rwork(LRW)
    real(rt)   :: time
    ! real(rt)   :: dt4
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    real(rt) :: rpar
    integer          :: ipar

    EXTERNAL jac, f_rhs, f_rhs_split
    
    logical, save :: firstCall = .true.

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    NR_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! We want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = 0.d0
    iwork(:) = 0
    
    ! Set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 2000

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    y(1) = e_in
    y(2) = rho_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    rtol(1) = 1.d-4

    atol(2) = 1.d-4 * rho_in
    rtol(2) = 1.d-4

    ! call the integration routine
    call dvode(f_rhs_split, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    e_out  = y(1)
    rho_out = y(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    T_out  = T_vode
    ne_out = ne_vode

    if (istate < 0) then
       print *, 'istate = ', istate, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call amrex_error("ERROR in vode_wrapper: integration failed")
    endif

end subroutine vode_wrapper_with_source
*/
int Nyx::integrate_state_struct
  (amrex::MultiFab &S_old,
   amrex::MultiFab &S_new,
   amrex::MultiFab &D_old,
   amrex::MultiFab &hydro_src,
   amrex::MultiFab &IR,
   amrex::MultiFab &reset_src,
   const Real& a, const amrex::Real& a_end,
   const Real& delta_time,
   const int sdc_iter)
{
    // time = starting time in the simulation
    f_rhs_data = (RhsData*)The_Arena()->alloc(sizeof(RhsData));
    Real gamma_minus_1 = gamma - 1.0;
    ode_eos_setup(f_rhs_data, gamma_minus_1, h_species);

  amrex::Gpu::LaunchSafeGuard lsg(true);
  long int store_steps=new_max_sundials_steps;

  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
  for ( MFIter mfi(S_old, false); mfi.isValid(); ++mfi )
#else
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
#endif
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();

      S_old[mfi].prefetchToDevice();
      D_old[mfi].prefetchToDevice();
      S_new[mfi].prefetchToDevice();
      hydro_src[mfi].prefetchToDevice();
      reset_src[mfi].prefetchToDevice();
      IR[mfi].prefetchToDevice();
      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);
      Array4<Real> const& state_n4 = S_new.array(mfi);
      Array4<Real> const& hydro_src4 = hydro_src.array(mfi);
      Array4<Real> const& reset_src4 = reset_src.array(mfi);
      Array4<Real> const& IR4 = IR.array(mfi);

      integrate_state_struct_mfin(state4,diag_eos4,state_n4,hydro_src4,reset_src4,IR4,tbx,a,a_end,delta_time,store_steps,new_max_sundials_steps,sdc_iter);
    }
      return 0;
}

int Nyx::integrate_state_struct_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   amrex::Array4<Real> const& state_n4,
   amrex::Array4<Real> const& hydro_src4,
   amrex::Array4<Real> const& reset_src4,
   amrex::Array4<Real> const& IR4,
   const Box& tbx,
   const Real& a, const amrex::Real& a_end, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter)
{

  realtype reltol, abstol;
  int flag;
    
  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  int one_in = 1;
  

      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box

      N_Vector u;
      N_Vector e_orig;
      N_Vector Data;
      N_Vector abstol_vec;
      N_Vector T_vec;
      N_Vector ne_vec;
      N_Vector rho_vec;
      N_Vector rho_init_vec;
      N_Vector rho_src_vec;
      N_Vector rhoe_src_vec;
      N_Vector e_src_vec;
      N_Vector IR_vec;

      void *cvode_mem;
      double *dptr, *eptr, *rpar, *rparh, *abstol_ptr;
      realtype t=0.0;
                                
      u = NULL;
      e_orig = NULL;
      Data = NULL;
      abstol_vec = NULL;
      cvode_mem = NULL;

      long int neq = len.x*len.y*len.z;
      amrex::Gpu::streamSynchronize();
      int loop = 1;

#ifdef AMREX_USE_CUDA
                        cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();  
                        if(sundials_alloc_type%2==0)
                        {
                          if(sundials_alloc_type==0)
                            u = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
                          else
                            u = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          dptr=N_VGetDeviceArrayPointer_Cuda(u);

                          if(sundials_alloc_type==0)
                            e_orig = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
                          else
                            e_orig = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
                          N_VSetCudaStream_Cuda(e_orig, &currentStream);
                          N_VSetCudaStream_Cuda(u, &currentStream);

                          if(sundials_alloc_type==0)
                            Data = N_VMakeWithManagedAllocator_Cuda(4*neq,sunalloc,sunfree);  // Allocate u vector 
                          else
                            Data = N_VNewManaged_Cuda(4*neq);  /* Allocate u vector */

                          rparh = N_VGetDeviceArrayPointer_Cuda(Data);
                          N_VSetCudaStream_Cuda(Data, &currentStream);
                          // shouldn't need to initialize 
                          //N_VConst(0.0,Data);

                          if(sundials_alloc_type==0)
                            abstol_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  
                          else
                            abstol_vec = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          abstol_ptr = N_VGetDeviceArrayPointer_Cuda(abstol_vec);
                          N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
                          amrex::Gpu::Device::streamSynchronize();
                        }
                        else
                        {
                          dptr=(double*) The_Managed_Arena()->alloc(neq*sizeof(double));
                          u = N_VMakeManaged_Cuda(neq,dptr);  /* Allocate u vector */
                          eptr= (double*) The_Managed_Arena()->alloc(neq*sizeof(double));
                          e_orig = N_VMakeManaged_Cuda(neq,eptr);  /* Allocate u vector */
                          N_VSetCudaStream_Cuda(e_orig, &currentStream);
                          N_VSetCudaStream_Cuda(u, &currentStream);

                          rparh = (double*) The_Managed_Arena()->alloc(4*neq*sizeof(double));
                          Data = N_VMakeManaged_Cuda(4*neq,rparh);  // Allocate u vector 
                          N_VSetCudaStream_Cuda(Data, &currentStream);
                          // shouldn't need to initialize 
                          //N_VConst(0.0,Data);

                          abstol_ptr = (double*) The_Managed_Arena()->alloc(neq*sizeof(double));
                          abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr);
                          N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
                          amrex::Gpu::streamSynchronize();
                        }

#else
#ifdef _OPENMP
                        int nthreads=omp_get_max_threads();
                        u = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
                        e_orig = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_OpenMP(u);

                        Data = N_VNew_OpenMP(4*neq,nthreads);  // Allocate u vector 
                        N_VConst(0.0,Data);
                        rparh=N_VGetArrayPointer_OpenMP(Data);
                        abstol_vec = N_VNew_OpenMP(neq,nthreads);
                        abstol_ptr=N_VGetArrayPointer_OpenMP(abstol_vec);
#else
                        u = N_VNew_Serial(neq);  /* Allocate u vector */
                        e_orig = N_VNew_Serial(neq);  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_Serial(u);

                        Data = N_VNew_Serial(4*neq);  // Allocate u vector 
                        N_VConst(0.0,Data);
                        rparh=N_VGetArrayPointer_Serial(Data);
                        abstol_vec = N_VNew_Serial(neq);
                        abstol_ptr=N_VGetArrayPointer_Serial(abstol_vec);
                        T_vec = N_VNew_Serial(neq);
                        ne_vec = N_VNew_Serial(neq);
                        rho_vec = N_VNew_Serial(neq);
                        rho_init_vec = N_VNew_Serial(neq);
                        rho_src_vec = N_VNew_Serial(neq);
                        rhoe_src_vec = N_VNew_Serial(neq);
                        e_src_vec = N_VNew_Serial(neq);
                        IR_vec = N_VNew_Serial(neq);
                        amrex::Real* T_vode=N_VGetArrayPointer_Serial(T_vec);
                        amrex::Real* ne_vode=N_VGetArrayPointer_Serial(ne_vec);
                        amrex::Real* rho_vode=N_VGetArrayPointer_Serial(rho_vec);
                        amrex::Real* rho_init_vode=N_VGetArrayPointer_Serial(rho_init_vec);
                        amrex::Real* rho_src_vode=N_VGetArrayPointer_Serial(rho_src_vec);
                        amrex::Real* rhoe_src_vode=N_VGetArrayPointer_Serial(rhoe_src_vec);
                        amrex::Real* e_src_vode=N_VGetArrayPointer_Serial(e_src_vec);
                        amrex::Real* IR_vode=N_VGetArrayPointer_Serial(IR_vec);
#endif
#endif
    ode_eos_initialize_single(f_rhs_data, a, dptr, eptr, T_vode, ne_vode, rho_vode, rho_init_vode, rho_src_vode, rhoe_src_vode, e_src_vode, IR_vode);

#ifdef _OPENMP
      const Dim3 hi = amrex::ubound(tbx);
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                //Skip setup since parameters are hard-coded
#else
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {                                 
#endif
                                  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
                                  ode_eos_initialize_arrays(i, j, k, idx, f_rhs_data,
                                                            a_end, lo, len, state4, diag_eos4,
                                                            hydro_src4, reset_src4, dptr, eptr,
                                                            abstol_ptr, sdc_iter, delta_time);

                                  dptr[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
                                  eptr[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
                                  rparh[4*idx+0]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
                                  rparh[4*idx+1]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
                                  rparh[4*idx+2]= state4(i,j,k,Density); //    rpar(3)=rho_vode
                                  rparh[4*idx+3]=1/a-1;    //    rpar(4)=z_vode
                                  abstol_ptr[idx]= state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol;
                                  //                            }
#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
      amrex::Gpu::Device::streamSynchronize();
#endif
                        
#ifdef CV_NEWTON
                                cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
                                cvode_mem = CVodeCreate(CV_BDF);
#endif
                                flag = CVodeInit(cvode_mem, f, t, u);

                                N_VScale(abstol,u,abstol_vec);
                                //                              N_VConst(N_VMin(abstol_vec),abstol_vec);

                                flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

                                //                              flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
                                flag = CVDiag(cvode_mem);

                                CVodeSetMaxNumSteps(cvode_mem,2000);

                                if(use_typical_steps)
                                    CVodeSetMaxStep(cvode_mem,delta_time/(old_max_steps));

                                N_Vector constrain;
                                if(use_sundials_constraint)
                                  {
                                    constrain=N_VClone(u);
                                    N_VConst(2,constrain);            
                                    flag =CVodeSetConstraints(cvode_mem,constrain);
                                  }

#ifdef SUNDIALS_VERSION_MAJOR
#if SUNDIALS_VERSION_MAJOR >= 5
#if SUNDIALS_VERSION_MINOR >= 3
                                if(use_sundials_fused)
                                {
                                     flag = CVodeSetUseIntegratorFusedKernels(cvode_mem, SUNTRUE);
                                }
#endif
#endif
#endif
                                CVodeSetUserData(cvode_mem, &Data);
                                //                              CVodeSetMaxStep(cvode_mem, delta_time/10);
                                //                              BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
                                flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
                                if(use_typical_steps)
                                  {
                                    long int nst=0;
                                    flag = CVodeGetNumSteps(cvode_mem, &nst);
                                    new_max_steps=std::max(nst,new_max_steps);
                                  }
                                //                              amrex::Gpu::Device::streamSynchronize();
                                //                              BL_PROFILE_VAR_STOP(cvode_timer2);

#ifdef AMREX_DEBUG
                                PrintFinalStats(cvode_mem);
#endif

#ifdef _OPENMP
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
#else
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {                                 
#endif
                                  int  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
                                //                              for (int i= 0;i < neq; ++i) {
                                  ode_eos_finalize((dptr[idx*loop]), &(rparh[4*idx*loop]), one_in);
                                  diag_eos4(i,j,k,Temp_comp)=rparh[4*idx*loop+0];   //rpar(1)=T_vode
                                  diag_eos4(i,j,k,Ne_comp)=rparh[4*idx*loop+1];//    rpar(2)=ne_vode
                                
                                  state4(i,j,k,Eint)  += state4(i,j,k,Density) * (dptr[idx*loop]-eptr[idx]);
                                  state4(i,j,k,Eden)  += state4(i,j,k,Density) * (dptr[idx*loop]-eptr[idx]);
                                  //                            }
                                //PrintFinalStats(cvode_mem);
#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
      amrex::Gpu::Device::streamSynchronize();
#endif


#ifdef AMREX_USE_CUDA
      if(sundials_alloc_type%2!=0)
      {
        The_Managed_Arena()->free(dptr);
        The_Managed_Arena()->free(eptr);
        if(use_sundials_constraint)
          The_Managed_Arena()->free(constrain);
        The_Managed_Arena()->free(rparh);
        The_Managed_Arena()->free(abstol_ptr);
      }
#endif

                                N_VDestroy(u);          /* Free the u vector */
                                N_VDestroy(e_orig);          /* Free the e_orig vector */
                                if(use_sundials_constraint)
                                  N_VDestroy(constrain);          /* Free the constrain vector */
                                N_VDestroy(abstol_vec);          /* Free the u vector */
                                N_VDestroy(Data);          /* Free the userdata vector */
                                CVodeFree(&cvode_mem);  /* Free the integrator memory */
                              //);
                                /*                          }
                        
        }
        }*/
                                return 0;
}

#ifdef AMREX_USE_CUDA
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  Real* udot_ptr=N_VGetDeviceArrayPointer_Cuda(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer_Cuda(u);
  int neq=N_VGetLength_Cuda(udot);
  double*  rpar=N_VGetDeviceArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
  
  cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
  AMREX_LAUNCH_DEVICE_LAMBDA ( neq, idx, {
      //  f_rhs_test(t,u_ptr,udot_ptr, rpar, neq);
        f_rhs_rpar(t, (u_ptr[idx]),(udot_ptr[idx]),&(rpar[4*idx]));
  });
  cudaStreamSynchronize(currentStream);

  return 0;
}

#else
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{

  Real* udot_ptr=N_VGetArrayPointer_Serial(udot);
  Real* u_ptr=N_VGetArrayPointer_Serial(u);
  int neq=N_VGetLength_Serial(udot);
  double*  rpar=N_VGetArrayPointer_Serial(*(static_cast<N_Vector*>(user_data)));

  #pragma omp parallel for
  for(int tid=0;tid<neq;tid++)
    {
        f_rhs_rpar(t, (u_ptr[tid]),(udot_ptr[tid]),&(rpar[4*tid]));
    }

  return 0;
}
#endif
