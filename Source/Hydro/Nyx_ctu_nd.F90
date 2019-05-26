! advection routines in support of the CTU unsplit advection scheme

module ctu_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains


  !> @brief Compute the normal interface states by reconstructing
  !! the primitive variables using the piecewise parabolic method
  !! and doing characteristic tracing.  We do not apply the
  !! transverse terms here.
  !!
  !! @param[in] q            (const)  input state, primitives
  !! @param[in] qaux         (const)  auxiliary hydro data
  !! @param[in] flatn        (const)  flattening parameter
  !! @param[in] srcQ         (const)  primitive variable source
  !! @param[in] dx           (const)  grid spacing in X, Y, Z direction
  !! @param[in] dt           (const)  time stepsize
  !! @param[inout] flux1        (modify) flux in X direction on X edges
  !! @param[inout] flux2        (modify) flux in Y direction on Y edges
  !! @param[inout] flux3        (modify) flux in Z direction on Z edges
  !! @param[inout] q1           (modify) Godunov interface state in X
  !! @param[inout] q2           (modify) Godunov interface state in Y
  !! @param[inout] q3           (modify) Godunov interface state in Z
  !!
  subroutine ctu_ppm_states(lo, hi, &
                            vlo, vhi, &
                            q, qd_lo, qd_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            srcQ, src_lo, src_hi, &
                            shk, sk_lo, sk_hi, &
                            Ip, Ip_lo, Ip_hi, &
                            Im, Im_lo, Im_hi, &
                            Ip_src, Ips_lo, Ips_hi, &
                            Im_src, Ims_lo, Ims_hi, &
                            Ip_gc, Ipg_lo, Ipg_hi, &
                            Im_gc, Img_lo, Img_hi, &
                            sm, sm_lo, sm_hi, &
                            sp, sp_lo, sp_hi, &
                            qxm, qxm_lo, qxm_hi, &
                            qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM >= 2
                            qym, qym_lo, qym_hi, &
                            qyp, qyp_lo, qyp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            qzm, qzm_lo, qzm_hi, &
                            qzp, qzp_lo, qzp_hi, &
#endif
                            dx, dt, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_ppm_states")

    use meth_params_module, only : NQSRC, QVAR, NVAR, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   ppm_predict_gammae, ppm_temp_fix, &
                                   hybrid_riemann
    use ca_ppm_module, only : ca_ppm_reconstruct, ppm_int_profile, ppm_reconstruct_with_eos
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use ca_trace_ppm_module, only : trace_ppm
#endif
    use advection_module, only : ca_shock
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
    integer, intent(in) :: sm_lo(3), sm_hi(3)
    integer, intent(in) :: sp_lo(3), sp_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)

    real(rt), intent(inout) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:3,QVAR)
    real(rt), intent(inout) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:3,QVAR)
    real(rt), intent(inout) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:3,NQSRC)
    real(rt), intent(inout) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:3,NQSRC)
    real(rt), intent(inout) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(inout) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1), sm_lo(2):sm_hi(2), sm_lo(3):sm_hi(3))
    real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1), sp_lo(2):sp_hi(2), sp_lo(3):sp_hi(3))

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), QVAR)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), QVAR)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), QVAR)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), QVAR)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), QVAR)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), QVAR)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n, idir

    logical :: source_nonzero(NQSRC)
    logical :: reconstruct_state(QVAR)

    logical :: compute_shock

    !$gpu

    hdt = HALF*dt

    ! multidimensional shock detection

#ifdef SHOCK_VAR
    compute_shock = .true.
#else
    compute_shock = .false.
#endif

    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1 .or. compute_shock) then
       call ca_shock(lo, hi, &
                     q, qd_lo, qd_hi, &
                     shk, sk_lo, sk_hi, &
                     dx)
    else
       shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
    endif

    ! we don't need to reconstruct all of the NQ state variables,
    ! depending on how we are tracing
    reconstruct_state(:) = .true.
    reconstruct_state(QREINT) = .false.

    ! preprocess the sources -- we don't want to trace under a source
    ! that is empty.  Note, we need to do this check over the entire
    ! grid, to be sure, e.g., use vlo:vhi. On the GPU, this check is
    ! expensive and for now we just disable this optimization and eat
    ! the cost of processing the sources, even if they're zero.

#ifdef AMREX_USE_CUDA
    source_nonzero(:) = .true.
#else
    do n = 1, NQSRC
       if (minval(srcQ(vlo(1)-2:vhi(1)+2,vlo(2)-2*dg(2):vhi(2)+2*dg(2),vlo(3)-2*dg(3):vhi(3)+2*dg(3),n)) == ZERO .and. &
           maxval(srcQ(vlo(1)-2:vhi(1)+2,vlo(2)-2*dg(2):vhi(2)+2*dg(2),vlo(3)-2*dg(3):vhi(3)+2*dg(3),n)) == ZERO) then
          source_nonzero(n) = .false.
       else
          source_nonzero(n) = .true.
       endif
    enddo
#endif


    do idir = 1, AMREX_SPACEDIM

       ! Compute Ip and Im -- this does the parabolic reconstruction,
       ! limiting, and returns the integral of each profile under each
       ! wave to each interface
       do n = 1, QVAR
          if (.not. reconstruct_state(n)) cycle

          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  q, qd_lo, qd_hi, QVAR, n, n, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               q, qd_lo, qd_hi, QVAR, n, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip, Ip_lo, Ip_hi, &
                               Im, Im_lo, Im_hi, QVAR, n, &
                               dx, dt)
       end do


       if (ppm_temp_fix /= 1) then
          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  qaux, qa_lo, qa_hi, NQAUX, QGAMC, QGAMC, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_gc, Ipg_lo, Ipg_hi, &
                               Im_gc, Img_lo, Img_hi, 1, 1, &
                               dx, dt)
       else

          ! temperature-based PPM
          call ppm_reconstruct_with_eos(lo, hi, idir, &
                                        Ip, Ip_lo, Ip_hi, &
                                        Im, Im_lo, Im_hi, &
                                        Ip_gc, Ipg_lo, Ipg_hi, &
                                        Im_gc, Img_lo, Img_hi)

       end if


       ! source terms
       do n = 1, NQSRC
          if (source_nonzero(n)) then
             call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                     srcQ, src_lo, src_hi, NQSRC, n, n, &
                                     flatn, f_lo, f_hi, &
                                     sm, sm_lo, sm_hi, &
                                     sp, sp_lo, sp_hi, &
                                     1, 1, 1)

             call ppm_int_profile(lo, hi, idir, &
                                  srcQ, src_lo, src_hi, NQSRC, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  Ip_src, Ips_lo, Ips_hi, &
                                  Im_src, Ims_lo, Ims_hi, NQSRC, n, &
                                  dx, dt)
          else
             Ip_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,n) = ZERO
             Im_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,n) = ZERO
          endif

       enddo


       ! compute the interface states

#ifdef RADIATION
       if (idir == 1) then
          call trace_ppm_rad(lo, hi, &
                             1, q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip, Ip_lo, Ip_hi, &
                             Im, Im_lo, Im_hi, &
                             Ip_src, Ips_lo, Ips_hi, &
                             Im_src, Ims_lo, Ims_hi, &
                             qxm, qxm_lo, qxm_hi, &
                             qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm_rad(lo, hi, &
                             2, q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip, Ip_lo, Ip_hi, &
                             Im, Im_lo, Im_hi, &
                             Ip_src, Ips_lo, Ips_hi, &
                             Im_src, Ims_lo, Ims_hi, &
                             qym, qym_lo, qym_hi, &
                             qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm_rad(lo, hi, &
                             3, q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip, Ip_lo, Ip_hi, &
                             Im, Im_lo, Im_hi, &
                             Ip_src, Ips_lo, Ips_hi, &
                             Im_src, Ims_lo, Ims_hi, &
                             qzm, qzm_lo, qzm_hi, &
                             qzp, qzp_lo, qzp_hi, &
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif
       endif
#else
       ! hydro (no radiation)
       if (idir == 1) then
          call trace_ppm(lo, hi, &
                         1, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         flatn, f_lo, f_hi, &
                         Ip, Ip_lo, Ip_hi, &
                         Im, Im_lo, Im_hi, &
                         Ip_src, Ips_lo, Ips_hi, &
                         Im_src, Ims_lo, Ims_hi, &
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &
                         qxm, qxm_lo, qxm_hi, &
                         qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm(lo, hi, &
                         2, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         flatn, f_lo, f_hi, &
                         Ip, Ip_lo, Ip_hi, &
                         Im, Im_lo, Im_hi, &
                         Ip_src, Ips_lo, Ips_hi, &
                         Im_src, Ims_lo, Ims_hi, &
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &
                         qym, qym_lo, qym_hi, &
                         qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm(lo, hi, &
                         3, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         flatn, f_lo, f_hi, &
                         Ip, Ip_lo, Ip_hi, &
                         Im, Im_lo, Im_hi, &
                         Ip_src, Ips_lo, Ips_hi, &
                         Im_src, Ims_lo, Ims_hi, &
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &
                         qzm, qzm_lo, qzm_hi, &
                         qzp, qzp_lo, qzp_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif
       end if
#endif

    end do

  end subroutine ctu_ppm_states


  !> @brief Compute the normal interface states by reconstructing
  !! the primitive variables using piecewise linear slopes and doing
  !! characteristic tracing.  We do not apply the transverse terms here.
  !!
  !! @todo we can get rid of the the different temporary q Godunov
  !! state arrays
  !!
  !! @param[in] q            (const)  input state, primitives
  !! @param[in] qaux         (const)  auxiliary hydro data
  !! @param[in] flatn        (const)  flattening parameter
  !! @param[in] srcQ         (const)  primitive variable source
  !! @param[in] dx           (const)  grid spacing in X, Y, Z direction
  !! @param[in] dt           (const)  time stepsize
  !! @param[inout] flux1        (modify) flux in X direction on X edges
  !! @param[inout] flux2        (modify) flux in Y direction on Y edges
  !! @param[inout] flux3        (modify) flux in Z direction on Z edges
  !! @param[inout] q1           (modify) Godunov interface state in X
  !! @param[inout] q2           (modify) Godunov interface state in Y
  !! @param[inout] q3           (modify) Godunov interface state in Z
  !!
  AMREX_CUDA_FORT_DEVICE subroutine ctu_plm_states(lo, hi, &
                            idir, vlo, vhi, &
                            q, qd_lo, qd_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            srcQ, src_lo, src_hi, &
                            shk, sk_lo, sk_hi, &
                            dq, dq_lo, dq_hi, &
                            qm, qm_lo, qm_hi, &
                            qp, qp_lo, qp_hi, &
                            dx, dt, a_old, a_new, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_plm_states")

    use meth_params_module, only : NQSRC, QVAR, NVAR, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   iorder, use_pslope, hybrid_riemann
    use trace_plm_module, only : trace_plm, trace_plm_orig
    use ca_slope_module, only : uslope, pslope
#ifndef AMREX_USE_CUDA
    use advection_module, only : ca_shock
#endif
!    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)
    integer, intent(in), value :: idir
    real(rt), intent(in), value :: dt
    real(rt), intent(in), value :: a_old
    real(rt), intent(in), value :: a_new
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)

    real(rt), intent(inout) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: dq(dq_lo(1):dq_hi(1), dq_lo(2):dq_hi(2), dq_lo(3):dq_hi(3), QVAR)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), QVAR)
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n

    logical :: compute_shock

    !$gpu

    hdt = HALF*dt

    ! multidimensional shock detection
    shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
    call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif

       do n = 1, QVAR
          call uslope(lo, hi, idir, &
                      q, qd_lo, qd_hi, n, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi)
       end do

       if (use_pslope == 1) then
          call pslope(lo, hi, idir, &
                      q, qd_lo, qd_hi, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi, &
                      srcQ, src_lo, src_hi, &
                      dx)
       endif


       ! compute the interface states

          call trace_plm_orig(lo, hi, &
                         idir, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qm, qm_lo, qm_hi, &
                         qp, qp_lo, qp_hi, &
#if AMREX_SPACEDIM < 3
                         dloga, dloga_lo, dloga_hi, &
#endif
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt,a_old)

  end subroutine ctu_plm_states

#ifndef AMREX_USE_CUDA
  subroutine ctu_consup(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q, q_lo, q_hi, &
                        shk,  sk_lo, sk_hi, &
                        uout, uout_lo, uout_hi, &
                        update, updt_lo, updt_hi, &
                        flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                        flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                        Erin, Erin_lo, Erin_hi, &
                        Erout, Erout_lo, Erout_hi, &
                        radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                        radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        radflux3, radflux3_lo, radflux3_hi, &
#endif
                        nstep_fsp, &
#endif
                        qx, qx_lo, qx_hi, &
#if AMREX_SPACEDIM >= 2
                        qy, qy_lo, qy_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        qz, qz_lo, qz_hi, &
#endif
                        area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                        area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        area3, area3_lo, area3_hi, &
#endif
                        vol, vol_lo, vol_hi, &
                        pdivu, pdivu_lo, pdivu_hi, &
                        dx, dt) bind(C, name="ctu_consup")

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
#ifdef RADIATION
                                   fspace_type, comoving, &
                                   GDU, GDV, GDW, GDLAMS, GDERADS, &
#endif

                                   UEDEN, UEINT, UTEMP, NGDNV, QVAR
    use prob_params_module, only : dg
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
#endif
    use advection_module, only : calc_pdivu
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use amrex_constants_module, only : ZERO, ONE, TWO, FOURTH, HALF

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) ::  uout_lo(3),  uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
#endif
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)
    integer, intent(in) ::   pdivu_lo(3),   pdivu_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
    integer, intent(inout) :: nstep_fsp
#endif

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(in) :: shk(sk_lo(1):sk_hi(1),sk_lo(2):sk_hi(2),sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)

    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

#if AMREX_SPACEDIM >= 2
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
#endif

#if AMREX_SPACEDIM == 3
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
#endif

    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: pdivu(pdivu_lo(1):pdivu_hi(1),pdivu_lo(2):pdivu_hi(2),pdivu_lo(3):pdivu_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

#ifdef RADIATION
    real(rt), intent(in) :: Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    real(rt), intent(in) :: radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(in) :: radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(in) :: radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
#endif

#endif

    integer :: i, j, g, k, n
    integer :: domlo(3), domhi(3)
    real(rt) :: volInv

#ifdef RADIATION
    real(rt), dimension(0:ngroups-1) :: Erscale
    real(rt), dimension(0:ngroups-1) :: ustar, af
    real(rt) :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    real(rt) :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    real(rt) :: Gf1E(3)
    real(rt) :: ux, uy, uz, divu, lamc, Egdc
    real(rt) :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    real(rt) :: dprdx, dprdy, dprdz, ek1, ek2, dek, dpdx
    real(rt) :: urho_new
    real(rt) :: umx_new1, umy_new1, umz_new1
    real(rt) :: umx_new2, umy_new2, umz_new2
#endif

    !$gpu

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1(area1_lo(1),area1_lo(2),area1_lo(3)), &
#if AMREX_SPACEDIM >= 2
                    qy, qy_lo, qy_hi, &
                    area2(area2_lo(1),area2_lo(2),area2_lo(3)), &
#endif
#if AMREX_SPACEDIM == 3
                    qz, qz_lo, qz_hi, &
                    area3(area3_lo(1),area3_lo(2),area3_lo(3)), &
#endif
                    vol(vol_lo(1),vol_lo(2),vol_lo(3)), &
                    dx, pdivu, pdivu_lo, pdivu_hi)


    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
#if AMREX_SPACEDIM >= 2
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
#endif
#if AMREX_SPACEDIM == 3
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
#endif
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

#ifdef SHOCK_VAR
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k,USHK) = shk(i,j,k)
          end do
       end do
    end do
#endif

#ifdef HYBRID_MOMENTUM
    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, uout_lo, uout_hi, &
                                     qx, qx_lo, qx_hi, &
                                     qy, qy_lo, qy_hi, &
                                     qz, qz_lo, qz_hi)
#endif


  end subroutine ctu_consup
#endif
end module ctu_module

