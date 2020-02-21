! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module ca_trace_ppm_module

  use prob_params_module, only : dg
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : ZERO, HALF, ONE

  implicit none

  private

  public trace_ppm

contains

  subroutine trace_ppm(lo, hi, &
                       idir, q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       flatn, f_lo, f_hi, &
                       Ip, Ip_lo, Ip_hi, &
                       Im, Im_lo, Im_hi, &
                       Ip_src, Ips_lo, Ips_hi, &
                       Im_src, Ims_lo, Ims_hi, &
                       Ip_gc, Ipg_lo, Ipg_hi, &
                       Im_gc, Img_lo, Img_hi, &
                       qm, qm_lo, qm_hi, &
                       qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                       dloga, dloga_lo, dloga_hi, &
#endif
                       vlo, vhi, domlo, domhi, &
                       dx, dt)

    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NQAUX, NQSRC, ppm_predict_gammae, &
                            ppm_temp_fix, QU, QV, QW, npassive, qpass_map, &
                            fix_mass_flux, gamma_minus_1, &
                            ppm_flatten_before_integrals
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    
    real(rt), intent(in) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:3,QVAR)
    real(rt), intent(in) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:3,QVAR)

    real(rt), intent(in) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:3,NQSRC)
    real(rt), intent(in) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:3,NQSRC)

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)


    real(rt) :: un, xi
    integer :: ipassive, n, i, j, k

#if AMREX_SPACEDIM == 1
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    !$gpu

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (lo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (hi(1) == domhi(1))
#endif

    ! the passive stuff is the same regardless of the tracing
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  We will deal with all 3 velocity
       ! components below, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ppm_flatten_before_integrals == 0) then
                   xi = flatn(i,j,k)
                else
                   xi = ONE
                endif
                ! Plus state on face i
                if ((idir == 1 .and. i >= vlo(1)) .or. &
                    (idir == 2 .and. j >= vlo(2)) .or. &
                    (idir == 3 .and. k >= vlo(3))) then

                   un = q(i,j,k,QU-1+idir)

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1-wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   qp(i,j,k,n) = merge(q(i,j,k,n), xi*Im(i,j,k,2,n), un > ZERO)
                   if (n <= NQSRC) qp(i,j,k,n) = qp(i,j,k,n) + HALF*dt*xi*Im_src(i,j,k,2,n)

                end if

                ! Minus state on face i+1
                if (idir == 1 .and. i <= vhi(1)) then
                   un = q(i,j,k,QU-1+idir)
                   qm(i+1,j,k,n) = merge(xi*Ip(i,j,k,2,n), q(i,j,k,n), un > ZERO)
                   if (n <= NQSRC) qm(i+1,j,k,n) = qm(i+1,j,k,n) + HALF*dt*xi*Ip_src(i,j,k,2,n)
                else if (idir == 2 .and. j <= vhi(2)) then
                   un = q(i,j,k,QU-1+idir)
                   qm(i,j+1,k,n) = merge(xi*Ip(i,j,k,2,n), q(i,j,k,n), un > ZERO)
                   if (n <= NQSRC) qm(i,j+1,k,n) = qm(i,j+1,k,n) + HALF*dt*xi*Ip_src(i,j,k,2,n)
                else if (idir == 3 .and. k <= vhi(3)) then
                   un = q(i,j,k,QU-1+idir)
                   qm(i,j,k+1,n) = merge(xi*Ip(i,j,k,2,n), q(i,j,k,n), un > ZERO)
                   if (n <= NQSRC) qm(i,j,k+1,n) = qm(i,j,k+1,n) + HALF*dt*xi*Ip_src(i,j,k,2,n)
                end if

             end do

#if AMREX_SPACEDIM == 1
             if (fix_mass_flux_hi) qp(vhi(1)+1,j,k,n) = q(vhi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(vlo(1),j,k,n) = q(vlo(1)-1,j,k,n)
#endif
          end do
       end do

    end do

    if (ppm_temp_fix < 3) then
       if (ppm_predict_gammae == 0) then
          call trace_ppm_rhoe(lo, hi, &
                              idir, q, qd_lo, qd_hi, &
                              qaux, qa_lo, qa_hi, &
                              flatn, f_lo, f_hi, &
                              Ip, Ip_lo, Ip_hi, &
                              Im, Im_lo, Im_hi, &
                              Ip_src, Ips_lo, Ips_hi, &
                              Im_src, Ims_lo, Ims_hi, &
                              Ip_gc, Ipg_lo, Ipg_hi, &
                              Im_gc, Img_lo, Img_hi, &
                              qm, qm_lo, qm_hi, &
                              qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                              dloga, dloga_lo, dloga_hi, &
#endif
                              vlo, vhi, domlo, domhi, &
                              dx, dt)
       else
          call trace_ppm_gammae(lo, hi, &
                                idir, q, qd_lo, qd_hi, &
                                qaux, qa_lo, qa_hi, &
                                Ip, Ip_lo, Ip_hi, &
                                Im, Im_lo, Im_hi, &
                                Ip_src, Ips_lo, Ips_hi, &
                                Im_src, Ims_lo, Ims_hi, &
                                Ip_gc, Ipg_lo, Ipg_hi, &
                                Im_gc, Img_lo, Img_hi, &
                                qm, qm_lo, qm_hi, &
                                qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                                dloga, dloga_lo, dloga_hi, &
#endif
                                vlo, vhi, domlo, domhi, &
                                dx, dt)
       end if
    else
       call trace_ppm_temp(lo, hi, &
                           idir, q, qd_lo, qd_hi, &
                           qaux, qa_lo, qa_hi, &
                           Ip, Ip_lo, Ip_hi, &
                           Im, Im_lo, Im_hi, &
                           Ip_src, Ips_lo, Ips_hi, &
                           Im_src, Ims_lo, Ims_hi, &
                           Ip_gc, Ipg_lo, Ipg_hi, &
                           Im_gc, Img_lo, Img_hi, &
                           qm, qm_lo, qm_hi, &
                           qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vlo, vhi, domlo, domhi, &
                           dx, dt)
    end if

  end subroutine trace_ppm



  subroutine trace_ppm_rhoe(lo, hi, &
                            idir, q, qd_lo, qd_hi, &
                            qaux, qa_lo, qa_hi, &
                            flatn, f_lo, f_hi, &
                            Ip, Ip_lo, Ip_hi, &
                            Im, Im_lo, Im_hi, &
                            Ip_src, Ips_lo, Ips_hi, &
                            Im_src, Ims_lo, Ims_hi, &
                            Ip_gc, Ipg_lo, Ipg_hi, &
                            Im_gc, Img_lo, Img_hi, &
                            qm, qm_lo, qm_hi, &
                            qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                            dloga, dloga_lo, dloga_hi, &
#endif
                            vlo, vhi, domlo, domhi, &
                            dx, dt)

    use network, only : nspec, naux

    use meth_params_module, only : QVAR, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QGAME, QC, QGAMC, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, &
                                   fix_mass_flux, ppm_flatten_before_integrals
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    
    real(rt), intent(in) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:3,QVAR)
    real(rt), intent(in) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:3,QVAR)

    real(rt), intent(in) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:3,NQSRC)
    real(rt), intent(in) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:3,NQSRC)

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g
    real(rt) :: xi, xi1

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref, h_g_ref

    real(rt) :: cc_ref, csq_ref, gam_g_ref
    real(rt) :: cc_ev, csq_ev, rho_ev, h_g_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt

#if AMREX_SPACEDIM == 1
    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (vlo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (vhi(1) == domhi(1))
#endif

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2

             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p + rhoe_g)/rho)/csq

             gam_g = qaux(i,j,k,QGAMC)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,1,QRHO)
                un_ref    = Im(i,j,k,1,QUN)

                p_ref    = Im(i,j,k,1,QPRES)
                rhoe_g_ref = Im(i,j,k,1,QREINT)

                gam_g_ref  = Im_gc(i,j,k,1,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by un-c
                ! *p are the jumps carried by un+c

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(i,j,k,1,QUN) - hdt*Im_src(i,j,k,1,QUN)
                dptotm = p_ref - Im(i,j,k,1,QPRES) - hdt*Im_src(i,j,k,1,QPRES)

                drho = rho_ref - Im(i,j,k,2,QRHO) - hdt*Im_src(i,j,k,2,QRHO)
                dptot = p_ref - Im(i,j,k,2,QPRES) - hdt*Im_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,2,QREINT) - hdt*Im_src(i,j,k,2,QREINT)

                dup = un_ref - Im(i,j,k,3,QUN) - hdt*Im_src(i,j,k,3,QUN)
                dptotp = p_ref - Im(i,j,k,3,QPRES) - hdt*Im_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   h_g_ev = h_g
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   h_g_ev = h_g_ref
                endif

                ! (rho, u, p, (rho e) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e_g = drhoe_g - dptot*h_g_ev  ! note h_g has a 1/c**2 in it

                alpham = merge(ZERO, -alpham, un-cc > ZERO)
                alphap = merge(ZERO, -alphap, un+cc > ZERO)
                alpha0r = merge(ZERO, -alpha0r, un > ZERO)
                alpha0e_g = merge(ZERO, -alpha0e_g, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                qp(i,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                qp(i,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(i,j,k,2,QUT) + hdt*Im_src(i,j,k,2,QUT)
                qp(i,j,k,QUTT) = Im(i,j,k,2,QUTT) + hdt*Im_src(i,j,k,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,3,QRHO)
                un_ref    = Ip(i,j,k,3,QUN)

                p_ref    = Ip(i,j,k,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,3,QREINT)

                gam_g_ref  = Ip_gc(i,j,k,3,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(i,j,k,1,QUN) - hdt*Ip_src(i,j,k,1,QUN)
                dptotm  = p_ref - Ip(i,j,k,1,QPRES) - hdt*Ip_src(i,j,k,1,QPRES)

                drho = rho_ref - Ip(i,j,k,2,QRHO) - hdt*Ip_src(i,j,k,2,QRHO)
                dptot = p_ref - Ip(i,j,k,2,QPRES) - hdt*Ip_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,2,QREINT) - hdt*Ip_src(i,j,k,2,QREINT)

                dup = un_ref - Ip(i,j,k,3,QUN) - hdt*Ip_src(i,j,k,3,QUN)
                dptotp = p_ref - Ip(i,j,k,3,QPRES) - hdt*Ip_src(i,j,k,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   h_g_ev = h_g
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   h_g_ev = h_g_ref
                endif

                ! (rho, u, p, (rho e)) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e_g = drhoe_g - dptot*h_g_ev  ! h_g has a 1/c**2 in it

                alpham = merge(-alpham, ZERO, un-cc > ZERO)
                alphap = merge(-alphap, ZERO, un+cc > ZERO)
                alpha0r = merge(-alpha0r, ZERO, un > ZERO)
                alpha0e_g = merge(-alpha0e_g, ZERO, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   qm(i+1,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   if (ppm_flatten_before_integrals == 0) then
                      xi  = flatn(i,j,k)
                      xi1 = ONE - flatn(i,j,k)
 
                      qm(i+1,j,k,QRHO  ) = xi1*rho  + xi*qm(i+1,j,k,QRHO  )
                      qm(i+1,j,k,QUN   ) = xi1*un   + xi*qm(i+1,j,k,QUN   )
                      qm(i+1,j,k,QUT   ) = xi1*ut   + xi*qm(i+1,j,k,QUT   )
                      qm(i+1,j,k,QUTT  ) = xi1*utt  + xi*qm(i+1,j,k,QUTT  )
                      qm(i+1,j,k,QREINT) = xi1*rhoe_g + xi*qm(i+1,j,k,QREINT)
                      qm(i+1,j,k,QPRES ) = xi1*p    + xi*qm(i+1,j,k,QPRES )
                   endif
                   
                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   qm(i,j+1,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   if (ppm_flatten_before_integrals == 0) then
                      xi  = flatn(i,j,k)
                      xi1 = ONE - flatn(i,j,k)
 
                      qm(i,j+1,k,QRHO  ) = xi1*rho  + xi*qm(i,j+1,k,QRHO  )
                      qm(i,j+1,k,QUN   ) = xi1*un   + xi*qm(i,j+1,k,QUN   )
                      qm(i,j+1,k,QUT   ) = xi1*ut   + xi*qm(i,j+1,k,QUT   )
                      qm(i,j+1,k,QUTT  ) = xi1*utt  + xi*qm(i,j+1,k,QUTT  )
                      qm(i,j+1,k,QREINT) = xi1*rhoe_g + xi*qm(i,j+1,k,QREINT)
                      qm(i,j+1,k,QPRES ) = xi1*p    + xi*qm(i,j+1,k,QPRES )
                   endif
                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   qm(i,j,k+1,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   if (ppm_flatten_before_integrals == 0) then
                      xi  = flatn(i,j,k)
                      xi1 = ONE - flatn(i,j,k)
 
                      qm(i,j,k+1,QRHO  ) = xi1*rho  + xi*qm(i,j,k+1,QRHO  )
                      qm(i,j,k+1,QUN   ) = xi1*un   + xi*qm(i,j,k+1,QUN   )
                      qm(i,j,k+1,QUT   ) = xi1*ut   + xi*qm(i,j,k+1,QUT   )
                      qm(i,j,k+1,QUTT  ) = xi1*utt  + xi*qm(i,j,k+1,QUTT  )
                      qm(i,j,k+1,QREINT) = xi1*rhoe_g + xi*qm(i,j,k+1,QREINT)
                      qm(i,j,k+1,QPRES ) = xi1*p    + xi*qm(i,j,k+1,QPRES )
                   endif
                   
                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)
                endif

             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             ! these only apply for x states (idir = 1)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= vhi(1)) then

                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then

                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             end if
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(vlo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(vlo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(vlo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(vlo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(vhi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(vhi(1)+1,j,k,QUN   ) = q(domhi(1)+1,j,k,QUN  )
                qp(vhi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(vhi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
             end if
#endif
          end do
       end do
    end do


  end subroutine trace_ppm_rhoe

  subroutine trace_ppm_gammae(lo, hi, &
                              idir, q, qd_lo, qd_hi, &
                              qaux, qa_lo, qa_hi, &
                              Ip, Ip_lo, Ip_hi, &
                              Im, Im_lo, Im_hi, &
                              Ip_src, Ips_lo, Ips_hi, &
                              Im_src, Ims_lo, Ims_hi, &
                              Ip_gc, Ipg_lo, Ipg_hi, &
                              Im_gc, Img_lo, Img_hi, &
                              qm, qm_lo, qm_hi, &
                              qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                              dloga, dloga_lo, dloga_hi, &
#endif
                              vlo, vhi, domlo, domhi, &
                              dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QC, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, &
                                   fix_mass_flux, &
                                   gamma_minus_1
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:3,QVAR)
    real(rt), intent(in) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:3,QVAR)

    real(rt), intent(in) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:3,NQSRC)
    real(rt), intent(in) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:3,NQSRC)

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, Clag
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g, game

    real(rt) :: dptot
    real(rt) :: dge, dtau, dtaum, dtaup
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref
    real(rt) :: tau_ref

    real(rt) :: Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: Clag_ev, tau_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt

#if AMREX_SPACEDIM == 1
    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (vlo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (vhi(1) == domhi(1))
#endif

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)

             gam_g = (gamma_minus_1+ONE)
             game = (gamma_minus_1+ONE)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,1,QRHO)
                un_ref    = Im(i,j,k,1,QUN)

                p_ref    = Im(i,j,k,1,QPRES)
                rhoe_g_ref = Im(i,j,k,1,QREINT)

                tau_ref  = ONE/Im(i,j,k,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,1,1)
                game_ref = (gamma_minus_1+ONE)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(i,j,k,1,QUN) - hdt*Im_src(i,j,k,1,QUN)
                dptotm = p_ref - Im(i,j,k,1,QPRES) - hdt*Im_src(i,j,k,1,QPRES)

                dptot = p_ref - Im(i,j,k,2,QPRES) - hdt*Im_src(i,j,k,2,QPRES)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,1,QRHO) + hdt*Im_src(i,j,k,1,QRHO)/Im(i,j,k,1,QRHO)**2
                dtau  = tau_ref - ONE/Im(i,j,k,2,QRHO) + hdt*Im_src(i,j,k,2,QRHO)/Im(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,3,QRHO) + hdt*Im_src(i,j,k,3,QRHO)/Im(i,j,k,3,QRHO)**2

                dup = un_ref - Im(i,j,k,3,QUN) - hdt*Im_src(i,j,k,3,QUN)
                dptotp = p_ref - Im(i,j,k,3,QPRES) - hdt*Im_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                endif

                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system

                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                dge   = game_ref - (gamma_minus_1+ONE)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                alpham = merge(ZERO, -alpham, un-cc > ZERO)
                alphap = merge(ZERO, -alphap, un+cc > ZERO)
                alpha0r = merge(ZERO, -alpha0r, un > ZERO)
                alpha0e_g = merge(ZERO, -alpha0e_g, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qp(i,j,k,QRHO) = max(small_dens, ONE/tau_s)

                qp(i,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qp(i,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

!                qp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                qp(i,j,k,QREINT) = qp(i,j,k,QPRES )/gamma_minus_1


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(i,j,k,2,QUT) + hdt*Im_src(i,j,k,2,QUT)
                qp(i,j,k,QUTT) = Im(i,j,k,2,QUTT) + hdt*Im_src(i,j,k,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,3,QRHO)
                un_ref    = Ip(i,j,k,3,QUN)

                p_ref    = Ip(i,j,k,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,3,QREINT)

                tau_ref  = ONE/Ip(i,j,k,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,3,1)
                game_ref = (gamma_minus_1+ONE)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(i,j,k,1,QUN) - hdt*Ip_src(i,j,k,1,QUN)
                dptotm  = p_ref - Ip(i,j,k,1,QPRES) - hdt*Ip_src(i,j,k,1,QPRES)

                dptot = p_ref - Ip(i,j,k,2,QPRES) - hdt*Ip_src(i,j,k,2,QPRES)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,1,QRHO) + hdt*Ip_src(i,j,k,1,QRHO)/Ip(i,j,k,1,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,2,QRHO) + hdt*Ip_src(i,j,k,2,QRHO)/Ip(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,3,QRHO) + hdt*Ip_src(i,j,k,3,QRHO)/Ip(i,j,k,3,QRHO)**2

                dup = un_ref - Ip(i,j,k,3,QUN) - hdt*Ip_src(i,j,k,3,QUN)
                dptotp = p_ref - Ip(i,j,k,3,QPRES) - hdt*Ip_src(i,j,k,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                endif

                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                dge = game_ref - (gamma_minus_1+ONE)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                alpham = merge(-alpham, ZERO, un-cc > ZERO)
                alphap = merge(-alphap, ZERO, un+cc > ZERO)
                alpha0r = merge(-alpha0r, ZERO, un > ZERO)
                alpha0e_g = merge(-alpha0e_g, ZERO, un > ZERO)


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i+1,j,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i+1,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

!                   qm(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QPRES )/gamma_minus_1

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j+1,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j+1,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

!                   qm(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i,j+1,k,QREINT) = qm(i,j+1,k,QPRES )/gamma_minus_1

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j,k+1,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j,k+1,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

!                   qm(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i,j,k+1,QREINT) = qm(i,j,k+1,QPRES )/gamma_minus_1

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(i,j,k,2,QUT) + hdt*Ip_src(i,j,k,2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(i,j,k,2,QUTT) + hdt*Ip_src(i,j,k,2,QUTT)

                end if
             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             ! these only apply for x states (dim = 1)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                h_g = ( (p + rhoe_g)/rho)/csq
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= vhi(1)) then
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             endif
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(vlo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(vlo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(vlo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(vlo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(vhi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(vhi(1)+1,j,k,QUN   ) = q(domhi(1)+1,j,k,QUN  )
                qp(vhi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(vhi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
             end if
#endif
          end do
       end do
    end do

  end subroutine trace_ppm_gammae


  subroutine trace_ppm_temp(lo, hi, &
                            idir, q, qd_lo, qd_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Ip_lo, Ip_hi, &
                            Im, Im_lo, Im_hi, &
                            Ip_src, Ips_lo, Ips_hi, &
                            Im_src, Ims_lo, Ims_hi, &
                            Ip_gc, Ipg_lo, Ipg_hi, &
                            Im_gc, Img_lo, Img_hi, &
                            qm, qm_lo, qm_hi, &
                            qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                            dloga, dloga_lo, dloga_hi, &
#endif
                            vlo, vhi, domlo, domhi, &
                            dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QC, QFS, QFX, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, &
                                   fix_mass_flux
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)

#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:3,QVAR)
    real(rt), intent(in) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:3,QVAR)

    real(rt), intent(in) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:3,NQSRC)
    real(rt), intent(in) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:3,NQSRC)

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, Clag
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g, temp
    real(rt) :: gam_g, game

    real(rt) :: drho, dptot
    real(rt) :: dtau, dtaum, dtaup
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm
    real(rt) :: dT0, dTp, dTm
    real(rt) :: p_r, p_T

    real(rt) :: rho_ref, un_ref, p_ref, temp_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, tau_ev, temp_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt

#if AMREX_SPACEDIM == 1
    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (vlo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (vhi(1) == domhi(1))
#endif


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    print*, "eos fail"
    stop

  end subroutine trace_ppm_temp

end module ca_trace_ppm_module
