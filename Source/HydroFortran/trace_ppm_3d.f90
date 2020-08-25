module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

    subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_g,Im_g, &
                           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QREINT, QPRES, &
                                   npassive, qpass_map, ppm_type, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) qxm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qxp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qym (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qyp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

    real(rt) dt, a_old

    ! Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq, rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref

    real(rt) drho, du, dv, dw, dp, drhoe
    real(rt) dup, dvp, dpp
    real(rt) dum, dvm, dpm

    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) halfdt
    real(rt) csqref,cref,enthref

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = HALF * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).

    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    ! Version 2 includes the source terms in the jump in 
    ! the velocities that goes through the characteristic projection.

    ! *********************************************************************************************
    ! x-direction
    ! *********************************************************************************************

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho  = q(i,j,k3d,QRHO)
          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)
          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
       !  csq = cc**2
       !  enth = ( (rhoe+p)/rho )/csq

          ! ******************************************************************************

          if (i .ge. ilo1) then

             ! Plus state on face i

             ! Set the reference state
                 ! This will be the fastest moving state to the left
                 rho_ref = Im(i,j,kc,1,1,QRHO)
                   u_ref = Im(i,j,kc,1,1,QU)
                   v_ref = Im(i,j,kc,1,1,QV)
                   w_ref = Im(i,j,kc,1,1,QW)
                   p_ref = Im(i,j,kc,1,1,QPRES)
                rhoe_ref = Im(i,j,kc,1,1,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)

             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum   =   u_ref - Im(i,j,kc,1,1,QU)    - halfdt*Im_g(i,j,kc,1,1,QU)/a_old
             dpm   =   p_ref - Im(i,j,kc,1,1,QPRES) - halfdt*Im_g(i,j,kc,1,1,QPRES)/a_old
   
             drho  =  rho_ref - Im(i,j,kc,1,2,QRHO)   - halfdt*Im_g(i,j,kc,1,2,QRHO)/a_old
             dp    =    p_ref - Im(i,j,kc,1,2,QPRES)  - halfdt*Im_g(i,j,kc,1,2,QPRES)/a_old
             drhoe = rhoe_ref - Im(i,j,kc,1,2,QREINT) - halfdt*Im_g(i,j,kc,1,2,QREINT)/a_old
   
             dup   =    u_ref - Im(i,j,kc,1,3,QU)     - halfdt*Im_g(i,j,kc,1,3,QU)/a_old
             dpp   =    p_ref - Im(i,j,kc,1,3,QPRES)  - halfdt*Im_g(i,j,kc,1,3,QPRES)/a_old
   

            ! These are analogous to the beta's from the original PPM
            ! paper (except we work with rho instead of tau).  This is
            ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dpm/(rho_ref*cref) - dum)*rho_ref/cref
             alphap = HALF*(dpp/(rho_ref*cref) + dup)*rho_ref/cref
             alpha0r = drho - dp/csqref
             alpha0e = drhoe - dp*enthref

             if (u-cc .gt. ZERO) then
                amright = ZERO
             else
                amright = -alpham
             endif

             if (u+cc .gt. ZERO) then
                apright = ZERO
             else
                apright = -alphap
             endif

             if (u .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else
                azrright = -alpha0r
                azeright = -alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             qxp(i,j,kc,QRHO  ) = max(small_dens,rho_ref +  apright + amright + azrright)
             qxp(i,j,kc,QU    ) =    u_ref + (apright - amright)*cref/rho_ref
             qxp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enthref*csqref + azeright
             qxp(i,j,kc,QPRES ) = max(small_pres, p_ref + (apright + amright)*csqref)

             ! Transverse velocities -- there's no projection here, so we don't
             ! need a reference state.  We only care about the state traced under
             ! the middle wave
             dv = Im(i,j,kc,1,2,QV) + halfdt*Im_g(i,j,kc,1,2,QV)/a_old
             dw = Im(i,j,kc,1,2,QW) + halfdt*Im_g(i,j,kc,1,2,QW)/a_old

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface

             qxp(i,j,kc,QV    ) = dv
             qxp(i,j,kc,QW    ) = dw

             qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES) / gamma_minus_1

          end if

          ! ******************************************************************************

          if (i .le. ihi1) then

             ! Minus state on face i+1

             ! Set the reference state
                 ! This will be the fastest moving state to the right
                 rho_ref = Ip(i,j,kc,1,3,QRHO)
                   u_ref = Ip(i,j,kc,1,3,QU)
                   v_ref = Ip(i,j,kc,1,3,QV)
                   w_ref = Ip(i,j,kc,1,3,QW)
                   p_ref = Ip(i,j,kc,1,3,QPRES)
                rhoe_ref = Ip(i,j,kc,1,3,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)

             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum   =    u_ref - Ip(i,j,kc,1,1,QU   ) - halfdt*Ip_g(i,j,kc,1,1,QU   )/a_old
             dpm   =    p_ref - Ip(i,j,kc,1,1,QPRES) - halfdt*Ip_g(i,j,kc,1,1,QPRES)/a_old
   
             drho  =  rho_ref - Ip(i,j,kc,1,2,QRHO  ) - halfdt*Ip_g(i,j,kc,1,2,QRHO)/a_old
             dp    =    p_ref - Ip(i,j,kc,1,2,QPRES ) - halfdt*Ip_g(i,j,kc,1,2,QPRES)/a_old
             drhoe = rhoe_ref - Ip(i,j,kc,1,2,QREINT) - halfdt*Ip_g(i,j,kc,1,2,QREINT)/a_old

             dup   =    u_ref - Ip(i,j,kc,1,3,QU   ) - halfdt*Ip_g(i,j,kc,1,3,QU)/a_old
             dpp   =    p_ref - Ip(i,j,kc,1,3,QPRES) - halfdt*Ip_g(i,j,kc,1,3,QPRES)/a_old
   
             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This is
             ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dpm/(rho_ref*cref) - dum)*rho_ref/cref
             alphap = HALF*(dpp/(rho_ref*cref) + dup)*rho_ref/cref
             alpha0r = drho - dp/csqref
             alpha0e = drhoe - dp*enthref

             if (u-cc .gt. ZERO) then
                amleft = -alpham
             else
                amleft = ZERO
             endif

             if (u+cc .gt. ZERO) then
                apleft = -alphap
             else
                apleft = ZERO 
             endif
   
             if (u .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else
                azrleft = ZERO
                azeleft = ZERO
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qxm(i+1,j,kc,QRHO  ) = max(small_dens, rho_ref + apleft + amleft + azrleft)
             qxm(i+1,j,kc,QU    ) =    u_ref + (apleft - amleft)*cref/rho_ref
             qxm(i+1,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enthref*csqref + azeleft
             qxm(i+1,j,kc,QPRES ) = max(small_pres, p_ref + (apleft + amleft)*csqref)

             ! Transverse velocities
             dv    = Ip(i,j,kc,1,2,QV) + halfdt*Ip_g(i,j,kc,1,2,QV)/a_old
             dw    = Ip(i,j,kc,1,2,QW) + halfdt*Ip_g(i,j,kc,1,2,QW)/a_old

             qxm(i+1,j,kc,QV    ) = dv
             qxm(i+1,j,kc,QW    ) = dw

             qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES) / gamma_minus_1
          end if

       end do
    end do

    ! ******************************************************************************
    ! Passively advected quantities 
    ! ******************************************************************************

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Plus state on face i
          do i = ilo1, ihi1+1
 
                qxp(i,j,kc,n) = Im(i,j,kc,1,2,n) + halfdt*Im_g(i,j,kc,1,2,n)/a_old
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
                qxm(i+1,j,kc,n) = Ip(i,j,kc,1,2,n) + halfdt*Ip_g(i,j,kc,1,2,n)/a_old 
          enddo

       enddo
    enddo

    ! *********************************************************************************************
    ! y-direction
    ! *********************************************************************************************

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          if (j .ge. ilo2) then

             ! Plus state on face j

             ! Set the reference state
                 ! This will be the fastest moving state to the left
                 rho_ref = Im(i,j,kc,2,1,QRHO)
                   u_ref = Im(i,j,kc,2,1,QU)
                   v_ref = Im(i,j,kc,2,1,QV)
                   w_ref = Im(i,j,kc,2,1,QW)
                   p_ref = Im(i,j,kc,2,1,QPRES)
                rhoe_ref = Im(i,j,kc,2,1,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)

             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - Im(i,j,kc,2,1,QV) - halfdt*Im_g(i,j,kc,2,1,QV)/a_old
             dpm   =    p_ref - Im(i,j,kc,2,1,QPRES) - halfdt*Im_g(i,j,kc,2,1,QPRES)/a_old
   
             drho  =  rho_ref - Im(i,j,kc,2,2,QRHO) - halfdt*Im_g(i,j,kc,2,2,QRHO)/a_old
             dp    =    p_ref - Im(i,j,kc,2,2,QPRES) - halfdt*Im_g(i,j,kc,2,2,QPRES)/a_old
             drhoe = rhoe_ref - Im(i,j,kc,2,2,QREINT) - halfdt*Im_g(i,j,kc,2,2,QREINT)/a_old

             dvp   =    v_ref - Im(i,j,kc,2,3,QV) - halfdt*Im_g(i,j,kc,2,3,QV)/a_old
             dpp   =    p_ref - Im(i,j,kc,2,3,QPRES) - halfdt*Im_g(i,j,kc,2,3,QPRES)/a_old
   

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This
             ! is simply (l . dq), where dq = qref - I(q)
 
             alpham = HALF*(dpm/(rho_ref*cref) - dvm)*rho_ref/cref
             alphap = HALF*(dpp/(rho_ref*cref) + dvp)*rho_ref/cref
             alpha0r = drho - dp/csqref
             alpha0e = drhoe - dp*enthref
 
             if (v-cc .gt. ZERO) then
                amright = ZERO
             else
                amright = -alpham
             endif
 
             if (v+cc .gt. ZERO) then
                apright = ZERO
             else
                apright = -alphap
             endif
 
             if (v .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else
                azrright = -alpha0r
                azeright = -alpha0e
             endif
 
             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qyp(i,j,kc,QRHO  ) = max(small_dens, rho_ref +  apright + amright + azrright)
             qyp(i,j,kc,QV    ) =    v_ref + (apright - amright)*cref/rho_ref
             qyp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enthref*csqref + azeright
             qyp(i,j,kc,QPRES ) = max(small_pres, p_ref + (apright + amright)*csqref)

             ! Transverse velocities
             du    = Im(i,j,kc,2,2,QU)
             dw    = Im(i,j,kc,2,2,QW)
   
             du  = du  + halfdt*Im_g(i,j,kc,2,2,QU)/a_old
             dw  = dw  + halfdt*Im_g(i,j,kc,2,2,QW)/a_old
 
             qyp(i,j,kc,QU    ) = du
             qyp(i,j,kc,QW    ) = dw

             qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES) / gamma_minus_1
          end if

          ! ******************************************************************************

          if (j .le. ihi2) then

             ! Minus state on face j+1

             ! Set the reference state
                 ! This will be the fastest moving state to the right
                 rho_ref = Ip(i,j,kc,2,3,QRHO)
                   u_ref = Ip(i,j,kc,2,3,QU)
                   v_ref = Ip(i,j,kc,2,3,QV)
                   w_ref = Ip(i,j,kc,2,3,QW)
                   p_ref = Ip(i,j,kc,2,3,QPRES)
                rhoe_ref = Ip(i,j,kc,2,3,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)

             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm   =    v_ref - Ip(i,j,kc,2,1,QV) - halfdt*Ip_g(i,j,kc,2,1,QV)/a_old
             dpm   =    p_ref - Ip(i,j,kc,2,1,QPRES) - halfdt*Ip_g(i,j,kc,2,1,QPRES)/a_old
   
             drho  =  rho_ref - Ip(i,j,kc,2,2,QRHO) - halfdt*Ip_g(i,j,kc,2,2,QRHO)/a_old
             dp    =    p_ref - Ip(i,j,kc,2,2,QPRES) - halfdt*Ip_g(i,j,kc,2,2,QPRES)/a_old
             drhoe = rhoe_ref - Ip(i,j,kc,2,2,QREINT) - halfdt*Ip_g(i,j,kc,2,2,QREINT)/a_old

             dvp   =    v_ref - Ip(i,j,kc,2,3,QV) - halfdt*Ip_g(i,j,kc,2,3,QV)/a_old
             dpp   =    p_ref - Ip(i,j,kc,2,3,QPRES) - halfdt*Ip_g(i,j,kc,2,3,QPRES)/a_old


             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
 
             alpham = HALF*(dpm/(rho_ref*cref) - dvm)*rho_ref/cref
             alphap = HALF*(dpp/(rho_ref*cref) + dvp)*rho_ref/cref
             alpha0r = drho - dp/csqref
             alpha0e = drhoe - dp*enthref
 
             if (v-cc .gt. ZERO) then
                amleft = -alpham
             else
                amleft = ZERO
             endif

             if (v+cc .gt. ZERO) then
                apleft = -alphap
             else
                apleft = ZERO
             endif

             if (v .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else
                azrleft = ZERO
                azeleft = ZERO
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             qym(i,j+1,kc,QRHO  ) = max(small_dens, rho_ref +  apleft + amleft + azrleft)
             qym(i,j+1,kc,QV    ) =    v_ref + (apleft - amleft)*cref/rho_ref
             qym(i,j+1,kc,QREINT) = rhoe_ref + (apleft + amleft)*enthref*csqref + azeleft
             qym(i,j+1,kc,QPRES ) = max(small_pres, p_ref + (apleft + amleft)*csqref)

             ! Transverse velocities
             du    = Ip(i,j,kc,2,2,QU)
             dw    = Ip(i,j,kc,2,2,QW)

             du  = du  + halfdt*Ip_g(i,j,kc,2,2,QU)/a_old
             dw  = dw  + halfdt*Ip_g(i,j,kc,2,2,QW)/a_old
 
             qym(i,j+1,kc,QU    ) = du
             qym(i,j+1,kc,QW    ) = dw

             qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES) / gamma_minus_1
          end if

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1

          ! Plus state on face j
          do j = ilo2, ihi2+1
                qyp(i,j,kc,n) = Im(i,j,kc,2,2,n) + halfdt*Im_g(i,j,kc,2,2,n)/a_old
          enddo
          
          ! Minus state on face j+1
          do j = ilo2-1, ihi2
                qym(i,j+1,kc,n) = Ip(i,j,kc,2,2,n) + halfdt*Ip_g(i,j,kc,2,2,n)/a_old
          enddo
          
       enddo
    enddo

    end subroutine tracexy_ppm

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

    subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          Ip,Im,Ip_g,Im_g, &
                          qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                          ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QREINT, QPRES, &
                                   npassive, qpass_map, ppm_type, &
                                   npassive, qpass_map, ppm_type, &
                                   small_dens, small_pres, gamma_minus_1
    use amrex_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt)   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt)   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    real(rt) Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    real(rt) qzm (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) qzp (qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    real(rt) dt, a_old

    !     Local variables
    integer i, j
    integer n, ipassive

    real(rt) cc, csq
    real(rt) rho, u, v, w, p, rhoe
    real(rt) rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    real(rt) dwp, dpp
    real(rt) dwm, dpm

    real(rt) drho, du, dv, dp, drhoe
    real(rt) enth, alpham, alphap, alpha0r, alpha0e
    real(rt) apright, amright, azrright, azeright
    real(rt) apleft, amleft, azrleft, azeleft
    real(rt) halfdt
    real(rt) cref,csqref,enthref

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    halfdt = HALF * dt

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM

    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! **************************************************************************
          ! This is all for qzp
          ! **************************************************************************

          rho  = q(i,j,k3d,QRHO)
          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)
          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc   = c(i,j,k3d)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Plus state on face kc

          ! Set the reference state
                 ! This will be the fastest moving state to the left
              rho_ref = Im(i,j,kc,3,1,QRHO)
                u_ref = Im(i,j,kc,3,1,QU)
                v_ref = Im(i,j,kc,3,1,QV)
                w_ref = Im(i,j,kc,3,1,QW)
                p_ref = Im(i,j,kc,3,1,QPRES)
             rhoe_ref = Im(i,j,kc,3,1,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)

             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   =    w_ref - Im(i,j,kc,3,1,QW) - halfdt*Im_g(i,j,kc,3,1,QW)/a_old
          dpm   =    p_ref - Im(i,j,kc,3,1,QPRES) - halfdt*Im_g(i,j,kc,3,1,QPRES)/a_old

          drho  =  rho_ref - Im(i,j,kc,3,2,QRHO) - halfdt*Im_g(i,j,kc,3,2,QRHO)/a_old
          dp    =    p_ref - Im(i,j,kc,3,2,QPRES) - halfdt*Im_g(i,j,kc,3,2,QPRES)/a_old
          drhoe = rhoe_ref - Im(i,j,kc,3,2,QREINT) - halfdt*Im_g(i,j,kc,3,2,QREINT)/a_old

          dwp   =    w_ref - Im(i,j,kc,3,3,QW) - halfdt*Im_g(i,j,kc,3,3,QW)/a_old
          dpp   =    p_ref - Im(i,j,kc,3,3,QPRES) - halfdt*Im_g(i,j,kc,3,3,QPRES)/a_old

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dpm/(rho_ref*cref) - dwm)*rho_ref/cref
          alphap = HALF*(dpp/(rho_ref*cref) + dwp)*rho_ref/cref
          alpha0r = drho - dp/csqref
          alpha0e = drhoe - dp*enthref

          if (w-cc .gt. ZERO) then
             amright = ZERO
          else
             amright = -alpham
          endif
          if (w+cc .gt. ZERO) then
             apright = ZERO
          else
             apright = -alphap
          endif
          if (w .gt. ZERO) then
             azrright = ZERO
             azeright = ZERO
          else
             azrright = -alpha0r
             azeright = -alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzp(i,j,kc,QRHO  ) = max(small_dens, rho_ref +  apright + amright + azrright)
          qzp(i,j,kc,QW    ) =    w_ref + (apright - amright)*cref/rho_ref
          qzp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enthref*csqref + azeright
          qzp(i,j,kc,QPRES ) = max(small_pres, p_ref + (apright + amright)*csqref)

          ! Transverse velocities
          du    = Im(i,j,kc,3,2,QU)
          dv    = Im(i,j,kc,3,2,QV)

          du  = du  + halfdt*Im_g(i,j,kc,3,2,QU)/a_old
          dv  = dv  + halfdt*Im_g(i,j,kc,3,2,QV)/a_old

          qzp(i,j,kc,QU    ) = du
          qzp(i,j,kc,QV    ) = dv

          qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES) / gamma_minus_1

          ! **************************************************************************
          ! This is all for qzm
          ! **************************************************************************

          ! Minus state on face kc

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)
          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)
          p    = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Set the reference state
              ! This will be the fastest moving state to the right
              rho_ref = Ip(i,j,km,3,3,QRHO)
                u_ref = Ip(i,j,km,3,3,QU)
                v_ref = Ip(i,j,km,3,3,QV)
                w_ref = Ip(i,j,km,3,3,QW)
                p_ref = Ip(i,j,km,3,3,QPRES)
             rhoe_ref = Ip(i,j,km,3,3,QREINT)

             rho_ref = max(rho_ref,small_dens)
               p_ref = max(  p_ref,small_pres)
 
             csqref = (1.d0+gamma_minus_1)*p_ref/rho_ref
             cref = sqrt(csqref)
             enthref = (rhoe_ref+p_ref)/(rho_ref*csqref)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = (   w_ref - Ip(i,j,km,3,1,QW)) - halfdt*Ip_g(i,j,km,3,1,QW)/a_old
          dpm   = (   p_ref - Ip(i,j,km,3,1,QPRES)) - halfdt*Ip_g(i,j,km,3,1,QPRES)/a_old

          drho  = ( rho_ref - Ip(i,j,km,3,2,QRHO)) - halfdt*Ip_g(i,j,km,3,2,QRHO)/a_old
          dp    = (   p_ref - Ip(i,j,km,3,2,QPRES)) - halfdt*Ip_g(i,j,km,3,2,QPRES)/a_old
          drhoe = (rhoe_ref - Ip(i,j,km,3,2,QREINT)) - halfdt*Ip_g(i,j,km,3,2,QREINT)/a_old

          dwp   = (   w_ref - Ip(i,j,km,3,3,QW)) - halfdt*Ip_g(i,j,km,3,3,QW)/a_old
          dpp   = (   p_ref - Ip(i,j,km,3,3,QPRES)) - halfdt*Ip_g(i,j,km,3,3,QPRES)/a_old


          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)

          alpham = HALF*(dpm/(rho_ref*cref) - dwm)*rho_ref/cref
          alphap = HALF*(dpp/(rho_ref*cref) + dwp)*rho_ref/cref
          alpha0r = drho - dp/csqref
          alpha0e = drhoe - dp*enthref
             
          if (w-cc .gt. ZERO) then
             amleft = -alpham
          else
             amleft = ZERO
          endif
          if (w+cc .gt. ZERO) then
             apleft = -alphap
          else
             apleft = ZERO
          endif
          if (w .gt. ZERO) then
             azrleft = -alpha0r
             azeleft = -alpha0e
          else
             azrleft = ZERO
             azeleft = ZERO
          endif
          
          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          qzm(i,j,kc,QRHO  ) = max(small_dens, rho_ref +  apleft + amleft + azrleft)
          qzm(i,j,kc,QW    ) =    w_ref + (apleft - amleft)*cref/rho_ref
          qzm(i,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enthref*csqref + azeleft
          qzm(i,j,kc,QPRES ) = max(small_pres, p_ref + (apleft + amleft)*csqref)

          ! Transverse velocity
          du = Ip(i,j,km,3,2,QU)
          dv = Ip(i,j,km,3,2,QV)

          du  = du  + halfdt*Ip_g(i,j,km,3,2,QU)/a_old
          dv  = dv  + halfdt*Ip_g(i,j,km,3,2,QV)/a_old
 
          qzm(i,j,kc,QU    ) = du
          qzm(i,j,kc,QV    ) = dv

          qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES) / gamma_minus_1

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Plus state on face kc
                     qzp(i,j,kc,n) = Im(i,j,kc,3,2,n) + halfdt*Im_g(i,j,kc,3,2,n)/a_old 
                     qzm(i,j,kc,n) = Ip(i,j,km,3,2,n) + halfdt*Ip_g(i,j,km,3,2,n)/a_old 

               enddo
         enddo
    enddo

    end subroutine tracez_ppm

end module trace_ppm_module
