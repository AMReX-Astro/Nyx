#include <Godunov.H>
#include <PPM.H>

using namespace amrex;

void
trace_ppm(const Box& bx,
          const int idir,
          Array4<Real const> const& q_arr,
          Array4<Real const> const& srcQ,
          Array4<Real> const& qm,
          Array4<Real> const& qp,
          const Box& vbx,
          const Real dt, 
          const Real* dx,
          const Real gamma,
          const Real small_dens, const Real small_pres,
#ifndef CONST_SPECIES
          const int NumSpec,
#endif
          const Real a_old)
{

  // here, lo and hi are the range we loop over -- this can include ghost cells
  // vlo and vhi are the bounds of the valid box (no ghost cells)

  //
  // rho : mass density
  // u, v, w : velocities
  // p : gas (hydro) pressure
  // ptot : total pressure (note for pure hydro, this is
  //        just the gas pressure)
  // cgas : sound speed for just the gas contribution
  // cc : total sound speed 
  //
  // for pure hydro, we will only consider:
  //    rho, u, v, w, ptot, cc

  Real hdtovera = 0.5_rt * dt / a_old;
  Real dtdx = dt / dx[idir];
  Real dtdxovera = dtdx / a_old;

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();

  // This does the characteristic tracing to build the interface
  // states using the normal predictor only (no transverse terms).
  //
  // For each zone, we construct Im and Ip arrays -- these are the averages
  // of the various primitive state variables under the parabolic
  // interpolant over the region swept out by one of the 3 different
  // characteristic waves.
  //
  // Im is integrating to the left interface of the current zone
  // (which will be used to build the right ("p") state at that interface)
  // and Ip is integrating to the right interface of the current zone
  // (which will be used to build the left ("m") state at that interface).
  //
  //
  // The choice of reference state is designed to minimize the
  // effects of the characteristic projection.  We subtract the I's
  // off of the reference state, project the quantity such that it is
  // in terms of the characteristic varaibles, and then add all the
  // jumps that are moving toward the interface to the reference
  // state to get the full state on that interface.

  int QUN, QUT, QUTT;

  if (idir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (idir == 1) {
    QUN = QV;
    QUT = QW;
    QUTT = QU;
  } else if (idir == 2) {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  Real lsmall_dens = small_dens;
  Real lsmall_pres = small_pres;

  // Trace to left and right edges using upwind PPM
  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real cc = std::sqrt(gamma * q_arr(i,j,k,QPRES)/q_arr(i,j,k,QRHO));

    Real un = q_arr(i,j,k,QUN);

    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[5];

    Real flat = 1.0;
    //Calculate flattening in-place
    for(int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++)
    {
        flat = amrex::min(flat,flatten(i, j, k, dir_flat, q_arr, small_pres));
    }

    Real sm;
    Real sp;

    Real Ip[QVAR][3];
    Real Im[QVAR][3];

    for (int n = 0; n < QVAR; n++) {

      if (idir == 0) {
        s[im2] = q_arr(i-2,j,k,n);
        s[im1] = q_arr(i-1,j,k,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i+1,j,k,n);
        s[ip2] = q_arr(i+2,j,k,n);

      } else if (idir == 1) {
        s[im2] = q_arr(i,j-2,k,n);
        s[im1] = q_arr(i,j-1,k,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i,j+1,k,n);
        s[ip2] = q_arr(i,j+2,k,n);

      } else {
        s[im2] = q_arr(i,j,k-2,n);
        s[im1] = q_arr(i,j,k-1,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i,j,k+1,n);
        s[ip2] = q_arr(i,j,k+2,n);

      }

      ppm_reconstruct(s, flat, sm, sp);
      ppm_int_profile(sm, sp, s[i0], un, cc, dtdxovera, Ip[n], Im[n]);

    }

    // source terms
    Real Ip_src[QVAR][3];
    Real Im_src[QVAR][3];

    for (int n = 0; n < QVAR; n++) 
    {
        if (idir == 0) {
          s[im2] = srcQ(i-2,j,k,n);
          s[im1] = srcQ(i-1,j,k,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i+1,j,k,n);
          s[ip2] = srcQ(i+2,j,k,n);

        } else if (idir == 1) {
          s[im2] = srcQ(i,j-2,k,n);
          s[im1] = srcQ(i,j-1,k,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i,j+1,k,n);
          s[ip2] = srcQ(i,j+2,k,n);

        } else {
          s[im2] = srcQ(i,j,k-2,n);
          s[im1] = srcQ(i,j,k-1,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i,j,k+1,n);
          s[ip2] = srcQ(i,j,k+2,n);
        }

        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdxovera, Ip_src[n], Im_src[n]);
    }

#ifndef CONST_SPECIES
    for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {

      // Plus state on face i
      if ((idir == 0 && i >= vlo[0]) ||
          (idir == 1 && j >= vlo[1]) ||
          (idir == 2 && k >= vlo[2])) {

        // We have
        //
        // q_l = q_ref - Proj{(q_ref - I)}
        //
        // and Proj{} represents the characteristic projection.
        // But for these, there is only 1-wave that matters, the u
        // wave, so no projection is needed.  Since we are not
        // projecting, the reference state doesn't matter

        qp(i,j,k,n) = Im[n][1];
      }

      // Minus state on face i+1
      if (idir == 0 && i <= vhi[0]) {
        qm(i+1,j,k,n) = Ip[n][1];

      } else if (idir == 1 && j <= vhi[1]) {
        qm(i,j+1,k,n) = Ip[n][1];

      } else if (idir == 2 && k <= vhi[2]) {
        qm(i,j,k+1,n) = Ip[n][1];
      }
    }
#endif

    // plus state on face i

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

      // Set the reference state
      // This will be the fastest moving state to the left --
      // this is the method that Miller & Colella and Colella &
      // Woodward use
      Real rho_ref = Im[QRHO][0];
      Real un_ref  = Im[QUN][0];

      Real p_ref      = Im[QPRES][0];

      Real gam_g_ref = gamma;

      rho_ref = amrex::max(rho_ref, lsmall_dens);
        p_ref = amrex::max(  p_ref, lsmall_pres);

      Real rho_ref_inv = 1.0_rt/rho_ref;

      // For tracing
      Real csq_ref = gam_g_ref*p_ref*rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;

      // *m are the jumps carried by un-c
      // *p are the jumps carried by un+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)

      // we add the sources here so they participate in the tracing
      Real dum    = un_ref - Im[QUN][0]   - hdtovera*Im_src[QUN][0];
      Real dptotm =  p_ref - Im[QPRES][0] - hdtovera*Im_src[QPRES][0];

      Real drho    = rho_ref    - Im[QRHO][1]   - hdtovera*Im_src[QRHO][1];
      Real dptot   = p_ref      - Im[QPRES][1]  - hdtovera*Im_src[QPRES][1];

      Real dup     = un_ref     - Im[QUN][2]    - hdtovera*Im_src[QUN][2];
      Real dptotp  = p_ref      - Im[QPRES][2]  - hdtovera*Im_src[QPRES][2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv;
      Real alphap = 0.5_rt*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv;
      Real alpha0r = drho - dptot/csq_ref;

      // These terms only have the velocity source terms
      Real dum_src = -hdtovera*Im_src[QUN][0];
      Real dup_src = -hdtovera*Im_src[QUN][2];

      // These terms only have the energy/pressure source terms
      Real dptotm_src = - hdtovera*Im_src[QPRES][0];
      Real dptot_src  = - hdtovera*Im_src[QPRES][1];
      Real dptotp_src = - hdtovera*Im_src[QPRES][2];

      // These are the definitions that will replace zero in the upwinding
      Real alpham_src  = 0.5_rt*(dptotm_src*rho_ref_inv*cc_ref_inv - dum_src)*rho_ref*cc_ref_inv;
      Real alphap_src  = 0.5_rt*(dptotp_src*rho_ref_inv*cc_ref_inv + dup_src)*rho_ref*cc_ref_inv;
      Real alpha0r_src =       - dptot_src/csq_ref;
          
      alpham    = un-cc > 0.0_rt ? -alpham_src  : -alpham;
      alphap    = un+cc > 0.0_rt ? -alphap_src  : -alphap;
      alpha0r   = un    > 0.0_rt ? -alpha0r_src : -alpha0r;
      
      // The final interface states are just
      // q_s = q_ref - sum(l . dq) r
      // note that the a{mpz}right as defined above have the minus already

      Real rho_pred = rho_ref +  alphap + alpham + alpha0r;

      qp(i,j,k,QRHO ) = amrex::max(lsmall_dens, rho_pred);
      qp(i,j,k,QUN  ) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
      qp(i,j,k,QPRES) = amrex::max(lsmall_pres, p_ref + (alphap + alpham)*csq_ref);

      // Transverse velocities -- there's no projection here, so we
      // don't need a reference state.  We only care about the state
      // traced under the middle wave

      // Recall that I already takes the limit of the parabola
      // in the event that the wave is not moving toward the
      // interface
      qp(i,j,k,QUT)  = Im[QUT][1]  + hdtovera*Im_src[QUT][1];
      qp(i,j,k,QUTT) = Im[QUTT][1] + hdtovera*Im_src[QUTT][1];

      // This allows the (rho e) to take advantage of (pressure > small_pres)
      qp(i,j,k,QREINT) = qp(i,j,k,QPRES) / (gamma - 1.0);
    }

    // minus state on face i + 1

    if ((idir == 0 && i <= vhi[0]) ||
        (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {

      // Set the reference state
      // This will be the fastest moving state to the right
      Real rho_ref = Ip[QRHO][2];
      Real un_ref = Ip[QUN][2];

      Real p_ref = Ip[QPRES][2];

      Real gam_g_ref = gamma;

      rho_ref = amrex::max(rho_ref, lsmall_dens);
        p_ref = amrex::max(  p_ref, lsmall_pres);

      Real rho_ref_inv = 1.0_rt/rho_ref;

      // For tracing
      Real csq_ref = gam_g_ref*p_ref*rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      Real dum     = un_ref - Ip[QUN][0]   - hdtovera*Ip_src[QUN][0];
      Real dptotm  =  p_ref - Ip[QPRES][0] - hdtovera*Ip_src[QPRES][0];

      Real drho    = rho_ref -    Ip[QRHO][1]   - hdtovera*Ip_src[QRHO][1];
      Real dptot   =   p_ref    - Ip[QPRES][1]  - hdtovera*Ip_src[QPRES][1];

      Real dup    = un_ref - Ip[QUN][2]   - hdtovera*Ip_src[QUN][2];
      Real dptotp =  p_ref - Ip[QPRES][2] - hdtovera*Ip_src[QPRES][2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv;
      Real alphap = 0.5_rt*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv;
      Real alpha0r = drho - dptot/csq_ref;

      // These terms only have the velocity source terms
      Real dum_src = -hdtovera*Ip_src[QUN][0];
      Real dup_src = -hdtovera*Ip_src[QUN][2];

      // These terms only have the energy/pressure source terms
      Real dptotm_src = - hdtovera*Ip_src[QPRES][0];
      Real dptot_src  = - hdtovera*Ip_src[QPRES][1];
      Real dptotp_src = - hdtovera*Ip_src[QPRES][2];

      // These are the definitions that will replace zero in the upwinding
      Real alpham_src  = 0.5_rt*(dptotm_src*rho_ref_inv*cc_ref_inv - dum_src)*rho_ref*cc_ref_inv;
      Real alphap_src  = 0.5_rt*(dptotp_src*rho_ref_inv*cc_ref_inv + dup_src)*rho_ref*cc_ref_inv;
      Real alpha0r_src =       - dptot_src/csq_ref;

      alpham  = un-cc > 0.0_rt ? -alpham  : -alpham_src;
      alphap  = un+cc > 0.0_rt ? -alphap  : -alphap_src;
      alpha0r = un    > 0.0_rt ? -alpha0r : -alpha0r_src;

      Real rho_pred = rho_ref +  alphap + alpham + alpha0r;

      // The final interface states are just
      // q_s = q_ref - sum (l . dq) r
      // note that the a{mpz}left as defined above have the minus already

      if (idir == 0) {

        qm(i+1,j,k,QRHO ) = amrex::max(lsmall_dens, rho_pred);
        qm(i+1,j,k,QUN  ) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i+1,j,k,QPRES) = amrex::max(lsmall_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i+1,j,k,QUT) = Ip[QUT][1]   + hdtovera*Ip_src[QUT][1];
        qm(i+1,j,k,QUTT) = Ip[QUTT][1] + hdtovera*Ip_src[QUTT][1];

        // This allows the (rho e) to take advantage of (pressure > small_pres)
        qm(i+1,j,k,QREINT) = qm(i+1,j,k,QPRES) / (gamma - 1.0);

      } else if (idir == 1) {

        qm(i,j+1,k,QRHO ) = amrex::max(lsmall_dens, rho_pred);
        qm(i,j+1,k,QUN  ) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j+1,k,QPRES) = amrex::max(lsmall_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j+1,k,QUT) = Ip[QUT][1]   + hdtovera*Ip_src[QUT][1];
        qm(i,j+1,k,QUTT) = Ip[QUTT][1] + hdtovera*Ip_src[QUTT][1];

        // This allows the (rho e) to take advantage of (pressure > small_pres)
        qm(i,j+1,k,QREINT) = qm(i,j+1,k,QPRES) / (gamma - 1.0);

      } else if (idir == 2) {

        qm(i,j,k+1,QRHO ) = amrex::max(lsmall_dens, rho_pred);
        qm(i,j,k+1,QUN  ) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j,k+1,QPRES) = amrex::max(lsmall_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j,k+1,QUT) = Ip[QUT][1]   + hdtovera*Ip_src[QUT][1];
        qm(i,j,k+1,QUTT) = Ip[QUTT][1] + hdtovera*Ip_src[QUTT][1];

        // This allows the (rho e) to take advantage of (pressure > small_pres)
        qm(i,j,k+1,QREINT) = qm(i,j,k+1,QPRES) / (gamma - 1.0);
      }
    }
  });
}



