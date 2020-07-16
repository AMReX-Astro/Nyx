#ifndef _RIEMANN_H_
#define _RIEMANN_H_
#include "Nyx.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
riemann(
  const amrex::Real rl,
  const amrex::Real ul,
  const amrex::Real vl,
  const amrex::Real v2l,
  const amrex::Real pl,
  const amrex::Real rel,
  //  const amrex::Real spl[NumSpec],
  const amrex::Real rr,
  const amrex::Real ur,
  const amrex::Real vr,
  const amrex::Real v2r,
  const amrex::Real pr,
  const amrex::Real rer,
  //  const amrex::Real spr[NumSpec],
  const int bc_test_val,
  const amrex::Real c,
  const amrex::Real cav,
  amrex::Real& ustar,
  amrex::Real& uflx_rho,
  amrex::Real& uflx_u,
  amrex::Real& uflx_v,
  amrex::Real& uflx_w,
  amrex::Real& uflx_eden,
  amrex::Real& uflx_eint,
  amrex::Real& qint_iu,
  amrex::Real& qint_iv1,
  amrex::Real& qint_iv2,
  amrex::Real& qint_gdpres,
  amrex::Real& qint_gdgame,
  const amrex::Real& small_dens,
  const amrex::Real& small_pres,
  const amrex::Real& small_vel,
  const amrex::Real& small,
  const amrex::Real& gamma,
  const int NumSpec)
{
  const amrex::Real csmall = small * c;
  const amrex::Real wsmall = small_dens * csmall;

  amrex::Real gdnv_state_rho = rl;
  amrex::Real gdnv_state_p = pl;
  //  amrex::Real gdnv_state_massfrac[NumSpec];
  //  for (int n = 0; n < NumSpec; n++) {
  //     gdnv_state_massfrac[n] = spl[n];
  //  }

  amrex::Real gdnv_state_cs;
  // Put soundspeed here:
  //  EOS::RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, gdnv_state_cs);
  //  const amrex::Real csl = gdnv_state_cs;

  gdnv_state_rho = rr;
  gdnv_state_p = pr;
  //  for (int n = 0; n < NumSpec; n++) {
  //    gdnv_state_massfrac[n] = spr[n];
  //  }
  // Put soundspeed here:
  //  EOS::RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, gdnv_state_cs);
  //  const amrex::Real csr = gdnv_state_cs;

  amrex::Print() << "WSMALL " << wsmall << std::endl; 
  amrex::Print() << "PL PR " << pl  << " " << pr << std::endl; 
  amrex::Print() << "GAMMA " << gamma  << std::endl; 

  amrex::Print() << "NUMSPEC " << NumSpec  << std::endl; 

  const amrex::Real wl =
    amrex::max(wsmall, std::sqrt(amrex::Math::abs(gamma * pl * rl)));
  const amrex::Real wr =
    amrex::max(wsmall, std::sqrt(amrex::Math::abs(gamma * pr * rr)));

  const amrex::Real csl = wl/rl;
  const amrex::Real csr = wr/rr;

  const amrex::Real pstar = amrex::max(
    small_pres, ((wr * pl + wl * pr) + wl * wr * (ul - ur)) / (wl + wr));
  ustar = ((wl * ul + wr * ur) + (pl - pr)) / (wl + wr);

  bool mask = ustar > 0.0;
  amrex::Real ro = mask ? rl : rr;
  amrex::Real uo = mask ? ul : ur;
  amrex::Real po = mask ? pl : pr;
  const amrex::Real  gamcr = gamma;
  const amrex::Real  gamcl = gamma;
  // Make data types/vars match:
  //             rel = ql(i,j,k,QREINT,comp)
  //             rer = qr(i,j,k,QREINT,comp)
  //reo = mask ? rel : rer
  //gamco = mask ? gamcl : gamcr
  amrex::Real sp[NumSpec];
  for (int n = 0; n < NumSpec; n++) {
  //   sp[n] = mask ? spl[n] : spr[n];
  }
  //0 case and also small ustar case
  //Figure out how to set small_vel to be consistent with:
  //          ro = max(small_dens, ro)
  mask = ( ( amrex::Math::abs(ustar) < small_vel * 0.5 * 
            (amrex::Math::abs(ul) + amrex::Math::abs(ur)) ) || ustar == 0.0 );
  ustar = mask ? 0.0 : ustar;
  ro = mask ? 0.5 * (rl + rr) : ro;
  uo = mask ? 0.5 * (ul + ur) : uo;
  po = mask ? 0.5 * (pl + pr) : po;
  //reo = mask ? 0.5 * (rel + rer) : reo
  //gamco = mask ? 0.5 * (gamcl + gamcr) : gamco
  for (int n = 0; n < NumSpec; n++) {
//     sp[n] = mask ? 0.5 * (spl[n] + spr[n]) : sp[n];
  }
  /*
  exit(0); 

  gdnv_state_rho = ro;
  gdnv_state_p = po;
  for (int n = 0; n < NumSpec; n++) {
    gdnv_state_massfrac[n] = sp[n];
  }
  
  amrex::Real gdnv_state_e;
  //Put nyx_eos_given_rp ???
  //  EOS::RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  const amrex::Real reo = gdnv_state_rho * gdnv_state_e;
  // Put soundspeed here:
  //  EOS::RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, gdnv_state_cs);
  const amrex::Real co = gdnv_state_cs;
  */

  //From riemannus, looks equivalent:
  /*
             ro = max(small_dens, ro)

             roinv = ONE/ro

             co = sqrt(abs(gamco*po*roinv))
             co = max(csmall,co)
  */
  const amrex::Real co = std::sqrt(std::abs(gamco * po / ro));

  //Is a max or limiter similar to the mask above needed for co?

  const amrex::Real drho = (pstar - po) / (co * co);
  const amrex::Real rstar = amrex::max(small_dens, ro + drho);
  /*
  gdnv_state_rho = rstar;
  gdnv_state_p = pstar;
  for (int n = 0; n < NumSpec; n++) {
    gdnv_state_massfrac[n] = sp[n];
  }
  //Put nyx_eos_given_rp ???
  //  EOS::RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  const amrex::Real estar = gdnv_state_rho * gdnv_state_e;
  // Put soundspeed here:
  //  EOS::RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, gdnv_state_cs);
  const amrex::Real cstar = gdnv_state_cs;
*/
  //Fix data type:
             entho = (reo + po)/(ro*co*co)
             estar = reo + (pstar - po)*entho

             cstar = sqrt(abs(gamco*pstar/rstar))
             cstar = max(cstar, csmall)
  const amrex::Real sgnm = amrex::Math::copysign(1.0, ustar);

  amrex::Real spout = co - sgnm * uo;
  amrex::Real spin = cstar - sgnm * ustar;
  const amrex::Real ushock = 0.5 * (spin + spout);

  mask = pstar < po;
  spout = mask ? spout : ushock;
  spin = mask ? spin : ushock;

  /////////////////// Continue here
  const amrex::Real scr = (spout == spin) ? small * cav : spout - spin;
  const amrex::Real frac =
    amrex::max(0.0, amrex::min(1.0, (1.0 + (spout + spin) / scr) * 0.5));

  mask = ustar > 0.0;
  qint_iv1 = mask ? vl : vr;
  qint_iv2 = mask ? v2l : v2r;

  mask = (ustar == 0.0);
  qint_iv1 = mask ? 0.5 * (vl + vr) : qint_iv1;
  qint_iv2 = mask ? 0.5 * (v2l + v2r) : qint_iv2;
  amrex::Real rgd = frac * rstar + (1.0 - frac) * ro;
  qint_iu = frac * ustar + (1.0 - frac) * uo;
  qint_gdpres = frac * pstar + (1.0 - frac) * po;

  gdnv_state_rho = rgd;
  gdnv_state_p = qint_gdpres;
  for (int n = 0; n < NumSpec; n++) {
    gdnv_state_massfrac[n] = sp[n];
  }
  //Put nyx_eos_given_rp ???
  //  EOS::RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  amrex::Real regd = gdnv_state_rho * gdnv_state_e;

  mask = (spout < 0.0);
  rgd = mask ? ro : rgd;
  qint_iu = mask ? uo : qint_iu;
  qint_gdpres = mask ? po : qint_gdpres;
  regd = mask ? reo : regd;

  mask = (spin >= 0.0);
  rgd = mask ? rstar : rgd;
  qint_iu = mask ? ustar : qint_iu;
  qint_gdpres = mask ? pstar : qint_gdpres;
  regd = mask ? estar : regd;

  gdnv_state_rho = rgd;
  gdnv_state_p = qint_gdpres;
  for (int n = 0; n < NumSpec; n++) {
    gdnv_state_massfrac[n] = sp[n];
  }
  //Put nyx_eos_given_rp ???
  //  EOS::RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  regd = gdnv_state_rho * gdnv_state_e;

  qint_gdgame = qint_gdpres / regd + 1.0;
  qint_iu = bc_test_val * qint_iu;
  uflx_rho = rgd * qint_iu;
  uflx_u = uflx_rho * qint_iu + qint_gdpres;
  uflx_v = uflx_rho * qint_iv1;
  uflx_w = uflx_rho * qint_iv2;
  const amrex::Real rhoetot =
    regd +
    0.5 * rgd * (qint_iu * qint_iu + qint_iv1 * qint_iv1 + qint_iv2 * qint_iv2);
  uflx_eden = qint_iu * (rhoetot + qint_gdpres);
  uflx_eint = qint_iu * regd;
  amrex::Print() << "LEAVING ROUTINE " << std::endl;
}

#endif
