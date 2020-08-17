#ifndef _RIEMANN_H_
#define _RIEMANN_H_
#include "IndexDefines.H"

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
  const amrex::Real rr,
  const amrex::Real ur,
  const amrex::Real vr,
  const amrex::Real v2r,
  const amrex::Real pr,
  const amrex::Real rer,
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
  const amrex::Real& small_in,
  const amrex::Real& csmall,
  const amrex::Real& gamma,
  const int NumSpec)
{
  const amrex::Real wsmall = small_dens * csmall;

  amrex::Real gdnv_state_rho = rl;
  amrex::Real gdnv_state_p = pl;

  amrex::Real gdnv_state_cs;

  gdnv_state_rho = rr;
  gdnv_state_p = pr;

  amrex::Real pstar;

#if 0
  const amrex::Real wl =
    amrex::max(wsmall, std::sqrt(amrex::Math::abs(gamma * pl * rl)));
  const amrex::Real wr =
    amrex::max(wsmall, std::sqrt(amrex::Math::abs(gamma * pr * rr)));

  const amrex::Real csl = wl/rl;
  const amrex::Real csr = wr/rr;

  const amrex::Real pstar = amrex::max(
    small_pres, ((wr * pl + wl * pr) + wl * wr * (ul - ur)) / (wl + wr));
  ustar = ((wl * ul + wr * ur) + (pl - pr)) / (wl + wr);
#else
{
      const amrex::Real weakwv  = 1.e-3;
      const amrex::Real  small  = 1.e-6;

      amrex::Real wl = amrex::Math::sqrt(amrex::Math::abs(gamma*pl*amrex::max(rl,small_dens)));
      amrex::Real wr = amrex::Math::sqrt(amrex::Math::abs(gamma*pr*amrex::max(rr,small_dens)));

      amrex::Real cleft  = wl/amrex::max(rl,small_dens);
      amrex::Real cright = wr/amrex::max(rr,small_dens);

      pstar = amrex::max( (wl*pr+wr*pl-wr*wl*(ur-ul))/(wl+wr), small_pres );

      amrex::Real pstnm1 = pstar;

      amrex::Real wlsq = (.5*(gamma-1.0)*(pstar+pl)+pstar) * amrex::max(rl,small_dens);
      amrex::Real wrsq = (.5*(gamma-1.0)*(pstar+pr)+pstar) * amrex::max(rr,small_dens);

      wl = amrex::Math::sqrt(amrex::Math::abs(wlsq));
      wr = amrex::Math::sqrt(amrex::Math::abs(wrsq));

      amrex::Real ustarp = ul - (pstar-pl)/wl;
      amrex::Real ustarm = ur + (pstar-pr)/wr;

      pstar = amrex::max( (wl*pr+wr*pl-wr*wl*(ur-ul))/(wl+wr), small_pres );

      for (int iter = 0; iter < 3; iter++)
      {
        wlsq = (.5*(gamma-1.0)*(pstar+pl)+pstar) * amrex::max(rl,small_dens);
        wrsq = (.5*(gamma-1.0)*(pstar+pr)+pstar) * amrex::max(rr,small_dens);

        wl = 1.0/sqrt(amrex::Math::abs(wlsq));
        wr = 1.0/sqrt(amrex::Math::abs(wrsq));

        amrex::Real ustnm1 = ustarm;
        amrex::Real ustnp1 = ustarp;

        ustarm = ur - (pr-pstar)*wr;
        ustarp = ul + (pl-pstar)*wl;

        amrex::Real dpditer = amrex::Math::abs(pstnm1-pstar);
        amrex::Real zp      = amrex::Math::abs(ustarp-ustnp1);

        if (zp-weakwv*cleft < 0.)
            zp = dpditer*wl;

        amrex::Real zm = amrex::Math::abs(ustarm-ustnm1);

        if (zm-weakwv*cright < 0.) 
           zm = dpditer*wr;

        amrex::Real denom  = dpditer / amrex::max(zp+zm,small*(cleft+cright));
        pstnm1 = pstar;

        pstar = pstar - denom*(ustarm-ustarp);

        pstar = amrex::max(pstar,small_pres);

        ustar = 0.5*(ustarm+ustarp);
     }
}
#endif

  bool mask = ustar > 0.0;
  amrex::Real ro = mask ? amrex::max(rl,small_dens) : amrex::max(rr,small_dens);
  amrex::Real uo = mask ? ul : ur;
  amrex::Real po = mask ? pl : pr;
  const amrex::Real  gamcr = gamma;
  const amrex::Real  gamcl = gamma;
  // Make data types/vars match:
  //             rel = ql(i,j,k,QREINT,comp)
  //             rer = qr(i,j,k,QREINT,comp)
  amrex::Real reo = mask ? rel : rer;
  amrex::Real gamco = mask ? gamcl : gamcr;
  //0 case and also small ustar case
  //Figure out how to set small_vel to be consistent with:
  //          ro = max(small_dens, ro)
  mask = ( ( amrex::Math::abs(ustar) < small_vel * 0.5 * 
            (amrex::Math::abs(ul) + amrex::Math::abs(ur)) ) || ustar == 0.0 );
  ustar = mask ? 0.0 : ustar;
  ro = mask ? 0.5 * (amrex::max(rl,small_dens) + amrex::max(rr,small_dens)) : ro;
  uo = mask ? 0.5 * (ul + ur) : uo;
  po = mask ? 0.5 * (pl + pr) : po;
  reo = mask ? 0.5 * (rel + rer) : reo;
  gamco = mask ? 0.5 * (gamcl + gamcr) : gamco;

  const amrex::Real co = amrex::max(csmall,std::sqrt(std::abs(gamco * po / ro)));

  //Is a max or limiter similar to the mask above needed for co?
  const amrex::Real drho = (pstar - po) / (co * co);
  const amrex::Real rstar = amrex::max(small_dens, ro + drho);

  //Fix data type:
  const amrex::Real entho = (reo + po)/(ro*co*co);
  const amrex::Real estar = reo + (pstar - po)*entho;

  amrex::Real cstar = std::sqrt(amrex::Math::abs(gamco*pstar/rstar));
  cstar = amrex::max(cstar, csmall);
  const amrex::Real sgnm = amrex::Math::copysign(1.0, ustar);

  amrex::Real spout = co - sgnm * uo;
  amrex::Real spin = cstar - sgnm * ustar;
  const amrex::Real ushock = 0.5 * (spin + spout);

  mask = pstar < po;
  spout = mask ? spout : ushock;
  spin = mask ? spin : ushock;

  /////////////////// Continue here
  const amrex::Real scr = (spout == spin) ? small_in * cav : spout - spin;
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

  // Might need rho factor
  amrex::Real regd = gdnv_state_p / (gamma - 1.0);

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

  // Might need rho factor
  regd = gdnv_state_p / (gamma - 1.0);

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
}

#endif
