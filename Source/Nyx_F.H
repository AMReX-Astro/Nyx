#ifndef _Nyx_F_H_
#define _Nyx_F_H_

#include <AMReX_BLFort.H>
#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void fort_alloc_cuda_managed();
  void fort_dealloc_cuda_managed();

  void fort_alloc_simd_vec();
  void fort_dealloc_simd_vec();

  void fort_integrate_time_given_a
    (amrex::Real* a0, amrex::Real* a1, amrex::Real* dt);

  void fort_integrate_comoving_a
    (amrex::Real* old_a, amrex::Real* new_a, amrex::Real* dt);

  void fort_integrate_comoving_a_to_z 
    (amrex::Real* old_a, amrex::Real* z_value, amrex::Real* dt);

  void fort_integrate_comoving_a_to_a 
    (amrex::Real* old_a, amrex::Real* a_value, amrex::Real* dt);

  void set_simd(const int *simd_width);

  //  void fort_get_omm    (amrex::Real* omm);
  //  void fort_get_omb    (amrex::Real* frac);
  void fort_get_hubble (amrex::Real* hubble);

  void fort_set_omm    (const amrex::Real& omm);
  void fort_set_omb    (const amrex::Real& frac);
  void fort_set_omr    (const amrex::Real& omr);
  void fort_set_hubble (const amrex::Real& hubble);

  AMREX_GPU_DEVICE void fort_estdt
  (const BL_FORT_FAB_ARG(state),
     const int* lo, const int* hi,
     const amrex::Real* dx, amrex::Real* dt, amrex::Real comoving_a);

  void fort_estdt_comoving_a
    (amrex::Real* old_a, amrex::Real* new_dummy_a, amrex::Real* dt, amrex::Real* change_allowed,
     amrex::Real* fixed_da_interval, amrex::Real* final_a, int* dt_modified);

  void fort_get_method_params(int* HYP_GROW);

  void fort_set_method_params
    (const int& dm, const int& NumAdv, const int& Ndiag, const int& do_hydro,
     const int& ppm_type, const int& ppm_ref,
     const int& ppm_flatten_before_integrals,
     const int& use_colglaz, const int& use_flattening,
     const int& use_analriem,
     const int& version_2, const int& use_const_species,
     const amrex::Real& gamma_in, const int& normalize_species,
     const int& heat_cool_type, const int& inhomo_reion);

  void fort_tabulate_rates();

  void filcc
    (const amrex::Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
     const int * domlo, const int * domhi,
     const amrex::Real * dx_crse, const amrex::Real * xlo,
     const int * bc);

  void fort_network_init();

  void fort_set_xhydrogen(amrex::Real& xhydrogen_in);

  void fort_get_num_spec(int* nspec);

  void fort_get_num_aux(int* naux);

  void fort_get_spec_names(int* spec_names, int* ispec, int* len);

  void fort_get_aux_names(int* aux_names, int* iaux, int* len);

  void fort_set_eos_params(const amrex::Real& h_species_in, const amrex::Real& he_species_in);

  void fort_set_small_values
    (const amrex::Real* average_dens, const amrex::Real* average_temp,
     const amrex::Real* comoving_a,
     const amrex::Real* small_dens,
     const amrex::Real* small_temp, const amrex::Real* small_pres);

  void fort_set_problem_params
    (const int& dm, const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, const int& Symmetry_value,
     const int& coord_type);

  void fort_initdata
    (const int& level, const amrex::Real& time, const int* lo, const int* hi,
     const int& num_state, BL_FORT_FAB_ARG(state),
     const int& num_diag,  BL_FORT_FAB_ARG(diag_eos), const amrex::Real dx[],
     const amrex::Real xlo[], const amrex::Real xhi[]);

  AMREX_GPU_DEVICE void ca_fort_initdata
    (const int& level, const amrex::Real& time, const int* lo, const int* hi,
     const int& num_state, BL_FORT_FAB_ARG(state),
     const int& num_diag,  BL_FORT_FAB_ARG(diag_eos), const amrex::Real dx[],
     const amrex::Real* z_in);

  void fort_init_zhi
    (const int* lo, const int* hi, 
     const int& num_diag, BL_FORT_FAB_ARG(diag_eos),
     const int& ratio, BL_FORT_FAB_ARG(zhi));

  void fort_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

  AMREX_GPU_DEVICE void fort_init_e_from_t
     (BL_FORT_FAB_ARG(state), const int* num_state,
      BL_FORT_FAB_ARG( diag), const int* num_diag,
      const int* lo, const int* hi, const amrex::Real* comoving_a);

  void fort_init_e_from_rhoe
     (BL_FORT_FAB_ARG(state), const int* num_state,
      const int* lo, const int* hi, amrex::Real* comoving_a);

  void fort_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

  void fort_enforce_nonnegative_species
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

  AMREX_GPU_DEVICE void ca_enforce_nonnegative_species
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

  AMREX_GPU_DEVICE void reset_internal_e
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(D_new),
     BL_FORT_FAB_ARG(reset_e_src),
     const int* print_fortran_warnings,
     const amrex::Real* comoving_a,
     const int* interp);

  AMREX_GPU_DEVICE void reset_internal_e_nostore
    (const int* lo, const int* hi,
     amrex::Real *S_new, const int, const int, const int, const int, const int, const int, 
     amrex::Real *D_new, const int, const int, const int, const int, const int, const int, 
     const int print_fortran_warnings,
     const amrex::Real comoving_a);

  void hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void generic_fill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void fort_make_hydro_sources
    (const amrex::Real* time,
     const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(ugdnvx),
     const BL_FORT_FAB_ARG(ugdnvy),
     const BL_FORT_FAB_ARG(ugdnvz),
     const BL_FORT_FAB_ARG(src),
           BL_FORT_FAB_ARG(hydro_src),
           BL_FORT_FAB_ARG(divu_cc),
     const BL_FORT_FAB_ARG(grav),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(const BL_FORT_FAB_ARG(xflux),
            const BL_FORT_FAB_ARG(yflux),
            const BL_FORT_FAB_ARG(zflux)),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const int* print_fortran_warnings);

  AMREX_GPU_DEVICE void ca_ctoprim(const int* lo, const int* hi,
	          const amrex::Real* u, const int* u_lo, const int* u_hi,
#ifdef RADIATION
		  const amrex::Real* Er, const int* Er_lo, const int* Er_hi,
		  const amrex::Real* lam, const int* lam_lo, const int* lam_hi,
#endif
		  amrex::Real* q, const int* q_lo, const int* q_hi,
		  amrex::Real* qaux, const int* qa_lo, const int* qa_hi );

  AMREX_GPU_DEVICE void ca_srctoprim(const int* lo,   const int* hi,
		    const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		    const amrex::Real* grav, const int* ga_lo, const int* ga_hi,
		    const amrex::Real* src, const int* src_lo, const int* src_hi,
		    amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi,
		    const amrex::Real a_old, const amrex::Real a_new, const amrex::Real dt);

  void ca_compute_cfl(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_ANYD(q),
                      const BL_FORT_FAB_ARG_ANYD(qaux),
                      const amrex::Real dt, const amrex::Real* dx,
                      amrex::Real* courno, const int verbose);

  AMREX_GPU_DEVICE void divu
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div));

  AMREX_GPU_DEVICE void ca_uflatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn), int pres_comp);

  void ca_rad_flatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(flatg));

  AMREX_GPU_DEVICE void ca_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(Sborder),
     BL_FORT_FAB_ARG(hydro_source),
      D_DECL(const BL_FORT_FAB_ARG(xflux),
             const BL_FORT_FAB_ARG(yflux),
             const BL_FORT_FAB_ARG(zflux)),
             const BL_FORT_FAB_ARG_3D(xqe),
             const BL_FORT_FAB_ARG_3D(yqe),
             const BL_FORT_FAB_ARG_3D(zqe),
     /*
     BL_FORT_FAB_ARG(flux0),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG(flux1),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG(flux2),
#endif
     */
     const BL_FORT_FAB_ARG_3D(divu_cc),
     const amrex::Real* dx, const amrex::Real dt,
     const amrex::Real a_old, const amrex::Real a_new);
  
#ifndef AMREX_USE_GPU
  void ctu_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(Sborder),
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(hydro_source),
     const BL_FORT_FAB_ARG_3D(flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(flux2),
#endif
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(Erborder),
     BL_FORT_FAB_ARG_3D(Er_out),
     const BL_FORT_FAB_ARG_3D(rad_flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(rad_flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(rad_flux2),
#endif
     int* nstep_fsp,
#endif
     const BL_FORT_FAB_ARG_3D(qe0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(qe1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(qe2),
#endif
     const BL_FORT_FAB_ARG_3D(area0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area2),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(pdivu),
     const amrex::Real* dx, const amrex::Real dt);
#endif

  AMREX_GPU_DEVICE void scale_flux
    (const int* lo, const int* hi,
#if AMREX_SPACEDIM == 1
     const BL_FORT_FAB_ARG_3D(qe0),
#endif
     BL_FORT_FAB_ARG_3D(flux),
     const amrex::Real area,
     const amrex::Real dt, const amrex::Real a_old, const amrex::Real a_new);

  void scale_rad_flux
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rad_flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);

  void ca_track_grid_losses
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(flux2),
#endif
     amrex::Real* mass_lost, amrex::Real* xmom_lost, amrex::Real* ymom_lost, amrex::Real* zmom_lost,
     amrex::Real* eden_lost, amrex::Real* xang_lost, amrex::Real* yang_lost, amrex::Real* zang_lost);

  AMREX_GPU_DEVICE void cmpflx_plus_godunov
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(ql),
     const BL_FORT_FAB_ARG_3D(qr), const int nc, const int comp,
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(q_int),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(radflux),
     BL_FORT_FAB_ARG_3D(lambda_int),
#endif
     BL_FORT_FAB_ARG_3D(qe),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(shk),
     const int idir, const int* domlo, const int* domhi);

#ifdef AMREX_USE_GPU
  void ca_ppm_reconstruct
    (const int* lo, const int* hi, int put_on_edges, int idir,
     BL_FORT_FAB_ARG_3D(q), int nc, int nstart, int nend,
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp), int nqc, int nqstart, int nqend);

  void ca_shock
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(shk),
     const amrex::Real* dx);

  AMREX_GPU_DEVICE void ca_store_godunov_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(qint),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(lambda_int),
#endif
     BL_FORT_FAB_ARG_3D(qgdnv));

  void ca_construct_flux_cuda
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real dt,
     const int idir,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     BL_FORT_FAB_ARG_3D(qe),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area));

  void ca_construct_hydro_update_cuda
    (const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(qx),
     BL_FORT_FAB_ARG_3D(qy),
     BL_FORT_FAB_ARG_3D(qz),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     BL_FORT_FAB_ARG_3D(xarea),
     BL_FORT_FAB_ARG_3D(yarea),
     BL_FORT_FAB_ARG_3D(zarea),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update));
#endif

  AMREX_GPU_DEVICE void apply_av
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(div),
     const BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(flux), 
     const int idir, const amrex::Real* dx, 
     const amrex::Real dt);

#ifdef RADIATION
  void apply_av_rad
    (const int* lo, const int* hi,
     int idir, const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(radflux));
#endif

  void limit_hydro_fluxes_on_small_dens
    (const int* lo, const int* hi,
     const int idir,
     BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(volume),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt, const amrex::Real* dx);

  AMREX_GPU_DEVICE void normalize_species_fluxes
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux));

  void ca_ctu_update
    (const int* lo, const int* hi,
     const int* is_finest_level,
     const amrex::Real time,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(Erout),
#endif
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(src_q),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(xflux),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(yflux),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(zflux),
#endif
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(rxflux),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(ryflux),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(rzflux),
#endif
#endif
     const BL_FORT_FAB_ARG_3D(xarea),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(yarea),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(zarea),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
     const BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
#ifdef RADIATION
     const int* priv_nstep_fsp,
#endif
     amrex::Real& mass_lost,
     amrex::Real& xmom_lost,
     amrex::Real& ymom_lost,
     amrex::Real& zmom_lost,
     amrex::Real& eden_lost,
     amrex::Real& xang_lost,
     amrex::Real& yang_lost,
     amrex::Real& zang_lost);

  void ctu_ppm_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(srcQ),
      const BL_FORT_FAB_ARG_3D(shk),
      BL_FORT_FAB_ARG_3D(Ip),
      BL_FORT_FAB_ARG_3D(Im),
      BL_FORT_FAB_ARG_3D(Ip_src),
      BL_FORT_FAB_ARG_3D(Im_src),
      BL_FORT_FAB_ARG_3D(Ip_gc),
      BL_FORT_FAB_ARG_3D(Im_gc),
      BL_FORT_FAB_ARG_3D(sm),
      BL_FORT_FAB_ARG_3D(sp),
      BL_FORT_FAB_ARG_3D(qxm),
      BL_FORT_FAB_ARG_3D(qxp),
#if AMREX_SPACEDIM >= 2
      BL_FORT_FAB_ARG_3D(qym),
      BL_FORT_FAB_ARG_3D(qyp),
#endif
#if AMREX_SPACEDIM == 3
      BL_FORT_FAB_ARG_3D(qzm),
      BL_FORT_FAB_ARG_3D(qzp),
#endif
      const amrex::Real* dx, const amrex::Real dt,
#if AMREX_SPACEDIM < 3
      const BL_FORT_FAB_ARG_3D(dloga),
#endif
      const int* domlo, const int* domhi);

  AMREX_GPU_DEVICE void ctu_plm_states
     (const int* lo, const int* hi,
      const int idir,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(srcQ),
      BL_FORT_FAB_ARG_3D(shk),
      BL_FORT_FAB_ARG_3D(dq),
      BL_FORT_FAB_ARG_3D(qm),
      BL_FORT_FAB_ARG_3D(qp),
      const amrex::Real* dx, const amrex::Real dt,
      const amrex::Real a_old, const amrex::Real a_new,
#if AMREX_SPACEDIM < 3
      const BL_FORT_FAB_ARG_3D(dloga),
#endif
      const int* domlo, const int* domhi);

  AMREX_GPU_DEVICE void transx_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
#if AMREX_SPACEDIM == 2
                         const BL_FORT_FAB_ARG_3D(area1),
                         const BL_FORT_FAB_ARG_3D(vol),
#endif
                         const amrex::Real hdt, const amrex::Real cdtdx);

  AMREX_GPU_DEVICE void transx_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
                         const amrex::Real hdt, const amrex::Real cdtdx);

  AMREX_GPU_DEVICE void transy_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  AMREX_GPU_DEVICE void transy_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  AMREX_GPU_DEVICE void transz_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  AMREX_GPU_DEVICE void transz_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  AMREX_GPU_DEVICE void transyz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fyz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyz),
#endif
               const BL_FORT_FAB_ARG_3D(fzy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzy),
#endif
               const BL_FORT_FAB_ARG_3D(qy),
               const BL_FORT_FAB_ARG_3D(qz),
	       const BL_FORT_FAB_ARG_3D(srcQ),
               const amrex::Real hdt, const amrex::Real cdtdy, const amrex::Real cdtdz,
	       const amrex::Real a_old, const amrex::Real a_new);

  AMREX_GPU_DEVICE void transxz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxz),
#endif
               const BL_FORT_FAB_ARG_3D(fzx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qz),
	       const BL_FORT_FAB_ARG_3D(srcQ),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdz,
	       const amrex::Real a_old, const amrex::Real a_new);

  AMREX_GPU_DEVICE void transxy(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxy),
#endif
               const BL_FORT_FAB_ARG_3D(fyx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qy),
	       const BL_FORT_FAB_ARG_3D(srcQ),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdy,
	       const amrex::Real a_old, const amrex::Real a_new);

  void ca_mol_single_stage
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const int* domlo, const int* domhi,
     const amrex::Real *stage_weight,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     BL_FORT_FAB_ARG_3D(update_flux),
     const amrex::Real* dx, const amrex::Real* dt,
     BL_FORT_FAB_ARG_3D(flux1),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(flux2),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(flux3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
#endif
     const BL_FORT_FAB_ARG_3D(area1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     const int& verbose);

  void ca_fourth_single_stage
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const int* domlo, const int* domhi,
     const amrex::Real *stage_weight,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(q_bar),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qaux_bar),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     BL_FORT_FAB_ARG_3D(update_flux),
     const amrex::Real* dx, const amrex::Real* dt,
     BL_FORT_FAB_ARG_3D(flux1),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(flux2),
#endif
#if AMREX_SPACEDIM >= 3
     BL_FORT_FAB_ARG_3D(flux3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
#endif
     const BL_FORT_FAB_ARG_3D(area1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     const int&  verbose);

  void ca_sdc_update_advection_o2
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     const int* m_start);

  void ca_sdc_update_advection_o4
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

#ifdef REACTIONS
  void ca_sdc_update_o2
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(R_0_old),
     BL_FORT_FAB_ARG_3D(R_1_old),
     const int* sdc_iteration,
     const int* m_start);

  void ca_sdc_update_centers_o4
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_center),
     BL_FORT_FAB_ARG_3D(U_new_center),
     const BL_FORT_FAB_ARG_3D(C_center),
     const int* sdc_iteration);

  void ca_sdc_conservative_update
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_old),
     BL_FORT_FAB_ARG_3D(U_new),
     const BL_FORT_FAB_ARG_3D(C),
     const BL_FORT_FAB_ARG_3D(R_new));

  void ca_sdc_compute_C4(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(A_m),
                         const BL_FORT_FAB_ARG_3D(A_0),
                         const BL_FORT_FAB_ARG_3D(A_1),
                         const BL_FORT_FAB_ARG_3D(A_2),
                         const BL_FORT_FAB_ARG_3D(R_0),
                         const BL_FORT_FAB_ARG_3D(R_1),
                         const BL_FORT_FAB_ARG_3D(R_2),
                         BL_FORT_FAB_ARG_3D(C),
                         const int* m_start);
#endif

  AMREX_GPU_DEVICE void ca_fort_update_state
    ( const int* lo, const int* hi,
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(src),
      const BL_FORT_FAB_ARG(hydro_src),
      const BL_FORT_FAB_ARG(divu_cc),
            BL_FORT_FAB_ARG_3D(sum_state),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new,
      const int* print_fortran_warnings);

  void fort_clean_state
    ( const int* lo, const int* hi,
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(src),
      const BL_FORT_FAB_ARG(hydro_src),
      const BL_FORT_FAB_ARG(divu_cc),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new,
      const int* print_fortran_warnings);

  AMREX_GPU_DEVICE void ca_fort_add_grav_source
    ( const int* lo, const int* hi,
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(grav),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new);

  void fort_update_state
    ( const int lo[], const int hi[],
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(src),
      const BL_FORT_FAB_ARG(hydro_src),
      const BL_FORT_FAB_ARG(divu_cc),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new,
      const int* print_fortran_warnings);

  void fort_add_grav_source
    ( const int lo[], const int hi[],
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(grav),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new);

  void time_center_sources
    (const int lo[], const int hi[], BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt, const int* print_fortran_warnings);

  void adjust_heat_cool
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt);

  AMREX_GPU_DEVICE void fort_correct_gsrc
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(grav_src_old),
     const BL_FORT_FAB_ARG(grad_phi_cc),
     const BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     const amrex::Real a_old, const amrex::Real a_new,
     const amrex::Real dt);

  AMREX_GPU_DEVICE void fort_syncgsrc
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(gPhi),
     const BL_FORT_FAB_ARG(gdPhi),
     const BL_FORT_FAB_ARG(S),
     const BL_FORT_FAB_ARG(dS),
     BL_FORT_FAB_ARG(src),
     const amrex::Real a_new,
     const amrex::Real dt);

  void sum_over_level
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const amrex::Real dx[], amrex::Real* sum);

  void sum_product
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2),
     const int lo[], const int hi[], const amrex::Real dx[], amrex::Real* sum);

  void sum_prod_prod
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2),
     BL_FORT_FAB_ARG(fab3),
     const int lo[], const int hi[], const amrex::Real dx[], amrex::Real* sum);

  void fort_avgdown
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);

#ifdef FORCING
  void fort_alloc_spect
    (const int* length);

  void fort_set_wavevector
    (const int k[], const int* indx);

  void fort_set_modes
    (const amrex::Real even[], const amrex::Real odd[],
     const int* length, const int* comp);
#endif

#ifdef AGN
  void fort_ext_src
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* particle_locs_and_mass,
     const amrex::Real* particle_data,
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* z,
     const amrex::Real* dt);
#else
  void fort_ext_src
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* z,
     const amrex::Real* dt);
#endif

  void integrate_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  void integrate_state_with_source
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     BL_FORT_FAB_ARG(diag_eos),
     BL_FORT_FAB_ARG(hydro_src),
     BL_FORT_FAB_ARG(reset_src),
     BL_FORT_FAB_ARG(IR),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  void integrate_state_fcvode_with_source
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     BL_FORT_FAB_ARG(diag_eos),
     BL_FORT_FAB_ARG(hydro_src),
     BL_FORT_FAB_ARG(reset_src),
     BL_FORT_FAB_ARG(IR),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  AMREX_GPU_DEVICE int RhsFnReal(amrex::Real t, amrex::Real* u, amrex::Real* udot, amrex::Real* rpar, int neq);
  AMREX_GPU_DEVICE void fort_ode_eos_finalize(double* e_out, double* rpar, int neq);
  void fort_update_eos(double dt, double* u, double* uout, double* rpar);
  void fort_ode_eos_setup(const amrex::Real& a,const amrex::Real& half_dt);

  //#ifdef HEATCOOL
  //  AMREX_GPU_DEVICE void fort_nyx_eos_T_given_Re_device(int JH, int JHe, amrex::Real* T, amrex::Real* Ne, amrex::Real R,amrex::Real e,amrex::Real comoving_a,amrex::Real*);
  //#else
  AMREX_GPU_DEVICE void fort_nyx_eos_T_given_Re_device(int JH, int JHe, amrex::Real* T, amrex::Real* Ne, amrex::Real R,amrex::Real e,amrex::Real comoving_a);
  //#endif
  AMREX_GPU_DEVICE void fort_nyx_eos_given_RT(amrex::Real* e, amrex::Real* P, amrex::Real R, amrex::Real T, amrex::Real Ne,amrex::Real comoving_a);
  
  AMREX_GPU_DEVICE void fort_compute_temp
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real comoving_a,
     const int print_fortran_warnings);

  void fort_compute_temp_host
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real* comoving_a,
     const int* print_fortran_warnings);

  void fort_compute_temp_vec
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* comoving_a,
     const int* print_fortran_warnings);

  void fort_interp_to_this_z
    (const amrex::Real* z);

  void fort_setup_eos_params
    (amrex::Real* eos_nr_eps,
     amrex::Real* vode_rtol,
     amrex::Real* vode_atol_scaled);

  void fort_compute_max_temp_loc
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real* max_temp, const amrex::Real* den_maxt,
     const int* imax, const int* jmax, const int* kmax);

  void fort_compute_rho_temp
    (const int lo[], const int hi[], const amrex::Real dx[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* rho_ave, amrex::Real* rho_T_sum, amrex::Real* T_sum,
     amrex::Real* Tinv_sum, amrex::Real* T_meanrho_sum, amrex::Real* rho_sum,
     amrex::Real* vol_sum, amrex::Real* vol_mn_sum);

  void fort_compute_gas_frac
    (const int lo[], const int hi[], const amrex::Real dx[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* rho_ave, amrex::Real* T_cut, amrex::Real* rho_cut,
     amrex::Real* whim_mass, amrex::Real* whim_vol,
     amrex::Real* hh_mass, amrex::Real* hh_vol,
     amrex::Real* igm_mass, amrex::Real* igm_vol,
     amrex::Real* mass_sum, amrex::Real* vol_sum);

#ifdef AUX_UPDATE
  void auxupdate
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const amrex::Real * dt);
#endif

  void get_rhoe
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));
#ifdef __cplusplus
}
#endif

#endif
