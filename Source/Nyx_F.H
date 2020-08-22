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

  void fort_set_method_params
    (const int& dm, const int& NumAdv, const int& Ndiag, const int& do_hydro,
     const int& ppm_type, 
     const int& use_flattening,
     const int& use_analriem,
     const int& use_const_species,
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

  void get_grav_source_type(int* grav_source_type);

  void fort_get_num_spec(int* nspec);

  void fort_get_num_aux(int* naux);

  void fort_get_spec_names(int* spec_names, int* ispec, int* len);

  void fort_get_aux_names(int* aux_names, int* iaux, int* len);

  void fort_get_UFS(int* ufs);

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

  void fort_enforce_nonnegative_species
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

  AMREX_GPU_DEVICE void ca_enforce_nonnegative_species
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

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

  AMREX_GPU_DEVICE void apply_av
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(div),
     const BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(flux), 
     const int idir, const amrex::Real* dx, 
     const amrex::Real dt);

  AMREX_GPU_DEVICE void normalize_species_fluxes
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux));

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

  void adjust_heat_cool
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt);

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
