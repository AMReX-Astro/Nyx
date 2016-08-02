#ifndef _Nyx_F_H_
#define _Nyx_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FILCC, filcc)
    (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
     const int * domlo, const int * domhi,
     const Real * dx_crse, const Real * xlo,
     const int * bc);

BL_FORT_PROC_DECL(NYX_NETWORK_INIT, nyx_network_init)();

BL_FORT_PROC_DECL(SET_XHYDROGEN, set_xhydrogen)
    (Real& xhydrogen_in);

BL_FORT_PROC_DECL(FORT_INTEGRATE_COMOVING_A, fort_integrate_comoving_a)
    (Real* old_a, Real* new_a, Real* dt);

BL_FORT_PROC_DECL(FORT_INTEGRATE_COMOVING_A_TO_Z, fort_integrate_comoving_a_to_z)
    (Real* old_a, Real* z_value, Real* dt);

BL_FORT_PROC_DECL(FORT_ESTDT, fort_estdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const Real dx[], Real* dt, Real* comoving_a);

BL_FORT_PROC_DECL(FORT_ESTDT_COMOVING_A, fort_estdt_comoving_a)
    (Real* old_a, Real* new_dummy_a, Real* dt, Real* change_allowed,
     Real* final_a, int* dt_modified);

BL_FORT_PROC_DECL(GET_NUM_SPEC, get_num_spec)(int* nspec);

BL_FORT_PROC_DECL(GET_NUM_AUX, get_num_aux)(int* naux);

BL_FORT_PROC_DECL(GET_SPEC_NAMES, get_spec_names)
    (int* spec_names, int* ispec, int* len);

BL_FORT_PROC_DECL(GET_AUX_NAMES, get_aux_names)
    (int* aux_names, int* iaux, int* len);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS, get_method_params)(int* HYP_GROW);

BL_FORT_PROC_DECL(SET_METHOD_PARAMS, set_method_params)
    (const int& dm, const int& NumAdv, const int& do_hydro,
     const int& ppm_type, const int& ppm_ref,
     const int& ppm_flatten_before_integrals,
     const int& use_colglaz, const int& use_flattening, 
     const int& corner_coupling, 
     const int& version_2, const int& use_const_species,
     const Real& gamma_in, const int& normalize_species,
     const int& heat_cool_type); 

BL_FORT_PROC_DECL(SET_EOS_PARAMS, set_eos_params)
     (const Real& h_species_in, const Real& he_species_in);

BL_FORT_PROC_DECL(SET_SMALL_VALUES, set_small_values)
    (const Real* average_dens, const Real* average_temp, 
     const Real* comoving_a, 
     const Real* small_dens, 
     const Real* small_temp, const Real* small_pres);

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS, set_problem_params)
    (const int& dm, const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, const int& Symmetry_value,
     const int& coord_type);

BL_FORT_PROC_DECL(CA_INITDATA, ca_initdata)
    (const int& level, const Real& time, const int* lo, const int* hi,
     const int& num_state, BL_FORT_FAB_ARG(state), 
     const int& num_diag,  BL_FORT_FAB_ARG(diag_eos), const Real dx[],
     const Real xlo[], const Real xhi[]);

BL_FORT_PROC_DECL(FORT_CHECK_INITIAL_SPECIES, fort_check_initial_species)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(INIT_E_FROM_T, init_e_from_t)
     (BL_FORT_FAB_ARG(state), const int* num_state,
      BL_FORT_FAB_ARG( diag), const int* num_diag, 
      const int* lo, const int* hi, Real* comoving_a);

BL_FORT_PROC_DECL(INIT_E_FROM_RHOE, init_e_from_rhoe)
     (BL_FORT_FAB_ARG(state), const int* num_state,
      const int* lo, const int* hi, Real* comoving_a);

BL_FORT_PROC_DECL(GET_OMB   , get_omb   )(Real* frac);
BL_FORT_PROC_DECL(GET_OMM   , get_omm   )(Real* omm);
BL_FORT_PROC_DECL(GET_HUBBLE, get_hubble)(Real* hubble);

BL_FORT_PROC_DECL(ENFORCE_CONSISTENT_E, enforce_consistent_e)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(HYPFILL, hypfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(DENFILL, denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(DENFILL, denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(GENERIC_FILL, generic_fill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(FORT_ADVANCE_GAS, fort_advance_gas)
    (const Real* time,
     const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state), BL_FORT_FAB_ARG(stateout),
     const BL_FORT_FAB_ARG(ugdnvx),
     const BL_FORT_FAB_ARG(ugdnvy),
     const BL_FORT_FAB_ARG(ugdnvz),
     const BL_FORT_FAB_ARG(src),
     const BL_FORT_FAB_ARG(grav),
     const Real dx[], const Real* dt,
     D_DECL(const BL_FORT_FAB_ARG(xflux),
            const BL_FORT_FAB_ARG(yflux),
            const BL_FORT_FAB_ARG(zflux)),
     const Real* cflLoc,
     const Real* a_old, const Real* a_new,
     const Real* e_added, const Real* ke_added, 
     const int* print_fortran_warnings,
     const int* do_gas);

BL_FORT_PROC_DECL(TIME_CENTER_SOURCES, time_center_sources)
    (const int lo[], const int hi[], BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const Real* a_old, const Real* a_new,
     const Real* dt, const int* print_fortran_warnings);

BL_FORT_PROC_DECL(ADJUST_HEAT_COOL, adjust_heat_cool)
    (const int lo[], const int hi[], 
     BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const Real* a_old, const Real* a_new,
     const Real* dt);

BL_FORT_PROC_DECL(CORRECT_GSRC, correct_gsrc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(grav_src_old),
     const BL_FORT_FAB_ARG(grad_phi_cc),
     const BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     const Real* a_old, const Real* a_new,
     const Real* dt, const Real* e_added, const Real* ke_added);

BL_FORT_PROC_DECL(RESET_INTERNAL_E, reset_internal_e)
    (BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(D_new),
     const int lo[], const int hi[],
     const int* print_fortran_warnings, 
     Real* comoving_a, Real* sum_energy_added,
     Real* sum_energy_total);

BL_FORT_PROC_DECL(FORT_SYNCGSRC, fort_syncgsrc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(gPhi),
     const BL_FORT_FAB_ARG(gdPhi),
     const BL_FORT_FAB_ARG(S),
     const BL_FORT_FAB_ARG(dS),
     BL_FORT_FAB_ARG(src),
     const Real* a_new,
     const Real& dt);

BL_FORT_PROC_DECL(SUM_OVER_LEVEL, sum_over_level)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* sum);

BL_FORT_PROC_DECL(SUM_PRODUCT, sum_product)
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2), 
     const int lo[], const int hi[], const Real dx[], Real* sum);

BL_FORT_PROC_DECL(SUM_PROD_PROD, sum_prod_prod)
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2), 
     BL_FORT_FAB_ARG(fab3), 
     const int lo[], const int hi[], const Real dx[], Real* sum);

BL_FORT_PROC_DECL(ENFORCE_NONNEGATIVE_SPECIES,
                  enforce_nonnegative_species)
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

BL_FORT_PROC_DECL(FORT_AVGDOWN, fort_avgdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);

#ifdef AGN
BL_FORT_PROC_DECL(CA_EXT_SRC, ca_ext_src)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const Real* particle_locs_and_mass,
     const Real* particle_data,
     const Real* prob_lo, const Real* dx,
     const Real* time, const Real* z, 
     const Real* dt);
#else
BL_FORT_PROC_DECL(CA_EXT_SRC, ca_ext_src)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const Real* prob_lo, const Real* dx,
     const Real* time, const Real* z, 
     const Real* dt);
#endif

BL_FORT_PROC_DECL(INTEGRATE_STATE, integrate_state)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(diag_eos),
     const Real* z, const Real* dt,
     const int* min_iter, const int* max_iter);

BL_FORT_PROC_DECL(COMPUTE_TEMP, compute_temp)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     Real* comoving_a,
     const int* print_fortran_warnings);

BL_FORT_PROC_DECL(INIT_THIS_Z, init_this_z)
    (Real* comoving_a);

BL_FORT_PROC_DECL(COMPUTE_MAX_TEMP_LOC, compute_max_temp_loc)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     const Real* max_temp, const Real* den_maxt,
     const int* imax, const int* jmax, const int* kmax);

BL_FORT_PROC_DECL(COMPUTE_RHO_TEMP, compute_rho_temp)
    (const int lo[], const int hi[], const Real dx[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     Real* rho_ave, Real* rho_T_sum, Real* T_sum,
     Real* T_meanrho_sum, Real* rho_sum,
     Real* vol_sum, Real* vol_mn_sum);

#ifdef AUX_UPDATE
BL_FORT_PROC_DECL(AUXUPDATE, auxupdate)
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const Real * dt);
#endif

BL_FORT_PROC_DECL(GET_RHOE, ca_get_rhoe)
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));
#endif
