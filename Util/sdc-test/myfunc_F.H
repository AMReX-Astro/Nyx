#ifndef MYFUNC_F_H_
#define MYFUNC_F_H_

#ifdef __cplusplus
extern "C"
{
#endif


  void FSetIC(amrex_real* phi,
                     const int* lo,
                     const int* hi,
	             double* dptr);

  void FSetIC_mfab(amrex_real* phi,
                     const int* lo,
		   const int* hi);

  void FSetInternalEnergy_mfab(amrex_real* e,
                     const int* lo,
		   const int* hi);

  void FSetSol(amrex_real* mf,
                     const int* lo,
                     const int* hi,
	             double* dptr);


  AMREX_GPU_DEVICE int RhsFnReal(amrex::Real t, amrex::Real* u, amrex::Real* udot, amrex::Real* rpar, int neq);
  AMREX_GPU_DEVICE void fort_ode_eos_finalize(double* e_out, double* rpar, int neq);
  void fort_ode_eos_setup(const amrex::Real& a,const amrex::Real& half_dt);

  void fort_setup_eos_params
  (amrex::Real* eos_nr_eps,
   amrex::Real* vode_rtol,
   amrex::Real* vode_atol_scaled);

  void fort_alloc_cuda_managed();

  void fort_init_tables_eos_params();

  void fort_dealloc_cuda_managed();

    void fort_fab_copytoreal(const int* lo, const int* hi,
                                  const int* bx_lo, const int* bx_hi,
                                  amrex_real* dst,
                                  const amrex_real* src, const int* slo, const int* shi,
                                  const int ncomp);

    void fort_fab_copyfromreal(const int* lo, const int* hi,
                                    const int* bx_lo, const int* bx_hi,
                                    const amrex_real* dst, const int* dlo, const int* dhi, const int ncomp,
                                    const amrex_real* src);

  /*  void FSetIC_mfab(amrex_real* mf,
                     const int* lo,
		   const int* hi);//, const int* ncomp, const int* src_comp);
  */

#ifdef __cplusplus
}
#endif

#endif
