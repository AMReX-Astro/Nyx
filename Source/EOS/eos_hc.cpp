
#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;
/*
#ifdef __cplusplus
extern "C"
{
#endif

  AMREX_GPU_DEVICE void fort_nyx_eos_T_given_Re_device(int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a);
  AMREX_GPU_DEVICE void fort_nyx_eos_given_RT(Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a);

#ifdef __cplusplus
}
#endif
*/


AMREX_GPU_DEVICE void Nyx::nyx_eos_T_given_Re_device(Real gamma_minus_1, Real h_species, int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a)
{

   //   Real species[5];
   //   Real* species_ptr=&species[0];
   fort_nyx_eos_T_given_Re_device(JH,JHe,T,Ne,R,e,comoving_a);//,species_ptr);

}

AMREX_GPU_DEVICE void Nyx::nyx_eos_given_RT(Real gamma_minus_1, Real h_species, Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a)
{
  
  Real mu = (1.0+4.0*YHELIUM) / (1.0+YHELIUM+Ne);
  *e = T / (gamma_minus_1 * mp_over_kb * mu);
  *P  = gamma_minus_1 * (R) * (*e);
  //  fort_nyx_eos_given_RT(e,P,R,T,Ne,comoving_a);
  ////  printf("etmp: %g Ptmp: %g\ne:    %g P:    %g\n",e_tmp,P_tmp,*e,*P);
	 //  if((*e!=e_tmp)&&(*P!=P_tmp))
  ////    amrex::Abort("found diff");//+e_tmp+P_tmp+(*e)+(*P));

}
