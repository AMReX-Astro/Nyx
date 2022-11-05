#ifndef __f_rhs_struct__
#define __f_rhs_struct__


#include <atomic_rates_data.H>
#include <constants_cosmo.H>
#include <eos_hc.H>
#include <IndexDefines.H>

struct RhsData {
    // Set once in init
    amrex::Real uvb_density_A;
    amrex::Real uvb_density_B;
    int         flash_h;
    int         flash_he;
    amrex::Real zhi_flash;
    amrex::Real zheii_flash;
    amrex::Real T_zhi;
    amrex::Real T_zheii;
    amrex::Real H_reion_z;
    amrex::Real He_reion_z;
    int         inhomogeneous_on;

    // Filled from Nyx.cpp
    amrex::Real gamma_minus_1;
    amrex::Real h_species;

    // Change at most twice per Nyx coarse time-step
    amrex::Real a;
    int         JH_vode;
    int*        JH_vode_arr;
    int         JHe_vode;
    amrex::Real* T_vode;
    amrex::Real* ne_vode;
    amrex::Real* rho_vode;
    amrex::Real* rho_init_vode;
    amrex::Real* rho_src_vode;
    amrex::Real* rhoe_src_vode;
    amrex::Real* e_src_vode;
    amrex::Real* IR;
    Dim3 lo;
    Dim3 len;
};

void
ode_eos_setup(RhsData* f_rhs_data, amrex::Real gamma_minus_1, amrex::Real h_species)
{

    RhsData f_rhs_data_host;
    f_rhs_data_host.gamma_minus_1 = gamma_minus_1;
    f_rhs_data_host.h_species = h_species;
    f_rhs_data_host.inhomogeneous_on=0;
    ParmParse pp_nyx("nyx");
    pp_nyx.query("inhomo_reion", f_rhs_data_host.inhomogeneous_on);

    f_rhs_data_host.uvb_density_A=1.0;
    f_rhs_data_host.uvb_density_B=0.0;
    f_rhs_data_host.zhi_flash=-1.0;
    f_rhs_data_host.zheii_flash=-1.0;
    f_rhs_data_host.T_zhi=0.0;
    f_rhs_data_host.T_zheii=0.0;
    f_rhs_data_host.flash_h=0;
    f_rhs_data_host.flash_he=0;
    pp_nyx.query("uvb_density_A"            , f_rhs_data_host.uvb_density_A);
    pp_nyx.query("uvb_density_B"            , f_rhs_data_host.uvb_density_B);
    pp_nyx.query("reionization_zHI_flash"   , f_rhs_data_host.zhi_flash);
    pp_nyx.query("reionization_zHeII_flash" , f_rhs_data_host.zheii_flash);
    pp_nyx.query("reionization_T_zHI"       , f_rhs_data_host.T_zhi);
    pp_nyx.query("reionization_T_zHeII"     , f_rhs_data_host.T_zheii);

    // Set options in reion_aux_module                                                                                                                                                                                                        
    //   Hydrogen reionization                                                                                                                                                                                                                
    if (f_rhs_data_host.zhi_flash > 0.0)
    {
        if (f_rhs_data_host.inhomogeneous_on > 0)
        {
            f_rhs_data_host.flash_h = 0;
        }
        else
        {
            f_rhs_data_host.flash_h = 1;
        }
    }
    else
    {
        f_rhs_data_host.flash_h = 0;
    }

    //   Helium reionization                                                                                                                                                                                                                  
    if (f_rhs_data_host.zheii_flash > 0.0)
        f_rhs_data_host.flash_he = 1;
    else
        f_rhs_data_host.flash_he = 0;

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(f_rhs_data,&f_rhs_data_host,sizeof(RhsData));
#else
    std::memcpy(f_rhs_data,&f_rhs_data_host,sizeof(RhsData));
#endif

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_get_idx(RhsData* f_rhs_data, const int i, const int j, const int k, long int& idx)
{
    idx = i+j*f_rhs_data->len.x+k*f_rhs_data->len.x*f_rhs_data->len.y-(f_rhs_data->lo.x+f_rhs_data->lo.y*f_rhs_data->len.x+f_rhs_data->lo.z*f_rhs_data->len.x*f_rhs_data->len.y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_initialize_single(RhsData* f_rhs_data,
                          const amrex::Real a,
                          amrex::Real* dptr, amrex::Real* eptr,
                          amrex::Real* T_vode,
                          amrex::Real* ne_vode,
                          amrex::Real* rho_vode,
                          amrex::Real* rho_init_vode,
                          amrex::Real* rho_src_vode,
                          amrex::Real* rhoe_src_vode,
                          amrex::Real* e_src_vode,
                          amrex::Real* IR,
                          int*         JH_vode_arr)
{

    f_rhs_data->a = a;
    amrex::Real z = 1/(f_rhs_data->a)-1.0;
    // Flash reionization?
    if ((f_rhs_data->flash_h == 1) && (z > f_rhs_data->zhi_flash))
       f_rhs_data->JH_vode = 0;
    else
       f_rhs_data->JH_vode = 1;

    if ((f_rhs_data->flash_he == 1) && (z > f_rhs_data->zheii_flash))
        f_rhs_data->JHe_vode = 0;
    else
        f_rhs_data->JHe_vode = 1;

    if (f_rhs_data->flash_h ) f_rhs_data->H_reion_z  = f_rhs_data->zhi_flash;
    if (f_rhs_data->flash_he) f_rhs_data->He_reion_z = f_rhs_data->zheii_flash;

    f_rhs_data->T_vode = T_vode;
    f_rhs_data->ne_vode = ne_vode;
    f_rhs_data->rho_vode = rho_vode;
    f_rhs_data->rho_init_vode = rho_init_vode;
    f_rhs_data->rho_src_vode = rho_src_vode;
    f_rhs_data->rhoe_src_vode = rhoe_src_vode;
    f_rhs_data->e_src_vode = e_src_vode;
    f_rhs_data->IR = IR;
    f_rhs_data->JH_vode_arr = JH_vode_arr;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_initialize_arrays(const int i, const int j, const int k, const int idx, RhsData* f_rhs_data,
                          const amrex::Real a_end, amrex::Dim3 lo, amrex::Dim3 len,
                          amrex::Array4<Real> const& state4,
                          amrex::Array4<Real> const& diag_eos4,
                          amrex::Array4<Real> const& hydro_src,
                          amrex::Array4<Real> const& reset_src,
                          amrex::Real* dptr, amrex::Real* eptr, amrex::Real* abstol_ptr,
                          const int sdc_iter, const amrex::Real delta_time)
{
    amrex::Real asq, aendsq, z, z_end;
    const amrex::Real abstol = 1e-4;

    z = 1.0/(f_rhs_data->a) - 1.0;
    z_end = 1.0/a_end - 1.0;
 
    asq = f_rhs_data->a*f_rhs_data->a;
    aendsq = a_end*a_end;

    f_rhs_data->lo = lo;
    f_rhs_data->len = len;

    //  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
    dptr[idx]=state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp);
    eptr[idx]=state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp);
    f_rhs_data->T_vode[idx]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
    f_rhs_data->ne_vode[idx]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
    f_rhs_data->rho_vode[idx]= state4(i,j,k,Density_comp); //    rpar(3)=rho_vode
    abstol_ptr[idx]= state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp)*abstol;
    if(sdc_iter >=0)
    {
    f_rhs_data->rho_src_vode[idx]= hydro_src(i,j,k,Density_comp) / delta_time;
    f_rhs_data->rhoe_src_vode[idx]= hydro_src(i,j,k,Eint_comp) / delta_time;
    f_rhs_data->e_src_vode[idx] = ( ((asq*state4(i,j,k,Eint_comp) + delta_time * f_rhs_data->rhoe_src_vode[idx] ) / aendsq + reset_src(i,j,k,0))/ 
                          (      state4(i,j,k,Density_comp ) + delta_time * f_rhs_data->rho_src_vode[idx]  ) - eptr[idx]) / delta_time;
    f_rhs_data->rho_init_vode[idx]= state4(i,j,k,Density_comp); //    rpar(3)=rho_vode
    }
    else
    {
    f_rhs_data->rho_src_vode=NULL;
    f_rhs_data->rhoe_src_vode=NULL;
    f_rhs_data->e_src_vode=NULL;
    f_rhs_data->rho_init_vode=NULL;
    }

    if (f_rhs_data->inhomogeneous_on)
    {
        f_rhs_data->H_reion_z = diag_eos4(i,j,k,Zhi_comp);
        if (z > f_rhs_data->H_reion_z)
            f_rhs_data->JH_vode_arr[idx] = 0;
        else
            f_rhs_data->JH_vode_arr[idx] = 1;
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_save_react_arrays(const int i, const int j, const int k, const int idx, AtomicRates* atomic_rates, RhsData* f_rhs_data,
                          const amrex::Real a_end,
			  amrex::Array4<const Real> const& state4, amrex::Array4<const Real> const& state_n4,
                          amrex::Array4<const Real> const& reset_src, amrex::Array4<const Real> const& diag_eos4,
			  amrex::Array4<const Real> const& I_R,
                          amrex::Array4<Real> const& react_in_arr,
                          amrex::Array4<Real> const& react_out_arr,
                          amrex::Array4<Real> const& react_out_work_arr,
                          amrex::Real* dptr, amrex::Real* eptr, amrex::Real* abstol_ptr, amrex::Real* abstol_achieve_ptr,
                          const int sdc_iter, const amrex::Real delta_time, const amrex::Real time_in,
			  long int nst, long int netf, long int nfe,
			  long int nni, long int ncfn, long int nsetups, long int nje, long int ncfl, long int nfeLS)
{
    /*
    e_out = dptr[idx];
    e_orig = eptr[idx];
    rho_out = state_n4(i,j,k,Density_comp);
    rho_orig = f_rhs_data->rho_init_vode[idx];
    rho_src = f_rhs_data->rho_src_vode[idx]; // used to construct rho during ode evolve, where time is given by sundials
    //rho_vode=f_rhs_data->rho_init_vode[idx]+time*f_rhs_data->rho_src_vode[idx];
    rhoe_src = f_rhs_data->rhoe_src_vode[idx]; // used to construct IR, this is A_{rhoe}
    e_src = f_rhs_data->e_src_vode[idx]; // used as forcing in ode solve, where lambda_HC is modified with A_{e}
    //energy  = energy / rho_vode * (1.0e0+amrex::Math::abs(z_vode)) + f_rhs_data->e_src_vode[idx];
    //if called immediately after ode_eos_finalize_struct these are outputs, if called after ode_eos_initialize_struct these are inputs
    T_out = f_rhs_data->T_vode[idx];
    ne_out = f_rhs_data->ne_vode[idx];
    abstol_request = abstol_ptr[idx];
    abstol_achieve_est = abstol_achieve_ptr[idx];
    */
    int comp1=0;
    int comp2=0;
    int comp3=0;
    react_in_arr(i,j,k,comp1++) = eptr[idx];
    react_out_arr(i,j,k,comp2++) = dptr[idx];
    react_in_arr(i,j,k,comp1++) = f_rhs_data->rho_init_vode[idx];
    react_out_arr(i,j,k,comp2++) = f_rhs_data->rho_init_vode[idx] + delta_time*f_rhs_data->rho_src_vode[idx];
    react_in_arr(i,j,k,comp1++) = f_rhs_data->rhoe_src_vode[idx];
    react_in_arr(i,j,k,comp1++) = f_rhs_data->e_src_vode[idx];
    react_out_arr(i,j,k,comp2++) = f_rhs_data->T_vode[idx];
    react_out_arr(i,j,k,comp2++) = f_rhs_data->ne_vode[idx];
    react_in_arr(i,j,k,comp1++) = abstol_ptr[idx];
    react_out_arr(i,j,k,comp2++) = abstol_achieve_ptr[idx];
    react_in_arr(i,j,k,comp1++) = f_rhs_data->a;
    react_out_arr(i,j,k,comp2++) = a_end;
    react_in_arr(i,j,k,comp1++) = time_in;
    react_out_arr(i,j,k,comp2++) = delta_time;
    react_out_work_arr(i,j,k,comp3++) = nst;
    react_out_work_arr(i,j,k,comp3++) = netf;
    react_out_work_arr(i,j,k,comp3++) = nfe;
    react_out_work_arr(i,j,k,comp3++) = nni;
    react_out_work_arr(i,j,k,comp3++) = ncfn;
    react_out_work_arr(i,j,k,comp3++) = nsetups;
    react_out_work_arr(i,j,k,comp3++) = nje;
    react_out_work_arr(i,j,k,comp3++) = ncfl;
    react_out_work_arr(i,j,k,comp3++) = nfeLS;

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_finalize_struct(const int i, const int j, const int k, const int idx, AtomicRates* atomic_rates, RhsData* f_rhs_data,
                        const amrex::Real a_end, amrex::Array4<Real> const& state4, amrex::Array4<Real> const& state_n4,
                        amrex::Array4<Real> const& reset_src, amrex::Array4<Real> const& diag_eos4, amrex::Array4<Real> const& I_R,
                        amrex::Real* dptr, amrex::Real* eptr, const amrex::Real delta_time)
{
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;
    amrex::Real rho, U, a, mu, e_out, e_orig, rho_out, rho_orig, asq, aendsq, ahalf, ahalf_inv;
    e_out = dptr[idx];
    e_orig = eptr[idx];
    rho_out = state_n4(i,j,k,Density_comp);
    if (f_rhs_data->inhomogeneous_on)
        f_rhs_data->H_reion_z = diag_eos4(i,j,k,Zhi_comp);
    //  int  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
    /////////////////   //update nyx_eos_T_given_Re
    diag_eos4(i,j,k,Temp_comp)=f_rhs_data->T_vode[idx];   //rpar(1)=T_vode
    diag_eos4(i,j,k,Ne_comp)=f_rhs_data->ne_vode[idx];//    rpar(2)=ne_vode
#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
    {
        a = f_rhs_data->a;
        asq = a*a;
        aendsq = a_end*a_end;
        ahalf     = 0.5 * (a + a_end);
        ahalf_inv  = 1.0 / ahalf;
        rho_orig = f_rhs_data->rho_init_vode[idx];
        //                             
        // I_R satisfies the equation anewsq * (rho_out  e_out ) = 
        //                            AOLDSQ * (rho_orig e_orig) + dt * a_half * I_R + dt * H_{rho e} + anewsq * reset_src
        I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*f_rhs_data->rhoe_src_vode[idx]) ) ) / (delta_time * ahalf) - aendsq * reset_src(i,j,k,0) / (delta_time * ahalf);
    }

    T_vode=f_rhs_data->T_vode[idx];
    ne_vode=f_rhs_data->ne_vode[idx];
    rho_vode=f_rhs_data->rho_vode[idx];
    z_vode=1/(f_rhs_data->a)-1.0;

#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
    {
        if ((state_n4(i,j,k,Eint_comp) + delta_time * ahalf * I_R(i,j,k) / aendsq) / state_n4(i,j,k,Density_comp) < 0.e0)
        {
            Real YHELIUM = (1.0-f_rhs_data->h_species)/(4.0*f_rhs_data->h_species);
            f_rhs_data->T_vode[idx]  = 10.0;
            f_rhs_data->ne_vode[idx] = 0.0;
            mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+f_rhs_data->ne_vode[idx]);
            e_out  = f_rhs_data->T_vode[idx] / (f_rhs_data->gamma_minus_1 * mp_over_kb * mu);
            I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*f_rhs_data->rhoe_src_vode[idx]) ) ) / (delta_time * ahalf) - aendsq * reset_src(i,j,k,0) / (delta_time * ahalf);
        }
    }
    else
    {
    if(e_out < 0.e0)
    {
        Real YHELIUM = (1.0-f_rhs_data->h_species)/(4.0*f_rhs_data->h_species);
        f_rhs_data->T_vode[idx]  = 10.0;
        f_rhs_data->ne_vode[idx] = 0.0;
        mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+f_rhs_data->ne_vode[idx]);
        e_out  = f_rhs_data->T_vode[idx] / (f_rhs_data->gamma_minus_1 * mp_over_kb * mu);
    }
    }

    Real nh0, nhep, nhp, nhe0, nhepp;
    /////////////////   //update nyx_eos_T_given_Re
    if (f_rhs_data->inhomogeneous_on)
      nyx_eos_T_given_Re_device(atomic_rates, f_rhs_data->gamma_minus_1, f_rhs_data->h_species, f_rhs_data->JH_vode_arr[idx], f_rhs_data->JHe_vode, &(f_rhs_data->T_vode[idx]), &(f_rhs_data->ne_vode[idx]), f_rhs_data->rho_vode[idx], e_out, f_rhs_data->a, nh0, nhep, nhp, nhe0, nhepp);
    else
      nyx_eos_T_given_Re_device(atomic_rates, f_rhs_data->gamma_minus_1, f_rhs_data->h_species, f_rhs_data->JH_vode, f_rhs_data->JHe_vode, &(f_rhs_data->T_vode[idx]), &(f_rhs_data->ne_vode[idx]), f_rhs_data->rho_vode[idx], e_out, f_rhs_data->a, nh0, nhep, nhp, nhe0, nhepp);

    // Instanteneous heating from reionization:
    Real z = 1/(f_rhs_data->a)-1.0;
    Real z_end = 1/(a_end)-1.0;

    Real T_H = 0.0e0;
    if (f_rhs_data->inhomogeneous_on || f_rhs_data->flash_h)
    {
      if ((f_rhs_data->H_reion_z  < z) && (f_rhs_data->H_reion_z  >= z_end))
        T_H  = (1.0e0 - nhp)*amrex::max((f_rhs_data->T_zhi-f_rhs_data->T_vode[idx]), 0.0e0);
    }
    Real T_He = 0.0e0;
    if (f_rhs_data->flash_he)
    {
      if ((f_rhs_data->He_reion_z < z) && (f_rhs_data->He_reion_z >= z_end)) 
        {
        T_He = (1.0e0 - nhepp)*amrex::max((f_rhs_data->T_zheii-f_rhs_data->T_vode[idx]), 0.0e0);
        }
    }

    if ((T_H > 0.0e0) || (T_He > 0.0e0))
    {
        Real YHELIUM = (1.0-f_rhs_data->h_species)/(4.0*f_rhs_data->h_species);
        f_rhs_data->T_vode[idx] = f_rhs_data->T_vode[idx] + T_H + T_He;                     // For simplicity, we assume
        f_rhs_data->ne_vode[idx] = 1.0e0 + YHELIUM;                                         //    completely ionized medium at
        if (T_He > 0.0e0) f_rhs_data->ne_vode[idx] = f_rhs_data->ne_vode[idx] + YHELIUM;    //    this point.  It's a very minor
        mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+f_rhs_data->ne_vode[idx]);          //    detail compared to the overall approximation.
        e_out  = f_rhs_data->T_vode[idx] / (f_rhs_data->gamma_minus_1 * mp_over_kb * mu);
#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
    {
        //
        // I_R satisfies the equation anewsq * (rho_out  e_out ) =
        //                            AOLDSQ * (rho_orig e_orig) + dt * a_half * I_R + dt * H_{rho e} + anewsq * reset_src
        I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*f_rhs_data->rhoe_src_vode[idx]) ) ) / (delta_time * ahalf) - aendsq * reset_src(i,j,k,0) / (delta_time * ahalf);
    }

    //Check for negative e_out
#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
    {
        if ((state_n4(i,j,k,Eint_comp) + delta_time * ahalf * I_R(i,j,k) / aendsq) / state_n4(i,j,k,Density_comp) < 0.e0)
        {
            Real YHELIUM = (1.0-f_rhs_data->h_species)/(4.0*f_rhs_data->h_species);
            f_rhs_data->T_vode[idx]  = 10.0;
            f_rhs_data->ne_vode[idx] = 0.0;
            mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+f_rhs_data->ne_vode[idx]);
            e_out  = f_rhs_data->T_vode[idx] / (f_rhs_data->gamma_minus_1 * mp_over_kb * mu);
            I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*f_rhs_data->rhoe_src_vode[idx]) ) ) / (delta_time * ahalf) - aendsq * reset_src(i,j,k,0) / (delta_time * ahalf);
        }
    }
    else
    {
    if(e_out < 0.e0)
    {
        Real YHELIUM = (1.0-f_rhs_data->h_species)/(4.0*f_rhs_data->h_species);
        f_rhs_data->T_vode[idx]  = 10.0;
        f_rhs_data->ne_vode[idx] = 0.0;
        mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+f_rhs_data->ne_vode[idx]);
        e_out  = f_rhs_data->T_vode[idx] / (f_rhs_data->gamma_minus_1 * mp_over_kb * mu);
    }
    }
    if (f_rhs_data->inhomogeneous_on)
      nyx_eos_T_given_Re_device(atomic_rates, f_rhs_data->gamma_minus_1, f_rhs_data->h_species, f_rhs_data->JH_vode_arr[idx], f_rhs_data->JHe_vode, &(f_rhs_data->T_vode[idx]), &(f_rhs_data->ne_vode[idx]), f_rhs_data->rho_vode[idx], e_out, f_rhs_data->a);
    else
      nyx_eos_T_given_Re_device(atomic_rates, f_rhs_data->gamma_minus_1, f_rhs_data->h_species, f_rhs_data->JH_vode, f_rhs_data->JHe_vode, &(f_rhs_data->T_vode[idx]), &(f_rhs_data->ne_vode[idx]), f_rhs_data->rho_vode[idx], e_out, f_rhs_data->a);
    }


#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
    {
        // Update (rho e) and (rho E)
        // Note that we add to state_n because those already have hydro_source in them
        state_n4(i,j,k,Eint_comp) = state_n4(i,j,k,Eint_comp) + delta_time * ahalf * I_R(i,j,k) / aendsq;
        state_n4(i,j,k,Eden_comp) = state_n4(i,j,k,Eden_comp) + delta_time * ahalf * I_R(i,j,k) / aendsq;
    }
    else
    {
        state4(i,j,k,Eint_comp)  += state4(i,j,k,Density_comp) * (e_out-e_orig);
        state4(i,j,k,Eden_comp)  += state4(i,j,k,Density_comp) * (e_out-e_orig);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int
f_rhs_struct(amrex::Real const time, amrex::Real& e_in, amrex::Real& energy, AtomicRates* atomic_rates, RhsData* f_rhs_data, const int idx)
{

    const amrex::Real compt_c = 1.01765467e-37;
    const amrex::Real T_cmb = 2.725e0;
    amrex::Real deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB;

    amrex::Real logT, tmp, fhi, flo;
    amrex::Real ggh0, gghe0, gghep, eh0, ehe0, ehep;
    amrex::Real ahp, ahep, ahepp, ad, geh0, gehe0, gehep;
    amrex::Real bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp;
    amrex::Real lambda_c, lambda_ff, lambda, heat;
    amrex::Real rho, U, a, rho_heat;
    amrex::Real nh, nh0, nhp, nhe0, nhep, nhepp;
    int j;
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;
    amrex::Real gamma_minus_1;

    T_vode=f_rhs_data->T_vode[idx];
    ne_vode=f_rhs_data->ne_vode[idx];
#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
        rho_vode=f_rhs_data->rho_init_vode[idx]+time*f_rhs_data->rho_src_vode[idx];
    else
        rho_vode=f_rhs_data->rho_vode[idx];
    z_vode=1/(f_rhs_data->a)-1.0;;
    a=f_rhs_data->a;
    gamma_minus_1=f_rhs_data->gamma_minus_1;
    if(e_in<=0 || amrex::isnan(e_in))     e_in=std::numeric_limits<amrex::Real>::min();

    // Converts from code units to CGS
    rho = rho_vode * density_to_cgs * (1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode));
    U = e_in * e_to_cgs;
    nh  = rho*f_rhs_data->h_species/MPROTON;

    // Get gas temperature and individual ionization species
    // testing different memory structures
    if (f_rhs_data->inhomogeneous_on)
        iterate_ne_device(atomic_rates, f_rhs_data->JH_vode_arr[idx], f_rhs_data->JHe_vode, z_vode, U, T_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1, f_rhs_data->h_species);
    else
        iterate_ne_device(atomic_rates, f_rhs_data->JH_vode, f_rhs_data->JHe_vode, z_vode, U, T_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1, f_rhs_data->h_species);
    // Convert species to CGS units: 
    ne_vode = nh * ne_vode;
    nh0   = nh * nh0;
    nhp   = nh * nhp;
    nhe0  = nh * nhe0;
    nhep  = nh * nhep;
    nhepp = nh * nhepp;

    logT = std::log10(T_vode);
    if (logT >= TCOOLMAX) // Only free-free and Compton cooling are relevant
    {
        lambda_ff = 1.42e-27 * std::sqrt(T_vode) * (1.1e0 + 0.34e0*std::exp(-(5.5e0 - logT)*(5.5e0 - logT) / 3.0e0))
            * (nhp + 4.0e0*nhepp)*ne_vode;
        lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));

        energy  = (-lambda_ff -lambda_c) * heat_from_cgs/((1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode)));

        // Convert to the actual term to be used in e_out = e_in + dt*energy
#ifdef AMREX_USE_HIP
        if((f_rhs_data->e_src_vode!=NULL)||true)
#else
        if((f_rhs_data->e_src_vode!=NULL))
#endif
               energy  = energy / rho_vode * (1.0e0+amrex::Math::abs(z_vode)) + f_rhs_data->e_src_vode[idx];
        else
               energy  = energy / rho_vode * (1.0e0+amrex::Math::abs(z_vode));
        ne_vode = ne_vode / nh;
        f_rhs_data->T_vode[idx]=T_vode;
        f_rhs_data->ne_vode[idx]=ne_vode;
        f_rhs_data->rho_vode[idx]=rho_vode;
        return 0;
    }

    // Temperature floor
    if (logT <= TCOOLMIN)  logT = TCOOLMIN + 0.5e0*deltaT;

    // Interpolate rates
    // Interpolate rates
    tmp = (logT-TCOOLMIN)/deltaT;
    j = amrex::Math::floor(tmp); // replacing int(tmp)
    fhi = tmp - j;
    flo = 1.0e0 - fhi;
    j = j + 1; // F90 arrays start with 1, using Array1D<Real,1,size>

    bh0   = flo*atomic_rates->BetaH0   (j) + fhi*atomic_rates->BetaH0   (j+1);
    bhe0  = flo*atomic_rates->BetaHe0  (j) + fhi*atomic_rates->BetaHe0  (j+1);
    bhep  = flo*atomic_rates->BetaHep  (j) + fhi*atomic_rates->BetaHep  (j+1);
    bff1  = flo*atomic_rates->Betaff1  (j) + fhi*atomic_rates->Betaff1  (j+1);
    bff4  = flo*atomic_rates->Betaff4  (j) + fhi*atomic_rates->Betaff4  (j+1);
    rhp   = flo*atomic_rates->RecHp    (j) + fhi*atomic_rates->RecHp    (j+1);
    rhep  = flo*atomic_rates->RecHep   (j) + fhi*atomic_rates->RecHep   (j+1);
    rhepp = flo*atomic_rates->RecHepp  (j) + fhi*atomic_rates->RecHepp  (j+1);

    // Cooling: 
    lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + 
               rhp*nhp + rhep*nhep + rhepp*nhepp + 
               bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode;

    lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));
    lambda = lambda + lambda_c;

    interp_to_this_z(atomic_rates, z_vode, ggh0, gghe0, gghep, eh0, ehe0, ehep);
    // Heating terms
    if (f_rhs_data->inhomogeneous_on)
        heat = f_rhs_data->JH_vode_arr[idx]*nh0*eh0 + f_rhs_data->JH_vode_arr[idx]*nhe0*ehe0 + f_rhs_data->JHe_vode*nhep*ehep;
    else
        heat = f_rhs_data->JH_vode*nh0*eh0 + f_rhs_data->JH_vode*nhe0*ehe0 + f_rhs_data->JHe_vode*nhep*ehep;

        rho_heat = f_rhs_data->uvb_density_A * std::pow((rho_vode/atomic_rates->mean_rhob),f_rhs_data->uvb_density_B);
    heat = rho_heat*heat;

    // Convert back to code units
    ne_vode     = ne_vode / nh;
    energy  = (heat - lambda)*heat_from_cgs / ((1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode)));

    // Convert to the actual term to be used in e_out = e_in + dt*energy
    a = 1.e0 / (1.e0 + amrex::Math::abs(z_vode));

#ifdef AMREX_USE_HIP
    if((f_rhs_data->e_src_vode!=NULL)||true)
#else
    if((f_rhs_data->e_src_vode!=NULL))
#endif
      energy  = energy / rho_vode / a + f_rhs_data->e_src_vode[idx];
    else
        energy = energy / rho_vode / a;
    f_rhs_data->T_vode[idx]=T_vode;
    f_rhs_data->ne_vode[idx]=ne_vode;
    f_rhs_data->rho_vode[idx]=rho_vode;
    return 0;
}

void sdc_writeOn(amrex::MultiFab &S_old,
                 amrex::MultiFab &S_new,
                 amrex::MultiFab &D_old,
                 amrex::MultiFab &hydro_src,
                 amrex::MultiFab &IR,
                 amrex::MultiFab &reset_src,
                 MFItInfo tiling,
                 const Real& a, const amrex::Real& a_end, const Real& delta_time,
                 long int& old_max_steps, long int& new_max_steps,
                 const int sdc_iter, const int loc_nStep,
                 std::string filename_inputs, std::string filename, std::string filename_chunk_prefix)
{

    //    if(ParallelDescriptor::IOProcessor())
        {
            {
            std::ofstream ofs_inputs(filename_inputs.c_str());
            ParmParse::dumpTable(ofs_inputs, true);
            //      ofs_inputs << "nyx.initial_a = "<<a<<std::endl;
            ofs_inputs << "nyx.initial_z = "<<1/a-1<<std::endl;
            //              ofs_inputs << "nyx.final_a = "<<a_end<<std::endl;
            ofs_inputs << "nyx.final_z = "<<1/a_end-1<<std::endl;
            ofs_inputs << "nyx.fixed_dt = "<<delta_time<<std::endl;
            ofs_inputs << "nyx.hctest_filename_inputs = "<<filename_inputs<<std::endl;
            ofs_inputs << "nyx.hctest_filename_badmap = "<<filename<<std::endl;
            ofs_inputs << "nyx.hctest_filename_chunk = "<<filename_chunk_prefix<<std::endl;
            ofs_inputs << "nyx.hctest_endIndex = "<<(MFIter(S_old,tiling).length())<<std::endl;
            ofs_inputs << "nyx.hctest_example_write = 0"<<std::endl;
            ofs_inputs << "nyx.hctest_example_read = 1"<<std::endl;
            ofs_inputs << "nyx.do_dm_particles = 0"<<std::endl;
            ofs_inputs << "nyx.do_hydro = 0"<<std::endl;
            ofs_inputs << "nyx.hctest_example_index = "<<loc_nStep<<std::endl;
            }
            {
                std::ofstream ofs(filename.c_str());
                (S_old.boxArray()).writeOn(ofs);
                (S_old.DistributionMap()).writeOn(ofs);
            }
        }

    // Since this is writing as a fab, set fab.format=ASCII if you want ascii output
    // The order is: old state data, old derived data, new state data, source terms from hydro, a correction term from enforcing conservation, the lagged IR term from old reactions
#ifdef AMREX_USE_OMP
#ifndef AMREX_USE_GPU
#pragma omp parallel
#endif
#endif
    for ( MFIter mfi(S_old, tiling); mfi.isValid(); ++mfi)
    {
      std::string filename_chunk = filename_chunk_prefix + std::to_string(mfi.index());
      std::ofstream ofs(filename_chunk.c_str());

      S_old[mfi].writeOn(ofs);
      D_old[mfi].writeOn(ofs);
      S_new[mfi].writeOn(ofs);
      hydro_src[mfi].writeOn(ofs);
      reset_src[mfi].writeOn(ofs);
      IR[mfi].writeOn(ofs);

      /*
      S_old[mfi].setVal(0);                                                                           
      D_old[mfi].setVal(0);                                                                           
      S_new[mfi].setVal(0);                                                                           
      hydro_src[mfi].setVal(0);                                                                       
      reset_src[mfi].setVal(0);                                                                       
      IR[mfi].setVal(0);
      */
    }
}

void sdc_readFrom(amrex::MultiFab &S_old,
                  amrex::MultiFab &S_new,
                  amrex::MultiFab &D_old,
                  amrex::MultiFab &hydro_src,
                  amrex::MultiFab &IR,
                  amrex::MultiFab &reset_src,
                  MFItInfo tiling,
                  const Real& a, const amrex::Real& a_end, const Real& delta_time,
                  long int& old_max_steps, long int& new_max_steps,
                  const int sdc_iter, const int loc_nStep,
                  std::string filename_inputs, std::string filename, std::string filename_chunk_prefix)
{

    BoxArray grids;
    DistributionMapping dmap;

    std::ifstream ifs(filename.c_str());
    grids.readFrom(ifs);
    dmap.readFrom(ifs);

    //Sanity check that we're getting the same BoxArray and DistributionMapping
    AMREX_ALWAYS_ASSERT(S_old.boxArray().CellEqual(grids) && dmap == S_old.DistributionMap());

#ifdef AMREX_USE_OMP
#ifndef AMREX_USE_GPU
#pragma omp parallel
#endif
#endif
    for ( MFIter mfi(S_old, tiling); mfi.isValid(); ++mfi)
    {
      std::string filename_chunk = filename_chunk_prefix + std::to_string(mfi.index());
      std::ifstream ifs(filename_chunk.c_str());
    
      S_old[mfi].readFrom(ifs);
      D_old[mfi].readFrom(ifs);
      S_new[mfi].readFrom(ifs);
      hydro_src[mfi].readFrom(ifs);
      reset_src[mfi].readFrom(ifs);
      IR[mfi].readFrom(ifs);
    }
}
        
#endif
