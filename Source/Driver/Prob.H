#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags);

void prob_param_special_fill(amrex::GpuArray<amrex::Real,max_prob_param>& /*prob_param*/);

#ifndef NO_HYDRO
void prob_initdata(const int i,
                   const int j,
                   const int k,
                   amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<amrex::Real> const& diag_eos,
                   amrex::GeometryData const& geomdata,
                   const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param);

void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param);

void prob_initdata_on_box(const Box& bx,
                          Array4<amrex::Real> const& state,
                          Array4<amrex::Real> const& diag_eos,
                          GeometryData const& geomdata,
                          const GpuArray<Real,max_prob_param>& prob_param);

void prob_initdata_state_on_box(const Box& bx,
                                Array4<amrex::Real> const& state,
                                GeometryData const& geomdata,
                                const GpuArray<Real,max_prob_param>& prob_param);
#endif

enum Prob_Type_Index {
  h_species_comp = 0,
  he_species_comp = 1,
  z_in_comp = 2,
  gamma_comp = 3,
  comoving_type_comp = 4,
  prob_type_comp = 5 // Note that if we add more components here, we must adjust the indices
                     // in HydroTests/Prob.H
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
    // Set to default values from Nyx.H
    amrex::Real  initial_z_loc = Nyx::initial_z;
#ifndef NO_HYDRO
    amrex::Real      gamma_loc = Nyx::gamma;
    amrex::Real  h_species_loc = Nyx::h_species;
    amrex::Real he_species_loc = Nyx::he_species;
#endif

    // Potentially override from values in inputs
    amrex::ParmParse pp_nyx("nyx");
    pp_nyx.query("initial_z" , initial_z_loc);
#ifndef NO_HYDRO
    pp_nyx.query("gamma"     , gamma_loc);
    pp_nyx.query("h_species" , h_species_loc);
    pp_nyx.query("he_species", he_species_loc);
#endif

#ifndef NO_HYDRO
    // Store in prob_param
    prob_param[h_species_comp ] =  h_species_loc;
    prob_param[he_species_comp] = he_species_loc;
    prob_param[z_in_comp]       =  initial_z_loc;
    prob_param[gamma_comp ]     =      gamma_loc;
#endif

    amrex::ParmParse pp("prob");
    prob_param[prob_type_comp] = 0;
    prob_param[comoving_type_comp] = 1.0;
    //    pp.query("comoving_type",prob_param[comoving_type_comp]);
    pp.query("prob_type",prob_param[prob_type_comp]);
}

#endif
