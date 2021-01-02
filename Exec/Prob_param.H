#ifndef _Prob_param_H_
#define _Prob_param_H_

#include "Nyx.H"

enum Prob_Type_Index {
  h_species_comp = 0,
  he_species_comp = 1,
  z_in_comp = 2,
  gamma_comp = 3, 
  prob_type_comp = 4 // Note that if we add more components here, we must adjust the indices
                     // in HydroTests/Prob.H
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
    // Set to default values from Nyx.H
    Real  initial_z_loc = Nyx::initial_z;
    Real      gamma_loc = Nyx::gamma;
    Real  h_species_loc = Nyx::h_species;
    Real he_species_loc = Nyx::he_species;

    // Potentially override from values in inputs
    amrex::ParmParse pp_nyx("nyx");
    pp_nyx.query("initial_z" , initial_z_loc);
    pp_nyx.query("gamma"     , gamma_loc);
    pp_nyx.query("h_species" , h_species_loc);
    pp_nyx.query("he_species", he_species_loc);

    // Store in prob_param
    prob_param[h_species_comp ] =  h_species_loc;
    prob_param[he_species_comp] = he_species_loc;
    prob_param[z_in_comp]       =  initial_z_loc;
    prob_param[gamma_comp ]     =      gamma_loc;

    amrex::ParmParse pp("prob");
    prob_param[prob_type_comp] = 0;
    pp.query("prob_type",prob_param[prob_type_comp]);
}

#endif
