#include "Nyx.H"
#include "Prob.H"

void prob_param_special_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
    prob_param[comoving_type_comp] = 0.0;
}

