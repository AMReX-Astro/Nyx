#include <atomic_rates_data.H>

AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> lzr;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rggh0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> reh0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehep;

AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHepp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Alphad;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeH0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeHe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaH0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaHe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Betaff1;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Betaff4;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHepp;
