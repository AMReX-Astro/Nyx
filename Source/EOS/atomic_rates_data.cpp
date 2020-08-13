#include <atomic_rates_data.H>

AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> lzr;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rggh0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> reh0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehep;

AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHepp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Alphad;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeH0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeHe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaH0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaHe0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Betaff1;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Betaff4;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHp;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHep;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHepp;
