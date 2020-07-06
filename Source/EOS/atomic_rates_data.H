#ifndef _Atomic_Rates_H_
#define _Atomic_Rates_H_


#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

static constexpr int  NCOOLFILE = 301;    // TREECOOL_middle
static constexpr int  NCOOLTAB  = 2000 + 1;

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> lzr;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rggh0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> reh0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehep;

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> AlphaHepp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Alphad;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeH0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeHe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> GammaeHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaH0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaHe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> BetaHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Betaff1;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> Betaff4;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB> RecHepp;


#endif
