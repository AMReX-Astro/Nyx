#ifndef _Atomic_Rates_H_
#define _Atomic_Rates_H_


#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

static constexpr int  NCOOLFILE = 301;    // TREECOOL_middle
static constexpr int  NCOOLTAB  = 2000 + 1;
static constexpr amrex::Real TCOOLMAX  = 9.0;
static constexpr amrex::Real TCOOLMIN  = 0.0;
static constexpr amrex::Real xacc  = 1e-6;
//These are ~1e-5 rel diff from m_proton and k_B in constants_cosmo.H
/*
static constexpr amrex::Real MPROTON = 1.6726231e-24;
static constexpr amrex::Real BOLTZMANN = 1.3806e-16;*/
//Fortran added noise version:
static constexpr amrex::Real MPROTON = 1.6726230999999999E-024;
static constexpr amrex::Real BOLTZMANN = 1.3806000442045675E-016;

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> lzr;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rggh0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rgghep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> reh0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLFILE> rehep;

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> AlphaHepp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Alphad;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeH0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeHe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> GammaeHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaH0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaHe0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> BetaHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Betaff1;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> Betaff4;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHp;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHep;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,NCOOLTAB+1> RecHepp;


#endif
