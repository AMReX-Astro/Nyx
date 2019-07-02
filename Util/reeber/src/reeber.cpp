#include "Nyx.H"
#include "reeber.H"

using namespace amrex;

/* // Global components of state
int Nyx::Density = -1;
int Nyx::Eden = -1;
int Nyx::Eint = -1;
int Nyx::Xmom = -1;
int Nyx::Ymom = -1;
int Nyx::Zmom = -1;*/

/* // Global density values
   average_gas_density;
   average_dm_density;
   average_neutr_density;
   average_total_density;*/

void Nyx::runReeberAnalysis(const amrex::MultiFab & new_state, amrex::MultiFab & particle_mf, const amrex::Geometry Geom, int nStep, bool do_analysis, std::vector<Halo> reeber_halos)
{

  if(verbose)
    amrex::Print()<<"Running Reeber anaylsis"<<std::endl;
  for ( MFIter mfi(new_state, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      const Box& tbx = mfi.tilebox();
      Array4<const Real> state = new_state.array(mfi);
      Array4<Real>  particle = particle_mf.array(mfi);
      const Dim3 lo = amrex::lbound(tbx);
      const Dim3 hi = amrex::ubound(tbx);
      int ncomp = 4;

      for (int n = 0; n < ncomp; ++n) {
	for (int z = lo.z; z <= hi.z; ++z) {
	  for (int y = lo.y; y <= hi.y; ++y) {
	            AMREX_PRAGMA_SIMD
		      for (int x = lo.x; x <= hi.x; ++x) {
			state(x,y,z,n);
			particle(x,y,z,n);
		      }
	  }
	}
      }


    }
  return;
}
