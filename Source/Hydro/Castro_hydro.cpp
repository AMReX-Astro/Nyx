#include "Nyx.H"
#include "Nyx_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Nyx::cons_to_prim(MultiFab& Sborder, MultiFab& q, MultiFab& qaux, MultiFab& grav, MultiFab& sources_for_hydro, MultiFab& src_q, MultiFab& csml, Real a_old, Real a_new, Real dt)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);
	//	const Box& qbx = mfi.tilebox();

	FArrayBox* fab_Sborder = Sborder.fabPtr(mfi);
	FArrayBox* fab_sources_for_hydro = sources_for_hydro.fabPtr(mfi);
	FArrayBox* fab_grav = grav.fabPtr(mfi);

	FArrayBox* fab_q = q.fabPtr(mfi);
	FArrayBox* fab_qaux = qaux.fabPtr(mfi);
	FArrayBox* fab_src_q = src_q.fabPtr(mfi);
	FArrayBox* fab_csml = csml.fabPtr(mfi);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.
  
#pragma gpu
	amrex::Cuda::setLaunchRegion(true);
	AMREX_LAUNCH_DEVICE_LAMBDA(qbx, tqbx,
	{
        ca_ctoprim(AMREX_INT_ANYD(tqbx.loVect()), AMREX_INT_ANYD(tqbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(*fab_Sborder),
                   BL_TO_FORTRAN_ANYD(*fab_q),
                   BL_TO_FORTRAN_ANYD(*fab_qaux),
		   BL_TO_FORTRAN_ANYD(*fab_csml));
	
        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
#pragma gpu
            ca_srctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                         BL_TO_FORTRAN_ANYD(*fab_q),
                         BL_TO_FORTRAN_ANYD(*fab_qaux),
                         BL_TO_FORTRAN_ANYD(*fab_grav),
			 BL_TO_FORTRAN_ANYD(*fab_sources_for_hydro),
                         BL_TO_FORTRAN_ANYD(*fab_src_q),
			 &a_old, &a_new, &dt);
	});
	amrex::Cuda::setLaunchRegion(false);


    }

}
