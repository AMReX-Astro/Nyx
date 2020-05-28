#include "Nyx.H"
#include "Nyx_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Nyx::cons_to_prim(MultiFab& Sborder, MultiFab& q, MultiFab& qaux, MultiFab& grav, MultiFab& sources_for_hydro, MultiFab& src_q, Real a_old, Real a_new, Real dt)
{

  BL_PROFILE("Nyx::cons_to_prim()");
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);
	//	const Box& qbx = mfi.tilebox();

	const auto fab_Sborder = Sborder.array(mfi);
	const auto fab_sources_for_hydro = sources_for_hydro.array(mfi);
	const auto fab_grav = grav.array(mfi);

	const auto fab_q = q.array(mfi);
	const auto fab_qaux = qaux.array(mfi);
	const auto fab_src_q = src_q.array(mfi);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.
	Sborder[mfi].prefetchToDevice();
	q[mfi].prefetchToDevice();
	qaux[mfi].prefetchToDevice();
	AMREX_LAUNCH_DEVICE_LAMBDA(qbx, tqbx,
	{
        ca_ctoprim(AMREX_INT_ANYD(tqbx.loVect()), AMREX_INT_ANYD(tqbx.hiVect()),
                   BL_ARR4_TO_FORTRAN_3D(fab_Sborder),
                   BL_ARR4_TO_FORTRAN_3D(fab_q),
                   BL_ARR4_TO_FORTRAN_3D(fab_qaux));
	});
        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
	q[mfi].prefetchToDevice();
	qaux[mfi].prefetchToDevice();
	grav[mfi].prefetchToDevice();
	sources_for_hydro[mfi].prefetchToDevice();
	src_q[mfi].prefetchToDevice();

	AMREX_LAUNCH_DEVICE_LAMBDA(qbx, tqbx,
	{
	ca_srctoprim(AMREX_INT_ANYD(tqbx.loVect()), AMREX_INT_ANYD(tqbx.hiVect()),
		     BL_ARR4_TO_FORTRAN_3D(fab_q),
		     BL_ARR4_TO_FORTRAN_3D(fab_qaux),
		     BL_ARR4_TO_FORTRAN_3D(fab_grav),
		     BL_ARR4_TO_FORTRAN_3D(fab_sources_for_hydro),
		     BL_ARR4_TO_FORTRAN_3D(fab_src_q),
		     a_old, a_new, dt);
	});

    }

}
