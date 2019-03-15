#include "Nyx.H"
#include "Nyx_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Nyx::cons_to_prim(MultiFab& Sborder, MultiFab& q, MultiFab& qaux, MultiFab& grav, MultiFab& sources_for_hydro, MultiFab& src_q, MultiFab& csml, Real& a_old, Real& a_new, Real& dt)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);
	//	const Box& qbx = mfi.tilebox();

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.
  
#pragma gpu
        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]),
		   BL_TO_FORTRAN_ANYD(csml[mfi]));

        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
#pragma gpu
            ca_srctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                         BL_TO_FORTRAN_ANYD(q[mfi]),
                         BL_TO_FORTRAN_ANYD(qaux[mfi]),
                         BL_TO_FORTRAN_ANYD(grav[mfi]),
			 BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                         BL_TO_FORTRAN_ANYD(src_q[mfi]),
			 &a_old, &a_new, &dt);

	    /*	    ca_srctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                         BL_TO_FORTRAN_ANYD(q[mfi]),
                         BL_TO_FORTRAN_ANYD(qaux[mfi]),
                         BL_TO_FORTRAN_ANYD(grav[mfi]),
			 BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                         BL_TO_FORTRAN_ANYD(src_q[mfi]),
			 &a_old, &a_new, &dt);*/
			 

    }

}

// Convert a MultiFab with conservative state data u to a primitive MultiFab q.
void
Nyx::cons_to_prim(MultiFab& u, MultiFab& q, MultiFab& qaux, Real time)
{

    BL_PROFILE("Nyx::cons_to_prim()");

    BL_ASSERT(u.nComp() == NUM_STATE);
    BL_ASSERT(q.nComp() == QVAR);
    BL_ASSERT(u.nGrow() >= q.nGrow());

    int ng = q.nGrow();

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(u, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng);

	ca_ctoprim(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		   BL_TO_FORTRAN_ANYD(u[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
		   BL_TO_FORTRAN_ANYD(q[mfi]),
		   BL_TO_FORTRAN_ANYD(qaux[mfi]),
		   BL_TO_FORTRAN_ANYD(u[mfi]));

    }

}
/*
void
Nyx::check_for_cfl_violation(const Real dt)
{

    Real courno = -1.0e+200;

    const Real *dx = geom.CellSize();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_compute_cfl(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       dt, AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MAX(&courno), print_fortran_warnings);

    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
      amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << "no cfl/_violation flag stored" << std::endl << std::endl;

	//        cfl_violation = 1;
    }

    }*/
