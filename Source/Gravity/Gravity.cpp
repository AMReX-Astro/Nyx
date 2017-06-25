#include <cmath>

#include <AMReX_ParmParse.H>
#include "Gravity.H"
#include "Nyx.H"
#include <Gravity_F.H>
#include <Nyx_F.H>

#include <AMReX_MultiGrid.H>
#include <AMReX_Laplacian.H>
#include <AMReX_MacBndry.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>
#include <mg_cpp_f.h>

#ifdef USEHPGMG
#include <BL_HPGMG.H>
#endif

using namespace amrex;

// MAX_LEV defines the maximum number of AMR levels allowed by the parent "Amr" object
#define MAX_LEV 15

// Give this a bogus default value to force user to define in inputs file
std::string Gravity::gravity_type = "fill_me_please";
int  Gravity::verbose       = 0;
int  Gravity::show_timings  = 0;
int  Gravity::no_sync       = 0;
int  Gravity::no_composite  = 0;
int  Gravity::dirichlet_bcs = 0;
int  Gravity::monopole_bcs  = 0;
int  Gravity::solve_with_cpp= 0;
int  Gravity::solve_with_hpgmg = 0;
Real Gravity::sl_tol        = 1.e-12;
Real Gravity::ml_tol        = 1.e-12;
Real Gravity::delta_tol     = 1.e-12;
Real Gravity::mass_offset   = 0;
int  Gravity::stencil_type  = CC_CROSS_STENCIL;

extern "C"
{void fort_get_grav_const(Real* Gconst);}

// ************************************************************************** //

// Ggravity is defined as -4 * pi * G, where G is the gravitational constant.
// G is defined as Gconst in `fParallel/extern/constants/nyx_constants.f90` if
// NYX is defined in the GNUmakefile. G is defined as Gconst in
// `fParallel/extern/constants/constants.f90` if NYX is not defined in the
// GNUmakefile

// In CGS, this constant is currently
//      Gconst   =  6.67428e-8           cm^3/g/s^2 , which results in
//      Ggravity = -83.8503442814844e-8  cm^3/g/s^2

// In cosmological units, this constant is currently
//      Gconst   =   4.3   e-9  (km/s)^2 Mpc / Msun, which results in
//      Ggravity = -54.0354e-9  (km/s)^2 Mpc / Msun

// ************************************************************************** //

static Real Ggravity = 0;

Gravity::Gravity (Amr*   Parent,
                  int    _finest_level,
                  BCRec* _phys_bc,
                  int    _density)
  :
    parent(Parent),
    LevelData(MAX_LEV),
    grad_phi_curr(MAX_LEV),
    grad_phi_prev(MAX_LEV),
    phi_flux_reg(MAX_LEV),
    grids(Parent->boxArray()),
    dmap(Parent->DistributionMap()),
    level_solver_resnorm(MAX_LEV),
    phys_bc(_phys_bc)
{
     density = _density;
     read_params();
     finest_level_allocated = -1;
#ifdef CGRAV
     if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav" || gravity_type == "StaticGrav")
          make_mg_bc();
#else
     if(gravity_type == "PoissonGrav") make_mg_bc();
#endif
}


Gravity::~Gravity ()
{
    // nothing to see here.
}

void
Gravity::read_params ()
{
    static bool done = false;

    if (!done)
    {
        const Real strt = ParallelDescriptor::second();

        ParmParse pp("gravity");
        pp.get("gravity_type", gravity_type);

#ifdef CGRAV
        if (gravity_type != "PoissonGrav" && gravity_type != "CompositeGrav" && gravity_type != "StaticGrav")
        {
            std::cout << "Sorry -- dont know this gravity type" << std::endl;
            amrex::Abort("Options are PoissonGrav, CompositeGrav and StaticGrav");
        }
#else
        if (gravity_type != "PoissonGrav")
        {
            std::cout << "Sorry -- dont know this gravity type" << std::endl;
            amrex::Abort("Options are PoissonGrav");
        }
#endif

        pp.query("v", verbose);
        pp.query("show_timings", show_timings);
        pp.query("no_sync", no_sync);
        pp.query("no_composite", no_composite);

        pp.query("dirichlet_bcs", dirichlet_bcs);
        pp.query("monopole_bcs"  , monopole_bcs);

        pp.query("solve_with_cpp", solve_with_cpp);
        pp.query("solve_with_hpgmg", solve_with_hpgmg);

        if (solve_with_cpp && solve_with_hpgmg)
          amrex::Error("Multiple gravity solvers selected.");

#ifndef USEHPGMG
        if (solve_with_hpgmg)
          amrex::Error("To use the HPGMG solver you must compile with USE_HPGMG = TRUE");
#endif

        // Allow run-time input of solver tolerances
        pp.query("ml_tol", ml_tol);
        pp.query("sl_tol", sl_tol);
        pp.query("delta_tol", delta_tol);

        Real Gconst;
        fort_get_grav_const(&Gconst);
        Ggravity = -4.0 * M_PI * Gconst;
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Getting Gconst from nyx_constants: " << Gconst
                      << '\n';
            std::cout << "Using " << Ggravity << " for 4 pi G in Gravity.cpp "
                      << '\n';
        }

        done = true;

        if (show_timings)
        {
            const int IOProc = ParallelDescriptor::IOProcessorNumber();
            Real end = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
            Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealMax(end,IOProc);
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Gravity::read_params() time = " << end << '\n';
#ifdef BL_LAZY
        });
#endif
        }
    }
}

void
Gravity::install_level (int       level,
                        AmrLevel* level_data_to_install)
{
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << "Installing Gravity level " << level << '\n';

    LevelData[level] = level_data_to_install;

    level_solver_resnorm[level] = 0;

#ifdef CGRAV
    if (gravity_type != "StaticGrav")
    {
#endif

    const auto& dm = level_data_to_install->DistributionMap();

    grad_phi_prev[level].resize(BL_SPACEDIM);
    for (int n=0; n<BL_SPACEDIM; ++n)
    {
        grad_phi_prev[level][n].reset(new MultiFab(level_data_to_install->getEdgeBoxArray(n),dm,1,1));
        grad_phi_prev[level][n]->setVal(0.);
    }

    grad_phi_curr[level].resize(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        grad_phi_curr[level][n].reset(new MultiFab(level_data_to_install->getEdgeBoxArray(n),dm,1,1));
        grad_phi_curr[level][n]->setVal(0.);
    }

    if (level > 0)
    {
        IntVect crse_ratio = parent->refRatio(level-1);
        phi_flux_reg[level].reset(new FluxRegister(level_data_to_install->boxArray(),
						   dm, crse_ratio, level, 1));
    }

#ifdef CGRAV
    }
#endif

    finest_level_allocated = level;
}

std::string
Gravity::get_gravity_type ()
{
    return gravity_type;
}

int
Gravity::get_no_sync ()
{
    return no_sync;
}

int
Gravity::get_no_composite ()
{
    return no_composite;
}

Array<MultiFab*>
Gravity::get_grad_phi_prev (int level)
{
    return amrex::GetArrOfPtrs(grad_phi_prev[level]);
}

Array<MultiFab*>
Gravity::get_grad_phi_curr (int level)
{
    return amrex::GetArrOfPtrs(grad_phi_curr[level]);
}

void
Gravity::plus_grad_phi_curr (int level, const Array<MultiFab*>& addend)
{
    for (int n = 0; n < BL_SPACEDIM; n++)
        grad_phi_curr[level][n]->plus(*addend[n], 0, 1, 0);
}

void
Gravity::swap_time_levels (int level)
{

#ifdef CGRAV
    if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav")
#else
    if (gravity_type == "PoissonGrav")
#endif
    {
        for (int n=0; n < BL_SPACEDIM; n++)
        {
	    std::swap(grad_phi_prev[level][n], grad_phi_curr[level][n]);
            grad_phi_curr[level][n].reset(new MultiFab(BoxArray(grids[level]).surroundingNodes(n), 
						       dmap[level], 1, 1));
            grad_phi_curr[level][n]->setVal(1.e50);
        }
    }
}

void
Gravity::zero_phi_flux_reg (int level)
{
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    phi_flux_reg[level]->setVal(0);
}

void
Gravity::solve_for_old_phi (int               level,
                            MultiFab&         phi,
                            const Array<MultiFab*>& grad_phi,
                            int               fill_interior)
{
    BL_PROFILE("Gravity::solve_for_old_phi()");
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Gravity ... single level solve for old phi at level "
                  << level << std::endl;
    MultiFab Rhs(grids[level], dmap[level], 1, 0);
    Rhs.setVal(0.0);

#ifndef NO_HYDRO
    if (Nyx::Do_Hydro() == 1)
    {
       MultiFab&  S_old = LevelData[level]->get_old_data(State_Type);
       MultiFab::Copy(Rhs, S_old, density, 0, 1, 0);
    }
#endif

    AddParticlesToRhs(level,Rhs,1);

    // We shouldn't need to use virtual or ghost particles for old phi solves.

    const Real time  = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Gravity::solve_for_new_phi (int               level,
                            MultiFab&         phi,
                            const Array<MultiFab*>& grad_phi,
                            int               fill_interior,
                            int               grav_n_grow)
{
    BL_PROFILE("Gravity::solve_for_new_phi()");
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Gravity ... single level solve for new phi at level "
                  << level << std::endl;

    MultiFab Rhs(grids[level], dmap[level], 1, 0);
    Rhs.setVal(0.0);

#ifndef NO_HYDRO
    if (Nyx::Do_Hydro() == 1)
    {
       MultiFab& S_new = LevelData[level]->get_new_data(State_Type);
       MultiFab::Copy(Rhs, S_new, density, 0, 1, 0);
    }
#endif

    AddParticlesToRhs(level,Rhs,grav_n_grow);
    AddVirtualParticlesToRhs(level,Rhs,grav_n_grow);
    AddGhostParticlesToRhs(level,Rhs);

    const Real time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Gravity::solve_for_phi (int               level,
                        MultiFab&         Rhs,
                        MultiFab&         phi,
                        const Array<MultiFab*>& grad_phi,
                        Real              time,
                        int               fill_interior)

{
    BL_PROFILE("Gravity::solve_for_phi()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << " ... solve for phi at level " << level << '\n';

    const Real strt = ParallelDescriptor::second();

    // This is a correction for fully periodic domains only
    if (Geometry::isAllPeriodic())
        CorrectRhsUsingOffset(level,Rhs);

    Rhs.mult(Ggravity);

    Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(level));

    BL_ASSERT(cs != 0);

    // Here we divide by a for the Poisson solve.
    Rhs.mult(1 / cs->get_comoving_a(time));

    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level], dmap[level], 1, geom);

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                     : IntVect::TheZeroVector();
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

#ifndef NDEBUG
    if (Rhs.contains_nan(0,1,0))
    {
        std::cout << "Rhs in solve_for_phi at level " << level << " has NaNs" << std::endl;
        amrex::Abort("");
    }
#endif

    // Need to set the boundary values here so they can get copied into "bndry"
    if (dirichlet_bcs) set_dirichlet_bcs(level,&phi);

    if (level == 0)
    {
        bndry.setBndryValues(phi, src_comp, dest_comp, num_comp, *phys_bc);
    }
    else
    {
        MultiFab c_phi;
        get_crse_phi(level, c_phi, time);
        BoxArray crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes, dmap[level],
			      in_rad, out_rad, extent_rad, num_comp);
        crse_br.copyFrom(c_phi, c_phi.nGrow(), src_comp, dest_comp, num_comp);
        bndry.setBndryValues(crse_br, src_comp, phi, src_comp, dest_comp,
                             num_comp, crse_ratio, *phys_bc);
    }

    Array<BoxArray> bav(1);
    bav[0] = phi.boxArray();
    Array<DistributionMapping> dmv(1);
    dmv[0] = Rhs.DistributionMap();
    Array<Geometry> fgeom(1);
    fgeom[0] = geom;

    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0;
            xb[0][i] = 0;
        }
    }
    else
    {
        const Real* dx_crse = parent->Geom(level-1).CellSize();
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0.5 * dx_crse[i];
            xb[0][i] = 0.5 * dx_crse[i];
        }
    }

    if ( Geometry::isAllPeriodic() )
    {
        if (grids[level].contains(parent->Geom(level).Domain()))
        {
            Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(level)));

            Real sum = nyx_level->vol_weight_sum(Rhs,false);

            Rhs.plus(-sum, 0, 1, 0);

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << " ... subtracting "
                          << sum
                          << " to ensure solvability "
                          << '\n';
        }
    }

    Array<MultiFab*> phi_p = { &phi };
    Array<MultiFab*> Rhs_p = { &Rhs };

    const Real  tol     = sl_tol;
    const Real  abs_tol = 0.;

    if (solve_with_cpp)
    {
        solve_with_Cpp(level, phi, grad_phi, Rhs, tol, abs_tol);
    }
    else if (solve_with_hpgmg)
    {
#ifdef USEHPGMG
        solve_with_HPGMG(level, phi, grad_phi, Rhs, tol, abs_tol);
#endif
    }
    else
    {
        MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);
        mgt_solver.set_const_gravity_coeffs(xa, xb);
        const int   mglev   = 0;
        const Real* dx      = geom.CellSize();

	int always_use_bnorm = 0;
	int need_grad_phi = 1;
        mgt_solver.solve(phi_p, Rhs_p, bndry, tol, abs_tol, always_use_bnorm,
			 level_solver_resnorm[level], need_grad_phi);

        mgt_solver.get_fluxes(mglev, grad_phi, dx);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::solve_for_phi() time = " << end << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

void
Gravity::solve_for_delta_phi (int                        crse_level,
                              int                        fine_level,
                              MultiFab&                  crse_rhs,
                              const Array<MultiFab*>&         delta_phi,
                              const Array<Array<MultiFab*> >& grad_delta_phi)
{
    BL_PROFILE("Gravity::solve_for_delta_phi()");
    const int num_levels = fine_level - crse_level + 1;
    const Box& crse_domain = (parent->Geom(crse_level)).Domain();

    BL_ASSERT(grad_delta_phi.size() == num_levels);
    BL_ASSERT(delta_phi.size() == num_levels);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "... solving for delta_phi at crse_level = " << crse_level
                  << '\n';
        std::cout << "...                    up to fine_level = " << fine_level
                  << '\n';
    }

    const Geometry& geom = parent->Geom(crse_level);
    MacBndry bndry(grids[crse_level], dmap[crse_level], 1, geom);

    IntVect crse_ratio = crse_level > 0 ? parent->refRatio(crse_level-1)
                                          : IntVect::TheZeroVector();

    // Set homogeneous Dirichlet values for the solve.
    bndry.setHomogValues(*phys_bc, crse_ratio);

    Array<BoxArray>            bav(num_levels);
    Array<DistributionMapping> dmv(num_levels);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        bav[lev-crse_level] = grids[lev];
        MultiFab& phi_new = LevelData[lev]->get_new_data(PhiGrav_Type);
        dmv[lev-crse_level] = phi_new.DistributionMap();
    }
    Array<Geometry> fgeom(num_levels);
    for (int lev = crse_level; lev <= fine_level; lev++)
        fgeom[lev-crse_level] = parent->Geom(lev);

    MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);

    Array< Array<Real> > xa(num_levels);
    Array< Array<Real> > xb(num_levels);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        xa[lev-crse_level].resize(BL_SPACEDIM);
        xb[lev-crse_level].resize(BL_SPACEDIM);
        if (lev == 0)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev-crse_level][i] = 0;
                xb[lev-crse_level][i] = 0;
            }
        }
        else
        {
            const Real* dx_crse = parent->Geom(lev-1).CellSize();
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev-crse_level][i] = 0.5 * dx_crse[i];
                xb[lev-crse_level][i] = 0.5 * dx_crse[i];
            }
        }
    }

    Array<std::unique_ptr<MultiFab> > raii;
    Array<MultiFab*> Rhs_p(num_levels);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        delta_phi[lev-crse_level]->setVal(0);

        if (lev == crse_level)
        {
            Rhs_p[0] = &crse_rhs;
        }
        else
        {
	    raii.push_back(std::unique_ptr<MultiFab>(new MultiFab(grids[lev], dmap[lev], 1, 0)));
	    Rhs_p[lev-crse_level] = raii.back().get();
            Rhs_p[lev-crse_level]->setVal(0);
        }

    }

    // If at coarsest level, subtract off average of RHS from all levels to ensure solvability
    if (Geometry::isAllPeriodic() &&
           (grids[crse_level].numPts() == crse_domain.numPts())) {
       Real local_correction = 0.0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:local_correction)
#endif
       for (MFIter mfi(crse_rhs,true); mfi.isValid(); ++mfi) {
           local_correction += crse_rhs[mfi].sum(mfi.tilebox(),0,1);
       }
       ParallelDescriptor::ReduceRealSum(local_correction);

       local_correction = local_correction / grids[crse_level].numPts();

       if (verbose && ParallelDescriptor::IOProcessor())
          std::cout << "WARNING: Adjusting RHS in solve_for_delta_phi by " << local_correction << std::endl;

       for (int lev = crse_level; lev <= fine_level; lev++) {
         Rhs_p[lev-crse_level]->plus(-local_correction,0,1,0);
       }
    }

    mgt_solver.set_const_gravity_coeffs(xa, xb);

    const Real tol     = delta_tol;
    Real       abs_tol = level_solver_resnorm[crse_level];
    for (int lev = crse_level + 1; lev < fine_level; lev++)
        abs_tol = std::max(abs_tol,level_solver_resnorm[lev]);

    Real final_resnorm;
    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    mgt_solver.solve(delta_phi, Rhs_p, bndry, tol, abs_tol, always_use_bnorm, final_resnorm, need_grad_phi);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        auto& gdphi = grad_delta_phi[lev-crse_level];
        const Real* dx = parent->Geom(lev).CellSize();
        mgt_solver.get_fluxes(lev-crse_level, gdphi, dx);
    }
}

void
Gravity::gravity_sync (int crse_level, int fine_level, int iteration, int ncycle,
                       const MultiFab& drho_and_drhoU, const MultiFab& dphi,
                       const Array<MultiFab*>& grad_delta_phi_cc)
{
    BL_PROFILE("Gravity::gravity_sync()");
    BL_ASSERT(parent->finestLevel()>crse_level);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << " ... gravity_sync at crse_level " << crse_level << '\n';
        std::cout << " ...     up to finest_level     " << fine_level << '\n';
    }

    // Build Rhs for solve for delta_phi
    MultiFab crse_rhs(grids[crse_level], dmap[crse_level], 1, 0);
    MultiFab::Copy(crse_rhs, drho_and_drhoU, 0, 0, 1, 0);
    crse_rhs.mult(Ggravity);
    crse_rhs.plus(dphi, 0, 1, 0);

    const Geometry& crse_geom   = parent->Geom(crse_level);
    const Box&      crse_domain = crse_geom.Domain();

    // In the all-periodic case we enforce that CrseRhsSync sums to zero.
    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()))
    {
        Real local_correction = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:local_correction)
#endif
        for (MFIter mfi(crse_rhs,true); mfi.isValid(); ++mfi)
            local_correction += crse_rhs[mfi].sum(mfi.tilebox(), 0, 1);
        ParallelDescriptor::ReduceRealSum(local_correction);

        local_correction /= grids[crse_level].numPts();

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "WARNING: Adjusting RHS in gravity_sync solve by " << local_correction << '\n';

        crse_rhs.plus(-local_correction,0,1,0);
    }

    // delta_phi needs a ghost cell for the solve
    Array<std::unique_ptr<MultiFab> >  delta_phi(fine_level - crse_level + 1);
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        delta_phi[lev-crse_level].reset(new MultiFab(grids[lev], dmap[lev], 1, 1));
        delta_phi[lev-crse_level]->setVal(0);
    }

    Array <Array<std::unique_ptr<MultiFab> > > ec_gdPhi(fine_level - crse_level + 1);
    for (int lev = crse_level; lev <= fine_level; lev++) {
        Nyx* Nyx_lev = dynamic_cast<Nyx*>(&parent->getLevel(lev));
        ec_gdPhi[lev-crse_level].resize(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
           ec_gdPhi[lev-crse_level][n].reset(new MultiFab(Nyx_lev->getEdgeBoxArray(n),
							  Nyx_lev->DistributionMap(),
							  1,0));
    }

    // Do multi-level solve for delta_phi
    solve_for_delta_phi(crse_level, fine_level, crse_rhs, 
			amrex::GetArrOfPtrs(delta_phi),
			amrex::GetArrOfArrOfPtrs(ec_gdPhi));

    crse_rhs.clear();

    // In the all-periodic case we enforce that delta_phi averages to zero.
    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()) ) {
       Real local_correction = 0.0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:local_correction)
#endif
       for (MFIter mfi(*delta_phi[0],true); mfi.isValid(); ++mfi) {
           local_correction += (*delta_phi[0])[mfi].sum(mfi.tilebox(),0,1);
       }
       ParallelDescriptor::ReduceRealSum(local_correction);

       local_correction = local_correction / grids[crse_level].numPts();

       for (int lev = crse_level; lev <= fine_level; lev++) {
           delta_phi[lev-crse_level]->plus(-local_correction,0,1,1);
       }
    }

    // Add delta_phi to phi_new, and grad(delta_phi) to grad(delta_phi_curr) on each level
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        LevelData[lev]->get_new_data(PhiGrav_Type).plus(*delta_phi[lev-crse_level], 0, 1, 0);
        for (int n = 0; n < BL_SPACEDIM; n++)
            grad_phi_curr[lev][n]->plus(*ec_gdPhi[lev-crse_level][n], 0, 1, 0);
    }

    int is_new = 1;

    // Average phi_new from fine to coarse level
    for (int lev = fine_level - 1; lev >= crse_level; lev--)
    {
        const IntVect& ratio = parent->refRatio(lev);
        amrex::average_down(LevelData[lev+1]->get_new_data(PhiGrav_Type),
                             LevelData[lev  ]->get_new_data(PhiGrav_Type),
                             0, 1, ratio);
    }

    // Average the edge-based grad_phi from finer to coarser level
    for (int lev = fine_level-1; lev >= crse_level; lev--)
        average_fine_ec_onto_crse_ec(lev, is_new);

    // Add the contribution of grad(delta_phi) to the flux register below if necessary.
    if (crse_level > 0 && iteration == ncycle)
    {
        for (MFIter mfi(*delta_phi[0]); mfi.isValid(); ++mfi)
            for (int n = 0; n < BL_SPACEDIM; ++n)
                phi_flux_reg[crse_level]->FineAdd((*ec_gdPhi[0][n])[mfi], n,
						  mfi.index(), 0, 0, 1, 1);
    }

    for (int lev = crse_level; lev <= fine_level; lev++) {
        grad_delta_phi_cc[lev-crse_level]->setVal(0.0);
        const Geometry& geom = parent->Geom(lev);
        amrex::average_face_to_cellcenter(*grad_delta_phi_cc[lev-crse_level],
                                           amrex::GetArrOfConstPtrs(ec_gdPhi[lev-crse_level]),
                                           geom);
    }
}

void
Gravity::get_crse_phi (int       level,
                       MultiFab& phi_crse,
                       Real      time)
{
    BL_PROFILE("Gravity::get_crse_phi()");
    BL_ASSERT(level != 0);

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);

    phi_crse.clear();
    phi_crse.define(grids[level-1], dmap[level-1], 1, 1);

    // BUT NOTE we don't trust phi's ghost cells.
    FArrayBox phi_crse_temp;

    // Note that we must do these cases separately because it's possible to do a
    //   new solve after a regrid when the old data on the coarse grid may not yet
    //   be defined.
    for (MFIter mfi(phi_crse,true); mfi.isValid(); ++mfi)
    {
        const Box& gtbx = mfi.growntilebox();

        phi_crse_temp.resize(gtbx,1);

        if (fabs(alpha-1.0) < 1.e-15)
        {
            phi_crse[mfi].copy(LevelData[level-1]->get_new_data(PhiGrav_Type)[mfi]);
        }
        else if (fabs(alpha) < 1.e-15)
        {
            phi_crse[mfi].copy(LevelData[level-1]->get_old_data(PhiGrav_Type)[mfi]);
        }
        else
        {
            phi_crse_temp.copy(LevelData[level-1]->get_old_data(PhiGrav_Type)[mfi]);
            Real omalpha = 1.0 - alpha;
            phi_crse_temp.mult(omalpha);

            phi_crse[mfi].copy(LevelData[level-1]->get_new_data(PhiGrav_Type)[mfi],gtbx);
            phi_crse[mfi].mult(alpha,gtbx);
            phi_crse[mfi].plus(phi_crse_temp);
        }
    }

    const Geometry& geom = parent->Geom(level-1);
    phi_crse.FillBoundary(geom.periodicity());
}

void
Gravity::get_crse_grad_phi (int               level,
                            Array<std::unique_ptr<MultiFab> >& grad_phi_crse,
                            Real              time)
{
    BL_PROFILE("Gravity::get_crse_grad_phi()");
    BL_ASSERT(level!=0);

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);
    const Real omalpha = 1.0 - alpha;

    Nyx* Nyx_crse_lev = dynamic_cast<Nyx*>(&parent->getLevel(level-1));

    BL_ASSERT(grad_phi_crse.size() == BL_SPACEDIM);

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        BL_ASSERT(!grad_phi_crse[i]);
        grad_phi_crse[i].reset(new MultiFab(Nyx_crse_lev->getEdgeBoxArray(i), 
					    Nyx_crse_lev->DistributionMap(),
					    1, 0));

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox grad_phi_crse_temp;
            for (MFIter mfi(*grad_phi_crse[i],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                grad_phi_crse_temp.resize(tbx,1);

                grad_phi_crse_temp.copy((*grad_phi_prev[level-1][i])[mfi]);
                grad_phi_crse_temp.mult(omalpha);

                (*grad_phi_crse[i])[mfi].copy((*grad_phi_curr[level-1][i])[mfi],tbx);
                (*grad_phi_crse[i])[mfi].mult(alpha,tbx);
                (*grad_phi_crse[i])[mfi].plus(grad_phi_crse_temp);
            }
        }
    }
}

void
Gravity::multilevel_solve_for_new_phi (int level,
                                       int finest_level,
                                       int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::multilevel_solve_for_new_phi()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Gravity ... multilevel solve for new phi at base level " << level
                  << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_curr[lev].size()==BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_curr[lev][n].reset(new MultiFab(eba, dmap[lev], 1, 1));
        }
    }

    int is_new = 1;
    actual_multilevel_solve(level, finest_level, 
			    amrex::GetArrOfArrOfPtrs(grad_phi_curr),
                            is_new, use_previous_phi_as_guess);
}

void
Gravity::multilevel_solve_for_old_phi (int level,
                                       int finest_level,
                                       int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::multilevel_solve_for_old_phi()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Gravity ... multilevel solve for old phi at base level " << level
                  << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_prev[lev].size() == BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_prev[lev][n].reset(new MultiFab(eba, dmap[lev], 1, 1));
        }
    }

    int is_new = 0;
    actual_multilevel_solve(level, finest_level,
			    amrex::GetArrOfArrOfPtrs(grad_phi_prev),
                            is_new, use_previous_phi_as_guess);
}

void
Gravity::multilevel_solve_for_phi(int level, int finest_level,
                                  int use_previous_phi_as_guess)
{
    multilevel_solve_for_new_phi(level, finest_level, use_previous_phi_as_guess);
}

void
Gravity::actual_multilevel_solve (int                       level,
                                  int                       finest_level,
                                  const Array<Array<MultiFab*> >& grad_phi,
                                  int                       is_new,
                                  int                       use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::actual_multilevel_solve()");
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const Real      strt = ParallelDescriptor::second();

    const int num_levels = finest_level - level + 1;

    Array<BoxArray> bav(num_levels);
    Array<DistributionMapping> dmv(num_levels);

    // Ok to use phi_new here because phi_new and phi_old have the same DistributionMap
    for (int lev = 0; lev < num_levels; lev++)
    {
        bav[lev]          = grids[level+lev];
        if (is_new == 1)
        {
           MultiFab& phi_new = LevelData[level+lev]->get_new_data(PhiGrav_Type);
           dmv[lev]          = phi_new.DistributionMap();
        } else {
           MultiFab& phi_old = LevelData[level+lev]->get_old_data(PhiGrav_Type);
           dmv[lev]          = phi_old.DistributionMap();
        }
    }
    Array<Geometry> fgeom(num_levels);
    for (int i = 0; i < num_levels; i++)
        fgeom[i] = parent->Geom(level+i);

    // FOR TIMINGS
    if (show_timings)
        ParallelDescriptor::Barrier();

    const Real strt_setup = ParallelDescriptor::second();

    // FOR TIMINGS
    if (show_timings)
        ParallelDescriptor::Barrier();

    if (show_timings)
    {
        Real    end_setup = ParallelDescriptor::second() - strt_setup;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end_setup,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity:: time in mgt_solver setup = " << end_setup << '\n';
#ifdef BL_LAZY
        });
#endif
    }

    Array< Array<Real> > xa(num_levels);
    Array< Array<Real> > xb(num_levels);

    for (int lev = 0; lev < num_levels; lev++)
    {
        xa[lev].resize(BL_SPACEDIM);
        xb[lev].resize(BL_SPACEDIM);
        if (level + lev == 0)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0;
                xb[lev][i] = 0;
            }
        }
        else
        {
            const Real* dx_crse = parent->Geom(level + lev - 1).CellSize();
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0.5 * dx_crse[i];
                xb[lev][i] = 0.5 * dx_crse[i];
            }
        }
    }

    // FOR TIMINGS
    if (show_timings)
        ParallelDescriptor::Barrier();

    const Real strt_rhs = ParallelDescriptor::second();

    Array<MultiFab*> phi_p(num_levels);
    Array<std::unique_ptr<MultiFab> > Rhs_p(num_levels);

    Array<std::unique_ptr<MultiFab> > Rhs_particles(num_levels);
    for (int lev = 0; lev < num_levels; lev++)
    {
	Rhs_particles[lev].reset(new MultiFab(grids[level+lev], dmap[level+lev], 1, 0));
       Rhs_particles[lev]->setVal(0.);
    }

    const auto& rpp = amrex::GetArrOfPtrs(Rhs_particles);
    AddParticlesToRhs(level,finest_level,rpp);
    AddGhostParticlesToRhs(level,rpp);
    AddVirtualParticlesToRhs(finest_level,rpp);

    Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(level));

    BL_ASSERT(cs != 0);

    Real time = 0;

    if (is_new == 1)
    {
        time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    }
    else if (is_new == 0)
    {
        time = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    }

    // Here we get comoving_a b/c the RHS should be 4 * pi * G * density / a
    const Real a_inverse = 1. / (cs->get_comoving_a(time));

// *****************************************************************************

    for (int lev = 0; lev < num_levels; lev++)
    {
        if (is_new == 0)
        {
           // Working in result data structure directly
           phi_p[lev] = &LevelData[level+lev]->get_old_data(PhiGrav_Type);
        }
        else
        {
           // Working in result data structure directly
           phi_p[lev] = &LevelData[level+lev]->get_new_data(PhiGrav_Type);
        }

        if (!use_previous_phi_as_guess)
            phi_p[lev]->setVal(0);

        // Need to set the boundary values before "bndry" is defined so they get copied in
        if (dirichlet_bcs) set_dirichlet_bcs(level+lev,phi_p[lev]);

        Rhs_p[lev].reset(new MultiFab(grids[level+lev], dmap[level+lev], 1, 0));
        Rhs_p[lev]->setVal(0.0);

#ifndef NO_HYDRO
        if (Nyx::Do_Hydro() == 1)
        {
            if (is_new == 1)
            {
                MultiFab::Copy(*Rhs_p[lev], LevelData[level+lev]->get_new_data(State_Type), density, 0, 1, 0);
            }
            else if (is_new == 0)
            {
                MultiFab::Copy(*Rhs_p[lev], LevelData[level+lev]->get_old_data(State_Type), density, 0, 1, 0);
            }
        }
#endif
        MultiFab::Add(*Rhs_p[lev], *Rhs_particles[lev], 0, 0, 1, 0);
    }

    // Average phi from fine to coarse level before the solve.
    for (int lev = num_levels-1; lev > 0; lev--)
    {
        amrex::average_down(*phi_p[lev], *phi_p[lev-1],
                             0, 1, parent->refRatio(level+lev-1));
    }

// *****************************************************************************

    // This correction is for fully periodic domains only.
    if (Geometry::isAllPeriodic())
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << " ... subtracting average density " << mass_offset
                      << " from RHS at each level " << '\n';

        for (int lev = 0; lev < num_levels; lev++)
            for (MFIter mfi(*Rhs_p[lev]); mfi.isValid(); ++mfi)
                (*Rhs_p[lev])[mfi].plus(-mass_offset);

       // This is used to enforce solvability if appropriate.
       if ( grids[level].contains(parent->Geom(level).Domain()) )
       {
           Real sum = 0;
           for (int lev = 0; lev < num_levels; lev++)
           {
               Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(level+lev)));
               sum += nyx_level->vol_weight_sum(*Rhs_p[lev],true);
           }

            sum /= parent->Geom(0).ProbSize();

            const Real eps = 1.e-10 * std::abs(mass_offset);
            if (std::abs(sum) > eps)
            {
                 if (ParallelDescriptor::IOProcessor())
                 {
                     std::cout << " ... current avg differs from mass_offset by " << sum << " " << '\n';
                     std::cout << " ... Gravity::actual_multilevel_solve -- total mass has changed!"
                               << '\n';;
                 }
            }

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << " ... subtracting " << sum << " to ensure solvability "
                          << '\n';

            for (int lev = 0; lev < num_levels; lev++)
                (*Rhs_p[lev]).plus(-sum, 0, 1, 0);
       }
    }

// *****************************************************************************

    for (int lev = 0; lev < num_levels; lev++)
    {
        Rhs_p[lev]->mult(Ggravity, 0, 1);
        Rhs_p[lev]->mult(a_inverse);
    }

    // FOR TIMINGS
    if (show_timings)
        ParallelDescriptor::Barrier();

    if (show_timings)
    {
        Real      end    = ParallelDescriptor::second() - strt_rhs;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity:: time in making rhs        = " << end << '\n';
#ifdef BL_LAZY
        });
#endif
    }

// *****************************************************************************

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                     : IntVect::TheZeroVector();

    //
    // Store the Dirichlet boundary condition for phi in bndry.
    //
    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level], dmap[level], 1, geom);
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;
    //
    // Build the homogeneous boundary conditions.  One could setVal
    // the bndry fabsets directly, but we instead do things as if
    // we had a fill-patched mf with grows--in that case the bndry
    // object knows how to grab grow data from the mf on physical
    // boundaries.  Here we create an mf, setVal, and pass that to
    // the bndry object.
    //
    if (level == 0)
    {
        bndry.setBndryValues(*phi_p[0], src_comp, dest_comp, num_comp,*phys_bc);
    }
    else
    {
        MultiFab CPhi;
        get_crse_phi(level, CPhi, time);
        BoxArray  crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes, dmap[level],
			      in_rad, out_rad, extent_rad, num_comp);
        crse_br.copyFrom(CPhi, CPhi.nGrow(), src_comp, dest_comp, num_comp);
        bndry.setBndryValues(crse_br, src_comp, LevelData[level]->get_new_data(PhiGrav_Type),
                             src_comp, dest_comp, num_comp, crse_ratio, *phys_bc);
    }

    // FOR TIMINGS
    if (show_timings)
        ParallelDescriptor::Barrier();
    const Real strt_solve = ParallelDescriptor::second();

    Real tol           = ml_tol;
    Real abs_tol       = 0;

    //
    // Can only use the C++ solvers if single-level
    //
    if (solve_with_cpp && (level == finest_level))
    {
        // We can only use the C++ solvers for a single level solve, but it's ok if level > 0
        solve_with_Cpp(level, *(phi_p[0]), grad_phi[0], *(Rhs_p[0]), tol, abs_tol);
    }
    else if ( solve_with_hpgmg && (level == finest_level) && (level == 0) )
    {
#ifdef USEHPGMG
        // Right now we can only use HPGMG for a single level = 0 solve
        solve_with_HPGMG(level, *(phi_p[0]), grad_phi[0], *(Rhs_p[0]), tol, abs_tol);
#endif
    }
    else
    {
        MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);
        mgt_solver.set_const_gravity_coeffs(xa, xb);

        Real final_resnorm;
	int always_use_bnorm = 0;
	int need_grad_phi = 1;

        //
        // Call the solver
        //
        mgt_solver.solve(phi_p, amrex::GetArrOfPtrs(Rhs_p),
			 bndry, tol, abs_tol, always_use_bnorm, final_resnorm, need_grad_phi);

        for (int lev = 0; lev < num_levels; lev++)
        {
            const Real* dx = parent->Geom(level+lev).CellSize();
            mgt_solver.get_fluxes(lev, grad_phi[level+lev], dx);
        }
    }
    Real end_solve = ParallelDescriptor::second() - strt_solve;

    if (show_timings)
    {
        ParallelDescriptor::ReduceRealMax(end_solve,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity:: time in solve          = " << end_solve << '\n';
    }

    // Average phi from fine to coarse level
    for (int lev = finest_level; lev > level; lev--)
    {
        if (is_new == 1)
        {
            amrex::average_down(LevelData[lev  ]->get_new_data(PhiGrav_Type),
                                 LevelData[lev-1]->get_new_data(PhiGrav_Type),
                                 0, 1, parent->refRatio(lev-1));

        }
        else if (is_new == 0)
        {
            amrex::average_down(LevelData[lev  ]->get_old_data(PhiGrav_Type),
                                 LevelData[lev-1]->get_old_data(PhiGrav_Type),
                                 0, 1, parent->refRatio(lev-1));
        }
    }

    // Average grad_phi from fine to coarse level
    for (int lev = finest_level; lev > level; lev--)
        average_fine_ec_onto_crse_ec(lev-1,is_new);

    if (show_timings)
    {
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Gravity:: all                  : time = " << end << '\n';
            std::cout << "Gravity:: all but solve        : time = " << end - end_solve << '\n';
        }
#ifdef BL_LAZY
        });
#endif
    }
}

void
Gravity::get_old_grav_vector (int       level,
                              MultiFab& grav_vector,
                              Real      time)
{
    BL_PROFILE("Gravity::get_old_grav_vector()");
    // Set to zero to fill ghost cells.
    grav_vector.setVal(0);

#ifdef CGRAV
    if (gravity_type == "StaticGrav")
    {
        make_prescribed_grav(level,time,grav_vector,0);
        grav_vector.FillBoundary();
    }
    else
    {
#endif

    // Fill grow cells in grad_phi, will need to compute grad_phi_cc in 1 grow cell
    const Geometry& geom = parent->Geom(level);
    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        {
            grad_phi_prev[level][i]->setBndry(0);
            grad_phi_prev[level][i]->FillBoundary(geom.periodicity());
        }
    }
    else
    {
        Array<std::unique_ptr<MultiFab> > crse_grad_phi(BL_SPACEDIM);
        get_crse_grad_phi(level, crse_grad_phi, time);
        fill_ec_grow(level, amrex::GetArrOfPtrs(grad_phi_prev[level]),
	                    amrex::GetArrOfPtrs(crse_grad_phi));
    }

    // Average edge-centered gradients to cell centers.
    amrex::average_face_to_cellcenter(grav_vector, 
				       amrex::GetArrOfConstPtrs(grad_phi_prev[level]),
				       geom);

#ifdef CGRAV
    if (gravity_type == "CompositeGrav")
    {
        make_prescribed_grav(level,time,grav_vector,1);
    }
#endif

    grav_vector.FillBoundary(geom.periodicity());

#ifdef CGRAV
    }
#endif

    MultiFab& G_old = LevelData[level]->get_old_data(Gravity_Type);

    // Fill G_old from grav_vector
    MultiFab::Copy(G_old, grav_vector, 0, 0, BL_SPACEDIM, 0);

    // This is a hack-y way to fill the ghost cell values of grav_vector
    //   before returning it
    AmrLevel* amrlev = &parent->getLevel(level);
    int ng = grav_vector.nGrow();
    AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
}

void
Gravity::get_new_grav_vector (int       level,
                              MultiFab& grav_vector,
                              Real      time)
{
    BL_PROFILE("Gravity::get_new_grav_vector()");
#ifdef CGRAV
    if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav")
#else
    if (gravity_type == "PoissonGrav")
#endif
    {
        // Set to zero to fill ghost cells
        grav_vector.setVal(0);

        // Fill grow cells in `grad_phi`, will need to compute `grad_phi_cc` in
        // 1 grow cell
        const Geometry& geom = parent->Geom(level);
        if (level == 0)
        {
            for (int i = 0; i < BL_SPACEDIM ; i++)
            {
                grad_phi_curr[level][i]->setBndry(0);
                grad_phi_curr[level][i]->FillBoundary(geom.periodicity());
            }
        }
        else
        {
            Array<std::unique_ptr<MultiFab> > crse_grad_phi(BL_SPACEDIM);
            get_crse_grad_phi(level, crse_grad_phi, time);
            fill_ec_grow(level, amrex::GetArrOfPtrs(grad_phi_curr[level]),
                                amrex::GetArrOfPtrs(crse_grad_phi));
        }

        // Average edge-centered gradients to cell centers, excluding grow cells
        amrex::average_face_to_cellcenter(grav_vector,
					   amrex::GetArrOfConstPtrs(grad_phi_curr[level]),
					   geom);

#ifdef CGRAV
        if (gravity_type == "CompositeGrav")
        {
          make_prescribed_grav(level,time,grav_vector,1);
        }
#endif

        grav_vector.FillBoundary(geom.periodicity());
    }

#ifdef CGRAV
    else if ( gravity_type == "StaticGrav")
    {
      make_prescribed_grav(level,time,grav_vector,0);
      grav_vector.FillBoundary();
    }
#endif

    MultiFab& G_new = LevelData[level]->get_new_data(Gravity_Type);

    // Fill G_new from grav_vector
    MultiFab::Copy(G_new, grav_vector, 0, 0, BL_SPACEDIM, 0);

    // This is a hack-y way to fill the ghost cell values of grav_vector
    //   before returning it
    AmrLevel* amrlev = &parent->getLevel(level) ;
    int ng = grav_vector.nGrow();
    AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
}

void
Gravity::add_to_fluxes(int level, int iteration, int ncycle)
{
    BL_PROFILE("Gravity::add_to_fluxes()");
    const int       finest_level      = parent->finestLevel();
    FluxRegister*   phi_fine          = (level<finest_level ? phi_flux_reg[level+1].get() : nullptr);
    FluxRegister*   phi_current       = (level>0 ? phi_flux_reg[level].get() : nullptr);
    const Geometry& geom              = parent->Geom(level);
    const Real*     dx                = geom.CellSize();
    const Real      area[BL_SPACEDIM] = { dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1] };

#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    int ngrow;
    fort_get_method_params(&ngrow);

    if (phi_fine)
    {
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            BoxArray ba = grids[level];
            ba.surroundingNodes(n);
            MultiFab fluxes(ba, dmap[level], 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(fluxes,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                FArrayBox& gphi_flux = fluxes[mfi];
                gphi_flux.copy((*grad_phi_curr[level][n])[mfi],tbx);
                gphi_flux.mult(area[n],tbx,0,1);
            }
            phi_fine->CrseInit(fluxes, n, 0, 0, 1, -1);
        }
    }

    if (phi_current && (iteration == ncycle))
    {
      MultiFab& phi_curr = LevelData[level]->get_new_data(PhiGrav_Type);
      for (MFIter mfi(phi_curr); mfi.isValid(); ++mfi)
      {
         for (int n=0; n<BL_SPACEDIM; ++n)
	     phi_current->FineAdd((*grad_phi_curr[level][n])[mfi],n,mfi.index(),0,0,1,area[n]);
      }
    }
}

void
Gravity::average_fine_ec_onto_crse_ec(int level, int is_new)
{
    BL_PROFILE("Gravity::average_fine_ec_to_crse_ec()");
    //
    // NOTE: this is called with level == the coarser of the two levels involved.
    //
    if (level == parent->finestLevel()) return;
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_gphi_fine_BA(grids[level+1].size());

    IntVect fine_ratio = parent->refRatio(level);

    for (int i = 0; i < crse_gphi_fine_BA.size(); ++i)
        crse_gphi_fine_BA.set(i, amrex::coarsen(grids[level+1][i],
                                                 fine_ratio));

    Array<std::unique_ptr<MultiFab> > crse_gphi_fine(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        const BoxArray eba = BoxArray(crse_gphi_fine_BA).surroundingNodes(n);
        crse_gphi_fine[n].reset(new MultiFab(eba, dmap[level+1], 1, 0));
    }

    auto& grad_phi = (is_new) ? grad_phi_curr : grad_phi_prev;

    amrex::average_down_faces(amrex::GetArrOfConstPtrs(grad_phi[level+1]),
			       amrex::GetArrOfPtrs(crse_gphi_fine), fine_ratio);

    const Geometry& cgeom = parent->Geom(level);

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
	grad_phi[level][n]->copy(*crse_gphi_fine[n], cgeom.periodicity());
    }
}

void
Gravity::reflux_phi (int       level,
                     MultiFab& dphi)
{
    BL_PROFILE("Gravity::reflux()");
    const Geometry& geom = parent->Geom(level);
    dphi.setVal(0);
    phi_flux_reg[level+1]->Reflux(dphi, 1.0, 0, 0, 1, geom);
}

void
Gravity::fill_ec_grow (int                     level,
                       const Array<MultiFab*>& ecF,
                       const Array<MultiFab*>& ecC) const
{
    BL_PROFILE("Gravity::fill_ec_grow()");
    //
    // Fill grow cells of the edge-centered mfs.  Assume
    // ecF built on edges of grids at this amr level, and ecC
    // is build on edges of the grids at amr level-1
    //
    BL_ASSERT(ecF.size() == BL_SPACEDIM);

    const int nGrow = ecF[0]->nGrow();

    if (nGrow == 0 || level == 0) return;

    BL_ASSERT(nGrow == ecF[1]->nGrow());
    BL_ASSERT(nGrow == ecF[2]->nGrow());

    const BoxArray& fgrids = grids[level];
    const Geometry& fgeom  = parent->Geom(level);

    BoxList bl = amrex::GetBndryCells(fgrids, 1);

    BoxArray f_bnd_ba(bl);

    bl.clear();

    BoxArray c_bnd_ba   = BoxArray(f_bnd_ba.size());
    IntVect  crse_ratio = parent->refRatio(level-1);

    for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
        c_bnd_ba.set(i, Box(f_bnd_ba[i]).coarsen(crse_ratio));
        f_bnd_ba.set(i, Box(c_bnd_ba[i]).refine(crse_ratio));
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        //
        // crse_src & fine_src must have same parallel distribution.
        // We'll use the KnapSack distribution for the fine_src_ba.
        // Since fine_src_ba should contain more points, this'll lead
        // to a better distribution.
        //
        BoxArray crse_src_ba(c_bnd_ba);
        BoxArray fine_src_ba(f_bnd_ba);

        crse_src_ba.surroundingNodes(n);
        fine_src_ba.surroundingNodes(n);

        Array<long> wgts(fine_src_ba.size());

        for (unsigned int i = 0; i < wgts.size(); i++)
        {
            wgts[i] = fine_src_ba[i].numPts();
        }
        DistributionMapping dm;
        //
        // This call doesn't invoke the MinimizeCommCosts() stuff.
        // There's very little to gain with these types of coverings
        // of trying to use SFC or anything else.
        // This also guarantees that these DMs won't be put into the
        // cache, as it's not representative of that used for more
        // usual MultiFabs.
        //
        dm.KnapSackProcessorMap(wgts, ParallelDescriptor::NProcs());

        MultiFab crse_src; crse_src.define(crse_src_ba, dm, 1, 0);
        MultiFab fine_src; fine_src.define(fine_src_ba, dm, 1, 0);

        crse_src.setVal(1.e200);
        fine_src.setVal(1.e200);
        //
        // We want to fill crse_src from ecC[n].
        //
	crse_src.copy(*ecC[n]);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
        {
            const int nComp = 1;
            const Box box = crse_src[mfi].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_PC_EDGE_INTERP, fort_pc_edge_interp)
                (box.loVect(), box.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(crse_src[mfi]), BL_TO_FORTRAN(fine_src[mfi]));
        }

        crse_src.clear();
        //
        // Replace pc-interpd fine data with preferred u_mac data at
        // this level u_mac valid only on surrounding faces of valid
        // region - this op will not fill grow region.
        //
        fine_src.copy(*ecF[n]); // parallel copy

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
        {
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine edges overlaying coarse edges.
            //
            const int nComp = 1;
            const Box& fbox = fine_src[mfi.index()].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_EDGE_INTERP, fort_edge_interp)
                (fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(fine_src[mfi]));
        }
        //
        // Build a mf with no grow cells on ecF grown boxes, do parallel copy into.
        //
        BoxArray fgridsG = ecF[n]->boxArray();
        fgridsG.grow(ecF[n]->nGrow());

        MultiFab ecFG(fgridsG, ecF[n]->DistributionMap(), 1, 0);

        ecFG.copy(fine_src); // Parallel copy

        fine_src.clear();

        ecFG.copy(*ecF[n]);   // Parallel copy

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*ecF[n]); mfi.isValid(); ++mfi)
            (*ecF[n])[mfi].copy(ecFG[mfi]);
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        ecF[n]->FillBoundary();
	ecF[n]->EnforcePeriodicity(fgeom.periodicity());
    }
}

void
Gravity::make_mg_bc ()
{
    BL_PROFILE("Gravity::make_mg_bc()");
    const Geometry& geom = parent->Geom(0);
    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        if (geom.isPeriodic(dir))
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
        else
        {
            if (phys_bc->lo(dir) == Symmetry)
            {
                mg_bc[2*dir + 0] = MGT_BC_NEU;
            }
            else if (phys_bc->lo(dir) == Outflow)
            {
                mg_bc[2*dir + 0] = MGT_BC_DIR;
            }
            else if (phys_bc->lo(dir) == Inflow)
            {
                mg_bc[2*dir + 0] = MGT_BC_DIR;
            }
            else
            {
                amrex::Abort("Unknown lo bc in make_mg_bc");
            }

            if (phys_bc->hi(dir) == Symmetry)
            {
                mg_bc[2*dir + 1] = MGT_BC_NEU;
            }
            else if (phys_bc->hi(dir) == Outflow)
            {
                mg_bc[2*dir + 1] = MGT_BC_DIR;
            }
            else if (phys_bc->hi(dir) == Inflow)
            {
                mg_bc[2*dir + 1] = MGT_BC_DIR;
            }
            else
            {
                amrex::Abort("Unknown hi bc in make_mg_bc");
            }
        }
    }
}

void
Gravity::set_mass_offset (Real time)
{
    BL_PROFILE("Gravity::set_mass_offset()");
    Real old_mass_offset = 0;

    int flev = parent->finestLevel();

    while (parent->getAmrLevels()[flev] == nullptr)
        flev--;

    if (flev > 0) old_mass_offset = mass_offset;

    mass_offset = 0;

    const Geometry& geom = parent->Geom(0);

    if (geom.isAllPeriodic())
    {
        //
        // Note: we must loop over levels because the volWgtSum routine zeros out
        //       crse regions under fine regions.
        //
        for (int lev = 0; lev <= flev; lev++)
        {
            Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(lev));

            BL_ASSERT(cs != 0);

#ifndef NO_HYDRO
            if (Nyx::Do_Hydro() == 1)
                mass_offset += cs->vol_weight_sum("density",time,true);
#endif
            for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                mass_offset += Nyx::theActiveParticles()[i]->sumParticleMass(lev);
        }

        mass_offset /= geom.ProbSize();

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "Gravity ... defining average density in Gravity::set_mass_offset to be " << mass_offset
                      << '\n';
    }

    const Real diff = std::abs(mass_offset - old_mass_offset);
    const Real eps  = 1.e-10 * std::abs(old_mass_offset);
    if (diff > eps && old_mass_offset > 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << " ... new vs old mass_offset " << mass_offset << " "
                      << old_mass_offset
                      << " ... diff is " << diff <<  '\n';
            std::cout << " ... Gravity::set_mass_offset -- total mass has changed!"
                      << '\n';;
        }
    }
}

void
Gravity::set_dirichlet_bcs (int       level,
                            MultiFab* phi)
{
    BL_PROFILE("Gravity::set_dirichlet_bcs()");
    const Real* dx        = parent->Geom(level).CellSize();
    const int*  domain_lo = parent->Geom(level).Domain().loVect();
    const int*  domain_hi = parent->Geom(level).Domain().hiVect();

    // Set phi to zero on all the ghost cells outside the domain.
    // If homogeneous bc's then we stop here; if not we add monopole bc's from each particle
    for (MFIter mfi(*phi); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        BL_FORT_PROC_CALL(FORT_SET_HOMOG_BCS, fort_set_homog_bcs)
            (lo, hi, domain_lo, domain_hi, BL_TO_FORTRAN((*phi)[mfi]), dx);
    }

    if (monopole_bcs && Nyx::theDMPC())
    {
        Array<Real> part_locs;
        Nyx::theDMPC()->GetParticleLocations(part_locs);

        Array<Real> part_mass;
        int start_comp = 0;
        int   num_comp = 1;
        Nyx::theDMPC()->GetParticleData(part_mass,start_comp,num_comp);

        // (x,y,z)
        int npart = part_locs.size()/3;

        for (MFIter mfi(*phi); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox(); const int* lo  = box.loVect(); const int* hi  = box.hiVect();
            BL_FORT_PROC_CALL(FORT_ADD_MONOPOLE_BCS, fort_add_monopole_bcs)
                (lo, hi, domain_lo, domain_hi, &npart,
                 part_locs.dataPtr(), part_mass.dataPtr(), BL_TO_FORTRAN((*phi)[mfi]), dx);
        }
    }
}

#ifdef CGRAV
void
Gravity::make_prescribed_grav (int       level,
                               Real      time,
                               MultiFab& grav_vector,
                               int       addToExisting)
{
    BL_PROFILE("Gravity::make_prescribed_grav()");
    const Real      strt = ParallelDescriptor::second();
    const Geometry& geom = parent->Geom(level);
    const Real*     dx   = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();

       BL_FORT_PROC_CALL(FORT_PRESCRIBE_GRAV,fort_prescribe_grav)
           (bx.loVect(), bx.hiVect(), dx,
            BL_TO_FORTRAN(grav_vector[mfi]),
            geom.ProbLo(), &addToExisting);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::make_prescribed_grav() time = " << end << '\n';
#ifdef BL_LAZY
        });
#endif
    }
}
#endif

void
Gravity::AddParticlesToRhs (int               level,
                            MultiFab&         Rhs,
                            int               ngrow)
{
    BL_PROFILE("Gravity::AddParticlesToRhs()");
    // Use the same multifab for all particle types
    MultiFab particle_mf(grids[level], dmap[level], 1, ngrow);

    for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
    {
        particle_mf.setVal(0.);
        Nyx::theActiveParticles()[i]->AssignDensitySingleLevel(particle_mf, level);
        MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
    }
}

void
Gravity::AddParticlesToRhs(int base_level, int finest_level, const Array<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddParticlesToRhsML()");
    const int num_levels = finest_level - base_level + 1;
    for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
    {
        Array<std::unique_ptr<MultiFab> > PartMF;
        Nyx::theActiveParticles()[i]->AssignDensity(PartMF, base_level, 1, finest_level);
        for (int lev = 0; lev < num_levels; lev++)
        {
            if (PartMF[lev]->contains_nan())
            {
                std::cout << "Testing particle density of type " << i << " at level " << base_level+lev << std::endl;
                amrex::Abort("...PartMF has NaNs in Gravity::actual_multilevel_solve()");
            }
        }

        for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        {
            amrex::average_down(*PartMF[lev+1], *PartMF[lev],
                                 0, 1, parent->refRatio(lev+base_level));
        }

        for (int lev = 0; lev < num_levels; lev++)
        {
            MultiFab::Add(*Rhs_particles[lev], *PartMF[lev], 0, 0, 1, 0);
        }
    }

}

void
Gravity::AddVirtualParticlesToRhs (int               level,
                                   MultiFab&         Rhs,
                                   int               ngrow)
{
    BL_PROFILE("Gravity::AddVirtualParticlesToRhs()");
    if (level <  parent->finestLevel())
    {
        // If we have virtual particles, add their density to the single level solve
        MultiFab particle_mf(grids[level], dmap[level], 1, ngrow);

        for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
        {
            particle_mf.setVal(0.);
            Nyx::theVirtualParticles()[i]->AssignDensitySingleLevel(particle_mf, level, 1, 1);
            MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddVirtualParticlesToRhs(int finest_level, const Array<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddVirtualParticlesToRhsML()");
    if (finest_level < parent->finestLevel())
    {
        // Should only need ghost cells for virtual particles if they're near
        // the simulation boundary and even then only maybe
        MultiFab VirtPartMF(grids[finest_level], dmap[finest_level], 1, 1);
        VirtPartMF.setVal(0.0);

        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            Nyx::theVirtualParticles()[i]->AssignDensitySingleLevel(VirtPartMF, finest_level, 1, 1);
            MultiFab::Add(*Rhs_particles[finest_level], VirtPartMF, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddGhostParticlesToRhs (int               level,
                                 MultiFab&         Rhs)
{
    BL_PROFILE("Gravity::AddGhostParticlesToRhs()");
    if (level > 0)
    {
        // If we have ghost particles, add their density to the single level solve
        MultiFab ghost_mf(grids[level], dmap[level], 1, 1);

        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            ghost_mf.setVal(0.);
            Nyx::theGhostParticles()[i]->AssignDensitySingleLevel(ghost_mf, level, 1, -1);
            MultiFab::Add(Rhs, ghost_mf, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddGhostParticlesToRhs(int level, const Array<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddGhostParticlesToRhsML()");
    if (level > 0)
    {
        // We require one ghost cell in GhostPartMF because that's how we handle
        // particles near fine-fine boundaries.  However we don't add any ghost
        // cells from GhostPartMF to the RHS.
        MultiFab GhostPartMF(grids[level], dmap[level], 1, 1);
        GhostPartMF.setVal(0.0);

        // Get the Ghost particle mass function. Note that Ghost particles should
        // only affect the coarsest level so we use a single level solve.  We pass in
        // -1 for the particle_lvl_offset because that makes the particles the size
        // of the coarse, not fine, dx.
        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            Nyx::theGhostParticles()[i]->AssignDensitySingleLevel(GhostPartMF, level, 1, -1);
            MultiFab::Add(*Rhs_particles[0], GhostPartMF, 0, 0, 1, 0);
        }
    }
}

void
Gravity::CorrectRhsUsingOffset(int level, MultiFab& Rhs)
{
    BL_PROFILE("Gravity::CorrectRhsUsingOffset()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << " ... subtracting average density from RHS in solve ... "
                  << mass_offset << '\n';

    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
        Rhs[mfi].plus(-mass_offset);
    // This checks if mass has been conserved--in particular if
    // virtual particles are correctly representing finer particles.
    if (level == 0)
    {
        Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(0)));
        Real sum = nyx_level->vol_weight_sum(Rhs,false);

        const Real eps = 1.e-10 * std::abs(mass_offset);
        if (std::abs(sum) > eps)
        {
            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << " ... current avg differs from mass_offset by " << sum << " " << '\n';
                std::cout << " ... Gravity : single_level solve for phi -- total mass has changed!"
                          << '\n';;
            }
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << " ... subtracting " << sum << " to ensure solvability "
                      << '\n';

        Rhs.plus(-sum, 0, 1, 0);
    }
}

void
Gravity::solve_with_Cpp(int level, MultiFab& soln, const Array<MultiFab*>& grad_phi,
                        MultiFab& rhs, Real tol, Real abs_tol)
{
  BL_PROFILE("Gravity::solve_with_Cpp()");
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();

  BndryData bd(grids[level], dmap[level], 1, geom);
  set_boundary(bd, rhs, dx);

  // Note that this actually solves Lap(phi) = RHS, not -Lap(phi) as in the F90 solve
  Laplacian lap_operator(bd, dx[0]);

  MultiGrid mg(lap_operator);
  mg.setVerbose(1);
  mg.solve(soln, rhs, tol, abs_tol);

  lap_operator.compFlux(*grad_phi[0],*grad_phi[1],*grad_phi[2],soln);

  // We have to multiply by -1 here because the compFlux routine returns
  // grad(phi), not -grad(phi) as in the F90 solver.
  grad_phi[0]->mult(-1.0);
  grad_phi[1]->mult(-1.0);
  grad_phi[2]->mult(-1.0);
}

#ifdef USEHPGMG
void
Gravity::solve_with_HPGMG(int level,
                          MultiFab& soln,
                          const Array<MultiFab*>& grad_phi,
                          MultiFab& rhs,
                          Real tol,
                          Real abs_tol)
{
  BL_PROFILE("Gravity::solve_with_HPGMG()");
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Box& domain = geom.Domain();

  BndryData bd(grids[level], dmap[level], 1, geom);
  set_boundary(bd, rhs, dx);
  const int n_cell = domain.length(0);

  // We'll get the max grid (box) size from the soln MultiFab already provided
  // to us. Just pick the first valid Box we can find and measure its size.
  // HPGMG requires the Boxes to be cubes, so if they're not then a sanity
  // check in MultiFab::CreateHPGMGLevel() will catch it and quit.
  MFIter mfi(soln);
  while (!mfi.isValid()) ++mfi;
  const Box& bx = mfi.validbox();
  const int max_grid_size = bx.length(0);

  const int my_rank = ParallelDescriptor::MyProc();
  const int num_ranks = ParallelDescriptor::NProcs();

  Laplacian lap_operator(bd, dx[0]);

  const Real a = 0.0; // coefficient in front of alpha in the Helmholtz operator
  // The canonical Helmholtz operator is a alpha u - b div (beta grad(u)) = f.
  // The self-gravity equation that we solve in Nyx is Lap(u) = f. So we need
  // either the betas or b to be -1. The other GMG solvers in BoxLib set the
  // b*beta to -1.
  const Real b = -1.0; // coefficient in front of beta in the Helmholtz operator

  const auto& ba = rhs.boxArray();
  const auto& dm = rhs.DistributionMap()
  MultiFab alpha(ba, dm, 1, 0);
  MultiFab beta_cc(ba, dm, 1, 1);
  alpha.setVal(0.0);
  beta_cc.setVal(1.0);

  const int domain_boundary_condition = BC_PERIODIC;

  const int minCoarseDim = 2; // avoid problems with black box calculation of D^{-1} for poisson with periodic BC's on a 1^3 grid

  level_type level_h;
  mg_type MG_h;

  int numVectors = 12;

  CreateHPGMGLevel(&level_h, rhs, n_cell, max_grid_size, my_rank, num_ranks, domain_boundary_condition, numVectors, *dx);
  SetupHPGMGCoefficients(a, b, alpha, beta_cc, &level_h);
  ConvertToHPGMGLevel(rhs, n_cell, max_grid_size, &level_h, VECTOR_F);

  if (level_h.boundary_condition.type == BC_PERIODIC)
  {
    Real average_value_of_f = mean (&level_h, VECTOR_F);
    if (average_value_of_f != 0.0)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        std::cerr << "WARNING: Periodic boundary conditions, but f does not sum to zero... mean(f)=" << average_value_of_f << std::endl;
        std::cerr << "Subtracting mean(f) from RHS ..." << std::endl;
      }
      shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rebuild_operator(&level_h,NULL,a,b);    // i.e. calculate Dinv and lambda_max
  MGBuild(&MG_h,&level_h,a,b,minCoarseDim,ParallelDescriptor::Communicator()); // build the Multigrid Hierarchy
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ParallelDescriptor::IOProcessor())
    std::cout << std::endl << std::endl << "===== STARTING SOLVE =====" << std::endl << std::flush;

  MGResetTimers (&MG_h);
  zero_vector (MG_h.levels[0], VECTOR_U);
  #ifdef USE_FCYCLES
  FMGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, abs_tol, tol);
  #else
  MGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, abs_tol, tol);
  #endif /* USE_FCYCLES */

  if (show_timings)
    MGPrintTiming (&MG_h, 0);

  // Now convert solution from HPGMG back to rhs MultiFab.
  ConvertFromHPGMGLevel(soln, &level_h, VECTOR_U);

  MGDestroy(&MG_h); // destroys all but the finest grid
  destroy_level(&level_h); // destroys the finest grid

  lap_operator.compFlux(*grad_phi[0],*grad_phi[1],*grad_phi[2],soln);

  // We have to multiply by -1 here because the compFlux routine returns
  // grad(phi), not -grad(phi) as in the F90 solver.
  grad_phi[0]->mult(-1.0);
  grad_phi[1]->mult(-1.0);
  grad_phi[2]->mult(-1.0);
}
#endif

void
Gravity::set_boundary(BndryData& bd, MultiFab& rhs, const Real* dx)
{
  BL_PROFILE("Gravity::set_boundary()");
  for (int n=0; n<BL_SPACEDIM; ++n) {
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index();

      // Our default will be that the face of this grid is either touching another grid
      //  across an interior boundary or a periodic boundary.
      {
        // Define the type of boundary conditions to be Dirichlet (even for periodic)
        bd.setBoundCond(Orientation(n, Orientation::low) ,i,0,LO_DIRICHLET);
        bd.setBoundCond(Orientation(n, Orientation::high),i,0,LO_DIRICHLET);

        // Set the boundary conditions to the cell centers outside the domain
        bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
        bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
      }
    }
  }
}



// Routine to duplicate Gravity class data onto sidecars
void
Gravity::AddProcsToComp(Amr *aptr, int level, AmrLevel *level_data_to_install,
                        int ioProcNumSCS, int ioProcNumAll, int scsMyId, MPI_Comm scsComm)
{
   parent = aptr;


   // ---- pack up the ints
   Array<int> allInts;

   if(scsMyId == ioProcNumSCS) {
     allInts.push_back(density);
     allInts.push_back(finest_level);
     allInts.push_back(finest_level_allocated);
     allInts.push_back(verbose);
     allInts.push_back(show_timings);
     allInts.push_back(no_sync);
     allInts.push_back(no_composite);
     allInts.push_back(dirichlet_bcs);
     allInts.push_back(monopole_bcs);
     allInts.push_back(solve_with_cpp);
     allInts.push_back(solve_with_hpgmg);
     allInts.push_back(stencil_type);
     for(int i(0); i < 2*BL_SPACEDIM; ++i)    { allInts.push_back(mg_bc[i]); }
   }

   amrex::BroadcastArray(allInts, scsMyId, ioProcNumSCS, scsComm);

   // ---- unpack the ints
   if(scsMyId != ioProcNumSCS) {
     int count(0);

     density = allInts[count++];
     finest_level = allInts[count++];
     finest_level_allocated = allInts[count++];
     verbose = allInts[count++];
     show_timings = allInts[count++];
     no_sync = allInts[count++];
     no_composite = allInts[count++];
     dirichlet_bcs = allInts[count++];
     monopole_bcs = allInts[count++];
     solve_with_cpp = allInts[count++];
     solve_with_hpgmg = allInts[count++];
     stencil_type = allInts[count++];
     for(int i(0); i < 2*BL_SPACEDIM; ++i)    { mg_bc[i] = allInts[count++]; }

     BL_ASSERT(count == allInts.size());
   }


   // ---- pack up the Reals
   Array<Real> allReals;
   if(scsMyId == ioProcNumSCS) {
     allReals.push_back(mass_offset);
     allReals.push_back(sl_tol);
     allReals.push_back(ml_tol);
     allReals.push_back(delta_tol);
     allReals.push_back(Ggravity);
   }

   amrex::BroadcastArray(allReals, scsMyId, ioProcNumSCS, scsComm);
   amrex::BroadcastArray(level_solver_resnorm, scsMyId, ioProcNumSCS, scsComm);

   // ---- unpack the Reals
   if(scsMyId != ioProcNumSCS) {
     int count(0);
     mass_offset = allReals[count++];
     sl_tol = allReals[count++];
     ml_tol = allReals[count++];
     delta_tol = allReals[count++];
     Ggravity = allReals[count++];

     BL_ASSERT(count == allReals.size());
   }


   // ---- pack up the strings
   Array<std::string> allStrings;
   Array<char> serialStrings;
   if(scsMyId == ioProcNumSCS) {
     allStrings.push_back(gravity_type);
     serialStrings = amrex::SerializeStringArray(allStrings);
   }

   amrex::BroadcastArray(serialStrings, scsMyId, ioProcNumSCS, scsComm);

   // ---- unpack the strings
   if(scsMyId != ioProcNumSCS) {
     int count(0);
     allStrings = amrex::UnSerializeStringArray(serialStrings);
     gravity_type = allStrings[count++];
   }


   // ---- BCRec
   Array<int> bcrLo(BL_SPACEDIM), bcrHi(BL_SPACEDIM);
   if(scsMyId == ioProcNumSCS) {
     for(int i(0); i < bcrLo.size(); ++i) { bcrLo[i] = phys_bc->lo(i); }
     for(int i(0); i < bcrHi.size(); ++i) { bcrHi[i] = phys_bc->hi(i); }
   }
   ParallelDescriptor::Bcast(bcrLo.dataPtr(), bcrLo.size(), ioProcNumSCS, scsComm);
   ParallelDescriptor::Bcast(bcrHi.dataPtr(), bcrHi.size(), ioProcNumSCS, scsComm);
   if(scsMyId != ioProcNumSCS) {
     for(int i(0); i < bcrLo.size(); ++i) { phys_bc->setLo(i, bcrLo[i]); }
     for(int i(0); i < bcrHi.size(); ++i) { phys_bc->setHi(i, bcrHi[i]); }
   }


   // ---- MultiFabs

   // ---- ---- grad_phi_curr :: Array< Array<std::unique_ptr<MultiFab> > > grad_phi_curr;
   if(scsMyId != ioProcNumSCS) {
     for(int j(0); j < grad_phi_curr.size(); ++j) {
       grad_phi_curr[j].clear();
     }
   }
   int gpcSize(grad_phi_curr.size());
   ParallelDescriptor::Bcast(&gpcSize, 1, ioProcNumSCS, scsComm);
   if(scsMyId != ioProcNumSCS) {
     grad_phi_curr.resize(gpcSize);
   }

   for(int j(0); j < grad_phi_curr.size(); ++j) {
     Array<int> isDefined;

     if(scsMyId == ioProcNumSCS) {
       isDefined.resize(grad_phi_curr[j].size());
       for(int i(0); i < grad_phi_curr[j].size(); ++i) {
	   isDefined[i] = (grad_phi_curr[j][i] != nullptr);
       }
     }
     amrex::BroadcastArray(isDefined, scsMyId, ioProcNumAll, scsComm);
     if(isDefined.size() > 0) {
       BL_ASSERT(isDefined.size() == BL_SPACEDIM);
       if(scsMyId != ioProcNumSCS) {
	   grad_phi_curr[j].resize(isDefined.size());
         for(int i(0); i < grad_phi_curr[j].size(); ++i) {
           if(isDefined[i]) {
             grad_phi_curr[j][i].reset(new MultiFab);
	   }
         }
       }

       for(int i(0); i < grad_phi_curr[j].size(); ++i) {
         if(grad_phi_curr[j][i]) {
           grad_phi_curr[j][i]->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
         }
       }
     }
   }


   // ---- ---- grad_phi_prev :: Array< Array<std::unique_ptr<MultiFab> > > grad_phi_prev;
   if(scsMyId != ioProcNumSCS) {
     for(int j(0); j < grad_phi_prev.size(); ++j) {
       grad_phi_prev[j].clear();
     }
   }
   int gppSize(grad_phi_prev.size());
   ParallelDescriptor::Bcast(&gppSize, 1, ioProcNumSCS, scsComm);
   if(scsMyId != ioProcNumSCS) {
     grad_phi_prev.resize(gppSize);
   }

   for(int j(0); j < grad_phi_prev.size(); ++j) {
     Array<int> isDefined;

     if(scsMyId == ioProcNumSCS) {
       isDefined.resize(grad_phi_prev[j].size());
       for(int i(0); i < grad_phi_prev[j].size(); ++i) {
	   isDefined[i] = (grad_phi_prev[j][i] != nullptr);
       }
     }
     amrex::BroadcastArray(isDefined, scsMyId, ioProcNumAll, scsComm);
     if(isDefined.size() > 0) {
       BL_ASSERT(isDefined.size() == BL_SPACEDIM);
       if(scsMyId != ioProcNumSCS) {
         grad_phi_prev[j].resize(isDefined.size());
         for(int i(0); i < grad_phi_prev[j].size(); ++i) {
           if(isDefined[i]) {
             grad_phi_prev[j][i].reset(new MultiFab);
	   }
         }
       }
       for(int i(0); i < grad_phi_prev[j].size(); ++i) {
         if(grad_phi_prev[j][i]) {
           grad_phi_prev[j][i]->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
         }
       }
     }
   }


   // ---- FluxRegisters :: Array<std::unique_ptr<FluxRegister> > phi_flux_reg;
   if(scsMyId != ioProcNumSCS) {
     phi_flux_reg.clear();
   }
   int pfrSize(phi_flux_reg.size());
   ParallelDescriptor::Bcast(&pfrSize, 1, ioProcNumSCS, scsComm);
   if(scsMyId != ioProcNumSCS) {
     phi_flux_reg.resize(pfrSize);
   }

   Array<int> isDefined;

   if(scsMyId == ioProcNumSCS) {
     isDefined.resize(phi_flux_reg.size());
     for(int i(0); i < phi_flux_reg.size(); ++i) {
	 isDefined[i] = (phi_flux_reg[i] != nullptr);
     }
   }
   amrex::BroadcastArray(isDefined, scsMyId, ioProcNumAll, scsComm);
   if(isDefined.size() > 0) {
     if(scsMyId != ioProcNumSCS) {
       phi_flux_reg.resize(isDefined.size());
       for(int i(0); i < phi_flux_reg.size(); ++i) {
         if(isDefined[i]) {
           phi_flux_reg[i].reset(new FluxRegister);
	 }
       }
     }

     for(int i(0); i < phi_flux_reg.size(); ++i) {
       if(phi_flux_reg[i]) {
         phi_flux_reg[i]->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
       }
     }
   }



   // ---- LevelData
    LevelData[level] = level_data_to_install;
}


