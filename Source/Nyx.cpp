#include <winstd.H>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <CONSTANTS.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <Particles_F.H>
#include <Utility.H>

#if BL_USE_MPI
#include "MemInfo.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef REEBER
#include <ReeberAnalysis.H>
#endif

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#endif

extern "C" {
  int get_comp_urho();
  int get_comp_temp();
  int get_comp_e_int();
}

const int NyxHaloFinderSignal = 42;
const int GimletSignal = 55;

static int sum_interval = -1;
static Real fixed_dt    = -1.0;
static Real initial_dt  = -1.0;
static Real dt_cutoff   =  0;

int Nyx::strict_subcycling = 0;

Real Nyx::old_a      = -1.0;
Real Nyx::new_a      = -1.0;
Real Nyx::old_a_time = -1.0;
Real Nyx::new_a_time = -1.0;

Array<Real> Nyx::plot_z_values;
Array<Real> Nyx::analysis_z_values;

bool Nyx::dump_old = false;
int Nyx::verbose      = 0;
int Nyx::show_timings = 0;

Real Nyx::cfl = 0.8;
Real Nyx::init_shrink = 1.0;
Real Nyx::change_max  = 1.1;

BCRec Nyx::phys_bc;
int Nyx::do_reflux = 1;
int Nyx::NUM_STATE = -1;
int Nyx::NUM_GROW = -1;

int Nyx::nsteps_from_plotfile = -1;

ErrorList Nyx::err_list;

int Nyx::Density = -1;
int Nyx::Eden = -1;
int Nyx::Eint = -1;
int Nyx::Xmom = -1;
int Nyx::Ymom = -1;
int Nyx::Zmom = -1;

int Nyx::Temp_comp = -1;
int Nyx::  Ne_comp = -1;

int Nyx::NumSpec  = 0;
int Nyx::NumAux   = 0;
int Nyx::NumAdv   = 0;

int Nyx::FirstSpec = -1;
int Nyx::FirstAux  = -1;
int Nyx::FirstAdv  = -1;

Real Nyx::small_dens = -1.e200;
Real Nyx::small_temp = -1.e200;
Real Nyx::gamma      =  0;

int Nyx::do_hydro = -1;
int Nyx::add_ext_src = 0;
int Nyx::heat_cool_type = 0;
int Nyx::strang_split = 0;

Real Nyx::average_gas_density = 0;
Real Nyx::average_dm_density = 0;
Real Nyx::average_neutr_density = 0;
Real Nyx::average_total_density = 0;

// Real Nyx::ave_lev_vorticity[10];
// Real Nyx::std_lev_vorticity[10];

#ifdef GRAVITY
Gravity* Nyx::gravity  =  0;
int Nyx::do_grav       = -1;
#else
int Nyx::do_grav       =  0;
#endif

int Nyx::allow_untagging    = 0;
int Nyx::use_const_species  = 0;
int Nyx::normalize_species  = 0;
int Nyx::do_special_tagging = 0;
int Nyx::ppm_type           = 1;
int Nyx::ppm_reference      = 1;
int Nyx::corner_coupling    = 1;

int Nyx::use_colglaz        = 0;
int Nyx::version_2          = 0;

int Nyx::use_flattening     = 1;
int Nyx::ppm_flatten_before_integrals = 0;

Real Nyx:: h_species        = 0.0;
Real Nyx::he_species        = 0.0;

int Nyx::use_exact_gravity  = 0;

#ifdef _OPENMP
#include <omp.h>
#endif

int Nyx::write_parameters_in_plotfile = true;
int Nyx::print_fortran_warnings       = true;

// Do we use separate SPH particles to initialize
//  the density and momentum on the grid?
int  Nyx::init_with_sph_particles = 0;

// Do we write the particles in single (IEEE32)
//  or doublue (NATIVE) precision?
#ifdef BL_SINGLE_PRECISION_PARTICLES
std::string Nyx::particle_plotfile_format = "IEEE32";
#else
std::string Nyx::particle_plotfile_format = "NATIVE";
#endif

// this will be reset upon restart
Real         Nyx::previousCPUTimeUsed = 0.0;

Real         Nyx::startCPUTime = 0.0;

int reeber_int(0);
int gimlet_int(0);

int Nyx::forceParticleRedist = false;
int Nyx::nSidecarProcs(0);

// Note: Nyx::variableSetUp is in Nyx_setup.cpp
void
Nyx::variable_cleanup ()
{
#ifdef GRAVITY
    if (gravity != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting gravity in variable_cleanup...\n";
        delete gravity;
        gravity = 0;
    }
#endif

    desc_lst.clear();
}

void
Nyx::read_params ()
{
    BL_PROFILE("Nyx::read_params()");
    static bool done = false;

    if (done) return;  // (caseywstark) when would this happen?

    done = true;  // ?

    ParmParse pp("nyx");

    pp.query("v", verbose);
    pp.query("show_timings", show_timings);
    //verbose = (verbose ? 1 : 0);
    pp.get("init_shrink", init_shrink);
    pp.get("cfl", cfl);
    pp.query("change_max", change_max);
    pp.query("fixed_dt", fixed_dt);
    pp.query("initial_dt", initial_dt);
    pp.query("sum_interval", sum_interval);
    pp.query("do_reflux", do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.get("dt_cutoff", dt_cutoff);

    pp.query("dump_old", dump_old);

    pp.query("small_dens", small_dens);
    pp.query("small_temp", small_temp);
    pp.query("gamma", gamma);

    pp.query("strict_subcycling",strict_subcycling);

    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, BL_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i, lo_bc[i]);
        phys_bc.setHi(i, hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
        }
    }

    pp.get("do_hydro", do_hydro);
#ifdef NO_HYDRO
    if (do_hydro == 1)
        BoxLib::Error("Cant have do_hydro == 1 when NO_HYDRO is true");
#endif

#ifdef NO_HYDRO
#ifndef GRAVITY
        BoxLib::Error("Dont know what to do with both hydro and gravity off");
#endif
#endif

    pp.query("add_ext_src", add_ext_src);
    pp.query("strang_split", strang_split);

    pp.query("heat_cool_type", heat_cool_type);

    pp.query("use_exact_gravity", use_exact_gravity);

#ifdef HEATCOOL
    if (heat_cool_type > 0 && add_ext_src == 0)
       BoxLib::Error("Nyx::must set add_ext_src to 1 if heat_cool_type > 0");
    if (heat_cool_type != 1 && heat_cool_type != 3)
       BoxLib::Error("Nyx:: nonzero heat_cool_type must equal 1 or 3");
    if (heat_cool_type == 0)
       BoxLib::Error("Nyx::contradiction -- HEATCOOL is defined but heat_cool_type == 0");
#else
    if (heat_cool_type > 0)
       BoxLib::Error("Nyx::you set heat_cool_type > 0 but forgot to set USE_HEATCOOL = TRUE");
#endif

    pp.query("allow_untagging", allow_untagging);
    pp.query("use_const_species", use_const_species);
    pp.query("normalize_species", normalize_species);
    pp.query("ppm_type", ppm_type);
    pp.query("ppm_reference", ppm_reference);
    pp.query("ppm_flatten_before_integrals", ppm_flatten_before_integrals);
    pp.query("use_flattening", use_flattening);
    pp.query("use_colglaz", use_colglaz);
    pp.query("version_2", version_2);
    pp.query("corner_coupling", corner_coupling);

    if (do_hydro == 1)
    {
        if (do_hydro == 1 && use_const_species == 1)
        {
           pp.get("h_species" ,  h_species);
           pp.get("he_species", he_species);
           BL_FORT_PROC_CALL(SET_XHYDROGEN,set_xhydrogen)(h_species);
           if (ParallelDescriptor::IOProcessor())
           {
               std::cout << "Nyx::setting species concentrations to "
                         << h_species << " and " << he_species
                         << " in the hydro and in the EOS " << std::endl;
           }
        }

        //
        if (use_colglaz == 1)
        {
           if (ppm_type == 0 && ParallelDescriptor::IOProcessor())
               std::cout << "Nyx::setting use_colglaz = 1 with ppm_type = 0 \n";
           if (ppm_type != 0)
               BoxLib::Error("Nyx::ppm_type must be 0 with use_colglaz = 1");
        }

        // ppm_flatten_before_integrals is only done for ppm_type != 0
        if (ppm_type == 0 && ppm_flatten_before_integrals > 0)
        {
            std::cerr << "ppm_flatten_before_integrals > 0 not implemented for ppm_type != 0 \n";
            BoxLib::Error();
        }

        if (version_2 > 0 && ppm_type == 0)
           BoxLib::Error("Nyx::version_2 only defined for ppm_type = 1");

        if (version_2 !=0 && version_2 != 1 && version_2 != 2 && version_2 != 3)
           BoxLib::Error("Nyx:: don't know what to do with version_2 flag");

        // Make sure ppm_type is set correctly.
        if (ppm_type != 0 && ppm_type != 1 && ppm_type != 2)
        {
           BoxLib::Error("Nyx::ppm_type must be 0, 1 or 2");
        }
    }

    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

#ifdef GRAVITY
    pp.get("do_grav", do_grav);
#endif

    read_particle_params();

    read_init_params();

    pp.query("write_parameter_file",write_parameters_in_plotfile);
    pp.query("print_fortran_warnings",print_fortran_warnings);

    read_comoving_params();

    if (pp.contains("plot_z_values"))
    {
      int num_z_values = pp.countval("plot_z_values");
      plot_z_values.resize(num_z_values);
      pp.queryarr("plot_z_values",plot_z_values,0,num_z_values);
    }

    if (pp.contains("analysis_z_values"))
    {
      int num_z_values = pp.countval("analysis_z_values");
      analysis_z_values.resize(num_z_values);
      pp.queryarr("analysis_z_values",analysis_z_values,0,num_z_values);
    }

    pp.query("gimlet_int", gimlet_int);
}

Nyx::Nyx ()
{
    BL_PROFILE("Nyx::Nyx()");
#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
    }
#endif
    fine_mask = 0;
}

Nyx::Nyx (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time)
{
    BL_PROFILE("Nyx::Nyx(Amr)");
    build_metrics();
    fine_mask = 0;

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    // Initialize to zero here in case we run with do_grav = false.
    MultiFab& new_grav_mf = get_new_data(Gravity_Type);
    new_grav_mf.setVal(0);

    if (do_grav)
    {
        // gravity is a static object, only alloc if not already there
        if (gravity == 0) {
          gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);
	}

        gravity->install_level(level, this);

        if (verbose && level == 0 && ParallelDescriptor::IOProcessor()) {
            std::cout << "Setting the gravity type to "
                      << gravity->get_gravity_type() << '\n';
	}
   }
#endif

    // Initialize the "a" variable
    if (level == 0 && time == 0.0 && old_a_time < 0.)
    {
       old_a_time = 0.0;
       new_a_time = 0.0;

       old_a = 1.0 / (1.0 + initial_z);
       new_a = old_a;
    }

     // Initialize "this_z" in the atomic_rates_module
     if (heat_cool_type == 1 || heat_cool_type == 3)
         BL_FORT_PROC_CALL(INIT_THIS_Z, init_this_z)(&old_a);

    // Set grav_n_grow to 3 on init. It'll be reset in advance.
    grav_n_grow = 3;
}

Nyx::~Nyx ()
{
#ifndef NO_HYDRO
    if (do_hydro == 1)
        delete flux_reg;
#endif
    delete fine_mask;
}

void
Nyx::restart (Amr&     papa,
              istream& is,
              bool     b_read_special)
{
    BL_PROFILE("Nyx::restart()");
    AmrLevel::restart(papa, is, b_read_special);

    build_metrics();


    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      // get ellapsed CPU time
      std::ifstream CPUFile;
      std::string FullPathCPUFile = parent->theRestartFile();
      FullPathCPUFile += "/CPUtime";
      CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

      CPUFile >> previousCPUTimeUsed;
      CPUFile.close();

      std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";
    }

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        BL_ASSERT(flux_reg == 0);
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    if (do_grav && level == 0)
    {
        BL_ASSERT(gravity == 0);
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);
    }
#endif
}

void
Nyx::build_metrics ()
{
}

void
Nyx::setTimeLevel (Real time,
                   Real dt_old,
                   Real dt_new)
{
    if (ParallelDescriptor::IOProcessor()) {
       std::cout << "Setting the current time in the state data to "
                 << parent->cumTime() << std::endl;
    }
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
Nyx::init (AmrLevel& old)
{
    BL_PROFILE("Nyx::init(old)");
    Nyx* old_level = (Nyx*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time = old_level->state[PhiGrav_Type].curTime();
    Real prev_time = old_level->state[PhiGrav_Type].prevTime();
#else
    Real cur_time = old_level->state[State_Type].curTime();
    Real prev_time = old_level->state[State_Type].prevTime();
#endif

    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);

        for (FillPatchIterator
                 fpi(old, S_new, 0, cur_time,   State_Type, 0, NUM_STATE),
                dfpi(old, D_new, 0, cur_time, DiagEOS_Type, 0, 2);
                fpi.isValid() && dfpi.isValid();
                ++fpi,++dfpi)
        {
            FArrayBox&  tmp =  fpi();
            FArrayBox& dtmp = dfpi();
            S_new[fpi].copy(tmp);
            D_new[fpi].copy(dtmp);
        }
    }
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    for (FillPatchIterator fpi(old, Phi_new, 0, cur_time, PhiGrav_Type, 0, 1);
         fpi.isValid(); ++fpi)
    {
        Phi_new[fpi].copy(fpi());
    }
#endif

    // Set E in terms of e + kinetic energy
    // if (do_hydro)
    // enforce_consistent_e(S_new);
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Nyx::init ()
{
    BL_PROFILE("Nyx::init()");
    Real dt        = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time  = get_level(level-1).state[PhiGrav_Type].curTime();
    Real prev_time = get_level(level-1).state[PhiGrav_Type].prevTime();
#else
    Real cur_time  = get_level(level-1).state[State_Type].curTime();
    Real prev_time = get_level(level-1).state[State_Type].prevTime();
#endif
    Real dt_old    = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time, dt_old, dt);

#ifndef NO_HYDRO
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    FillCoarsePatch(S_new, 0, cur_time,   State_Type, 0, S_new.nComp());
    FillCoarsePatch(D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    FillCoarsePatch(Phi_new, 0, cur_time, PhiGrav_Type, 0, Phi_new.nComp());
#endif

    // We set dt to be large for this new level to avoid screwing up
    // computeNewDt.
    parent->setDtLevel(1.e100, level);

    // Set E in terms of e + kinetic energy
    // if (do_hydro)
    // enforce_consistent_e(S_new);
}

Real
Nyx::initial_time_step ()
{
    BL_PROFILE("Nyx::initial_time_step()");
    Real dummy_dt = 0;
    Real init_dt = 0;

    if (initial_dt > 0)
    {
        init_dt = initial_dt;
    }
    else
    {
        init_dt = init_shrink * est_time_step(dummy_dt);
    }

    bool dt_changed = false;
    if (level == 0 && plot_z_values.size() > 0)
        plot_z_est_time_step(init_dt,dt_changed);

    if (level == 0 && analysis_z_values.size() > 0)
        analysis_z_est_time_step(init_dt,dt_changed);

    return init_dt;
}

Real
Nyx::est_time_step (Real dt_old)
{
    BL_PROFILE("Nyx::est_time_step()");
    if (fixed_dt > 0)
        return fixed_dt;

    // This is just a dummy value to start with
    Real est_dt = 1.0e+200;

#ifndef NO_HYDRO
    const MultiFab& stateMF = get_new_data(State_Type);
#endif

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

#ifndef NO_HYDRO
    if (do_hydro)
    {
        Real a = get_comoving_a(cur_time);
        const Real* dx = geom.CellSize();

	Real dt = est_dt;

#ifdef _OPENMP
#pragma omp parallel firstprivate(dt)
#endif
	{
	  for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
	    {
	      const Box& box = mfi.tilebox();

	      BL_FORT_PROC_CALL(FORT_ESTDT, fort_estdt)
                (BL_TO_FORTRAN(stateMF[mfi]), box.loVect(), box.hiVect(), dx,
                 &dt, &a);
	    }
#ifdef _OPENMP
#pragma omp critical (nyx_estdt)
#endif
	  {
	    est_dt = std::min(est_dt, dt);
	  }
	}

        // If in comoving coordinates, then scale dt (based on u and c) by a
        est_dt *= a;

        ParallelDescriptor::ReduceRealMin(est_dt);
        est_dt *= cfl;
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "...estdt from hydro at level "
                      << level << ": "
                      << est_dt << '\n';
    }
#endif

#ifdef GRAVITY
    particle_est_time_step(est_dt);
#endif

    if (level==0)
        comoving_est_time_step(cur_time,est_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Nyx::est_time_step at level "
                  << level
                  << ":  estdt = "
                  << est_dt << '\n';

    return est_dt;
}

void
Nyx::computeNewDt (int                   finest_level,
                   int                   sub_cycle,
                   Array<int>&           n_cycle,
                   const Array<IntVect>& ref_ratio,
                   Array<Real>&          dt_min,
                   Array<Real>&          dt_level,
                   Real                  stop_time,
                   int                   post_regrid_flag)
{
    BL_PROFILE("Nyx::computeNewDt()");
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Nyx& adv_level = get_level(i);
        dt_min[i] = adv_level.est_time_step(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
        if (post_regrid_flag == 1)
        {
            //
            // Limit dt's by pre-regrid dt
            //
            for (i = 0; i <= finest_level; i++)
            {
                dt_min[i] = std::min(dt_min[i], dt_level[i]);
            }
            //
            // Find the minimum over all levels
            //
            for (i = 0; i <= finest_level; i++)
            {
                n_factor *= n_cycle[i];
                dt_0 = std::min(dt_0, n_factor * dt_min[i]);
            }
        }
        else
        {
            bool sub_unchanged=true;
            if ((parent->maxLevel() > 0) && (level == 0) &&
                (parent->subcyclingMode() == "Optimal") &&
                (parent->okToRegrid(level) || parent->levelSteps(0) == 0) )
            {
                int new_cycle[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                    new_cycle[i] = n_cycle[i];
                // The max allowable dt
                Real dt_max[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    dt_max[i] = dt_min[i];
                }
                // find the maximum number of cycles allowed.
                int cycle_max[finest_level+1];
                cycle_max[0] = 1;
                for (i = 1; i <= finest_level; i++)
                {
                    cycle_max[i] = parent->MaxRefRatio(i-1);
                }
                // estimate the amout of work to advance each level.
                Real est_work[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    est_work[i] = parent->getLevel(i).estimateWork();
                }
                // this value will be used only if the subcycling pattern is changed.
                dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
                for (i = 0; i <= finest_level; i++)
                {
                    if (n_cycle[i] != new_cycle[i])
                    {
                        sub_unchanged = false;
                        n_cycle[i] = new_cycle[i];
                    }
                }

            }

            if (sub_unchanged)
            //
            // Limit dt's by change_max * old dt
            //
            {
                for (i = 0; i <= finest_level; i++)
                {
                    if (verbose && ParallelDescriptor::IOProcessor())
                    {
                        if (dt_min[i] > change_max*dt_level[i])
                        {
                            std::cout << "Nyx::compute_new_dt : limiting dt at level "
                                      << i << '\n';
                            std::cout << " ... new dt computed: " << dt_min[i]
                                      << '\n';
                            std::cout << " ... but limiting to: "
                                      << change_max * dt_level[i] << " = " << change_max
                                      << " * " << dt_level[i] << '\n';
                        }
                    }

                    dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
                }
                //
                // Find the minimum over all levels
                //
                for (i = 0; i <= finest_level; i++)
                {
                    n_factor *= n_cycle[i];
                    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
                }
            }
            else
            {
                if (verbose && ParallelDescriptor::IOProcessor())
                {
                   std::cout << "Nyx: Changing subcycling pattern. New pattern:\n";
                   for (i = 1; i <= finest_level; i++)
                    std::cout << "   Lev / n_cycle: " << i << " " << n_cycle[i] << '\n';
                }
            }
        }
    }
    else
    {
        dt_0 = fixed_dt;
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    // Shrink the time step if necessary in order to hit the next plot_z_value
    if (level == 0 && ( plot_z_values.size() > 0 || analysis_z_values.size() > 0 ) )
    {
        bool dt_changed_plot     = false;
        bool dt_changed_analysis = false;
   
        if (plot_z_values.size() > 0)
           plot_z_est_time_step(dt_0,dt_changed_plot);

        if (analysis_z_values.size() > 0)
           analysis_z_est_time_step(dt_0,dt_changed_analysis);

        // Update the value of a if we didn't change dt in the call to plot_z_est_time_step or analysis_z_est_time_step.
        // If we didn't change dt there, then we have already done the integration.
        // If we did    change dt there, then we need to re-integrate here.
        if (dt_changed_plot || dt_changed_analysis)
            integrate_comoving_a(cur_time,dt_0);
    }
    else 
    {
        integrate_comoving_a(cur_time,dt_0);
    }
     
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void
Nyx::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Array<int>&           n_cycle,
                       const Array<IntVect>& ref_ratio,
                       Array<Real>&          dt_level,
                       Real                  stop_time)
{
    BL_PROFILE("Nyx::computeInitialDt()");
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    if (parent->subcyclingMode() == "Optimal")
    {
        int new_cycle[finest_level+1];
        for (i = 0; i <= finest_level; i++)
            new_cycle[i] = n_cycle[i];
        Real dt_max[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            dt_max[i] = get_level(i).initial_time_step();
        }
        // Find the maximum number of cycles allowed
        int cycle_max[finest_level+1];
        cycle_max[0] = 1;
        for (i = 1; i <= finest_level; i++)
        {
            cycle_max[i] = parent->MaxRefRatio(i-1);
        }
        // estimate the amout of work to advance each level.
        Real est_work[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            est_work[i] = parent->getLevel(i).estimateWork();
        }
        dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
        for (i = 0; i <= finest_level; i++)
        {
            n_cycle[i] = new_cycle[i];
        }
        if (verbose && ParallelDescriptor::IOProcessor() && finest_level > 0)
        {
           std::cout << "Nyx: Initial subcycling pattern:\n";
           for (i = 0; i <= finest_level; i++)
               std::cout << "Level " << i << ": " << n_cycle[i] << '\n';
        }
    }
    else
    {
        for (i = 0; i <= finest_level; i++)
        {
            dt_level[i] = get_level(i).initial_time_step();
            n_factor *= n_cycle[i];
            dt_0 = std::min(dt_0, n_factor * dt_level[i]);
        }
    }
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }

    integrate_comoving_a(cur_time,dt_0);
}

bool
Nyx::writePlotNow ()
{
    BL_PROFILE("Nyx::writePlotNow()");
    if (level > 0)
        BoxLib::Error("Should only call writePlotNow at level 0!");

    bool found_one = false;

    if (plot_z_values.size() > 0)
    {

#ifdef NO_HYDRO
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
#else
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#endif
        Real a_old = get_comoving_a(prev_time);
        Real z_old = (1. / a_old) - 1.;

        Real a_new = get_comoving_a( cur_time);
        Real z_new = (1. / a_new) - 1.;

        for (int i = 0; i < plot_z_values.size(); i++)
        {
            if (std::abs(z_new - plot_z_values[i]) < (0.01 * (z_old - z_new)) )
                found_one = true;
        }
    }

    if (found_one) {
        return true;
    } else {
        return false;
    }
}

bool
Nyx::doAnalysisNow ()
{
    BL_PROFILE("Nyx::doAnalysisNow()");
    if (level > 0)
        BoxLib::Error("Should only call doAnalysisNow at level 0!");

    bool found_one = false;

    if (analysis_z_values.size() > 0)
    {

#ifdef NO_HYDRO
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
#else
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#endif
        Real a_old = get_comoving_a(prev_time);
        Real z_old = (1. / a_old) - 1.;

        Real a_new = get_comoving_a( cur_time);
        Real z_new = (1. / a_new) - 1.;

        for (int i = 0; i < analysis_z_values.size(); i++)
        {
            if (std::abs(z_new - analysis_z_values[i]) < (0.01 * (z_old - z_new)) )
                found_one = true;
        }
    }

    if (found_one) {
        return true;
    } else {
        return false;
    }
}

void
Nyx::do_energy_diagnostics ()
{
    // nothing to see here, folks
}

void
Nyx::post_timestep (int iteration)
{
    BL_PROFILE("Nyx::post_timestep()");
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    const int ncycle = parent->nCycle(level);

    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration
    //
    if (iteration == ncycle)
        remove_ghost_particles();

    //
    // Redistribute if it is not the last subiteration
    //
    if (iteration < ncycle || level == 0)
    {
         for (int i = 0; i < theActiveParticles().size(); i++)
         {
             theActiveParticles()[i]->Redistribute(false, true, level, grav_n_grow);
         }
    }

#ifndef NO_HYDRO
    if (do_reflux && level < finest_level)
    {
        MultiFab& S_new_crse = get_new_data(State_Type);
#ifdef GRAVITY
        MultiFab drho_and_drhoU;
#ifdef CGRAV
        if (do_grav &&
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav"))
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            // Define the update to rho and rhoU due to refluxing.
            drho_and_drhoU.define(grids, BL_SPACEDIM + 1, 0, Fab_allocate);
            MultiFab::Copy(drho_and_drhoU, S_new_crse, Density, 0,
                           BL_SPACEDIM + 1, 0);
            drho_and_drhoU.mult(-1.0);
        }
#endif // GRAVITY

        //We must reflux if the next finer level is subcycled relative to this level;
        //   otherwise the reflux was done as part of the multilevel advance
        if (parent->nCycle(level+1) != 1)
          reflux();

        // We need to do this before anything else because refluxing changes the
        // values of coarse cells underneath fine grids with the assumption
        // they'll be over-written by averaging down
        if (level < finest_level)
            average_down();

        // This needs to be done after any changes to the state from refluxing.
        enforce_nonnegative_species(S_new_crse);

#ifdef GRAVITY
#ifdef CGRAV
        if (do_grav &&
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
            && gravity->get_no_sync() == 0)
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->get_no_sync() == 0)
#endif
        {
            MultiFab::Add(drho_and_drhoU, S_new_crse, Density, 0, BL_SPACEDIM+1, 0);

            MultiFab dphi(grids, 1, 0);
            dphi.setVal(0);

            gravity->reflux_phi(level, dphi);

            // Compute (cross-level) gravity sync based on drho, dphi
            PArray<MultiFab> grad_delta_phi_cc(finest_level - level + 1,
                                               PArrayManage);
            for (int lev = level; lev <= finest_level; lev++)
            {
                grad_delta_phi_cc.set(lev - level,
                                      new MultiFab(get_level(lev).boxArray(),BL_SPACEDIM, 0));
                grad_delta_phi_cc[lev-level].setVal(0);
            }

            gravity->gravity_sync(level,finest_level,iteration,ncycle,drho_and_drhoU,dphi,grad_delta_phi_cc);
            dphi.clear();

            for (int lev = level; lev <= finest_level; lev++)
            {
                Real dt_lev = parent->dtLevel(lev);
                MultiFab&  S_new_lev = get_level(lev).get_new_data(State_Type);
                Real cur_time = state[State_Type].curTime();
                Real a_new = get_comoving_a(cur_time);

                const BoxArray& ba = get_level(lev).boxArray();
                MultiFab grad_phi_cc(ba, BL_SPACEDIM, 0);
                gravity->get_new_grav_vector(lev, grad_phi_cc, cur_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
		  FArrayBox sync_src;
		  FArrayBox dstate;

		  for (MFIter mfi(S_new_lev,true); mfi.isValid(); ++mfi)
                  {
                    const Box& bx = mfi.tilebox();
                    dstate.resize(bx, BL_SPACEDIM + 1);
                    if (lev == level)
                    {
		      dstate.copy(drho_and_drhoU[mfi]);
                    }
                    else
                    {
		      dstate.setVal(0);
                    }

                    // Compute sync source
                    sync_src.resize(bx, BL_SPACEDIM+1);
                    int i = mfi.index();
                    BL_FORT_PROC_CALL(FORT_SYNCGSRC,fort_syncgsrc)
                        (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grad_phi_cc[i]),
                         BL_TO_FORTRAN(grad_delta_phi_cc[lev-level][i]),
                         BL_TO_FORTRAN(S_new_lev[i]), BL_TO_FORTRAN(dstate),
                         BL_TO_FORTRAN(sync_src), &a_new, dt_lev);

                    sync_src.mult(0.5 * dt_lev);
                    S_new_lev[mfi].plus(sync_src, 0, Xmom, BL_SPACEDIM);
                    S_new_lev[mfi].plus(sync_src, BL_SPACEDIM, Eden, 1);
		  }
		}
	    }
        }
#endif
    }
#endif // end ifndef NO_HYDRO

    if (level < finest_level)
        average_down();

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
#ifndef NO_HYDRO
        if ( (do_hydro == 1) && (sum_interval > 0) && (nstep % sum_interval == 0))
        {
            sum_integrated_quantities();
        }
#endif
        write_info();

#if BL_USE_MPI
        // Memory monitoring:
        MemInfo* mInfo = MemInfo::GetInstance();
        char info[32];
        snprintf(info, sizeof(info), "Step %4d", nstep);
        mInfo->LogSummary(info);
#endif
    }

#ifndef NO_HYDRO
    if (do_hydro)
    {
       // Re-compute temperature after all the other updates.
       compute_new_temp();
    }
#endif
}

void
Nyx::post_restart ()
{
    BL_PROFILE("Nyx::post_restart()");
    if (level == 0)
        particle_post_restart(parent->theRestartFile());

    if (level == 0)
        comoving_a_post_restart(parent->theRestartFile());

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    // Update the value of a only if restarting from chk00000
    //   (special case for which computeNewDt is *not* called from Amr::coarseTimeStep)
    if (level == 0 && cur_time == 0.0)
        integrate_comoving_a(cur_time,parent->dtLevel(0));

#ifdef TISF
     int blub = parent->finestLevel();
     BL_FORT_PROC_CALL(FORT_SET_FINEST_LEVEL, fort_set_finest_level)(&blub);
#endif

#ifdef GRAVITY

    if (do_grav)
    {
        if (level == 0)
        {
            for (int lev = 0; lev <= parent->finestLevel(); lev++)
            {
                AmrLevel& this_level = get_level(lev);
                gravity->install_level(lev, &this_level);
            }

            gravity->set_mass_offset(cur_time);

            if (
#ifdef CGRAV
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
#else
	    gravity->get_gravity_type() == "PoissonGrav"
#endif
)
            {
                // Do multilevel solve here.  We now store phi in the checkpoint file so we can use it
                //  at restart.
                int use_previous_phi_as_guess = 1;
                gravity->multilevel_solve_for_phi(0,parent->finestLevel(),use_previous_phi_as_guess);

#ifndef AGN
                if (do_dm_particles)
#endif
                {
                    for (int k = 0; k <= parent->finestLevel(); k++)
                    {
                        const BoxArray& ba = get_level(k).boxArray();
                        MultiFab grav_vec_new(ba, BL_SPACEDIM, 0);
                        gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
                    }
                }
            }
        }
    }
#endif

#ifndef NO_HYDRO
    if (level == 0)
    {
       // Need to compute this *before* regridding in case this is needed
       compute_average_density();
       set_small_values();
    }
#endif
}

#ifndef NO_HYDRO
void
Nyx::set_small_values ()
{
       if (do_hydro == 0) {
          return;
       }

       Real small_pres;

       const Real cur_time = state[State_Type].curTime();
       Real a = get_comoving_a(cur_time);

       Real average_temperature;
       compute_average_temperature(average_temperature);
       //
       // Get the number of species from the network model.
       //
       BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);
       BL_FORT_PROC_CALL(GET_NUM_AUX , get_num_aux )(&NumAux);

       BL_FORT_PROC_CALL(SET_SMALL_VALUES, set_small_values)
            (&average_gas_density, &average_temperature,
             &a,  &small_dens, &small_temp, &small_pres);

       if (verbose && ParallelDescriptor::IOProcessor())
       {
          std::cout << "... setting small_dens to " << small_dens << '\n';
          std::cout << "... setting small_temp to " << small_temp << '\n';
          std::cout << "... setting small_pres to " << small_pres << '\n';
       }
}
#endif

void
Nyx::postCoarseTimeStep (Real cumtime)
{
   BL_PROFILE("Nyx::postCoarseTimeStep()");

   AmrLevel::postCoarseTimeStep(cumtime);

   const Real cur_time = state[State_Type].curTime();
   const int whichSidecar(0);

#ifdef REEBER
   const auto& reeber_density_var_list = getReeberHaloDensityVars();
   bool do_analysis(doAnalysisNow());
   if (do_analysis || (reeber_int > 0 && nStep() % reeber_int == 0)) {
     if (ParallelDescriptor::NProcsSidecar(0) <= 0) { // we have no sidecars, so do everything in situ
       const Real time1 = ParallelDescriptor::second();

       BoxArray ba;
       DistributionMapping dm;
       getAnalysisDecomposition(Geom(), ParallelDescriptor::NProcs(), ba, dm);
       MultiFab reeberMF(ba, reeber_density_var_list.size() + 1, 0, dm);
       int cnt = 1;
       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           MultiFab *derive_dat = particle_derive(*it, cur_time, 0); // FIXME: Is this the right way? 
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
           delete derive_dat;
           cnt++;
       }

       reeberMF.setVal(0, 0, 1, 0);
       for (int comp = 1; comp < reeberMF.nComp(); ++comp)
           MultiFab::Add(reeberMF, reeberMF, comp, 0, 1, 0);

       std::vector<Halo> reeber_halos;
       runReeberAnalysis(reeberMF, Geom(), nStep(), do_analysis, &reeber_halos);

       // Redistribute halos to "correct" processor for simulation
       // FIXME: This is a hack that maps things to MultiFabs and back. This should be
       // changed to use Nyx's particle redistribution infrastructure.
       reeberMF.setVal(0, 0, reeber_density_var_list.size() + 1, 0);
       // Deposit halo into MultiFab. This works because (i) There is only one box
       // per processors and halos returned by Reeber will always be in that box;
       // (ii) Halos have different positions; (iii) Halos have a mass that differs from
       // zero.
       BL_ASSERT(dm[ParallelDescriptor::MyProc()] == ParallelDescriptor::MyProc());
       FArrayBox& my_fab = reeberMF[ParallelDescriptor::MyProc()];
       for (const Halo& h : reeber_halos)
       {
           BL_ASSERT(reeberMF.fabbox(ParallelDescriptor::MyProc()).contains(h.position));
           my_fab(h.position, 0) = h.totalMass;
           for (int comp = 0; comp < reeber_density_var_list.size(); ++comp)
           {
               my_fab(h.position, comp + 1) = h.individualMasses[comp];
           }
       }
       // Actual redistribution
       //MultiFab redistributeFab(m_leveldata.boxArray(), reeber_density_var_list.size() + 1, 0, m_leveldata.DistributionMap());
       const MultiFab& simMF = get_new_data(State_Type);
       const BoxArray& simBA = simMF.boxArray();
       const DistributionMapping& simDM = simMF.DistributionMap();

       MultiFab redistributeFab(simBA, reeber_density_var_list.size() + 1, 0, simDM);
       redistributeFab.copy(reeberMF);
       // Re-extract halos
       reeber_halos.clear();
       for (MFIter mfi(redistributeFab); mfi.isValid(); ++mfi)
       {
           const Box& currBox = mfi.fabbox();
           for (IntVect iv = currBox.smallEnd(); iv <= currBox.bigEnd(); currBox.next(iv))
           {
               Real totalMass = redistributeFab[mfi](iv, 0);
               if (totalMass > 0)
               {
                   std::vector<Real> masses(reeber_density_var_list.size(), 0);
                   for (int comp = 0; comp < reeber_density_var_list.size(); ++comp)
                   {
                       masses[comp] = redistributeFab[mfi](iv, comp + 1);
                   }
                   reeber_halos.emplace_back(iv, totalMass, masses);
               }
           }
        }
       // NOTE: ZARIJA, GET YOUR FRESH HALOS HERE!!!
#if 0
       std::ofstream os(BoxLib::Concatenate(BoxLib::Concatenate("debug-halos-", nStep(), 5), ParallelDescriptor::MyProc(), 2));
       for (const Halo& h : reeber_halos)
       {
           os << h.position << " " << h.totalMass << std::endl;
       }
#endif

       const Real time2 = ParallelDescriptor::second();
       if (ParallelDescriptor::IOProcessor())
       {
         std::cout << std::endl << "===== Time to post-process: " << time2 - time1 << " sec" << std::endl;
       }
     } else { // we have sidecars, so do everything in-transit
       int sidecarSignal(NyxHaloFinderSignal);
       const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
       ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       Geometry geom(Geom());
       Geometry::SendGeometryToSidecar(&geom, whichSidecar);

       MultiFab reeberMF(grids, reeber_density_var_list.size(), 0);
       int cnt = 0;
       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           MultiFab *derive_dat = particle_derive(*it, cur_time, 0); // FIXME: Is this the right way?
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
           delete derive_dat;
           cnt++;
       }

       int time_step(nStep()), nComp(reeberMF.nComp());

       Real time1(ParallelDescriptor::second());
       ParallelDescriptor::Bcast(&nComp, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       MultiFab *mfSource = &reeberMF;
       MultiFab *mfDest = 0;
       int srcComp(0), destComp(1);
       int srcNGhost(0), destNGhost(0);
       MPI_Comm commSrc(ParallelDescriptor::CommunicatorComp());
       MPI_Comm commDest(ParallelDescriptor::CommunicatorSidecar());
       MPI_Comm commInter(ParallelDescriptor::CommunicatorInter(whichSidecar));
       MPI_Comm commBoth(ParallelDescriptor::CommunicatorBoth(whichSidecar));
       bool isSrc(true);

       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost,
                           commSrc, commDest, commInter, commBoth,
                           isSrc);


       ParallelDescriptor::Bcast(&time_step, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       int do_analysis_bcast(do_analysis);
       ParallelDescriptor::Bcast(&do_analysis_bcast, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       Real time2(ParallelDescriptor::second());
       if(ParallelDescriptor::IOProcessor()) {
         std::cout << "COMPUTE PROCESSES: time spent sending data to sidecars: " << time2 - time1 << std::endl;
       }
     }
   }
#endif /* REEBER */

#ifdef GIMLET
   if ( (nStep() % gimlet_int == 0) || doAnalysisNow() ) 
   {
     if (ParallelDescriptor::NProcsSidecar(0) <= 0) { // we have no sidecars, so do everything in situ
       const Real time1 = ParallelDescriptor::second();

       const MultiFab &state_data = get_data(State_Type, cur_time);
       const MultiFab &eos_data = get_data(DiagEOS_Type, cur_time);

       // Grab only components of state and EOS data that gimlet needs.
       MultiFab density (state_data.boxArray(), 1, 0);
       density.copy(state_data, get_comp_urho()-1, 0, 1);
       MultiFab temperature (state_data.boxArray(), 1, 0);
       temperature.copy(eos_data, get_comp_temp()-1, 0, 1);
       MultiFab e_int(state_data.boxArray(), 1, 0);
       e_int.copy(state_data, get_comp_e_int()-1, 0, 1);

       MultiFab *dm_density = particle_derive("particle_mass_density", cur_time, 0);
       MultiFab *xmom = particle_derive("xmom", cur_time, 0);
       MultiFab *ymom = particle_derive("ymom", cur_time, 0);
       MultiFab *zmom = particle_derive("zmom", cur_time, 0);
       Geometry geom(Geom());
       Real omega_m, omega_b, comoving_h;
       BL_FORT_PROC_CALL(GET_OMM, get_omm)(&omega_m);
       BL_FORT_PROC_CALL(GET_OMB, get_omb)(&omega_b);
       omega_b *= omega_m;
       BL_FORT_PROC_CALL(GET_HUBBLE, get_hubble)(&comoving_h);
       const double omega_l = 1.0 - omega_m;

       do_analysis(omega_b, omega_m, omega_l, comoving_h, new_a, density,
                   temperature, e_int, *dm_density, *xmom, *ymom, *zmom, geom, nStep());

       const Real time_to_postprocess = ParallelDescriptor::second() - time1;
       if (ParallelDescriptor::IOProcessor())
         std::cout << std::endl << "===== Time to post-process: " << time_to_postprocess << " sec" << std::endl;
     } else { // we have sidecars, so do everything in-transit
       int signal = GimletSignal;
       const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
       ParallelDescriptor::Bcast(&signal, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       const MultiFab &state_data = get_data(State_Type, cur_time);
       const MultiFab &eos_data = get_data(DiagEOS_Type, cur_time);

       // Grab only components of state and EOS data that gimlet needs.
       MultiFab density (state_data.boxArray(), 1, 0);
       density.copy(state_data, get_comp_urho()-1, 0, 1);
       MultiFab temperature (state_data.boxArray(), 1, 0);
       temperature.copy(eos_data, get_comp_temp()-1, 0, 1);
       MultiFab e_int(state_data.boxArray(), 1, 0);
       e_int.copy(state_data, get_comp_e_int()-1, 0, 1);

       MultiFab *dm_density = particle_derive("particle_mass_density", cur_time, 0);
       MultiFab *xmom = particle_derive("xmom", cur_time, 0);
       MultiFab *ymom = particle_derive("ymom", cur_time, 0);
       MultiFab *zmom = particle_derive("zmom", cur_time, 0);
       Geometry geom(Geom());
       int time_step = nStep();
       Real omega_m, omega_b, comoving_h;
       BL_FORT_PROC_CALL(GET_OMM, get_omm)(&omega_m);
       BL_FORT_PROC_CALL(GET_OMB, get_omb)(&omega_b);
       BL_FORT_PROC_CALL(GET_HUBBLE, get_hubble)(&comoving_h);
       const Real time1 = ParallelDescriptor::second();

       MultiFab *mfSource = 0;
       MultiFab *mfDest = 0;
       int srcComp(0), destComp(0), nComp(1);
       int srcNGhost(0), destNGhost(0);
       const MPI_Comm &commSrc = ParallelDescriptor::CommunicatorComp();
       const MPI_Comm &commDest = ParallelDescriptor::CommunicatorSidecar();
       const MPI_Comm &commInter = ParallelDescriptor::CommunicatorInter(whichSidecar);
       const MPI_Comm &commBoth = ParallelDescriptor::CommunicatorBoth(whichSidecar);
       bool isSrc(true);

       BL_ASSERT(density.boxArray() == dm_density->boxArray());  // ---- need to check them all
       BoxArray::SendBoxArray(density.boxArray(), whichSidecar);

       mfSource = &density;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = &temperature;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = &e_int;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = dm_density;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = xmom;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = ymom;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       mfSource = zmom;
       MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

       Geometry::SendGeometryToSidecar(&geom, whichSidecar);

       ParallelDescriptor::Bcast(&Nyx::new_a, 1, MPI_IntraGroup_Broadcast_Rank, commInter);
       ParallelDescriptor::Bcast(&omega_m, 1, MPI_IntraGroup_Broadcast_Rank, commInter);
       omega_b *= omega_m;
       ParallelDescriptor::Bcast(&omega_b, 1, MPI_IntraGroup_Broadcast_Rank, commInter);
       ParallelDescriptor::Bcast(&comoving_h, 1, MPI_IntraGroup_Broadcast_Rank, commInter);
       ParallelDescriptor::Bcast(&time_step, 1, MPI_IntraGroup_Broadcast_Rank, commInter);
       const Real time2 = ParallelDescriptor::second();
       if (ParallelDescriptor::IOProcessor()) {
         std::cout << "COMPUTE PROCESSES: time spent sending data to sidecars: " << time2 - time1 << std::endl;
       }

       delete dm_density;
       delete xmom;
       delete ymom;
       delete zmom;
     }
   }
#endif /* GIMLET */

    //
    // postCoarseTimeStep() is only called by level 0.
    //
    if (Nyx::theDMPC() && particle_move_type == "Random")
        particle_move_random();
}

void
Nyx::post_regrid (int lbase,
                  int new_finest)
{
    BL_PROFILE("Nyx::post_regrid()");
#ifndef NO_HYDRO
#ifdef TISF
     BL_FORT_PROC_CALL(FORT_SET_FINEST_LEVEL, fort_set_finest_level)(&new_finest);
#endif
#endif

    if (level == lbase) {
        particle_redistribute(lbase, forceParticleRedist);
    }

    int which_level_being_advanced = parent->level_being_advanced();

#ifdef GRAVITY
    bool do_grav_solve_here;
    if (which_level_being_advanced >= 0)
    {
        do_grav_solve_here = (level == which_level_being_advanced) && (lbase == which_level_being_advanced);
    } else {
        do_grav_solve_here = (level == lbase);
    }

    // Only do solve here if we will be using it in the timestep right after without re-solving,
    //      or if this is called from somewhere other than Amr::timeStep
    const Real cur_time = state[PhiGrav_Type].curTime();
    if (do_grav && (cur_time > 0) && do_grav_solve_here)
    {
#ifdef CGRAV
        if (gravity->get_gravity_type() == "PoissonGrav" || gravity->get_gravity_type() == "CompositeGrav")
#else
        if (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            int use_previous_phi_as_guess = 1;
            gravity->multilevel_solve_for_phi(level, new_finest, use_previous_phi_as_guess);
        }
    }
#endif
    delete fine_mask;
    fine_mask = 0;
}

void
Nyx::post_init (Real stop_time)
{
    BL_PROFILE("Nyx::post_init()");
    if (level > 0) {
        return;
    }

    // If we restarted from a plotfile, we need to reset the level_steps counter
    if ( ! parent->theRestartPlotFile().empty()) {
        parent->setLevelSteps(0,nsteps_from_plotfile);
    }

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; --k) {
        get_level(k).average_down();
    }

#ifdef GRAVITY
    if (do_grav)
    {
        const Real cur_time = state[PhiGrav_Type].curTime();
        if
#ifdef CGRAV
            (gravity->get_gravity_type() == "PoissonGrav" ||
             gravity->get_gravity_type() == "CompositeGrav")
#else
	    (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            //
            // Calculate offset before first multilevel solve.
            //
            gravity->set_mass_offset(cur_time);

            //
            // Solve on full multilevel hierarchy
            //
            gravity->multilevel_solve_for_phi(0, finest_level);
        }

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= finest_level; k++)
        {
            const BoxArray& ba = get_level(k).boxArray();
            MultiFab grav_vec_new(ba, BL_SPACEDIM, 0);
            gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
        }
    }
#endif

#ifndef NO_HYDRO
    if ( (do_hydro == 1) && (sum_interval > 0) && (parent->levelSteps(0) % sum_interval == 0) )
    {
        sum_integrated_quantities();
    }
    else
    {
        // Even if we don't call `sum_integrated_quantities` we need to compute
        // average_density before regridding
        compute_average_density();
    }

    if (do_hydro == 1)
    {
        set_small_values();
    }
#endif

    write_info();
}

int
Nyx::okToContinue ()
{
    if (level > 0) {
        return 1;
    }

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff) {
        test = 0;
    }

    if ((test == 1) && (final_a > 0))
    {
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
        Real a = get_comoving_a(cur_time);
        if (a >= final_a) test = 0;
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            if (test == 0) {
                std::cout << "...a " << a
                          << " is greater than or equal to final_a " << final_a
                          << '\n';
	    }
        }
    }
    return test;
}

#ifdef AUX_UPDATE
void
Nyx::advance_aux (Real time,
                  Real dt)
{
    BL_PROFILE("Nyx::advance_aux()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... special update for auxiliary variables \n";

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
        BL_FORT_PROC_CALL(FORT_AUXUPDATE, fort_auxupdate)
            (BL_TO_FORTRAN(old_fab), BL_TO_FORTRAN(new_fab), box.loVect(),
             box.hiVect(), &dt);
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::reflux ()
{
    BL_PROFILE("Nyx::reflux()");
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = ParallelDescriptor::second();

    get_flux_reg(level+1).Reflux(get_new_data(State_Type), 1.0, 0, 0, NUM_STATE,
                                 geom);

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end, IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Nyx::reflux() at level "
                      << level
                      << " : time = "
                      << end << '\n';
    }
}
#endif // NO_HYDRO

void
Nyx::average_down ()
{
    BL_PROFILE("Nyx::average_down()");
    if (level == parent->finestLevel()) return;

#ifndef NO_HYDRO
    // With State_Type we do DiagEOS_Type
    average_down(State_Type);
#endif

#ifdef GRAVITY
    average_down(PhiGrav_Type);
    average_down(Gravity_Type);
#endif
}

#ifndef NO_HYDRO
void
Nyx::enforce_nonnegative_species (MultiFab& S_new)
{
    BL_PROFILE("Nyx::enforce_nonnegative_species()");
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(ENFORCE_NONNEGATIVE_SPECIES,
			  enforce_nonnegative_species)
	  (BL_TO_FORTRAN(S_new[mfi]), bx.loVect(), bx.hiVect(),
	   &print_fortran_warnings);
    }
}

void
Nyx::enforce_consistent_e (MultiFab& S)
{
    BL_PROFILE("Nyx::enforce_consistent_e()");
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();
        BL_FORT_PROC_CALL(ENFORCE_CONSISTENT_E,
			  enforce_consistent_e)
	  (lo, hi, BL_TO_FORTRAN(S[mfi]));
    }
}
#endif

void
Nyx::average_down (int state_index)
{
    BL_PROFILE("Nyx::average_down(si)");
#ifndef NO_HYDRO
    // We average DiagEOS_Type when average_down is called with State_Type
    if (state_index == DiagEOS_Type) return;
#endif

    if (level == parent->finestLevel()) return;

    Nyx& fine_level = get_level(level + 1);

    MultiFab& S_crse = get_new_data(state_index);
    MultiFab& S_fine = fine_level.get_new_data(state_index);
    const int num_comps = S_fine.nComp();

#ifndef NO_HYDRO
    if (state_index == State_Type)
    {
        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA(S_fine.boxArray().size());

        for (int i = 0; i < S_fine.boxArray().size(); ++i)
        {
            crse_S_fine_BA.set(i, BoxLib::coarsen(S_fine.boxArray()[i],
                                                  fine_ratio));
        }
        MultiFab crse_S_fine(crse_S_fine_BA, num_comps, 0);

        MultiFab& D_crse =            get_new_data(DiagEOS_Type);
        MultiFab& D_fine = fine_level.get_new_data(DiagEOS_Type);
        const int num_D_comps = D_fine.nComp();
        MultiFab crse_D_fine(crse_S_fine_BA, num_D_comps, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
        {
 	    const Box&       overlap  = mfi.tilebox();

            const FArrayBox& fine_S_fab = S_fine[mfi];
            FArrayBox&       crse_S_fab = crse_S_fine[mfi];

            const FArrayBox& fine_D_fab = D_fine[mfi];
            FArrayBox&       crse_D_fab = crse_D_fine[mfi];

            BL_FORT_PROC_CALL(FORT_AVGDOWN, fort_avgdown)
                (BL_TO_FORTRAN(crse_S_fab), num_comps,
		 BL_TO_FORTRAN(fine_S_fab),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());

            BL_FORT_PROC_CALL(FORT_AVGDOWN, fort_avgdown)
                (BL_TO_FORTRAN(crse_D_fab), num_D_comps,
		 BL_TO_FORTRAN(fine_D_fab),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());
        }
        D_crse.copy(crse_D_fine);
        S_crse.copy(crse_S_fine);
    }
    else
#endif
    {
      const Geometry& fine_geom = parent->Geom(level+1);
      const Geometry& crse_geom = parent->Geom(level  );
      BoxLib::average_down(S_fine,S_crse,fine_geom,crse_geom,0,num_comps,fine_ratio);
    }
}

void
Nyx::errorEst (TagBoxArray& tags,
               int          clearval,
               int          tagval,
               Real         time,
               int          n_error_buf,
               int          ngrow)
{
    BL_PROFILE("Nyx::errorEst()");
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    Real avg;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        if (err_list[j].errType() == ErrorRec::UseAverage)
        {
            if (err_list[j].name() == "density")
            {
                avg = average_gas_density;
            }
            else if (err_list[j].name() == "particle_mass_density")
            {
                avg = average_dm_density;
#ifdef NEUTRINO_PARTICLES
                avg += average_neutr_density;
#endif
            }
            else if (err_list[j].name() == "total_density")
            {
                avg = average_total_density;
            }
#if 0
            else if (err_list[j].name() == "magvort")
            {
                avg = std::fabs(ave_lev_vorticity[level]);
                stddev = std_lev_vorticity[level];
                thresh = avg + std::max(stddev,avg);
                //std::cout << "errorEst, level " << level << ": magvort avg " << avg << ", stddev " << stddev
                //        << ", max " << std::max(stddev,avg) << ", thresh " << thresh << std::endl;
                thresh = std::max(thresh, 500.0);
                //std::cout << "errorEst, level " << level << ": thresh cut " << thresh << std::endl;
                avg = thresh;
            }
#endif
            else
            {
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "Dont know the average of this variable "
                              << err_list[j].name() << '\n';
                avg = 0;
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Array<int> itags;

            for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
            {
		// FABs
		FArrayBox&  datfab  = (*mf)[mfi];
		TagBox&     tagfab  = tags[mfi];

		// Box in physical space
		int         idx     = mfi.index();
		RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());

		// tile box
		const Box&  tilebx  = mfi.tilebox();

		//fab box
		const Box&  datbox  = datfab.box();

		// We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
		// So we are going to get a temporary integer array.
		tagfab.get_itags(itags, tilebx);

		// data pointer and index space
		int*        tptr    = itags.dataPtr();
		const int*  tlo     = tilebx.loVect();
		const int*  thi     = tilebx.hiVect();
		//
		const int*  lo      = tlo;
		const int*  hi      = thi;
		//
		const Real* xlo     = gridloc.lo();
		//
		Real*       dat     = datfab.dataPtr();
		const int*  dlo     = datbox.loVect();
		const int*  dhi     = datbox.hiVect();
		const int   ncomp   = datfab.nComp();

                if (err_list[j].errType() == ErrorRec::Standard)
                {
                    err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, xlo, prob_lo, &time, &level);
                }
                else if (err_list[j].errType() == ErrorRec::UseAverage)
                {
                   err_list[j].errFunc2()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, &level, &avg);
                }

                //
                // Don't forget to set the tags in the TagBox.
                //
                if (allow_untagging == 1)
                {
                    tagfab.tags_and_untags(itags, tilebx);
                }
                else
                {
                   tagfab.tags(itags, tilebx);
                }
            }
        }

        delete mf;
    }
}

MultiFab*
Nyx::derive (const std::string& name,
             Real               time,
             int                ngrow)
{
    BL_PROFILE("Nyx::derive()");
    if (name == "Rank")
    {
        MultiFab* derive_dat = new MultiFab(grids, 1, 0);
        for (MFIter mfi(*derive_dat); mfi.isValid(); ++mfi)
        {
           (*derive_dat)[mfi].setVal(ParallelDescriptor::MyProc());
        }
        return derive_dat;
    } else {
        return particle_derive(name, time, ngrow);
    }
}

void
Nyx::derive (const std::string& name,
             Real               time,
             MultiFab&          mf,
             int                dcomp)
{
    BL_PROFILE("Nyx::derive(mf)");
    AmrLevel::derive(name, time, mf, dcomp);
}

void
Nyx::network_init ()
{
    BL_FORT_PROC_CALL(NYX_NETWORK_INIT, nyx_network_init)();
}

#ifndef NO_HYDRO
void
Nyx::reset_internal_energy (MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::reset_internal_energy()");
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    const Real  cur_time = state[State_Type].curTime();
    Real        a        = get_comoving_a(cur_time);
    const Real* dx       = geom.CellSize();
    const Real  vol      = D_TERM(dx[0],*dx[1],*dx[2]);

    Real sum_energy_added = 0;
    Real sum_energy_total = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_energy_added,sum_energy_total)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Real s  = 0;
        Real se = 0;
        BL_FORT_PROC_CALL(RESET_INTERNAL_E, reset_internal_e)
            (BL_TO_FORTRAN(S_new[mfi]), BL_TO_FORTRAN(D_new[mfi]),
             bx.loVect(), bx.hiVect(),
             &print_fortran_warnings, &a, &s, &se);
        sum_energy_added += s;
        sum_energy_total += se;
    }

    if (verbose > 1)
    {
        Real sums[2] = {sum_energy_added,sum_energy_total};

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sums,2,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            sum_energy_added = vol*sums[0];
            sum_energy_total = vol*sums[1];

            if (sum_energy_added > (1.e-12)*sum_energy_total)
            {
                std::cout << "Adding to (rho E) "
                          << sum_energy_added
                          << " out of total (rho E) "
                          << sum_energy_total << '\n';
            }
        }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_new_temp ()
{
    BL_PROFILE("Nyx::compute_new_temp()");
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real cur_time   = state[State_Type].curTime();

    reset_internal_energy(S_new,D_new);

    Real a = get_comoving_a(cur_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        BL_FORT_PROC_CALL(COMPUTE_TEMP, compute_temp)
            (bx.loVect(), bx.hiVect(),
            BL_TO_FORTRAN(S_new[mfi]),
            BL_TO_FORTRAN(D_new[mfi]), &a,
             &print_fortran_warnings);
    }

    // Compute the maximum temperature
    Real max_temp = D_new.norm0(Temp_comp);

    int imax = -1;
    int jmax = -1;
    int kmax = -1;

    Real den_maxt;

    // Find the cell which has the maximum temp -- but only if not the first
    // time step because in the first time step too many points have the same
    // value.
    Real prev_time   = state[State_Type].prevTime();
    if (prev_time > 0.0 && verbose > 0)
    {
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            BL_FORT_PROC_CALL(COMPUTE_MAX_TEMP_LOC, compute_max_temp_loc)
                (bx.loVect(), bx.hiVect(),
                BL_TO_FORTRAN(S_new[mfi]),
                BL_TO_FORTRAN(D_new[mfi]),
                &max_temp,&den_maxt,&imax,&jmax,&kmax);
        }

        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            if (imax > -1 && jmax > -1 && kmax > -1)
            {
              std::cout << "Maximum temp. at level " << level << " is " << max_temp
                        << " at density " << den_maxt
                        << " at (i,j,k) " << imax << " " << jmax << " " << kmax << std::endl;
            }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_rho_temp (Real& rho_T_avg, Real& T_avg, Real& T_meanrho)
{
    BL_PROFILE("Nyx::compute_rho_temp()");
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real rho_T_sum=0.0,   T_sum=0.0, T_meanrho_sum=0.0;
    Real   rho_sum=0.0, vol_sum=0.0,    vol_mn_sum=0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:rho_T_sum, rho_sum, T_sum, T_meanrho_sum, vol_sum, vol_mn_sum)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        BL_FORT_PROC_CALL(COMPUTE_RHO_TEMP, compute_rho_temp)
            (bx.loVect(), bx.hiVect(), geom.CellSize(),
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(D_new[mfi]), &average_gas_density,
             &rho_T_sum, &T_sum, &T_meanrho_sum, &rho_sum, &vol_sum, &vol_mn_sum);
    }
    Real sums[6] = {rho_T_sum, rho_sum, T_sum, T_meanrho_sum, vol_sum, vol_mn_sum};
    ParallelDescriptor::ReduceRealSum(sums,6);

    rho_T_avg = sums[0] / sums[1];  // density weighted T
        T_avg = sums[2] / sums[4];  // volume weighted T
    if (sums[5] > 0) {
       T_meanrho = sums[3] / sums[5];  // T at mean density
       T_meanrho = pow(10.0, T_meanrho);
    }
}
#endif


Real
Nyx::getCPUTime()
{

  int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores*omp_get_max_threads();
#endif

  Real T = numCores*(ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}

void
Nyx::AddProcsToComp(Amr *aptr, int nSidecarProcs, int prevSidecarProcs,
                    int ioProcNumSCS, int ioProcNumAll, int scsMyId,
                    MPI_Comm scsComm)
{
      BL_PROFILE("Nyx::AddProcsToComp()");

      forceParticleRedist = true;

      AmrLevel::AddProcsToComp(aptr, nSidecarProcs, prevSidecarProcs,
                               ioProcNumSCS, ioProcNumAll, scsMyId,
                               scsComm);


      // ---- pack up the bools
      Array<int> allBools;  // ---- just use ints here
      if(scsMyId == ioProcNumSCS) {
        allBools.push_back(dump_old);
        allBools.push_back(do_dm_particles);
        allBools.push_back(particle_initrandom_serialize);
        allBools.push_back(FillPatchedOldState_ok);
      }

      BoxLib::BroadcastArray(allBools, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the bools
      if(scsMyId != ioProcNumSCS) {
        int count(0);
        dump_old = allBools[count++];
        do_dm_particles = allBools[count++];
        particle_initrandom_serialize = allBools[count++];
        FillPatchedOldState_ok = allBools[count++];
        BL_ASSERT(count == allBools.size());
      }


      // ---- pack up the ints
      Array<int> allInts;

      if(scsMyId == ioProcNumSCS) {
        allInts.push_back(write_parameters_in_plotfile);
        allInts.push_back(print_fortran_warnings);
        allInts.push_back(particle_verbose);
        allInts.push_back(write_particles_in_plotfile);
        allInts.push_back(write_particle_density_at_init);
        allInts.push_back(write_coarsened_particles);
        allInts.push_back(NUM_STATE);
        allInts.push_back(Density);
        allInts.push_back(Xmom);
        allInts.push_back(Ymom);
        allInts.push_back(Zmom);
        allInts.push_back(Eden);
        allInts.push_back(Eint);
        allInts.push_back(Temp_comp);
        allInts.push_back(Ne_comp);
        allInts.push_back(FirstSpec);
        allInts.push_back(FirstAux);
        allInts.push_back(FirstAdv);
        allInts.push_back(NumSpec);
        allInts.push_back(NumAux);
        allInts.push_back(NumAdv);
        allInts.push_back(strict_subcycling);
        allInts.push_back(init_with_sph_particles);
        allInts.push_back(verbose);
        allInts.push_back(show_timings);
        allInts.push_back(do_reflux);
        allInts.push_back(NUM_GROW);
        allInts.push_back(nsteps_from_plotfile);
        allInts.push_back(allow_untagging);
        allInts.push_back(use_const_species);
        allInts.push_back(normalize_species);
        allInts.push_back(do_special_tagging);
        allInts.push_back(ppm_type);
        allInts.push_back(ppm_reference);
        allInts.push_back(ppm_flatten_before_integrals);
        allInts.push_back(use_colglaz);
        allInts.push_back(use_flattening);
        allInts.push_back(corner_coupling);
        allInts.push_back(version_2);
        allInts.push_back(use_exact_gravity);
        allInts.push_back(particle_initrandom_iseed);
        allInts.push_back(do_hydro);
        allInts.push_back(do_grav);
        allInts.push_back(add_ext_src);
        allInts.push_back(heat_cool_type);
        allInts.push_back(strang_split);
        allInts.push_back(reeber_int);
        allInts.push_back(gimlet_int);
        allInts.push_back(grav_n_grow);
        allInts.push_back(forceParticleRedist);
      }

      BoxLib::BroadcastArray(allInts, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the ints
      if(scsMyId != ioProcNumSCS) {
        int count(0);

        write_parameters_in_plotfile = allInts[count++];
        print_fortran_warnings = allInts[count++];
        particle_verbose = allInts[count++];
        write_particles_in_plotfile = allInts[count++];
        write_particle_density_at_init = allInts[count++];
        write_coarsened_particles = allInts[count++];
        NUM_STATE = allInts[count++];
        Density = allInts[count++];
        Xmom = allInts[count++];
        Ymom = allInts[count++];
        Zmom = allInts[count++];
        Eden = allInts[count++];
        Eint = allInts[count++];
        Temp_comp = allInts[count++];
        Ne_comp = allInts[count++];
        FirstSpec = allInts[count++];
        FirstAux = allInts[count++];
        FirstAdv = allInts[count++];
        NumSpec = allInts[count++];
        NumAux = allInts[count++];
        NumAdv = allInts[count++];
        strict_subcycling = allInts[count++];
        init_with_sph_particles = allInts[count++];
        verbose = allInts[count++];
        show_timings = allInts[count++];
        do_reflux = allInts[count++];
        NUM_GROW = allInts[count++];
        nsteps_from_plotfile = allInts[count++];
        allow_untagging = allInts[count++];
        use_const_species = allInts[count++];
        normalize_species = allInts[count++];
        do_special_tagging = allInts[count++];
        ppm_type = allInts[count++];
        ppm_reference = allInts[count++];
        ppm_flatten_before_integrals = allInts[count++];
        use_colglaz = allInts[count++];
        use_flattening = allInts[count++];
        corner_coupling = allInts[count++];
        version_2 = allInts[count++];
        use_exact_gravity = allInts[count++];
        particle_initrandom_iseed = allInts[count++];
        do_hydro = allInts[count++];
        do_grav = allInts[count++];
        add_ext_src = allInts[count++];
        heat_cool_type = allInts[count++];
        strang_split = allInts[count++];
        reeber_int = allInts[count++];
        gimlet_int = allInts[count++];
        grav_n_grow = allInts[count++];
        forceParticleRedist = allInts[count++];

        BL_ASSERT(count == allInts.size());
      }


      // ---- longs
      ParallelDescriptor::Bcast(&particle_initrandom_count, 1, ioProcNumAll, scsComm);


      // ---- pack up the Reals
      Array<Real> allReals;
      if(scsMyId == ioProcNumSCS) {
        allReals.push_back(initial_z);
        allReals.push_back(final_a);
        allReals.push_back(final_z);
        allReals.push_back(relative_max_change_a);
        allReals.push_back(old_a_time);
        allReals.push_back(new_a_time);
        allReals.push_back(old_a);
        allReals.push_back(new_a);
        allReals.push_back(particle_cfl);
        allReals.push_back(cfl);
        allReals.push_back(init_shrink);
        allReals.push_back(change_max);
        allReals.push_back(small_dens);
        allReals.push_back(small_temp);
        allReals.push_back(gamma);
        allReals.push_back( h_species);
        allReals.push_back(he_species);
        allReals.push_back(particle_initrandom_mass);
        allReals.push_back(average_gas_density);
        allReals.push_back(average_dm_density);
        allReals.push_back(average_neutr_density);
        allReals.push_back(average_total_density);
#ifdef NEUTRINO_PARTICLES
        allReals.push_back(neutrino_cfl);
#endif
      }

      BoxLib::BroadcastArray(allReals, scsMyId, ioProcNumAll, scsComm);
      BoxLib::BroadcastArray(plot_z_values, scsMyId, ioProcNumAll, scsComm);
      BoxLib::BroadcastArray(analysis_z_values, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the Reals
      if(scsMyId != ioProcNumSCS) {
        int count(0);
        initial_z = allReals[count++];
        final_a = allReals[count++];
        final_z = allReals[count++];
        relative_max_change_a = allReals[count++];
        old_a_time = allReals[count++];
        new_a_time = allReals[count++];
        old_a = allReals[count++];
        new_a = allReals[count++];
        particle_cfl = allReals[count++];
        cfl = allReals[count++];
        init_shrink = allReals[count++];
        change_max = allReals[count++];
        small_dens = allReals[count++];
        small_temp = allReals[count++];
        gamma = allReals[count++];
         h_species = allReals[count++];
        he_species = allReals[count++];
        particle_initrandom_mass = allReals[count++];
        average_gas_density = allReals[count++];
        average_dm_density = allReals[count++];
        average_neutr_density = allReals[count++];
        average_total_density = allReals[count++];
#ifdef NEUTRINO_PARTICLES
        neutrino_cfl = allReals[count++];
#endif
        BL_ASSERT(count == allReals.size());
      }


      // ---- pack up the strings
      Array<std::string> allStrings;
      Array<char> serialStrings;
      if(scsMyId == ioProcNumSCS) {
        allStrings.push_back(particle_plotfile_format);
        allStrings.push_back(particle_init_type);
        allStrings.push_back(particle_move_type);

        serialStrings = BoxLib::SerializeStringArray(allStrings);
      }

      BoxLib::BroadcastArray(serialStrings, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the strings
      if(scsMyId != ioProcNumSCS) {
        int count(0);
        allStrings = BoxLib::UnSerializeStringArray(serialStrings);

        particle_plotfile_format = allStrings[count++];
        particle_init_type = allStrings[count++];
        particle_move_type = allStrings[count++];
      }


      // ---- maps
      std::cout << "_in AddProcsToComp:  fix maps." << std::endl;
      //std::map<std::string,MultiFab*> auxDiag;
      //static std::map<std::string,Array<std::string> > auxDiag_names;


      // ---- pack up the IntVects
      Array<int> allIntVects;
      if(scsMyId == ioProcNumSCS) {
        for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(Nrep[i]); }

        BL_ASSERT(allIntVects.size() == BL_SPACEDIM);
      }

      BoxLib::BroadcastArray(allIntVects, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the IntVects
      if(scsMyId != ioProcNumSCS) {
          int count(0);

          BL_ASSERT(allIntVects.size() == BL_SPACEDIM);
          for(int i(0); i < BL_SPACEDIM; ++i)    { Nrep[i] = allIntVects[count++]; }

          BL_ASSERT(allIntVects.size() == BL_SPACEDIM);
      }



      // ---- BCRec
      Array<int> bcrLo(BL_SPACEDIM), bcrHi(BL_SPACEDIM);
      if(scsMyId == ioProcNumSCS) {
        for(int i(0); i < bcrLo.size(); ++i) { bcrLo[i] = phys_bc.lo(i); }
        for(int i(0); i < bcrHi.size(); ++i) { bcrHi[i] = phys_bc.hi(i); }
      }
      ParallelDescriptor::Bcast(bcrLo.dataPtr(), bcrLo.size(), ioProcNumSCS, scsComm);
      ParallelDescriptor::Bcast(bcrHi.dataPtr(), bcrHi.size(), ioProcNumSCS, scsComm);
      if(scsMyId != ioProcNumSCS) {
        for(int i(0); i < bcrLo.size(); ++i) { phys_bc.setLo(i, bcrLo[i]); }
        for(int i(0); i < bcrHi.size(); ++i) { phys_bc.setHi(i, bcrHi[i]); }
      }




      // ---- ErrorList
      if(scsMyId != ioProcNumSCS) {
        InitErrorList();
      }

      // ---- DeriveList
      if(scsMyId != ioProcNumSCS) {
        InitDeriveList();
      }



      int isAllocated(0);
#ifndef NO_HYDRO
      // ---- FluxRegister
      if(scsMyId == ioProcNumSCS) {
        if(flux_reg == 0) {
          isAllocated = 0;
        } else {
          isAllocated = 1;
        }
      }
      ParallelDescriptor::Bcast(&isAllocated, 1, ioProcNumSCS, scsComm);
      if(isAllocated == 1) {
        if(scsMyId != ioProcNumSCS) {
          BL_ASSERT(flux_reg == 0);
          flux_reg = new FluxRegister;
        }
        flux_reg->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }
#endif


      // ---- fine_mask
      isAllocated = 0;
      if(scsMyId == ioProcNumSCS) {
        if(fine_mask == 0) {
          isAllocated = 0;
        } else {
          isAllocated = 1;
        }
      }
      ParallelDescriptor::Bcast(&isAllocated, 1, ioProcNumSCS, scsComm);
      if(isAllocated == 1) {
        if(scsMyId != ioProcNumSCS) {
          BL_ASSERT(fine_mask == 0);
          std::cout << "**** check fine_mask." << std::endl;
          fine_mask = build_fine_mask();
        }
        fine_mask->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }


#ifdef GRAVITY
      // ---- gravity
      if(do_grav) {
        if(gravity != 0) {
          gravity->AddProcsToComp(parent, level, this, ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
	}
   }
#endif


       NyxParticlesAddProcsToComp(parent, nSidecarProcs, prevSidecarProcs, ioProcNumSCS,
                                  ioProcNumAll, scsMyId, scsComm);

}


void
Nyx::InitErrorList() {
    //err_list.clear(true);
    //err_list.add("FULLSTATE",1,ErrorRec::Special,FORT_DENERROR);
}


//static Box the_same_box (const Box& b) { return b; }


void
Nyx::InitDeriveList() {

/*
    derive_lst.clear();
    //
    // Set number of state variables
    //
    int counter = (int)Density + 1;
    int Xmom = counter++;
    int Ymom = counter++;
#if(BL_SPACEDIM==3)
    int Zmom = counter++;
#endif
    int Eden = counter++;
#if(NADV>0)
    int Tracer = counter++;
#if(NADV>1)
    BoxLib::Error("Prob::variableSetUp: Only one Advected quantity allowed");
#endif
#endif
#ifdef BL_USE_CHEM
    NumSpec = getChemDriver().numSpecies();
    if (NumSpec > 0)
    {
        FirstSpec = counter++;
        counter += NumSpec - 2;
        LastSpec = counter++;
    }

    const int tmp = (int)Density;
    FORT_SETCOMPS (&tmp, &Xmom,&Ymom,&Zmom,&Eden, &FirstSpec,&Tracer);

    const Array<std::string>& names = getChemDriver().speciesNames();
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Chemical species interpreted:\n { ";
        for (int i = 0; i < names.size(); i++)
            std::cout << names[i] << ' ' << ' ';
        std::cout << '}' << '\n' << '\n';
    }
#endif

    NUM_STATE = counter;

    //
    // DEFINE DERIVED QUANTITIES
    //
    // log of Density
    //
    derive_lst.add("log_den",IndexType::TheCellType(),1,FORT_DERLOGS,the_same_box);
    derive_lst.addComponent("log_den",desc_lst,State_Type,Density,1);
    //
    // Pressure
    //
#ifdef BL_USE_CHEM
    const int nWorkPres = 4 + NumSpec;
#else
    const int nWorkPres = 4;
#endif
    derive_lst.add("pressure",IndexType::TheCellType(),nWorkPres,
                   FORT_DERPRES,the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

#ifdef  PRISCILLA
    // mach
    derive_lst.add("mach",IndexType::TheCellType(),4, FORT_DERMACH,the_same_box);
    derive_lst.addComponent("mach",desc_lst,State_Type,Density,NUM_STATE);
#endif

    //
    // Xvel
    //
    derive_lst.add("xvel",IndexType::TheCellType(),1,FORT_DERVEL,the_same_box);
    derive_lst.addComponent("xvel",desc_lst,State_Type,Density,2);
    //
    // Yvel
    //
    derive_lst.add("yvel",IndexType::TheCellType(),1,FORT_DERVEL,the_same_box);
    derive_lst.addComponent("yvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("yvel",desc_lst,State_Type,Ymom,1);

#if(BL_SPACEDIM==3)
    //
    // Zvel
    //
    derive_lst.add("zvel",IndexType::TheCellType(),1,FORT_DERVEL,the_same_box);
    derive_lst.addComponent("zvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("zvel",desc_lst,State_Type,Zmom,1);
#endif
    //
    // log of Eden
    //
    derive_lst.add("log_eden",IndexType::TheCellType(),1,FORT_DERLOGS,the_same_box);
    derive_lst.addComponent("log_eden",desc_lst,State_Type,Eden,1);
    //
    // A derived quantity equal to all the state variables.
    //
    derive_lst.add("FULLSTATE",IndexType::TheCellType(),NUM_STATE,FORT_DERCOPY,the_same_box);
    derive_lst.addComponent("FULLSTATE",desc_lst,State_Type,Density,NUM_STATE);
*/
}
