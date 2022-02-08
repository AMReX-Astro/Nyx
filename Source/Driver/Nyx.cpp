#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <AMReX_CONSTANTS.H>
#include <Nyx.H>
#include <Gravity.H>
#include <constants_cosmo.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

#ifdef BL_USE_MPI
#include <MemInfo.H>
#endif

#ifndef NO_HYDRO
#include <atomic_rates_data.H>
#include <Derive.H>
#endif

#include <Forcing.H>

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#endif

#ifdef AGN
#include <agn_F.H>
#endif

using namespace amrex;

extern "C" {
  int get_comp_urho();
  int get_comp_temp();
  int get_comp_e_int();
}

const int NyxHaloFinderSignal = 42;
const int GimletSignal = 55;

static int sum_interval = -1;
static int  max_temp_dt = -1;

static Real fixed_dt    = -1.0;
static Real initial_dt  = -1.0;
static Real dt_cutoff   =  0;

Real Nyx::old_a      = -1.0;
Real Nyx::new_a      = -1.0;
Real Nyx::old_a_time = -1.0;
Real Nyx::new_a_time = -1.0;

Vector<Real> Nyx::plot_z_values;
Vector<Real> Nyx::analysis_z_values;
int Nyx::insitu_start = 0;
int Nyx::insitu_int = 0;

int Nyx::load_balance_int = -1;
amrex::Real Nyx::load_balance_start_z = 15;
int Nyx::load_balance_wgt_strategy = 0;
int Nyx::load_balance_wgt_nmax = -1;
int Nyx::load_balance_strategy = DistributionMapping::SFC;

bool Nyx::dump_old = false;
int Nyx::verbose      = 0;

Real Nyx::cfl = 0.8;
Real Nyx::init_shrink = 1.0;
Real Nyx::change_max  = 1.1;

BCRec Nyx::phys_bc;
int Nyx::do_reflux = 1;
int Nyx::NUM_STATE = -1;

int Nyx::nsteps_from_plotfile = -1;

ErrorList Nyx::err_list;
Vector<AMRErrorTag> Nyx::errtags;

int Nyx:: Zhi_comp = -1;

int Nyx::NumSpec  = 0;

Real Nyx::small_dens = -1.e200;
Real Nyx::small_temp = -1.e200;
Real Nyx::small_pres = -1.e200;
Real Nyx::large_temp =  1.e9;
Real Nyx::gamma      =  5.0/3.0;

// This is no longer an optional input;
//   we use the value hard-wired here
Real Nyx::small      =  1.e-6;

Real Nyx::comoving_OmB;
Real Nyx::comoving_OmM;
Real Nyx::comoving_OmR = 0.0;
Real Nyx::comoving_h;
int  Nyx::comoving_type = 1;

int Nyx::do_hydro = 0;
int Nyx::add_ext_src = 0;
int Nyx::heat_cool_type = 0;
int Nyx::use_sundials_constraint = 0;
int Nyx::use_sundials_fused = 0;
int Nyx::use_typical_steps = 0;
#ifndef AMREX_USE_GPU
int Nyx::sundials_alloc_type = 0;
int Nyx::sundials_atomic_reductions = -1; // CUDA and HIP only
#else
#ifdef AMREX_USE_CUDA
int Nyx::sundials_atomic_reductions = 1;
#ifndef _OPENMP
int Nyx::sundials_alloc_type = 0; //consider changing to 5
#else
int Nyx::sundials_alloc_type = 5;
#endif
#endif
#ifdef AMREX_USE_HIP
int Nyx::sundials_atomic_reductions = 0;
int Nyx::sundials_alloc_type = 5;
#endif
#ifdef AMREX_USE_DPCPP
int Nyx::sundials_atomic_reductions = -1; // CUDA and HIP only
int Nyx::sundials_alloc_type = 5;
#endif
#endif

Real Nyx::sundials_reltol = 1e-4;
Real Nyx::sundials_abstol = 1e-4;

int Nyx::minimize_memory = 0;
int Nyx::shrink_to_fit = 0;

bool Nyx::sundials_use_tiling = true;

#ifndef AMREX_USE_GPU
IntVect      Nyx::hydro_tile_size(1024,16,16);
#else
IntVect      Nyx::hydro_tile_size(1048576,1048576,1048576);
#endif

#ifndef AMREX_USE_GPU
IntVect      Nyx::sundials_tile_size(1024,16,16);
#else
#ifdef AMREX_USE_OMP
IntVect      Nyx::sundials_tile_size(1024,16,16);
#else
IntVect      Nyx::sundials_tile_size(1048576,1048576,1048576);
#endif
#endif

int Nyx::strang_split = 1;
int Nyx::strang_grown_box = 1;
#ifdef SDC
int Nyx::sdc_split    = 0;
int Nyx::strang_restart_from_sdc    = 0;
#endif

Real Nyx::average_gas_density = 0.;
Real Nyx::average_dm_density = 0.;
Real Nyx::average_neutr_density = 0.;
Real Nyx::average_total_density = 0.;

long int Nyx::old_max_sundials_steps = 3;
long int Nyx::new_max_sundials_steps = 3;

int         Nyx::inhomo_reion = 0;
std::string Nyx::inhomo_zhi_file = "";
int         Nyx::inhomo_grid = -1;

Gravity* Nyx::gravity  =  0;
int Nyx::do_grav       = -1;

#ifndef NO_HYDRO
StochasticForcing* Nyx::forcing = 0;
#endif
int Nyx::do_forcing =  0;

int Nyx::nghost_state       = 1;
Real Nyx::tagging_base       = 8.0;
int Nyx::reuse_mlpoisson     = 0;
int Nyx::ppm_type           = 1;

// Options are "floor" or "conservative"
std::string Nyx::enforce_min_density_type = "floor";

Real Nyx:: h_species        = 0.76;
Real Nyx::he_species        = 0.24;

#ifdef REEBER
Real Nyx::mass_halo_min     = 1.e10;
Real Nyx::mass_seed         = 1.e5;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// The default for how many digits to use for each column in the runlog
// We have a second parameter (rlp_terse) for regression-testing those quantities 
//    which show more variation than others
int Nyx::runlog_precision = 6;
int Nyx::runlog_precision_terse = 6;

int Nyx::write_parameters_in_plotfile = true;
int Nyx::write_skip_prepost = 0;
int Nyx::write_hdf5 = 0;

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

#ifndef NO_HYDRO
#ifndef HEATCOOL
AtomicRates* atomic_rates_glob;
#endif
#endif

// this will be reset upon restart
Real         Nyx::previousCPUTimeUsed = 0.0;

Real         Nyx::startCPUTime = 0.0;

int reeber_int(0);
int gimlet_int(0);

// Note: Nyx::variableSetUp is in Nyx_setup.cpp
void
Nyx::variable_cleanup ()
{
if (do_grav)
{
    if (gravity != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting gravity in variable_cleanup...\n";
        delete gravity;
        gravity = 0;
    }
}

#ifndef NO_HYDRO
    if (forcing != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting forcing in variable_cleanup...\n";
        delete forcing;
        forcing = 0;
    }
#endif

    desc_lst.clear();
}

void
Nyx::read_params ()
{
    BL_PROFILE("Nyx::read_params()");

    ParmParse pp_nyx("nyx");

    pp_nyx.query("v", verbose);
    pp_nyx.get("init_shrink", init_shrink);
    pp_nyx.get("cfl", cfl);
    pp_nyx.query("change_max", change_max);
    pp_nyx.query("fixed_dt", fixed_dt);
    pp_nyx.query("initial_dt", initial_dt);
    pp_nyx.query("max_temp_dt", max_temp_dt);
    pp_nyx.query("sum_interval", sum_interval);
    pp_nyx.get("dt_cutoff", dt_cutoff);

    // Get boundary conditions
    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp_nyx.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp_nyx.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        phys_bc.setLo(i, lo_bc[i]);
        phys_bc.setHi(i, hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (DefaultGeometry().isAnyPeriodic() || (!do_dm_particles && !do_hydro))
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (DefaultGeometry().isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior" << std::endl;
                    amrex::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior" << std::endl;
                    amrex::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                amrex::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                amrex::Error();
            }
        }
    }

    read_comoving_params();

#ifdef AMREX_PARTICLES
    read_particle_params();
#endif

    if (do_grav)
       pp_nyx.get("do_grav", do_grav);

#ifndef NO_HYDRO
    read_hydro_params();
#else
    if (!do_grav)
        amrex::Error("Dont know what to do with both hydro and gravity off");
#endif
    if (do_dm_particles || do_hydro)
        read_init_params();

    pp_nyx.query("runlog_precision",runlog_precision);
    pp_nyx.query("runlog_precision_terse",runlog_precision_terse);

    pp_nyx.query("write_parameter_file",write_parameters_in_plotfile);
    if(pp_nyx.query("write_hdf5",write_hdf5))
        write_skip_prepost = write_hdf5;
    else
        pp_nyx.query("write_skip_prepost",write_skip_prepost);

#ifndef AMREX_USE_HDF5
    if(write_hdf5 == 1)
        amrex::Error("Must compile with USE_HDF5 = TRUE for write_hdf5 = 1");
#endif

#ifdef AMREX_USE_HDF5_ASYNC
    // Complete all previous async writes on amrex::Finalize()
    amrex::ExecOnFinalize(H5VLasync_waitall);
#endif

    if (pp_nyx.contains("plot_z_values"))
    {
      int num_z_values = pp_nyx.countval("plot_z_values");
      plot_z_values.resize(num_z_values);
      pp_nyx.queryarr("plot_z_values",plot_z_values,0,num_z_values);
    }

    if (pp_nyx.contains("analysis_z_values"))
    {
      int num_z_values = pp_nyx.countval("analysis_z_values");
      analysis_z_values.resize(num_z_values);
      pp_nyx.queryarr("analysis_z_values",analysis_z_values,0,num_z_values);
    }

    ParmParse pp_insitu("insitu");
    pp_insitu.query("int", insitu_int);
    pp_insitu.query("start", insitu_start);
    pp_insitu.query("reeber_int", reeber_int);

    pp_nyx.query("load_balance_int",          load_balance_int);
    pp_nyx.query("load_balance_start_z",          load_balance_start_z);
    pp_nyx.query("load_balance_wgt_strategy", load_balance_wgt_strategy);
    load_balance_wgt_nmax = amrex::ParallelDescriptor::NProcs();
    pp_nyx.query("load_balance_wgt_nmax",     load_balance_wgt_nmax);

    std::string theStrategy;

    if (pp_nyx.query("load_balance_strategy", theStrategy))
    {
        if (theStrategy == "SFC")
        {
            load_balance_strategy=DistributionMapping::Strategy::SFC;
        }
        else if (theStrategy == "KNAPSACK")
        {
            load_balance_strategy=DistributionMapping::Strategy::KNAPSACK;
        }
        else if (theStrategy == "ROUNDROBIN")
        {
            load_balance_strategy=DistributionMapping::Strategy::ROUNDROBIN;
        }
        else
        {
            std::string msg("Unknown strategy: ");
            msg += theStrategy;
            amrex::Warning(msg.c_str());
        }
    }

    pp_nyx.query("gimlet_int", gimlet_int);
#ifdef REEBER
    pp_nyx.query("mass_halo_min", mass_halo_min);
    pp_nyx.query("mass_seed", mass_seed);
#endif
}

#ifndef NO_HYDRO
void
Nyx::read_hydro_params ()
{
    ParmParse pp_nyx("nyx");

    pp_nyx.get("do_hydro", do_hydro);

    pp_nyx.query("small_dens", small_dens);
    pp_nyx.query("small_temp", small_temp);
    pp_nyx.query("small_pres", small_pres);
    pp_nyx.query("large_temp", large_temp);
    pp_nyx.query("gamma", gamma);

#ifndef NO_HYDRO
#ifdef HEATCOOL
    atomic_rates_glob = (AtomicRates*)The_Arena()->alloc(sizeof(AtomicRates));
#else
    atomic_rates_glob = NULL;
#endif
#endif

    pp_nyx.query("add_ext_src", add_ext_src);
    pp_nyx.query("strang_split", strang_split);
    pp_nyx.query("strang_grown_box", strang_grown_box);

#ifdef HEATCOOL
#ifdef SDC
    pp_nyx.query("sdc_split", sdc_split);
    pp_nyx.query("strang_restart_from_sdc", strang_restart_from_sdc);
    if (sdc_split == 1 && strang_split == 1)
        amrex::Error("Cant have strang_split == 1 and sdc_split == 1");
    if (sdc_split == 0 && strang_split == 0)
        amrex::Error("Cant have strang_split == 0 and sdc_split == 0");
    if (sdc_split != 1 && strang_split != 1)
        amrex::Error("Cant have strang_split != 1 and sdc_split != 1");
#else
    if (strang_split != 1)
        amrex::Error("Cant have strang_split != 1 with USE_SDC != TRUE");
#endif
#endif

    pp_nyx.query("do_forcing", do_forcing);
    if (do_forcing == 1 && add_ext_src == 0)
       amrex::Error("Nyx::must set add_ext_src to 1 if do_forcing = 1 ");

    pp_nyx.query("heat_cool_type", heat_cool_type);
    pp_nyx.query("inhomo_reion", inhomo_reion);

    if (inhomo_reion) {
        pp_nyx.get("inhomo_zhi_file", inhomo_zhi_file);
        pp_nyx.get("inhomo_grid", inhomo_grid);
    }

#ifdef HEATCOOL
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Integrating heating/cooling method with the following method: ";
      switch (heat_cool_type)
        {
        case 11:
          std::cout << "Vectorized CVODE";
          break;
      }
      std::cout << std::endl;
    }

#else
    if (inhomo_reion > 0)
       amrex::Error("Nyx::you set inhomo_reion > 0 but forgot to set USE_HEATCOOL = TRUE");
#endif  // HEATCOOL

    pp_nyx.query("use_sundials_constraint", use_sundials_constraint);
    pp_nyx.query("use_sundials_fused", use_sundials_fused);
    pp_nyx.query("nghost_state", nghost_state);
    pp_nyx.query("sundials_atomic_reductions", sundials_atomic_reductions);
    pp_nyx.query("sundials_alloc_type", sundials_alloc_type);
    pp_nyx.query("sundials_reltol", sundials_reltol);
    pp_nyx.query("sundials_abstol", sundials_abstol);
    pp_nyx.query("minimize_memory", minimize_memory);
    pp_nyx.query("shrink_to_fit", shrink_to_fit);
    pp_nyx.query("use_typical_steps", use_typical_steps);
    pp_nyx.query("tagging_base", tagging_base);
    pp_nyx.query("reuse_mlpoisson", reuse_mlpoisson);
    pp_nyx.query("ppm_type", ppm_type);
    pp_nyx.query("enforce_min_density_type", enforce_min_density_type);

    pp_nyx.query("sundials_use_tiling", sundials_use_tiling);

    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp_nyx.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) {
          hydro_tile_size[i] = tilesize[i];
        }
    }
    else
    {
        amrex::Print()<<"Nyx::hydro_tile_size unset, using fabarray.mfiter_tile_size default: "<<
            FabArrayBase::mfiter_tile_size<<
            "\nSuggested default for currently compiled CPU / GPU: nyx.hydro_tile_size="<<
            hydro_tile_size<<std::endl;
        hydro_tile_size = FabArrayBase::mfiter_tile_size;
    }

    if (pp_nyx.queryarr("sundials_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) {
          sundials_tile_size[i] = tilesize[i];
        }
    }
    else
    {
        amrex::Print()<<"Nyx::sundials_tile_size unset, using fabarray.mfiter_tile_size default: "<<
            FabArrayBase::mfiter_tile_size<<
            "\nSuggested default for currently compiled CPU / GPU: nyx.sundials_tile_size="<<
            sundials_tile_size<<std::endl;
        sundials_tile_size = FabArrayBase::mfiter_tile_size;
    }

    if (use_typical_steps != 0 && strang_grown_box == 0)
    {
          amrex::Error("Nyx::use_typical_steps must be 0 with strang_grown_box = 0");
    }

    if (do_hydro == 1)
    {
#ifdef CONST_SPECIES
        pp_nyx.get("h_species" ,  h_species);
        pp_nyx.get("he_species", he_species);

        amrex::Print() << "Nyx::setting species concentrations to "
                       << h_species << " and " << he_species
                       << " in the hydro and in the EOS " << std::endl;
#endif

        // Make sure ppm_type is set correctly.
        if (ppm_type != 0 && ppm_type != 1)
           amrex::Error("Nyx::ppm_type must be 0 or 1");
    }

    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;
}
#endif //NO_HYDRO

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
          const DistributionMapping& dm,
          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    BL_PROFILE("Nyx::Nyx(Amr)");

    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(lev));

    build_metrics();
    fine_mask = 0;

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
    }
#endif

    // Initialize to zero here in case we run with do_grav = false.
    MultiFab& new_grav_mf = get_new_data(Gravity_Type);
    new_grav_mf.setVal(0);

    MultiFab& new_phi_grav = get_new_data(PhiGrav_Type);
    new_phi_grav.setVal(0);

    if (do_grav)
    {
        // gravity is a static object, only alloc if not already there
        if (gravity == 0) {
          gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, 0);
        }

        gravity->install_level(level, this);
   }

#ifndef NO_HYDRO
    if (do_forcing)
    {
        // forcing is a static object, only alloc if not already there
        if (forcing == 0)
           forcing = new StochasticForcing();

        const Real* prob_lo = geom.ProbLo();
        const Real* prob_hi = geom.ProbHi();

        forcing->init(AMREX_SPACEDIM, prob_lo, prob_hi);
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
      // get elapsed CPU time
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
            flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
    }
#endif

    if (do_grav && level == 0)
    {
        BL_ASSERT(gravity == 0);
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, 0);
    }

#ifndef NO_HYDRO
    if (do_forcing)
    {
        // forcing is a static object, only alloc if not already there
        if (forcing == 0)
           forcing = new StochasticForcing();

        const Real* prob_lo = geom.ProbLo();
        const Real* prob_hi = geom.ProbHi();

        forcing->init(AMREX_SPACEDIM, prob_lo, prob_hi);
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
    if (verbose && ParallelDescriptor::IOProcessor()) {
       std::cout << "Setting the current time in the state data to "
                 << parent->cumTime() << std::endl;
    }
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
Nyx::init (AmrLevel& old)
{
    BL_PROFILE("Nyx::init(old)");

    MultiFab::RegionTag amrInit_tag("Init_" + std::to_string(level));
    Nyx* old_level = (Nyx*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);

#ifndef NO_HYDRO
    Real cur_time  = old_level->state[State_Type].curTime();
    Real prev_time = old_level->state[State_Type].prevTime();
#else
    Real cur_time  = old_level->state[PhiGrav_Type].curTime();
    Real prev_time = old_level->state[PhiGrav_Type].prevTime();
#endif

    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);

        FillPatch(old, S_new, 0, cur_time,   State_Type, 0, NUM_STATE);
        FillPatch(old, D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());

        MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
        reset_e_src.setVal(0.0);
        reset_internal_energy_interp(S_new,D_new,reset_e_src);

    }
#endif

    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    if (do_grav)
    {
        FillPatch(old, Phi_new, 0, cur_time, PhiGrav_Type, 0, 1);
    } else {
        // We need to initialize it otherwise we might write out NaNs in the checkpoint 
        Phi_new.setVal(0.);
    }

#ifdef SDC
    MultiFab& IR_new = get_new_data(SDC_IR_Type);
    FillPatch(old, IR_new, 0, cur_time, SDC_IR_Type, 0, 1);
#endif

    amrex::Gpu::Device::streamSynchronize();

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

#ifndef NO_HYDRO
    Real cur_time  = get_level(level-1).state[State_Type].curTime();
    Real prev_time = get_level(level-1).state[State_Type].prevTime();
#else
    Real cur_time  = get_level(level-1).state[PhiGrav_Type].curTime();
    Real prev_time = get_level(level-1).state[PhiGrav_Type].prevTime();
#endif

    Real dt_old    = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time, dt_old, dt);

#ifndef NO_HYDRO

    if(do_hydro)
    {

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    FillCoarsePatch(S_new, 0, cur_time,   State_Type, 0, S_new.nComp());
    FillCoarsePatch(D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());

    Real rho_E = 0;
    Real rho_e = 0;

    for (int lev = 0; lev <= 0; lev++)
      //      for (int lev = 0; lev <= finest_level; lev++)
    {
        Nyx& nyx_lev = get_level(lev);

        rho_E += nyx_lev.vol_weight_sum("rho_E", prev_time, true);
        rho_e += nyx_lev.vol_weight_sum("rho_e", prev_time, true);
    }
    amrex::Print()<<"total RHO*E    "<<rho_E/(geom.ProbSize())<<std::endl;

    MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
    reset_e_src.setVal(0.0);
    reset_internal_energy_interp(S_new,D_new,reset_e_src);
    }
#endif

    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    if (do_grav)
    {
        FillCoarsePatch(Phi_new, 0, cur_time, PhiGrav_Type, 0, Phi_new.nComp());
    } else {
        // We need to initialize it otherwise we might write out NaNs in the checkpoint 
        Phi_new.setVal(0.);
    }

#ifdef SDC
    MultiFab& IR_new = get_new_data(SDC_IR_Type);
    FillCoarsePatch(IR_new, 0, cur_time, SDC_IR_Type, 0, IR_new.nComp());
#endif

    // We set dt to be large for this new level to avoid screwing up
    // computeNewDt.
    parent->setDtLevel(1.e100, level);
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
Nyx::est_time_step (Real /*dt_old*/)
{
    BL_PROFILE("Nyx::est_time_step()");
    if (fixed_dt > 0)
        return fixed_dt;

    // This is just a dummy value to start with
    Real est_dt = 1.0e+200;

#ifndef NO_HYDRO
    const MultiFab& stateMF = get_new_data(State_Type);
#endif

#ifndef NO_HYDRO
    Real cur_time = state[State_Type].curTime();
#else
    Real cur_time = state[PhiGrav_Type].curTime();
#endif

#ifndef NO_HYDRO
    if (do_hydro)
    {
        Real a = get_comoving_a(cur_time);
        const auto dx = geom.CellSizeArray();

        {
          Real dt = 1.e200;
          Real sound_speed_factor=sqrt(gamma*(gamma-1));

          Real local_small_dens  = small_dens;
          int  local_max_temp_dt = max_temp_dt;

          dt = amrex::ReduceMin(stateMF, 0,
              [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& u) -> Real
              {
                  const auto lo = amrex::lbound(bx);
                  const auto hi = amrex::ubound(bx);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
                  amrex::Real dt_gpu = std::numeric_limits<amrex::Real>::max();
#else
                  amrex::Real dt_gpu = 1.e37;
#endif

                  for         (int k = lo.z; k <= hi.z; ++k) {
                    for     (int j = lo.y; j <= hi.y; ++j) {
                      for (int i = lo.x; i <= hi.x; ++i) {
                        if(u(i,j,k,Density_comp)<=1.1*local_small_dens && local_max_temp_dt==1)
                          continue;

                            Real rhoInv = 1.0 / u(i,j,k,Density_comp);
                            Real ux     = u(i,j,k,Xmom_comp)*rhoInv;
                            Real uy     = u(i,j,k,Ymom_comp)*rhoInv;
                            Real uz     = u(i,j,k,Zmom_comp)*rhoInv;

                            // Use internal energy for calculating dt
                            Real e  = u(i,j,k,Eint_comp)*rhoInv;

                            Real c;
                            // Protect against negative e
#ifdef HEATCOOL
                            if (e > 0.0)
                              c=sound_speed_factor*std::sqrt(e);
#else
                            if (e > 0.0)
                              c=sound_speed_factor*std::sqrt(u(i,j,k,Density_comp)*e/u(i,j,k,Density_comp));
#endif
                            else
                              c = 0.0;

                            Real dt1 = dx[0]/(c + amrex::Math::abs(ux));
                            Real dt2 = dx[1]/(c + amrex::Math::abs(uy));
                            Real dt3 = dx[2]/(c + amrex::Math::abs(uz));
                            dt_gpu = amrex::min(dt_gpu,amrex::min(dt1,amrex::min(dt2,dt3)));
                       }
                    }
                  }
                  return dt_gpu;
              });

          est_dt = std::min(est_dt, dt);

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

    if (do_grav)
        particle_est_time_step(est_dt);

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
Nyx::computeNewDt (int                      finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&             n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&            dt_min,
                   Vector<Real>&            dt_level,
                   Real                     stop_time,
                   int                      post_regrid_flag)
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
                Vector<int> new_cycle(finest_level+1);
                for (i = 0; i <= finest_level; i++)
                    new_cycle[i] = n_cycle[i];
                // The max allowable dt
                Vector<Real> dt_max(finest_level+1);
                for (i = 0; i <= finest_level; i++)
                {
                    dt_max[i] = dt_min[i];
                }
                // find the maximum number of cycles allowed.
                Vector<int> cycle_max(finest_level+1);
                cycle_max[0] = 1;
                for (i = 1; i <= finest_level; i++)
                {
                    cycle_max[i] = parent->MaxRefRatio(i-1);
                }
                // estimate the amout of work to advance each level.
                Vector<Real> est_work(finest_level+1);
                for (i = 0; i <= finest_level; i++)
                {
                    est_work[i] = parent->getLevel(i).estimateWork();
                }

                // This value will be used only if the subcycling pattern is changed.

                dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle.dataPtr(), dt_max.dataPtr(), 
                                                        est_work.dataPtr(), cycle_max.dataPtr());

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
#ifndef NO_HYDRO
    Real cur_time = state[State_Type].curTime();
#else
    Real cur_time = state[PhiGrav_Type].curTime();
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
Nyx::computeInitialDt (int                      finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&             n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&            dt_level,
                       Real                     stop_time)
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
        Vector<int> new_cycle(finest_level+1);
        for (i = 0; i <= finest_level; i++)
            new_cycle[i] = n_cycle[i];
        Vector<Real> dt_max(finest_level+1);
        for (i = 0; i <= finest_level; i++)
        {
            dt_max[i] = get_level(i).initial_time_step();
        }
        // Find the maximum number of cycles allowed
        Vector<int> cycle_max(finest_level+1);
        cycle_max[0] = 1;
        for (i = 1; i <= finest_level; i++)
        {
            cycle_max[i] = parent->MaxRefRatio(i-1);
        }
        // estimate the amout of work to advance each level.
        Vector<Real> est_work(finest_level+1);
        for (i = 0; i <= finest_level; i++)
        {
            est_work[i] = parent->getLevel(i).estimateWork();
        }

        dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle.dataPtr(), dt_max.dataPtr(), 
                                                est_work.dataPtr(), cycle_max.dataPtr());

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
#ifndef NO_HYDRO
    Real cur_time = state[State_Type].curTime();
#else
    Real cur_time = state[PhiGrav_Type].curTime();
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
        amrex::Error("Should only call writePlotNow at level 0!");

    bool found_one = false;

    if (plot_z_values.size() > 0)
    {
#ifndef NO_HYDRO
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#else
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
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
        amrex::Error("Should only call doAnalysisNow at level 0!");

    bool found_one = false;

    if (analysis_z_values.size() > 0)
    {

#ifndef NO_HYDRO
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#else
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
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

    MultiFab::RegionTag amrPost_tag("Post_" + std::to_string(level));

    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    const int ncycle = parent->nCycle(level);
#ifdef AMREX_PARTICLES
    BL_PROFILE_VAR("Nyx::post_timestep()::remove_virt_ghost",rm);
    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration
    //
    if (iteration == ncycle)
        remove_ghost_particles();
    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(rm);
    BL_PROFILE_VAR("Nyx::post_timestep()::redist",redist);

    if(load_balance_int < 0 || !(nStep() % load_balance_int == 0 && level == 0 && (new_a >= 1.0/(load_balance_start_z + 1.0))))
    {

    //
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level.
    //
    if ((iteration < ncycle && level < finest_level) || level == 0)
    {
        for (int i = 0; i < theActiveParticles().size(); i++)
        {
             if(finest_level == 0)
                 theActiveParticles()[i]->RedistributeLocal(level,
                                                  theActiveParticles()[i]->finestLevel(),
                                                  iteration);
             else
                 theActiveParticles()[i]->Redistribute(level,
                                                  theActiveParticles()[i]->finestLevel(),
                                                  iteration);

             if(shrink_to_fit)
                 theActiveParticles()[i]->ShrinkToFit();
        }
    }

    }
    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(redist);
#endif
    BL_PROFILE_VAR("Nyx::post_timestep()::do_reflux",do_reflux);

#ifndef NO_HYDRO
    if (do_reflux && level < finest_level)
    {
        MultiFab& S_new_crse = get_new_data(State_Type);
        MultiFab drho_and_drhoU;
        if (do_grav)
        {
            // Define the update to rho and rhoU due to refluxing.
            drho_and_drhoU.define(grids, dmap, AMREX_SPACEDIM + 1, 0);
            MultiFab::Copy(drho_and_drhoU, S_new_crse, Density_comp, 0,
                           AMREX_SPACEDIM + 1, 0);
            drho_and_drhoU.mult(-1.0);
        }

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
#ifndef CONST_SPECIES
        enforce_nonnegative_species(S_new_crse);
#endif

        if (do_grav && gravity->get_no_sync() == 0)
        {
            MultiFab::Add(drho_and_drhoU, S_new_crse, Density_comp, 0, AMREX_SPACEDIM+1, 0);

            MultiFab dphi(grids, dmap, 1, 0);
            dphi.setVal(0);

            gravity->reflux_phi(level, dphi);

            // Compute (cross-level) gravity sync based on drho, dphi
            Vector<std::unique_ptr<MultiFab> > grad_delta_phi_cc(finest_level - level + 1);
            for (int lev = level; lev <= finest_level; lev++)
            {
                grad_delta_phi_cc[lev-level].reset(
                                      new MultiFab(get_level(lev).boxArray(),
                                                   get_level(lev).DistributionMap(),
                                                   AMREX_SPACEDIM, 0));
                grad_delta_phi_cc[lev-level]->setVal(0);
            }

            gravity->gravity_sync(level,finest_level,iteration,ncycle,drho_and_drhoU,dphi,
                                  amrex::GetVecOfPtrs(grad_delta_phi_cc));
            dphi.clear();

            for (int lev = level; lev <= finest_level; lev++)
            {
                Real dt_lev = parent->dtLevel(lev);
                MultiFab&  S_new_lev = get_level(lev).get_new_data(State_Type);
                Real cur_time = state[State_Type].curTime();
                Real a_new = get_comoving_a(cur_time);

                const auto& ba = get_level(lev).boxArray();
                const auto& dm = get_level(lev).DistributionMap();
                MultiFab grad_phi_cc(ba, dm, AMREX_SPACEDIM, 0);
                gravity->get_new_grav_vector(lev, grad_phi_cc, cur_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                  {
                  FArrayBox dstate;

                  for (MFIter mfi(S_new_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                  {
                    const Box& bx = mfi.tilebox();

                    dstate.resize(bx, AMREX_SPACEDIM + 1);
                    Array4<Real> d_fab = dstate.array();

                    if (lev == level)
                    {
                      dstate.copy<RunOn::Device>(drho_and_drhoU[mfi]);
                    }
                    else
                    {
                      dstate.setVal<RunOn::Device>(0);
                    }

                    Array4<Real const> const&  gphi = grad_phi_cc.array(mfi);
                    Array4<Real const> const& gdphi = grad_delta_phi_cc[lev-level]->array(mfi);
                    Array4<Real      > const& s_fab = S_new_lev.array(mfi);

                    int iden  = Density_comp;
                    int ieden = Eden_comp;

                    amrex::ParallelFor(bx, [s_fab,d_fab,gphi,gdphi,a_new,dt_lev,iden,ieden]
                       AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                           Real rho_pre  = s_fab(i,j,k,iden  ) - d_fab(i,j,k,0);
                           Real rhoU_pre = s_fab(i,j,k,iden+1) - d_fab(i,j,k,1);
                           Real rhoV_pre = s_fab(i,j,k,iden+2) - d_fab(i,j,k,2);
                           Real rhoW_pre = s_fab(i,j,k,iden+3) - d_fab(i,j,k,3);

                           Real SrU = d_fab(i,j,k,1)*gphi(i,j,k,0) + rho_pre*gdphi(i,j,k,0);
                           Real SrV = d_fab(i,j,k,1)*gphi(i,j,k,1) + rho_pre*gdphi(i,j,k,1);
                           Real SrW = d_fab(i,j,k,1)*gphi(i,j,k,2) + rho_pre*gdphi(i,j,k,2);

                           Real SrE = ( SrU * (rhoU_pre + (0.5*dt_lev)*SrU) +
                                        SrV * (rhoV_pre + (0.5*dt_lev)*SrV) +
                                        SrW * (rhoW_pre + (0.5*dt_lev)*SrW) ) / rho_pre;

                           Real a_new_inv = 1. / a_new;

                           s_fab(i,j,k,iden+1) += SrU * a_new_inv * 0.5 * dt_lev;
                           s_fab(i,j,k,iden+2) += SrV * a_new_inv * 0.5 * dt_lev;
                           s_fab(i,j,k,iden+3) += SrW * a_new_inv * 0.5 * dt_lev;
                           s_fab(i,j,k,ieden ) += SrE * a_new_inv * 0.5 * dt_lev;
                       });

                  }
                }
            }
        }
    }
#endif // end ifndef NO_HYDRO

    if (level < finest_level)
        average_down();

    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(do_reflux);
    BL_PROFILE_VAR("Nyx::post_timestep()::sum_write",sum_write);

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
#ifndef NO_HYDRO
        if ( (do_hydro == 1) && (sum_interval > 0) && (nstep % sum_interval == 0))
        {
            sum_integrated_quantities();
        }
#endif
        bool do_insitu = ((nstep+1) >= insitu_start) &&
            (insitu_int > 0) && ((nstep+1) % insitu_int == 0);

    if(do_insitu || doAnalysisNow())
            updateInSitu();

        write_info();

#ifdef BL_USE_MPI
        // Memory monitoring:
        MemInfo* mInfo = MemInfo::GetInstance();
        char info[32];
        snprintf(info, sizeof(info), "Step %4d", nstep);
        mInfo->LogSummary(info);
#endif
    }

    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(sum_write);
    BL_PROFILE_VAR("Nyx::post_timestep()::compute_temp",compute_temp);

#ifndef NO_HYDRO
    if (do_hydro)
    {
       MultiFab& S_new = get_new_data(State_Type);
       MultiFab& D_new = get_new_data(DiagEOS_Type);

       // First reset internal energy before call to compute_temp
       reset_internal_energy_nostore(S_new,D_new);

       // Re-compute temperature after all the other updates.
       compute_new_temp(S_new,D_new);
    }
#endif

    amrex::Gpu::Device::streamSynchronize();
    BL_PROFILE_VAR_STOP(compute_temp);

}

void
Nyx::typical_values_post_restart (const std::string& restart_file)
{
    if (level > 0)
        return;

    if (use_typical_steps)
    {
      if (ParallelDescriptor::IOProcessor())
        {
          std::string FileName = restart_file + "/first_max_steps";
          std::ifstream File;
          File.open(FileName.c_str(),std::ios::in);
          if (!File.good())
            amrex::FileOpenFailed(FileName);
          File >> old_max_sundials_steps;
        }
      ParallelDescriptor::Bcast(&old_max_sundials_steps, 1, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor())
        {
          std::string FileName = restart_file + "/second_max_steps";
          std::ifstream File;
          File.open(FileName.c_str(),std::ios::in);
          if (!File.good())
            amrex::FileOpenFailed(FileName);
          File >> new_max_sundials_steps;
        }
      ParallelDescriptor::Bcast(&new_max_sundials_steps, 1, ParallelDescriptor::IOProcessorNumber());
    }
}

void
Nyx::post_restart ()
{
    BL_PROFILE("Nyx::post_restart()");
#ifdef AMREX_PARTICLES
    if (level == 0)
        particle_post_restart(parent->theRestartFile());
#endif
    if (level == 0)
        comoving_a_post_restart(parent->theRestartFile());

    if (level == 0)
        typical_values_post_restart(parent->theRestartFile());

    if (inhomo_reion) init_zhi();

#ifndef NO_HYDRO
    Real cur_time = state[State_Type].curTime();
#else
    Real cur_time = state[PhiGrav_Type].curTime();
#endif

    // Update the value of a only if restarting from chk00000
    //   (special case for which computeNewDt is *not* called from Amr::coarseTimeStep)
    if (level == 0 && cur_time == 0.0)
        integrate_comoving_a(cur_time,parent->dtLevel(0));

#ifdef TISF
     int blub = parent->finestLevel();
     fort_set_finest_level(&blub);
#endif

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

            // Do multilevel solve here.  We now store phi in the checkpoint file so we can use it
            //  at restart.
            int ngrow_for_solve = 1;
            int use_previous_phi_as_guess = 1;
            gravity->multilevel_solve_for_new_phi(0,parent->finestLevel(),ngrow_for_solve,use_previous_phi_as_guess);

#ifndef AGN
            if (do_dm_particles)
#endif
            {
                for (int k = 0; k <= parent->finestLevel(); k++)
                {
                    const auto& ba = get_level(k).boxArray();
                    const auto& dm = get_level(k).DistributionMap();
                    MultiFab grav_vec_new(ba, dm, AMREX_SPACEDIM, 0);
                    gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
                }
            }
        }
    }

#ifndef NO_HYDRO
    if (do_forcing)
    {
        if (level == 0)
           forcing_post_restart(parent->theRestartFile());
    }

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

       const Real cur_time = state[State_Type].curTime();
       Real a = get_comoving_a(cur_time);

       Real average_temperature;
       compute_average_temperature(average_temperature);
       //
       // Get the number of species from the network model.
       //
#ifndef CONST_SPECIES
       NumSpec = 2;
#endif
       Real gamma_minus_1 = gamma - 1.0;
       set_small_values_given_average
           (average_gas_density, average_temperature,
           a, small_dens, small_temp, small_pres, gamma_minus_1, h_species);

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
   MultiFab::RegionTag amrPost_tag("Post_" + std::to_string(level));

#ifdef AMREX_PARTICLES
    if(load_balance_int >= 0 && nStep() % load_balance_int == 0 && (new_a >= 1.0/(load_balance_start_z + 1.0)))
    {
      if(verbose>0)
        amrex::Print()<<"Load balancing since "<<nStep()<<" mod "<<load_balance_int<<" == 0"<<std::endl;

    for (int lev = 0; lev <= parent->finestLevel(); lev++)
    {

        Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(lev));
        Vector<long> wgts(parent->boxArray(lev).size());
        DistributionMapping dm;

        //Should these weights be constructed based on level=0, or lev from for?
        if(load_balance_wgt_strategy == 0)
        {
            for (unsigned int i = 0; i < wgts.size(); i++)
            {
                wgts[i] = parent->boxArray(lev)[i].numPts();
            }
            if(load_balance_strategy==DistributionMapping::Strategy::KNAPSACK)
                dm.KnapSackProcessorMap(wgts, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::SFC)
                dm.SFCProcessorMap(parent->boxArray(lev), wgts, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::ROUNDROBIN)
                dm.RoundRobinProcessorMap(wgts, load_balance_wgt_nmax);
        }
        else if(load_balance_wgt_strategy == 1)
        {
            wgts = cs->theDMPC()->NumberOfParticlesInGrid(lev,false,false);
            if(load_balance_strategy==DistributionMapping::Strategy::KNAPSACK)
                dm.KnapSackProcessorMap(wgts, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::SFC)
              dm.SFCProcessorMap(parent->boxArray(lev), wgts, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::ROUNDROBIN)
                dm.RoundRobinProcessorMap(wgts, load_balance_wgt_nmax);
        }
        else if(load_balance_wgt_strategy == 2)
        {
            MultiFab particle_mf(parent->boxArray(lev),theDMPC()->ParticleDistributionMap(lev),1,1);
            cs->theDMPC()->Increment(particle_mf, lev);
            if(load_balance_strategy==DistributionMapping::Strategy::KNAPSACK)
                dm = DistributionMapping::makeKnapSack(particle_mf, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::SFC)
                dm = DistributionMapping::makeSFC(particle_mf, load_balance_wgt_nmax);
            else if(load_balance_strategy==DistributionMapping::Strategy::ROUNDROBIN)
                dm = DistributionMapping::makeRoundRobin(particle_mf);
        }
        else
        {
            amrex::Abort("Selected load balancing strategy not implemented");
        }

        amrex::Gpu::Device::streamSynchronize();
        const DistributionMapping& newdmap = dm;

        if(verbose > 2)
          amrex::Print()<<"Using ba: "<<parent->boxArray(lev)<<"\nUsing dm: "<<newdmap<<std::endl;        
        for (int i = 0; i < theActiveParticles().size(); i++)
        {
             if(lev > 0)
                 amrex::Abort("Particle load balancing needs multilevel testing");
          /*
             cs->theActiveParticles()[i]->Redistribute(lev,
                                                       theActiveParticles()[i]->finestLevel(),
                                                       1);
          */
             cs->theActiveParticles()[i]->Regrid(newdmap, parent->boxArray(lev), lev);

             if(shrink_to_fit)
                 cs->theActiveParticles()[i]->ShrinkToFit();
        }

        if(cs->Nyx::theVirtPC() != 0)
        {
            cs->Nyx::theVirtPC()->Regrid(newdmap, parent->boxArray(lev), lev);
        }

        if(cs->Nyx::theGhostPC() != 0)
        {
            cs->Nyx::theGhostPC()->Regrid(newdmap, parent->boxArray(lev), lev);
        }

        amrex::Gpu::streamSynchronize();
    }

   }
#endif
   AmrLevel::postCoarseTimeStep(cumtime);

#ifdef AGN
   if (level == 0)
     {
       agn_halo_find(parent->dtLevel(level));
     }
#endif

#ifdef GIMLET
   LyA_statistics();
#endif

   int nstep = parent->levelSteps(0);

   if (verbose>1)
   {
       amrex::Print() << "End of postCoarseTimeStep, printing:" <<std::endl;
       MultiFab::printMemUsage();
       amrex::Arena::PrintUsage();
   }
}

void
Nyx::post_regrid (int lbase,
                  int new_finest)
{
    BL_PROFILE("Nyx::post_regrid()");
#ifndef NO_HYDRO
#ifdef TISF
     fort_set_finest_level(&new_finest);
#endif
#endif
#ifdef AMREX_PARTICLES
    if (level == lbase) {
        particle_redistribute(lbase, false);
    }
    amrex::Gpu::Device::streamSynchronize();
#endif

    if (do_grav)
    {
        int which_level_being_advanced = parent->level_being_advanced();
        bool do_grav_solve_here;

        if (which_level_being_advanced >= 0)
        {
            do_grav_solve_here = (level == which_level_being_advanced) && (lbase == which_level_being_advanced);
        } else {
            do_grav_solve_here = (level == lbase);
        }

        //        if(parent->maxLevel() == 0)
        if(reuse_mlpoisson != 0)
          gravity->setup_Poisson(level,new_finest);

        // Only do solve here if we will be using it in the timestep right after without re-solving,
        //      or if this is called from somewhere other than Amr::timeStep
#ifndef NO_HYDRO
    const Real cur_time = state[State_Type].curTime();

#else
    const Real cur_time = state[PhiGrav_Type].curTime();
#endif
        if ((cur_time > 0) && do_grav_solve_here)
        {
            // We have done a particle redistribute above so we shouldn't need to grow by more than 1
            // to capture the effect of all particles on this level
            int ngrow_for_solve = 1;
            int use_previous_phi_as_guess = 1;
        gravity->multilevel_solve_for_new_phi(level, new_finest, ngrow_for_solve, use_previous_phi_as_guess);

#ifndef AGN
            if (do_dm_particles)
#endif
            {
                for (int k = 0; k <= parent->finestLevel(); k++)
                {
                    const auto& ba = get_level(k).boxArray();
                    const auto& dm = get_level(k).DistributionMap();
                    MultiFab grav_vec_new(ba, dm, AMREX_SPACEDIM, 0);
                    gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
                }
            }
        }
    }
    delete fine_mask;
    fine_mask = 0;
}

void
Nyx::post_init (Real /*stop_time*/)
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

    if (do_grav)
    {
#ifndef NO_HYDRO
        const Real cur_time = state[State_Type].curTime();
#else
        const Real cur_time = state[PhiGrav_Type].curTime();
#endif

        //
        // Calculate offset before first multilevel solve.
        //
        gravity->set_mass_offset(cur_time);

        //
        // Solve on full multilevel hierarchy
        //
        int ngrow_for_solve = 1;
        gravity->multilevel_solve_for_new_phi(0, finest_level, ngrow_for_solve);

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= finest_level; k++)
        {
            const auto& ba = get_level(k).boxArray();
            const auto& dm = get_level(k).DistributionMap();
            MultiFab grav_vec_new(ba, dm, AMREX_SPACEDIM, 0);
            gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
        }
    }

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
#ifndef NO_HYDRO
        const Real cur_time = state[State_Type].curTime();
#else
        const Real cur_time = state[PhiGrav_Type].curTime();
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
        fort_auxupdate
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

    get_flux_reg(level+1).Reflux(get_new_data(State_Type), 1.0, 0, 0, NUM_STATE,
                                 geom);
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

    if (do_grav)
    {
        average_down(PhiGrav_Type);
        average_down(Gravity_Type);
    }
}

#ifndef NO_HYDRO
#ifndef CONST_SPECIES
void
Nyx::enforce_nonnegative_species (MultiFab& S_new)
{
    BL_PROFILE("Nyx::enforce_nonnegative_species()");
    Real eps = -1.0e-16;
    int NumSpec = S_new.nComp()-FirstSpec_comp;

    if (NumSpec > 0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto uout = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               bool any_negative=false;
               Real x;
               Real dom_spec;
               int int_dom_spec;
               for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; n++)
               {
                   if (uout(i,j,k,n) < 0.0)
                   {
                       x = uout(i,j,k,n)/uout(i,j,k,Density_comp);
                       if (x > eps)
                           uout(i,j,k,n) = 0.0;
                       else
                           any_negative = true;
                   }
               }

               //
               // We know there are one or more undershoots needing correction.
               //
               if (any_negative)
               {
               //
               // Find the dominant species.
               //
                   int_dom_spec = FirstSpec_comp;
                   dom_spec     = uout(i,j,k,int_dom_spec);

                   for (int n = 0; n < FirstSpec_comp + NumSpec; n++)
                   {
                       if (uout(i,j,k,n) > dom_spec)
                       {
                           dom_spec     = uout(i,j,k,n);
                           int_dom_spec = n;
                       }
                   }
                   //
                   // Now take care of undershoots greater in magnitude than 1e-16.
                   //
                   for (int n = 0; n < FirstSpec_comp + NumSpec; n++)
                   {
                       if (uout(i,j,k,n) < 0.0)
                       {
                           x = uout(i,j,k,n)/uout(i,j,k,Density_comp);
                           //
                           // Take enough from the dominant species to fill the negative one.
                           //
                           uout(i,j,k,int_dom_spec) = uout(i,j,k,int_dom_spec) + uout(i,j,k,n);
                           //
                           // Test that we didn't make the dominant species negative.
                           //
                           //
                           // Now set the negative species to zero.
                           //
                           uout(i,j,k,n) = 0.00;
                       }
                   }
               }
            });
        }
    }
}
#endif

void
Nyx::enforce_consistent_e (MultiFab& S)
{
  BL_PROFILE("Nyx::enforce_consistent_e()");

  //  Set (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
        const Box& bx = mfi.tilebox();
        auto const s_arr = S.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            s_arr(i,j,k,Eden_comp) = s_arr(i,j,k,Eint_comp) + 0.5 * (
                                s_arr(i,j,k,Xmom_comp)*s_arr(i,j,k,Xmom_comp) +
                                s_arr(i,j,k,Ymom_comp)*s_arr(i,j,k,Ymom_comp) +
                                s_arr(i,j,k,Zmom_comp)*s_arr(i,j,k,Zmom_comp) ) / s_arr(i,j,k,Density_comp);
        });
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

    Nyx& fine_lev = get_level(level+1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom =          geom;

#ifndef NO_HYDRO
    if (state_index == State_Type)
    {
        MultiFab& S_crse =          get_new_data(State_Type);
        MultiFab& S_fine = fine_lev.get_new_data(State_Type);

        amrex::average_down(S_fine, S_crse,
                            fgeom, cgeom,
                            0, S_fine.nComp(), fine_ratio);

        MultiFab& D_crse =          get_new_data(DiagEOS_Type);
        MultiFab& D_fine = fine_lev.get_new_data(DiagEOS_Type);

        amrex::average_down(D_fine, D_crse,
                            fgeom, cgeom,
                            0, D_fine.nComp(), fine_ratio);
    }
    else
#endif
    {
      MultiFab& S_crse = get_new_data(state_index);
      MultiFab& S_fine = fine_lev.get_new_data(state_index);

      const int num_comps = S_fine.nComp();

      amrex::average_down(S_fine,S_crse,fgeom,cgeom,0,num_comps,fine_ratio);
    }
}

void
Nyx::errorEst (TagBoxArray& tags,
               int          clearval,
               int          tagval,
               Real         time,
               int        /*n_error_buf*/,
               int        /*ngrow*/)
{
    BL_PROFILE("Nyx::errorEst()");

    for (int j=0; j<errtags.size(); ++j) {
      std::unique_ptr<MultiFab> mf;
      if (errtags[0].Field() != std::string()) {
        mf = std::unique_ptr<MultiFab>(derive(errtags[j].Field(), time, errtags[j].NGrow()));
      }
      errtags[j](tags,mf.get(),clearval,tagval,time,level,geom);
    }
}

std::unique_ptr<MultiFab>
Nyx::derive (const std::string& name,
             Real               time,
             int                ngrow)
{
    BL_PROFILE("Nyx::derive()");


    if (name == "Rank")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids, dmap, 1, 0));
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*derive_dat,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const auto fab_derive_dat=derive_dat->array(mfi);
          const Box& bx = mfi.tilebox();
          Real my_proc = ParallelDescriptor::MyProc();
          AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
                                   {
                                     fab_derive_dat(i,j,k,0)=my_proc;
                                     });
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
    if (name == "Rank")
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const auto fab_derive_dat=mf.array(mfi);
          const Box& bx = mfi.tilebox();
          Real my_proc = ParallelDescriptor::MyProc();
          AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
                                   {
                                     fab_derive_dat(i,j,k,dcomp)=my_proc;
                                     });
        }
    } else {
        const auto& derive_dat = particle_derive(name, time, mf.nGrow());
        MultiFab::Copy(mf, *derive_dat, 0, dcomp, 1, mf.nGrow());
    }
}

#ifndef NO_HYDRO
void
Nyx::reset_internal_energy (MultiFab& S_new, MultiFab& D_new, MultiFab& reset_e_src)
{
    BL_PROFILE("Nyx::reset_internal_energy()");
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    const Real  cur_time = state[State_Type].curTime();
    Real        a        = get_comoving_a(cur_time);
    int interp=false;
    Real gamma_minus_1 = gamma - 1.0;
    auto atomic_rates = atomic_rates_glob;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto fab = S_new.array(mfi);
        const auto fab_diag = D_new.array(mfi);
        const auto fab_reset = reset_e_src.array(mfi);
        Real h_species_in=h_species;
        Real small_temp_in=small_temp;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        reset_internal_e
            (i,j,k,fab,fab_diag,fab_reset,
             atomic_rates, a, gamma_minus_1,h_species_in, small_temp_in, interp);
          });
    }
}

void
Nyx::reset_internal_energy_interp (MultiFab& S_new, MultiFab& D_new, MultiFab& reset_e_src)
{
    BL_PROFILE("Nyx::reset_internal_energy()");
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    const Real  cur_time = state[State_Type].curTime();
    Real        a        = get_comoving_a(cur_time);
    int interp=true;
    Real gamma_minus_1 = gamma - 1.0;
    auto atomic_rates = atomic_rates_glob;


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
          const auto fab = S_new.array(mfi);
          const auto fab_diag = D_new.array(mfi);
          const auto fab_reset = reset_e_src.array(mfi);
          Real h_species_in=h_species;
          Real small_temp_in=small_temp;
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        reset_internal_e
            (i,j,k,fab,fab_diag,fab_reset,
             atomic_rates, a, gamma_minus_1,h_species_in, small_temp_in, interp);
          });
    }

}

void
Nyx::reset_internal_energy_nostore (MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::reset_internal_energy_nostore()");
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
    reset_e_src.setVal(0.0);

    reset_internal_energy (S_new, D_new,  reset_e_src);
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_new_temp (MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::compute_new_temp()");

    Real cur_time  = state[State_Type].curTime();
    Real a        = get_comoving_a(cur_time);

    amrex::Gpu::synchronize();
    FArrayBox test, test_d;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            const Box& bx = mfi.tilebox();

            /*
            test.resize(bx,S_new.nComp());
            test_d.resize(bx,S_new.nComp());
            test.copy(S_new[mfi],bx);
            test_d.copy(D_new[mfi],bx);
            //      test.copy(S_new[mfi],0,0,S_new.nComp());
            //test_d.copy(D_new[mfi],0,0,D_new.nComp());*/

            const auto state_fab    = S_new.array(mfi);
            const auto diag_eos_fab = D_new.array(mfi);

            //      AMREX_LAUNCH_DEVICE_LAMBDA


            Real local_small_temp  = small_temp;
            Real local_large_temp  = large_temp;
            int  local_max_temp_dt = max_temp_dt;

            Real h_species_in=h_species;
            Real gamma_minus_1_in=gamma - 1.0;
            auto atomic_rates = atomic_rates_glob;
            AMREX_PARALLEL_FOR_3D(bx, i, j ,k,
            {
              Real rhoInv = 1.0 / state_fab(i,j,k,Density_comp);
              Real eint = state_fab(i,j,k,Eint_comp) * rhoInv;

            if (state_fab(i,j,k,Eint_comp) > 0.0)
            {
                eint = state_fab(i,j,k,Eint_comp) * rhoInv;

                int JH = 1;
                int JHe = 1;

                nyx_eos_T_given_Re_device(atomic_rates, gamma_minus_1_in, h_species_in, 1, 1, 
                                          &diag_eos_fab(i,j,k,Temp_comp), &diag_eos_fab(i,j,k,Ne_comp),
                                          state_fab(i,j,k,Density_comp), state_fab(i,j,k,Eint_comp) * (1.0 / state_fab(i,j,k,Density_comp)), a);

                if(diag_eos_fab(i,j,k,Temp_comp)>=local_large_temp && local_max_temp_dt == 1)
                {
                  diag_eos_fab(i,j,k,Temp_comp) = local_large_temp;

                  Real dummy_pres=0.0;
                  // Set temp to small_temp and compute corresponding internal energy
                  nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, &dummy_pres, 
                                   state_fab(i,j,k,Density_comp), diag_eos_fab(i,j,k,Temp_comp),
                                   diag_eos_fab(i,j,k,Ne_comp), a);

                   Real ke = 0.5e0 * (state_fab(i,j,k,Xmom_comp) * state_fab(i,j,k,Xmom_comp) +
                                      state_fab(i,j,k,Ymom_comp) * state_fab(i,j,k,Ymom_comp) +
                                      state_fab(i,j,k,Zmom_comp) * state_fab(i,j,k,Zmom_comp)) * rhoInv;

                   state_fab(i,j,k,Eint_comp) = state_fab(i,j,k,Density_comp) * eint;
                   state_fab(i,j,k,Eden_comp) = state_fab(i,j,k,Eint_comp) + ke;

                }
              }
            else
              {
                Real dummy_pres=0.0;
                // Set temp to small_temp and compute corresponding internal energy
                nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, &dummy_pres, 
                                 state_fab(i,j,k,Density_comp), local_small_temp,
                                 diag_eos_fab(i,j,k,Ne_comp), a);

                Real ke = 0.5e0 * (state_fab(i,j,k,Xmom_comp) * state_fab(i,j,k,Xmom_comp) +
                                   state_fab(i,j,k,Ymom_comp) * state_fab(i,j,k,Ymom_comp) +
                                   state_fab(i,j,k,Zmom_comp) * state_fab(i,j,k,Zmom_comp)) * rhoInv;

                diag_eos_fab(i,j,k,Temp_comp) = local_small_temp;
                state_fab(i,j,k,Eint_comp) = state_fab(i,j,k,Density_comp) * eint;
                state_fab(i,j,k,Eden_comp) = state_fab(i,j,k,Eint_comp) + ke;
              }

            });
            amrex::Gpu::synchronize();
          }

    // Find the cell which has the maximum temp -- but only if not the first
    // time step because in the first time step too many points have the same
    // value.
    Real prev_time = state[State_Type].prevTime();
    if (prev_time > 0.0 && verbose > 1)
    {
        BL_PROFILE("Nyx::compute_new_temp()::max_temp");
        // Compute the maximum temperature
        Real max_temp = D_new.norm0(Temp_comp);
        IntVect max_temp_loc = D_new.maxIndex(Temp_comp);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            if (bx.contains(max_temp_loc))
            {
                const auto fab_state = S_new.array(mfi);
                Real den_maxt=fab_state(max_temp_loc,Density_comp);
                std::cout << "Maximum temp. at level " << level << " is " << max_temp
                            << " at density " << den_maxt
                            << " at (i,j,k) " << max_temp_loc << std::endl;
            }
        }
    }
    amrex::Gpu::synchronize();
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_rho_temp (Real& rho_T_avg, Real& T_avg, Real& Tinv_avg, Real& T_meanrho)
{
    BL_PROFILE("Nyx::compute_rho_temp()");

    {

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real rho_T_sum=0.0,   T_sum=0.0, Tinv_sum=0.0, T_meanrho_sum=0.0;
    Real   rho_sum=0.0, vol_sum=0.0,    vol_mn_sum=0.0;

    Real rho_hi = 1.1*average_gas_density;
    Real rho_lo = 0.9*average_gas_density;
    const auto dx= geom.CellSizeArray();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        BL_PROFILE("Nyx::compute_rho_temp()::ReduceOpsOnDevice");
        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
                  ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
        ReduceData<Real,Real,Real,Real,Real,Real,Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto state    = S_new.const_array(mfi);
            const auto diag_eos = D_new.const_array(mfi);
            Real vol = dx[0]*dx[1]*dx[2];
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real T_tmp = diag_eos(i,j,k,Temp_comp);
                Real rho_tmp = state(i,j,k,Density_comp);
                Real l_rho_T_sum = 0._rt;
                Real l_rho_sum = 0._rt;
                Real l_T_sum = 0._rt;
                Real l_Tinv_sum = 0._rt;
                Real l_T_meanrho_sum = 0._rt;
                Real l_vol_sum = 0._rt;
                Real l_vol_mn_sum = 0._rt;
                l_T_sum = vol*T_tmp;
                l_Tinv_sum = rho_tmp/T_tmp;
                l_rho_T_sum = rho_tmp*T_tmp;
                l_rho_sum = rho_tmp;
                if ( (rho_tmp < rho_hi) &&  (rho_tmp > rho_lo) && (T_tmp <= 1.0e5) ) {
                    l_T_meanrho_sum = vol*log10(T_tmp);
                    l_vol_mn_sum = vol;
                }
                l_vol_sum = vol;
                return {l_rho_T_sum, l_rho_sum, l_T_sum, l_Tinv_sum, l_T_meanrho_sum, l_vol_sum, l_vol_mn_sum};
            });
        }

        ReduceTuple hv = reduce_data.value();
        rho_T_sum      = amrex::get<0>(hv);
        rho_sum        = amrex::get<1>(hv);
        T_sum          = amrex::get<2>(hv);
        Tinv_sum       = amrex::get<3>(hv);
        T_meanrho_sum  = amrex::get<4>(hv);
        vol_sum        = amrex::get<5>(hv);
        vol_mn_sum     = amrex::get<6>(hv);
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion() && !system::regtest_reduction)               \
    reduction(+:rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum, vol_sum, vol_mn_sum)
#endif
    {
        BL_PROFILE("Nyx::compute_rho_temp()::OnHost");
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto s_arr = S_new.const_array(mfi);
            const auto d_arr = D_new.const_array(mfi);
            Real vol = dx[0]*dx[1]*dx[2];
            AMREX_LOOP_3D(bx, i, j, k,
            {
                Real T_tmp   = d_arr(i,j,k,Temp_comp);
                Real rho_tmp = s_arr(i,j,k,Density_comp);

                T_sum     += vol*T_tmp;
                Tinv_sum  += rho_tmp/T_tmp;
                rho_T_sum += rho_tmp*T_tmp;
                rho_sum   += rho_tmp;

                if ( (rho_tmp < rho_hi) &&  (rho_tmp > rho_lo) && (T_tmp <= 1.0e5) ) {
                    T_meanrho_sum += vol*log10(T_tmp);
                    vol_mn_sum += vol;
                }
                vol_sum += vol;
            });
        }
    }
    }

    Real sums[7] = {rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum, vol_sum, vol_mn_sum};

    ParallelDescriptor::ReduceRealSum(sums,7);

    rho_T_avg = sums[0] / sums[1];  // density weighted T
        T_avg = sums[2] / sums[5];  // volume weighted T
     Tinv_avg = sums[3] / sums[1];  // 21cm T

    if (sums[6] > 0) {
       T_meanrho = sums[4] / sums[6];  // T at mean density
       T_meanrho = pow(10.0, T_meanrho);
    }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_gas_fractions (Real T_cut, Real rho_cut,
                            Real& whim_mass_frac, Real& whim_vol_frac,
                            Real& hh_mass_frac,   Real& hh_vol_frac,
                            Real& igm_mass_frac,  Real& igm_vol_frac)
{
    BL_PROFILE("Nyx::compute_gas_fractions()");

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real whim_mass=0.0, whim_vol=0.0, hh_mass=0.0, hh_vol=0.0;
    Real igm_mass=0.0, igm_vol=0.0, mass_sum=0.0, vol_sum=0.0;

    const auto dx= geom.CellSizeArray();

#ifdef AMREX_USE_GPU
    Real local_average_gas_density = average_gas_density;
    if (Gpu::inLaunchRegion())
    {
        BL_PROFILE("Nyx::compute_gas_fractions()::ReduceOpsOnDevice");
        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
                  ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
        ReduceData<Real,Real,Real,Real,Real,Real,Real,Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto state    = S_new.const_array(mfi);
            const auto diag_eos = D_new.const_array(mfi);
            Real vol = dx[0]*dx[1]*dx[2];
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real T = diag_eos(i,j,k,Temp_comp);
                Real R = state(i,j,k,Density_comp) / local_average_gas_density;
                Real rho_vol = state(i,j,k,Density_comp)*vol;
                Real l_whim_mass = 0._rt;
                Real l_whim_vol = 0._rt;
                Real l_hh_mass = 0._rt;
                Real l_hh_vol = 0._rt;
                Real l_igm_mass = 0._rt;
                Real l_igm_vol = 0._rt;
                if ( (T > T_cut) && (R <= rho_cut) ) {
                    l_whim_mass = rho_vol;
                    l_whim_vol = vol;
                }
                else if ( (T > T_cut) && (R > rho_cut) ) {
                    l_hh_mass = rho_vol;
                    l_hh_vol = vol;
                }
                else if ( (T <= T_cut) && (R <= rho_cut) ) {
                    l_igm_mass = rho_vol;
                    l_igm_vol = vol;
                }
                Real l_mass_sum = rho_vol;
                Real l_vol_sum = vol;
                return {l_whim_mass, l_whim_vol, l_hh_mass, l_hh_vol,
                        l_igm_mass, l_igm_vol, l_mass_sum, l_vol_sum};
            });
        }

        ReduceTuple hv = reduce_data.value();
        whim_mass = amrex::get<0>(hv);
        whim_vol  = amrex::get<1>(hv);
        hh_mass   = amrex::get<2>(hv);
        hh_vol    = amrex::get<3>(hv);
        igm_mass  = amrex::get<4>(hv);
        igm_vol   = amrex::get<5>(hv);
        mass_sum  = amrex::get<6>(hv);
        vol_sum   = amrex::get<7>(hv);
    }
    else
#endif
    {
    BL_PROFILE("Nyx::compute_gas_fractions()::OnHost");
#ifdef _OPENMP
#pragma omp parallel  if (amrex::Gpu::notInLaunchRegion())               \
    reduction(+:whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol, mass_sum, vol_sum)
#endif
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto s_arr    = S_new.const_array(mfi);
            const auto diag_eos = D_new.const_array(mfi);
            Real vol = dx[0]*dx[1]*dx[2];
            AMREX_LOOP_3D(bx, i, j, k,
            {
                Real T = diag_eos(i,j,k,Temp_comp);
                Real R = s_arr(i,j,k,Density_comp) / average_gas_density;
                Real rho_vol = s_arr(i,j,k,Density_comp)*vol;
                if ( (T > T_cut) && (R <= rho_cut) ) {
                    whim_mass += rho_vol;
                    whim_vol += vol;
                }
                else if ( (T > T_cut) && (R > rho_cut) ) {
                    hh_mass += rho_vol;
                    hh_vol += vol;
                }
                else if ( (T <= T_cut) && (R <= rho_cut) ) {
                    igm_mass += rho_vol;
                    igm_vol += vol;
                }
                mass_sum += rho_vol;
                vol_sum += vol;
            });
        }
    }

    Real sums[8] = {whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol, mass_sum, vol_sum};

    ParallelDescriptor::ReduceRealSum(sums,8);

    whim_mass_frac = sums[0] / sums[6];
    whim_vol_frac  = sums[1] / sums[7];
    hh_mass_frac   = sums[2] / sums[6];
    hh_vol_frac    = sums[3] / sums[7];
    igm_mass_frac  = sums[4] / sums[6];
    igm_vol_frac   = sums[5] / sums[7];
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
Nyx::InitErrorList() {
    //err_list.clear(true);
    //err_list.add("FULLSTATE",1,ErrorRec::Special,FORT_DENERROR);
}


//static Box the_same_box (const Box& b) { return b; }

void
Nyx::InitDeriveList() {
}


void
Nyx::LevelDirectoryNames(const std::string &dir,
                         const std::string &secondDir,
                         std::string &LevelDir,
                         std::string &FullPath)
{
    LevelDir = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    FullPath = dir;
    if( ! FullPath.empty() && FullPath.back() != '/') {
        FullPath += '/';
    }
    FullPath += secondDir;
    FullPath += "/";
    FullPath += LevelDir;
}


void
Nyx::CreateLevelDirectory (const std::string &dir)
{
    AmrLevel::CreateLevelDirectory(dir);  // ---- this sets levelDirectoryCreated = true
#ifdef AMREX_PARTICLES
    std::string dm(dir + "/" + Nyx::retrieveDM());
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(dm, 0755)) {
        amrex::CreateDirectoryFailed(dm);
      }
    }

    std::string LevelDir, FullPath;
    LevelDirectoryNames(dir, Nyx::retrieveDM(), LevelDir, FullPath);
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }

#ifdef AGN
    std::string agn(dir + "/" + Nyx::retrieveAGN());
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(agn, 0755)) {
        amrex::CreateDirectoryFailed(agn);
      }
    }

    LevelDirectoryNames(dir, Nyx::retrieveAGN(), LevelDir, FullPath);
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
#endif

#ifdef NEUTRINO_DARK_MATTER
    std::string npc(dir + "/" + Nyx::retrieveNPC());
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(npc, 0755)) {
        amrex::CreateDirectoryFailed(npc);
      }
    }

    LevelDirectoryNames(dir, Nyx::retrieveNPC(), LevelDir, FullPath);
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
#endif

    if(parent->UsingPrecreateDirectories()) {
      if(Nyx::theDMPC()) {
        Nyx::theDMPC()->SetLevelDirectoriesCreated(true);
      }
#ifdef AGN
      if(Nyx::theAPC()) {
        Nyx::theAPC()->SetLevelDirectoriesCreated(true);
      }
#endif
    }
#endif
}
