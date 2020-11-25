#include <iomanip>
#include <Nyx.H>
#include <Prob.H>

#ifdef GRAVITY
#include <Gravity.H>
#endif

using namespace amrex;

namespace
{
    std::string ascii_particle_file;
    std::string binary_particle_file;
    std::string    sph_particle_file;

#ifdef AGN
    std::string agn_particle_file;
#endif

#ifdef NEUTRINO_PARTICLES
    std::string neutrino_particle_file;
#endif
}

static int  do_santa_barbara = 0;
static int  init_sb_vels     = 1;
static int  do_readin_ics    = 0;
static int  fix_random_seed  = 0;
std::string readin_ics_fname;

#ifndef NO_HYDRO
void Nyx::set_small_values_given_average (amrex::Real average_dens, amrex::Real average_temp, amrex::Real a,
                                          amrex::Real & small_dens_inout, amrex::Real & small_temp_inout, 
                                          amrex::Real &small_pres_inout, Real gamma_minus_1_in, 
                                          Real h_species_in, Real He_species_in)
{

    Real small_dens, small_temp, frac;
    frac = 1e-6;
    if (small_dens_inout <= 0.e0)
        small_dens = frac * average_dens;
    else
        small_dens = small_dens_inout;

    if (small_temp_inout <= 0.e0)
        small_temp = frac * average_temp;
    else
        small_temp = small_temp_inout;

    // HACK HACK HACK -- FIX ME!!!!!
    ParmParse pp_nyx("nyx");
    const amrex::Real typical_Ne = 1.e0;
    Real small_pres=0.0;
    Real dummy_small_pres=0.0;
    Real eint=0.0;
    auto atomic_rates = atomic_rates_glob;
    if (!(pp_nyx.query("small_pres", dummy_small_pres)))
        nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, &small_pres, small_dens, small_temp, typical_Ne, a);
    else
        small_pres = small_pres_inout;

    small_dens_inout = small_dens;
    small_temp_inout = small_temp;
    small_pres_inout = small_pres;
}
#endif

void
Nyx::read_init_params ()
{
    BL_PROFILE("Nyx::read_init_params()");
    ParmParse pp("nyx");

    pp.query("do_santa_barbara", do_santa_barbara);
    pp.query("init_sb_vels", init_sb_vels);

    pp.query("fix_random_seed", fix_random_seed);
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    //      same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed)
        amrex::InitRandom(1024UL);

    if (do_hydro == 0 && do_santa_barbara == 1)
           amrex::Error("Nyx::cant have do_hydro == 0 and do_santa_barbara == 1");
    if (do_santa_barbara == 0 && init_with_sph_particles == 1)
           amrex::Error("Nyx::cant have do_santa_barbara == 0 and init_with_sph_particles == 1");
    if (do_santa_barbara == 0 && init_sb_vels == 1)
    {
       init_sb_vels = 0;
       if (ParallelDescriptor::IOProcessor())
           std::cout << "Nyx::setting init_sb_vels to 0 since do_santa_barbara = 0\n";
    }

    pp.query("do_readin_ics",       do_readin_ics);
    pp.query("readin_ics_fname", readin_ics_fname);
#ifdef AMREX_PARTICLES	
    pp.query("particle_launch_ics",       particle_launch_ics);
    pp.query("ascii_particle_file", ascii_particle_file);

    // Input error check
    if (do_dm_particles && !ascii_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified ascii_particle_file" << std::endl;
        amrex::Error();
    }

    pp.query("sph_particle_file", sph_particle_file);

    // Input error check
    if (init_with_sph_particles != 1 && !sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is not 1 but you specified sph_particle_file" << std::endl;
        amrex::Error();
    }

    // Input error check
    if (init_with_sph_particles == 1 && sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is 1 but you did not specify sph_particle_file" << std::endl;
        amrex::Error();
    }

    pp.query("binary_particle_file", binary_particle_file);

    // Input error check
    if (!binary_particle_file.empty() && (particle_init_type != "BinaryFile" &&
                                          particle_init_type != "BinaryMetaFile" && 
                                          particle_init_type != "BinaryMortonFile"))
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not BinaryFile, BinaryMetaFile, or BinaryMortonFile but you specified binary_particle_file" << std::endl;
        amrex::Error();
    }

#ifdef AGN
    pp.query("agn_particle_file", agn_particle_file);
    if (!agn_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified agn_particle_file" << std::endl;
        amrex::Error();
    }
#endif

#ifdef NEUTRINO_PARTICLES
    pp.query("neutrino_particle_file", neutrino_particle_file);
    if (!neutrino_particle_file.empty() && particle_init_type != "AsciiFile" &&
                                           particle_init_type != "BinaryMetaFile" && 
                                           particle_init_type != "BinaryFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile or BinaryFile but you specified neutrino_particle_file" << std::endl;
        amrex::Error();
    }
#endif
#endif
}

void
Nyx::init_zhi ()
{
    BL_PROFILE("Nyx::init_zhi()");

    if (ParallelDescriptor::IOProcessor()) std::cout << "Reading z_HI from file...";

    const int file_res = inhomo_grid;
    const int prob_res = geom.Domain().longside();
    const int ratio = prob_res / file_res;

    BL_ASSERT(ratio >= 1);
    
#ifndef NO_HYDRO
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    int nd = D_new.nComp();

    const BoxArray& my_ba = D_new.boxArray();
    const DistributionMapping& my_dmap = D_new.DistributionMap();

    BL_ASSERT(my_ba.coarsenable(ratio));
    BoxArray coarse_ba = my_ba;
    coarse_ba.coarsen(ratio);
    MultiFab zhi(coarse_ba, my_dmap, 1, 0);

    MultiFab zhi_from_file;
    VisMF::Read(zhi_from_file, inhomo_zhi_file);
    zhi.copy(zhi_from_file, geom.periodicity());
    if(D_new.nComp()>2)
    {
        int l_Zhi_comp = Zhi_comp;
        for (MFIter mfi(D_new,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_zhi=zhi.array(mfi);
                const auto fab_D_new=D_new.array(mfi);

                amrex::ParallelFor(
                               bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                   fab_D_new(i, j ,k, l_Zhi_comp) = fab_zhi(i/ratio,j/ratio,k/ratio);
                               });

            }
    }
#endif
    if (ParallelDescriptor::IOProcessor()) std::cout << "done.\n";
}

void
Nyx::initData ()
{
    BL_PROFILE("Nyx::initData()");
    
    amrex::Gpu::LaunchSafeGuard lsg(true);
    // Here we initialize the grid data and the particles from a plotfile.
    if (!parent->theRestartPlotFile().empty())
    {
        amrex::Abort("AmrData requires fortran");
        return;
    }

#ifndef NO_HYDRO    
    MultiFab&   S_new    = get_new_data(State_Type);

    // We need this because otherwise we might operate on uninitialized data.
    S_new.setVal(0.0);
#endif

    // If you run a pure N-body simulation and Nyx segfaults here, then
    // please check Prob_3d.f90... You might set other variables then density...
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Initializing the data at level " << level << '\n';

    const auto dx = geom.CellSizeArray();
    const auto geomdata = geom.data();

    // Make sure dx = dy = dz -- that's all we guarantee to support
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL) || (fabs(dx[0] - dx[2]) > SMALL) )
        amrex::Abort("We don't support dx != dy != dz");

#ifndef NO_HYDRO    
    Real  cur_time = state[State_Type].curTime();
#ifdef AMREX_PARTICLES
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) && (particle_init_type != "Cosmological") )
#else
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) )
#endif
    {
        if (do_hydro == 1) 
        {
            MultiFab&   D_new    = get_new_data(DiagEOS_Type);
            D_new.setVal(0., Temp_comp);
            D_new.setVal(0.,   Ne_comp);

            for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_S_new=S_new.array(mfi);
                const auto fab_D_new=D_new.array(mfi);

                GpuArray<amrex::Real,max_prob_param> prob_param;
                prob_param_fill(prob_param, initial_z);

                amrex::ParallelFor(
                                   bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                   prob_initdata
                                       (i, j ,k, fab_S_new, fab_D_new, geomdata, prob_param);
                               });
            }

            if (inhomo_reion) init_zhi();

            // First reset internal energy before call to compute_temp
            MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
            reset_e_src.setVal(0.0);

            reset_internal_energy(S_new,D_new,reset_e_src);
            compute_new_temp     (S_new,D_new);

            // Define (rho E) given (rho e) and the momenta
            enforce_consistent_e(S_new);
        }
        else
        {
            for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_S_new=S_new.array(mfi);                

                GpuArray<amrex::Real,max_prob_param> prob_param;
                prob_param_fill(prob_param, initial_z);
                amrex::ParallelFor(
                                   bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                   prob_initdata_state
                                      (i, j ,k, fab_S_new, geomdata, prob_param);
                               });
            }
        }
    }
#else
        Real  cur_time = state[PhiGrav_Type].curTime();
#endif // end NO_HYDRO

#ifdef GRAVITY

    if (!do_grav)
    {
        //
        // Set these to zero so they're defined for the plotfile.
        //
        MultiFab& G_new = get_new_data(Gravity_Type);
        G_new.setVal(0);
    }
    else 
    {
        //
        // Initialize this to zero before first solve.
        //
        MultiFab& Phi_new = get_new_data(PhiGrav_Type);
        Phi_new.setVal(0.);
    }

#endif

#ifdef SDC
    //
    // Initialize this to zero before we use it in advance
    //
    MultiFab& IR_new = get_new_data(SDC_IR_Type);
    IR_new.setVal(0.0);
#endif

#ifndef NO_HYDRO
    //
    // Read in initial conditions from a file.
    // By now only for fixed grid ics.
    // Layout and units have to be as in \vec U.
    //
    if (do_readin_ics)
    {
        std::string mfDirName(readin_ics_fname);

        MultiFab mf;
        mfDirName.append("/Level_0/Cell");

        VisMF::Read(mf, mfDirName.c_str());

        MultiFab& S_new_crse = get_level(0).get_new_data(State_Type);
        
        S_new_crse.copy(mf, 0, 0, 6);
        S_new_crse.copy(mf, 0, FirstSpec, 1);

        if (do_hydro == 1) 
        {
            MultiFab&  D_new = get_new_data(DiagEOS_Type);
            D_new.setVal(0.);
        }

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Readin stuff...done\n";
    }
#endif


    // Add partially redundant setup for small_dens
    Real average_gas_density = -1e200;
    Real average_temperature = -1e200;

    // Make sure small values give non-negative small_pres                                                           
    Real small_dens_loc = amrex::max(1e-2,small_dens);
    Real small_temp_loc = amrex::max(1e-2,small_temp);
    Real a = get_comoving_a(cur_time);
    Real gamma_minus_1 = gamma - 1.0;

    amrex::Gpu::Device::synchronize();

#ifndef NO_HYDRO
    set_small_values_given_average
      (average_gas_density, average_temperature,
       a,  small_dens_loc, small_temp_loc, small_pres, gamma_minus_1, h_species, he_species);
#endif

    amrex::Gpu::Device::synchronize();

#ifdef AMREX_PARTICLES
    if (level == 0)
        init_particles();

    amrex::Gpu::Device::synchronize();

    if ( particle_init_type == "Cosmological")
        initcosmo();

    //
    // Must redistribute particles before calling `init_santa_barbara` so that
    // the particles already live on the higher level when we go to put some of
    // the mass onto the grid.
    //
    if (level > 0)
        particle_redistribute();
#endif

    //
    // With this call we define the initial data on the current level but we
    // also may need to modify the data on the level below since particles
    // previous at the coarser level may now live at the finer level and
    // distribute their mass differently.
    //
#ifndef NO_HYDRO
#ifdef GRAVITY
    if (do_santa_barbara == 1)
        init_santa_barbara(init_sb_vels);
#endif

    if (do_hydro)
        check_initial_species();

    //
    // Need to compute this in case we want to use overdensity for regridding.
    //
    if (level == 0) 
        compute_average_density();
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Done initializing the level " << level << " data\n";
}

#if 0
void
Nyx::init_from_plotfile ()
{
    BL_PROFILE("Nyx::init_from_plotfile()");
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << " " << std::endl; 
        std::cout << "Initializing the data from " << parent->theRestartPlotFile() << std::endl;
    }

    if (parent->maxLevel() > 0)
        amrex::Abort("We can only restart from single-level plotfiles");

    // Make sure to read in "a" before we call ReadPlotFile since we will use a 
    //      when we construct e from T.

    bool is_checkpoint = false;

    // Now read in the time as well as the grid and particle data
    bool first = true;
    bool rhoe_infile;
    ReadPlotFile(first,parent->theRestartPlotFile(),rhoe_infile);

    // This is just a dummy value so that we can set the current time of the StateData
    Real dummy_dt = 1.e100;
    setTimeLevel(parent->cumTime(), dummy_dt, dummy_dt);

    comoving_a_post_restart(parent->theRestartPlotFile());

#ifndef NO_HYDRO
    for (int lev = 0; lev <= parent->finestLevel(); ++lev)
    {
        Nyx& nyx_lev = get_level(lev);
        MultiFab& S_new = nyx_lev.get_new_data(State_Type);
        MultiFab& D_new = nyx_lev.get_new_data(DiagEOS_Type);
        int ns = S_new.nComp();
        int nd = D_new.nComp();

        // Construct internal energy given density, temperature and species
        if (not rhoe_infile)
           init_e_from_T(old_a);

        // Define (rho E) given (rho e) and the momenta
        nyx_lev.enforce_consistent_e(S_new);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Done initializing the grid data from " << parent->theRestartPlotFile() << std::endl;
#endif

    // Now read the particles from the plotfile
    particle_post_restart(parent->theRestartPlotFile(),is_checkpoint);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Done initializing the particles from the plotfile " << std::endl;
        std::cout << " " << std::endl; 
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::check_initial_species ()
{
    // Verify that the sum of (rho X)_i = rho at every cell.

    MultiFab&   S_new    = get_new_data(State_Type);

    int iden  = Density;

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef CONST_SPECIES
    if (amrex::Math::abs(1.0 - h_species - he_species) > 1.e-8)
        amrex::Abort("Error:: Failed check of initial species summing to 1");
#else
    if (FirstSpec > 0)
    {
        int iufs = FirstSpec;
        int nspec = NumSpec;
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const auto state = S_new.array(mfi);

          const Box& tbx = mfi.tilebox();

          reduce_op.eval(tbx, reduce_data,
          [state,nspec,iden,iufs]
          AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {
               Real sum = state(i,j,k,iufs);
               for (int n = 1; n < nspec; n++)
                  sum += state(i,j,k,iufs+n);

               sum /= state(i,j,k,iden);

               Real x = amrex::Math::abs(amrex::Math::abs(sum) - 1.);
               return x;
          });
        }

        ReduceTuple hv = reduce_data.value();
        ParallelDescriptor::ReduceRealMax(amrex::get<0>(hv));
        if (get<0>(hv) > 1.e-8)
            amrex::Abort("Error:: Failed check of initial species summing to 1");
    }
#endif
}

void
Nyx::init_e_from_T (const Real& a)
{

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    Real h_species_in=h_species;
    Real gamma_minus_1_in=gamma - 1.0;
    auto atomic_rates = atomic_rates_glob;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto state    = S_new.array(mfi);
        const auto diag_eos = D_new.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
             Real rho =  state(i,j,k,Density);
             Real T   = diag_eos(i,j,k,Temp_comp);
             Real ne  = diag_eos(i,j,k,  Ne_comp);

             Real eint;

             // Call EOS to get the internal energy
             Real dummy_pres=0.0;
             // Set temp to small_temp and compute corresponding internal energy
             nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, &dummy_pres, rho, T, ne, a);

             state(i,j,k,Eint) = state(i,j,k,Density) * eint;

             state(i,j,k,Eden) = state(i,j,k,Eint) + 0.5 * ( 
                state(i,j,k,Xmom)*state(i,j,k,Xmom) +
                state(i,j,k,Ymom)*state(i,j,k,Ymom) +
                state(i,j,k,Zmom)*state(i,j,k,Zmom) ) / state(i,j,k,Density);
        });
        amrex::Gpu::Device::streamSynchronize();
    }
}
#endif
