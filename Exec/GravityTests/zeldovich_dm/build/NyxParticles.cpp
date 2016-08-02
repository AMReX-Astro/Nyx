#include <iomanip>
#include "Nyx.H"

#ifdef GRAVITY
#include "Gravity.H"
#include "Gravity_F.H" //needed for get_grav_constant, but there might be a better way
#endif

#include "Nyx_F.H"
#include <Particles_F.H>

static std::string ascii_particle_file;

// There's really only one of these.
static DarkMatterParticleContainer* DMPC = 0;
DarkMatterParticleContainer* Nyx::theDMPC ()
{
    return DMPC;
}

bool Nyx::do_dm_particles = false;

std::string Nyx::particle_init_type = "";
std::string Nyx::particle_move_type = "";

bool Nyx::particle_initrandom_serialize = false;
Real Nyx::particle_initrandom_mass;
long Nyx::particle_initrandom_count;
int  Nyx::particle_initrandom_iseed;

int Nyx::particle_verbose            = 1;
int Nyx::write_particles_in_plotfile = 0;

Real Nyx::particle_cfl = 0.5;

static const std::string chk_particle_file("DM");

void
Nyx::read_particle_params()
{
    ParmParse pp("nyx");
    pp.query("do_dm_particles", do_dm_particles);

    if (do_dm_particles)
    {
        pp.get("particle_init_type", particle_init_type);
        pp.get("particle_move_type", particle_move_type);
    }

#ifdef GRAVITY
    if (!do_grav && particle_move_type == "Gravitational")
    {
        std::cerr << "ERROR:: doesnt make sense to have do_grav=false but move_type = Gravitational"
                  << std::endl;
        BoxLib::Error();
    }
#endif

    pp.query("particle_initrandom_serialize", particle_initrandom_serialize);
    pp.query("particle_initrandom_count", particle_initrandom_count);
    pp.query("particle_initrandom_mass", particle_initrandom_mass);
    pp.query("particle_initrandom_iseed", particle_initrandom_iseed);

    pp.query("ascii_particle_file", ascii_particle_file);

    if (!ascii_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified ascii_particle_file"
                  << std::endl;;
        BoxLib::Error();
    }

    //
    // Control the verbosity of the Particle class
    //
    ParmParse ppp("particles");
    ppp.query("v", particle_verbose);
    ppp.query("write_in_plotfile", write_particles_in_plotfile);

    //
    // Set the cfl for particle motion (fraction of cell that a particle can
    // move in a timestep).
    //
    ppp.query("cfl", particle_cfl);
}

void
Nyx::init_particles()
{
    if (level > 0)
        return;

    //
    // Need to initialize particles before defining gravity.
    //
    if (do_dm_particles)
    {
        BL_ASSERT(DMPC == 0);
        DMPC = new DarkMatterParticleContainer(parent);

        //
        // 2 gives more stuff than 1.
        //
        DMPC->SetVerbose(particle_verbose);

        if (particle_init_type == "Random")
        {
            if (particle_initrandom_count <= 0)
            {
                BoxLib::Abort("Nyx::init_particles(): particle_initrandom_count must be > 0");
            }
            if (particle_initrandom_iseed <= 0)
            {
                BoxLib::Abort("Nyx::init_particles(): particle_initrandom_iseed must be > 0");
            }

            if (verbose && ParallelDescriptor::IOProcessor())
            {
                std::cout << "\nInitializing DM with cloud of "
                          << particle_initrandom_count
                          << " random particles with initial seed: "
                          << particle_initrandom_iseed << "\n\n";
            }

            DMPC->InitRandom(particle_initrandom_count,
                             particle_initrandom_iseed,
                             particle_initrandom_mass,
                             particle_initrandom_serialize);
        }
#ifdef GRAVITY
        else if (particle_init_type == "Cosmological")
        {
            Real comoving_OmM, comoving_OmL, comoving_h;
            Real comoving_OmB, Gconst;
            const Real len[BL_SPACEDIM] = {geom.ProbLength(0),geom.ProbLength(1),geom.ProbLength(2)};

            Real particleMass;
            std::string mfDirName;
            //
            // Read the init directory name and particle mass from the inputs file
            //
            ParmParse pp("cosmo");
            pp.get("initDirName", mfDirName);
            pp.get("particle_mass", particleMass);

            //
            // Read from the directory into mf
            //
            MultiFab mf;
            mfDirName.append("/Level_0/Cell");
            VisMF::Read(mf,mfDirName.c_str());
            if (mf.contains_nan(Density, mf.nComp(), 0))
            {
                for (int i = 0; i < mf.nComp(); i++)
                {
                    if (mf.contains_nan(Density + i, 1, 0))
                    {
                        std::cout << "Found NaNs in component " << i << ". " << std::endl;
                        BoxLib::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
                    }
                }
            }

            BL_FORT_PROC_CALL(GET_OMM, get_omm)(&comoving_OmM);
            BL_FORT_PROC_CALL(GET_OML, get_oml)(&comoving_OmL);
            BL_FORT_PROC_CALL(GET_OMB, get_omb)(&comoving_OmB);
            BL_FORT_PROC_CALL(GET_HUBBLE, get_hubble)(&comoving_h);
	    BL_FORT_PROC_CALL(GET_GRAV_CONST, get_grav_const)(&Gconst);

	    //For whatever reason OmB is divided by OmM, which I need to undo here...
	    Real realOmB = comoving_OmB*comoving_OmM;

	    //compute \Omega_K - just for checking flatness...
	    Real comoving_OmK = 1 - comoving_OmM - comoving_OmL;
	    if ( std::abs(comoving_OmK) > 1e-3 ){
		    std::cout << "You assume a non-flat universe - \\Omega_K = "
			      << comoving_OmK << std::endl;
		    BoxLib::Abort();
	    }

	    //compute \rho_{baryon}=3H_0^2\Omega_{baryon}/(8\pi G)
	    Real rhoB = 3*comoving_h*100*comoving_h*100*realOmB
		        / (8*M_PI*Gconst);
	    std::cout << "Mean baryonic matter density is " << rhoB
		      << " M_sun/Mpc^3." << std::endl;

	    //compute \rho_{dark matter}=3H_0^2\Omega_{dark matter}/(8\pi G)
	    Real comoving_OmD = comoving_OmM - realOmB;
	    Real rhoD = 3*comoving_h*100*comoving_h*100*comoving_OmD
		        / (8*M_PI*Gconst);
	    std::cout << "Mean dark matter density is " << rhoD
		      << " M_sun/Mpc^3." << std::endl;

	    //Compute particle mass
	    Real simulationVolume  = len[0]*len[1]*len[2];
	    const Real* dx = geom.CellSize();
	    Real  numberOfParticles = (len[0]/dx[0]) * (len[1]/dx[1]) * (len[2]/dx[2]);	//FIXME: In future there might be more particles than cells...
	    particleMass = rhoD * simulationVolume / numberOfParticles;
	    std::cout << "Particle mass is " << particleMass
		      << " M_sun." << std::endl;

	    if(!do_hydro && realOmB>0){
		    std::cout << std::endl;
		    std::cout << std::endl;
		    std::cout << "You chose a baryonic matter content greater than 0, but chose not to do hydro" << std::endl;
		    std::cout << "I will be using \\rho_M for the dark matter density..." << std::endl;
		    rhoD += rhoB;
		    std::cout << std::endl;
		    std::cout << std::endl;
	    }


            // Velocities are proportional to displacements by
            // LBox [MPc] * a * H [km/s/MPc]
            // we have to calculate the initial a on our own
            // as there is no code path to get_comoving_a (Castro.H sources Particles.H)
            Real redshift=-1;
            ParmParse pp2("nyx");
            pp2.get("initial_z",redshift);

            Real comoving_a = 1/(1+redshift);

            Real vel_fac[BL_SPACEDIM];
            for(int n=0;n<BL_SPACEDIM;n++)
                vel_fac[n]=len[n]*comoving_a*std::sqrt(comoving_OmM/pow(comoving_a,3)+comoving_OmL)*comoving_h*100;

	    //initialize particles
            DMPC->InitCosmo(mf, vel_fac, particleMass);

	    //initialize hydro
	    //units for velocity should be okay
	    if(do_hydro){
            	MultiFab& S_new = get_level(0).get_new_data(State_Type);
    	     	//copy density
	     	S_new.copy(mf, 3, Density, 1);
	     	S_new.plus(1,     Density, 1, S_new.nGrow());
	     	S_new.mult(rhoB,  Density, 1, S_new.nGrow());

//		//This block assigns "the same" density for the baryons as for the dm.
//              PArray<MultiFab> particle_mf;
//              DMPC->AssignDensity(particle_mf);
//              particle_mf[0].mult(realOmB / comoving_OmD);
//              S_new.copy(particle_mf[0], 0, Density, 1);

	     	//copy velocities...
	     	S_new.copy(mf, 0, Xmom, 3);

		//...and "transform" to momentum
	     	S_new.mult(vel_fac[0], Xmom, 1, S_new.nGrow());
	     	MultiFab::Multiply(S_new, S_new, Density, Xmom, 1, S_new.nGrow());
	     	S_new.mult(vel_fac[1], Ymom, 1, S_new.nGrow());
	     	MultiFab::Multiply(S_new, S_new, Density, Ymom, 1, S_new.nGrow());
	     	S_new.mult(vel_fac[2], Zmom, 1, S_new.nGrow());
		MultiFab::Multiply(S_new, S_new, Density, Zmom, 1, S_new.nGrow());

             	for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
             	{
             	    const Box& box = mfi.validbox();
             	    const int* lo = box.loVect();
             	    const int* hi = box.hiVect();
		    Real a = old_a;
		    int	ns = NUM_STATE;

		    // Temp unused for GammaLaw, set it here so that pltfiles have
		    // defined numbers
              	    S_new[mfi].setVal(100., Temp);

		    //is this the correct function?!
             	    BL_FORT_PROC_CALL(INIT_E_FROM_T, init_e_from_t)
              	       (BL_TO_FORTRAN(S_new[mfi]), &ns, lo, hi, &a);
             	}
	    }

        }
#endif
        else if (particle_init_type == "AsciiFile")
        {
            if (verbose && ParallelDescriptor::IOProcessor())
            {
                std::cout << "\nInitializing DM particles from \""
                          << ascii_particle_file << "\" ...\n\n";
            }
            //
            // The second argument is how many Reals we read into `m_data[]`
            // after reading in `m_pos[]`. Here we're reading in the particle
            // mass and velocity.
            //
            DMPC->InitFromAsciiFile(ascii_particle_file, BL_SPACEDIM + 1);
        }
        else
        {
            BoxLib::Error("not a valid input for nyx.particle_init_type");
        }
    }
}

#ifdef GRAVITY
void
Nyx::init_santa_barbara()
{
    Real frac_for_hydro;
    BL_FORT_PROC_CALL(GET_OMB, get_omb)(&frac_for_hydro);

    Real omfrac = 1. - frac_for_hydro;

    if (level == 0 && frac_for_hydro != 1.0)
        DMPC->MultiplyParticleMass(level, omfrac);

    PArray<MultiFab> particle_mf;
    DMPC->AssignDensity(particle_mf);
    for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
    {
        const IntVect ratio = parent->refRatio(lev);
        Gravity::average_down(particle_mf[lev], particle_mf[lev+1], ratio);
    }

    if (frac_for_hydro == 1.0)
    {
       for (int lev = 0; lev <= level; lev++)
          particle_mf[lev].mult(0.0);
    }
    else
    {
       for (int lev = 0; lev <= level; lev++)
          particle_mf[lev].mult(frac_for_hydro / omfrac);
    }

    int ns = NUM_STATE;

    Real cur_time = state[State_Type].curTime();
    Real a = old_a;

    if (ParallelDescriptor::IOProcessor())
        std::cout << "... time and comoving a when data is initialized "
                  << cur_time << " " << a << std::endl;

    for (int lev = 0; lev <= level; lev++)
    {
        MultiFab& S_new = get_level(lev).get_new_data(State_Type);
        const Real * dx = get_level(lev).geom.CellSize();
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            RealBox gridloc = RealBox(get_level(lev).grids[mfi.index()],
                                      get_level(lev).geom.CellSize(),
                                      get_level(lev).geom.ProbLo());
            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();

            // Temp unused for GammaLaw, set it here so that pltfiles have
            // defined numbers
            S_new[mfi].setVal(0., Temp);

            BL_FORT_PROC_CALL(CA_INITDATA, ca_initdata)
                (lev, cur_time, lo, hi, ns, BL_TO_FORTRAN(S_new[mfi]), dx,
                 gridloc.lo(), gridloc.hi());
        }

        MultiFab::Add(S_new, particle_mf[lev], 0, Density, 1, S_new.nGrow());

        // Make sure we've finished initializing the density before calling this.
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const int* lo = box.loVect();
            const int* hi = box.hiVect();
            BL_FORT_PROC_CALL(INIT_E_FROM_T, init_e_from_t)
                (BL_TO_FORTRAN(S_new[mfi]), &ns, lo, hi, &a);
        }
    }
}
#endif

void
Nyx::particle_check_point(const std::string& dir)
{
    if (level == 0)
    {
        if (DMPC)
            DMPC->Checkpoint(dir, chk_particle_file);

        Real cur_time = state[State_Type].curTime();

        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if (!File.good())
                BoxLib::FileOpenFailed(FileName);
            File.precision(15);
            if (cur_time == 0.0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
        }
    }
}

void
Nyx::particle_plot_file(const std::string& dir)
{
    if (level == 0)
    {
        if (DMPC && write_particles_in_plotfile)
            DMPC->Checkpoint(dir, chk_particle_file);

        Real cur_time = state[State_Type].curTime();
        DMPC->WriteAsciiFile(dir + "/particles.ascii");

        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if (!File.good())
                BoxLib::FileOpenFailed(FileName);
            File.precision(15);
            if (cur_time == 0.0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
            File.close();
        }
    }
}

void
Nyx::particle_post_restart(const std::string& restart_file)
{
    if (level == 0)
    {
        if (do_dm_particles)
        {
            BL_ASSERT(DMPC == 0);
            DMPC = new DarkMatterParticleContainer(parent);

            //
            // 2 gives more stuff than 1.
            //
            DMPC->SetVerbose(particle_verbose);
            DMPC->Restart(restart_file, chk_particle_file);

            //
            // We want the ability to write the particles out to an ascii file.
            //
            ParmParse pp("particles");

            std::string particle_output_file;

            pp.query("particle_output_file", particle_output_file);

            if (!particle_output_file.empty())
            {
                DMPC->WriteAsciiFile(particle_output_file);
            }
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = parent->theRestartFile() + "/comoving_a";
            std::ifstream File;
            File.open(FileName.c_str(),std::ios::in);
            if (!File.good())
                BoxLib::FileOpenFailed(FileName);
            File >> old_a;
        }
        ParallelDescriptor::Bcast(&old_a, 1,
                                  ParallelDescriptor::IOProcessorNumber());
        new_a = old_a;
    }
}

#ifdef GRAVITY
void
Nyx::particle_est_time_step(Real& est_dt)
{
    if (DMPC && particle_move_type == "Gravitational")
    {
        const Real cur_time = state[State_Type].curTime();
        const Real a = get_comoving_a(cur_time);
        const MultiFab& grav = get_new_data(Gravity_Type);
        const Real est_dt_particle = DMPC->estTimestep(grav, a, level,
                                                        particle_cfl);

        if (est_dt_particle > 0.0)
            est_dt = std::min(est_dt, est_dt_particle);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            if (est_dt_particle > 0.0)
            {
                std::cout << "...estdt from particles at level " << level
                          << ": " << est_dt_particle << '\n';
            }
            else
            {
                std::cout << "...there are no particles at level " << level
                          << std::endl;
            }
        }
    }
}
#endif

void
Nyx::particle_redistribute()
{
    if (DMPC)
        DMPC->Redistribute();
}

void
Nyx::particle_move_random()
{
    if (DMPC && particle_move_type == "Random")
    {
        BL_ASSERT(level == 0);

        DMPC->MoveRandom();
    }
}

MultiFab*
Nyx::particle_derive(const std::string& name, Real time, int ngrow)
{
    if (DMPC && name == "particle_count")
    {
        MultiFab* derive_dat = new MultiFab(grids, 1, 0);
        MultiFab temp_dat(grids, 1, 0);
        temp_dat.setVal(0);
        DMPC->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
    else if (DMPC && name == "total_particle_count")
    {
        //
        // We want the total particle count at this level or higher.
        //
        MultiFab* derive_dat = particle_derive("particle_count", time, ngrow);
        IntVect trr(D_DECL(1, 1, 1));

        // @todo: level vs. lev
        for (int lev = level + 1; lev <= parent->finestLevel(); lev++)
        {
            BoxArray ba = parent->boxArray(lev);
            MultiFab temp_dat(ba, 1, 0);

            trr *= parent->refRatio(lev - 1);

            ba.coarsen(trr);
            MultiFab ctemp_dat(ba, 1, 0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            DMPC->Increment(temp_dat, lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab = temp_dat[mfi];
                FArrayBox& cfab = ctemp_dat[mfi];
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == BoxLib::coarsen(fbx, trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
                    const Real val = ffab(p);
                    if (val > 0)
                        cfab(BoxLib::coarsen(p, trr)) += val;
                }
            }

            temp_dat.clear();

            MultiFab dat(grids, 1, 0);
            dat.setVal(0);
            dat.copy(ctemp_dat);

            MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
        }

        return derive_dat;
    }
#ifdef GRAVITY
    else if (DMPC && name == "particle_mass_density")
    {
        MultiFab* derive_dat = new MultiFab(grids,1,0);

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        PArray<MultiFab> particle_mf;
        DMPC->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            const IntVect ratio = parent->refRatio(lev);
            Gravity::average_down(particle_mf[lev], particle_mf[lev+1], ratio);
        }

        MultiFab::Copy(*derive_dat, particle_mf[level], 0, 0, 1, 0);

        return derive_dat;
    }
#endif

    else if (name == "total_density")
    {
#ifdef GRAVITY
      if (DMPC)
      {
        MultiFab* derive_dat = new MultiFab(grids,1,0);

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        PArray<MultiFab> particle_mf;
        DMPC->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            const IntVect ratio = parent->refRatio(lev);
            Gravity::average_down(particle_mf[lev], particle_mf[lev+1], ratio);
        }

        MultiFab::Copy(*derive_dat, particle_mf[level], 0, 0, 1, 0);

        MultiFab* gas_density = derive("density",time,0);

        MultiFab::Add(*derive_dat,*gas_density, 0, 0, 1, 0);

        delete gas_density;
        return derive_dat;
      }
      else
      {
        MultiFab* derive_dat = derive("density",time,0);
        return derive_dat;
      }
#else
        MultiFab* derive_dat = derive("density",time,0);
        return derive_dat;
#endif
    }
    else
    {
        return AmrLevel::derive(name, time, ngrow);
    }
}
