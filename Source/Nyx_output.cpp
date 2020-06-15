#include <unistd.h>
#include <iomanip>
#include <Nyx.H>
#include <Nyx_F.H>
#include "Nyx_output.H"

#include "AMReX_buildInfo.H"

#ifdef FORCING
#include "Forcing.H"

void mt_write(std::ofstream& output);
#endif

using namespace amrex;

namespace
{
    const std::string dm_chk_particle_file("DM");
    const std::string dm_plt_particle_file("DM");

    const std::string agn_chk_particle_file("AGN");
    const std::string agn_plt_particle_file("AGN");

    const std::string npc_chk_particle_file("NPC");
    const std::string npc_plt_particle_file("NPC");

}

std::string
Nyx::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");
    return the_plot_file_type;
}

std::string
Nyx::retrieveDM () {
    return dm_chk_particle_file;
}

#ifdef AGN
std::string
Nyx::retrieveAGN () {
    return agn_chk_particle_file;
}
#endif

#ifdef NEUTRINO_DARK_PARTICLES
std::string
Nyx::retrieveNPC () {
    return npc_chk_particle_file;
}
#endif

void
Nyx::setPlotVariables ()
{
    AmrLevel::setPlotVariables();

    ParmParse pp("nyx");
    bool plot_X = false;
    bool plot_rank = false;
    if (pp.query("plot_rank", plot_rank))
    {
        if (plot_rank)
        {
            //
            // Write the processor ID for each grid into the plotfile
            //
            std::string proc_string = "Rank";
            parent->addDerivePlotVar(proc_string);
        }
    }
    if (pp.query("plot_X", plot_X))
    {
        if (plot_X)
        {
            //
            // Get the number of species from the network model.
            //
            fort_get_num_spec(&NumSpec);
            //
            // Get the species names from the network model.
            //
            for (int i = 0; i < NumSpec; i++)
            {
                int len = 20;
                Vector<int> int_spec_names(len);
                //
                // This call return the actual length of each string in "len"
                //
                fort_get_spec_names(int_spec_names.dataPtr(), &i, &len);
                char* spec_name = new char[len+1];

                for (int j = 0; j < len; j++)
                    spec_name[j] = int_spec_names[j];
                spec_name[len] = '\0';

                // @todo: better string ops
                std::string spec_string = "X(";
                spec_string += spec_name;
                spec_string += ')';
                parent->addDerivePlotVar(spec_string);
                delete [] spec_name;
            }
        }
    }
}

void
Nyx::writePlotFilePre (const std::string& dir, ostream& os)
{
    amrex::Gpu::LaunchSafeGuard lsg(true);
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePre();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePre();
  }
#endif

#ifdef NEUTRINO_DARK_PARTICLES
  if(Nyx::theNPC()) {
    Nyx::theNPC()->WritePlotFilePre();
  }
#endif

}

void
Nyx::writePlotFile (const std::string& dir,
                    ostream&           os,
                    VisMF::How         how)
{

    amrex::Gpu::LaunchSafeGuard lsg(true);
    AmrLevel::writePlotFile(dir, os, how);

}

void
Nyx::writePlotFilePost (const std::string& dir, ostream& os)
{

  amrex::Gpu::LaunchSafeGuard lsg(true);

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

        // Write comoving_a into its own file in the particle directory
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
            }
            File.precision(15);
            if (cur_time == 0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
            File.close();
        }

        if (ParallelDescriptor::IOProcessor() && use_typical_steps)
        {
            std::string FileName = dir + "/first_max_steps";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
            }
            File.precision(15);
            File << old_max_sundials_steps << '\n';
        }

        if (ParallelDescriptor::IOProcessor() && use_typical_steps)
        {
            std::string FileName = dir + "/second_max_steps";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
            }
            File.precision(15);
            File << old_max_sundials_steps << '\n';
        }

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      writeJobInfo(dir);
    }

    //
    // Write the particles and `comoving_a` in a plotfile directory. 
    //
    particle_plot_file(dir);

    // Write out all parameters into the plotfile
    if (write_parameters_in_plotfile) {
	write_parameter_file(dir);
    }

    if(Nyx::theDMPC()) {
      Nyx::theDMPC()->SetLevelDirectoriesCreated(false);
    }
#ifdef AGN
    if(Nyx::theAPC()) {
      Nyx::theAPC()->SetLevelDirectoriesCreated(false);
    }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
    if(Nyx::theNPC()) {
      Nyx::theNPC()->SetLevelDirectoriesCreated(false);
    }
#endif


  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePost();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePost();
  }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
  if(Nyx::theNPC()) {
    Nyx::theNPC()->WritePlotFilePost();
  }
#endif

  if(verbose) {

    if (level == 0)
    {
      if (cur_time == 0) {
	amrex::Print().SetPrecision(15) << "Output file " << dir << " at time " << std::to_string(old_a) << " and step " << std::to_string(nStep()) << std::endl;
      } else {
	amrex::Print().SetPrecision(15) << "Output file " << dir << " at time " << std::to_string(new_a) << " and step " << std::to_string(nStep()) << std::endl;
      }
    }
  }
}

void
Nyx::writeJobInfo (const std::string& dir)
{
        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	std::string PrettyLine = std::string(78, '=') + "\n";
	std::string OtherLine = std::string(78, '-') + "\n";
	std::string SkipSpace = std::string(8, ' ');

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Nyx Job Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "inputs file: " << inputs_name << "\n\n";

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
	jobInfoFile << "\n";
	jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
	  getCPUTime()/3600.0;

	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << "\n";
	}

	jobInfoFile << "\n\n";


	// cosmology information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Cosmology Information\n";
	jobInfoFile << PrettyLine;

        //	Real comoving_OmM, comoving_OmL, comoving_h;
        //	fort_get_omm(&comoving_OmM);
	// Omega lambda is defined algebraically
	Real comoving_OmL = 1. - comoving_OmM;

	// fort_get_hubble(&comoving_h);

	jobInfoFile << "Omega_m (comoving):      " << comoving_OmM << "\n";
	jobInfoFile << "Omega_lambda (comoving): " << comoving_OmL << "\n";
	jobInfoFile << "h (comoving):            " << comoving_h << "\n";

	jobInfoFile << "\n\n";

        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
	jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
	jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
	jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "Nyx    git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "AMReX git hash:  " << githash2 << "\n";
	}

	jobInfoFile << "\n\n";

	// grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        int f_lev = parent->finestLevel();

        for (int i = 0; i <= f_lev; i++)
          {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
              {
                jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
                //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
              }
            jobInfoFile << "\n\n";
          }

        jobInfoFile << " Boundary conditions\n";
        Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
        ParmParse pp("nyx");
        pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
        pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


        // these names correspond to the integer flags setup in the
        // Nyx_setup.cpp
        const char* names_bc[] =
          { "interior", "inflow", "outflow",
            "symmetry", "slipwall", "noslipwall" };


        jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
        jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
        if (BL_SPACEDIM >= 2) {
          jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
          jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
        }
        if (BL_SPACEDIM == 3) {
          jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
          jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
        }

        jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;

	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
}

void
Nyx::particle_plot_file (const std::string& dir)
{
    if (level == 0)
    {
        if (Nyx::theDMPC())
          {
            Nyx::theDMPC()->WriteNyxPlotFile(dir, dm_plt_particle_file);
	    ParmParse pp("particles");

	    int dm_particle_output_ascii = 0;
	    pp.query("dm_particle_output_ascii", dm_particle_output_ascii);
	    if (dm_particle_output_ascii != 0)
	    {
	      Nyx::theDMPC()->WriteAsciiFile(dir+"/particles.ascii");
	    }
          }

#ifdef AGN
        if (Nyx::theAPC())
          {
            Nyx::theAPC()->WriteNyxPlotFile(dir, agn_plt_particle_file);
          }
#endif

#ifdef NEUTRINO_DARK_PARTICLES
        if (Nyx::theNPC())
          {
            Nyx::theNPC()->WriteNyxPlotFile(dir, npc_plt_particle_file);
          }
#endif

#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif

        // Write particle_plotfile_format into its own file in the particle directory
        if (Nyx::theDMPC() && ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/" + dm_plt_particle_file + "/precision";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            File << particle_plotfile_format << '\n';
            File.close();
        }

#ifdef AGN
        // Write particle_plotfile_format into its own file in the particle directory
        if (Nyx::theAPC() && ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/" + agn_plt_particle_file + "/precision";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            File << particle_plotfile_format << '\n';
            File.close();
        }
#endif

#ifdef NEUTRINO_DARK_PARTICLES
        // Write particle_plotfile_format into its own file in the particle directory
        if (Nyx::theNPC() && ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/" + npc_plt_particle_file + "/precision";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            File << particle_plotfile_format << '\n';
            File.close();
        }
#endif
    }
}

void
Nyx::particle_check_point (const std::string& dir)
{
  BL_PROFILE("Nyx::particle_check_point");
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (level == 0)
    {
      if (Nyx::theDMPC())
        {
          Nyx::theDMPC()->NyxCheckpoint(dir, dm_chk_particle_file);
        }
#ifdef AGN
      if (Nyx::theAPC())
        {
          Nyx::theAPC()->NyxCheckpoint(dir, agn_chk_particle_file);
        }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
      if (Nyx::theNPC())
        {
          Nyx::theNPC()->NyxCheckpoint(dir, npc_chk_particle_file);
        }
#endif

#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif

        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/comoving_a";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            if (cur_time == 0)
            {
               File << old_a << '\n';
            } else {
               File << new_a << '\n';
            }
        }

	if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/first_max_steps";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
	    File << old_max_sundials_steps << '\n';
        }

	if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/second_max_steps";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
	    File << old_max_sundials_steps << '\n';
        }
    }
}

void
Nyx::write_parameter_file (const std::string& dir)
{
    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/the_parameters";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.precision(15);
            ParmParse::dumpTable(File,true);
            File.close();
        }
    }
}

void
Nyx::writeMultiFabAsPlotFile(const std::string& pltfile,
                             const MultiFab&    mf,
                             std::string        componentName)
{
    std::ofstream os;
    if (ParallelDescriptor::IOProcessor())
    {
        if( ! amrex::UtilCreateDirectory(pltfile, 0755)) {
          amrex::CreateDirectoryFailed(pltfile);
	}
        std::string HeaderFileName = pltfile + "/Header";
        os.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        // The first thing we write out is the plotfile type.
        os << thePlotFileType() << '\n';
        // Just one component ...
        os << 1 << '\n';
        // ... with name
        os << componentName << '\n';
        // Dimension
        os << BL_SPACEDIM << '\n';
        // Time
        os << "0\n";
        // One level
        os << "0\n";
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        // Only one level -> no refinement ratios
        os << '\n';
        // Geom
        os << parent->Geom(0).Domain() << ' ';
        os << '\n';
        os << parent->levelSteps(0) << ' ';
        os << '\n';
        for (int k = 0; k < BL_SPACEDIM; k++)
            os << parent->Geom(0).CellSize()[k] << ' ';
        os << '\n';
        os << (int) Geom().Coord() << '\n';
        os << "0\n"; // Write bndry data.
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    std::string Level = "Level_0";
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = pltfile;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
        if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
	}
    }
    //
    // Force other processors to wait until directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        std::string PathNameInHeader = Level;
        PathNameInHeader += BaseName;
        os << PathNameInHeader << '\n';
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(mf, TheFullPath);
    ParallelDescriptor::Barrier();
}

void
Nyx::checkPoint (const std::string& dir,
                 std::ostream&      os,
                 VisMF::How         how,
                 bool               dump_old_default)
{

  for (int s = 0; s < desc_lst.size(); ++s) {
      if (dump_old && state[s].hasOldData()) {
          MultiFab& old_MF = get_old_data(s);
          amrex::prefetchToHost(old_MF);
      }
      MultiFab& new_MF = get_new_data(s);
      amrex::prefetchToHost(new_MF);
  }

  AmrLevel::checkPoint(dir, os, how, dump_old);

  for (int s = 0; s < desc_lst.size(); ++s) {
      if (dump_old && state[s].hasOldData()) {
          MultiFab& old_MF = get_old_data(s);
          amrex::prefetchToDevice(old_MF);
      }
      MultiFab& new_MF = get_new_data(s);
      amrex::prefetchToDevice(new_MF);
  }

  particle_check_point(dir);

  if (level == 0 && ParallelDescriptor::IOProcessor())
  {
      writeJobInfo(dir);
  }

#ifdef FORCING
  forcing_check_point(dir);
#endif

  if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      {
	// store elapsed CPU time
	std::ofstream CPUFile;
	std::string FullPathCPUFile = dir;
	FullPathCPUFile += "/CPUtime";
	CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

	CPUFile << std::setprecision(15) << getCPUTime();
	CPUFile.close();
      }
    }

    if(Nyx::theDMPC()) {
      Nyx::theDMPC()->SetLevelDirectoriesCreated(false);
    }
#ifdef AGN
    if(Nyx::theAPC()) {
      Nyx::theAPC()->SetLevelDirectoriesCreated(false);
    }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
    if(Nyx::theNPC()) {
      Nyx::theNPC()->SetLevelDirectoriesCreated(false);
    }
#endif

}

void
Nyx::checkPointPre (const std::string& dir,
                    std::ostream&      os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->CheckpointPre();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->CheckpointPre();
  }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
  if(Nyx::theNPC()) {
    Nyx::theNPC()->CheckpointPre();
  }
#endif

}


void
Nyx::checkPointPost (const std::string& dir,
                 std::ostream&      os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->CheckpointPost();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->CheckpointPost();
  }
#endif
#ifdef NEUTRINO_DARK_PARTICLES
  if(Nyx::theNPC()) {
    Nyx::theNPC()->CheckpointPost();
  }
#endif
}


#ifdef FORCING
void
Nyx::forcing_check_point (const std::string& dir)
{
    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::string FileName = dir + "/forcing";
            std::ofstream File;
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            File.setf(std::ios::scientific, std::ios::floatfield);
            File.precision(16);
            forcing->write_Spectrum(File);
            File.close();

            FileName = dir + "/mt";
            File.open(FileName.c_str(), std::ios::out|std::ios::trunc);
            if ( ! File.good()) {
                amrex::FileOpenFailed(FileName);
	    }
            mt_write(File);
        }
    }
}
#endif

int
Nyx::updateInSitu ()
{
#if defined(BL_USE_SENSEI_INSITU) || defined(AMREX_USE_ASCENT)
    BL_PROFILE("Nyx::UpdateInSitu()");

#if defined(BL_USE_SENSEI_INSITU)
    if (insitu_bridge && insitu_bridge->update(this))
    {
        amrex::ErrorStream() << "Amr::updateInSitu : Failed to update." << std::endl;
        amrex::Abort();
    }
#endif

#ifdef AMREX_USE_CONDUIT
    blueprint_check_point();
#endif

#ifdef REEBER
    halo_find(parent->dtLevel(level));
#endif

#endif
    return 0;
}


#ifdef AMREX_USE_CONDUIT
void
Nyx::blueprint_check_point ()
{
  //    MultiFab& S_new = get_new_data(State_Type);
    //MultiFab S_new_tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0 NUM_GROW);

    const Real cur_time = state[State_Type].curTime();
    int cycle = nStep();

    Vector<std::string> varnames;

    varnames.push_back("Density");
    varnames.push_back("Xmom");
    varnames.push_back("Ymom");
    varnames.push_back("Zmom");
    varnames.push_back("Eden");
    varnames.push_back("Eint");

    if(!use_const_species)
    {
        varnames.push_back("H");
        varnames.push_back("He");
    }

    ///////////////////////////////////////////////////////////////////////////
    // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
    ///////////////////////////////////////////////////////////////////////////
    conduit::Node bp_mesh;
    SingleLevelToBlueprint(get_new_data(State_Type),
                           varnames,
                           geom,
                           cur_time,
                           cycle,
                           bp_mesh);

    //conduit::Node bp_particles;
    
    Vector<std::string> particle_varnames;
    particle_varnames.push_back("particle_mass");
    particle_varnames.push_back("particle_xvel");
    particle_varnames.push_back("particle_yvel");
    particle_varnames.push_back("particle_zvel");

    Vector<std::string> particle_int_varnames;
    amrex::ParticleContainerToBlueprint(*(Nyx::theDMPC()),
                                        particle_varnames,
                                        particle_int_varnames,
                                        bp_mesh);
    
    // very helpful for debugging when we actual try
    // to pull the varnames list from amrex, vs hand initing
    //
    // amrex::Print()<<varnames.size()<<S_new.nComp()<<std::endl;
    // amrex::Print()<<particle_varnames.size()<<4<<std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Uncomment below to: 
    // Save the Blueprint Mesh to a set of files that we can 
    // view in VisIt. 
    // (For debugging and to demonstrate how to do this w/o Ascent)
    ///////////////////////////////////////////////////////////////////////////
    //
    // amrex::Print()<< "Exporting Conduit Blueprint HDF5 files (cycle="
    //               << cycle <<")"
    //               << std::endl;
    //
    // WriteBlueprintFiles(bp_mesh,"bp_example_",cycle,"hdf5");

    ///////////////////////////////////////////////////////////////////
    // Render with Ascent
    ///////////////////////////////////////////////////////////////////

    if(verbose)
        amrex::Print()<< "Executing Ascent (cycle="
                      << cycle <<")"
                      << std::endl;

    Ascent ascent;
    conduit::Node open_opts;
    // tell ascent to use the ghost_indicator field to exclude ghosts

    open_opts["ghost_field_name"] = "ghost_indicator";
    
#ifdef BL_USE_MPI
    // if mpi, we need to provide the mpi comm to ascent
    open_opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

    ascent.open(open_opts);
    // publish structured mesh to ascent
    ascent.publish(bp_mesh);
    
    // call ascent, with empty actions.
    // actions below will be overridden by those in
    // ascent_actions.yaml
    Node actions;
    ascent.execute(actions);
    ascent.close();



}
#endif
