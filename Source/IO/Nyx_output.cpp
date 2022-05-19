#include <unistd.h>
#include <iomanip>
#include <Nyx.H>
#ifdef AMREX_USE_SUNDIALS
#include <sundials/sundials_version.h>
#endif
#include <AMReX_PlotFileUtil.H>
#include <Nyx_output.H>
#include <AMReX_buildInfo.H>
#include <Forcing.H>

void mt_write(std::ofstream& output);

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
}

void
Nyx::writePlotFilePre (const std::string& /*dir*/, ostream& /*os*/)
{
    if(write_hdf5 == 1 && (parent->maxLevel() > 0))
        amrex::Error("Calling single-level hdf5 interface for multilevel code (max_level > 0)");
    if(write_skip_prepost == 1)
    {
        amrex::Print()<<"Skip writePlotFilePre"<<std::endl;
    }
#ifdef AMREX_PARTICLES
    else
    {
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
#endif
}

void
Nyx::writePlotFile (const std::string& dir,
                    ostream&           os,
                    VisMF::How         how)
{

#ifdef AMREX_USE_HDF5
    if(write_hdf5==1 && parent->finestLevel() == 0)
    {
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    std::string dir_final = dir;
    if(!amrex::AsyncOut::UseAsyncOut())
    {
        auto start_position_to_erase = dir_final.find(".temp");
        dir_final.erase(start_position_to_erase,5);
    }

    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
    {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
        {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp))
                && desc_lst[typ].getType() == IndexType::TheCellType())
            {
                plot_var_map.push_back(std::pair<int,int>(typ, comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;

    for (std::list<DeriveRec>::const_iterator it = derive_lst.dlist().begin();
         it != derive_lst.dlist().end(); ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
        {
#ifdef AMREX_PARTICLES
            if (it->name() == "particle_count" ||
                it->name() == "total_particle_count" ||
                it->name() == "particle_mass_density" ||
                it->name() == "total_density" || 
                it->name() == "particle_x_velocity" ||
                it->name() == "particle_y_velocity" ||
                it->name() == "particle_z_velocity" )
            {
                if (Nyx::theDMPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#ifdef AGN
            } else if (it->name() == "agn_particle_count" ||
                       it->name() == "agn_mass_density")
            {
                if (Nyx::theAPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#endif
#ifdef NEUTRINO_PARTICLES
            } else if (it->name() == "neutrino_particle_count" ||
                       it->name() == "neutrino_mass_density" ||
                       it->name() == "neutrino_x_velocity" ||
                       it->name() == "neutrino_y_velocity" ||
                       it->name() == "neutrino_z_velocity" )
            {
                if (Nyx::theNPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#endif
            } else
#endif
                        if (it->name() == "Rank") {
                derive_names.push_back(it->name());
                num_derive++;
            } else {
                derive_names.push_back(it->name());
                num_derive++;
            }
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int cnt = 0;
    const int nGrow = 0;
    MultiFab plotMF(grids, dmap, n_data_items, nGrow);
    MultiFab* this_dat = 0;
    Vector<std::string> varnames;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        varnames.push_back(desc_lst[typ].name(comp));
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
        cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (derive_names.size() > 0)
    {
      for (std::list<DeriveRec>::const_iterator it = derive_lst.dlist().begin();
           it != derive_lst.dlist().end(); ++it)
        {
            varnames.push_back(it->name());
            const auto& derive_dat = derive(it->name(), cur_time, nGrow);
            MultiFab::Copy(plotMF, *derive_dat, 0, cnt, 1, nGrow);
            cnt++;
        }
    }

    WriteSingleLevelPlotfileHDF5(dir_final,
                          plotMF, varnames,
                          Geom(), cur_time, nStep());
//                          const std::string &versionName,
//                          const std::string &levelPrefix,
//                          const std::string &mfPrefix,
//                          const Vector<std::string>& extra_dirs)
    }
    else
#endif
        AmrLevel::writePlotFile(dir, os, how);

}

void
Nyx::writePlotFilePost (const std::string& dir, ostream& /*os*/)
{

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    if(write_skip_prepost == 1)
    {
        amrex::Print()<<"Skip writePlotFilePost"<<std::endl;
    }
    else
    {

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

#ifdef AMREX_PARTICLES
    //
    // Write the particles and `comoving_a` in a plotfile directory. 
    //
    particle_plot_file(dir);
#endif

    // Write out all parameters into the plotfile
    if (write_parameters_in_plotfile) {
        write_parameter_file(dir);
    }
#ifdef AMREX_PARTICLES
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
#endif
        }
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
Nyx::writeBuildInfo ()
{
  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // build information
  std::cout << PrettyLine;
  std::cout << " Nyx Build Information\n";
  std::cout << PrettyLine;

  std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
  std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
  std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
  std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  std::cout << "\n";

  std::cout << "COMP:          " << buildInfoGetComp() << "\n";
  std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  std::cout << "\n";

  std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  std::cout << "\n";

  std::cout << "Fortran comp:  " << buildInfoGetFName() << "\n";
  std::cout << "Fortran flags: " << buildInfoGetFFlags() << "\n";

  std::cout << "\n";

  std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

  std::cout << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  std::cout << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);
  int major, minor, patch;
  int len=256;
  char* label = new char[len];
#ifdef AMREX_USE_SUNDIALS
  SUNDIALSGetVersionNumber(&major, &minor, &patch, label, len);
#endif
#ifdef AMREX_USE_CONDUIT
    conduit::Node about;
    ascent::about(about);
    const std::string ascent_version = about["version"].to_json();
    const std::string ascent_githash = about["git_sha1"].to_json();
    //    about.print_detailed();
#endif
  if (strlen(githash1) > 0) {
    std::cout << "Nyx          git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    std::cout << "AMReX        git describe: " << githash2 << "\n";
  }
#ifdef AMREX_USE_SUNDIALS
  std::cout << "Sundials     git version:  " << SUNDIALS_GIT_VERSION << "\n";
  std::cout << "Sundials     int version:  v" <<major<<"."<<minor<<"."<<patch<< "\n";
#endif

#ifdef AMREX_USE_CONDUIT
  std::cout << "Ascent       git version:   " << ascent_githash << "\n";
  std::cout << "Ascent    string version:   " << ascent_version << "\n";
#endif
  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0){
    std::cout << buildgitname << " git describe: " << buildgithash << "\n";
  }

  std::cout << "\n\n";
  /*
// Possibly put this into job_info
#ifdef AMREX_USE_CONDUIT
  std::cout << PrettyLine;
  std::cout << " Ascent Build Information\n";
  std::cout << PrettyLine;

  std::cout << about.to_json();

  std::cout << "\n\n";
#endif
  */
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

        //      Real comoving_OmM, comoving_OmL, comoving_h;
        //      fort_get_omm(&comoving_OmM);
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
            for (int n = 0; n < AMREX_SPACEDIM; n++)
              {
                jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
                //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
              }
            jobInfoFile << "\n\n";
          }

        jobInfoFile << " Boundary conditions\n";
        Vector<int> lo_bc_out(AMREX_SPACEDIM), hi_bc_out(AMREX_SPACEDIM);
        ParmParse pp("nyx");
        pp.getarr("lo_bc",lo_bc_out,0,AMREX_SPACEDIM);
        pp.getarr("hi_bc",hi_bc_out,0,AMREX_SPACEDIM);


        // these names correspond to the integer flags setup in the
        // Nyx_setup.cpp
        const char* names_bc[] =
          { "interior", "inflow", "outflow",
            "symmetry", "slipwall", "noslipwall" };


        jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
        jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
        if (AMREX_SPACEDIM >= 2) {
          jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
          jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
        }
        if (AMREX_SPACEDIM == 3) {
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

#ifdef AMREX_PARTICLES
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
#endif
void
Nyx::particle_check_point (const std::string& dir)
{
  BL_PROFILE("Nyx::particle_check_point");

  if (level == 0)
    {
#ifdef AMREX_PARTICLES
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
#endif
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
        os << AMREX_SPACEDIM << '\n';
        // Time
        os << "0\n";
        // One level
        os << "0\n";
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        // Only one level -> no refinement ratios
        os << '\n';
        // Geom
        os << parent->Geom(0).Domain() << ' ';
        os << '\n';
        os << parent->levelSteps(0) << ' ';
        os << '\n';
        for (int k = 0; k < AMREX_SPACEDIM; k++)
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
            for (int n = 0; n < AMREX_SPACEDIM; n++)
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
                 bool               /*dump_old_default*/)
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

#ifndef NO_HYDRO
  if (do_forcing)
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
#ifdef AMREX_PARTICLES
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
#endif

}

void
Nyx::checkPointPre (const std::string& /*dir*/,
                    std::ostream&      /*os*/)
{
#ifdef AMREX_PARTICLES
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
#endif
}


void
Nyx::checkPointPost (const std::string& dir,
                     std::ostream&      /*os*/)
{
#ifdef AMREX_PARTICLES
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
#endif
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
  if(level==0)
    {
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

#ifndef NO_HYDRO
void
Nyx::forcing_check_point (const std::string& dir)
{
    if (level == 0 && do_forcing)
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
    amrex::Vector<Halo>& reeber_halos);
    halo_find(parent->dtLevel(level), reeber_halos);
    halo_print(reeber_halos);
#endif

#endif
    return 0;
}

#ifdef REEBER
void Nyx::halo_print(amrex::Vector<Halo>& reeber_halos)
{

       amrex::Real    halo_mass;
       amrex::IntVect halo_pos ;
       const auto dx = geom.CellSizeArray();

       for (const Halo& h : reeber_halos)
       {
           halo_mass = h.total_mass;
           halo_pos  = h.position;

           if (halo_mass > mass_halo_min)
           {
                amrex::Real x = (halo_pos[0]+0.5) * dx[0];
                amrex::Real y = (halo_pos[1]+0.5) * dx[1];
                amrex::Real z = (halo_pos[2]+0.5) * dx[2];

                amrex::Real mass = mass_seed;

                int lev = 0;
                int grid = 0;
                int tile = 0;

                amrex::AllPrintToFile("reeber_halos")<<h.id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<h.total_mass<<"\t"<<h.n_cells<<"\t"<<h.n_cells<<std::endl;

           }
       } // end of loop over creating new particles from halos

}
#endif

#ifdef AMREX_USE_CONDUIT
void
Nyx::blueprint_check_point ()
{
  //    MultiFab& S_new = get_new_data(State_Type);
    //MultiFab S_new_tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0 NUM_GROW);
#ifdef NO_HYDRO
    const Real cur_time = state[PhiGrav_Type].curTime();
#else
    const Real cur_time = state[State_Type].curTime();
#endif
    int cycle = nStep();

    conduit::Node bp_mesh;
#ifndef NO_HYDRO
    Vector<std::string> varnames;

    varnames.push_back("Density");
    varnames.push_back("Xmom");
    varnames.push_back("Ymom");
    varnames.push_back("Zmom");
    varnames.push_back("Eden");
    varnames.push_back("Eint");

#ifndef CONST_SPECIES
    varnames.push_back("H");
    varnames.push_back("He");
#endif

    ///////////////////////////////////////////////////////////////////////////
    // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
    ///////////////////////////////////////////////////////////////////////////
    SingleLevelToBlueprint(get_new_data(State_Type),
                           varnames,
                           geom,
                           cur_time,
                           cycle,
                           bp_mesh);
#else
    Vector<std::string> varnames;
    const int n_data_items = 5;
    const int nGrow = 0;
    MultiFab plotMF(grids, dmap, n_data_items, nGrow);
    varnames.push_back("particle_mass_density");
    varnames.push_back("particle_count");
    varnames.push_back("particle_x_velocity");
    varnames.push_back("particle_y_velocity");
    varnames.push_back("particle_z_velocity");

    for(int cnt = 0; cnt < n_data_items; cnt++)
    {
        const auto& derive_dat = derive(varnames[cnt], cur_time, nGrow);
        MultiFab::Copy(plotMF, *derive_dat, 0, cnt, 1, nGrow);
    }
    SingleLevelToBlueprint(plotMF,
                           varnames,
                           geom,
                           cur_time,
                           cycle,
                           bp_mesh);
#endif
    //conduit::Node bp_particles;
#ifdef AMREX_PARTICLES    
    Vector<std::string> particle_varnames;
    particle_varnames.push_back("particle_mass");
    particle_varnames.push_back("particle_xvel");
    particle_varnames.push_back("particle_yvel");
    particle_varnames.push_back("particle_zvel");

    Vector<std::string> particle_int_varnames;
    amrex::ParticleContainerToBlueprint(*(Nyx::theDMPC()),
                                        particle_varnames,
                                        particle_int_varnames,
                                        bp_mesh,dm_plt_particle_file);

#ifdef NEUTRINO_PARTICLES

#ifdef NEUTRINO_DARK_PARTICLES    
    Vector<std::string> neutrino_varnames;
    neutrino_varnames.push_back("neutrino_mass");
    neutrino_varnames.push_back("neutrino_xvel");
    neutrino_varnames.push_back("neutrino_yvel");
    neutrino_varnames.push_back("neutrino_zvel");

    Vector<std::string> neutrino_int_varnames;
    amrex::ParticleContainerToBlueprint(*(Nyx::theNPC()),
                                        neutrino_varnames,
                                        neutrino_int_varnames,
                                        bp_mesh,npc_plt_particle_file);
#endif
#endif
#endif
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
