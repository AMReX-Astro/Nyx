#include <unistd.h>
#include <iomanip>
#include <Nyx.H>
#include <Nyx_F.H>
#include "Nyx_output.H"

#include "AMReX_buildInfo.H"

#ifdef FORCING
#include "Forcing.H"

#ifdef BL_HDF5
#include <vector>
using std::endl;
#endif


void mt_write(std::ofstream& output);
#endif

using namespace amrex;

namespace
{
    const std::string dm_chk_particle_file("DM");
    const std::string dm_plt_particle_file("DM");

    const std::string agn_chk_particle_file("AGN");
    const std::string agn_plt_particle_file("AGN");
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

void
Nyx::setPlotVariables ()
{
    AmrLevel::setPlotVariables();

    ParmParse pp("nyx");
    bool plot_X, plot_rank;
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
Nyx::writePlotFile (const std::string& dir,
                    ostream&           os,
                    VisMF::How         how)
{

#ifdef BL_HDF5
  std::string dirNoTemp(dir.substr(0, dir.length() - 5));
  writePlotFileHDF5(dirNoTemp, os, how);
  return;
#endif

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
            if (it->name() == "particle_count" ||
                it->name() == "total_particle_count" ||
                it->name() == "particle_mass_density" ||
                it->name() == "total_density")
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
                       it->name() == "neutrino_mass_density")
            {
                if (Nyx::theNPC())
                {
                    derive_names.push_back(it->name());
                    num_derive++;
                }
#endif
            } else if (it->name() == "Rank") {
                derive_names.push_back(it->name());
                num_derive++;
            } else {
                derive_names.push_back(it->name());
                num_derive++;
            }
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0) {
            amrex::Error("Must specify at least one valid data item to plot");
	}

        os << n_data_items << '\n';
        //
        // Names of variables -- first state, then derived
        //
        for (i = 0; i < plot_var_map.size(); i++)
        {
            int typ = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            os << desc_lst[typ].name(comp) << '\n';
        }

        for (std::list<std::string>::iterator it = derive_names.begin();
             it != derive_names.end(); ++it)
        {
            const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geom().Coord() << '\n';
        os << "0\n"; // Write bndry data.

        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	std::string PrettyLine = std::string(78, '=') + "\n";
	std::string OtherLine = std::string(78, '-') + "\n";
	std::string SkipSpace = std::string(8, ' ') + "\n";

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

        for (i = 0; i <= f_lev; i++)
          {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (n = 0; n < BL_SPACEDIM; n++)
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
        // Castro_setup.cpp
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
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if( ! levelDirectoryCreated) {
      amrex::Print() << "IOIOIOIO:CD  Nyx::writePlotFile:  ! ldc:  creating directory:  "
                     << FullPath << '\n';
      if (ParallelDescriptor::IOProcessor()) {
        if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
	}
      }
      //
      // Force other processors to wait until directory is built.
      //
      ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int cnt = 0;
    const int nGrow = 0;
    MultiFab plotMF(grids, dmap, n_data_items, nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
        cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (derive_names.size() > 0)
    {
        for (std::list<std::string>::iterator it = derive_names.begin();
             it != derive_names.end(); ++it)
        {
            const auto& derive_dat = derive(*it, cur_time, nGrow);
            MultiFab::Copy(plotMF, *derive_dat, 0, cnt, 1, nGrow);
            cnt++;
        }
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF, TheFullPath, how, true);

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

}


#ifdef BL_HDF5

template <class T>
void VWriteData(T &vData, hid_t vLevelGroup, long H5Ttype, const std::string &vName) {
  hid_t aid = H5Screate(H5S_SCALAR);
  hid_t attr = H5Acreate2(vLevelGroup, vName.c_str(), H5Ttype, aid, H5P_DEFAULT, H5P_DEFAULT);
  if(attr < 0) {
    std::cerr << " Problem writing attribute " << vName.c_str() << std::endl;
  }
  H5Awrite(attr, H5Ttype, &vData);
  H5Sclose(aid);
  H5Aclose(attr);
}


herr_t VWriteToLocation(hid_t loc_id,
                        std::map<std::string, int>  &m_int,
                        std::map<std::string, Real> &m_real,
                        std::map<std::string, std::string> &m_string)
{
  H5E_auto_t efunc; void* edata;
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
  herr_t  ret;

#define INSERT2(Ttype, mapName, H5Ttype)                                   \
  for (std::map<std::string, Ttype>::const_iterator p = mapName.begin();        \
      p!= mapName.end(); ++p)                                             \
    {                                                                     \
      hid_t aid  = H5Screate(H5S_SCALAR);                                 \
      H5Eset_auto2(H5E_DEFAULT, NULL, NULL);                            \
      hid_t attr = H5Acreate2(loc_id, p->first.c_str(), H5Ttype, aid, H5P_DEFAULT, H5P_DEFAULT); \
      if (attr < 0) {                                                      \
          H5Adelete(loc_id, p->first.c_str());                            \
          attr = H5Acreate2(loc_id, p->first.c_str(), H5Ttype,             \
                            aid, H5P_DEFAULT, H5P_DEFAULT);               \
	  if (attr < 0) {                                                  \
            std::cerr << " Problem writing attribute " << p->first.c_str() << std::endl;  \
          }                                                             \
        }                                                                 \
      H5Eset_auto2(H5E_DEFAULT, efunc, edata);                          \
      Ttype tmp = p->second;                                              \
      ret = H5Awrite(attr, H5Ttype, &tmp);                                \
      if (ret < 0) return ret;                                             \
      H5Sclose(aid);                                                      \
      H5Aclose(attr);                                                     \
    }
  INSERT2(Real, m_real, H5T_NATIVE_DOUBLE);
  INSERT2(int, m_int, H5T_NATIVE_INT);

    // string is different, of course
    for (std::map<std::string, std::string>::const_iterator p = m_string.begin();
        p!= m_string.end(); ++p)
    {
      hid_t s_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(s_type, p->second.length()); //extra requirement for strings
      hid_t aid  = H5Screate(H5S_SCALAR);
      H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
      hid_t attr = H5Acreate2(loc_id, p->first.c_str(), s_type,
                              aid, H5P_DEFAULT, H5P_DEFAULT);
      if (attr < 0) {
          H5Adelete(loc_id, p->first.c_str());
          attr = H5Acreate2(loc_id, p->first.c_str(), s_type,
                            aid, H5P_DEFAULT,H5P_DEFAULT);
          if (attr < 0) {
              std::cerr << " Problem writing attribute " << p->first.c_str() << std::endl;
            }
        }
      H5Eset_auto2(H5E_DEFAULT, efunc, edata);
      char* tmp = (char*)p->second.c_str();
      ret = H5Awrite(attr, s_type, tmp);
      if (ret < 0) return ret;
      H5Sclose(aid);
      H5Aclose(attr);
      H5Tclose(s_type);
    }

    return 0;
}


void
Nyx::writePlotFileHDF5 (const std::string &dir,
                          std::ostream      &os,
                          VisMF::How        how)
{
    BL_PROFILE("HyperCLaw::writePlotFileHDF5");

    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());

    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    Vector<std::string> derive_namesV;
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ(0); typ < desc_lst.size(); ++typ) {
        for (int comp(0); comp < desc_lst[typ].nComp(); ++comp) {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
            {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
                derive_namesV.push_back(desc_lst[typ].name(comp));
            }
        }
    }

    int num_derive(0);
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
         it != dlist.end();
         ++it)
    {
        if (parent->isDerivePlotVar(it->name())) {
            derive_names.push_back(it->name());
            derive_namesV.push_back(it->name());
            ++num_derive;
        }
    }

    int num_diag(0);

    int nComp = plot_var_map.size() + num_derive + num_diag;
    int a_numLevels(parent->finestLevel() + 1);
    std::string filename(dir + ".hdf5");
    Real a_time(parent->cumTime());

    // ---- get the data
    int cnt(0);
    const int nGrow(0);
    MultiFab  plotMF(grids, dmap, nComp, nGrow);
    MultiFab *this_dat = nullptr;
    for(int i(0); i < plot_var_map.size(); ++i) {
      int typ  = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      this_dat = &state[typ].newData();
      MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
      ++cnt;
    }
    if(derive_names.size() > 0) {
      for (std::list<std::string>::iterator it = derive_names.begin();
           it != derive_names.end(); ++it)
      {
        auto derive_dat = derive(*it, a_time, nGrow);
        MultiFab::Copy(plotMF, *derive_dat, 0, cnt, 1, nGrow);
        ++cnt;
      }
    }

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime80(ParallelDescriptor::second());

  // ---- start hdf5 part
  herr_t  ret;
  IntVect iv1;
  Box b3(geom.Domain());

  int  b3int[2 * BL_SPACEDIM];
  for(int i(0); i < BL_SPACEDIM; ++i) {
    b3int[i] = b3.smallEnd(i);
    b3int[i + BL_SPACEDIM] = b3.bigEnd(i);
  }
  hid_t box_id = H5Tcreate (H5T_COMPOUND, sizeof(Box));

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime741(ParallelDescriptor::second());
  double dPlotFileTime742(dPlotFileTime741 - dPlotFileTime80);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime742);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_1_0 = " << dPlotFileTime742 << "  seconds." << std::endl;
  }

#if BL_SPACEDIM == 1
  amrex::Abort("writePlotFileHDF5 not implemented in 1d.");
#elif BL_SPACEDIM == 2
  amrex::Abort("writePlotFileHDF5 not implemented in 2d.");
#elif BL_SPACEDIM == 3
  H5Tinsert (box_id, "lo_i", b3int[0] + 0 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", b3int[0] + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_k", b3int[0] + 2 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", b3int[0] + 3 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", b3int[0] + 4 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_k", b3int[0] + 5 * sizeof(int), H5T_NATIVE_INT);
#endif

  // ASim@lbl.gov 6/15/2016
  double dPlotFileTime711(ParallelDescriptor::second());
  double dPlotFileTime712(dPlotFileTime711 - dPlotFileTime741);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime712);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_1_1 = " << dPlotFileTime712 << "  seconds." << std::endl;
  }

  std::string vGroupName = "/";
  hid_t vFile;

  if(level == 0) {
    std::string filedescriptor("VanillaAMRFileType");
    hid_t file_access = H5Pcreate (H5P_FILE_ACCESS);
    // ASim@lbl.gov CollectiveMetaData
	// commented out 10/3/2017
    ret = H5Pset_all_coll_metadata_ops(file_access, true);
    ret = H5Pset_coll_metadata_write(file_access, true);

    H5Pset_fapl_mpio(file_access,  ParallelDescriptor::Communicator(), MPI_INFO_NULL);

    vFile = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_access);
    std::string vGlobalGroupName = "Chombo_global";
    H5Pclose(file_access);
    hid_t vCurrentGroup = H5Gopen2(vFile, vGroupName.c_str(),H5P_DEFAULT);

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime713(ParallelDescriptor::second());
    double dPlotFileTime714(dPlotFileTime713 - dPlotFileTime711);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime714);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write_H5M_time_1_2_0 = " << dPlotFileTime714 << "  seconds." << std::endl;
    }

    std::map<std::string, int>  vMInt;
    std::map<std::string, Real> vMReal;
    std::map<std::string, std::string> vMString;

    vMInt["SpaceDim"] = amrex::SpaceDim;
    vMReal["testReal"] = 0.0;
    vMString["testString"] = "vMString::testString";
    hid_t vGroup = H5Gcreate2(vFile, vGlobalGroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    ret = VWriteToLocation(vGroup, vMInt, vMReal, vMString);
    if(ret < 0) {
      std::cout << myProc << "**** Error 0:  ret = " << ret << std::endl;
    }
    H5Gclose(vGroup);

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime7131(ParallelDescriptor::second());
    double dPlotFileTime7141(dPlotFileTime713 - dPlotFileTime713);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime7141);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write_H5M_time_1_2 = " << dPlotFileTime7141 << "  seconds." << std::endl;
    }

    vMInt.clear();
    vMReal.clear();

    vMString ["filetype"]    = filedescriptor;
    vMInt ["num_levels"]     = a_numLevels;
    vMInt ["num_components"] = nComp;
    for(int ivar(0); ivar < nComp; ++ivar) {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      std::string label(labelChSt);
      vMString[label] = derive_namesV[ivar];
    }

    ret = VWriteToLocation(vCurrentGroup, vMInt, vMReal, vMString);
    if(ret < 0) {
      std::cout << myProc << "**** Error 1:  ret = " << ret << std::endl;
    }

    H5Gclose(vCurrentGroup);

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime715(ParallelDescriptor::second());
    double dPlotFileTime716(dPlotFileTime715 - dPlotFileTime7131);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime716);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write_H5M_time_1_3 = " << dPlotFileTime716 << "  seconds." << std::endl;
    }

  } else {
    hid_t file_access = H5Pcreate (H5P_FILE_ACCESS);
    // ASim@lbl.gov CollectiveMetaData
	// commented out 10/3/2017
    ret = H5Pset_all_coll_metadata_ops(file_access, true);
    ret = H5Pset_coll_metadata_write(file_access, true);
    H5Pset_fapl_mpio(file_access,  ParallelDescriptor::Communicator(), MPI_INFO_NULL);

    vFile = H5Fopen(filename.c_str(), H5F_ACC_RDWR, file_access);
    H5Pclose(file_access);

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime713(ParallelDescriptor::second());
    double dPlotFileTime714(dPlotFileTime713 - dPlotFileTime711);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime714);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write_H5M_time_1_4 = " << dPlotFileTime714 << "  seconds." << std::endl;
    }

  }

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime71(ParallelDescriptor::second());
  double dPlotFileTime72(dPlotFileTime71 - dPlotFileTime80);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime72);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_1 = " << dPlotFileTime72 << "  seconds." << std::endl;
  }

/*
  // ASim@lbl.gov 4/10/2017 metadata sizing test : get size
  // https://support.hdfgroup.org/HDF5/doc/RM/H5F/H5Fget_mdc_config.htm
  H5AC_cache_config_t config_ptr;
  config_ptr.version = H5AC__CURR_CACHE_CONFIG_VERSION;
  H5Fget_mdc_config(vFile, &config_ptr);
  if(ParallelDescriptor::IOProcessor()) {
      size_t max_size, min_clean_size, cur_size;
      int cur_num_entries;
      H5Fget_mdc_size(vFile, &max_size, &min_clean_size, &cur_size, &cur_num_entries);
      std::cout << "GET_MDC_SIZE = " << cur_size << std::endl;
      std::cout << "GET_MDC_SIZE_2 = " << max_size << " : " << min_clean_size << " : " << cur_num_entries << std::endl;
      std::cout << "GET_MDC_CONFIG = " << config_ptr.initial_size << std::endl;
  }
  // ASim@lbl.gov 4/11/2017 for setting larger meta cache
  // https://support.hdfgroup.org/HDF5/doc/RM/H5F/H5Fset_mdc_config.htm
  //config_ptr.set_initial_size = true;
  //config_ptr.initial_size = 24*1048576; // change the default 2097152
  config_ptr.evictions_enabled = false;
  config_ptr.incr_mode = H5C_incr__off;
  config_ptr.decr_mode = H5C_decr__off;
  config_ptr.flash_incr_mode = H5C_flash_incr__off;
  H5Fset_mdc_config(vFile, &config_ptr);
*/
/*
  H5AC_cache_config_t config_ptr;
  config_ptr.version = H5AC__CURR_CACHE_CONFIG_VERSION;
  config_ptr.rpt_fcn_enabled = false;
  config_ptr.set_initial_size = true;
  config_ptr.open_trace_file = false;
  config_ptr.close_trace_file = false;
  config_ptr.evictions_enabled = true;
  config_ptr.incr_mode = H5C_incr__threshold;
  config_ptr.decr_mode = H5C_decr__age_out;
  config_ptr.flash_incr_mode = H5C_flash_incr__add_space;
  config_ptr.increment = 2.0;
  config_ptr.decrement = 0.9;
  config_ptr.lower_hr_threshold = 0.9;
  config_ptr.upper_hr_threshold = 0.99995;
  config_ptr.apply_max_increment = false;
  config_ptr.apply_max_decrement = false;
  config_ptr.flash_multiple = 0.5;
  config_ptr.flash_threshold = 0.2;
  config_ptr.apply_empty_reserve = true;
  config_ptr.dirty_bytes_threshold = 524288;
  config_ptr.min_size = 1048576;
  config_ptr.max_size = 64*1048576;
  config_ptr.epoch_length = 262144;
  config_ptr.epochs_before_eviction = 3;

  config_ptr.initial_size = 24*1048576; // change the default 2097152
  H5Fset_mdc_config(vFile, &config_ptr);
*/

  char levelName[10];
  sprintf(levelName, "/level_%i", level);
  std::string gL(vGroupName + levelName);
  hid_t vLevelGroup = H5Gcreate2(vFile, gL.c_str(), H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

  std::string gLDA(gL + "/data_attributes");
  hid_t vLevelGroupDA = H5Gcreate2(vFile, gLDA.c_str(), H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
  VWriteData(nComp, vLevelGroupDA, H5T_NATIVE_INT, "comps");

  IntVect giv;
  int gint[BL_SPACEDIM];
  for(int gi(0); gi < BL_SPACEDIM; ++gi) {
    giv[gi] = 0;
    gint[gi] = 0;
  }
  hid_t gintvect_id = H5Tcreate (H5T_COMPOUND, sizeof(IntVect));
  H5Tinsert (gintvect_id, "intvecti", gint[0] + 0 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (gintvect_id, "intvectj", gint[0] + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (gintvect_id, "intvectk", gint[0] + 2 * sizeof(int), H5T_NATIVE_INT);
  VWriteData(gintvect_id, vLevelGroupDA, gintvect_id, "ghost");
  // the following does not seem to work
  //VWriteData(gintvect_id, vLevelGroupDA, gintvect_id, "outputGhost");


  Real a_dt(parent->dtLevel(level));
  const Real *a_dx = geom.CellSize();
  Real vData(a_dt);
  std::string vName("dt");
  long H5Ttype(H5T_NATIVE_DOUBLE);

  VWriteData(vData, vLevelGroup, H5Ttype, vName);
  VWriteData(a_dx[level], vLevelGroup, H5T_NATIVE_DOUBLE, "dx");
  VWriteData(a_time, vLevelGroup, H5T_NATIVE_DOUBLE, "time");
  VWriteData(b3, vLevelGroup, box_id, "prob_domain");

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime73(ParallelDescriptor::second());
  double dPlotFileTime74(dPlotFileTime73 - dPlotFileTime71);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime74);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_2 = " << dPlotFileTime74 << "  seconds." << std::endl;
  }


  // ---- "boxes" and "Processors" data
  Vector<int> procMap = plotMF.DistributionMap().ProcessorMap();
  hid_t procdataset, procdataspace;
  hid_t boxdataset, boxdataspace;
  hid_t offsetdataset, offsetdataspace;
  std::string pdsname("Processors");
  std::string bdsname("boxes");
  std::string odsname("data:offsets=0");
  hsize_t  flatdims[1], count[1], ocount[1];
  flatdims[0] = grids.size();
  H5E_auto_t efunc; void *edata; // turn auto error messaging off
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(vLevelGroup, pdsname.c_str(), H5P_DEFAULT);  // removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
  procdataspace = H5Screate_simple(1, flatdims, NULL);
  procdataset   = H5Dcreate2(vLevelGroup, pdsname.c_str(),  H5T_NATIVE_INT,
                             procdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  flatdims[0] = grids.size();

  boxdataspace = H5Screate_simple(1, flatdims, NULL);

  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(vLevelGroup, bdsname.c_str(),H5P_DEFAULT);  // removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);

  int bab3int[2 * BL_SPACEDIM];
  for(int i(0); i < BL_SPACEDIM; ++i) {
    bab3int[i] = 0;
    bab3int[i + BL_SPACEDIM] = 1;
  }
  int  boxSize(2 * BL_SPACEDIM);
  hid_t babox_id;
  babox_id = H5Tcreate (H5T_COMPOUND, boxSize * sizeof(int));
  H5Tinsert (babox_id, "lo_i", bab3int[0] + 0 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (babox_id, "lo_j", bab3int[0] + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (babox_id, "lo_k", bab3int[0] + 2 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (babox_id, "hi_i", bab3int[0] + 3 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (babox_id, "hi_j", bab3int[0] + 4 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (babox_id, "hi_k", bab3int[0] + 5 * sizeof(int), H5T_NATIVE_INT);

  boxdataset = H5Dcreate2(vLevelGroup, bdsname.c_str(),
                          babox_id, boxdataspace, H5P_DEFAULT,
                          H5P_DEFAULT,H5P_DEFAULT);

  int iRefRatio(1);
  if(level < parent->finestLevel()) {
    iRefRatio = parent->refRatio(level)[0];
  }
  VWriteData(iRefRatio, vLevelGroup, H5T_NATIVE_INT, "ref_ratio");

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime75(ParallelDescriptor::second());
  double dPlotFileTime76(dPlotFileTime75 - dPlotFileTime73);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime76);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_3 = " << dPlotFileTime76 << "  seconds." << std::endl;
  }

  // ---- create a boxarray sorted by rank
  std::map<int, Vector<Box> > gridMap;
  for(int i(0); i < grids.size(); ++i) {
    int gridProc(procMap[i]);
    Vector<Box> &boxesAtProc = gridMap[gridProc];
    boxesAtProc.push_back(grids[i]);
  }
  BoxArray sortedGrids(grids.size());
  Vector<int> sortedProcs(grids.size());
  int bIndex(0);
  for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
    int proc = it->first;
    Vector<Box> &boxesAtProc = it->second;
    for(int ii(0); ii < boxesAtProc.size(); ++ii) {
      sortedGrids.set(bIndex, boxesAtProc[ii]);
      sortedProcs[bIndex] = proc;
      ++bIndex;
    }
  }

  hsize_t  oflatdims[1];
  oflatdims[0] = sortedGrids.size() + 1;
  offsetdataspace = H5Screate_simple(1, oflatdims, NULL);
  offsetdataset   = H5Dcreate2(vLevelGroup, odsname.c_str(),  H5T_NATIVE_LLONG,
                                offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  Vector<unsigned long long> offsets(sortedGrids.size() + 1);
  unsigned long long currentOffset(0L);
  ocount[0] = sortedGrids.size() + 1;
  hid_t omemdataspace = H5Screate_simple(1, ocount, NULL);
  for(int b(0); b < sortedGrids.size(); ++b) {
    offsets[b] = currentOffset;
    currentOffset += sortedGrids[b].numPts() * nComp;
  }
  offsets[sortedGrids.size()] = currentOffset;

  Vector<unsigned long long> procOffsets(nProcs);
  int posCount(0);
  Vector<unsigned long long> procBufferSize(nProcs);
  unsigned long long totalOffset(0);
  for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
    int proc = it->first;
    Vector<Box> &boxesAtProc = it->second;
    BL_ASSERT(posCount == proc);
    procOffsets[posCount] = totalOffset;
    ++posCount;
    procBufferSize[proc] = 0L;
    for(int b(0); b < boxesAtProc.size(); ++b) {
      procBufferSize[proc] += boxesAtProc[b].numPts() * nComp;
    }
    totalOffset += procBufferSize[proc];
  }

  if(ParallelDescriptor::IOProcessor()) {
    int vbCount(0);
    Vector<int> vbox(sortedGrids.size() * boxSize);
    Vector<int> pid(sortedGrids.size());
    count[0] = sortedGrids.size();
    hid_t bmemdataspace = H5Screate_simple(1, count, NULL);
    hid_t pmemdataspace = H5Screate_simple(1, count, NULL);
    for(int b(0); b < sortedGrids.size(); ++b) {
      for(int i(0); i < BL_SPACEDIM; ++i) {
        vbox[(vbCount * boxSize) + i] = sortedGrids[b].smallEnd(i);
        vbox[(vbCount * boxSize) + i + BL_SPACEDIM] = sortedGrids[b].bigEnd(i);
      }
      ++vbCount;
      pid[b] = sortedProcs[b];
    }
    // ASim@lbl.gov CollectiveMetaData
	// commented out 10/3/2017
    //ret = H5Pset_all_coll_metadata_ops(bmemdataspace, true);
    //ret = H5Pset_coll_metadata_write(bmemdataspace, true);
    //ret = H5Pset_all_coll_metadata_ops(pmemdataspace, true);
    //ret = H5Pset_coll_metadata_write(pmemdataspace, true);


/*
    // ASim@lbl.gov 03/20/2017 for collective io setting: H5FD_MPIO_COLLECTIVE
    // Collective IO does not work here. H5FD_MPIO_INDEPENDENT = H5P_DEFAULT 
    // Is there a reason why this array needs to be written 
    // on this particular IOProcessor?
    // leaving along dxfer_template = H5P_DEFAULT;
    hid_t dxfer_template;
    dxfer_template = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
    H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
*/

    if(vbox.size() > 0) {
	hid_t dxfer_template;
	dxfer_template = H5P_DEFAULT;
        ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, omemdataspace, offsetdataspace,
                 dxfer_template, &(offsets[0]));
	if(ret < 0) { std::cout << "_here 0:  ret = " << ret << std::endl; }
        ret = H5Dwrite(boxdataset, babox_id, bmemdataspace, boxdataspace,
                 dxfer_template, &(vbox[0]));
	if(ret < 0) { std::cout << "_here 1:  ret = " << ret << std::endl; }
        ret = H5Dwrite(procdataset, H5T_NATIVE_INT, pmemdataspace, procdataspace,
                 dxfer_template, &(pid[0]));
	if(ret < 0) { std::cout << "_here 2:  ret = " << ret << std::endl; }
    } else {
	hid_t dxfer_template;
	dxfer_template = H5P_DEFAULT;
        ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, omemdataspace, offsetdataspace,
                 dxfer_template, NULL);
	if(ret < 0) { std::cout << "_here 3:  ret = " << ret << std::endl; }
        ret = H5Dwrite(boxdataset, babox_id, bmemdataspace, boxdataspace,
                 dxfer_template, NULL);
	if(ret < 0) { std::cout << "_here 4:  ret = " << ret << std::endl; }
        ret = H5Dwrite(procdataset, H5T_NATIVE_INT, pmemdataspace, procdataspace,
                 dxfer_template, NULL);
	if(ret < 0) { std::cout << "_here 5:  ret = " << ret << std::endl; }
    }

/*
    // ASim@lbl.gov 03/20/2017 for closing collective io
    H5Pclose(dxfer_template);
*/


    H5Sclose(bmemdataspace);
    H5Sclose(pmemdataspace);
  }

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime77(ParallelDescriptor::second());
  double dPlotFileTime78(dPlotFileTime77 - dPlotFileTime75);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime78);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_4 = " << dPlotFileTime78 << "  seconds." << std::endl;
  }

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime85(ParallelDescriptor::second());
  double dPlotFileTime86(dPlotFileTime85 - dPlotFileTime80);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime86);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5MD_time = " << dPlotFileTime86 << "  seconds." << std::endl;
  }


  {  // ---- data write
    BL_PROFILE_VAR("H5Dwritedata", h5dwd);
    hsize_t hs_procsize[1], hs_allprocsize[1], ch_offset[1];

    ch_offset[0]      = procOffsets[myProc];          // ---- offset on this proc
    hs_procsize[0]    = procBufferSize[myProc];       // ---- size of buffer on this proc
    hs_allprocsize[0] = offsets[sortedGrids.size()];  // ---- size of buffer on all procs

    char dataname[1024];
    sprintf(dataname, "data:datatype=0");

    hid_t dataspace    = H5Screate_simple(1, hs_allprocsize, NULL);
    hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

    hid_t dataset = H5Dcreate(vLevelGroup, dataname, H5T_NATIVE_DOUBLE, dataspace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // ASim@lbl.gov CollectiveMetaData
	// commented out 10/3/2017
    //ret = H5Pset_all_coll_metadata_ops(dataset, true);
    //ret = H5Pset_coll_metadata_write(dataset, true);

    //select where in the file it will be written
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, ch_offset, NULL,
                        hs_procsize, NULL);

    Vector<Real> a_buffer(procBufferSize[myProc], -1.0);
    long dataCount(0);
    for(MFIter mfi(plotMF); mfi.isValid(); ++mfi) {
      const Box &vbox    = mfi.validbox();
      const Real *dataPtr = plotMF[mfi].dataPtr();
      for(int i(0); i < vbox.numPts() * nComp; ++i) {
        a_buffer[dataCount++] = dataPtr[i];
      }
    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "::---- calling H5Dwrite for the grid data on level " << level << std::endl;
    }

/*
    // ASim@lbl.gov 4/10/2017 metadata sizing test 
    if(ParallelDescriptor::IOProcessor()) {
        // H5Fget_mdc_size(hid_t file_id, size_t *max_size_ptr, size_t *min_clean_size_ptr, size_t *cur_size_ptr, int *cur_num_entries_ptr)
        size_t max_size, min_clean_size, cur_size;
        int cur_num_entries;
        H5Fget_mdc_size(vFile, &max_size, &min_clean_size, &cur_size, &cur_num_entries);
        std::cout << "GET_MDC_SIZE = " << cur_size << std::endl;
        std::cout << "GET_MDC_SIZE_2 = " << max_size << " : " << min_clean_size << " : " << cur_num_entries << std::endl;
    }
*/

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime771(ParallelDescriptor::second());
  double dPlotFileTime781(dPlotFileTime771 - dPlotFileTime77);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime781);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_5 = " << dPlotFileTime781 << "  seconds." << std::endl;
  }

    hid_t dxfer_template;
    dxfer_template = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5INDEP
    ret = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
#else
    ret = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
#endif

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime00(ParallelDescriptor::second());

    BL_PROFILE_VAR("H5DwriteGrids", h5dwg);
    ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxfer_template,  a_buffer.dataPtr());
    BL_PROFILE_VAR_STOP(h5dwg);
    if(ret < 0) {
      std::cout << ParallelDescriptor::MyProc() << "_here 6:  ret = " << ret << std::endl;
    }

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime11(ParallelDescriptor::second());
    double dPlotFileTime22(dPlotFileTime11 - dPlotFileTime00);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime22);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write_H5Dwrite_time = " << dPlotFileTime22 << "  seconds." << std::endl;
      std::cout << "Write_H5Dwrite_time_since = " << ParallelDescriptor::second() << std::endl;
    }

	// ASim@lbl.gov 6/15/2017 for closing collective io
    H5Pclose(dxfer_template);

    H5Sclose(memdataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    BL_PROFILE_VAR_STOP(h5dwd);
  }

  H5Sclose(omemdataspace);
  H5Sclose(offsetdataspace);
  H5Dclose(offsetdataset);
  H5Sclose(boxdataspace);
  H5Dclose(boxdataset);
  H5Sclose(procdataspace);
  H5Dclose(procdataset);

  H5Gclose(vLevelGroupDA);
  H5Gclose(vLevelGroup);

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime791(ParallelDescriptor::second());
  H5Fclose(vFile);

  // ASim@lbl.gov 6/15/2017
  double dPlotFileTime81(ParallelDescriptor::second());
  double dPlotFileTime82(dPlotFileTime81 - dPlotFileTime80);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime82);
  double dPlotFileTime792(dPlotFileTime81 - dPlotFileTime791);
  ParallelDescriptor::ReduceRealMax(dPlotFileTime791);
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Write_H5M_time_7_closing = " << dPlotFileTime792 << "  seconds." << std::endl;
    std::cout << "Write_HDF5_time = " << dPlotFileTime82 << "  seconds." << std::endl;
  }


}
#endif





void
Nyx::writePlotFilePre (const std::string& dir, ostream& os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePre();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePre();
  }
#endif

}


void
Nyx::writePlotFilePost (const std::string& dir, ostream& os)
{
  if(Nyx::theDMPC()) {
    Nyx::theDMPC()->WritePlotFilePost();
  }
#ifdef AGN
  if(Nyx::theAPC()) {
    Nyx::theAPC()->WritePlotFilePost();
  }
#endif
}


void
Nyx::particle_plot_file (const std::string& dir)
{
    if (level == 0)
    {
        if (Nyx::theDMPC())
          {
            Nyx::theDMPC()->WriteNyxPlotFile(dir, dm_plt_particle_file);
          }

#ifdef AGN
        if (Nyx::theAPC())
          {
            Nyx::theAPC()->WriteNyxPlotFile(dir, agn_plt_particle_file);
          }
#endif

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
    }
}

void
Nyx::particle_check_point (const std::string& dir)
{
  BL_PROFILE("Nyx::particle_check_point");
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
  AmrLevel::checkPoint(dir, os, how, dump_old);
  particle_check_point(dir);
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
