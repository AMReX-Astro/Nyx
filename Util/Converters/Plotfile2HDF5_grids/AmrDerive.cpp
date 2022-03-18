/*

AMReX plotfile conversion to HDF5.  This converts only grid data.  There is a separate
code (in python) for particle data conversion.

The main routine first loads everything with amrex routines. It opens the HDF5
file, and adds the required groups and metadata in prepare_output. It loads
each field into the field_buffer array in c-order and writes it to the file
in write_field.

*/

#include <stdint.h>
#include <string>

// amrex includes
#include <new>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"
// end amrex includes

#include <hdf5.h>

#ifdef MEM_CHECK
#define CHK(x) x
#else
#define CHK(x)
#endif

using namespace amrex;


static void print_usage(int argc, char **argv) {
    printf("Usage: %s input_path=... output_path=...\n", argv[0]);
    printf("  Converts an AMReX plotfile to Lyman alpha HDF5 format.\n");
    exit(1);
}

bool dmvars_on(std::string& fns) {
    // Looks for nyx.do_hydro = 1 in the jobinfo file, if it is there, Hydro is on
    std::ifstream fileInput;
    std::string line;
    bool dmvars = false;
    char* search = "particle_x_velocity"; // test variable to search in file
    // open file to search
    fileInput.open(fns.c_str());
    while(getline(fileInput, line)) {
        if (line.find(search, 0) != std::string::npos) {
            dmvars = true;
        }
    }
    return dmvars;
}

bool neutrinos_on(std::string& fns) {
    // Looks for the directory NPC in the plotfile, if this exists neutrinos are on
    struct stat sbuf;
    const char *fn = fns.c_str();
    if( lstat( fn, &sbuf ) == -1 ) {
        Print()<< "Neutrinos are off\n" << endl;
        return false;
    } else if( S_ISDIR( sbuf.st_mode ) ) {
        Print()<<"Neutrinos are on \n" << endl;
        return true;
    } else {
        Print() << "What?\n" << endl;;
        return false;
    }
}


bool hydro_on(std::string& fns) {
    // Looks for nyx.do_hydro = 1 in the jobinfo file, if it is there, Hydro is on
    std::ifstream fileInput;
    std::string line;
    bool hydro = false;
    char* search = "nyx.do_hydro = 1"; // test variable to search in file
    // open file to search
    fileInput.open(fns.c_str());
    while(getline(fileInput, line)) {
        if (line.find(search, 0) != std::string::npos) {
            hydro = true;
        }   
    }
    return hydro;
}


bool append_on(std::string& fns) {
    // Looks for nyx.do_hydro = 1 in the jobinfo file, if it is there, Hydro is on
    std::ifstream fileInput;
    std::string line;
    bool append = false;
    char* search = "nyx.particle_init_type = Restart"; // test variable to search in file
    // open file to search
    fileInput.open(fns.c_str());
    while(getline(fileInput, line)) {
        if (line.find(search, 0) != std::string::npos) {
            append = true;
        }
    }
    return append;
}


double parseByName(std::string file_name, std::string var_name) {
    std::ifstream fin;
    std::string line = "";

    fin.open(file_name.c_str());
    while(std::getline(fin, line)) {
        if(line.find("//") == 0 || line.empty())
            continue;
        else {
            int pos = line.find("=");
            std::string name   = line.substr(0,pos-1);
            std::string svalue = line.substr(pos+1);
            if(name == var_name) {
                std::istringstream instr(svalue);
                double value;
                instr >> value;
                return value;
            }
        }
    }
    return -9999.0;
}


/*
Creates the HDF5 file in truncate mode and closes it.
Should be run only by the root process.
*/
void output_create(const std::string file_path) {
    hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        printf("Error: could not create file at %s.", file_path.c_str());
        exit(1);
    }
    H5Fclose(file);
}


/*
Opens the output file and writes all of metadata attributes.
Should be run only by the root process.
*/
void output_write_metadata(const std::string file_path,
 const int nx, const int ny, const int nz, const double domain_size,
 const double omega_b, const double omega_m, const double omega_l,
 const double h, const double z) {

    hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    //
    // Output format root attribute
    //

    const std::string format_text = "nyx-lyaf";

    hid_t str_type     = H5Tcopy(H5T_C_S1);
    hid_t scalar_space = H5Screate(H5S_SCALAR);

    // Fix the str_type length for the format string.
    H5Tset_size(str_type, strlen(format_text.c_str()));

    hid_t attr = H5Acreate(file, "format", str_type, scalar_space, H5P_DEFAULT,
       H5P_DEFAULT);
    H5Awrite(attr, str_type, format_text.c_str());

    //
    // Domain metadata group
    //

    hid_t group = H5Gcreate(file, "domain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Dataspace for domain attributes
    hsize_t num_dims[1] = {3};
    hid_t grid_attr_space = H5Screate_simple(1, num_dims, NULL);

    // Grid shape
    int shape[3] = {nx, ny, nz};
    attr = H5Acreate(group, "shape", H5T_STD_I32LE, grid_attr_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, shape);

    // Grid size
    double size[3] = {domain_size, domain_size, domain_size};
    attr = H5Acreate(group, "size", H5T_IEEE_F64LE, grid_attr_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, size);

    //
    // Universe metadata group
    //

    group = H5Gcreate(file, "universe", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);

    attr = H5Acreate(group, "omega_b", H5T_IEEE_F64LE, scalar_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_b);
    attr = H5Acreate(group, "omega_m", H5T_IEEE_F64LE, scalar_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_m);
    attr = H5Acreate(group, "omega_l", H5T_IEEE_F64LE, scalar_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_l);
    attr = H5Acreate(group, "hubble", H5T_IEEE_F64LE, scalar_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &h);
    attr = H5Acreate(group, "redshift", H5T_IEEE_F64LE, scalar_space,
     H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &z);

    //
    // Field groups
    //

    group = H5Gcreate(file, "native_fields", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);

    group = H5Gcreate(file, "derived_fields", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);

    // Close all resources.
    H5Sclose(grid_attr_space);
    H5Gclose(group);
    H5Aclose(attr);
    H5Sclose(scalar_space);
    H5Tclose(str_type);
    H5Fclose(file);
    H5close();
}


/*
Creates a dataset with the given cell dimensions, at the path
"/native_fields/(field_name)".
Should be run only by the master rank.
*/
void output_create_field(const std::string file_path, const std::string field_path,
 const std::string units, const int nx, const int ny, const int nz) {
    // Open the output.
    hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Create a 3D, nx x ny x nz dataspace.
    hsize_t dims[3] = {nx, ny, nz};
    hid_t grid_space = H5Screate_simple(3, dims, NULL);
    // Create the dataset.
    hid_t dataset = H5Dcreate(file, field_path.c_str(), H5T_IEEE_F32LE,
      grid_space, H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);

    if (dataset < 0) {
        printf("Error: could not create dataset. H5 returned %i.\n", dataset);
        exit(1);
    }

    // Don't forget to attach units attr.
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, strlen(units.c_str()));
    hid_t attr = H5Acreate(dataset, "units", str_type, scalar_space,
       H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, str_type, units.c_str());

    // Close resources.
    H5Aclose(attr);
    H5Tclose(str_type);
    H5Sclose(scalar_space);
    H5Dclose(dataset);
    H5Sclose(grid_space);
    H5Fclose(file);
}


/**
Write the only component in the multifab to the dataset given by field_name.
Uses hdf5-parallel.
*/
void output_write_field(const std::string file_path, const std::string field_path,
 MultiFab &mf) {

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // Create the file access prop list.
    hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(pa_plist, comm, info);

    // Open the file, and the group.
    hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);
    // Open the field dataset.
    hid_t dataset = H5Dopen(file, field_path.c_str(), H5P_DEFAULT);

    // Make sure the dataset is there.
    if (dataset < 0) {
        printf("Error on rank %i: Could not find dataset %s.\n", mpi_rank,
           field_path.c_str());
        exit(1);
    }

    // Grab the dataspace of the field dataset from file.
    hid_t file_dataspace = H5Dget_space(dataset);

    // Create collective io prop list.
    hid_t collective_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(collective_plist, H5FD_MPIO_COLLECTIVE);

    // Iterate over Fabs, select matching hyperslab and write.
    hid_t status;
    // slab lo index and shape.
    hsize_t slab_offsets[3], slab_dims[3];
    hid_t slab_dataspace;

    int write_count = 0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        // Get fab, lo and hi vectors.
        const Box& box = mfi.validbox();
        const int *lo_vec = box.loVect();
        const int *hi_vec = box.hiVect();

        // this is in fortran order.
        double *fab_data = mf[mfi].dataPtr();

        // We must make a copy of the fab data in c-order
        int nx = hi_vec[0] - lo_vec[0] + 1;
        int ny = hi_vec[1] - lo_vec[1] + 1;
        int nz = hi_vec[2] - lo_vec[2] + 1;
        size_t num_cells = nx * ny * nz;

        //size_t block_size = num_cells;
        std::vector<double> h5_data_vec(num_cells);
        double *h5_data = h5_data_vec.data();
        //double* h5_data = (double*) malloc(block_size * sizeof(double));
        if (h5_data == NULL) {
            printf("Error allocating h5 data block.\n");
            exit(1);
        }

        //int ix, iy, iz, c_index, f_index;
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    size_t c_index = (ix * ny + iy) * nz + iz;
                    size_t f_index = (iz * ny + iy) * nx + ix;
                    h5_data[c_index] = fab_data[f_index];
                }
            }
        }

        // Data is preped. Now write it.

        // Set slab offset and shape.
        slab_offsets[0] = lo_vec[0];
        slab_offsets[1] = lo_vec[1];
        slab_offsets[2] = lo_vec[2];
        slab_dims[0] = nx;
        slab_dims[1] = ny;
        slab_dims[2] = nz;

        // Create the slab space.
        slab_dataspace = H5Screate_simple(3, slab_dims, NULL);

        // Select the hyperslab matching this fab.
        status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET,
         slab_offsets, NULL, slab_dims, NULL);
        if (status < 0) {
            printf("Error on rank %i: could not select hyperslab.\n", mpi_rank);
            exit(1);
        }

        // Write this pencil.
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, slab_dataspace,
          file_dataspace, collective_plist, h5_data);
        if (status < 0) {
            printf("Error on rank %i: could not write hyperslab.\n", mpi_rank);
            exit(1);
        }

        H5Sclose(slab_dataspace);
        write_count++;
    }

    // Close HDF5 resources.
    H5Pclose(collective_plist);
    H5Sclose(file_dataspace);
    H5Dclose(dataset);
    H5Fclose(file);
    H5Pclose(pa_plist);
}


size_t proc_status_value(const std::string& field) {
    std::ifstream in("/proc/self/status");
    std::string line;
    while (in) {
        std::getline(in, line);
        if (line.compare(0, field.length(), field) == 0) {
            std::istringstream iss(line);
            std::string f; size_t res;
            iss >> f >> res;
            return res;
        }
    }
    return 0;
}


void printHWM(const char* info, MPI_Comm comm) {
    int io_rank=0, mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    //size_t hwm = proc_status_value("VmHWM");
    size_t hwm = proc_status_value("VmRSS");
    size_t avg_hwm = 0, min_hwm = 0, max_hwm = 0;
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &min_hwm, 1, MPI_LONG, MPI_MIN, io_rank, comm));
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &max_hwm, 1, MPI_LONG, MPI_MAX, io_rank, comm));
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &avg_hwm, 1, MPI_LONG, MPI_SUM, io_rank, comm));

    if (mpi_rank == io_rank) {
        printf("%s: Max = %zu, Min = %zu, Average = %zu\n", info, max_hwm, min_hwm, avg_hwm/mpi_size);
        fflush(stdout);
    }
}


int main(int argc, char **argv) {
    amrex::Initialize(argc, argv);

    // MPI info...
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    CHK(printHWM("Start", comm));

    //
    // Parameter parsing
    //

    if (argc > 4 || argc < 3) {
        print_usage(argc, argv);
    }

    ParmParse pp;
    if (pp.contains("help")) {
        print_usage(argc, argv);
    }
    bool verbose = false;
    if (pp.contains("verbose")) {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    std::string input_path;
    pp.get("input_path", input_path);
    std::string output_path;
    pp.get("output_path", output_path);
    std::string params_file = input_path + "/the_parameters";
    pp.query("params_file", params_file);

    //
    // Open input plotfile
    //

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl << "Opening " << input_path << " ...";
    }

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(input_path, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << " done." << std::endl;
    }

    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel();
    int Nlev = finestLevel + 1;

    //
    // Read metadata.
    //

    // Read comoving_a file
    Real comoving_a;
    if (ParallelDescriptor::IOProcessor()) {
        std::string a_file_path = input_path + "/comoving_a";
        std::ifstream a_file;
        a_file.open(a_file_path.c_str(), std::ios::in);
        if (!a_file.good()) {
            amrex::FileOpenFailed(a_file_path);
        }
        a_file >> comoving_a;
    }

    double z = 1.0/comoving_a - 1.0;

    // Read the_parameters
    double h = parseByName(params_file, "nyx.comoving_h");
    double omega_b = parseByName(params_file, "nyx.comoving_OmB");
    double omega_m = parseByName(params_file, "nyx.comoving_OmM");
    double omega_l = 1.0-omega_m;
    double domain_size = parseByName(params_file, "geometry.prob_hi") * h;

    Print() << "Read file: " << params_file << std::endl;
    Print() << "Parameters:  ";
    Print() << "Om = " << omega_m << ";  " << "Ob = " << omega_b << ";  " 
    << "h_0 = " << h << ";  " << "L = " << domain_size << std::endl;

    if (h < 0.0 || omega_b < 0.0 || omega_m < 0.0 || domain_size < 0.0) {
        Print() << "Error reading the_parameters file!" << std::endl;
        exit(-1);
    }

    // Parallel info... We might not need it.
    int MyPE = ParallelDescriptor::MyProc();
    int NumPEs = ParallelDescriptor::NProcs();

    //
    // Specify the data we need from the plotfile.
    //

    // For all lyaf analysis, we will need:
    //   - dark matter density
    //   - gas density
    //   - gas momemtum x, y, z
    //   - gas temperature

    // check if neutrinos are there:
    // NPC directory = neutrinos on
    std::string fn=input_path+"/NPC";
    bool neutrinos=neutrinos_on(fn);

    fn=input_path+"/job_info";
    bool append=append_on(fn);
    bool hydro=append ? false : hydro_on(fn);
    neutrinos=append ? false : neutrinos;
    bool dmvars=append ? true : dmvars_on(fn);

    const int nGhost = 0;
    int nComp = 0;
    if (dmvars) {
        Print() << "DMvars on" << endl;
        nComp+=4;
    }

    if (hydro) {
        Print() << "Hydro on" << endl;
        nComp+=5;
    }

    if (neutrinos) {
        Print() << "Neutrinos on" << endl;
        nComp+=1;
    }

    if (append) {
        Print() << "Append on" << endl;
        nComp+=0;
    }

    Vector<int> comps(nComp);
    int val = -1;
    if (dmvars) {
    //Dark matter fields
    Print() << "Converting (dm) particle quantities \n";
    int i_dm_density(amrData.StateNumber("particle_mass_density"));
    comps[0] = i_dm_density;
    Print() << "   particle mass density \n";
    int i_dm_vx(amrData.StateNumber("particle_x_velocity"));
    comps[1] = i_dm_vx;
    Print() << "   particle x velocity \n";
    int i_dm_vy(amrData.StateNumber("particle_y_velocity"));
    comps[2]= i_dm_vy;
    Print() << "   particle y velocity \n";
    int i_dm_vz(amrData.StateNumber("particle_z_velocity"));
    comps[3]= i_dm_vz;
    Print() << "   particle z velocity \n";

    val = 3;
    }

    //Hydro fields
    if (hydro) {
        Print() << "Converting hydro quantities \n";
        int i_gas_density(amrData.StateNumber("density"));
        comps[val+1] = i_gas_density;
        Print() << "   density \n";
        int i_xmom(amrData.StateNumber("xmom"));
        comps[val+2] = i_xmom;
        Print() << "   x momentum \n";
        int i_ymom(amrData.StateNumber("ymom"));
        comps[val+3] = i_ymom;
        Print() << "   y momentum \n";
        int i_zmom(amrData.StateNumber("zmom"));
        comps[val+4] = i_zmom;
        Print() << "   z momentum \n";
        int i_temp(amrData.StateNumber("Temp"));
        comps[val+5] = i_temp;
        Print() << "   temperature \n";
        val += 5;
    }

    //Neutrino fields
    if (neutrinos) {
        Print() << "Converting neutrino quantities \n";
        int i_neu_density(amrData.StateNumber("neutrino_mass_density"));
        comps[val+1] = i_neu_density;
        Print() << "   neutrino mass density \n";
        val += 1;
    }

    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    Vector<std::string> inVarNames(nComp);
    Vector<int> destFillComps(nComp);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl << "Converting the following states: "
        << std::endl;
    }
    for (int i = 0; i < nComp; ++i) {
        inVarNames[i] = plotVarNames[comps[i]];
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "    " << amrData.StateNumber(inVarNames[i])
            << " (" << inVarNames[i] << ")" << std::endl;
        }
        destFillComps[i] = i;
    }

    //
    // Make boxes and boxarray.
    //

    // Grab the number of cells in each direction (should be the same, but
    // we're keeping things general).
    Box pd(amrData.ProbDomain()[0]);
    int64_t grid_nx = pd.bigEnd(0) - pd.smallEnd(0) + 1;
    int64_t grid_ny = pd.bigEnd(1) - pd.smallEnd(1) + 1;
    int64_t grid_nz = pd.bigEnd(2) - pd.smallEnd(2) + 1;
    int64_t num_cells = grid_nx * grid_ny * grid_nz;

    if (ParallelDescriptor::IOProcessor()) {
        printf("The grid has a shape of (%i, %i, %i), %.3e cells total.\n\n",
           (int)grid_nx, (int)grid_ny, (int)grid_nz, (double)num_cells);
        printf("Making skewer chunks.\n");
        fflush(stdout);
    }

    // Check if we can split along x evenly.
    if (grid_nx % mpi_size != 0) {
        if (ParallelDescriptor::IOProcessor()) {
            printf("ERROR: domain decomposition.\n");
            printf("The number of MPI ranks must fit evenly in the number of cells along x.\n");
        }
        MPI_Barrier(comm);
        amrex::Finalize();
    }

    int chunk_size = grid_nx / mpi_size;

    // The box for z-pencils
    Box bx_pencil(amrData.ProbDomain()[0]);
    // List of pencil boxes.
    BoxList box_list;
    // indexes
    int i, ix_lo, ix_hi;
    for (i = 0; i < mpi_size; ++i) {
        ix_lo = chunk_size * i;
        ix_hi = ix_lo + chunk_size - 1;

        Box skewer;
        skewer.setSmall(0, ix_lo);
        skewer.setBig(0, ix_hi);
        skewer.setSmall(1, 0);
        skewer.setBig(1, grid_ny - 1);
        skewer.setSmall(2, 0);
        skewer.setBig(2, grid_nz - 1);

        box_list.push_back(skewer);
    }

    // Make box array from the skewer box list.
    BoxArray ba(box_list);

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Creating multifabs." << std::endl;
    }

    // Make the multifabs.
    int num_comps = 1;  // keeping it simple...
    int ng = 0;         // don't add any ghost zones.
    int level = 0;
    int comp_start = 0;

    DistributionMapping dmap(ba, ParallelDescriptor::NProcs());

    MultiFab mf1(ba, dmap, num_comps, ng);
    MultiFab mf2(ba, dmap, num_comps, ng);

    //
    // Start outputting.
    //

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Preparing output file." << std::endl;
    }

    // Create on master process.
    if (ParallelDescriptor::IOProcessor() && !append) {
        // create the file.
        output_create(output_path.c_str());

        // write metadata.
        output_write_metadata(output_path.c_str(),
            grid_nx, grid_ny, grid_nz, domain_size,
            omega_b, omega_m, omega_l, h, z);
    }

    ParallelDescriptor::Barrier();
    Print() << "Before reads allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;
    CHK(printHWM("Before read/write", comm));

    //
    // Dark matter particles, always present
    //

    std::string field_path;
    if (dmvars)
    {
    if (hydro || neutrinos || append) {
        field_path = "native_fields/dm_density";
    } else {
        field_path = "native_fields/matter_density";
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nReading particle_mass_density." << std::endl;
    }

    amrData.FillVar(mf1, level, "particle_mass_density", comp_start);
    amrData.FlushGrids(amrData.StateNumber("particle_mass_density"));

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Computing mean." << std::endl;
    }

    double mean_dm_density = mf1.norm1() / num_cells;

    if (ParallelDescriptor::IOProcessor()) {
        printf("Mean DM density: %.14e Msun/Mpc**3\n", mean_dm_density);
        fflush(stdout);
    }

    // Convert to mean units.
    double mean_dm_density_inv = 1.0/mean_dm_density;
    mf1.mult(mean_dm_density_inv);

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Writing to file." << std::endl;
        output_create_field(output_path, field_path, "(mean)",
            grid_nx, grid_ny, grid_nz);
    }
    ParallelDescriptor::Barrier();
    output_write_field(output_path, field_path, mf1);
    ParallelDescriptor::Barrier();

    // Particle velocities

    field_path = "native_fields/particle_vx";

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nReading particle_x_velocity." << std::endl;
    }

    amrData.FillVar(mf1, level, "particle_x_velocity", comp_start);
    amrData.FlushGrids(amrData.StateNumber("particle_x_velocity"));

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Writing to file." << std::endl;
        output_create_field(output_path, field_path, "cm/s",
            grid_nx, grid_ny, grid_nz);
    }
    // Fix units, km/s -> cm/s.
    mf1.mult(1.0e5);
    ParallelDescriptor::Barrier();
    output_write_field(output_path, field_path, mf1);
    ParallelDescriptor::Barrier();

    field_path = "native_fields/particle_vy";

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nReading particle_y_velocity." << std::endl;
    }

    amrData.FillVar(mf1, level, "particle_y_velocity", comp_start);
    amrData.FlushGrids(amrData.StateNumber("particle_y_velocity"));

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Writing to file." << std::endl;
        output_create_field(output_path, field_path, "cm/s",
            grid_nx, grid_ny, grid_nz);
    }
    // Fix units, km/s -> cm/s.
    mf1.mult(1.0e5);
    ParallelDescriptor::Barrier();
    output_write_field(output_path, field_path, mf1);
    ParallelDescriptor::Barrier();

    field_path = "native_fields/particle_vz";

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nReading particle_z_velocity." << std::endl;
    }

    amrData.FillVar(mf1, level, "particle_z_velocity", comp_start);
    amrData.FlushGrids(amrData.StateNumber("particle_z_velocity"));

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Writing to file." << std::endl;
        output_create_field(output_path, field_path, "cm/s",
            grid_nx, grid_ny, grid_nz);
    }
    // Fix units, km/s -> cm/s.
    mf1.mult(1.0e5);
    ParallelDescriptor::Barrier();
    output_write_field(output_path, field_path, mf1);
    ParallelDescriptor::Barrier();
    }
    //
    // Baryonic gas, optional
    //

    if(hydro) {
        field_path = "native_fields/baryon_density";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading density." << std::endl;
        }

        amrData.FillVar(mf2, level, "density", comp_start);
        amrData.FlushGrids(amrData.StateNumber("density"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Computing mean." << std::endl;
        }

        double mean_baryon_density = mf2.norm1() / num_cells;

        if (ParallelDescriptor::IOProcessor()) {
            printf("Mean baryon density: %.14e Msun/Mpc**3\n", mean_baryon_density);
            fflush(stdout);
        }

        // Convert to mean units.
        double mean_baryon_density_inv = 1.0/mean_baryon_density;
        mf2.mult(mean_baryon_density_inv);

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "(mean)",
                grid_nx, grid_ny, grid_nz);
        }
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf2);
        ParallelDescriptor::Barrier();

        // Density is in multifab 2, but in the wrong units.
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Restoring density..." << std::endl;
        }

        mf2.mult(mean_baryon_density);

        // Gas velocity componens

        field_path = "native_fields/velocity_x";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading xmom." << std::endl;
        }

        amrData.FillVar(mf1, level, "xmom", comp_start);
        amrData.FlushGrids(amrData.StateNumber("xmom"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "cm/s",
                grid_nx, grid_ny, grid_nz);
        }
        // Divide out density to get velocity.
        mf1.divide(mf2,comp_start,num_comps,ng);
        // Fix units, km/s -> cm/s.
        mf1.mult(1.0e5);
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf1);
        ParallelDescriptor::Barrier();

        field_path = "native_fields/velocity_y";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading ymom." << std::endl;
        }

        amrData.FillVar(mf1, level, "ymom", comp_start);
        amrData.FlushGrids(amrData.StateNumber("ymom"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "cm/s",
                grid_nx, grid_ny, grid_nz);
        }
        // Divide out density to get velocity.
        mf1.divide(mf2,comp_start,num_comps,ng);
        // Fix units, km/s -> cm/s.
        mf1.mult(1.0e5);
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf1);
        ParallelDescriptor::Barrier();

        field_path = "native_fields/velocity_z";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading zmom." << std::endl;
        }

        amrData.FillVar(mf1, level, "zmom", comp_start);
        amrData.FlushGrids(amrData.StateNumber("zmom"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "cm/s",
                grid_nx, grid_ny, grid_nz);
        }
        // Divide out density to get velocity.
        mf1.divide(mf2,comp_start,num_comps,ng);
        // Fix units, km/s -> cm/s.
        mf1.mult(1.0e5);
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf1);
        ParallelDescriptor::Barrier();

        field_path = "native_fields/temperature";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading Temp." << std::endl;
        }

        amrData.FillVar(mf1, level, "Temp", comp_start);
        amrData.FlushGrids(amrData.StateNumber("Temp"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "K",
                grid_nx, grid_ny, grid_nz);
        }
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf1);
        ParallelDescriptor::Barrier();
    }

    //
    // Neutrinos, optional
    //

    if(neutrinos) {
        field_path = "native_fields/neutrino_density";

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "\nReading neutrino_mass_density." << std::endl;
        }
        amrData.FillVar(mf1, level, "neutrino_mass_density", comp_start);
        amrData.FlushGrids(amrData.StateNumber("neutrino_mass_density"));

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Computing mean." << std::endl;
        }

        double mean_neu_density = mf1.norm1() / num_cells;

        if (ParallelDescriptor::IOProcessor()) {
            printf("Mean NPC density: %.14e Msun/Mpc**3\n", mean_neu_density);
            fflush(stdout);
        }

        // Convert to mean units.
        double mean_neu_density_inv = 1.0/mean_neu_density;
        mf1.mult(mean_neu_density_inv);

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing to file." << std::endl;
            output_create_field(output_path, field_path, "(mean)",
                grid_nx, grid_ny, grid_nz);
        }
        ParallelDescriptor::Barrier();
        output_write_field(output_path, field_path, mf1);
        ParallelDescriptor::Barrier();
    }

    Print() << "End allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;

    amrex::Finalize();
    return 0;
}
