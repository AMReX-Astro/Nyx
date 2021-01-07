#include <iostream>
#include <stdint.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

#ifdef MEM_CHECK
#define CHK(x) x
#else
#define CHK(x)
#endif

using namespace amrex;


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


void printHWM(const char* info) {
    int io_rank=0, mpi_size=1, mpi_rank=0;
    //size_t hwm = proc_status_value("VmHWM");
    size_t hwm = proc_status_value("VmRSS");
    size_t avg_hwm = 0, min_hwm = 0, max_hwm = 0;

#ifdef AMREX_USE_MPI
    MPI_Comm comm = ParallelDescriptor::Communicator();
    BL_MPI_REQUIRE(MPI_Comm_size(comm, &mpi_size));
    BL_MPI_REQUIRE(MPI_Comm_rank(comm, &mpi_rank));
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &min_hwm, 1, MPI_LONG, MPI_MIN, io_rank, comm));
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &max_hwm, 1, MPI_LONG, MPI_MAX, io_rank, comm));
    BL_MPI_REQUIRE(MPI_Reduce(&hwm, &avg_hwm, 1, MPI_LONG, MPI_SUM, io_rank, comm));
#else
    avg_hwm = hwm;
    min_hwm = hwm;
    max_hwm = hwm;
#endif

    if (mpi_rank == io_rank) {
        printf("%s: Max = %zu, Min = %zu, Average = %zu\n", info, max_hwm, min_hwm, avg_hwm/mpi_size);
        fflush(stdout);
    }
}


int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);

    Print() << "Starting, allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;
    CHK(printHWM("Start"));
    
    uint64_t num_cells = 2048;
    int max_grid_size = 256;
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(num_cells - 1,
                                   num_cells - 1,
                                   num_cells - 1));
    const Box domain(domain_lo, domain_hi);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dmap(ba);
    
    MultiFab zhi(ba, dmap, 1, 0);

    std::string file_name = "zhi.bin";
    std::ifstream ifs;
    ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
    if (!ifs ) {
        amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
        amrex::Abort();
    }
    Vector<float> values(num_cells*num_cells*num_cells);
    ifs.read((char*) &values[0], num_cells*num_cells*num_cells*sizeof(float));    

    Print() << "After file read, allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;
    CHK(printHWM("After file read"));

    for (MFIter mfi(zhi); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        for (unsigned x = box.smallEnd(0); x <= box.bigEnd(0); ++x) {
            for (unsigned y = box.smallEnd(1); y <= box.bigEnd(1); ++y) {
                for (unsigned z = box.smallEnd(2); z <= box.bigEnd(2); ++z) {
                    IntVect iv(x, y, z);
                    uint64_t index = x*num_cells*num_cells + y*num_cells + z;
                    zhi[mfi](iv) = values[index];
                }
            }
        }
    }

    Print() << "Before write, allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;
    CHK(printHWM("Before write"));

    amrex::VisMF::Write(zhi, "zhi/zhi");

    Print() << "End, allocation in Fabs: " << TotalBytesAllocatedInFabs()/(1024*1024.) << " Mb." << std::endl;
    CHK(printHWM("End"));

    amrex::Finalize();
}
