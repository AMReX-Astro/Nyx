#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <AMReX_ParallelDescriptor.H>

#include <LightConeParticle.H>

void SwapEnd(float& val) {
    // Swap endianess if necessary
    char* bytes = reinterpret_cast<char*>(&val);
    std::swap(bytes[0], bytes[3]);
    std::swap(bytes[1], bytes[2]);
}

void writeBinaryVTK(const std::string& filename, const std::vector<LightConeParticle>& particles) {

	int rank = amrex::ParallelDescriptor::MyProc();
	int size = amrex::ParallelDescriptor::NProcs();

    long int local_num_particles = particles.size();
    long int total_num_particles = local_num_particles;

    // Get total particles across all ranks
	amrex::ParallelDescriptor::ReduceLongSum(total_num_particles);

    // Compute offset for this rank's data
    size_t offset = 0;
    MPI_Exscan(&local_num_particles, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    // Header handling
    size_t header_size = 0;

    if (rank == 0) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            amrex::Abort("Error: Could not open file " + filename + "\n");
        }

        // Write the header
        file << "# vtk DataFile Version 2.0\n";
        file << "Particle Cloud Data\n";
        file << "BINARY\n";
        file << "DATASET POLYDATA\n";
        file << "POINTS " << total_num_particles << " float\n";

        // Determine header size
        file.seekp(0, std::ios::end);
        header_size = file.tellp();
        file.close();
    }

    // Broadcast the header size to all ranks
	amrex::ParallelDescriptor::Bcast(&header_size, 1, 0);

    // Use MPI collective I/O for binary data
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_file);

    // Compute byte offset for this rank
    size_t byte_offset = header_size + sizeof(float) * 3 * offset;

    // Prepare local data
    std::vector<float> data(3 * local_num_particles);
    for (size_t i = 0; i < local_num_particles; ++i) {
        data[3 * i] = particles[i].x;
        data[3 * i + 1] = particles[i].y;
        data[3 * i + 2] = particles[i].z;

        // Convert to big-endian if needed
        SwapEnd(data[3 * i]);
        SwapEnd(data[3 * i + 1]);
        SwapEnd(data[3 * i + 2]);
    }

    // Write particle data collectively
    MPI_File_write_at_all(mpi_file, byte_offset, data.data(), data.size(), MPI_FLOAT, MPI_STATUS_IGNORE);

    MPI_File_close(&mpi_file);

    if (rank == 0) {
        std::cout << "Successfully wrote VTK file: " << filename << "\n";
    }
}

void writeBinarySimple(const std::string& filename, const std::vector<LightConeParticle>& particles) {

	int rank = amrex::ParallelDescriptor::MyProc();
	int size = amrex::ParallelDescriptor::NProcs();

    long int local_num_particles = particles.size();

    // Compute offset for this rank's data
    size_t offset = 0;
    MPI_Exscan(&local_num_particles, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    // Header handling
    unsigned long header_size = 0;

    if (rank == 0) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            amrex::Abort("Error: Could not open file " + filename + "\n");
        }
        file.seekp(0, std::ios::end);
        header_size = file.tellp();
        file.close();
    }

    // Broadcast the header size to all ranks
	amrex::ParallelDescriptor::Bcast(&header_size, 1, 0);

    // Use MPI collective I/O for binary data
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_file);

    // Compute byte offset for this rank
    size_t byte_offset = header_size + sizeof(float) * 6 * offset;

    // Prepare local data
    std::vector<float> data(6 * local_num_particles);
    for (size_t i = 0; i < local_num_particles; ++i) {
        data[6 * i] = particles[i].x;
        data[6 * i + 1] = particles[i].y;
        data[6 * i + 2] = particles[i].z;
        data[6 * i + 3] = particles[i].vx;
        data[6 * i + 4] = particles[i].vy;
        data[6 * i + 5] = particles[i].vz;

        // Convert to big-endian if needed
        SwapEnd(data[6 * i]);
        SwapEnd(data[6 * i + 1]);
        SwapEnd(data[6 * i + 2]);
        SwapEnd(data[6 * i + 3]);
        SwapEnd(data[6 * i + 4]);
        SwapEnd(data[6 * i + 5]);
    }

    // Write particle data collectively
    MPI_File_write_at_all(mpi_file, byte_offset, data.data(), data.size(), MPI_FLOAT, MPI_STATUS_IGNORE);

    MPI_File_close(&mpi_file);

    if (rank == 0) {
        std::cout << "Successfully wrote VTK file: " << filename << "\n";
    }
}

