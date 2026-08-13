#include <GadgetFileWriter.H>

void write_vtk_binary(const std::string& filename,
                      const std::vector<float>& pos,
                      const std::vector<float>& vel,
                      const std::vector<int32_t>& ids)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Failed to open VTK output file: " << filename << "\n";
        return;
    }

    uint64_t ntot = ids.size();

    // 1. Write the ASCII Header text
    out << "# vtk DataFile Version 3.0\n";
    out << "Gadget ICs Snapshot Visualization\n";
    out << "BINARY\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // 2. Write Points (Positions)
    out << "POINTS " << ntot << " float\n";
    for (size_t i = 0; i < pos.size(); ++i) {
        float swapped = to_big_endian(pos[i]);
        out.write(reinterpret_cast<const char*>(&swapped), sizeof(float));
    }
    out << "\n"; // VTK requires a newline after binary blocks

    // 3. Write Topology (Treating every particle as a VTK_VERTEX / Cell Type 1)
    // VTK expects: CELLS [num_cells] [size of cell list array]
    // Each cell list entries are: [num_points_in_cell, point_id_0, point_id_1...]
    out << "CELLS " << ntot << " " << (ntot * 2) << "\n";
    for (uint32_t i = 0; i < ntot; ++i) {
        int32_t cell_size = to_big_endian(static_cast<int32_t>(1));
        int32_t pt_idx = to_big_endian(static_cast<int32_t>(i));

        out.write(reinterpret_cast<const char*>(&cell_size), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&pt_idx), sizeof(int32_t));
    }
    out << "\n";

    // Cell Types (1 = VTK_VERTEX)
    out << "CELL_TYPES " << ntot << "\n";
    int32_t cell_type = to_big_endian(static_cast<int32_t>(1));
    for (uint32_t i = 0; i < ntot; ++i) {
        out.write(reinterpret_cast<const char*>(&cell_type), sizeof(int32_t));
    }
    out << "\n";

    // 4. Write Point Data (Velocities and Particle IDs)
    out << "POINT_DATA " << ntot << "\n";

    // Velocity vectors
    out << "VECTORS velocity float\n";
    for (size_t i = 0; i < vel.size(); ++i) {
        float swapped = to_big_endian(vel[i]);
        out.write(reinterpret_cast<const char*>(&swapped), sizeof(float));
    }
    out << "\n";

    // ID scalars
    out << "SCALARS particle_id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < ids.size(); ++i) {
        int32_t swapped = to_big_endian(ids[i]);
        out.write(reinterpret_cast<const char*>(&swapped), sizeof(int32_t));
    }
    out << "\n";

    out.close();
    std::cout << "Successfully generated VTK file: " << filename << "\n";
}

// ------------------------------------------------------------
// Write a Gadget block
// ------------------------------------------------------------
void write_block(std::ofstream& out,
                 const void* data,
                 uint32_t nbytes)
{
    out.write(reinterpret_cast<const char*>(&nbytes), sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(data), nbytes);
    out.write(reinterpret_cast<const char*>(&nbytes), sizeof(uint32_t));
}

void WriteGadgetFileBlock(const double comoving_OmM,
                          const double comoving_h,
                          const double comoving_a,
                          const std::vector<float>& vec_pos,
                          const std::vector<float>& vec_vel,
                          const std::vector<int64_t>& vec_id,
                          const int total_num_blocks, 
                          uint64_t num_particles_in_block,
                          uint64_t total_num_particles,
                          const double domain_size,
                          const std::string filename_str)
{
    GadgetHeader hdr{};
    std::memset(&hdr, 0, sizeof(hdr));

    if (num_particles_in_block > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error(
            "Too many particles in a single Gadget file: " +
            std::to_string(num_particles_in_block) +
            " (maximum supported is " +
            std::to_string(std::numeric_limits<uint32_t>::max()) + ").");
    }

    if (vec_pos.size() != 3 * num_particles_in_block ||
        vec_vel.size() != 3 * num_particles_in_block ||
        vec_id.size()  != num_particles_in_block) {
        throw std::runtime_error(
            "Inconsistent particle data sizes when writing Gadget file.");
    }

    hdr.num_particles[1] =
    static_cast<uint32_t>(num_particles_in_block);

    hdr.num_total_particles[1] =
    static_cast<uint32_t>(total_num_particles & 0xFFFFFFFFULL);

    hdr.num_total_particles_hw[1] =
    static_cast<uint32_t>(total_num_particles >> 32);

    // Constant particle mass
    hdr.particle_masses[1] = 1.0;

    hdr.scale_factor = comoving_a;
    hdr.redshift     = 1.0/comoving_a - 1.0;

    hdr.flag_sfr      = 0;
    hdr.flag_feedback = 0;
    hdr.flag_cooling  = 0;

    hdr.num_files_per_snapshot = total_num_blocks;

    //hdr.box_size      = domain_size * comoving_h * 1.0/comoving_a;
    hdr.box_size      = domain_size;
    hdr.omega_0       = comoving_OmM;
    hdr.omega_lambda  = 1.0 - comoving_OmM;
    hdr.h_0           = comoving_h;

    hdr.flag_stellarage   = 0;
    hdr.flag_metals       = 0;
    hdr.flag_entropy_ics  = 0;

    std::ofstream out(filename_str, std::ios::binary);
    if (!out) {
        throw std::runtime_error(
            "Failed to open Gadget file '" + filename_str + "' for writing.");
    }

    write_block(out,
                &hdr,
                sizeof(GadgetHeader));

    write_block(out,
                vec_pos.data(),
                static_cast<uint32_t>(
                vec_pos.size()*sizeof(float)));

    write_block(out,
                vec_vel.data(),
                static_cast<uint32_t>(
                vec_vel.size()*sizeof(float)));

    write_block(out,
                vec_id.data(),
                static_cast<uint64_t>(
                vec_id.size()*sizeof(int64_t)));

    out.close();
}
