#ifndef Nyx_output_H_
#define Nyx_output_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <array>
#include <vector>

#define GADGET_HEADER_SIZE 256

#include <fstream>
#include <iostream>
#include <vector>
#include <bit> // Required for std::byteswap (C++20)
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>

extern std::string inputs_name;

// ------------------------------------------------------------
// Gadget-2 header
// ------------------------------------------------------------
struct GadgetHeader
{
    uint32_t num_particles[6];
    double   particle_masses[6];
    double   scale_factor;
    double   redshift;
    int32_t  flag_sfr;
    int32_t  flag_feedback;
    uint32_t num_total_particles[6];
    int32_t  flag_cooling;
    int32_t  num_files_per_snapshot;
    double   box_size;
    double   omega_0;
    double   omega_lambda;
    double   h_0;
    int32_t  flag_stellarage;
    int32_t  flag_metals;
    int32_t  num_total_particles_hw[6];
    int32_t  flag_entropy_ics;
    char     unused[60];
};

// Runtime check for little endianness
static inline bool is_little_endian()
{
    uint16_t number = 0x1;
    char* numPtr = reinterpret_cast<char*>(&number);
    return (numPtr[0] == 1);
}

// Manual byte-swapping that safely handles both integers and floats
template <typename T>
static inline T to_big_endian(T val)
{
    if (!is_little_endian()) {
        return val;
    }

    T swapped;
    const char* src = reinterpret_cast<const char*>(&val);
    char* dest = reinterpret_cast<char*>(&swapped);
    size_t size = sizeof(T);

    for (size_t i = 0; i < size; ++i) {
        dest[i] = src[size - 1 - i];
    }

    return swapped;
}

void write_vtk_binary(const std::string& filename,
                      const std::vector<float>& pos,
                      const std::vector<float>& vel,
                      const std::vector<int32_t>& ids);

void write_block(std::ofstream& out,
                 const void* data,
                 uint32_t nbytes);

void WriteGadgetFileBlock(const double comoving_OmM,
                          const double comoving_h,
                          const double comoving_a,
                          const std::vector<float>& vec_pos,
                          const std::vector<float>& vec_vel,
                          const std::vector<int64_t>& vec_id,
                          const int total_num_blocks,
                          const uint64_t num_partcles_in_block,
                          const uint64_t total_num_particles,
                          const double domain_size,
                          const std::string filename_str);
#endif
