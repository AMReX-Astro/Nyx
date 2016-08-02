// Zeldovich pancake problem IC generator
// --------------------------------------
//
// The `main` function generates an ascii file containing the particle data for
// the Zeldovich pancake problem with the parameters defined there.
//
// Compile and run this code on your own with something simple like:
//
//     $ g++ generate_zeldovich_ics.cpp -o generate_zeldovich_ics.ex
//     $ ./generate_zeldovich_ics.ex
//
// Writing the data to an ascii file is incredibly wasteful, but it's the
// easiest way for Nyx to read the data for now.
//
// Written by: Casey W. Stark
// Copyright 2011
// Date: September 2011
// Included by Casey W. Stark, September 2, 2011
//
// Modified: date by person.
// Comments:
//

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// because we aren't using boost.
const double PI = acos(-1.0);

double perturbed_x(double q, double current_z, double final_z, double k)
{
    // Takes the unperturbed position of the particle on a regular grid, and
    // returns the "perturbed" position, in comoving units.

    return (q - (1.0 + final_z) / (1.0 + current_z) * sin(k * q) / k);
}

double perturbed_v_x(double q, double current_z, double final_z,
                     double k, double H_0)
{
    // Returns the "perturbed" velocity, in comoving units.
    return (-H_0 * (1.0 + final_z) * sqrt(1.0 + current_z) * sin(k * q) / k);
}

double perturbed_rho(double q, double average_rho, double current_z,
                     double final_z, double k)
{
    // Returns the "perturbed" density, in proper units.
    return (average_rho / (1.0 - (1.0 + final_z) / (1.0 + current_z)
            * cos(k * q)));
}

double get_scale_factor(double current_z)
{
    // Take in the redshift, `z`, return the scale factor, `a`.
    return 1.0 / (1.0 + current_z);
}

void write_particle(ofstream& p_file, double x, double y, double z,
                    double mass, double v_x, double v_y, double v_z)
{
    // Write the particle data, (x, y, z, mass, v_x, v_y, v_z), to the given
    // file. This is just an abstraction for the format Nyx expects.
    p_file << x << " " << y << " " << z << " " << mass << " " << v_x << " "
           << v_y << " " << v_z << std::endl;
}

int main(int argc, char *argv[])
{
    // these are params you can adjust.
    double initial_z = 10.0;
    double final_z = 5.0;

    // required cosmology
    double h = 0.50;  // the weird little h used in the def. of Hubble rate.
    double H_0 = 100.0 * h;  // the Hubble rate in the standard km Mpc^-1 s^-1

    // Lengths. Note these are all in Mpc, not Mpc h^-1.
    double box_length = 15.0;
    double perturbation_wavelength = 5.0;
    double k = 2 * PI / perturbation_wavelength;

    // particle number and mass
    long num_particles_1d = 64;
    double crit_density = 5e11;  // critcal density in M_sun / Mpc^3

    // You should never have to adjust anything below this. I hope.

    // extra vars
    long num_particles = pow(num_particles_1d, 3);
    double particle_mass = crit_density * pow(box_length, 3) / num_particles;

    double grid_x, zel_x, zel_v_x;
    double grid_y[num_particles_1d];
    double v_y = 0.0;
    // no need to define z's -- just a copy of y

    // distance between particles in one dimension
    double p_spacing = (float) box_length / num_particles_1d;

    // the current scale factor
    double current_a = get_scale_factor(initial_z);

    string p_filename = "simple_particles.ascii";
    ofstream particle_file;  // for the output file

    // Give the user that warm, fuzzy feeling.
    std::cout << std::endl;
    std::cout << "Zeldovich pancake IC generator" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << std::endl;
    std::cout << "Number of particles:\t" << num_particles_1d << "^3 = "
              << num_particles << std::endl;
    std::cout << "Particle mass:\t\t" << particle_mass << " M_sun" << std::endl;
    std::cout << "Box Length:\t\t" << box_length << " Mpc" << std::endl;
    std::cout << std::endl;

    // compute particle positions for the non-perturbed dims
    for (int y = 0; y != num_particles_1d; ++y) {
        grid_y[y] = p_spacing * y;
    }

    // open particle file
    particle_file.open(p_filename.c_str());  // need a char* here

    // first line has to be total number of particles
    particle_file << num_particles << std::endl;

    std::cout << "Computing and writings cells at x = ..." << std::endl;

    // loop through particles and set the perturbed 1D position and velocity
    for (int x = 0; x != num_particles_1d; ++x) {
        // figure out the correct position on the grid
        grid_x = p_spacing * x;

        // give the user something to stare at
        std::cout << grid_x << "\t";
        std::cout.flush();  // need to flush the buffer to print without `endl`

        zel_x = perturbed_x(grid_x, initial_z, final_z, k);

        // multiply by current_a because the function returns comoving and nyx
        // uses proper velocity.
        zel_v_x = perturbed_v_x(grid_x, initial_z, final_z, k, H_0) * current_a;

        for (int y = 0; y != num_particles_1d; ++y) {
            for (int z = 0; z != num_particles_1d; ++z) {
                write_particle(particle_file, zel_x, grid_y[y], grid_y[z],
                               particle_mass, zel_v_x, v_y, v_y);
            }
        }
    }

    std::cout << std::endl << std::endl << "Done!" << std::endl;
    std::cout << "The output is in the ascii file `" << p_filename << "`"
              << std::endl << std::endl;

    // close it and exit
    particle_file.close();

    return 0;
}
