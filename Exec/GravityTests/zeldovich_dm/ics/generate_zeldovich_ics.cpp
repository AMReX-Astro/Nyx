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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

// because we aren't using boost.
const double PI = acos(-1.0);
const double grav_const = 6.673e-8;  // in cm^3 g^-1 s^-2
// units factor for critical density (rho_c), km^2 Mpc^-2 cm^-3 g -> M_sun Mpc^-3
// terms are (km^2 / cm^2) * (Mpc / cm) * (g / M_sun)
const double rho_c_units_factor = 1.0e10 * 3.08568025e24 * 5.02785431e-34;

// This is like our "header", so we can put `main` first.
// Define a point in 3d.
class Point {
private:

public:
    double x, y, z;

    Point(double a=0.0, double b=0.0, double c=0.0) {
        x = a;
        y = b;
        z = c;
    }
};

// A particle is just two points in 3d: position and velocity.
class Particle {
private:

public:
    Point position, velocity;

    Particle(Point pos=Point(0.0, 0.0, 0.0), Point vel=Point(0.0, 0.0, 0.0)) {
        position = pos;
        velocity = vel;
    }
};

double perturbed_x(double q, double current_z, double final_z, double k);
double perturbed_v_x(double q, double current_z, double final_z, double k,
                     double H_0);
double perturbed_rho(double q, double average_rho, double current_z,
                     double final_z, double k);
double get_scale_factor(double current_z);
Point get_unit_vector(Point vec);
void write_particle(ofstream& p_file, Particle p, double mass);
// end header.

// ================ //
// Params to adjust
// ================ //

// Redshifts
double initial_z = 100.0;
double caustic_z = 10.0;

// Required cosmology
double h = 0.675;  // the weird little h used in the def. of Hubble rate.
double H_0 = 100.0 * h;  // the Hubble rate in the standard km Mpc^-1 s^-1.

// The size of the box. Note that all lengths are in Mpc, not Mpc h^-1.
double box_length = 1.0/0.675;

// How many sheets should fit in the box. This determines the perturbation
// wavelength (k mode).
double num_sheets = 1;
double offset = 0.75*box_length;//PI/2;  // wave offset

// Particle number (mass determined by critical density)
long num_particles_1d = 128;

// Normal vector to the sheets
Point sheet_normal = get_unit_vector(Point(1.0, 0.0, 0.0));

// =========== //
// End params. You shouldn't have to adjust anything below this. I hope.
// =========== //

int main(int argc, char *argv[]) {
    // Particle things.
    long num_particles = pow(num_particles_1d, 3);
    double crit_density = (3.0 * pow(H_0, 2) / (8.0 * PI * grav_const)
                           * rho_c_units_factor);  // in M_sun / Mpc^3

    // Lengths in the problem.
    double half_cell = (box_length / num_particles_1d) / 2.0;
    double spacing = (float) box_length / num_particles_1d;  // distance between particles in one dimension
    double initial_a = get_scale_factor(initial_z);  // the current scale factor

    // Fit an integer number of sheets in box
    // The dumb way to do min of >2 floats, but it works fine for 3...
    double normal_length = min(box_length / sheet_normal.x,
                               min(box_length / sheet_normal.y,
                                   box_length / sheet_normal.z));
    double wavelength = normal_length / num_sheets;
    double k = 2 * PI / wavelength;

    // overfill control
    int p_fill_1d = 1.2 * num_particles_1d;
    int p_fill = pow(p_fill_1d, 3);
    double overflow_edge = box_length / 10.0;

    // dummy stuff
    double q, perturbation, x, y, z;
    Point p;
    Point * pp;

    // the big vector with all positions
    vector<Particle> particles, final_particles;

    // Give the user that warm, fuzzy feeling.
    std::cout << std::endl;
    std::cout << "Zeldovich pancake IC generator" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << std::endl;
    std::cout << "Initial redshift:\t" << initial_z << std::endl;
    std::cout << "Caustic redshift:\t" << caustic_z << std::endl;
    std::cout << "Box length:\t\t" << box_length << " Mpc" << std::endl;
    std::cout << "Perturbation length:\t" << wavelength << " Mpc" << std::endl;
    std::cout << "Critical density:\t" << crit_density << " M_sun Mpc^-3"
              << std::endl;
    std::cout << "Number of particles:\t" << num_particles_1d << "^3 = "
              << num_particles << std::endl << std::endl;

    std::cout << "Filling grid." << std::endl;

    // Fill particle_pos with regular grid. Note we haven't added elements to
    // the vector yet, so use push_back.
    for (int i = 0; i != p_fill_1d; ++i) {
        x = -overflow_edge + (i * spacing) + half_cell;
        for (int j = 0; j != p_fill_1d; ++j) {
            y = -overflow_edge + (j * spacing) + half_cell;
            for (int k = 0; k != p_fill_1d; ++k) {
                z = -overflow_edge + (k * spacing) + half_cell;
                particles.push_back(Particle(Point(x, y, z),
                                             Point(0.0, 0.0, 0.0)));
            }
        }
    }

    std::cout << "Perturbing grid." << std::endl;

    // Now "perturb" the x positions and velocities
    for (int i = 0; i != particles.size(); ++i) {
        pp = &(particles[i].position);

        // zeldovich grid position
        q = (pp->x * sheet_normal.x + pp->y * sheet_normal.y
             + pp->z * sheet_normal.z) + offset;

        // zeldovich position perturbation
        perturbation = perturbed_x(q, initial_z, caustic_z, k) - q;

        pp->x += perturbation * sheet_normal.x;
        pp->y += perturbation * sheet_normal.y;
        pp->z += perturbation * sheet_normal.z;

        pp = &(particles[i].velocity);

        // note perturbed_v_x will be in comoving km s^-1 and we need it in
        // proper km s^-1, so multiply by ``initial_a``.
        perturbation = perturbed_v_x(q, initial_z, caustic_z, k, H_0) * initial_a;

        // don't need ``+=`` because velocity components are already 0.0
        pp->x = perturbation * sheet_normal.x;
        pp->y = perturbation * sheet_normal.y;
        pp->z = perturbation * sheet_normal.z;
    }

    std::cout << "Filtering particles outside of domain." << std::endl;

    //double fudged_box_length = box_length + 0.001;

    // Run through the particles and add the particles inside the box to
    // ``final_positions``. This is cheaper than Vector::erase, O(N).
    for (int i = 0; i != particles.size(); ++i) {
        p = particles[i].position;

        if (p.x > 0.0 && p.x < box_length && p.y > 0.0 && p.y < box_length
            && p.z > 0.0 && p.z < box_length)
        {
            // this particle is in the box, so add it to ``final_positions``
            final_particles.push_back(particles[i]);
        }
    }

    unsigned int final_num_particles = final_particles.size();
    double particle_mass = crit_density * pow(box_length, 3) / final_num_particles;

    std::cout << std::endl << "Saving " << final_num_particles
              << " particles to file." << std::endl;
    std::cout << "Mass density:\t\t"
              << final_num_particles * particle_mass / pow(box_length, 3)
              << " M_sun Mpc^-3" << std::endl;
    std::cout << "Particle mass:\t\t" << particle_mass << " M_sun" << std::endl;

    // params output
    string params_filename = "zeldovich_params.csv";
    ofstream params_file;
    params_file.open(params_filename.c_str());  // need a char* here

    // write headers
    params_file << "initial_z,caustic_z,box_length,k,offset,H_0,"
                << "normal_x,normal_y,normal_z" << std::endl;
    // write param values
    params_file << std::scientific;
    params_file.precision(15);
    params_file << initial_z << "," << caustic_z << "," << box_length << ","
                << k << "," << offset << "," << H_0 << "," << sheet_normal.x
                << "," << sheet_normal.y << "," << sheet_normal.z << ","
                << std::endl;

    // particles output
    string p_filename = "zeldovich_particles.ascii";
    ofstream particle_file;
    particle_file.open(p_filename.c_str());  // need a char* here

    // first line has to be total number of particles
    particle_file << final_num_particles << std::endl;

    particle_file << std::scientific;
    particle_file.precision(15);  // just in case

    for (int i = 0; i != final_num_particles; ++i) {
        write_particle(particle_file, final_particles[i], particle_mass);
    }

    std::cout << std::endl << "Done!" << std::endl;
    std::cout << "The particle data is in `" << p_filename << "`" << std::endl;
    std::cout << "The params are in `" << params_filename << "`"
              << std::endl << std::endl;

    // close it and exit
    particle_file.close();

    return 0;
}


// =========================== //
// Utils for Zeldovich problem
// =========================== //

// Zeldovich problem setup functions
double perturbed_x(double q, double current_z, double final_z, double k) {
    // Takes the unperturbed position of the particle on a regular grid, and
    // returns the "perturbed" position, in comoving units.

    return (q - (1.0 + final_z) / (1.0 + current_z) * sin(k * q) / k);
}

double perturbed_v_x(double q, double current_z, double final_z, double k,
                     double H_0) {
    // Returns the "perturbed" velocity, in comoving units.
    return (-H_0 * (1.0 + final_z) * sqrt(1.0 + current_z)
            * sin(k * q) / k);
}

double perturbed_rho(double q, double average_rho, double current_z,
                     double final_z, double k) {
    // Returns the "perturbed" density, in proper units.
    return (average_rho / (1.0 - (1.0 + final_z) / (1.0 + current_z)
            * cos(k * q)));
}
// end Zeldovich problem functions

Point get_unit_vector(Point vec) {
    // take `vec` and divide it by its norm, making it a unit vecotr in the same
    // direction.
    double norm = sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
    vec.x /= norm;
    vec.y /= norm;
    vec.z /= norm;

    return vec;
}

double get_scale_factor(double current_z) {
    // Take in the redshift, `z`, return the scale factor, `a`.
    return 1.0 / (1.0 + current_z);
}

void write_particle(ofstream& p_file, Particle p, double mass) {
    // Write the particle data, (x, y, z, mass, v_x, v_y, v_z), to the given
    // file. This is just an abstraction for the format Nyx expects.
    p_file << p.position.x << " " << p.position.y << " " << p.position.z << " "
           << mass << " " << p.velocity.x << " " << p.velocity.y << " "
           << p.velocity.z << std::endl;
}
