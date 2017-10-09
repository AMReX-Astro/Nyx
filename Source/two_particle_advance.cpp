#ifdef  GRAVITY

#include "Nyx.H"
#include "Nyx_F.H"
#include "Gravity.H"
#include <AMReX_Particles_F.H>

using namespace amrex;

extern "C"
{void fort_get_grav_const(Real* Gconst);}

using std::string;

void
Nyx::moveKickDriftExact (Real dt)
{
    std::cout  << "Doing exact moveKickDrift " << std::endl;

    // Find the current particle locations
    Vector<Real> part_locs;
    Nyx::theDMPC()->GetParticleLocations(part_locs);

    // Find the current particle masses
    Vector<Real> part_mass;
    int start_comp = 0;
    int   num_comp = 1;
    Nyx::theDMPC()->GetParticleData(part_mass,start_comp,num_comp);
 
    // Find the current particle velocity components 
    Vector<Real> part_vels;
    start_comp = 1;
      num_comp = BL_SPACEDIM;
    Nyx::theDMPC()->GetParticleData(part_vels,start_comp,num_comp);

    // Sanity check
    if (part_locs.size() != 2*BL_SPACEDIM)  
    {
        std::cout << "part_locs.size() is " << part_locs.size() << std::endl;
        amrex::Abort("moveKickDriftExact: we only call the exact solver for two particles");
    }

    // Sanity check
    if (part_vels.size() != 2*BL_SPACEDIM)  
    {
        std::cout << "part_vels.size() is " << part_vels.size() << std::endl;
        amrex::Abort("moveKickDriftExact: we only call the exact solver for two particles");
    }

    // These define the vector from the first to the second particle
    Real x = (part_locs[3]-part_locs[0]);
    Real y = (part_locs[4]-part_locs[1]);
    Real z = (part_locs[5]-part_locs[2]);

    // Distance between the particles
    Real r = std::sqrt( x*x + y*y + z*z); 

    // Gravitation acceleration = G m / r^2
    Real Gconst;
    fort_get_grav_const(&Gconst);

    // First particle 
    part_vels[0] += 0.5 * dt * Gconst * part_mass[1] * (x/std::pow(r,3));
    part_vels[1] += 0.5 * dt * Gconst * part_mass[1] * (y/std::pow(r,3));
    part_vels[2] += 0.5 * dt * Gconst * part_mass[1] * (z/std::pow(r,3));

    // Second particle 
    part_vels[3] += 0.5 * dt * Gconst * part_mass[0] * (-x/std::pow(r,3));
    part_vels[4] += 0.5 * dt * Gconst * part_mass[0] * (-y/std::pow(r,3));
    part_vels[5] += 0.5 * dt * Gconst * part_mass[0] * (-z/std::pow(r,3));

    // First particle -- position
    part_locs[0] += dt * part_vels[0];
    part_locs[1] += dt * part_vels[1];
    part_locs[2] += dt * part_vels[2];

    // Second particle -- position
    part_locs[3] += dt * part_vels[3];
    part_locs[4] += dt * part_vels[4];
    part_locs[5] += dt * part_vels[5];

    Nyx::theDMPC()->SetParticleLocations (part_locs);
    Nyx::theDMPC()->SetParticleVelocities(part_vels);
}

void
Nyx::moveKickExact (Real dt)
{
    std::cout  << "Doing exact moveKick " << std::endl;

    // Find the current particle locations
    Vector<Real> part_locs;
    Nyx::theDMPC()->GetParticleLocations(part_locs);

    // Find the current particle masses
    Vector<Real> part_mass;
    int start_comp = 0;
    int   num_comp = 1;
    Nyx::theDMPC()->GetParticleData(part_mass,start_comp,num_comp);

    // Find the current velocity components
    Vector<Real> part_vels;
    start_comp = 1;
      num_comp = BL_SPACEDIM;
    Nyx::theDMPC()->GetParticleData(part_vels,start_comp,num_comp);

    // Sanity check
    if (part_locs.size() != 2*BL_SPACEDIM)  
    {
        std::cout << "part_locs.size() is " << part_locs.size() << std::endl;
        amrex::Abort("moveKickExact: we only call the exact solver for two particles");
    }

    // Sanity check
    if (part_vels.size() != 2*BL_SPACEDIM)  
    {
        std::cout << "part_vels.size() is " << part_vels.size() << std::endl;
        amrex::Abort("moveKickExact: we only call the exact solver for two particles");
    }

    // Gravitation acceleration = G m / r^2 (mass included below)
    Real Gconst;
    fort_get_grav_const(&Gconst);

    // These define the vector from the first to the second particle
    Real x = (part_locs[3]-part_locs[0]);
    Real y = (part_locs[4]-part_locs[1]);
    Real z = (part_locs[5]-part_locs[2]);

    // Distance between the particles
    Real r = std::sqrt( x*x + y*y + z*z); 

    // First particle 
    part_vels[0] += 0.5 * dt * Gconst * part_mass[1] * (x/std::pow(r,3));
    part_vels[1] += 0.5 * dt * Gconst * part_mass[1] * (y/std::pow(r,3));
    part_vels[2] += 0.5 * dt * Gconst * part_mass[1] * (z/std::pow(r,3));

    // Second particle 
    part_vels[3] += 0.5 * dt * Gconst * part_mass[0] * (-x/std::pow(r,3));
    part_vels[4] += 0.5 * dt * Gconst * part_mass[0] * (-y/std::pow(r,3));
    part_vels[5] += 0.5 * dt * Gconst * part_mass[0] * (-z/std::pow(r,3));

    Nyx::theDMPC()->SetParticleVelocities(part_vels);
}
#endif
