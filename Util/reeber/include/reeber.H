#ifndef Halo_H_
#define Halo_H_

struct Halo {
  long long int id;
  amrex::Real total_mass;
  amrex::Real dm_mass;
  amrex::Real gas_mass;
  amrex::Real n_cells;  // can be fractional in multi-level case
  amrex::Real volume;
  amrex::IntVect position;
  amrex::Real real_position[AMREX_SPACEDIM];
  amrex::Real gas_vel_x;
  amrex::Real gas_vel_y;
  amrex::Real gas_vel_z;
};


#endif
