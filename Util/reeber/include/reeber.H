#ifndef Halo_H_
#define Halo_H_

struct Halo {
  long long int id;
  amrex::Real total_mass;
  amrex::Real dm_mass;
  amrex::Real gas_mass;
  long long int int_volume;
  amrex::Real volume;
  amrex::IntVect position;
  amrex::Real real_position[AMREX_SPACEDIM];
};


#endif
