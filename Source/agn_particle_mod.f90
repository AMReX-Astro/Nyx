module particle_mod

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t, ghost_t
  
  type, bind(C)  :: particle_t
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: mass       !< Particle mass
     real(c_real)    :: vel(3)     !< Particle velocity
     real(c_real)    :: energy     !< Particle energy
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type particle_t
  
  type, bind(C)  :: ghost_t
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: mass       !< Particle mass
     real(c_real)    :: vel(3)     !< Particle velocity
     real(c_real)    :: energy     !< Particle energy
  end type ghost_t


end module
