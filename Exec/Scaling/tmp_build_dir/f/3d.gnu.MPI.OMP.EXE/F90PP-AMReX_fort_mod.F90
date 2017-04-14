












module amrex_fort_module

  use iso_c_binding, only : c_float, c_double

  implicit none

  integer, parameter ::    bl_spacedim = 3
  integer, parameter :: amrex_spacedim = 3

  integer, parameter :: amrex_real = c_double

end module amrex_fort_module
