module reion_aux_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! Global variables (re)set on inputs
  real(rt), allocatable :: zhi_flash, zheii_flash, T_zhi, T_zheii
  logical, allocatable  :: flash_h, flash_he, inhomogeneous_on
#ifdef AMREX_USE_CUDA
  attributes(managed) :: zhi_flash, zheii_flash, T_zhi, T_zheii
  attributes(managed)  :: flash_h, flash_he, inhomogeneous_on
#endif

end module reion_aux_module
