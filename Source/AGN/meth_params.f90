module meth_params_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none
  integer         , save :: NTHERM, NVAR, NDIAG
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  integer         , save :: TEMP_COMP, NE_COMP, ZHI_COMP

end module meth_params_module
