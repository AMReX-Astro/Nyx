
! This module stores the extra parameters for the VODE calls.

module vode_aux_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt), save :: z_vode
  real(rt), save :: rho_vode, T_vode, ne_vode
  integer , save :: i_vode, j_vode, k_vode
  !$OMP THREADPRIVATE (rho_vode, T_vode, ne_vode, i_vode, j_vode, k_vode)

end module vode_aux_module
