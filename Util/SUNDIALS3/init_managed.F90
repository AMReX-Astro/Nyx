module init_managed

contains

attributes(host) subroutine init_allocations() &
     bind(C,name="fort_init_allocations")
use vode_aux_module
use atomic_rates_module
use meth_params_module, only: gamma_minus_1
allocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
    allocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
    allocate(uvb_density_A, uvb_density_B, mean_rhob)
! parameters need to be allocated or not?
!    allocate(MPROTON,XHYDROGEN,YHELIUM, BOLTZMANN)
    allocate(gamma_minus_1)
  TCOOLMIN = 0.0d0
  TCOOLMAX = 9.0d0
  TCOOLMIN_R = 10.0d0**TCOOLMIN 
  TCOOLMAX_R = 10.0d0**TCOOLMAX
  deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB
  uvb_density_A = 1.0d0
  uvb_density_B = 0.0d0
  gamma_minus_1 = 2.d0/3.d0

end subroutine init_allocations

attributes(host) subroutine init_tables_eos_params() &
     bind(C,name="fort_init_tables_eos_params")
use atomic_rates_module

call fort_tabulate_rates()

end subroutine init_tables_eos_params

attributes(host) subroutine fin_allocations() &
     bind(C,name="fort_fin_allocations")
use vode_aux_module
use atomic_rates_module
use meth_params_module, only: gamma_minus_1
deallocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
    deallocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
    deallocate(uvb_density_A, uvb_density_B, mean_rhob)
! parameters need to be allocated or not?
!    allocate(MPROTON,XHYDROGEN,YHELIUM, BOLTZMANN)
    deallocate(gamma_minus_1)

end subroutine fin_allocations
end module init_managed
