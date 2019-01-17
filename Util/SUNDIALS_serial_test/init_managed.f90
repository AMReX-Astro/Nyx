module init_managed

contains

subroutine init_allocations() &
     bind(C,name="fort_init_allocations")
use vode_aux_module
use atomic_rates_module
!!allocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
!!    allocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
!!    allocate(uvb_density_A, uvb_density_B, mean_rhob)
! parameters need to be allocated or not?
!    allocate(MPROTON,XHYDROGEN,YHELIUM, BOLTZMANN)
!!    allocate(gamma_minus_1)
!!  TCOOLMIN = 0.0d0
!!  TCOOLMAX = 9.0d0
!!  TCOOLMIN_R = 10.0d0**TCOOLMIN 
!!  TCOOLMAX_R = 10.0d0**TCOOLMAX
!!  deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB
!!  uvb_density_A = 1.0d0
!!  uvb_density_B = 0.0d0

end subroutine init_allocations

subroutine init_tables_eos_params(a) &
     bind(C,name="fort_init_tables_eos_params")
  use constants_module, only : rt => type_real, M_PI
  use atomic_rates_module
  use eos_module, only: fort_setup_eos_params
  use fundamental_constants_module, only: Gconst, mp_over_kb
  use comoving_module, only: comoving_h,comoving_OmB
  use comoving_nd_module, only: fort_integrate_comoving_a
  use reion_aux_module, only: zhi_flash, zheii_flash, T_zhi, T_zheii, &
       flash_h, flash_he, inhomogeneous_on
  use vode_aux_module, only: fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
  use meth_params_module, only: gamma_minus_1
implicit none
  real(c_double) :: vode_atol_scaled_in, vode_rtol_in, xacc_in
  integer :: simd_width
  real(rt), intent(inout) :: a
  real(rt) :: half_dt
  integer :: i, j, k
  real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z

  gamma_minus_1 = 2.d0/3.d0
  
  call fort_tabulate_rates()

    simd_width = 1    
    vode_atol_scaled_in = 1e-4
    vode_rtol_in = 1e-4
    xacc_in = 1e-6

    call fort_setup_eos_params(xacc_in, vode_rtol_in, vode_atol_scaled_in)

    print*,"Finished reading table"

    fn_vode = 0
    NR_vode = 0
!    print*,"Read parameters"

!     FMT="(A6,I1,/,ES21.15,/,ES21.15E2,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15)"
!    read(1,FMT) string, STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig

!    print*,"Finished reading parameters:"
!    print(FMT), string,STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig
!    a=1.635780036449432E-01

    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, half_dt)
    z_end = 1.0d0/a_end - 1.0d0
    !Added z_vode arbitrarily to be z, since it was set to 0
    call fort_interp_to_this_z(z)
    z_vode = z

    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    call SetCellInit()
end subroutine init_tables_eos_params

!!!!!!!!!!!!!!!!!!!! ASSUME inhomogeneous_on .eq. .false.
subroutine SetCellInit()

  use vode_aux_module, only: fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
print*,   fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
!                if (inhomogeneous_on) then
!                   H_reion_z = 1*H_reion_z!diag_eos(i,j,k,ZHI_COMP)
!                   if (z .gt. H_reion_z) then
!                      JH_vode = 0
!                   else
!                      JH_vode = 1
!                   endif
!                endif

!                if (e_orig .lt. 0.d0) then
!!!                    !$OMP CRITICAL
!                    print *,'negative e entering strang integration ',z, i,j,k, rho/mean_rhob, e_orig
!                    call bl_abort('bad e in strang')
!!!                    !$OMP END CRITICAL
!                end if
end subroutine SetCellInit

!!!!!!!!!!!!!!!!!!!! ASSUME inhomogeneous_on .eq. .false.
subroutine SetCellFin()

  use vode_aux_module, only: fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
print*,   fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
!                if (inhomogeneous_on) then
!                   H_reion_z = 1*H_reion_z!diag_eos(i,j,k,ZHI_COMP)
!                   if (z .gt. H_reion_z) then
!                      JH_vode = 0
!                   else
!                      JH_vode = 1
!                   endif
!                endif

!                if (e_orig .lt. 0.d0) then
!!!                    !$OMP CRITICAL
!                    print *,'negative e entering strang integration ',z, i,j,k, rho/mean_rhob, e_orig
!                    call bl_abort('bad e in strang')
!!!                    !$OMP END CRITICAL
!                end if
end subroutine SetCellFin

subroutine fin_allocations() &
     bind(C,name="fort_fin_allocations")
use vode_aux_module
use atomic_rates_module
use meth_params_module, only: gamma_minus_1
!!deallocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
!!    deallocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
!!    deallocate(uvb_density_A, uvb_density_B, mean_rhob)
! parameters need to be allocated or not?
!    allocate(MPROTON,XHYDROGEN,YHELIUM, BOLTZMANN)
!!    deallocate(gamma_minus_1)

end subroutine fin_allocations
end module init_managed
