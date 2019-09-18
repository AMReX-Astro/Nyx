module init_managed
use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
contains

subroutine init_allocations() &
     bind(C,name="fort_alloc_cuda_managed")
#ifdef HEATCOOL
  use vode_aux_module
#endif
use atomic_rates_module
#ifdef HEATCOOL
use eos_module, only: xacc, vode_rtol, vode_atol_scaled
#endif
use meth_params_module
use reion_aux_module

allocate(this_z)
!! meth_params
allocate(gamma_minus_1,iorder)!,gamma_const)
allocate(URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX)
allocate(TEMP_COMP, NE_COMP, ZHI_COMP, NTHERM, NVAR, NDIAG, small_temp, heat_cool_type)

allocate(nadv, small_pres, small_dens)
allocate(npassive)
allocate(QTHERM, NQAUX, QVAR, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, NGDNV, QTEMP, QFX,  QGC)
allocate(QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO, GDPRES, GDU, GDV, GDW)

allocate(difmag)
!allocate(NGDNV, GDPRES, GDU, GDV, GDW,QVAR,use_pressure_law_pdivu,use_area_dt_scale_apply)
allocate(ppm_type,ppm_reference,ppm_flatten_before_integrals,use_colglaz,use_flattening,version_2,use_const_species,normalize_species,inhomo_reion,grav_source_type)
allocate(cg_maxiter,cg_tol,cg_blend,fix_mass_flux,ppm_predict_gammae,ppm_temp_fix,ppm_reference_eigenvectors, &
     hybrid_riemann,riemann_solver,use_pslope,transverse_reset_density,transverse_reset_rhoe,use_pressure_law_pdivu, &
     use_analriem,use_srcQ_in_trace,use_csmall_gamma,use_reset_state,use_gamma_minus,use_area_dt_scale_apply)
cg_maxiter=12
cg_tol=1.0d-5
cg_blend=0
fix_mass_flux=0
ppm_predict_gammae=1
ppm_temp_fix=0
ppm_reference_eigenvectors=0
hybrid_riemann=0
riemann_solver=0
use_pslope=0
transverse_reset_density=1
transverse_reset_rhoe=0
use_pressure_law_pdivu=1
use_analriem=0
use_srcQ_in_trace=0
use_csmall_gamma=1
use_reset_state=0
use_gamma_minus=1
use_area_dt_scale_apply = 1
!!eos
allocate(XHYDROGEN, YHELIUM)
#ifdef HEATCOOL
allocate(xacc,vode_rtol,vode_atol_scaled)
!! vode
allocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
#endif
!! atomic
allocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
allocate(uvb_density_A, uvb_density_B, mean_rhob)
!!Reion
allocate(zhi_flash, zheii_flash, T_zhi, T_zheii)
allocate(flash_h, flash_he, inhomogeneous_on)
TCOOLMIN = 0.0d0
TCOOLMAX = 9.0d0
TCOOLMIN_R = 10.0d0**TCOOLMIN 
TCOOLMAX_R = 10.0d0**TCOOLMAX
deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB
uvb_density_A = 1.0d0
uvb_density_B = 0.0d0
XHYDROGEN = 0.76d0
YHELIUM   = 7.8947368421d-2
gamma_minus_1 = 2.d0/3.d0

!reion
zhi_flash=-1.0
zheii_flash=-1.0
T_zhi=0.0
T_zheii=0.0
flash_h=.false.
flash_he=.false.
inhomogeneous_on=.false.
print*,xacc
end subroutine init_allocations

subroutine init_tables_eos_params() &
     bind(C,name="fort_init_tables_eos_params")
  use amrex_constants_module, only : rt => amrex_real, M_PI
  use atomic_rates_module
  use fundamental_constants_module, only: Gconst, mp_over_kb
  use comoving_module, only: comoving_h,comoving_OmB
  use comoving_nd_module, only: fort_integrate_comoving_a
  use reion_aux_module, only: zhi_flash, zheii_flash, T_zhi, T_zheii, &
       flash_h, flash_he, inhomogeneous_on
#ifdef HEATCOOL
  use eos_module, only: fort_setup_eos_params
  use vode_aux_module, only: fn_vode, NR_vode, z_vode, JH_vode, JHe_vode
#endif
  implicit none
  real(rt) :: vode_atol_scaled_in, vode_rtol_in, xacc_in
  integer :: simd_width
  real(rt) :: a, half_dt
  integer :: i, j, k
  real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
  
  call fort_tabulate_rates()

#ifdef HEATCOOL
    simd_width = 1    
    vode_atol_scaled_in = 1e-4
    vode_rtol_in = 1e-4
    xacc_in = 1e-6

    print*,"start reading table"

    call fort_setup_eos_params(xacc_in, vode_rtol_in, vode_atol_scaled_in)

    print*,"Finished reading table"
    fn_vode = 0
    NR_vode = 0
!    print*,"Read parameters"

!     FMT="(A6,I1,/,ES21.15,/,ES21.15E2,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15)"
!    read(1,FMT) string, STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig

!    print*,"Finished reading parameters:"
!    print(FMT), string,STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig
    a=0.0688707121 !1.635780036449432E-01

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
#endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

end subroutine init_tables_eos_params

subroutine fin_allocations() &
     bind(C,name="fort_dealloc_cuda_managed")
#ifdef HEATCOOL
  use vode_aux_module
  use eos_module, only: xacc, vode_rtol, vode_atol_scaled
  use atomic_rates_module
#endif
use meth_params_module
use reion_aux_module

#ifdef HEATCOOL
deallocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
deallocate(xacc,vode_rtol,vode_atol_scaled)
#endif
#ifdef HEATCOOL
deallocate(this_z)
#endif
deallocate(gamma_minus_1,iorder)!, gamma_const)
deallocate(URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX)
deallocate(TEMP_COMP, NE_COMP, ZHI_COMP, NTHERM, NVAR, NDIAG, small_temp, heat_cool_type)

deallocate(nadv, small_pres, small_dens)
deallocate(npassive)
deallocate(QTHERM, NQAUX, QVAR, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, NGDNV, QTEMP, QFX,  QGC)
deallocate(QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO, GDPRES, GDU, GDV, GDW)

deallocate(difmag)
!deallocate(NGDNV, GDPRES, GDU, GDV, GDW, QVAR,use_pressure_law_pdivu,use_area_dt_scale_apply)
deallocate(ppm_type,ppm_reference,ppm_flatten_before_integrals,use_colglaz,use_flattening,version_2, &
     use_const_species,normalize_species,inhomo_reion,grav_source_type)
deallocate(cg_maxiter,cg_tol,cg_blend,fix_mass_flux,ppm_predict_gammae,ppm_temp_fix,ppm_reference_eigenvectors, &
     hybrid_riemann,riemann_solver,use_pslope,transverse_reset_density,transverse_reset_rhoe,use_pressure_law_pdivu, &
     use_analriem,use_srcQ_in_trace,use_csmall_gamma,use_reset_state,use_gamma_minus,use_area_dt_scale_apply)
#ifdef HEATCOOL
deallocate(XHYDROGEN,YHELIUM)
deallocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
deallocate(uvb_density_A, uvb_density_B, mean_rhob)
deallocate(zhi_flash, zheii_flash, T_zhi, T_zheii)
deallocate(flash_h, flash_he, inhomogeneous_on)

! Initially allocated when the table is read
if(allocated(ggh0)) deallocate(ggh0)
if(allocated(gghe0)) deallocate(gghe0)
if(allocated(gghep)) deallocate(gghep)
if(allocated(eh0)) deallocate(eh0)
if(allocated(ehe0)) deallocate(ehe0)
if(allocated(ehep)) deallocate(ehep)
if(allocated(lzr)) deallocate( lzr)
if(allocated(rggh0)) deallocate(rggh0)
if(allocated(rgghe0)) deallocate( rgghe0)
if(allocated(rgghep)) deallocate( rgghep)

if(allocated(reh0)) deallocate(reh0)
if(allocated(rehe0)) deallocate(rehe0)
if(allocated(rehep)) deallocate(rehep)
if(allocated(NCOOLFILE)) deallocate(NCOOLFILE)
if(allocated(AlphaHp)) deallocate(AlphaHp)
if(allocated(AlphaHep)) deallocate(AlphaHep)
if(allocated(AlphaHepp)) deallocate(AlphaHepp)
if(allocated(Alphad)) deallocate(Alphad)
if(allocated(GammaeH0)) deallocate(GammaeH0)
if(allocated(GammaeHe0)) deallocate(GammaeHe0)
if(allocated(GammaeHep)) deallocate(GammaeHep)
if(allocated(BetaH0)) deallocate(BetaH0)
if(allocated(BetaHe0)) deallocate(BetaHe0)
if(allocated(BetaHep)) deallocate(BetaHep)
if(allocated(Betaff1)) deallocate(Betaff1)
if(allocated(Betaff4)) deallocate(Betaff4)
if(allocated(RecHp)) deallocate(RecHp)
if(allocated(RecHep)) deallocate(RecHep)
if(allocated(RecHepp)) deallocate(RecHepp)
#endif
end subroutine fin_allocations

!!! subroutine dummy()
!!! allocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
!!!     allocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
!!!     allocate(uvb_density_A, uvb_density_B, mean_rhob)
!!! ! parameters need to be allocated or not?
!!!     allocate(XHYDROGEN,YHELIUM)
!!! !    allocate(gamma_minus_1)
!!!   TCOOLMIN = 0.0d0
!!!   TCOOLMAX = 9.0d0
!!!   TCOOLMIN_R = 10.0d0**TCOOLMIN 
!!!   TCOOLMAX_R = 10.0d0**TCOOLMAX
!!!   deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB
!!!   uvb_density_A = 1.0d0
!!!   uvb_density_B = 0.0d0
!!!   XHYDROGEN = 0.76d0
!!!   YHELIUM   = 7.8947368421d-2  ! (1.0d0-XHYDROGEN)/(4.0d0*XHYDROGEN)
!!! !  gamma_minus_1 = 2.d0/3.d0
!!! deallocate(z_vode, rho_vode, T_vode, ne_vode, JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, firstcall)
!!!     deallocate(TCOOLMIN, TCOOLMAX, TCOOLMAX_R, TCOOLMIN_R, deltaT)
!!!     deallocate(uvb_density_A, uvb_density_B, mean_rhob)
!!! ! parameters need to be allocated or not?
!!!     deallocate(XHYDROGEN,YHELIUM)
!!!     deallocate(gamma_minus_1)
!!! end subroutine dummy
end module init_managed
