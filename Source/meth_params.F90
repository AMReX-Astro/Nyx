
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  use amrex_fort_module, only : rt => amrex_real
!  use cudafor

  implicit none

  integer , allocatable :: iorder        ! used only in uslope and uflaten

  real(rt), save, public :: gamma_const
  real(rt), allocatable, public  :: gamma_minus_1

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 5

  integer         , allocatable :: NVAR, NDIAG
  integer         , allocatable :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
  integer         , allocatable :: TEMP_COMP, NE_COMP, ZHI_COMP

  real(rt)        , allocatable :: small_dens, small_pres  
  real(rt)        , allocatable :: small_temp

<<<<<<< HEAD
=======
  integer,allocatable::ppm_type
  integer,allocatable::use_flattening
  integer,allocatable::use_const_species
>>>>>>> development
  integer,allocatable::heat_cool_type
  integer,allocatable::inhomo_reion
  integer,allocatable::grav_source_type
  integer,allocatable::cg_maxiter
  integer,allocatable::cg_tol
  integer,allocatable::cg_blend
  integer,allocatable::fix_mass_flux
  integer,allocatable::use_reset_state
  
#ifdef AMREX_USE_CUDA_FORTRAN
<<<<<<< HEAD
  attributes(managed) :: gamma_minus_1, iorder
  attributes(managed) :: URHO, UMX, UMY, UMZ, UEDEN, UEINT
  attributes(managed) :: TEMP_COMP, NE_COMP, ZHI_COMP, NVAR, NDIAG, small_temp, heat_cool_type
  attributes(managed) :: small_pres, small_dens
  attributes(managed) :: inhomo_reion,grav_source_type
  attributes(managed) :: cg_maxiter,cg_tol,cg_blend,fix_mass_flux,use_reset_state
  attributes(managed) :: use_reset_state
=======
  attributes(managed) :: gamma_minus_1, iorder!, gamma_const
  attributes(managed) :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  attributes(managed) :: TEMP_COMP, NE_COMP, ZHI_COMP, NTHERM, NVAR, NDIAG, small_temp, heat_cool_type
!  attributes(managed) :: NGDNV, GDPRES, GDU, GDV, GDW, QVAR
  attributes(managed) ::  nadv, small_pres, small_dens
  attributes(managed) ::  npassive, qpass_map, upass_map
  attributes(managed) ::  QTHERM, NQAUX, QVAR, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, NGDNV, QTEMP, QFX,  QGC
  attributes(managed) :: QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO, GDPRES, GDU, GDV, GDW
  attributes(managed) :: difmag
 attributes(managed) :: ppm_type,use_flattening,use_const_species,inhomo_reion,grav_source_type
 attributes(managed) :: cg_maxiter,cg_tol,cg_blend,fix_mass_flux,ppm_predict_gammae,ppm_temp_fix,&
      hybrid_riemann,riemann_solver,use_pslope,transverse_reset_density,transverse_reset_rhoe,use_pressure_law_pdivu, &
      use_analriem,use_srcQ_in_trace,use_csmall_gamma,use_reset_state,use_gamma_minus, use_area_dt_scale_apply
  attributes(managed) :: transverse_reset_density, transverse_reset_rhoe, ppm_predict_gammae, use_srcQ_in_trace, use_reset_state
>>>>>>> development
#endif

end module meth_params_module
