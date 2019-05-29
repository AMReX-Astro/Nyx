
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  use amrex_fort_module, only : rt => amrex_real
!  use cudafor

  implicit none

  real(rt), allocatable :: difmag        ! used only in consup to weight the divu contributin
  integer , allocatable :: iorder        ! used only in uslope and uflaten

  real(rt), save, public :: gamma_const
  real(rt), allocatable, public  :: gamma_minus_1

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 5

  ! NTHERM: number of thermodynamic variables
  integer         , allocatable :: NTHERM, NVAR, NDIAG
  integer         , allocatable :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  integer         , allocatable :: TEMP_COMP, NE_COMP, ZHI_COMP
!  integer         , allocatable :: NGDNV, GDPRES, GDU, GDV, GDW, QVAR

  ! QTHERM: number of primitive variables
  integer         , allocatable :: QTHERM, NQAUX, QVAR, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, NGDNV, QTEMP, QFX,  QGC
  integer         , allocatable :: QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO, GDPRES, GDU, GDV, GDW
  
  integer         , allocatable :: nadv

  real(rt)        , allocatable :: small_dens, small_pres  
  real(rt)        , allocatable :: small_temp

  integer,allocatable::ppm_type
  integer,allocatable::ppm_reference
  integer,allocatable::ppm_flatten_before_integrals
  integer,allocatable::use_colglaz
  integer,allocatable::use_flattening
  integer,allocatable::version_2
  integer,allocatable::use_const_species
  integer,allocatable::normalize_species
  integer,allocatable::heat_cool_type
  integer,allocatable::inhomo_reion
  integer,allocatable::grav_source_type
  integer,allocatable::cg_maxiter
  integer,allocatable::cg_tol
  integer,allocatable::cg_blend
  integer,allocatable::fix_mass_flux
  integer,allocatable::ppm_predict_gammae
  integer,allocatable::ppm_temp_fix
  integer,allocatable::ppm_reference_eigenvectors
  integer,allocatable::hybrid_riemann
  integer,allocatable::riemann_solver
  integer,allocatable::use_pslope
  integer,allocatable::transverse_reset_density
  integer,allocatable::transverse_reset_rhoe
  integer,allocatable::use_pressure_law_pdivu
  integer,allocatable::use_analriem
  integer,allocatable::use_srcQ_in_trace
  integer,allocatable::use_csmall_gamma
  integer,allocatable::use_reset_state
  integer,allocatable::use_gamma_minus
  integer         , allocatable :: use_area_dt_scale_apply
  
  integer, allocatable :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)
#ifdef AMREX_USE_CUDA
  attributes(managed) :: gamma_minus_1, iorder!, gamma_const
  attributes(managed) :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  attributes(managed) :: TEMP_COMP, NE_COMP, ZHI_COMP, NTHERM, NVAR, NDIAG, small_temp, heat_cool_type
!  attributes(managed) :: NGDNV, GDPRES, GDU, GDV, GDW, QVAR
  attributes(managed) ::  nadv, small_pres, small_dens
  attributes(managed) ::  npassive, qpass_map, upass_map
  attributes(managed) ::  QTHERM, NQAUX, QVAR, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, NGDNV, QTEMP, QFX,  QGC
  attributes(managed) :: QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO, GDPRES, GDU, GDV, GDW
  attributes(managed) :: difmag
 attributes(managed) :: ppm_type,ppm_reference,ppm_flatten_before_integrals,use_colglaz,use_flattening,version_2,use_const_species,normalize_species,inhomo_reion,grav_source_type
 attributes(managed) :: cg_maxiter,cg_tol,cg_blend,fix_mass_flux,ppm_predict_gammae,ppm_temp_fix,ppm_reference_eigenvectors, &
      hybrid_riemann,riemann_solver,use_pslope,transverse_reset_density,transverse_reset_rhoe,use_pressure_law_pdivu, &
      use_analriem,use_srcQ_in_trace,use_csmall_gamma,use_reset_state,use_gamma_minus, use_area_dt_scale_apply
  attributes(managed) :: transverse_reset_density, transverse_reset_rhoe, ppm_predict_gammae, use_srcQ_in_trace, use_reset_state
#endif

end module meth_params_module
