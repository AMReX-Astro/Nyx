
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  use amrex_fort_module, only : rt => amrex_real
!  use cudafor

  implicit none

  real(rt), save :: difmag        ! used only in consup to weight the divu contributin
  integer , save :: iorder        ! used only in uslope and uflaten

  real(rt), save, public  :: gamma_const
  real(rt), allocatable, public  :: gamma_minus_1

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 5

  ! NTHERM: number of thermodynamic variables
  integer         , allocatable :: NTHERM, NVAR, NDIAG
  integer         , allocatable :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  integer         , allocatable :: TEMP_COMP, NE_COMP, ZHI_COMP
  integer         , allocatable :: NGDNV, GDPRES, GDU, GDV, GDW, QVAR

  ! QTHERM: number of primitive variables
  integer         , save :: QTHERM, NQAUX, QC, NQSRC, NQ, UTEMP, QGAME, QGAMC, &!NGDNV, 
QTEMP, QFX,  QGC
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, GDGAME, GDRHO !GDPRES, GDRHO, GDU, GDV, GDW
  
  integer         , save :: nadv

  real(rt)        , save :: small_dens, small_pres  
  real(rt)        , allocatable :: small_temp

  integer         , save :: ppm_type
  integer         , save :: ppm_reference
  integer         , save :: ppm_flatten_before_integrals
  integer         , save :: use_colglaz
  integer         , save :: use_flattening
  integer         , save :: version_2
  integer         , save :: use_const_species
  integer         , save :: normalize_species
  integer         , allocatable :: heat_cool_type
  integer         , save :: inhomo_reion
  integer         , save :: grav_source_type
  integer         , save :: cg_maxiter = 12
  integer         , save :: cg_tol = 1.0d-5
  integer         , save :: cg_blend = 0
  integer         , save :: fix_mass_flux = 0
  !ppm ppm_predict_gammae=0 matches castro
  integer         , save :: ppm_predict_gammae = 1
  integer         , save :: ppm_temp_fix = 0
  integer         , save :: ppm_reference_eigenvectors = 0
  integer         , save :: hybrid_riemann = 0
  integer         , save :: riemann_solver = 0
  integer         , save :: use_pslope = 0
  integer         , save :: transverse_reset_density = 1
  integer         , save :: transverse_reset_rhoe = 0
  !use_pressure_law_pdivu = 1 matches old Nyx
  integer         , allocatable :: use_pressure_law_pdivu
  !use_analriem = 0 matches castro
  integer         , save :: use_analriem = 0
  !use_srcQ_in_trace = 1 matches castro
  integer         , save :: use_srcQ_in_trace = 0
  integer         , save :: use_csmall_gamma = 1
  integer         , save :: use_reset_state = 0
  integer         , save :: use_gamma_minus = 1
  integer         , allocatable :: use_area_dt_scale_apply
  
  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)
#ifdef AMREX_USE_CUDA
  attributes(managed) :: gamma_minus_1
  attributes(managed) :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  attributes(managed) :: TEMP_COMP, NE_COMP, ZHI_COMP, NTHERM, NVAR, NDIAG, small_temp, heat_cool_type
  attributes(managed) :: NGDNV, GDPRES, GDU, GDV, GDW, QVAR
  attributes(managed) :: use_pressure_law_pdivu, use_area_dt_scale_apply
#endif

end module meth_params_module
