
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
  integer         , allocatable :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  integer         , allocatable :: TEMP_COMP, NE_COMP, ZHI_COMP

  integer         , allocatable :: QC, NQ, UTEMP, QTEMP, QFX,  QGC
  
  integer         , allocatable :: nadv

  real(rt)        , allocatable :: small_dens, small_pres  
  real(rt)        , allocatable :: small_temp

  integer,allocatable::use_const_species
  integer,allocatable::normalize_species
  integer,allocatable::heat_cool_type
  integer,allocatable::inhomo_reion
  integer,allocatable::grav_source_type
  integer,allocatable::cg_maxiter
  integer,allocatable::cg_tol
  integer,allocatable::cg_blend
  integer,allocatable::fix_mass_flux
  integer,allocatable::use_reset_state
  
#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(managed) :: gamma_minus_1, iorder
  attributes(managed) :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX
  attributes(managed) :: TEMP_COMP, NE_COMP, ZHI_COMP, NVAR, NDIAG, small_temp, heat_cool_type
  attributes(managed) :: nadv, small_pres, small_dens
  attributes(managed) :: QC, NQ, UTEMP, QTEMP, QFX,  QGC
  attributes(managed) :: use_const_species,normalize_species,inhomo_reion,grav_source_type
  attributes(managed) :: cg_maxiter,cg_tol,cg_blend,fix_mass_flux,use_reset_state
  attributes(managed) :: use_reset_state
#endif

end module meth_params_module
