module probdata_module

!     These determine the refinement criteria
      double precision, save :: denerr,  dengrad
      double precision, save :: presserr,pressgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev

!     Sod variables
      double precision, save ::  p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, frac
      double precision, save :: center(3)

!     These help specify which specific problem
      integer        , save ::  probtype,idir

end module probdata_module
