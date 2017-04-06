module probdata_module

!     These determine the refinement criteria
      double precision, save ::  denerr, dengrad
      integer         , save ::  max_denerr_lev, max_dengrad_lev

      integer         , save ::  radiative_cooling_type

      double precision, save ::  alpha, rho0, temp0


!     Use these in add_turb_forcing.
      double precision, save ::  prob_lo(3), prob_hi(3)
      
!because of harald... ;)
      double precision, save :: center(3)   
end module probdata_module
