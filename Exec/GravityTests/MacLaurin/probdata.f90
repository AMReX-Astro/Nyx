module probdata_module


      use amrex_fort_module, only : rt => amrex_real
!     Tagging variables
      integer, save :: max_num_part

!     Use these to define the spheroid
      real(rt), save :: center(3)
      real(rt), save :: a1, a3
      
end module probdata_module
