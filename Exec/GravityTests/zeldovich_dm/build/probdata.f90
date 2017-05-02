module probdata_module

    use amrex_fort_module, only : rt => amrex_real

!   Tagging variables
    integer, save :: max_num_part

!   Subcluster variables
    real(rt), save :: rho_c, r_c

!   Residual variables
    real(rt), save :: center(3)

end module probdata_module
