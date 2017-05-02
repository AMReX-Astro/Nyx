module mg_smoother_module
 
  use amrex_fort_module, only : rt => amrex_real
  use multifab_module
  use cc_stencil_module
  use mg_tower_module
  use bl_timer_module
 
  implicit none

contains

  subroutine mg_tower_smoother(mgt, lev, ss, uu, ff, mm)

    use bl_prof_module
    use cc_smoothers_module, only: gs_rb_smoother_3d

    integer        , intent(in   ) :: lev
    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: ff
    type( multifab), intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm

    real(rt), pointer :: fp(:,:,:,:)
    real(rt), pointer :: up(:,:,:,:)
    real(rt), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, k, n, ng, nn, stat, npts
    integer :: lo(mgt%dim)
    type(bl_prof_timer), save :: bpt
    logical :: pmask(mgt%dim), singular_test
    real(rt) :: local_eps

    if (.not.nodal_q(ff)) then
       singular_test =  mgt%bottom_singular .and. mgt%coeffs_sum_to_zero
    end if

    pmask = get_pmask(get_layout(uu))
 
    ! Make sure to define this here so we don't assume a certain number of ghost cells for uu
    ng = nghost(uu)

    call build(bpt, "mgt_smoother")

    if (mgt%skewed_not_set(lev)) then 
       do i = 1, nfabs(mm)
          mp => dataptr(mm, i)
          mgt%skewed(lev,i) = skewed_q(mp)
       end do
       mgt%skewed_not_set(lev) = .false.
    end if

    do nn = 0, 1
       call multifab_fill_boundary(uu, cross = mgt%lcross)

       do i = 1, nfabs(ff)
          up => dataptr(uu, i)
          fp => dataptr(ff, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          lo =  lwb(get_box(ss, i))
          call gs_rb_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,1), &
                                 fp(:,:,:,1), mp(:,:,:,1), lo, ng, nn, &
                                 mgt%skewed(lev,i))
       end do

    end do

    call destroy(bpt)

  end subroutine mg_tower_smoother

end module mg_smoother_module
