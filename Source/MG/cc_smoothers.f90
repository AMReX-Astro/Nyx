module cc_smoothers_module

  use bl_constants_module
  use cc_stencil_module

  implicit none

contains

  subroutine gs_rb_smoother_3d(omega, ss, uu, ff, mm, lo, ng, n, skwd)
    use bl_prof_module
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in) :: ss(0:,lo(1):, lo(2):, lo(3):)
    integer            ,intent(in) :: mm(lo(1):,lo(2):,lo(3):)
    logical, intent(in), optional :: skwd
    integer :: i, j, k, ioff
    integer :: hi(size(lo))
    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd
    real(dp_t) :: dd, dhsq_inv, ss0, ss0_inv

    type(bl_prof_timer), save :: bpt

    call build(bpt, "gs_rb_smoother_3d")

    hi = ubound(ff)

    ss0_inv  =  1.d0/ss(0,lo(1),lo(2),lo(3))
    dhsq_inv = -ss(0,lo(1),lo(2),lo(3))/6.d0
        
    !$OMP PARALLEL DO PRIVATE(k,j,i,ioff,dd) IF((hi(3)-lo(3)).ge.3)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          ioff = 0; if ( mod (lo(1) + j + k, 2) /= n ) ioff = 1
          do i = lo(1)+ioff, hi(1), 2

             dd = (-6.d0*uu(i,j,k)   + &
                   uu(i+1,j,k) + uu(i-1,j,k) + &
                   uu(i,j+1,k) + uu(i,j-1,k) + &
                   uu(i,j,k+1) + uu(i,j,k-1) ) * dhsq_inv
 
             uu(i,j,k) = uu(i,j,k) + (ff(i,j,k) - dd)*ss0_inv
!            uu(i,j,k) = uu(i,j,k) + (ff(i,j,k) - dd)* 1.d0/ss(0,i,j,k)

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine gs_rb_smoother_3d

end module cc_smoothers_module
