
  subroutine fill_slice_3d(full_data, flo, fhi, fstart, nfull, slice_data, slo, shi, tlo, thi, ncomp) &
       bind(C, name="fill_slice")

    use amrex_fort_module, only : amrex_real

    integer       ,   intent(in)    :: ncomp, fstart, nfull
    integer       ,   intent(in)    :: flo(3), fhi(3)
    integer       ,   intent(in)    :: slo(3), shi(3)
    integer       ,   intent(in)    :: tlo(3), thi(3)
    real(amrex_real), intent(inout) ::  full_data(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3), nfull)
    real(amrex_real), intent(inout) :: slice_data(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3), ncomp)

    integer n, i, j, k

    do n = 1, ncomp
       do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
             do i = tlo(1), thi(1)
                slice_data(i, j, k, n) = full_data(i, j, k, fstart+n)
             end do
          end do
       end do
    end do

  end subroutine fill_slice_3d
