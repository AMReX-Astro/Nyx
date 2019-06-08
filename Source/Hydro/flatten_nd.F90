module ca_flatten_module

  use amrex_constants_module, only : ZERO
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public :: ca_uflatten, uflaten

contains

  AMREX_CUDA_FORT_DEVICE subroutine ca_uflatten(lo, hi, &
                         q, q_lo, q_hi, &
                         flatn, f_lo, f_hi, pres_comp) bind(c,name='ca_uflatten')

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: QVAR, QU, QV, QW

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer, intent(in), value :: pres_comp

    integer :: i, j, k, ishft

    real(rt), parameter :: small_pres = 1.e-200_rt

    real(rt) :: denom, zeta, tst, tmp, ftmp
    real(rt) :: dp, z, z2, chi, chi2

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    !$gpu

    ! x-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i+1,j,k,pres_comp) - q(i-1,j,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i+2,j,k,pres_comp)-q(i-2,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1,j,k,QU)-q(i+1,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1,j,k,pres_comp), q(i-1,j,k,pres_comp))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i+1-ishft,j,k,pres_comp) - q(i-1-ishft,j,k,pres_comp)

             denom = max(small_pres, abs(q(i+2-ishft,j,k,pres_comp)-q(i-2-ishft,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1-ishft,j,k,QU)-q(i+1-ishft,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1-ishft,j,k,pres_comp), q(i-1-ishft,j,k,pres_comp))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             flatn(i,j,k) = ONE - max(chi2 * z2, chi * z)

          end do
       end do
    end do

#if AMREX_SPACEDIM >= 2
    ! y-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i,j+1,k,pres_comp) - q(i,j-1,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j+2,k,pres_comp)-q(i,j-2,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1,k,QV)-q(i,j+1,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1,k,pres_comp), q(i,j-1,k,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j+1-ishft,k,pres_comp) - q(i,j-1-ishft,k,pres_comp)

             denom = max(small_pres, abs(q(i,j+2-ishft,k,pres_comp)-q(i,j-2-ishft,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1-ishft,k,QV)-q(i,j+1-ishft,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1-ishft,k,pres_comp), q(i,j-1-ishft,k,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             ftmp = ONE - max(chi2 * z2, chi * z)
             flatn(i,j,k) = min(flatn(i,j,k), ftmp)

          end do
       end do
    end do
#endif

#if AMREX_SPACEDIM == 3
    ! z-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i,j,k+1,pres_comp) - q(i,j,k-1,pres_comp)

             if(dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j,k+2,pres_comp)-q(i,j,k-2,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1,QW)-q(i,j,k+1,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1,pres_comp),q(i,j,k-1,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j,k+1-ishft,pres_comp) - q(i,j,k-1-ishft,pres_comp)

             denom = max(small_pres, abs(q(i,j,k+2-ishft,pres_comp)-q(i,j,k-2-ishft,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1-ishft,QW)-q(i,j,k+1-ishft,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1-ishft,pres_comp),q(i,j,k-1-ishft,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             ftmp = ONE - max(chi2 * z2, chi * z)
             flatn(i,j,k) = min(flatn(i,j,k), ftmp)
          enddo
       enddo
    enddo
#endif

  end subroutine ca_uflatten

  ! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,w,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : iorder, small_pres
    use amrex_constants_module

    implicit none

    integer lo(3),hi(3)
    integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    
    real(rt) p(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    real(rt) u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    real(rt) v(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    real(rt) w(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    real(rt) flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

    integer i, j, k, idx, ishft
    integer nx,ny,nz,nmax

    real(rt) denom, zeta, tst, tmp, ftmp

    ! Local arrays
    real(rt), allocatable :: dp(:,:,:), z(:,:,:), chi(:,:,:)
    
    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33d0, zcut1 = 0.75d0, zcut2 = 0.85d0, dzcut = ONE/(zcut2-zcut1)

    nx = hi(1)-lo(1)+3
    ny = hi(2)-lo(2)+3
    nz = hi(3)-lo(3)+3

    nmax = max(nx,ny,nz)

    if (iorder .eq. 3) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                flatn(i,j,k) = ONE
             enddo
          enddo
       enddo
       return
    endif

    ! x-direction flattening coef
    allocate(dp (0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
    allocate(z  (0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
    allocate(chi(0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
    do k = lo(3),hi(3)
       do j = lo(2),hi(2) 
          do i = lo(1)-1,hi(1)+1
             idx = i-lo(1)+1
             dp(idx,j,k) = p(i+1,j,k) - p(i-1,j,k)
             denom = max(small_pres,abs(p(i+2,j,k)-p(i-2,j,k)))
             zeta = abs(dp(idx,j,k))/denom
             z(idx,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (u(i-1,j,k)-u(i+1,j,k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i+1,j,k),p(i-1,j,k))
             if ((abs(dp(idx,j,k))/tmp).gt.shktst) then
                chi(idx,j,k) = tst
             else
                chi(idx,j,k) = ZERO
             endif
          enddo
          do i = lo(1),hi(1)
             idx = i-lo(1)+1
             if(dp(idx,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             flatn(i,j,k) = ONE - &
                  max(chi(idx-ishft,j,k)*z(idx-ishft,j,k),chi(idx,j,k)*z(idx,j,k))
          enddo
       enddo
    enddo

    deallocate(dp,z,chi)

    ! y-direction flattening coef
    allocate(dp (lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
    allocate(z  (lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
    allocate(chi(lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
    do k = lo(3),hi(3)
       do i = lo(1),hi(1)
          do j = lo(2)-1,hi(2)+1
             idx = j-lo(2)+1
             dp(i,idx,k) = p(i,j+1,k) - p(i,j-1,k)
             denom = max(small_pres,abs(p(i,j+2,k)-p(i,j-2,k)))
             zeta = abs(dp(i,idx,k))/denom
             z(i,idx,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (v(i,j-1,k)-v(i,j+1,k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j+1,k),p(i,j-1,k))
             if ((abs(dp(i,idx,k))/tmp).gt.shktst) then
                chi(i,idx,k) = tst
             else
                chi(i,idx,k) = ZERO
             endif
          enddo
          do j = lo(2),hi(2)
             idx = j-lo(2)+1
             if(dp(i,idx,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,idx-ishft,k)*z(i,idx-ishft,k),chi(i,idx,k)*z(i,idx,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

    deallocate(dp,z,chi)

    ! z-direction flattening coef
    allocate(dp (lo(1):hi(1),lo(2):hi(2),0:nmax-1))
    allocate(z  (lo(1):hi(1),lo(2):hi(2),0:nmax-1))
    allocate(chi(lo(1):hi(1),lo(2):hi(2),0:nmax-1))
    do j = lo(2),hi(2) 
       do i = lo(1),hi(1)
          do k = lo(3)-1,hi(3)+1
             idx = k-lo(3)+1
             dp(i,j,idx) = p(i,j,k+1) - p(i,j,k-1)
             denom = max(small_pres,abs(p(i,j,k+2)-p(i,j,k-2)))
             zeta = abs(dp(i,j,idx))/denom
             z(i,j,idx) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (w(i,j,k-1)-w(i,j,k+1) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j,k+1),p(i,j,k-1))
             if ((abs(dp(i,j,idx))/tmp).gt.shktst) then
                chi(i,j,idx) = tst
             else
                chi(i,j,idx) = ZERO
             endif
          enddo
          do k = lo(3),hi(3)
             idx = k-lo(3)+1
             if(dp(i,j,idx).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j,idx-ishft)*z(i,j,idx-ishft),chi(i,j,idx)*z(i,j,idx))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo
    
    deallocate(dp,z,chi)

  end subroutine uflaten
  
end module ca_flatten_module
