module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                 u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                 Ip,Im, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old)

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type

    implicit none

    integer         , intent(in   ) ::   s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in   ) ::  qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer         , intent(in   ) ::   f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer         , intent(in   ) ::  ilo1,ilo2,ihi1,ihi2
 
    real(rt), intent(in   ) ::      s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in   ) ::      u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in   ) ::   cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in   ) ::  flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(inout) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(inout) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    real(rt), intent(in   ) ::  dx,dy,dz,dt,a_old

    real(rt) :: dt_over_a
    integer          :: k3d,kc

    dt_over_a = dt / a_old
   
    call ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
    
  subroutine ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt_over_a,k3d,kc)

    use amrex_error_module
    use amrex_mempool_module, only: amrex_allocate, amrex_deallocate
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use amrex_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    real(rt), intent(in) :: s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in) :: u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in) :: cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    real(rt), intent(in) :: flatn(f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    real(rt), intent(out) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(out) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    ! Note that dt_over_a = dt / a_old
    real(rt), intent(in) :: dx,dy,dz,dt_over_a
    integer, intent(in)          :: k3d,kc

    ! local
    integer i,j,k

    real(rt) :: dxinv,dyinv,dzinv

    real(rt), pointer :: dsl(:), dsr(:), dsc(:)
    real(rt), pointer :: sigma(:), s6(:)

    ! s_{\ib,+}, s_{\ib,-}
    real(rt), pointer :: sp(:)
    real(rt), pointer :: sm(:)

    ! \delta s_{\ib}^{vL}
    real(rt), pointer :: dsvl(:,:)
    real(rt), pointer :: dsvlm(:,:)
    real(rt), pointer :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    real(rt), pointer :: sedge(:,:)
    real(rt), pointer :: sedgez(:,:,:)

    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy
    dzinv = 1.0d0/dz

    ! cell-centered indexing
    call amrex_allocate(sp,ilo1-1,ihi1+1)
    call amrex_allocate(sm,ilo1-1,ihi1+1)

    call amrex_allocate(sigma,ilo1-1,ihi1+1)
    call amrex_allocate(s6,ilo1-1,ihi1+1)

    if (ppm_type .ne. 1) &
         call amrex_error("Should have ppm_type = 1 in ppm_type1")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    call amrex_allocate(dsvl,ilo1-2,ihi1+2,ilo2-1,ihi2+1)

    ! edge-centered indexing for x-faces -- ppm_type = 1 only
    call amrex_allocate(sedge,ilo1-1,ihi1+2,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-2,ihi1+2)
    call amrex_allocate(dsl,ilo1-2,ihi1+2)
    call amrex_allocate(dsr,ilo1-2,ihi1+2)

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    dsvl = ZERO
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc(i) = HALF * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl(i) = TWO  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr(i) = TWO  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate s to x-edges
       do i=ilo1-1,ihi1+2
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do

       ! copy sedge into sp and sm
       sp(ilo1-1:ihi1+1) = sedge(ilo1:ihi1+2,j)
       sm(ilo1-1:ihi1+1) = sedge(ilo1-1:ihi1+1,j)

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters -- note this version of the limiting comes
       ! from Colella and Sekora (2008), not the original PPM paper.
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute x-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! Ip/m is the integral under the parabola for the extent
       ! that a wave can travel over a timestep
       !
       ! Ip integrates to the right edge of a cell
       ! Im integrates to the left edge of a cell

       ! u-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,1) = sp(i)
          else
             Ip(i,j,kc,1,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,1) = sm(i)
          else
             Im(i,j,kc,1,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! u wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) <= ZERO) then
             Ip(i,j,kc,1,2) = sp(i)
          else
             Ip(i,j,kc,1,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1) >= ZERO) then
             Im(i,j,kc,1,2) = sm(i)
          else
             Im(i,j,kc,1,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! u+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,1)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dxinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,3) = sp(i)
          else
             Ip(i,j,kc,1,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,3) = sm(i)
          else
             Im(i,j,kc,1,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(sedge)
    call amrex_deallocate(dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-2,ihi2+2)

    ! edge-centered indexing for y-faces
    call amrex_allocate(sedge,ilo1-1,ihi1+1,ilo2-1,ihi2+2)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    dsvl = ZERO
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc(i) = HALF * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl(i) = TWO  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr(i) = TWO  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       do i=ilo1-1,ihi1+1
          sedge(i,j) = HALF*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    do j=ilo2-1,ihi2+1

       ! copy sedge into sp and sm
       sp = sedge(ilo1-1:ihi1+1,j+1)
       sm = sedge(ilo1-1:ihi1+1,j  )

       if (ppm_flatten_before_integrals == 1) then
          ! Flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! Modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! Flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute y-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! v-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,1) = sp(i)
          else
             Ip(i,j,kc,2,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,1) = sm(i)
          else
             Im(i,j,kc,2,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! v wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) <= ZERO) then
             Ip(i,j,kc,2,2) = sp(i)
          else
             Ip(i,j,kc,2,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2) >= ZERO) then
             Im(i,j,kc,2,2) = sm(i)
          else
             Im(i,j,kc,2,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! v+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,2)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dyinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,3) = sp(i)
          else
             Ip(i,j,kc,2,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,3) = sm(i)
          else
             Im(i,j,kc,2,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    call amrex_allocate( dsvl,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlm,ilo1-1,ihi1+1,ilo2-1,ihi2+1)
    call amrex_allocate(dsvlp,ilo1-1,ihi1+1,ilo2-1,ihi2+1)

    ! cell-centered indexing
    call amrex_allocate(dsc,ilo1-1,ihi1+1)
    call amrex_allocate(dsl,ilo1-1,ihi1+1)
    call amrex_allocate(dsr,ilo1-1,ihi1+1)

    call amrex_allocate(sedgez,ilo1-1,ihi1+1,ilo2-2,ihi2+3,k3d-1,k3d+2)

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    dsvl  = ZERO
    dsvlm = ZERO
    dsvlp = ZERO

    do j=ilo2-1,ihi2+1

       ! compute on slab below
       k = k3d-1
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvlm(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on slab above
       k = k3d+1
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvlp(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! compute on current slab
       k = k3d
       dsc = HALF * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k-1))
       dsl = TWO  * (s(ilo1-1:ihi1+1,j,k  ) - s(ilo1-1:ihi1+1,j,k-1))
       dsr = TWO  * (s(ilo1-1:ihi1+1,j,k+1) - s(ilo1-1:ihi1+1,j,k  ))
       do i = ilo1-1, ihi1+1
          if (dsl(i)*dsr(i) .gt. ZERO) &
               dsvl(i,j) = sign(ONE,dsc(i))*min(abs(dsc(i)),abs(dsl(i)),abs(dsr(i)))
       end do

       ! interpolate to lo face
       k = k3d
       sm = HALF*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - SIXTH*(dsvl(ilo1-1:ihi1+1,j)-dsvlm(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sm = max(sm,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sm = min(sm,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       ! interpolate to hi face
       k = k3d+1
       sp = HALF*(s(ilo1-1:ihi1+1,j,k)+s(ilo1-1:ihi1+1,j,k-1)) - SIXTH*(dsvlp(ilo1-1:ihi1+1,j)-dsvl(ilo1-1:ihi1+1,j))
       ! make sure sedge lies in between adjacent cell-centered values
       sp = max(sp,min(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))
       sp = min(sp,max(s(ilo1-1:ihi1+1,j,k),s(ilo1-1:ihi1+1,j,k-1)))

       if (ppm_flatten_before_integrals == 1) then
          ! flatten the parabola BEFORE doing the other
          ! monotonization -- this is the method that Flash does
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! modify using quadratic limiters
       do i = ilo1-1, ihi1+1
          if ((sp(i)-s(i,j,k3d))*(s(i,j,k3d)-sm(i)) .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(sp(i)-s(i,j,k3d)) .ge. TWO*abs(sm(i)-s(i,j,k3d))) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(sm(i)-s(i,j,k3d)) .ge. TWO*abs(sp(i)-s(i,j,k3d))) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       if (ppm_flatten_before_integrals == 2) then
          ! flatten the parabola AFTER doing the monotonization --
          ! this is the method that Miller & Colella do
          sm = flatn(ilo1-1:ihi1+1,j,k3d)*sm + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
          sp = flatn(ilo1-1:ihi1+1,j,k3d)*sp + (ONE-flatn(ilo1-1:ihi1+1,j,k3d))*s(ilo1-1:ihi1+1,j,k3d)
       endif

       ! compute z-component of Ip and Im
       s6 = SIX*s(ilo1-1:ihi1+1,j,k3d) - THREE*(sm+sp)

       ! w-c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)-cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,1) = sp(i)
          else
             Ip(i,j,kc,3,1) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,1) = sm(i)
          else
             Im(i,j,kc,3,1) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! w wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) <= ZERO) then
             Ip(i,j,kc,3,2) = sp(i)
          else
             Ip(i,j,kc,3,2) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3) >= ZERO) then
             Im(i,j,kc,3,2) = sm(i)
          else
             Im(i,j,kc,3,2) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       ! w+c wave
       sigma = abs(u(ilo1-1:ihi1+1,j,k3d,3)+cspd(ilo1-1:ihi1+1,j,k3d))*dt_over_a*dzinv

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,3) = sp(i)
          else
             Ip(i,j,kc,3,3) = sp(i) - &
               HALF*sigma(i)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

       do i = ilo1-1, ihi1+1
          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,3) = sm(i)
          else
             Im(i,j,kc,3,3) = sm(i) + &
               HALF*sigma(i)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma(i))*s6(i))
          endif
       end do

    end do

    call amrex_deallocate(dsc)
    call amrex_deallocate(dsl)
    call amrex_deallocate(dsr)
    call amrex_deallocate(dsvl)
    call amrex_deallocate(dsvlm)
    call amrex_deallocate(dsvlp)
    call amrex_deallocate(sp)
    call amrex_deallocate(sm)
    call amrex_deallocate(sedgez)
    call amrex_deallocate(sigma)
    call amrex_deallocate(s6)

  end subroutine ppm_type1

end module ppm_module

