

c ::: -----------------------------------------------------------
c ::: Add fine grid flux to flux register.  Flux array is a fine grid
c ::: edge based object, Register is a coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3  => index limits for flx
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------
      subroutine frfineadd(reg,reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, r
     &eg_h3,flx,flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3,numcomp
     &,dir,ratio,mult)
      implicit none
      integer    reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3
      integer    flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3
      integer    ratio(3), dir, numcomp
      DOUBLE PRECISION     mult
      DOUBLE PRECISION reg(reg_l1:reg_h1, reg_l2:reg_h2, reg_l3:reg_h3,n
     &umcomp)
      DOUBLE PRECISION flx(flx_l1:flx_h1, flx_l2:flx_h2, flx_l3:flx_h3,n
     &umcomp)

      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

      if (dir .eq. 0) then
c
c        ::::: flux normal to X direction
c
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = reg_l3, reg_h3
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = reg_l2, reg_h2
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j
     &,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else if (dir .eq. 1) then
c        ::::: flux normal to Y direction
         jc = reg_l2
         j = jc*ratioy
         if (reg_l2 .ne. reg_h2) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (j .lt. flx_l2 .or. j .gt. flx_h2) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = reg_l3, reg_h3
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = reg_l1, reg_h1
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j
     &,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else
c
c        ::::: flux normal to Z direction
c
         kc = reg_l3
         k = kc*ratioz
         if (reg_l3 .ne. reg_h3) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (k .lt. flx_l3 .or. k .gt. flx_h3) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

         do joff = 0, ratioy-1
            do n = 1, numcomp
               do jc = reg_l2, reg_h2
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = reg_l1, reg_h1
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j
     &,k,n)
                     end do
                  end do
               end do
            end do
         end do

      end if

      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux times area to flux register.  
c ::: Flux array is a fine grid edge based object, Register is a 
c ::: coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3  => index limits for flx
c ::: area       => edge centered area array
c ::: area_l1, area_l2, area_l3, area_h1, area_h2, area_h3 => index limits for area
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------
      subroutine frfaadd(reg,reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg
     &_h3,flx,flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3,area,area
     &_l1, area_l2, area_l3, area_h1, area_h2, area_h3,numcomp,dir,rat
     &io,mult)
      implicit none
      integer    reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3
      integer    flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3
      integer    area_l1, area_l2, area_l3, area_h1, area_h2, area_h3
      integer    ratio(3), dir, numcomp
      DOUBLE PRECISION     mult
      DOUBLE PRECISION reg(reg_l1:reg_h1, reg_l2:reg_h2, reg_l3:reg_h3,n
     &umcomp)
      DOUBLE PRECISION flx(flx_l1:flx_h1, flx_l2:flx_h2, flx_l3:flx_h3,n
     &umcomp)
      DOUBLE PRECISION area(area_l1:area_h1, area_l2:area_h2, area_l3:ar
     &ea_h3)

      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

      if (dir .eq. 0) then
c
c        ::::: flux normal to X direction
c
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = reg_l3, reg_h3
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = reg_l2, reg_h2
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*area(i,
     &j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else if (dir .eq. 1) then
c
c        ::::: flux normal to Y direction
c
         jc = reg_l2
         j = jc*ratioy
         if (reg_l2 .ne. reg_h2) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (j .lt. flx_l2 .or. j .gt. flx_h2) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = reg_l3, reg_h3
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = reg_l1, reg_h1
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*area(i,
     &j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else
c
c        ::::: flux normal to Z direction
c
         kc = reg_l3
         k = kc*ratioz
         if (reg_l3 .ne. reg_h3) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (k .lt. flx_l3 .or. k .gt. flx_h3) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

         do joff = 0, ratioy-1
            do n = 1, numcomp
               do jc = reg_l2, reg_h2
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = reg_l1, reg_h1
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*area(i,
     &j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      end if

      end

      subroutine frreflux (lo, hi, s, slo, shi, f, flo, fhi,v, vlo, vhi,
     & nc, mult, dir, isloface)
      implicit none
      integer, intent(in) :: lo(3), hi(3), slo(3), shi(3)
      integer, intent(in) :: flo(3), fhi(3), vlo(3), vhi(3)
      integer, intent(in) :: nc, dir, isloface
      DOUBLE PRECISION , intent(in) :: mult
      DOUBLE PRECISION , intent(inout) :: s(slo(1):shi(1),slo(2):shi(2),
     &slo(3):shi(3),nc)
      DOUBLE PRECISION , intent(in ) :: f(flo(1):fhi(1),flo(2):fhi(2),fl
     &o(3):fhi(3),nc)
      DOUBLE PRECISION , intent(in ) :: v(vlo(1):vhi(1),vlo(2):vhi(2),vl
     &o(3):vhi(3))
      !
      integer :: i, j, k, n
      if (isloface .eq. 1) then
         if (dir .eq. 0) then
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i+1,j,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         else if (dir .eq. 1) then
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i,j+1,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         else
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i,j,k+1,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         end if
      else
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)+mult*f(i,j,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
      end if
      end subroutine frreflux
