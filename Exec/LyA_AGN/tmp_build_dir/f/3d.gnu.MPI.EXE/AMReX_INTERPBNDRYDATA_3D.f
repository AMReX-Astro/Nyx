

c ---------------------------------------------------------------
c ::  bdintrpxlo : Interpolation on Xlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpxlo (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l2:cb_h2, cb_l3:cb_h3,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi, ratioy, ratioz

      ratioy = ratios(2)
      ratioz = ratios(3)

      kclo = cb_l3
      kchi = cb_h3
      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_l1-1
      i    = lo(1)-1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do koff = 0, ratioz - 1
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do joff = 0, ratioy - 1
                     do jc = jclo, jchi
                        j = ratioy*jc + joff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do kc = kclo, kchi
               k = ratioz*kc
               do jc = jclo, jchi
                  j = ratioy*jc

                  if ( mask(i,j-1,k) .eq. not_covered .and.mask(i,j+rati
     &oy,k) .eq. not_covered) then
                     derives(jc,kc,1) = 0.5D0*(crse(ic,jc+1,kc,n) - crse
     &(ic,jc-1,kc,n))
                     derives(jc,kc,3) = 0.5D0*(crse(ic,jc+1,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc-1,kc,n))
                  else if (mask(i,j-1,k) .eq. not_covered) then
                     derives(jc,kc,1) = crse(ic,jc,kc,n) - crse(ic,jc-1,
     &kc,n)
                     derives(jc,kc,3) = 0.0D0
                  else if (mask(i,j+ratioy,k) .eq. not_covered) then
                     derives(jc,kc,1) = crse(ic,jc+1,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(jc,kc,3) = 0.0D0
                  else
                     derives(jc,kc,1)  = 0.0D0
                     derives(jc,kc,3) = 0.0D0
                  end if

                  if ( mask(i,j,k-1) .eq. not_covered .and.mask(i,j,k+ra
     &tioz) .eq. not_covered) then
                     derives(jc,kc,2) = 0.5D0*(crse(ic,jc,kc+1,n) - crse
     &(ic,jc,kc-1,n))
                     derives(jc,kc,4) = 0.5D0*(crse(ic,jc,kc+1,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc,kc-1,n))
                  else if (mask(i,j,k-1) .eq. not_covered) then
                     derives(jc,kc,2) = crse(ic,jc,kc,n) - crse(ic,jc,kc
     &-1,n)
                     derives(jc,kc,4) = 0.0D0
                  else if (mask(i,j,k+ratioz) .eq. not_covered) then
                     derives(jc,kc,2) = crse(ic,jc,kc+1,n) - crse(ic,jc,
     &kc,n)
                     derives(jc,kc,4) = 0.0D0
                  else
                     derives(jc,kc,2)  = 0.0D0
                     derives(jc,kc,4)  = 0.0D0
                  end if

                  if ( ( mask(i,j+ratioy,k+ratioz) .ne. not_covered ) .o
     &r.( mask(i,j-1,k+ratioz) .ne. not_covered ) .or.( ma
     &sk(i,j+ratioy,k-1) .ne. not_covered ) .or.( mask(i,j
     &-1,k-1) .ne. not_covered ) ) then

                     derives(jc,kc,5) = 0.0D0
                  else
                     derives(jc,kc,5) = 0.25D0*(crse(ic,jc+1,kc+1,n) - c
     &rse(ic,jc-1,kc+1,n)+ crse(ic,jc-1,kc-1,n) - crse(
     &ic,jc+1,kc-1,n))
                  end if
               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               yy = (dble(koff - ratioz/2) + 0.5D0)/ratioz
               yysq = yy**2
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do joff = 0, ratioy - 1
                     xx = (dble(joff - ratioy/2) + 0.5D0)/ratioy
                     xxsq = xx**2
                     do jc = jclo, jchi
                        j = ratioy*jc + joff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(jc
     &,kc,1) + derives(jc,kc,3)*xxsq + yy*derives(jc
     &,kc,2) + derives(jc,kc,4)*yysq + xx*yy*derives
     &(jc,kc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

c ---------------------------------------------------------------
c ::  bdintrpxhi : Interpolation on Xhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpxhi (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l2:cb_h2, cb_l3:cb_h3,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi, ratioy, ratioz

      ratioy = ratios(2)
      ratioz = ratios(3)

      kclo = cb_l3
      kchi = cb_h3
      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_h1+1
      i    = hi(1)+1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do joff = 0, ratioy - 1
                     do jc = jclo, jchi
                        j = ratioy*jc + joff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do kc = kclo, kchi
               k = ratioz*kc
               do jc = jclo, jchi
                  j = ratioy*jc

                  if (mask(i,j-1,k) .eq. not_covered .and.mask(i,j+ratio
     &y,k) .eq. not_covered) then
                     derives(jc,kc,1) = 0.5D0*(crse(ic,jc+1,kc,n) - crse
     &(ic,jc-1,kc,n))
                     derives(jc,kc,3) = 0.5D0*(crse(ic,jc+1,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc-1,kc,n))
                  else if (mask(i,j-1,k) .eq. not_covered) then
                     derives(jc,kc,1) = crse(ic,jc,kc,n) - crse(ic,jc-1,
     &kc,n)
                     derives(jc,kc,3) = 0.0D0
                  else if (mask(i,j+ratioy,k) .eq. not_covered) then
                     derives(jc,kc,1) = crse(ic,jc+1,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(jc,kc,3) = 0.0D0
                  else
                     derives(jc,kc,1)  = 0.0D0
                     derives(jc,kc,3) = 0.0D0
                  end if

                  if (mask(i,j,k-1) .eq. not_covered .and.mask(i,j,k+rat
     &ioz) .eq. not_covered) then
                     derives(jc,kc,2) = 0.5D0*(crse(ic,jc,kc+1,n) - crse
     &(ic,jc,kc-1,n))
                     derives(jc,kc,4) = 0.5D0*(crse(ic,jc,kc+1,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc,kc-1,n))
                  else if (mask(i,j,k-1) .eq. not_covered) then
                     derives(jc,kc,2) = crse(ic,jc,kc,n) - crse(ic,jc,kc
     &-1,n)
                     derives(jc,kc,4) = 0.0D0
                  else if (mask(i,j,k+ratioz) .eq. not_covered) then
                     derives(jc,kc,2) = crse(ic,jc,kc+1,n) - crse(ic,jc,
     &kc,n)
                     derives(jc,kc,4) = 0.0D0
                  else
                     derives(jc,kc,2) = 0.0D0
                     derives(jc,kc,4) = 0.0D0
                  end if

                  if (( mask(i,j+ratioy,k+ratioz) .ne. not_covered ) .or
     &.( mask(i,j-1,k+ratioz) .ne. not_covered ) .or.( mas
     &k(i,j+ratioy,k-1) .ne. not_covered ) .or.( mask(i,j-
     &1,k-1) .ne. not_covered ) ) then

                     derives(jc,kc,5) = 0.0D0
                  else
                     derives(jc,kc,5) = 0.25D0*(crse(ic,jc+1,kc+1,n) - c
     &rse(ic,jc-1,kc+1,n)+ crse(ic,jc-1,kc-1,n) - crse(
     &ic,jc+1,kc-1,n))
                  end if

               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               yy = (dble(koff - ratioz/2) + 0.5D0)/ratioz
               yysq = yy**2
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do joff = 0, ratioy - 1
                     xx = (dble(joff - ratioy/2) + 0.5D0)/ratioy
                     xxsq = xx**2
                     do jc = jclo, jchi
                        j = ratioy*jc + joff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(jc
     &,kc,1) + derives(jc,kc,3)*xxsq + yy*derives(jc
     &,kc,2) + derives(jc,kc,4)*yysq + xx*yy*derives
     &(jc,kc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

c ---------------------------------------------------------------
c ::  bdintrpylo : Interpolation on Ylo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpylo (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1, cb_l3:cb_h3,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi, ratiox, ratioz

      ratiox = ratios(1)
      ratioz = ratios(3)

      kclo = cb_l3
      kchi = cb_h3
      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_l2-1
      j    = lo(2)-1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox - 1
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do kc = kclo, kchi
               k = ratioz*kc
               do ic = iclo, ichi
                  i = ratiox*ic

                  if (mask(i-1,j,k) .eq. not_covered .and.mask(i+ratiox,
     &j,k) .eq. not_covered) then
                     derives(ic,kc,1) = 0.5D0*(crse(ic+1,jc,kc,n) - crse
     &(ic-1,jc,kc,n))
                     derives(ic,kc,3) = 0.5D0*(crse(ic+1,jc,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic-1,jc,kc,n))
                  else if (mask(i-1,j,k) .eq. not_covered) then
                     derives(ic,kc,1) = crse(ic,jc,kc,n) - crse(ic-1,jc,
     &kc,n)
                     derives(ic,kc,3) = 0.0D0
                  else if (mask(i+ratiox,j,k) .eq. not_covered) then
                     derives(ic,kc,1) = crse(ic+1,jc,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,kc,3) = 0.0D0
                  else
                     derives(ic,kc,1)  = 0.0D0
                     derives(ic,kc,3)  = 0.0D0
                  end if

                  if (mask(i,j,k-1) .eq. not_covered .and.mask(i,j,k+rat
     &ioz) .eq. not_covered) then
                     derives(ic,kc,2) = 0.5D0*(crse(ic,jc,kc+1,n) - crse
     &(ic,jc,kc-1,n))
                     derives(ic,kc,4) = 0.5D0*(crse(ic,jc,kc+1,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc,kc-1,n))
                  else if (mask(i,j,k-1) .eq. not_covered) then
                     derives(ic,kc,2) = crse(ic,jc,kc,n) - crse(ic,jc,kc
     &-1,n)
                     derives(ic,kc,4) = 0.0D0
                  else if (mask(i,j,k+ratioz) .eq. not_covered) then
                     derives(ic,kc,2) = crse(ic,jc,kc+1,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,kc,4) = 0.0D0
                  else
                     derives(ic,kc,2) = 0.0D0
                     derives(ic,kc,4) = 0.0D0
                  end if

                  if (( mask(i+ratiox,j,k+ratioz) .ne. not_covered ) .or
     &.( mask(i-1,j,k+ratioz) .ne. not_covered ) .or.( mas
     &k(i+ratiox,j,k-1) .ne. not_covered ) .or.( mask(i-1,
     &j,k-1) .ne. not_covered ) ) then

                     derives(ic,kc,5) = 0.0D0
                  else
                     derives(ic,kc,5) = 0.25D0*(crse(ic+1,jc,kc+1,n) - c
     &rse(ic-1,jc,kc+1,n)+ crse(ic-1,jc,kc-1,n) - crse(
     &ic+1,jc,kc-1,n))
                  end if

               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               yy = (dble(koff - ratioz/2) + 0.5D0)/ratioz
               yysq = yy**2
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox - 1
                     xx = (dble(ioff - ratiox/2) + 0.5D0)/ratiox
                     xxsq = xx**2
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic
     &,kc,1) + derives(ic,kc,3)*xxsq + yy*derives(ic
     &,kc,2) + derives(ic,kc,4)*yysq + xx*yy*derives
     &(ic,kc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

c ---------------------------------------------------------------
c ::  bdintrpyhi : Interpolation on Yhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpyhi (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1, cb_l3:cb_h3,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi, ratiox, ratioz

      ratiox = ratios(1)
      ratioz = ratios(3)

      kclo = cb_l3
      kchi = cb_h3
      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_h2+1
      j    = hi(2)+1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox - 1
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do kc = kclo, kchi
               k = ratioz*kc
               do ic = iclo, ichi
                  i = ratiox*ic

                  if (mask(i-1,j,k) .eq. not_covered .and.mask(i+ratiox,
     &j,k) .eq. not_covered) then
                     derives(ic,kc,1) = 0.5D0*(crse(ic+1,jc,kc,n) - crse
     &(ic-1,jc,kc,n))
                     derives(ic,kc,3) = 0.5D0*(crse(ic+1,jc,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic-1,jc,kc,n))
                  else if (mask(i-1,j,k) .eq. not_covered) then
                     derives(ic,kc,1) = crse(ic,jc,kc,n) - crse(ic-1,jc,
     &kc,n)
                     derives(ic,kc,3) = 0.0D0
                  else if (mask(i+ratiox,j,k) .eq. not_covered) then
                     derives(ic,kc,1) = crse(ic+1,jc,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,kc,3) = 0.0D0
                  else
                     derives(ic,kc,1) = 0.0D0
                     derives(ic,kc,3) = 0.0D0
                  end if

                  if (mask(i,j,k-1) .eq. not_covered .and.mask(i,j,k+rat
     &ioz) .eq. not_covered) then
                     derives(ic,kc,2) = 0.5D0*(crse(ic,jc,kc+1,n) - crse
     &(ic,jc,kc-1,n))
                     derives(ic,kc,4) = 0.5D0*(crse(ic,jc,kc+1,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc,kc-1,n))
                  else if (mask(i,j,k-1) .eq. not_covered) then
                     derives(ic,kc,2) = crse(ic,jc,kc,n) - crse(ic,jc,kc
     &-1,n)
                     derives(ic,kc,4) = 0.0D0
                  else if (mask(i,j,k+ratioz) .eq. not_covered) then
                     derives(ic,kc,2) = crse(ic,jc,kc+1,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,kc,4) = 0.0D0
                  else
                     derives(ic,kc,2)  = 0.0D0
                     derives(ic,kc,4)  = 0.0D0
                  end if

                  if ( ( mask(i+ratiox,j,k+ratioz) .ne. not_covered ) .o
     &r.( mask(i-1,j,k+ratioz) .ne. not_covered ) .or.( ma
     &sk(i+ratiox,j,k-1) .ne. not_covered ) .or.( mask(i-1
     &,j,k-1) .ne. not_covered ) ) then

                     derives(ic,kc,5) = 0.0D0
                  else
                     derives(ic,kc,5) = 0.25D0*(crse(ic+1,jc,kc+1,n) - c
     &rse(ic-1,jc,kc+1,n)+ crse(ic-1,jc,kc-1,n) - crse(
     &ic+1,jc,kc-1,n))
                  end if
               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do koff = 0, ratioz - 1
               yy = (dble(koff - ratioz/2) + 0.5D0)/ratioz
               yysq = yy**2
               do kc = kclo,kchi
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox - 1
                     xx = (dble(ioff - ratiox/2) + 0.5D0)/ratiox
                     xxsq = xx**2
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic
     &,kc,1) + derives(ic,kc,3)*xxsq + yy*derives(ic
     &,kc,2) + derives(ic,kc,4)*yysq + xx*yy*derives
     &(ic,kc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

c ---------------------------------------------------------------
c ::  bdintrpzlo : Interpolation on Zlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpzlo (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1, cb_l2:cb_h2,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi, ratiox, ratioy

      ratiox = ratios(1)
      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      iclo = cb_l1
      ichi = cb_h1
      kc   = cb_l3-1
      k    = lo(3)-1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: interpolate to fine grid
            !
            do joff = 0, ratioy - 1
               do jc = jclo,jchi
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox - 1
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do jc = jclo, jchi
               j = ratioy*jc
               do ic = iclo, ichi
                  i = ratiox*ic

                  if (mask(i-1,j,k) .eq. not_covered .and. mask(i+ratiox
     &,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = 0.5D0*(crse(ic+1,jc,kc,n) - crse
     &(ic-1,jc,kc,n) )
                     derives(ic,jc,3) = 0.5D0*(crse(ic+1,jc,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n) + crse(ic-1,jc,kc,n) )
                  else if (mask(i-1,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = crse(ic,jc,kc,n) - crse(ic-1,jc,
     &kc,n)
                     derives(ic,jc,3) = 0.0D0
                  else if (mask(i+ratiox,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = crse(ic+1,jc,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,jc,3) = 0.0D0                     
                  else
                     derives(ic,jc,1)  = 0.0D0
                     derives(ic,jc,3)  = 0.0D0
                  end if

                  if (mask(i,j-1,k) .eq. not_covered .and. mask(i,j+rati
     &oy,k) .eq. not_covered) then
                     derives(ic,jc,2) = 0.5D0*(crse(ic,jc+1,kc,n) - crse
     &(ic,jc-1,kc,n) )
                     derives(ic,jc,4) = 0.5D0*(crse(ic,jc+1,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n) + crse(ic,jc-1,kc,n) )
                  else if (mask(i,j-1,k) .eq. not_covered) then
                     derives(ic,jc,2) = crse(ic,jc,kc,n) - crse(ic,jc-1,
     &kc,n)
                     derives(ic,jc,4) = 0.0D0
                  else if (mask(i,j+ratioy,k) .eq. not_covered) then
                     derives(ic,jc,2) = crse(ic,jc+1,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,jc,4) = 0.0D0
                  else
                     derives(ic,jc,2)  = 0.0D0
                     derives(ic,jc,4)  = 0.0D0
                  end if

                  if (( mask(i+ratiox,j+ratioy,k) .ne. not_covered ) .or
     &.( mask(i-1,j+ratioy,k) .ne. not_covered ) .or.( mas
     &k(i+ratiox,j-1,k) .ne. not_covered ) .or.( mask(i-1,
     &j-1,k) .ne. not_covered ) ) then

                     derives(ic,jc,5) = 0.0D0
                  else
                     derives(ic,jc,5) = 0.25D0*(crse(ic+1,jc+1,kc,n) - c
     &rse(ic-1,jc+1,kc,n)+ crse(ic-1,jc-1,kc,n) - crse(
     &ic+1,jc-1,kc,n))
                  end if
               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do joff = 0, ratioy - 1
               yy = (dble(joff - ratioy/2) + 0.5D0)/ratioy
               yysq = yy**2
               do jc = jclo,jchi
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox - 1
                     xx = (dble(ioff - ratiox/2) + 0.5D0)/ratiox
                     xxsq = xx**2
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic
     &,jc,1) + derives(ic,jc,3)*xxsq + yy*derives(ic
     &,jc,2) + derives(ic,jc,4)*yysq + xx*yy*derives
     &(ic,jc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

c ---------------------------------------------------------------
c ::  bdintrpzhi : Interpolation on Zhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------
      subroutine bdintrpzhi (bdry,bdry_l1, bdry_l2, bdry_l3, bdry_h1, bd
     &ry_h2, bdry_h3,lo,hi,cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,nv
     &ar,ratios,not_covered,mask,mask_l1, mask_l2, mask_l3, mask_h1, m
     &ask_h2, mask_h3,crse,crse_l1, crse_l2, crse_l3, crse_h1, crse_h2
     &, crse_h3,derives,max_order)
      implicit none
      integer  nvar, ratios(3), not_covered,max_order
      integer  lo(3), hi(3)
      integer  bdry_l1, bdry_l2, bdry_l3, bdry_h1, bdry_h2, bdry_h3
      integer  cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer  mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3
      integer  crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      DOUBLE PRECISION bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2, bdry_l3:bd
     &ry_h3,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1, cb_l2:cb_h2,5)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2, mask_l3:mask_h3)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)

      DOUBLE PRECISION   xx, yy, xxsq, yysq
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi, ratiox, ratioy

      ratiox = ratios(1)
      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      iclo = cb_l1
      ichi = cb_h1
      kc   = cb_h3+1
      k    = hi(3)+1

      if (max_order.eq.1) then

         do n = 1, nvar
            !
            ! ::::: interpolate to fine grid
            !
            do joff = 0, ratioy - 1
               do jc = jclo,jchi
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox - 1
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         do n = 1, nvar
            !
            ! ::::: define interp coefs
            !
            do jc = jclo, jchi
               j = ratioy*jc
               do ic = iclo, ichi
                  i = ratiox*ic

                  if (mask(i-1,j,k) .eq. not_covered .and. mask(i+ratiox
     &,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = 0.5D0*(crse(ic+1,jc,kc,n) - crse
     &(ic-1,jc,kc,n))
                     derives(ic,jc,3) = 0.5D0*(crse(ic+1,jc,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic-1,jc,kc,n))
                  else if (mask(i-1,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = crse(ic,jc,kc,n) - crse(ic-1,jc,
     &kc,n)
                     derives(ic,jc,3) = 0.0D0
                  else if (mask(i+ratiox,j,k) .eq. not_covered) then
                     derives(ic,jc,1) = crse(ic+1,jc,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,jc,3) = 0.0D0
                  else
                     derives(ic,jc,1) = 0.0D0
                     derives(ic,jc,3) = 0.0D0
                  end if

                  if (mask(i,j-1,k) .eq. not_covered .and. mask(i,j+rati
     &oy,k) .eq. not_covered) then
                     derives(ic,jc,2) = 0.5D0*(crse(ic,jc+1,kc,n) - crse
     &(ic,jc-1,kc,n))
                     derives(ic,jc,4) = 0.5D0*(crse(ic,jc+1,kc,n) - 2.0D
     &0*crse(ic,jc,kc,n)+ crse(ic,jc-1,kc,n))
                  else if (mask(i,j-1,k) .eq. not_covered) then
                     derives(ic,jc,2) = crse(ic,jc,kc,n) - crse(ic,jc-1,
     &kc,n)
                     derives(ic,jc,4) = 0.0D0
                  else if (mask(i,j+ratioy,k) .eq. not_covered) then
                     derives(ic,jc,2) = crse(ic,jc+1,kc,n) - crse(ic,jc,
     &kc,n)
                     derives(ic,jc,4) = 0.0D0
                  else 
                     derives(ic,jc,2)  = 0.0D0
                     derives(ic,jc,4)  = 0.0D0
                  end if

                  if (( mask(i+ratiox,j+ratioy,k) .ne. not_covered ) .or
     &.( mask(i-1,j+ratioy,k) .ne. not_covered ) .or.( mas
     &k(i+ratiox,j-1,k) .ne. not_covered ) .or.( mask(i-1,
     &j-1,k) .ne. not_covered ) ) then

                     derives(ic,jc,5) = 0.0D0
                  else
                     derives(ic,jc,5) = 0.25D0*(crse(ic+1,jc+1,kc,n) - c
     &rse(ic-1,jc+1,kc,n)+ crse(ic-1,jc-1,kc,n) - crse(
     &ic+1,jc-1,kc,n))
                  end if
               end do
            end do
            !
            ! ::::: interpolate to fine grid
            !
            do joff = 0, ratioy - 1
               yy = (dble(joff - ratioy/2) + 0.5D0)/ratioy
               yysq = yy**2
               do jc = jclo,jchi
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox - 1
                     xx = (dble(ioff - ratiox/2) + 0.5D0)/ratiox
                     xxsq = xx**2
                     do ic = iclo, ichi
                        i = ratiox*ic + ioff
                        bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic
     &,jc,1) + derives(ic,jc,3)*xxsq + yy*derives(ic
     &,jc,2) + derives(ic,jc,4)*yysq + xx*yy*derives
     &(ic,jc,5) 
                     end do
                  end do
               end do
            end do
         end do

      endif

      end

