

c ::: --------------------------------------------------------------
c ::: nbinterp:  node based bilinear interpolation
c :::
c ::: INPUTS/OUTPUTS
c ::: fine        <=>  (modify) fine grid array
c ::: fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c :::
c ::: crse         =>  (const)  coarse grid data widened by 1 zone
c ::: crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3   =>  (const)  index limits of coarse grid
c :::
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: num_slp      =>  (const)  number of types of slopes
c :::
c ::: TEMPORARY ARRAYS
c ::: sl           =>  num_slp 1-D slope arrays
c ::: --------------------------------------------------------------
c :::
      subroutine nbinterp (crse, crse_l1, crse_l2, crse_l3, crse_h1, crs
     &e_h2, crse_h3, cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,fine, fi
     &ne_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3, fb_l1, fb_l2
     &, fb_l3, fb_h1, fb_h2, fb_h3,lratiox, lratioy, lratioz, nvar,sl,
     & num_slp,actual_comp, actual_state)
      implicit none
      integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      integer cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer fb_l1, fb_l2, fb_l3, fb_h1, fb_h2, fb_h3
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer actual_comp,actual_state
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3,nvar)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)
      DOUBLE PRECISION sl(cb_l1:cb_h1,num_slp)

      integer ioff, joff, koff
      integer i, j, k, ic, jc, kc, n
      integer ilo, ihi, jlo, jhi, klo, khi
      integer iratio, jratio, kratio
      integer kstrt, kstop, jstrt, jstop, istrt, istop

      DOUBLE PRECISION fx, fy,fz
      DOUBLE PRECISION sx, sy, sz, sxy, sxz, syz, sxyz
      DOUBLE PRECISION RX, RY, RZ, RXY, RXZ, RYZ, RXYZ
      DOUBLE PRECISION dx00, d0x0, d00x, dx10, dx01, d0x1, dx11

      RX   = 1.0D0/dfloat(lratiox)
      RY   = 1.0D0/dfloat(lratioy)
      RZ   = 1.0D0/dfloat(lratioz)
      RXY  = RX*RY
      RXZ  = RX*RZ
      RYZ  = RY*RZ
      RXYZ = RX*RY*RZ

      do n = 1, nvar

         do kc = cb_l3, cb_h3-1

            kratio = lratioz-1
            if (kc .eq. cb_h3-1) kratio = lratioz
            kstrt = kc*lratioz
            kstop = kstrt + kratio
            klo = max(fine_l3,kstrt) - kstrt
            khi = min(fine_h3,kstop) - kstrt
            do jc = cb_l2, cb_h2-1

               jratio = lratioy-1
               if (jc .eq. cb_h2-1) jratio = lratioy
               jstrt = jc*lratioy
               jstop = jstrt + jratio
               jlo = max(fine_l2,jstrt) - jstrt
               jhi = min(fine_h2,jstop) - jstrt

               do ic = cb_l1, cb_h1-1

                  iratio = lratiox-1
                  if (ic .eq. cb_h1-1) iratio = lratiox
                  istrt = ic*lratiox
                  istop = istrt + iratio
                  ilo = max(fine_l1,istrt) - istrt
                  ihi = min(fine_h1,istop) - istrt
                  !
                  ! ::::: compute slopes
                  !
                  dx00 = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
                  d0x0 = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
                  d00x = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)

                  dx10 = crse(ic+1,jc+1,kc,n) - crse(ic,jc+1,kc,n)
                  dx01 = crse(ic+1,jc,kc+1,n) - crse(ic,jc,kc+1,n)
                  d0x1 = crse(ic,jc+1,kc+1,n) - crse(ic,jc,kc+1,n)

                  dx11 = crse(ic+1,jc+1,kc+1,n) - crse(ic,jc+1,kc+1,n)

                  sx   = RX*dx00
                  sy   = RY*d0x0
                  sz   = RZ*d00x
                  sxy  = RXY*(dx10 - dx00)
                  sxz  = RXZ*(dx01 - dx00)
                  syz  = RYZ*(d0x1 - d0x0)
                  sxyz = RXYZ*(dx11 - dx01 - dx10 + dx00)
                  !
                  ! ::::: interpolate to fine grid
                  !
                  do koff = klo, khi
                     k = lratioz*kc + koff
                     fz = dfloat(koff)
                     do joff = jlo, jhi
                        j = lratioy*jc + joff
                        fy = dfloat(joff)
                        do ioff = ilo, ihi
                           i = lratiox*ic + ioff
                           fx = dfloat(ioff)
                           fine(i,j,k,n) = crse(ic,jc,kc,n) +fx*sx + fy*
     &sy + fz*sz +fx*fy*sxy + fx*fz*sxz + fy*fz*s
     &yz +fx*fy*fz*sxyz
                        end do
                     end do
                  end do           

               end do
            end do
         end do
      end do

      end

c ::: 
c ::: --------------------------------------------------------------
c ::: cbinterp:  cell centered bilinear interpolation
c ::: 
c ::: NOTE: it is assumed that the coarse grid array is
c ::: large enough to define interpolated values
c ::: in the region fblo:fbhi on the fine grid
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3   =>  (const)  index limits of fine grid
c ::: fb_l1, fb_l2, fb_l3, fb_h1, fb_h2, fb_h3     =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3   =>  (const)  index limits of coarse grid
c ::: 
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: slx,sly,slxy =>  1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine cbinterp (crse, crse_l1, crse_l2, crse_l3, crse_h1, crs
     &e_h2, crse_h3, cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3,fine, fi
     &ne_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3, fb_l1, fb_l2
     &, fb_l3, fb_h1, fb_h2, fb_h3,lratiox, lratioy, lratioz, nvar,sl,
     & num_slp, strip, strip_lo, strip_hi,actual_comp, actual_state)

      implicit none

      integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      integer cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer fb_l1, fb_l2, fb_l3, fb_h1, fb_h2, fb_h3
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer strip_lo, strip_hi
      integer actual_comp,actual_state
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3,nvar)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3,nvar)
      DOUBLE PRECISION strip(strip_lo:strip_hi)
      DOUBLE PRECISION sl(cb_l1:cb_h1, num_slp)

      call bl_abort("FORT_CBINTERP not implemented")

      end
c ::: 
c ::: --------------------------------------------------------------
c ::: linccinterp:   linear conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: The interpolation is linear in that it uses a
c ::: a limiting scheme that preserves the value of 
c ::: any linear combination of the
c ::: coarse grid data components--e.g.,
c ::: if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
c ::: sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
c ::: in all fine cells if,jf covering coarse cell ic,jc.
c ::: 
c ::: If lin_limit = 0, the interpolation scheme is identical to
c ::: that used in ccinterp for limslope=1; the results should
c ::: be exactly the same -- difference = hard 0.
c :::
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: clo,chi      =>  (const)  index limits of crse grid
c ::: cslo,cshi    =>  (const)  coarse grid index limits where
c :::				slopes are to be defined. This is
c :::				the projection of (fblo,fbhi) down
c :::				to the coarse level 
c ::: ucslope      =>  (modify) temp array of unlimited coarse grid slopes
c ::: lcslope      =>  (modify) temp array of limited coarse grid slopes
c ::: slope_factor =>  (modify) temp array of slope limiting factors
c ::: lin_limit    =>  (const)  != 0 => do linear slope limiting scheme
c :::
c ::: --------------------------------------------------------------
c ::: 
      subroutine linccinterp (fine, fine_l1, fine_l2, fine_l3, fine_h1, 
     &fine_h2, fine_h3, fblo, fbhi, fvcb_l1, fvcb_l2, fvcb_l3, fvcb_h1
     &, fvcb_h2, fvcb_h3, crse, crse_l1, crse_l2, crse_l3, crse_h1, cr
     &se_h2, crse_h3, cvcb_l1, cvcb_l2, cvcb_l3, cvcb_h1, cvcb_h2, cvc
     &b_h3,uc_xslope, lc_xslope, xslope_factor,uc_yslope, lc_yslope, y
     &slope_factor,uc_zslope, lc_zslope, zslope_factor,cslope_l1, cslo
     &pe_l2, cslope_l3, cslope_h1, cslope_h2, cslope_h3,cslopelo, cslo
     &pehi,nvar, lratiox, lratioy, lratioz, bc, lim_slope, lin_limit,f
     &vcx, fvcy, fvcz, cvcx, cvcy, cvcz,voffx,voffy,voffz, alpha, cmax
     &, cmin,actual_comp, actual_state)
      implicit none

      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      integer fvcb_l1, fvcb_l2, fvcb_l3, fvcb_h1, fvcb_h2, fvcb_h3
      integer cvcb_l1, cvcb_l2, cvcb_l3, cvcb_h1, cvcb_h2, cvcb_h3
      integer cslope_l1, cslope_l2, cslope_l3, cslope_h1, cslope_h2, csl
     &ope_h3
      integer fblo(3), fbhi(3)
      integer cslopelo(3), cslopehi(3)
      integer lratiox, lratioy, lratioz, nvar
      integer lim_slope, lin_limit
      integer bc(3,2,nvar)
      integer actual_comp,actual_state
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3,nvar)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3, nvar)
      DOUBLE PRECISION uc_xslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION lc_xslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION xslope_factor(cslope_l1:cslope_h1, cslope_l2:cslo
     &pe_h2, cslope_l3:cslope_h3)
      DOUBLE PRECISION uc_yslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION lc_yslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION yslope_factor(cslope_l1:cslope_h1, cslope_l2:cslo
     &pe_h2, cslope_l3:cslope_h3)
      DOUBLE PRECISION uc_zslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION lc_zslope(cslope_l1:cslope_h1, cslope_l2:cslope_h
     &2, cslope_l3:cslope_h3,nvar)
      DOUBLE PRECISION zslope_factor(cslope_l1:cslope_h1, cslope_l2:cslo
     &pe_h2, cslope_l3:cslope_h3)
      DOUBLE PRECISION alpha(cslope_l1:cslope_h1, cslope_l2:cslope_h2, c
     &slope_l3:cslope_h3,nvar)
      DOUBLE PRECISION cmax(cslope_l1:cslope_h1, cslope_l2:cslope_h2, cs
     &lope_l3:cslope_h3,nvar)
      DOUBLE PRECISION cmin(cslope_l1:cslope_h1, cslope_l2:cslope_h2, cs
     &lope_l3:cslope_h3,nvar)
      DOUBLE PRECISION fvcx(fvcb_l1:fvcb_h1)
      DOUBLE PRECISION fvcy(fvcb_l2:fvcb_h2)
      DOUBLE PRECISION fvcz(fvcb_l3:fvcb_h3)
      DOUBLE PRECISION voffx(fvcb_l1:fvcb_h1)
      DOUBLE PRECISION voffy(fvcb_l2:fvcb_h2)
      DOUBLE PRECISION voffz(fvcb_l3:fvcb_h3)       
      DOUBLE PRECISION cvcx(cvcb_l1:cvcb_h1)
      DOUBLE PRECISION cvcy(cvcb_l2:cvcb_h2)
      DOUBLE PRECISION cvcz(cvcb_l3:cvcb_h3)

      integer n 
      integer i, ic
      integer j, jc
      integer k, kc
      DOUBLE PRECISION cen, forw, back, slp
      DOUBLE PRECISION factorn, denom
      DOUBLE PRECISION fxcen, cxcen, fycen, cycen, fzcen, czcen
      DOUBLE PRECISION corr_fact,orig_corr_fact
      DOUBLE PRECISION dummy_fine
      logical xok, yok, zok
      integer ncbx, ncby, ncbz
      integer ioff,joff,koff

      integer voff_lo(3), voff_hi(3)

      ncbx = cslopehi(1)-cslopelo(1)+1
      ncby = cslopehi(2)-cslopelo(2)+1
      ncbz = cslopehi(3)-cslopelo(3)+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      zok = (ncbz .ge. 2)

      voff_lo(1) = cslopelo(1) * lratiox
      voff_lo(2) = cslopelo(2) * lratioy
      voff_lo(3) = cslopelo(3) * lratioz
      voff_hi(1) = (cslopehi(1)+1) * lratiox - 1
      voff_hi(2) = (cslopehi(2)+1) * lratioy - 1
      voff_hi(3) = (cslopehi(3)+1) * lratioz - 1

      do k = voff_lo(3),voff_hi(3)
        kc = (k+lratioz*iabs(k))/lratioz-iabs(k)
        fzcen = 0.5D0*(fvcz(k)+fvcz(k+1))
        czcen = 0.5D0*(cvcz(kc)+cvcz(kc+1))
        voffz(k) = (fzcen-czcen)/(cvcz(kc+1)-cvcz(kc))
      end do
      do j = voff_lo(2),voff_hi(2)
        jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
        fycen = 0.5D0*(fvcy(j)+fvcy(j+1))
        cycen = 0.5D0*(cvcy(jc)+cvcy(jc+1))
        voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
      end do
      do i = voff_lo(1),voff_hi(1)
         ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
         fxcen = 0.5D0*(fvcx(i)+fvcx(i+1))
         cxcen = 0.5D0*(cvcx(ic)+cvcx(ic+1))
         voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
      end do

      do n = 1, nvar 
c
c     Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
c
          do k = cslopelo(3),cslopehi(3)
            do j = cslopelo(2),cslopehi(2)
              do i = cslopelo(1), cslopehi(1)
                alpha(i,j,k,n) = 1.d0
                cmax(i,j,k,n) = crse(i,j,k,n)
                cmin(i,j,k,n) = crse(i,j,k,n)
                do koff = -1,1
                do joff = -1,1
                do ioff = -1,1
                  cmax(i,j,k,n) = max(cmax(i,j,k,n),crse(i+ioff,j+joff,k
     &+koff,n))
                  cmin(i,j,k,n) = min(cmin(i,j,k,n),crse(i+ioff,j+joff,k
     &+koff,n))
                end do
                end do
                end do
              end do
            end do
          end do
      end do
c
c     Computed unlimited and limited slopes
c
      do n = 1, nvar 

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_xslope(i,j,k,n) = 0.5D0*(crse(i+1,j,k,n)-crse(i-1,j,k
     &,n))
                cen  = uc_xslope(i,j,k,n)
                forw = 2.0D0*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_xslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
            end do
          end do

          if (bc(1,1,n) .eq. 3 .or. bc(1,1,n).eq.4) then
            i = cslopelo(1)
            if (xok) then
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = -16.0D0/15.0D0*crse(i-1,j,k,n) + 
     &0.5D0*crse(i,j,k,n)+ 0.66666666666666667D0*crse(i+1,
     &j,k,n) - 0.1D0*crse(i+2,j,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = 0.25D0 * (crse(i+1,j,k,n) + 5.0D0
     &*crse(i,j,k,n) - 6.0D0*crse(i-1,j,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                cen  = uc_xslope(i,j,k,n)
                forw = 2.0D0*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_xslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bc(1,2,n) .eq. 3 .or. bc(1,2,n).eq.4) then
            i = cslopehi(1)
            if (xok) then
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = 16.0D0/15.0D0*crse(i+1,j,k,n) - 0
     &.5D0*crse(i,j,k,n)- 0.66666666666666667D0*crse(i-1,j
     &,k,n) + 0.1D0*crse(i-2,j,k,n)
                end do
              end do
            else 
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = -0.25D0 * (crse(i-1,j,k,n) + 5.0D
     &0*crse(i,j,k,n) - 6.0D0*crse(i+1,j,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                cen  = uc_xslope(i,j,k,n)
                forw = 2.0D0*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_xslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,k,n) = 0.5D0*(crse(i,j+1,k,n)-crse(i,j-1,k
     &,n))
                cen  = uc_yslope(i,j,k,n)
                forw = 2.0D0*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_yslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
               end do
            end do
          end do

          if (bc(2,1,n) .eq. 3 .or. bc(2,1,n).eq.4) then
            j = cslopelo(2)
            if (yok) then
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = -16.0D0/15.0D0*crse(i,j-1,k,n) + 
     &0.5D0*crse(i,j,k,n)+ 0.66666666666666667D0*crse(i,j+
     &1,k,n) - 0.1D0*crse(i,j+2,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = 0.25D0 * (crse(i,j+1,k,n) + 5.0D0
     &*crse(i,j,k,n) - 6.0D0*crse(i,j-1,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,k,n)
                forw = 2.0D0*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_yslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bc(2,2,n) .eq. 3 .or. bc(2,2,n).eq.4) then
            j = cslopehi(2)
            if (yok) then
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = 16.0D0/15.0D0*crse(i,j+1,k,n) - 0
     &.5D0*crse(i,j,k,n)- 0.66666666666666667D0*crse(i,j-1
     &,k,n) + 0.1D0*crse(i,j-2,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = -0.25D0 * (crse(i,j-1,k,n) + 5.0D
     &0*crse(i,j,k,n) - 6.0D0*crse(i,j+1,k,n) )
                end do
              end do
            endif
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  cen  = uc_yslope(i,j,k,n)
                  forw = 2.0D0*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = 2.0D0*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                  lc_yslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
                end do
              end do
          end if

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_zslope(i,j,k,n) = 0.5D0*(crse(i,j,k+1,n)-crse(i,j,k-1
     &,n))
                cen  = uc_zslope(i,j,k,n)
                forw = 2.0D0*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_zslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
              end do
            end do
          end do

          if (bc(3,1,n) .eq. 3 .or. bc(3,1,n).eq.4) then
            k = cslopelo(3)
            if (zok) then
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = -16.0D0/15.0D0*crse(i,j,k-1,n) + 
     &0.5D0*crse(i,j,k,n)+ 0.66666666666666667D0*crse(i,j,
     &k+1,n) - 0.1D0*crse(i,j,k+2,n)
                end do
              end do
            else
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = 0.25D0 * (crse(i,j,k+1,n) + 5.0D0
     &*crse(i,j,k,n) - 6.0D0*crse(i,j,k-1,n) )
                end do
              end do
            endif
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_zslope(i,j,k,n)
                forw = 2.0D0*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_zslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bc(3,2,n) .eq. 3 .or. bc(3,2,n).eq.4) then
            k = cslopehi(3)
            if (zok) then
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = 16.0D0/15.0D0*crse(i,j,k+1,n) - 0
     &.5D0*crse(i,j,k,n)- 0.66666666666666667D0*crse(i,j,k
     &-1,n) + 0.1D0*crse(i,j,k-2,n)
               end do
              end do
            else
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = -0.25D0 * (crse(i,j,k-1,n) + 5.0D
     &0*crse(i,j,k,n) - 6.0D0*crse(i,j,k+1,n) )
               end do
              end do
            endif
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_zslope(i,j,k,n)
                forw = 2.0D0*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = 2.0D0*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_zslope(i,j,k,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
            end do
          end if

       end do

       if (lim_slope.eq.0) then
c
c        Do the interpolation using unlimited slopes.
c
          do n = 1, nvar
             do k = fblo(3), fbhi(3)
               kc = (k+lratioz*iabs(k))/lratioz-iabs(k)
               do j = fblo(2), fbhi(2)
                  jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
                  do i = fblo(1), fbhi(1)
                     ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
                     fine(i,j,k,n) = crse(ic,jc,kc,n) + voffx(i)*uc_xslo
     &pe(ic,jc,kc,n)+ voffy(j)*uc_yslope(ic,jc,kc,n)+ v
     &offz(k)*uc_zslope(ic,jc,kc,n)
                  end do
               end do
             end do
          end do

       else

         if (lin_limit.eq.1)then
c
c     compute linear limited slopes
c     Note that the limited and the unlimited slopes
c     have the same sign, and it is assumed that they do.
c
c     compute slope factors
c
          xslope_factor = 1.0D0
          yslope_factor = 1.0D0
          zslope_factor = 1.0D0

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  denom = uc_xslope(i,j,k,n)
                  denom = merge(denom,1.0D0,denom.ne.0.0D0)
                  factorn = lc_xslope(i,j,k,n)/denom
                  factorn = merge(1.0D0,factorn,denom.eq.0.0D0)
                  xslope_factor(i,j,k) = min(xslope_factor(i,j,k),factor
     &n)

                  denom = uc_yslope(i,j,k,n)
                  denom = merge(denom,1.0D0,denom.ne.0.0D0)
                  factorn = lc_yslope(i,j,k,n)/denom
                  factorn = merge(1.0D0,factorn,denom.eq.0.0D0)
                  yslope_factor(i,j,k) = min(yslope_factor(i,j,k),factor
     &n)

                  denom = uc_zslope(i,j,k,n)
                  denom = merge(denom,1.0D0,denom.ne.0.0D0)
                  factorn = lc_zslope(i,j,k,n)/denom
                  factorn = merge(1.0D0,factorn,denom.eq.0.0D0)
                  zslope_factor(i,j,k) = min(zslope_factor(i,j,k),factor
     &n)
                end do
              end do
            end do
          end do
c
c         Compute linear limited slopes
c
          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  lc_xslope(i,j,k,n) = xslope_factor(i,j,k)*uc_xslope(i,
     &j,k,n)
                  lc_yslope(i,j,k,n) = yslope_factor(i,j,k)*uc_yslope(i,
     &j,k,n)
                  lc_zslope(i,j,k,n) = zslope_factor(i,j,k)*uc_zslope(i,
     &j,k,n)
                end do
              end do
            end do
          end do

         else
c
c         Limit slopes so as to not introduce new maxs or mins.
c
            do n = 1,nvar
               do k = voff_lo(3),voff_hi(3)
                  kc = (k+lratioz*iabs(k))/lratioz-iabs(k)
                  do j = voff_lo(2),voff_hi(2)
                     jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
                     do i = voff_lo(1),voff_hi(1)
                        ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
                        orig_corr_fact = voffx(i)*lc_xslope(ic,jc,kc,n)+
     & voffy(j)*lc_yslope(ic,jc,kc,n)+ voffz(k)*lc_z
     &slope(ic,jc,kc,n)
                        dummy_fine = crse(ic,jc,kc,n) + orig_corr_fact
                        if ((dummy_fine .gt. cmax(ic,jc,kc,n)) .and.(abs
     &(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,kc
     &,n)))) then
                           corr_fact = (cmax(ic,jc,kc,n) - crse(ic,jc,kc
     &,n)) / orig_corr_fact
                           alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),cor
     &r_fact)
                        endif

                        if ((dummy_fine .lt. cmin(ic,jc,kc,n)) .and.(abs
     &(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,kc
     &,n)))) then
                           corr_fact = (cmin(ic,jc,kc,n) - crse(ic,jc,kc
     &,n)) / orig_corr_fact
                           alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),cor
     &r_fact)
                        endif

                     end do
                  end do
               end do
            end do
         end if
c
c       Do the interpolation with limited slopes.
c
        do n = 1, nvar
          do k = fblo(3), fbhi(3)
            kc = (k+lratioz*iabs(k))/lratioz-iabs(k)
            do j = fblo(2), fbhi(2)
              jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
              do i = fblo(1), fbhi(1)
                ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
                fine(i,j,k,n) = crse(ic,jc,kc,n) + alpha(ic,jc,kc,n) *( 
     &voffx(i)*lc_xslope(ic,jc,kc,n)+voffy(j)*lc_yslope(ic,j
     &c,kc,n)+voffz(k)*lc_zslope(ic,jc,kc,n) )
              end do
            end do
          end do
        end do

      end if

      end

      subroutine cqinterp (fine, fine_l1, fine_l2, fine_l3, fine_h1, fin
     &e_h2, fine_h3, fb_l1, fb_l2, fb_l3, fb_h1, fb_h2, fb_h3,nvar, lr
     &atiox, lratioy, lratioz, crse,clo, chi, cb_l1, cb_l2, cb_l3, cb_
     &h1, cb_h2, cb_h3,fslo, fshi, cslope, clen, fslope, fdat,flen, vo
     &ff, bc, limslope,fvcx, fvcy, fvcz, cvcx, cvcy, cvcz,actual_comp,
     & actual_state)

      implicit none

      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer fb_l1, fb_l2, fb_l3, fb_h1, fb_h2, fb_h3
      integer cb_l1, cb_l2, cb_l3, cb_h1, cb_h2, cb_h3
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      integer actual_comp,actual_state
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3,nvar)
      DOUBLE PRECISION crse(clo:chi, nvar)
      DOUBLE PRECISION cslope(clo:chi, 3)
      DOUBLE PRECISION fslope(flen, 3)
      DOUBLE PRECISION fdat(flen)
      DOUBLE PRECISION voff(flen)
      DOUBLE PRECISION fvcx(fb_l1:fb_h1+1)
      DOUBLE PRECISION fvcy(fb_l2:fb_h2+1)
      DOUBLE PRECISION fvcz(fb_l3:fb_h3+1)
      DOUBLE PRECISION cvcx(cb_l1:cb_h1+1)
      DOUBLE PRECISION cvcy(cb_l2:cb_h2+1)
      DOUBLE PRECISION cvcz(cb_l3:cb_h3+1)

      call bl_abort('QUADRATIC INTERP NOT IMPLEMEMNTED IN 3-D')

      end

c ::: 
c ::: --------------------------------------------------------------
c ::: pcinterp:  cell centered piecewise constant interpolation
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3   =>  (const)  index limits of coarse grid
c ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
c ::: 
c ::: longdir      =>  (const)  which index direction is longest (1 or 2)
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine pcinterp (crse,crse_l1, crse_l2, crse_l3, crse_h1, crse
     &_h2, crse_h3,cblo,cbhi,fine,fine_l1, fine_l2, fine_l3, fine_h1, 
     &fine_h2, fine_h3,fblo,fbhi,longdir,lratiox,lratioy,lratioz,nvar,
     &ftmp, ftmp_lo, ftmp_hi,actual_comp, actual_state)

      implicit none

      integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      integer cblo(3), cbhi(3)
      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer fblo(3), fbhi(3)
      integer nvar, lratiox, lratioy, lratioz, longdir
      integer ftmp_lo, ftmp_hi
      integer actual_comp,actual_state
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3, nvar)
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3, nvar)
      DOUBLE PRECISION  ftmp(ftmp_lo:ftmp_hi)

      integer i, j, k, ic, jc, kc, ioff, joff, koff, n
      integer istrt, istop, jstrt, jstop, kstrt, kstop
      integer ilo, ihi, jlo, jhi, klo, khi

      if (longdir .eq. 1) then
         do n = 1, nvar
         do kc = cblo(3), cbhi(3)
            kstrt = kc*lratioz
            kstop = (kc+1)*lratioz - 1
            klo = max(fblo(3),kstrt)
            khi = min(fbhi(3),kstop)
            do jc = cblo(2), cbhi(2)

c	       ::::: fill strip in i direction
               do ioff = 0, lratiox-1
                  do ic = cblo(1), cbhi(1)
                     i = lratiox*ic + ioff
                     ftmp(i) = crse(ic,jc,kc,n)
                  end do
               end do

c	       ::::: stuff into fine array
               jstrt = jc*lratioy
               jstop = (jc+1)*lratioy - 1
               jlo = max(fblo(2),jstrt)
               jhi = min(fbhi(2),jstop)
               do k = klo, khi
               do j = jlo, jhi
               do i = fblo(1), fbhi(1)
                  fine(i,j,k,n) = ftmp(i)
               end do
               end do
               end do
            end do
         end do
         end do
      else if (longdir.eq.2) then
         do n = 1, nvar
         do kc = cblo(3), cbhi(3)
            kstrt = kc*lratioz
            kstop = (kc+1)*lratioz - 1
            klo = max(fblo(3),kstrt)
            khi = min(fbhi(3),kstop)
            do ic = cblo(1), cbhi(1)

c	       ::::: fill strip in j direction
               do joff = 0, lratioy-1
                  do jc = cblo(2), cbhi(2)
                     j = lratioy*jc + joff
                     ftmp(j) = crse(ic,jc,kc,n)
                  end do
               end do

c	       ::::: stuff into fine array
               istrt = ic*lratiox
               istop = (ic+1)*lratiox - 1
               ilo = max(fblo(1),istrt)
               ihi = min(fbhi(1),istop)
               do k = klo, khi
               do i = ilo, ihi
               do j = fblo(2), fbhi(2)
                  fine(i,j,k,n) = ftmp(j)
               end do
               end do
               end do
            end do
         end do
         end do
      else
         do n = 1, nvar
         do ic = cblo(1), cbhi(1)
            istrt = ic*lratiox
            istop = (ic+1)*lratiox - 1
            ilo = max(fblo(1),istrt)
            ihi = min(fbhi(1),istop)
            do jc = cblo(2), cbhi(2)

c	       ::::: fill strip in k direction
               do koff = 0, lratioz-1
                  do kc = cblo(3), cbhi(3)
                     k = lratioz*kc + koff
                     ftmp(k) = crse(ic,jc,kc,n)
                  end do
               end do

c	       ::::: stuff into fine array
               jstrt = jc*lratioy
               jstop = (jc+1)*lratioy - 1
               jlo = max(fblo(2),jstrt)
               jhi = min(fbhi(2),jstop)
               do i = ilo, ihi
               do j = jlo, jhi
               do k = fblo(3), fbhi(3)
                  fine(i,j,k,n) = ftmp(k)
               end do
               end do
               end do
            end do
         end do
         end do
      end if

      end

c ::: 
c ::: --------------------------------------------------------------
c ::: protect_interp:   redo interpolation if the result of linccinterp
c ::: generates under- or overshoots.
c ::: 
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: cblo,cbhi    =>  (const)  coarse equivalent of fblo,fbhi
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: clo,chi      =>  (const)  index limits of crse grid
c :::
c ::: --------------------------------------------------------------
c ::: 
      subroutine printerp (fine, fine_l1, fine_l2, fine_l3, fine_h1, fin
     &e_h2, fine_h3, fblo, fbhi, crse, crse_l1, crse_l2, crse_l3, crse
     &_h1, crse_h2, crse_h3, cblo, cbhi,fine_state, state_l1, state_l2
     &, state_l3, state_h1, state_h2, state_h3, nvar, lratiox,lratioy,
     &lratioz, bc)

      implicit none

      integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
      integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
      integer state_l1, state_l2, state_l3, state_h1, state_h2, state_h3
      integer fblo(3), fbhi(3)
      integer cblo(3), cbhi(3)
      integer lratiox, lratioy, lratioz, nvar
      integer bc(3,2,nvar)
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:fi
     &ne_h3,nvar)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:cr
     &se_h3, nvar)
      DOUBLE PRECISION fine_state(state_l1:state_h1, state_l2:state_h2, 
     &state_l3:state_h3, nvar)

      integer rMAX
      parameter (rMAX = 16)
      DOUBLE PRECISION alpha, sumN, sumP, crseTot, negVal, posVal
      DOUBLE PRECISION sum_fine_new,sum_fine_old
      DOUBLE PRECISION orig_fine(0:rMAX-1,0:rMAX-1,0:rMAX-1)
      integer redo_me
      integer ilo,ihi,jlo,jhi,klo,khi
      integer i,j,k,ic,jc,kc,n
      integer numFineCells
      integer icase

      if (MAX(lratiox,lratioy,lratioz).gt.rMAX) then
         print *,'rMAX in INTERP_3D::FORT_PROTECT_INTERP must be >= ',MA
     &X(lratiox,lratioy,lratioz)
         call bl_abort(" ")
      endif

      do kc = cblo(3), cbhi(3)
      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)
         klo = max(lratioz*kc            ,fine_l3)
         khi = min(lratioz*kc+(lratioz-1),fine_h3)

         do n = 2, nvar-1

            redo_me = 0
            do k = klo,khi
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,k,n)+fine(i,j,k,n)) .lt. 0.d0) redo_m
     &e = 1
            enddo
            enddo
            enddo
c
c ****************************************************************************************
c
c           If all the fine values are non-negative after the original interpolated 
c            correction, then we do nothing here.
c
c           If any of the fine values are negative after the original interpolated
c            correction, then we do our best.
c
c           Special cases:
c
c             1) Coarse correction > 0, and fine_state has some cells with 
c                negative values which will be filled before adding to the other cells.
c                Use the correction to bring negative cells to 0.0D0, then
c                distribute the remaining positive proportionally.
c
c             2) Coarse correction > 0, and correction can not make them all
c                positive.  Add correction only to the negative cells, in proportion
c                to their magnitude.
c
c             3) Coarse correction < 0, and fine_state DOES NOT have enough
c                  have enough positive state to absorb it.  Here we bring
c                  all the positive fine cells to 0.0D0 then distribute the remaining
c                  negative amount in such a way as to make them all as close to the
c                  same negative value as possible.
c
c             4) Coarse correction < 0, fine_state has enough
c                  positive state to absorb it without making any fine 
c                  cells negative, BUT fine_state+fine is currently negative
c                  in at least 1.0D0 fine cell.  Here just take a constant percentage
c                  away from each positive and don't touch the negatives.
c
c             crseTot = sum of all interpolated values of the correction,
c                       which is equivalent to the coarse correction * ratio**3
c             SumN = sum of all negative values of fine_state
c             SumP = sum of all positive values of fine_state
c
c ****************************************************************************************
c

            if (redo_me .eq. 1) then

               icase = 0
               sum_fine_old = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_old = sum_fine_old + fine(i,j,k,n)
                  orig_fine(i-ilo,j-jlo,k-klo) = fine(i,j,k,n)
               enddo
               enddo
               enddo

               crseTot = sum_fine_old
               numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1)

               sumN = 0.0D0
               sumP = 0.0D0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  if (fine_state(i,j,k,n) .le. 0.d0) then
                    sumN = SumN + fine_state(i,j,k,n)
                  else
                    sumP = sumP + fine_state(i,j,k,n)
                  endif
               enddo
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
c              Here we want to fill in the negative values first, then add
c                the remaining positive proportionally.

                   icase = 1
                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .le. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

                   if (sumP .gt. 0.d0) then

                    alpha = (crseTot - abs(sumN)) / sumP

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       if (fine_state(i,j,k,n) .ge. 0.d0) then
                         fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                       endif
                    enddo
                    enddo
                    enddo

                  else

                    posVal = (crseTot - abs(sumN)) / float(numFineCells)

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       fine(i,j,k,n) = fine(i,j,k,n) + posVal
                    enddo
                    enddo
                    enddo

                  endif

               endif

               if (crseTot .gt. 0.d0. and. crseTot .lt. abs(sumN)) then
c              Here we don't have enough positive correction to fill all the
c                negative values of state, so we just try to fill them proportionally
c                and don't add any correction to the states already positive.

                   icase = 2
                   alpha = crseTot / abs(sumN)

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = alpha * abs(fine_state(i,j,k,n))
                      else 
                        fine(i,j,k,n) = 0.d0
                      endif
                   enddo
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0. and. abs(crseTot) .gt. sumP) then
c              Here we don't have enough positive states to absorb all the
c                negative correction, so we want to end up with all the fine
c                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/float(numFineCells)

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,k,n) = negVal - fine_state(i,j,k,n)
                   enddo
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP.and. (
     &sumP+sumN+crseTot) .gt. 0.d0) then
c              Here we have enough positive states to absorb all the
c                negative correction *and* redistribute to make negative cells
c                positive. 

                   icase = 4
                   alpha = (crseTot + sumN) / sumP

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif  
                   enddo
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0. and. abs(crseTot) .lt. sumP.and. (
     &sumP+sumN+crseTot) .le. 0.d0) then
c              Here we have enough positive states to absorb all the
c                negative correction, but not to fix the states already negative. 
c                We bring all the positive states to 0.0D0, and use whatever 
c                remaining positiveness from the states to help the negative states.

                   icase = 5
                   alpha = (crseTot + sumP) / sumN

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .gt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else 
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

               endif

               sum_fine_new = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_new = sum_fine_new + fine(i,j,k,n)
               enddo
               enddo
               enddo

               if (abs(sum_fine_new - sum_fine_old) .gt. 1.e-8) then
                  print *,' '
                  print *,'PROTECT_INTERP: BLEW CONSERVATION with ICASE 
     &= ',icase
                  print *,'AT COARSE CELL ',ic,jc,kc,' AND COMPONENT ',n
                  print *,'NEW SUM / OLD SUM ',sum_fine_new, sum_fine_ol
     &d
                  print *,'CRSETOT ',crseTot
                  print *,'SUMP SUMN ',sumP,sumN
                  do k = klo,khi
                  do j = jlo,jhi
                  do i = ilo,ihi
                     print *,'FINE OLD NEW ',i,j,k,orig_fine(i-ilo,j-jlo
     &,k-klo),fine(i,j,k,n), fine_state(i,j,k,n)
                  enddo
                  enddo
                  enddo
               endif

c              do k = klo,khi
c              do j = jlo,jhi
c              do i = ilo,ihi
c                 if ((fine_state(i,j,k,n) + fine(i,j,k,n)) .lt. 0.d0) then
c                    print *,'STILL NEGATIVE AT ',i,j,k,n
c                    print *,'AT COARSE CELL ',ic,jc,kc
c                    print *,'FINE STATE ',fine_state(i,j,k,n)
c                    print *,'FINE CORRECTION ',fine(i,j,k,n)
c                    print *,'CRSETOT ',crseTot
c                    print *,'SUMN / SUMP ',sumN, sumP
c                    print *,' '
c                 endif
c              enddo
c              enddo
c              enddo
c           End (if redo .eq. 1)
            endif

         enddo

c     Set sync for density (n=1) to sum of spec sync (2:nvar-1)
         do k = klo,khi
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,k,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,k,1) = fine(i,j,k,1) + fine(i,j,k,n)
            enddo
         enddo
         enddo
         enddo

c     End of coarse index loops
      enddo
      enddo
      enddo
      end

c ::: 
c ::: --------------------------------------------------------------
c ::: quartinterp: quartic conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio[xyz]  =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data
c ::: clo,chi      =>  (const)  index limits of crse grid
c ::: cblo,cbhi    =>  (const)  coarse grid region containing fblo,fbhi and widen by 2 or 4 cells
c :::
c ::: cb2lo,cb2hi  =>  (const)  coarse grid region containing fblo,fbhi
c ::: fb2lo,fb2hi  =>  (const)  fine version of cb2. It could be wider than fb
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: ctmp         =>  2-D temp array
c ::: ctmp2        =>  2-D temp array
c ::: --------------------------------------------------------------
c ::: 
       subroutine quartinterp (fine, fine_l1, fine_l2, fine_l3, fine_h1,
     & fine_h2, fine_h3, fblo, fbhi, fb2lo, fb2hi,crse, crse_l1, crse
     &_l2, crse_l3, crse_h1, crse_h2, crse_h3, cblo, cbhi, cb2lo, cb2
     &hi,nvar, lratiox, lratioy, lratioz,ftmp, ctmp, ctmp2,bc,actual_
     &comp,actual_state)

       implicit none

       integer fine_l1, fine_l2, fine_l3, fine_h1, fine_h2, fine_h3
       integer crse_l1, crse_l2, crse_l3, crse_h1, crse_h2, crse_h3
       integer fblo(3), fbhi(3), fb2lo(3), fb2hi(3)
       integer cblo(3), cbhi(3), cb2lo(3), cb2hi(3)
       integer nvar,lratiox,lratioy,lratioz
       integer bc(3,2,nvar)
       integer actual_comp,actual_state
       DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2, fine_l3:f
     &ine_h3,nvar)
       DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, crse_l3:c
     &rse_h3,nvar)
       DOUBLE PRECISION ftmp(fb2lo(1):fb2hi(1))
       DOUBLE PRECISION ctmp(cblo(1):cbhi(1),0:lratioy-1)
       DOUBLE PRECISION ctmp2(cblo(1):cbhi(1),cblo(2):cbhi(2),0:lratioz-
     &1)

c      Local variables
       integer i,j,k,ii,jj,kk,n,iry,irz
       DOUBLE PRECISION cL(-2:2)
c       DOUBLE PRECISION cR(-2:2)
       data cL/ -0.01171875D0, 0.0859375D0, 0.5d0, -0.0859375D0,0.011718
     &75D0 /
c       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0,
c     $          -0.01171875D0 /

       if (lratiox.eq.2 .and. lratioy.eq.2 .and. lratioz.eq.2) then

          do n = 1, nvar
          do k = cb2lo(3), cb2hi(3)

             do j = cblo(2), cbhi(2)
             do i = cblo(1), cbhi(1)
                ctmp2(i,j,0) = 2.d0*(cL(-2)*crse(i,j,k-2,n) + cL(-1)*crs
     &e(i,j,k-1,n)+ cL( 0)*crse(i,j,k ,n)+ cL( 1)*crse(i,j,k
     &+1,n)+ cL( 2)*crse(i,j,k+2,n))
                ctmp2(i,j,1) = 2.d0*crse(i,j,k,n) - ctmp2(i,j,0)
c$$$                ctmp2(i,j,1) = 2.d0*(cR(-2)*crse(i,j,k-2,n) 
c$$$     $               +               cR(-1)*crse(i,j,k-1,n)
c$$$     $               +               cR( 0)*crse(i,j,k  ,n)
c$$$     $               +               cR( 1)*crse(i,j,k+1,n)
c$$$     $               +               cR( 2)*crse(i,j,k+2,n))
             enddo
             enddo

             do irz = 0, 1
                kk = k*2+irz
                if (kk.ge.fblo(3) .and. kk.le.fbhi(3)) then

                   do j = cb2lo(2), cb2hi(2)

                      do i = cblo(1), cbhi(1)
                         ctmp(i,0) = 2.d0*(cL(-2)*ctmp2(i,j-2,irz) + cL(
     &-1)*ctmp2(i,j-1,irz)+ cL( 0)*ctmp2(i,j ,irz)+
     & cL( 1)*ctmp2(i,j+1,irz)+ cL( 2)*ctmp2(i,j+2,
     &irz))
                         ctmp(i,1) = 2.d0*ctmp2(i,j,irz) - ctmp(i,0)
c$$$                         ctmp(i,1) = 2.d0*(cR(-2)*ctmp2(i,j-2,irz) 
c$$$     $                        +            cR(-1)*ctmp2(i,j-1,irz)
c$$$     $                        +            cR( 0)*ctmp2(i,j  ,irz)
c$$$     $                        +            cR( 1)*ctmp2(i,j+1,irz)
c$$$     $                        +            cR( 2)*ctmp2(i,j+2,irz))
                      enddo

                      do iry = 0, 1
                         jj = j*2+iry

                         if (jj.ge.fblo(2).and.jj.le.fbhi(2)) then
                            do i = cb2lo(1), cb2hi(1)
                               ii = 2*i
                               ftmp(ii ) = 2.d0*(cL(-2)*ctmp(i-2,iry) + 
     &cL(-1)*ctmp(i-1,iry)+ cL( 0)*ctmp(i ,ir
     &y)+ cL( 1)*ctmp(i+1,iry)+ cL( 2)*ctmp(i
     &+2,iry))
                               ftmp(ii+1) = 2.d0*ctmp(i,iry) - ftmp(ii)
c$$$                               ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) 
c$$$     $                              +             cR(-1)*ctmp(i-1,iry)
c$$$     $                              +             cR( 0)*ctmp(i  ,iry)
c$$$     $                              +             cR( 1)*ctmp(i+1,iry)
c$$$     $                              +             cR( 2)*ctmp(i+2,iry))
                            enddo
                            do ii = fblo(1), fbhi(1)
                               fine(ii,jj,kk,n) = ftmp(ii)
                            enddo
                         endif  ! if (jj.ge.......
                      enddo  ! do iry

                   enddo  ! do j

                endif  ! if (kk.ge.......
             enddo  ! do irz

          enddo  ! do k
          enddo  ! do n

       else if (lratiox.eq.4 .and. lratioy.eq.4 .and. lratioz.eq.4) then
c      todo
          write(6,*) 'FORT_QUARTINTERP: refinement ratio = 4 TODO'
          stop
       else
          write(6,*) 'FORT_QUARTINTERP: unsupported refinement ratio'
          stop
       endif

       end
