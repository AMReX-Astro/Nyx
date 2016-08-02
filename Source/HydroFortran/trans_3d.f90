module transverse_module
 
  use bl_constants_module
 
  implicit none
 
contains

      !===========================================================================
      ! transx1 -- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transx1(qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,k3d)

      ! Note that what we call ilo here is ilo = lo(1)
      ! Note that what we call ihi here is ihi = hi(1)
      ! Note that what we call jlo here is jlo = lo(2) - 1
      ! Note that what we call jhi here is jhi = hi(2) + 1

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      double precision cdtdx

      integer i, j
      integer n, nq

      double precision rrnew, rr
      double precision rrry, rrly
      double precision rury, ruly
      double precision rvry, rvly
      double precision rwry, rwly
      double precision ekenry, ekenly
      double precision pnewry, pnewly
      double precision rery, rely
      double precision rrnewry, rrnewly
      double precision runewry, runewly
      double precision rvnewry, rvnewly
      double precision rwnewry, rwnewly
      double precision renewry, renewly
      double precision rhoekenry, rhoekenly
      double precision compn, compu
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n)) 

               if (j.ge.jlo+1) then
                  rr = qyp(i,j,kc,QRHO)
                  rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                  compu = rr*qyp(i,j,kc,nq) - compn
                  qypo(i,j,kc,nq) = compu/rrnew
               end if

               if (j.le.jhi-1) then
                  rr = qym(i,j+1,kc,QRHO)
                  rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                  compu = rr*qym(i,j+1,kc,nq) - compn
                  qymo(i,j+1,kc,nq) = compu/rrnew
               end if

            enddo
         enddo
      enddo

      ! NOTE: it is better *not* to protect against small density in this routine

      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvx(i+1,j,kc)
            pgm = pgdnvx(i,j,kc)
            ugp = ugdnvx(i+1,j,kc)
            ugm = ugdnvx(i,j,kc)

            if (j.ge.jlo+1) then
               ! Convert to conservation form
               rrry = qyp(i,j,kc,QRHO)
               rury = rrry*qyp(i,j,kc,QU)
               rvry = rrry*qyp(i,j,kc,QV)
               rwry = rrry*qyp(i,j,kc,QW)
               ekenry = HALF*rrry*(qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
               rery = qyp(i,j,kc,QREINT) + ekenry

               ! Add transverse terms
               rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO ) - fx(i,j,kc,URHO ))
               runewry = rury - cdtdx*(fx(i+1,j,kc,UMX  ) - fx(i,j,kc,UMX  ))
               rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY  ) - fx(i,j,kc,UMY  ))
               rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ  ) - fx(i,j,kc,UMZ  ))
               renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewry .lt. ZERO) then
                 rrnewry = rrry
                 runewry = rury
                 rvnewry = rvry
                 rwnewry = rwry
                 renewry = rery
               end if
            end if

            if (j.le.jhi-1) then
               rrly = qym(i,j+1,kc,QRHO)
               ruly = rrly*qym(i,j+1,kc,QU)
               rvly = rrly*qym(i,j+1,kc,QV)
               rwly = rrly*qym(i,j+1,kc,QW)
               ekenly = HALF*rrly* &
                    (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
               rely = qym(i,j+1,kc,QREINT) + ekenly

               ! Add transverse terms
               rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO ) - fx(i,j,kc,URHO ))
               runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX  ) - fx(i,j,kc,UMX  ))
               rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY  ) - fx(i,j,kc,UMY  ))
               rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ  ) - fx(i,j,kc,UMZ  ))
               renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewly .lt. ZERO) then
                 rrnewly = rrly
                 runewly = ruly
                 rvnewly = rvly
                 rwnewly = rwly
                 renewly = rely
               end if
            end if

            dup = pgp*ugp - pgm*ugm
            pav = HALF*(pgp+pgm)
            du = ugp-ugm

            ! Convert back to non-conservation form
            if (j.ge.jlo+1) then
               qypo(i,j,kc,QRHO) = rrnewry
               qypo(i,j,kc,QU) = runewry/qypo(i,j,kc,QRHO)
               qypo(i,j,kc,QV) = rvnewry/qypo(i,j,kc,QRHO)
               qypo(i,j,kc,QW) = rwnewry/qypo(i,j,kc,QRHO)
               rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,QRHO)

               qypo(i,j,kc,QREINT) = renewry - rhoekenry
               qypo(i,j,kc,QPRES) = qypo(i,j,kc,QREINT) * gamma_minus_1

               if (qypo(i,j,kc,QPRES) .lt. small_pres) then
                   pnewry = qyp(i,j  ,kc,QPRES) - cdtdx*(dup + pav*du*gamma_minus_1)
                   qypo(i,j,kc,QPRES ) = pnewry
                   qypo(i,j,kc,QREINT) = qypo(i,j,kc,QPRES) / gamma_minus_1
               end if
            end if

            if (j.le.jhi-1) then
               qymo(i,j+1,kc,QRHO) = rrnewly
               qymo(i,j+1,kc,QU) = runewly/qymo(i,j+1,kc,QRHO)
               qymo(i,j+1,kc,QV) = rvnewly/qymo(i,j+1,kc,QRHO)
               qymo(i,j+1,kc,QW) = rwnewly/qymo(i,j+1,kc,QRHO)
               rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,QRHO)

               qymo(i,j+1,kc,QREINT) = renewly - rhoekenly
               qymo(i,j+1,kc,QPRES) = qymo(i,j+1,kc,QREINT) * gamma_minus_1

               if (qymo(i,j+1,kc,QPRES) .lt. small_pres) then
                   pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*gamma_minus_1)
                   qymo(i,j+1,kc,QPRES ) = pnewly
                   qymo(i,j+1,kc,QREINT) = qymo(i,j+1,kc,QPRES) / gamma_minus_1
               end if
            end if

         enddo
      enddo

      end subroutine transx1

      !===========================================================================
      ! transx2 -- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transx2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      double precision cdtdx

      integer i, j
      integer n, nq

      double precision rrnew, rr
      double precision rrrz, rrlz
      double precision rurz, rulz
      double precision rvrz, rvlz
      double precision rwrz, rwlz
      double precision ekenrz, ekenlz
      double precision rerz, relz
      double precision pnewrz, pnewlz
      double precision rrnewrz, rrnewlz
      double precision runewrz, runewlz
      double precision rvnewrz, rvnewlz
      double precision rwnewrz, rwnewlz
      double precision renewrz, renewlz
      double precision rhoekenrz, rhoekenlz
      double precision compn, compu
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive
      
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

                compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

                rr = qzp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qzp(i,j,kc,nq) - compn
                qzpo(i,j,kc,nq) = compu/rrnew

                compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

                rr = qzm(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
                compu = rr*qzm(i,j,kc,nq) - compn
                qzmo(i,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo

      do j = jlo, jhi
          do i = ilo, ihi

             ! ************************************************************************
             ! Convert to conservation form
             rrrz =      qzp(i,j,kc,QRHO)
             rurz = rrrz*qzp(i,j,kc,QU)
             rvrz = rrrz*qzp(i,j,kc,QV)
             rwrz = rrrz*qzp(i,j,kc,QW)
             ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
             rerz = qzp(i,j,kc,QREINT) + ekenrz

             ! Add transverse terms
             rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO ) - fx(i,j,kc,URHO ))
             runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX  ) - fx(i,j,kc,UMX  ))
             rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY  ) - fx(i,j,kc,UMY  ))
             rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ  ) - fx(i,j,kc,UMZ  ))
             renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewrz .lt. ZERO) then
                 rrnewrz = rrrz
                 runewrz = rurz
                 rvnewrz = rvrz
                 rwnewrz = rwrz
                 renewrz = rerz
            end if

             ! Convert back to non-conservation form
             qzpo(i,j,kc,QRHO) = rrnewrz
             qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
             qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
             qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)

             qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz
             qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT) * gamma_minus_1

             if (qzpo(i,j,kc,QPRES) .lt. small_pres) then
                 pgp = pgdnvx(i+1,j,kc)
                 pgm = pgdnvx(i,j,kc)
                 ugp = ugdnvx(i+1,j,kc)
                 ugm = ugdnvx(i,j,kc)
                 dup = pgp*ugp - pgm*ugm
                 pav = HALF*(pgp+pgm)
                 du = ugp-ugm
                 pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*gamma_minus_1)
                 qzpo(i,j,kc,QPRES ) = pnewrz
                 qzpo(i,j,kc,QREINT) = qzpo(i,j,kc,QPRES) / gamma_minus_1
             end if
             ! ************************************************************************

             ! ************************************************************************
             ! Convert to conservation form
             rrlz =      qzm(i,j,kc,QRHO)
             rulz = rrlz*qzm(i,j,kc,QU)
             rvlz = rrlz*qzm(i,j,kc,QV)
             rwlz = rrlz*qzm(i,j,kc,QW)
             ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
             relz = qzm(i,j,kc,QREINT) + ekenlz

             ! Add transverse terms
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO ) - fx(i,j,km,URHO ))
             runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX  ) - fx(i,j,km,UMX  ))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY  ) - fx(i,j,km,UMY  ))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ  ) - fx(i,j,km,UMZ  ))
             renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewlz .lt. ZERO) then
                 rrnewlz = rrlz
                 runewlz = rulz
                 rvnewlz = rvlz
                 rwnewlz = rwlz
                 renewlz = relz
            end if

             ! Convert back to non-conservation form
             qzmo(i,j,kc,QRHO) = rrnewlz
             qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
             qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
             qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)

             qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz
             qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT) * gamma_minus_1

             if (qzmo(i,j,kc,QPRES) .lt. small_pres) then
                 pgp = pgdnvx(i+1,j,km)
                 pgm = pgdnvx(i,j,km)
                 ugp = ugdnvx(i+1,j,km)
                 ugm = ugdnvx(i,j,km)
                 dup = pgp*ugp - pgm*ugm
                 pav = HALF*(pgp+pgm)
                 du = ugp-ugm
                 pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*gamma_minus_1)
                 qzmo(i,j,kc,QPRES ) = pnewlz
                 qzmo(i,j,kc,QREINT) = qzmo(i,j,kc,QPRES) / gamma_minus_1
             end if
             ! ************************************************************************

          enddo
      enddo

      end subroutine transx2

      !===========================================================================
      ! transy1 -- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transy1(qxm,qxmo,qxp,qxpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,k3d)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)

      ! Note that cdtdx = dtdx/3.d0/a_half
      double precision cdtdy

      integer i, j
      integer n, nq

      double precision rrnew, rr
      double precision compn, compu
      double precision rrrx, rrlx
      double precision rurx, rulx
      double precision rvrx, rvlx
      double precision rwrx, rwlx
      double precision ekenrx, ekenlx
      double precision rerx, relx
      double precision rrnewrx, rrnewlx
      double precision runewrx, runewlx
      double precision rvnewrx, rvnewlx
      double precision rwnewrx, rwnewlx
      double precision renewrx, renewlx
      double precision pnewrx, pnewlx
      double precision rhoekenrx, rhoekenlx
      double precision pgp, pgm, ugp, ugm
      double precision du,dup,pav

      integer          :: ipassive
      
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               if (i.ge.ilo+1) then
                  rr = qxp(i,j,kc,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                  compu = rr*qxp(i,j,kc,nq) - compn
                  qxpo(i,j,kc,nq) = compu/rrnew
               end if

               if (i.le.ihi-1) then
                  rr = qxm(i+1,j,kc,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                  compu = rr*qxm(i+1,j,kc,nq) - compn
                  qxmo(i+1,j,kc,nq) = compu/rrnew
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvy(i,j+1,kc)
            pgm = pgdnvy(i,j,kc)
            ugp = ugdnvy(i,j+1,kc)
            ugm = ugdnvy(i,j,kc)

            ! Convert to conservation form
            if (i.ge.ilo+1) then
               rrrx = qxp(i,j,kc,QRHO)
               rurx = rrrx*qxp(i,j,kc,QU)
               rvrx = rrrx*qxp(i,j,kc,QV)
               rwrx = rrrx*qxp(i,j,kc,QW)
               ekenrx = HALF*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
                    + qxp(i,j,kc,QW)**2)
               rerx = qxp(i,j,kc,QREINT) + ekenrx

               ! Add transverse terms
               rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO ) - fy(i,j,kc,URHO ))
               runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX  ) - fy(i,j,kc,UMX  ))
               rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY  ) - fy(i,j,kc,UMY  ))
               rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ  ) - fy(i,j,kc,UMZ  ))
               renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewrx .lt. ZERO) then
                 rrnewrx = rrrx
                 runewrx = rurx
                 rvnewrx = rvrx
                 rwnewrx = rwrx
                 renewrx = rerx
               end if
            end if
   
            if (i.le.ihi-1) then
               rrlx = qxm(i+1,j,kc,QRHO)
               rulx = rrlx*qxm(i+1,j,kc,QU)
               rvlx = rrlx*qxm(i+1,j,kc,QV)
               rwlx = rrlx*qxm(i+1,j,kc,QW)
               ekenlx = HALF*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
                    + qxm(i+1,j,kc,QW)**2)
               relx = qxm(i+1,j,kc,QREINT) + ekenlx

               ! Add transverse terms
               rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO ) - fy(i,j,kc,URHO ))
               runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX  ) - fy(i,j,kc,UMX  ))
               rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY  ) - fy(i,j,kc,UMY  ))
               rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ  ) - fy(i,j,kc,UMZ  ))
               renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewlx .lt. ZERO) then
                    rrnewlx = rrlx
                    runewlx = rulx
                    rvnewlx = rvlx
                    rwnewlx = rwlx
                    renewlx = relx
               end if
            end if

            dup = pgp*ugp - pgm*ugm
            pav = HALF*(pgp+pgm)
            du = ugp-ugm

            ! Convert back to non-conservation form

            ! ************************************************************************
            if (i.ge.ilo+1) then
               qxpo(i,j,kc,QRHO) = rrnewrx
               qxpo(i,j,kc,QU) = runewrx/qxpo(i,j,kc,QRHO)
               qxpo(i,j,kc,QV) = rvnewrx/qxpo(i,j,kc,QRHO)
               qxpo(i,j,kc,QW) = rwnewrx/qxpo(i,j,kc,QRHO)
               rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,QRHO)

               qxpo(i,j,kc,QREINT)= renewrx - rhoekenrx
               qxpo(i,j,kc,QPRES) = qxpo(i,j,kc,QREINT) * gamma_minus_1

               if (qxpo(i,j,kc,QPRES) .lt. small_pres) then
                   pnewrx = qxp(i  ,j,kc,QPRES) - cdtdy*(dup + pav*du*gamma_minus_1)
                   qxpo(i,j,kc,QPRES) = pnewrx
                   qxpo(i,j,kc,QREINT) = qxpo(i,j,kc,QPRES) / gamma_minus_1
               end if
            end if
            ! ************************************************************************

            ! ************************************************************************
            if (i.le.ihi-1) then
               qxmo(i+1,j,kc,QRHO) = rrnewlx
               qxmo(i+1,j,kc,QU) = runewlx/qxmo(i+1,j,kc,QRHO)
               qxmo(i+1,j,kc,QV) = rvnewlx/qxmo(i+1,j,kc,QRHO)
               qxmo(i+1,j,kc,QW) = rwnewlx/qxmo(i+1,j,kc,QRHO)
               rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,QRHO)

               qxmo(i+1,j,kc,QREINT)= renewlx - rhoekenlx
               qxmo(i+1,j,kc,QPRES) = qxmo(i+1,j,kc,QREINT) * gamma_minus_1

               if (qxmo(i+1,j,kc,QPRES) .lt. small_pres) then
                   pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*gamma_minus_1)
                   qxmo(i+1,j,kc,QPRES ) = pnewlx
                   qxmo(i+1,j,kc,QREINT) = qxmo(i+1,j,kc,QPRES) / gamma_minus_1
               end if
            end if
            ! ************************************************************************

         enddo
      enddo

      end subroutine transy1

      !===========================================================================
      ! transy2 -- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transy2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)

      ! Note that cdtdy = dtdy/3.d0/a_half
      double precision cdtdy

      integer i, j
      integer n, nq

      double precision rrnew, rr
      double precision compn, compu
      double precision rrrz, rrlz
      double precision rurz, rulz
      double precision rvrz, rvlz
      double precision rwrz, rwlz
      double precision ekenrz, ekenlz
      double precision rerz, relz
      double precision rrnewrz, rrnewlz
      double precision runewrz, runewlz
      double precision rvnewrz, rvnewlz
      double precision rwnewrz, rwnewlz
      double precision renewrz, renewlz
      double precision pnewrz , pnewlz
      double precision rhoekenrz, rhoekenlz
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            !DIR$ vector always
            do i = ilo, ihi

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qzp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qzp(i,j,kc,nq) - compn
               qzpo(i,j,kc,nq) = compu/rrnew

               compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

               rr = qzm(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
               compu = rr*qzm(i,j,kc,nq) - compn
               qzmo(i,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            ! Convert to conservation form

            ! ************************************************************************
            rrrz = qzp(i,j,kc,QRHO)
            rurz = rrrz*qzp(i,j,kc,QU)
            rvrz = rrrz*qzp(i,j,kc,QV)
            rwrz = rrrz*qzp(i,j,kc,QW)
            ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
                 + qzp(i,j,kc,QW)**2)
            rerz = qzp(i,j,kc,QREINT) + ekenrz

            ! Add transverse terms
            rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO ) - fy(i,j,kc,URHO ))
            runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX  ) - fy(i,j,kc,UMX  ))
            rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY  ) - fy(i,j,kc,UMY  ))
            rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ  ) - fy(i,j,kc,UMZ  ))
            renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewrz .lt. ZERO) then
                 rrnewrz = rrrz
                 runewrz = rurz
                 rvnewrz = rvrz
                 rwnewrz = rwrz
                 renewrz = rerz
            end if

            ! Convert back to non-conservation form
            qzpo(i,j,kc,QRHO) = rrnewrz
            qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
            qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
            qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
            rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)

            qzpo(i,j,kc,QREINT)= renewrz - rhoekenrz
            qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT) * gamma_minus_1

            if (qzpo(i,j,kc,QPRES) .lt. small_pres) then
                pgp = pgdnvy(i,j+1,kc)
                pgm = pgdnvy(i,j,kc)
                ugp = ugdnvy(i,j+1,kc)
                ugm = ugdnvy(i,j,kc)
                dup = pgp*ugp - pgm*ugm
                pav = HALF*(pgp+pgm)
                du = ugp-ugm
                pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*gamma_minus_1)
                qzpo(i,j,kc,QPRES ) = pnewrz
                qzpo(i,j,kc,QREINT) = qzpo(i,j,kc,QPRES) / gamma_minus_1
            end if

            ! ************************************************************************

            ! ************************************************************************
            rrlz = qzm(i,j,kc,QRHO)
            rulz = rrlz*qzm(i,j,kc,QU)
            rvlz = rrlz*qzm(i,j,kc,QV)
            rwlz = rrlz*qzm(i,j,kc,QW)
            ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
                 + qzm(i,j,kc,QW)**2)
            relz = qzm(i,j,kc,QREINT) + ekenlz

            ! Add transverse terms
            rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO ) - fy(i,j,km,URHO ))
            runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX  ) - fy(i,j,km,UMX  ))
            rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY  ) - fy(i,j,km,UMY  ))
            rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ  ) - fy(i,j,km,UMZ  ))
            renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN) - fy(i,j,km,UEDEN))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewlz .lt. ZERO) then
                 rrnewlz = rrlz
                 runewlz = rulz
                 rvnewlz = rvlz
                 rwnewlz = rwlz
                 renewlz = relz
            end if

            qzmo(i,j,kc,QRHO) = rrnewlz
            qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
            qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
            qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
            rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)

            qzmo(i,j,kc,QREINT)= renewlz - rhoekenlz
            qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT) * gamma_minus_1

            if (qzmo(i,j,kc,QPRES) .lt. small_pres) then
                pgp = pgdnvy(i,j+1,km)
                pgm = pgdnvy(i,j,km)
                ugp = ugdnvy(i,j+1,km)
                ugm = ugdnvy(i,j,km)
                dup = pgp*ugp - pgm*ugm
                pav = HALF*(pgp+pgm)
                du = ugp-ugm
                pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*gamma_minus_1)
                qzmo(i,j,kc,QPRES ) = pnewlz
                qzmo(i,j,kc,QREINT) = qzmo(i,j,kc,QPRES) / gamma_minus_1
            end if
            ! ************************************************************************

         enddo
      enddo

      end subroutine transy2

      !===========================================================================
      ! transz -- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transz(qxm,qxmo,qxp,qxpo, &
                        qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        fz,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                        ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                        cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     npassive, upass_map, qpass_map, & 
                                     small_pres, gamma_minus_1
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)

      ! Note that cdtdz = dtdz/3.d0/a_half
      double precision cdtdz

      integer n, nq
      integer i, j

      double precision rrnew, rr
      double precision compn, compu
      double precision, dimension(ilo:ihi,jlo:jhi) :: rrrx, rrry, rrlx, rrly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rurx, rury, rulx, ruly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rvrx, rvry, rvlx, rvly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rwrx, rwry, rwlx, rwly
      double precision, dimension(ilo:ihi,jlo:jhi) :: ekenrx, ekenry, ekenlx, ekenly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rerx, rery, relx, rely
      double precision, dimension(ilo:ihi,jlo:jhi) :: rrnewrx, rrnewry, rrnewlx, rrnewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: runewrx, runewry, runewlx, runewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rvnewrx, rvnewry, rvnewlx, rvnewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rwnewrx, rwnewry, rwnewlx, rwnewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: renewrx, renewry, renewlx, renewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: pnewrx,  pnewry,  pnewlx,  pnewly
      double precision, dimension(ilo:ihi,jlo:jhi) :: rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
      double precision, dimension(ilo:ihi,jlo:jhi) :: pgp, pgm, ugp, ugm, dup, pav, du

      integer          :: ipassive

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

                 compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

                 if (i.ge.ilo+1) then
                    rr = qxp(i,j,km,QRHO)
                    rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                    compu = rr*qxp(i,j,km,nq) - compn
                    qxpo(i,j,km,nq) = compu/rrnew
                 end if

                 if (j.ge.jlo+1) then
                    rr = qyp(i,j,km,QRHO)
                    rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                    compu = rr*qyp(i,j,km,nq) - compn
                    qypo(i,j,km,nq) = compu/rrnew
                 end if

                 if (i.le.ihi-1) then
                    rr = qxm(i+1,j,km,QRHO)
                    rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                    compu = rr*qxm(i+1,j,km,nq) - compn
                    qxmo(i+1,j,km,nq) = compu/rrnew
                 end if

                 if (j.le.jhi-1) then
                    rr = qym(i,j+1,km,QRHO)
                    rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                    compu = rr*qym(i,j+1,km,nq) - compn
                    qymo(i,j+1,km,nq) = compu/rrnew
                 end if

            enddo
         enddo
      enddo

      pgp(ilo:ihi,jlo:jhi) = pgdnvz(ilo:ihi,jlo:jhi,kc)
      pgm(ilo:ihi,jlo:jhi) = pgdnvz(ilo:ihi,jlo:jhi,km)
      ugp(ilo:ihi,jlo:jhi) = ugdnvz(ilo:ihi,jlo:jhi,kc)
      ugm(ilo:ihi,jlo:jhi) = ugdnvz(ilo:ihi,jlo:jhi,km)

      do i = ilo, ihi

          if (i.ge.ilo+1) then
            ! Convert to conservation form
            rrrx(i,jlo:jhi) = qxp(i,jlo:jhi,km,QRHO)
            rurx(i,jlo:jhi) = rrrx(i,jlo:jhi)*qxp(i,jlo:jhi,km,QU)
            rvrx(i,jlo:jhi) = rrrx(i,jlo:jhi)*qxp(i,jlo:jhi,km,QV)
            rwrx(i,jlo:jhi) = rrrx(i,jlo:jhi)*qxp(i,jlo:jhi,km,QW)
            ekenrx(i,jlo:jhi) = HALF*rrrx(i,jlo:jhi)*(qxp(i,jlo:jhi,km,QU)**2 + qxp(i,jlo:jhi,km,QV)**2 &
                 + qxp(i,jlo:jhi,km,QW)**2)
            rerx(i,jlo:jhi) = qxp(i,jlo:jhi,km,QREINT) + ekenrx(i,jlo:jhi)

            ! Add transverse terms
            rrnewrx(i,jlo:jhi) = rrrx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,URHO ) - fz(i,jlo:jhi,km,URHO ))
            runewrx(i,jlo:jhi) = rurx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMX  ) - fz(i,jlo:jhi,km,UMX  ))
            rvnewrx(i,jlo:jhi) = rvrx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMY  ) - fz(i,jlo:jhi,km,UMY  ))
            rwnewrx(i,jlo:jhi) = rwrx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMZ  ) - fz(i,jlo:jhi,km,UMZ  ))
            renewrx(i,jlo:jhi) = rerx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UEDEN) - fz(i,jlo:jhi,km,UEDEN))
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (i.ge.ilo+1) then
                ! Reset to original value if adding transverse terms made density negative
                if (rrnewrx(i,j) .lt. ZERO) then
                   rrnewrx(i,j) = rrrx(i,j)
                   runewrx(i,j) = rurx(i,j)
                   rvnewrx(i,j) = rvrx(i,j)
                   rwnewrx(i,j) = rwrx(i,j)
                   renewrx(i,j) = rerx(i,j)
                end if
              end if

          enddo
      enddo

      do j = jlo, jhi

          if (j.ge.jlo+1) then
            rrry(ilo:ihi,j) = qyp(ilo:ihi,j,km,QRHO)
            rury(ilo:ihi,j) = rrry(ilo:ihi,j)*qyp(ilo:ihi,j,km,QU)
            rvry(ilo:ihi,j) = rrry(ilo:ihi,j)*qyp(ilo:ihi,j,km,QV)
            rwry(ilo:ihi,j) = rrry(ilo:ihi,j)*qyp(ilo:ihi,j,km,QW)
            ekenry(ilo:ihi,j) = HALF*rrry(ilo:ihi,j)*(qyp(ilo:ihi,j,km,QU)**2 + qyp(ilo:ihi,j,km,QV)**2 &
                 + qyp(ilo:ihi,j,km,QW)**2)
            rery(ilo:ihi,j) = qyp(ilo:ihi,j,km,QREINT) + ekenry(ilo:ihi,j)

            ! Add transverse terms
            rrnewry(ilo:ihi,j) = rrry(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,URHO ) - fz(ilo:ihi,j,km,URHO ))
            runewry(ilo:ihi,j) = rury(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMX  ) - fz(ilo:ihi,j,km,UMX  ))
            rvnewry(ilo:ihi,j) = rvry(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMY  ) - fz(ilo:ihi,j,km,UMY  ))
            rwnewry(ilo:ihi,j) = rwry(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMZ  ) - fz(ilo:ihi,j,km,UMZ  ))
            renewry(ilo:ihi,j) = rery(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UEDEN) - fz(ilo:ihi,j,km,UEDEN))
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.ge.jlo+1) then
                ! Reset to original value if adding transverse terms made density negative
                if (rrnewry(i,j) .lt. ZERO) then
                   rrnewry(i,j) = rrry(i,j)
                   runewry(i,j) = rury(i,j)
                   rvnewry(i,j) = rvry(i,j)
                   rwnewry(i,j) = rwry(i,j)
                   renewry(i,j) = rery(i,j)
                end if
              end if

          enddo
      enddo

      do i = ilo, ihi

          if (i.le.ihi-1) then
            rrlx(i,jlo:jhi) = qxm(i+1,jlo:jhi,km,QRHO)
            rulx(i,jlo:jhi) = rrlx(i,jlo:jhi)*qxm(i+1,jlo:jhi,km,QU)
            rvlx(i,jlo:jhi) = rrlx(i,jlo:jhi)*qxm(i+1,jlo:jhi,km,QV)
            rwlx(i,jlo:jhi) = rrlx(i,jlo:jhi)*qxm(i+1,jlo:jhi,km,QW)
            ekenlx(i,jlo:jhi) = HALF*rrlx(i,jlo:jhi)*(qxm(i+1,jlo:jhi,km,QU)**2 + qxm(i+1,jlo:jhi,km,QV)**2 &
                 + qxm(i+1,jlo:jhi,km,QW)**2)
            relx(i,jlo:jhi) = qxm(i+1,jlo:jhi,km,QREINT) + ekenlx(i,jlo:jhi)

            ! Add transverse terms
            rrnewlx(i,jlo:jhi) = rrlx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,URHO ) - fz(i,jlo:jhi,km,URHO ))
            runewlx(i,jlo:jhi) = rulx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMX  ) - fz(i,jlo:jhi,km,UMX  ))
            rvnewlx(i,jlo:jhi) = rvlx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMY  ) - fz(i,jlo:jhi,km,UMY  ))
            rwnewlx(i,jlo:jhi) = rwlx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UMZ  ) - fz(i,jlo:jhi,km,UMZ  ))
            renewlx(i,jlo:jhi) = relx(i,jlo:jhi) - cdtdz*(fz(i,jlo:jhi,kc,UEDEN) - fz(i,jlo:jhi,km,UEDEN))
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (i.le.ihi-1) then
                ! Reset to original value if adding transverse terms made density negative
                if (rrnewlx(i,j) .lt. ZERO) then
                   rrnewlx(i,j) = rrlx(i,j)
                   runewlx(i,j) = rulx(i,j)
                   rvnewlx(i,j) = rvlx(i,j)
                   rwnewlx(i,j) = rwlx(i,j)
                   renewlx(i,j) = relx(i,j)
                end if
              end if

          enddo
      enddo

      do j = jlo, jhi

          if (j.le.jhi-1) then
            rrly(ilo:ihi,j) = qym(ilo:ihi,j+1,km,QRHO)
            ruly(ilo:ihi,j) = rrly(ilo:ihi,j)*qym(ilo:ihi,j+1,km,QU)
            rvly(ilo:ihi,j) = rrly(ilo:ihi,j)*qym(ilo:ihi,j+1,km,QV)
            rwly(ilo:ihi,j) = rrly(ilo:ihi,j)*qym(ilo:ihi,j+1,km,QW)
            ekenly(ilo:ihi,j) = HALF*rrly(ilo:ihi,j)*(qym(ilo:ihi,j+1,km,QU)**2 + qym(ilo:ihi,j+1,km,QV)**2 &
                 + qym(ilo:ihi,j+1,km,QW)**2)
            rely(ilo:ihi,j) = qym(ilo:ihi,j+1,km,QREINT) + ekenly(ilo:ihi,j)

            ! Add transverse terms
            rrnewly(ilo:ihi,j) = rrly(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,URHO ) - fz(ilo:ihi,j,km,URHO )) 
            runewly(ilo:ihi,j) = ruly(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMX  ) - fz(ilo:ihi,j,km,UMX  )) 
            rvnewly(ilo:ihi,j) = rvly(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMY  ) - fz(ilo:ihi,j,km,UMY  ))
            rwnewly(ilo:ihi,j) = rwly(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UMZ  ) - fz(ilo:ihi,j,km,UMZ  ))
            renewly(ilo:ihi,j) = rely(ilo:ihi,j) - cdtdz*(fz(ilo:ihi,j,kc,UEDEN) - fz(ilo:ihi,j,km,UEDEN))
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.le.jhi-1) then
                ! Reset to original value if adding transverse terms made density negative
                if (rrnewly(i,j) .lt. ZERO) then
                   rrnewly(i,j) = rrly(i,j)
                   runewly(i,j) = ruly(i,j)
                   rvnewly(i,j) = rvly(i,j)
                   rwnewly(i,j) = rwly(i,j)
                   renewly(i,j) = rely(i,j)
                end if
              end if

          enddo
      enddo

      dup(ilo:ihi,jlo:jhi) = pgp(ilo:ihi,jlo:jhi)*ugp(ilo:ihi,jlo:jhi) - pgm(ilo:ihi,jlo:jhi)*ugm(ilo:ihi,jlo:jhi)
      pav(ilo:ihi,jlo:jhi) = HALF*(pgp(ilo:ihi,jlo:jhi)+pgm(ilo:ihi,jlo:jhi))
      du(ilo:ihi,jlo:jhi) = ugp(ilo:ihi,jlo:jhi)-ugm(ilo:ihi,jlo:jhi)


      do i = ilo, ihi

          ! Convert back to non-conservation form
          if (i.ge.ilo+1) then
             qxpo(i,jlo:jhi,km,QRHO) = rrnewrx(i,jlo:jhi)
             qxpo(i,jlo:jhi,km,QU) = runewrx(i,jlo:jhi)/qxpo(i,jlo:jhi,km,QRHO)
             qxpo(i,jlo:jhi,km,QV) = rvnewrx(i,jlo:jhi)/qxpo(i,jlo:jhi,km,QRHO)
             qxpo(i,jlo:jhi,km,QW) = rwnewrx(i,jlo:jhi)/qxpo(i,jlo:jhi,km,QRHO)
             rhoekenrx(i,jlo:jhi) = HALF*(runewrx(i,jlo:jhi)**2 + rvnewrx(i,jlo:jhi)**2 + rwnewrx(i,jlo:jhi)**2) / &
                                    qxpo(i,jlo:jhi,km,QRHO)

             qxpo(i,jlo:jhi,km,QREINT)= renewrx(i,jlo:jhi) - rhoekenrx(i,jlo:jhi)
             qxpo(i,jlo:jhi,km,QPRES) = qxpo(i,jlo:jhi,km,QREINT) * gamma_minus_1
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (i.ge.ilo+1) then
                 if (qxpo(i,j,km,QPRES) .lt. small_pres) then
                     pnewrx = qxp(i  ,j,km,QPRES) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qxpo(i,j,km,QPRES ) = pnewrx(i,j)
                     qxpo(i,j,km,QREINT) = qxpo(i,j,km,QPRES) / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      do j = jlo, jhi

          if (j.ge.jlo+1) then
             qypo(ilo:ihi,j,km,QRHO) = rrnewry(ilo:ihi,j)
             qypo(ilo:ihi,j,km,QU) = runewry(ilo:ihi,j)/qypo(ilo:ihi,j,km,QRHO)
             qypo(ilo:ihi,j,km,QV) = rvnewry(ilo:ihi,j)/qypo(ilo:ihi,j,km,QRHO)
             qypo(ilo:ihi,j,km,QW) = rwnewry(ilo:ihi,j)/qypo(ilo:ihi,j,km,QRHO)
             rhoekenry(ilo:ihi,j) = HALF*(runewry(ilo:ihi,j)**2 + rvnewry(ilo:ihi,j)**2 + rwnewry(ilo:ihi,j)**2) / &
                                    qypo(ilo:ihi,j,km,QRHO)

             qypo(ilo:ihi,j,km,QREINT)= renewry(ilo:ihi,j) - rhoekenry(ilo:ihi,j)
             qypo(ilo:ihi,j,km,QPRES) = qypo(ilo:ihi,j,km,QREINT) * gamma_minus_1
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.ge.jlo+1) then
                 if (qypo(i,j,km,QPRES) .lt. small_pres) then
                     pnewry(i,j) = qyp(i,j  ,km,QPRES) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qypo(i,j,km,QPRES ) = pnewry(i,j)
                     qypo(i,j,km,QREINT) = qypo(i,j,km,QPRES) / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      do i = ilo, ihi

          if (i.le.ihi-1) then
             qxmo(i+1,jlo:jhi,km,QRHO) = rrnewlx(i,jlo:jhi)
             qxmo(i+1,jlo:jhi,km,QU) = runewlx(i,jlo:jhi)/qxmo(i+1,jlo:jhi,km,QRHO)
             qxmo(i+1,jlo:jhi,km,QV) = rvnewlx(i,jlo:jhi)/qxmo(i+1,jlo:jhi,km,QRHO)
             qxmo(i+1,jlo:jhi,km,QW) = rwnewlx(i,jlo:jhi)/qxmo(i+1,jlo:jhi,km,QRHO)
             rhoekenlx(i,jlo:jhi) = HALF*(runewlx(i,jlo:jhi)**2 + rvnewlx(i,jlo:jhi)**2 + rwnewlx(i,jlo:jhi)**2) / &
                                    qxmo(i+1,jlo:jhi,km,QRHO)

             qxmo(i+1,jlo:jhi,km,QREINT)= renewlx(i,jlo:jhi) - rhoekenlx(i,jlo:jhi)
             qxmo(i+1,jlo:jhi,km,QPRES) = qxmo(i+1,jlo:jhi,km,QREINT) * gamma_minus_1
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (i.le.ihi-1) then
                 if (qxmo(i+1,j,km,QPRES) .lt. small_pres) then
                     pnewlx(i,j) = qxm(i+1,j,km,QPRES) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qxmo(i+1,j,km,QPRES ) = pnewlx(i,j)
                     qxmo(i+1,j,km,QREINT) = qxmo(i+1,j,km,QPRES)  / gamma_minus_1
                 end if
              end if

          enddo
      enddo

      do j = jlo, jhi

          if (j.le.jhi-1) then
             qymo(ilo:ihi,j+1,km,QRHO) = rrnewly(ilo:ihi,j)
             qymo(ilo:ihi,j+1,km,QU) = runewly(ilo:ihi,j)/qymo(ilo:ihi,j+1,km,QRHO)
             qymo(ilo:ihi,j+1,km,QV) = rvnewly(ilo:ihi,j)/qymo(ilo:ihi,j+1,km,QRHO)
             qymo(ilo:ihi,j+1,km,QW) = rwnewly(ilo:ihi,j)/qymo(ilo:ihi,j+1,km,QRHO)
             rhoekenly(ilo:ihi,j) = HALF*(runewly(ilo:ihi,j)**2 + rvnewly(ilo:ihi,j)**2 + rwnewly(ilo:ihi,j)**2) / &
                                    qymo(ilo:ihi,j+1,km,QRHO)

             qymo(ilo:ihi,j+1,km,QREINT)= renewly(ilo:ihi,j) - rhoekenly(ilo:ihi,j)
             qymo(ilo:ihi,j+1,km,QPRES) = qymo(ilo:ihi,j+1,km,QREINT) * gamma_minus_1
          end if

      enddo

      do j = jlo, jhi
          do i = ilo, ihi

              if (j.le.jhi-1) then
                 if (qymo(i,j+1,km,QPRES) .lt. small_pres) then
                     pnewly(i,j) = qym(i,j+1,km,QPRES) - cdtdz*(dup(i,j) + pav(i,j)*du(i,j)*gamma_minus_1)
                     qymo(i,j+1,km,QPRES ) = pnewly(i,j)
                     qymo(i,j+1,km,QREINT) = qymo(i,j+1,km,QPRES) / gamma_minus_1
                 endif
              end if

          enddo
      enddo

      end subroutine transz

      !===========================================================================
      ! transxy-- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transxy(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxy,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fyx,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         hdt,hdtdx,hdtdy,ilo,ihi,jlo,jhi,kc,km,k3d,a_old,a_new)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt,a_old,a_new

      ! Note that hdtdx = dtdx/2.d0/a_half
      !       and hdtdy = dtdy/2.d0/a_half
      double precision hdtdx,hdtdy

      integer i, j
      integer n , nq

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      double precision pgxp, pgxm, ugxp, ugxm
      double precision pgyp, pgym, ugyp, ugym
      double precision pgxpm, pgxmm, ugxpm, ugxmm
      double precision pgypm, pgymm, ugypm, ugymm
      double precision compr, compl, compnr, compnl
      double precision a_half

      double precision :: dux,duxm,duxp,duxpm,duy,duym,duyp,duypm
      double precision :: pxav,pxavm,pxnew,pxnewm,pyav,pyavm,pynew,pynewm
      integer          :: ipassive

      a_half = HALF * (a_old + a_new)

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               rrr = qp(i,j,kc,QRHO)
               rrl = qm(i,j,kc,QRHO)

               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)

               rrnewr = rrr - hdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                            - hdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
               rrnewl = rrl - hdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                            - hdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

               compnr = compr - hdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                              - hdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - hdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                              - hdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d  ,nq) / a_half
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq) / a_half

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

             pgxp = pgdnvx(i+1,j,kc)
             pgxm = pgdnvx(i,j,kc)
             ugxp = ugdnvx(i+1,j,kc)
             ugxm = ugdnvx(i,j,kc)
 
             pgyp = pgdnvy(i,j+1,kc)
             pgym = pgdnvy(i,j,kc)
             ugyp = ugdnvy(i,j+1,kc)
             ugym = ugdnvy(i,j,kc)
 
             pgxpm = pgdnvx(i+1,j,km)
             pgxmm = pgdnvx(i,j,km)
             ugxpm = ugdnvx(i+1,j,km)
             ugxmm = ugdnvx(i,j,km)
 
             pgypm = pgdnvy(i,j+1,km)
             pgymm = pgdnvy(i,j,km)
             ugypm = ugdnvy(i,j+1,km)
             ugymm = ugdnvy(i,j,km)
 
             ! Convert to conservation form
             rrr = qp(i,j,kc,QRHO)
             rur = rrr*qp(i,j,kc,QU)
             rvr = rrr*qp(i,j,kc,QV)
             rwr = rrr*qp(i,j,kc,QW)
             ekenr = HALF*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
                  qp(i,j,kc,QW)**2)
             rer = qp(i,j,kc,QREINT) + ekenr
 
             rrl = qm(i,j,kc,QRHO)
             rul = rrl*qm(i,j,kc,QU)
             rvl = rrl*qm(i,j,kc,QV)
             rwl = rrl*qm(i,j,kc,QW)
             ekenl = HALF*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
                  qm(i,j,kc,QW)**2)
             rel = qm(i,j,kc,QREINT) + ekenl

             rrnewr = rrr - hdtdx*(fxy(i+1,j,kc,URHO ) - fxy(i,j,kc,URHO )) &
                          - hdtdy*(fyx(i,j+1,kc,URHO ) - fyx(i,j,kc,URHO ))
             runewr = rur - hdtdx*(fxy(i+1,j,kc,UMX  ) - fxy(i,j,kc,UMX  )) &
                          - hdtdy*(fyx(i,j+1,kc,UMX  ) - fyx(i,j,kc,UMX  ))
             rvnewr = rvr - hdtdx*(fxy(i+1,j,kc,UMY  ) - fxy(i,j,kc,UMY  )) &
                          - hdtdy*(fyx(i,j+1,kc,UMY  ) - fyx(i,j,kc,UMY  ))
             rwnewr = rwr - hdtdx*(fxy(i+1,j,kc,UMZ  ) - fxy(i,j,kc,UMZ  )) &
                          - hdtdy*(fyx(i,j+1,kc,UMZ  ) - fyx(i,j,kc,UMZ  ))
             renewr = rer - hdtdx*(fxy(i+1,j,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
                          - hdtdy*(fyx(i,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))

             rrnewl = rrl - hdtdx*(fxy(i+1,j,km,URHO ) - fxy(i,j,km,URHO )) &
                          - hdtdy*(fyx(i,j+1,km,URHO ) - fyx(i,j,km,URHO ))
             runewl = rul - hdtdx*(fxy(i+1,j,km,UMX  ) - fxy(i,j,km,UMX  )) &
                          - hdtdy*(fyx(i,j+1,km,UMX  ) - fyx(i,j,km,UMX  ))
             rvnewl = rvl - hdtdx*(fxy(i+1,j,km,UMY  ) - fxy(i,j,km,UMY  )) &
                          - hdtdy*(fyx(i,j+1,km,UMY  ) - fyx(i,j,km,UMY  ))
             rwnewl = rwl - hdtdx*(fxy(i+1,j,km,UMZ  ) - fxy(i,j,km,UMZ  )) &
                          - hdtdy*(fyx(i,j+1,km,UMZ  ) - fyx(i,j,km,UMZ  ))
             renewl = rel - hdtdx*(fxy(i+1,j,km,UEDEN) - fxy(i,j,km,UEDEN)) &
                          - hdtdy*(fyx(i,j+1,km,UEDEN) - fyx(i,j,km,UEDEN))

            ! Reset to original value if adding transverse terms made density negative
            if (rrnewr .lt. ZERO) then
                 rrnewr = rrr
                 runewr = rur
                 rvnewr = rvr
                 rwnewr = rwr
                 renewr = rer
            end if
            if (rrnewl .lt. ZERO) then
                 rrnewl = rrl
                 runewl = rul
                 rvnewl = rvl
                 rwnewl = rwl
                 renewl = rel
            end if

            rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
            rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            qpo(i,j,kc,QRHO  ) = rrnewr

            qpo(i,j,kc,QU    ) = runewr/rrnewr
            qpo(i,j,kc,QV    ) = rvnewr/rrnewr
            qpo(i,j,kc,QW    ) = rwnewr/rrnewr

            qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt* srcQ(i,j,k3d,QREINT) / a_old
            qpo(i,j,kc,QPRES) = qpo(i,j,kc,QREINT) * gamma_minus_1

            if (qpo(i,j,kc,QPRES) .lt. small_pres) then
                duxp = pgxp*ugxp - pgxm*ugxm
                pxav = HALF*(pgxp+pgxm)
                dux = ugxp-ugxm
                pxnew = hdtdx*(duxp + pxav*dux*gamma_minus_1)

                duyp = pgyp*ugyp - pgym*ugym
                pyav = HALF*(pgyp+pgym)
                duy = ugyp-ugym
                pynew = hdtdy*(duyp + pyav*duy*gamma_minus_1)

                pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
                qpo(i,j,kc,QPRES ) = pnewr + hdt* srcQ(i,j,k3d,QPRES ) / a_old
                qpo(i,j,kc,QREINT) = qpo(i,j,kc,QPRES) / gamma_minus_1
            end if

            ! ************************************************************************* 
            qmo(i,j,kc,QRHO  ) = rrnewl

            qmo(i,j,kc,QU    ) = runewl/rrnewl
            qmo(i,j,kc,QV    ) = rvnewl/rrnewl
            qmo(i,j,kc,QW    ) = rwnewl/rrnewl

            qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt* srcQ(i,j,k3d-1,QREINT) / a_old
            qmo(i,j,kc,QPRES) = qmo(i,j,kc,QREINT) * gamma_minus_1

            if (qmo(i,j,kc,QPRES) .lt. small_pres) then

                duxpm = pgxpm*ugxpm - pgxmm*ugxmm
                pxavm = HALF*(pgxpm+pgxmm)
                duxm = ugxpm-ugxmm
                pxnewm = hdtdx*(duxpm + pxavm*duxm*gamma_minus_1)

                duypm = pgypm*ugypm - pgymm*ugymm
                pyavm = HALF*(pgypm+pgymm)
                duym = ugypm-ugymm
                pynewm = hdtdy*(duypm + pyavm*duym*gamma_minus_1)

                pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm
                qmo(i,j,kc,QPRES ) = pnewl + hdt* srcQ(i,j,k3d-1,QPRES ) / a_old
                qmo(i,j,kc,QREINT) = qmo(i,j,kc,QPRES) / gamma_minus_1
            end if
            ! ************************************************************************* 
         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               qpo(i,j,kc,QRHO) = qpo(i,j,kc,QRHO) + hdt*srcQ(i,j,k3d,QRHO) / a_old
               qpo(i,j,kc,QU  ) = qpo(i,j,kc,QU  ) + hdt*srcQ(i,j,k3d,QU  ) / a_old
               qpo(i,j,kc,QV  ) = qpo(i,j,kc,QV  ) + hdt*srcQ(i,j,k3d,QV  ) / a_old
               qpo(i,j,kc,QW  ) = qpo(i,j,kc,QW  ) + hdt*srcQ(i,j,k3d,QW  ) / a_old

               qmo(i,j,kc,QRHO) = qmo(i,j,kc,QRHO) + hdt*srcQ(i,j,k3d-1,QRHO) / a_old
               qmo(i,j,kc,QU  ) = qmo(i,j,kc,QU  ) + hdt*srcQ(i,j,k3d-1,QU  ) / a_old
               qmo(i,j,kc,QV  ) = qmo(i,j,kc,QV  ) + hdt*srcQ(i,j,k3d-1,QV  ) / a_old 
               qmo(i,j,kc,QW  ) = qmo(i,j,kc,QW  ) + hdt*srcQ(i,j,k3d-1,QW  ) / a_old
            enddo
         enddo
      endif

      end subroutine transxy

      !===========================================================================
      ! transz-- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transxz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxz,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fzx,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         hdt,hdtdx,hdtdz,ilo,ihi,jlo,jhi,km,kc,k3d,a_old,a_new)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt,a_old,a_new

      ! Note that hdtdx = dtdx/2.d0/a_half
      !       and hdtdz = dtdz/2.d0/a_half
      double precision hdtdx,hdtdz

      integer i, j
      integer n, nq

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      double precision pgxp, pgxm, ugxp, ugxm
      double precision pgzp, pgzm, ugzp, ugzm
      double precision compr, compl, compnr, compnl
      double precision a_half

      double precision :: drr, dcompn
      double precision :: dux, duxp, duz, duzp, pxav, pxnew, pzav, pznew
      integer          :: ipassive

      a_half = HALF * (a_old + a_new)
 
      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               drr    = - hdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                        - hdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
               dcompn = - hdtdx*(fxz(i+1,j,km,n   ) - fxz(i,j,km,n)) &
                        - hdtdz*(fzx(i  ,j,kc,n   ) - fzx(i,j,km,n))

               if (j.ge.jlo+1) then
                  rrr = qp(i,j,km,QRHO)
                  compr = rrr*qp(i,j,km,nq)

                  rrnewr = rrr   + drr
                  compnr = compr + dcompn

                  qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq) / a_half
               end if

               if (j.le.jhi-1) then
                  rrl = qm(i,j+1,km,QRHO)
                  compl = rrl*qm(i,j+1,km,nq)

                  rrnewl = rrl   + drr
                  compnl = compl + dcompn

                  qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq) / a_half
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgxp = pgdnvx(i+1,j,km)
            pgxm = pgdnvx(i,j,km)
            ugxp = ugdnvx(i+1,j,km)
            ugxm = ugdnvx(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            if (j.ge.jlo+1) then
               ! Convert to conservation form
               rrr = qp(i,j,km,QRHO)
               rur = rrr*qp(i,j,km,QU)
               rvr = rrr*qp(i,j,km,QV)
               rwr = rrr*qp(i,j,km,QW)
               ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
               rer = qp(i,j,km,QREINT) + ekenr

               ! Add transverse terms
               rrnewr = rrr - hdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO )) &
                            - hdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO ))
               runewr = rur - hdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX  )) &
                            - hdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX  ))
               rvnewr = rvr - hdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY  )) &
                            - hdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY  ))
               rwnewr = rwr - hdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ  )) &
                            - hdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ  ))
               renewr = rer - hdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                            - hdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewr .lt. ZERO) then
                    rrnewr = rrr
                    runewr = rur
                    rvnewr = rvr
                    rwnewr = rwr
                    renewr = rer
               end if

               rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            end if

            if (j.le.jhi-1) then
               rrl = qm(i,j+1,km,QRHO)
               rul = rrl*qm(i,j+1,km,QU)
               rvl = rrl*qm(i,j+1,km,QV)
               rwl = rrl*qm(i,j+1,km,QW)
               ekenl = HALF*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
               rel = qm(i,j+1,km,QREINT) + ekenl

              ! Add transverse terms
              rrnewl = rrl - hdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO )) &
                           - hdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO ))
              runewl = rul - hdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX  )) &
                           - hdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX  ))
              rvnewl = rvl - hdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY  )) &
                           - hdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY  ))
              rwnewl = rwl - hdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ  )) &
                           - hdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ  ))
              renewl = rel - hdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                           - hdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN))
               ! Reset to original value if adding transverse terms made density negative
               if (rrnewl .lt. ZERO) then
                    rrnewl = rrl
                    runewl = rul
                    rvnewl = rvl
                    rwnewl = rwl
                    renewl = rel
               end if

               rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            end if

            ! These are used only if the (rho e) version goes negative 
            duxp = pgxp*ugxp - pgxm*ugxm
            pxav = HALF*(pgxp+pgxm)
            dux = ugxp-ugxm
            pxnew = hdtdx*(duxp + pxav*dux*gamma_minus_1)

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = HALF*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = hdtdz*(duzp + pzav*duz*gamma_minus_1)
            ! End comment

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            if (j.ge.jlo+1) then
               qpo(i,j,km,QRHO  ) = rrnewr
               qpo(i,j,km,QU    ) = runewr/rrnewr
               qpo(i,j,km,QV    ) = rvnewr/rrnewr
               qpo(i,j,km,QW    ) = rwnewr/rrnewr

               qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt* srcQ(i,j,k3d,QREINT) / a_old
               qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT) * gamma_minus_1

               if (qpo(i,j,km,QPRES) .lt. small_pres) then
                   pnewr = qp(i,j  ,km,QPRES) - pxnew - pznew
                   qpo(i,j,km,QPRES ) = pnewr + hdt* srcQ(i,j,k3d,QPRES ) / a_old
                   qpo(i,j,km,QREINT) = qpo(i,j,km,QPRES) / gamma_minus_1
               end if
            end if

            ! ************************************************************************* 
            if (j.le.jhi-1) then
               qmo(i,j+1,km,QRHO  ) = rrnewl
               qmo(i,j+1,km,QU    ) = runewl/rrnewl
               qmo(i,j+1,km,QV    ) = rvnewl/rrnewl
               qmo(i,j+1,km,QW    ) = rwnewl/rrnewl

               qmo(i,j+1,km,QREINT) = renewl - rhoekenl + hdt* srcQ(i,j,k3d,QREINT) / a_old
               qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QREINT) * gamma_minus_1

               if (qmo(i,j+1,km,QPRES) .lt. small_pres) then
                   pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew
                   qmo(i,j+1,km,QPRES ) = pnewl + hdt* srcQ(i,j,k3d,QPRES ) / a_old
                   qmo(i,j+1,km,QREINT) = qmo(i,j+1,km,QPRES) / gamma_minus_1
               end if

            end if
            ! ************************************************************************* 

         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               if (j.ge.jlo+1) then
                  qpo(i,j,km,QRHO) = qpo(i,j,km,QRHO) + hdt*srcQ(i,j,k3d,QRHO  ) / a_old
                  qpo(i,j,km,QU  ) = qpo(i,j,km,QU  ) + hdt*srcQ(i,j,k3d,QU) / a_old
                  qpo(i,j,km,QV  ) = qpo(i,j,km,QV  ) + hdt*srcQ(i,j,k3d,QV) / a_old
                  qpo(i,j,km,QW  ) = qpo(i,j,km,QW  ) + hdt*srcQ(i,j,k3d,QW) / a_old
               end if

               if (j.le.jhi-1) then
                  qmo(i,j+1,km,QRHO) = qmo(i,j+1,km,QRHO) + hdt*srcQ(i,j,k3d,QRHO  ) / a_old
                  qmo(i,j+1,km,QU  ) = qmo(i,j+1,km,QU  ) + hdt*srcQ(i,j,k3d,QU) / a_old
                  qmo(i,j+1,km,QV  ) = qmo(i,j+1,km,QV  ) + hdt*srcQ(i,j,k3d,QV) / a_old
                  qmo(i,j+1,km,QW  ) = qmo(i,j+1,km,QW  ) + hdt*srcQ(i,j,k3d,QW) / a_old
               end if
            enddo
         enddo
      endif

      end subroutine transxz

      !===========================================================================
      ! transyz-- called from within threaded loops in advance_gas_tile so *no* OMP here ...
      !===========================================================================
      subroutine transyz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fyz,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         fzy,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         hdt,hdtdy,hdtdz,ilo,ihi,jlo,jhi,km,kc,k3d,a_old,a_new)

      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, &
                                     URHO, UMX, UMY, UMZ, UEDEN, &
                                     small_pres, gamma_minus_1, &
                                     npassive, upass_map, qpass_map, & 
                                     version_2
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision fzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt,a_old,a_new

      ! Note that hdtdy = dtdy/2.d0/a_half
      !       and hdtdz = dtdz/2.d0/a_half
      double precision hdtdy,hdtdz

      integer i, j
      integer n, nq

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr, pnewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl, pnewl
      double precision pgyp, pgym, ugyp, ugym
      double precision pgzp, pgzm, ugzp, ugzm
      double precision compr, compl, compnr, compnl
      double precision a_half

      double precision :: drr, dcompn
      double precision :: duy,duyp,duz,duzp,pyav,pynew,pzav,pznew
      integer          :: ipassive

      a_half = HALF * (a_old + a_new)

      do ipassive = 1,npassive
         n  = upass_map(ipassive)
         nq = qpass_map(ipassive)
         do j = jlo, jhi
            do i = ilo, ihi

               drr    = - hdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                        - hdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
               dcompn = - hdtdy*(fyz(i,j+1,km,n   ) - fyz(i,j,km,n)) &
                        - hdtdz*(fzy(i,j  ,kc,n   ) - fzy(i,j,km,n))

               if (i.ge.ilo+1) then
                  rrr = qp(i,j,km,QRHO)
                  compr  = rrr*qp(i,j,km,nq)
   
                  rrnewr = rrr   + drr
                  compnr = compr + dcompn

                  qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq) / a_half

               end if

               if (i.le.ihi-1) then
                  rrl = qm(i+1,j,km,QRHO)
                  compl  = rrl*qm(i+1,j,km,nq)

                  rrnewl = rrl + drr
                  compnl = compl + dcompn

                  qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq) / a_half
               end if

            enddo
         enddo
      enddo

      do j = jlo, jhi
         do i = ilo, ihi

            pgyp = pgdnvy(i,j+1,km)
            pgym = pgdnvy(i,j,km)
            ugyp = ugdnvy(i,j+1,km)
            ugym = ugdnvy(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            if (i.ge.ilo+1) then
               ! Convert to conservation form
               rrr =     qp(i,j,km,QRHO)
               rur = rrr*qp(i,j,km,QU)
               rvr = rrr*qp(i,j,km,QV)
               rwr = rrr*qp(i,j,km,QW)
               ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
                                  qp(i,j,km,QW)**2)
               rer = qp(i,j,km,QREINT) + ekenr

               ! Add transverse terms
               rrnewr = rrr - hdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO )) &
                            - hdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO ))
               runewr = rur - hdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX  ))
               rvnewr = rvr - hdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY  ))
               rwnewr = rwr - hdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ  ))
               renewr = rer - hdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                            - hdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewr .lt. ZERO) then
                    rrnewr = rrr
                    runewr = rur
                    rvnewr = rvr
                    rwnewr = rwr
                    renewr = rer
               end if

               rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            end if

            if (i.le.ihi-1) then
               rrl =     qm(i+1,j,km,QRHO)
               rul = rrl*qm(i+1,j,km,QU)
               rvl = rrl*qm(i+1,j,km,QV)
               rwl = rrl*qm(i+1,j,km,QW)
               ekenl = HALF*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
                                  qm(i+1,j,km,QW)**2)
               rel = qm(i+1,j,km,QREINT) + ekenl

            ! Add transverse terms
               rrnewl = rrl - hdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO )) &
                            - hdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO ))
               runewl = rul - hdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX  ))
               rvnewl = rvl - hdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY  ))
               rwnewl = rwl - hdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ  )) &
                            - hdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ  ))
               renewl = rel - hdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                         - hdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN))

               ! Reset to original value if adding transverse terms made density negative
               if (rrnewl .lt. ZERO) then
                    rrnewl = rrl
                    runewl = rul
                    rvnewl = rvl
                    rwnewl = rwl
                    renewl = rel
               end if

               rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            end if

            ! These are used only if the (rho e) version goes negative 
            duyp = pgyp*ugyp - pgym*ugym
            pyav = HALF*(pgyp+pgym)
            duy = ugyp-ugym
            pynew = hdtdy*(duyp + pyav*duy*gamma_minus_1)

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = HALF*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = hdtdz*(duzp + pzav*duz*gamma_minus_1)
            ! End section

            ! Convert back to non-conservation form
            ! ************************************************************************* 
            if (i.ge.ilo+1) then
               qpo(i,j,km,QRHO  ) = rrnewr
               qpo(i,j,km,QU    ) = runewr/rrnewr
               qpo(i,j,km,QV    ) = rvnewr/rrnewr
               qpo(i,j,km,QW    ) = rwnewr/rrnewr

               qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt* srcQ(i,j,k3d,QREINT) / a_old
               qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT) * gamma_minus_1

               if (qpo(i,j,km,QPRES) .lt. small_pres) then
                   pnewr = qp(i,j,km,QPRES) - pynew - pznew
                   qpo(i,j,km,QPRES ) = pnewr + hdt* srcQ(i,j,k3d,QPRES ) / a_old
                   qpo(i,j,km,QREINT) = qpo(i,j,km,QPRES) / gamma_minus_1
               end if
            end if

            ! ************************************************************************* 
            if (i.le.ihi-1) then
               qmo(i+1,j,km,QRHO  ) = rrnewl
               qmo(i+1,j,km,QU    ) = runewl/rrnewl
               qmo(i+1,j,km,QV    ) = rvnewl/rrnewl
               qmo(i+1,j,km,QW    ) = rwnewl/rrnewl

               qmo(i+1,j,km,QREINT) = renewl - rhoekenl + hdt* srcQ(i,j,k3d,QREINT) / a_old
               qmo(i+1,j,km,QPRES) = qmo(i+1,j,km,QREINT) * gamma_minus_1

               if (qmo(i+1,j,km,QPRES) .lt. small_pres) then
                   pnewl = qm(i+1,j,km,QPRES) - pynew - pznew
                   qmo(i+1,j,km,QPRES ) = pnewl + hdt* srcQ(i,j,k3d,QPRES ) / a_old
                   qmo(i+1,j,km,QREINT) = qmo(i+1,j,km,QPRES) / gamma_minus_1
               end if
            end if
            ! ************************************************************************* 

         enddo
      enddo

      ! Version_2 = 0: we add the full source term here
      ! Version_2 = 1: we already added the (piecewise constant) source term in the original tracing
      ! Version_2 = 2: we already added the (parabolic         ) source term in the original tracing
      ! Version_2 = 3: we project then add the source term in trace*_src_3d

      if (version_2 .eq. 0) then
         do j = jlo, jhi
            do i = ilo, ihi
               if (i.ge.ilo+1) then
                  qpo(i,j,km,QRHO) = qpo(i,j,km,QRHO) + hdt*srcQ(i,j,k3d,QRHO  ) / a_old
                  qpo(i,j,km,QU  ) = qpo(i,j,km,QU  ) + hdt*srcQ(i,j,k3d,QU) / a_old
                  qpo(i,j,km,QV  ) = qpo(i,j,km,QV  ) + hdt*srcQ(i,j,k3d,QV) / a_old
                  qpo(i,j,km,QW  ) = qpo(i,j,km,QW  ) + hdt*srcQ(i,j,k3d,QW) / a_old
               end if
   
               if (i.le.ihi-1) then
                  qmo(i+1,j,km,QRHO) = qmo(i+1,j,km,QRHO) + hdt*srcQ(i,j,k3d,QRHO  ) / a_old
                  qmo(i+1,j,km,QU  ) = qmo(i+1,j,km,QU  ) + hdt*srcQ(i,j,k3d,QU) / a_old
                  qmo(i+1,j,km,QV  ) = qmo(i+1,j,km,QV  ) + hdt*srcQ(i,j,k3d,QV) / a_old
                  qmo(i+1,j,km,QW  ) = qmo(i+1,j,km,QW  ) + hdt*srcQ(i,j,k3d,QW) / a_old
               end if
            enddo
         enddo
      endif

      end subroutine transyz

end module transverse_module
