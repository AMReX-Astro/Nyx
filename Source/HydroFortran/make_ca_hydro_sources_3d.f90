
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_ca_make_hydro_sources(time,lo,hi,&
           uin,ulo, uhi, &
           ugdnvx_out,ugdnvx_lo,ugdnvx_hi,&
           ugdnvy_out,ugdnvy_lo,ugdnvy_hi ,&
           ugdnvz_out,ugdnvz_lo,ugdnvz_hi ,&
           src,src_lo,src_hi ,&
           hydro_src,hsrc_lo,hsrc_hi ,&
           divu_cc,d_lo,d_hi ,&
           grav,gv_lo,gv_hi ,&
           delta,dt,&
           flux1,flux1_lo,flux1_hi ,&
           flux2,flux2_lo,flux2_hi ,&
           flux3,flux3_lo,flux3_hi ,&
           a_old,a_new,print_fortran_warnings) &
           bind(C, name="fort_ca_make_hydro_sources")

      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use meth_params_module, only : QVAR, NVAR, NHYP, normalize_species, &
           QC, QPRES, use_flattening, UFS
      use flatten_module
      use amrex_constants_module
      use advection_module

      implicit none

      integer lo(3),hi(3),print_fortran_warnings
      integer ulo(3),uhi(3)
      integer ugdnvx_lo(3),ugdnvx_hi(3)
      integer ugdnvy_lo(3),ugdnvy_hi(3)
      integer ugdnvz_lo(3),ugdnvz_hi(3)
      integer flux1_lo(3),flux1_hi(3)
      integer flux2_lo(3),flux2_hi(3)
      integer flux3_lo(3),flux3_hi(3)
      integer src_lo(3),src_hi(3)
      integer hsrc_lo(3),hsrc_hi(3)
      integer d_lo(3),d_hi(3)
      integer gv_lo(3),gv_hi(3)
      real(rt) uin(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),NVAR)
      real(rt) ugdnvx_out(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3))
      real(rt) ugdnvy_out(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3))
      real(rt) ugdnvz_out(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3))
      real(rt)   src(  src_lo(1):src_hi(1),    src_lo(2):src_hi(2),     src_lo(3):src_hi(3),  NVAR)
      real(rt) hydro_src(hsrc_lo(1):hsrc_hi(1),hsrc_lo(2):hsrc_hi(2),hsrc_lo(3):hsrc_hi(3),NVAR)
      real(rt)  divu_cc(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
      real(rt)  grav( gv_lo(1):gv_hi(1),  gv_lo(2):gv_hi(2),   gv_lo(3):gv_hi(3),    3)
      real(rt) flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),NVAR)
      real(rt) flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),NVAR)
      real(rt) flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),NVAR)
      real(rt) delta(3),dt,time
      real(rt) a_old, a_new

      ! Automatic arrays for workspace
      real(rt), pointer :: q(:,:,:,:)
      real(rt), pointer :: qaux(:,:,:,:)
      real(rt), pointer :: flatn(:,:,:)
      real(rt), pointer :: c(:,:,:)
      real(rt), pointer :: csml(:,:,:)
      real(rt), pointer :: divu_nd(:,:,:)
      real(rt), pointer :: srcQ(:,:,:,:)

      real(rt) dx,dy,dz
      integer ngq,ngf
      integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
      integer qlo(3), qhi(3), loq(3), hiq(3), n, tmp_hi(3), glo(3), ghi(3)
      ngq = NHYP
      ngf = 1

      q_l1 = lo(1)-NHYP
      q_l2 = lo(2)-NHYP
      q_l3 = lo(3)-NHYP
      q_h1 = hi(1)+NHYP
      q_h2 = hi(2)+NHYP
      q_h3 = hi(3)+NHYP

      srcq_l1 = lo(1)-NHYP
      srcq_l2 = lo(2)-NHYP
      srcq_l3 = lo(3)-NHYP
      srcq_h1 = hi(1)+NHYP
      srcq_h2 = hi(2)+NHYP
      srcq_h3 = hi(3)+NHYP

      call amrex_allocate(     q, lo-NHYP, hi+NHYP, QVAR)
      call amrex_allocate(  qaux, lo-NHYP, hi+NHYP,    1)
      call amrex_allocate( flatn, lo-NHYP, hi+NHYP      )
      call amrex_allocate(     c, lo-NHYP, hi+NHYP      )
      call amrex_allocate(  csml, lo-NHYP, hi+NHYP      )

      call amrex_allocate(   srcQ, lo-NHYP, hi+NHYP, QVAR)
      call amrex_allocate(divu_nd, lo  , hi+1)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! 1) Translate conserved variables (u) to primitive variables (q).
      ! 2) Compute sound speeds (c) 
      !    Note that (q,c,csml,flatn) are all dimensioned the same
      !    and set to correspond to coordinates of (lo:hi)
      ! 3) Translate source terms
      qlo(1)=q_l1
      qlo(2)=q_l2
      qlo(3)=q_l3
      qhi(1)=q_h1
      qhi(2)=q_h2
      qhi(3)=q_h3

      ! Note that for now, csml is passed seperately
      ! It's unclear whether qaux will bw generally useful
      call ca_ctoprim(qlo,qhi,uin,ulo, uhi, &
           q,qlo, qhi, &
           qaux,qlo, qhi,csml)

      c(:,:,:)=qaux(:,:,:,QC)
      if (use_flattening == 1) then
         call ca_uflatten(lo-1,hi+1,&
         q,qlo,qhi, &
         flatn, qlo, qhi, QPRES)
      else
         flatn = ONE
      endif

      call ca_srctoprim(qlo,qhi, &
           q,qlo, qhi, &
           qaux,qlo, qhi, &
           grav,gv_lo, gv_hi, &
           src, qlo, qhi, &
           srcQ, qlo, qhi, a_old, a_new, dt)

      ! Compute hyperbolic fluxes using unsplit Godunov
      call umeth3d(q,c,csml,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                   srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
                   flux1,flux1_lo(1),flux1_lo(2),flux1_lo(3),flux1_hi(1),flux1_hi(2),flux1_hi(3), &
                   flux2,flux2_lo(1),flux2_lo(2),flux2_lo(3),flux2_hi(1),flux2_hi(2),flux2_hi(3), &
                   flux3,flux3_lo(1),flux3_lo(2),flux3_lo(3),flux3_hi(1),flux3_hi(2),flux3_hi(3), &
                   grav,gv_lo(1),gv_lo(2),gv_lo(3),gv_hi(1),gv_hi(2),gv_hi(3), &
                   ugdnvx_out,ugdnvx_lo(1),ugdnvx_lo(2),ugdnvx_lo(3),ugdnvx_hi(1),ugdnvx_hi(2),ugdnvx_hi(3), &
                   ugdnvy_out,ugdnvy_lo(1),ugdnvy_lo(2),ugdnvy_lo(3),ugdnvy_hi(1),ugdnvy_hi(2),ugdnvy_hi(3), &
                   ugdnvz_out,ugdnvz_lo(1),ugdnvz_lo(2),ugdnvz_lo(3),ugdnvz_hi(1),ugdnvz_hi(2),ugdnvz_hi(3), &
                   divu_cc,d_lo(1),d_lo(2),d_lo(3),d_hi(1),d_hi(2),d_hi(3), &
                   a_old,a_new,print_fortran_warnings)

      ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi+1, &
           q,qlo, qhi, &
           delta, divu_nd, lo, hi+1)

      tmp_hi=hi
      tmp_hi(1)=tmp_hi(1)+1
     call ca_apply_av(lo, tmp_hi, 1, delta, &
          divu_nd, lo, hi+1, &
          uin, ulo, uhi, &
          flux1, flux1_lo, flux1_hi, dt)
     tmp_hi=hi
      tmp_hi(2)=tmp_hi(2)+1
     call ca_apply_av(lo, tmp_hi, 2, delta, &
          divu_nd, lo, hi+1, &
          uin, ulo, uhi, &
          flux2, flux2_lo, flux2_hi, dt)
     tmp_hi=hi
     tmp_hi(3)=tmp_hi(3)+1
     call ca_apply_av(lo, tmp_hi, 3, delta, &
          divu_nd, lo, hi+1, &
          uin, ulo, uhi, &
          flux3, flux3_lo, flux3_hi, dt)

     if (UFS .gt. 0 .and. normalize_species .eq. 1) then
     tmp_hi=hi
      tmp_hi(1)=tmp_hi(1)+1
      call ca_normalize_species_fluxes(lo, tmp_hi, &
               flux1, flux1_lo, flux1_hi)
     tmp_hi=hi
     tmp_hi(2)=tmp_hi(2)+1
     call ca_normalize_species_fluxes(lo, tmp_hi, &
          flux2, flux2_lo, flux2_hi)
     tmp_hi=hi
     tmp_hi(3)=tmp_hi(3)+1
     call ca_normalize_species_fluxes(lo, tmp_hi, &
          flux3, flux3_lo, flux3_hi)
     endif
      
      ! Conservative update to make hydro sources
      call ca_consup(uin,ulo, uhi, &
                  hydro_src , hsrc_lo, hsrc_hi, &
                  flux1,flux1_lo,flux1_hi, &
                  flux2,flux2_lo, flux2_hi, &
                  flux3,flux3_lo,flux3_hi, &
                  divu_nd,lo, hi+1, &
                  divu_cc,d_lo,d_hi, &
                  lo,hi,delta,dt,a_old,a_new)

      ! We are done with these here so can go ahead and free up the space.
      call amrex_deallocate(q)
      call amrex_deallocate(qaux)
      call amrex_deallocate(flatn)
      call amrex_deallocate(c)
      call amrex_deallocate(csml)
      call amrex_deallocate(divu_nd)
      call amrex_deallocate(srcQ)

      end subroutine fort_ca_make_hydro_sources

