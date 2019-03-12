
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_ca_make_hydro_sources(time,lo,hi,&
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
           divu_cc,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta,dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           a_old,a_new,print_fortran_warnings) &
           bind(C, name="fort_ca_make_hydro_sources")

      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use meth_params_module, only : QVAR, NVAR, NHYP, normalize_species, &
           QC, QPRES, use_flattening
      use flatten_module
      use amrex_constants_module
      use advection_module

      implicit none

      integer lo(3),hi(3),print_fortran_warnings
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      real(rt)   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,     uin_l3:uin_h3,  NVAR)
      real(rt) ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      real(rt) ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      real(rt) ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      real(rt)   src(  src_l1:src_h1,    src_l2:src_h2,     src_l3:src_h3,  NVAR)
      real(rt) hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt)  divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
      real(rt)  grav( gv_l1:gv_h1,  gv_l2:gv_h2,   gv_l3:gv_h3,    3)
      real(rt) flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
      real(rt) flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
      real(rt) flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
      real(rt) delta(3),dt,time
      real(rt) a_old, a_new

      ! Automatic arrays for workspace
      real(rt), pointer :: q(:,:,:,:)
      real(rt), pointer :: q2(:,:,:,:)
      real(rt), pointer :: qaux(:,:,:,:)
      real(rt), pointer :: flatn(:,:,:)
      real(rt), pointer :: c(:,:,:)
      real(rt), pointer :: csml(:,:,:)
      real(rt), pointer :: divu_nd(:,:,:)
      real(rt), pointer :: srcQ(:,:,:,:)
      real(rt), pointer :: srcQ2(:,:,:,:)

      real(rt) dx,dy,dz
      integer ngq,ngf
      integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
      integer ulo(3), uhi(3), qlo(3), qhi(3), loq(3), hiq(3), n, tmp_hi(3)
      integer flux1_lo(3), flux1_hi(3), flux2_lo(3), flux2_hi(3), flux3_lo(3), flux3_hi(3)
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
      call amrex_allocate(     q2, lo-NHYP, hi+NHYP, QVAR)
      call amrex_allocate(  qaux, lo-NHYP, hi+NHYP,    1)
      call amrex_allocate( flatn, lo-NHYP, hi+NHYP      )
      call amrex_allocate(     c, lo-NHYP, hi+NHYP      )
      call amrex_allocate(  csml, lo-NHYP, hi+NHYP      )

      call amrex_allocate(   srcQ, lo-NHYP, hi+NHYP, QVAR)
      call amrex_allocate(   srcQ2, lo-NHYP, hi+NHYP, QVAR)
      call amrex_allocate(divu_nd, lo  , hi+1)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! 1) Translate conserved variables (u) to primitive variables (q).
      ! 2) Compute sound speeds (c) 
      !    Note that (q,c,csml,flatn) are all dimensioned the same
      !    and set to correspond to coordinates of (lo:hi)
      ! 3) Translate source terms
      ulo(1)=uin_l1
      ulo(2)=uin_l2
      ulo(3)=uin_l3
      uhi(1)=uin_h1
      uhi(2)=uin_h2
      uhi(3)=uin_h3

      qlo(1)=q_l1
      qlo(2)=q_l2
      qlo(3)=q_l3
      qhi(1)=q_h1
      qhi(2)=q_h2
      qhi(3)=q_h3

      flux1_lo(1)=flux1_l1
      flux1_lo(2)=flux1_l2
      flux1_lo(3)=flux1_l3
      flux1_hi(1)=flux1_h1
      flux1_hi(2)=flux1_h2
      flux1_hi(3)=flux1_h3
      flux2_lo(1)=flux2_l1
      flux2_lo(2)=flux2_l2
      flux2_lo(3)=flux2_l3
      flux2_hi(1)=flux2_h1
      flux2_hi(2)=flux2_h2
      flux2_hi(3)=flux2_h3
      flux3_lo(1)=flux3_l1
      flux3_lo(2)=flux3_l2
      flux3_lo(3)=flux3_l3
      flux3_hi(1)=flux3_h1
      flux3_hi(2)=flux3_h2
      flux3_hi(3)=flux3_h3

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
           grav,qlo, qhi, &
           src, qlo, qhi, &
           srcQ, qlo, qhi, a_old, a_new, dt)

      ! Compute hyperbolic fluxes using unsplit Godunov
      call umeth3d(q,c,csml,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                   srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                   flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                   flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                   grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                   ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                   ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                   ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                   divu_cc,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                   a_old,a_new,print_fortran_warnings)

      ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call make_divu_nd(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                        dx,dy,dz,divu_nd,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

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
      
      ! Conservative update
      call ca_consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                  hydro_src , hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3, &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  divu_nd,divu_cc,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                  lo,hi,dx,dy,dz,dt,a_old,a_new)

      ! We are done with these here so can go ahead and free up the space.
      call amrex_deallocate(q)
      call amrex_deallocate(flatn)
      call amrex_deallocate(c)
      call amrex_deallocate(csml)
      call amrex_deallocate(divu_nd)
      call amrex_deallocate(srcQ)

      end subroutine fort_ca_make_hydro_sources

