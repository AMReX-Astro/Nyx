module advection_module
      use amrex_fort_module, only : rt => amrex_real
contains
  
subroutine ca_ctoprim(lo, hi, &
       uin, uin_lo, uin_hi, &
       q,     q_lo,   q_hi, &
       qaux, qa_lo,  qa_hi, csml) bind(c,name='ca_ctoprim')

    use network, only : nspec, naux
    use eos_module, only : nyx_eos_soundspeed
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, &
                                     UEDEN, UEINT, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, QFA, QFS, &!                                   UEDEN, UEINT, UTEMP, &
!                                   QRHO, QU, QV, QW, &
!                                   QREINT, QPRES, QTEMP, QGAME, QFS, QFX, &
                                     !                                   NQ, QC, QGAMC, QGC, QDPDR, QDPDE,
                                   QC, NQAUX, &
                                   npassive, upass_map, qpass_map, &
                                   small_dens, small_pres, &
                                   gamma_const, gamma_minus_1, use_flattening

    use amrex_constants_module, only: ZERO, HALF, ONE
    use amrex_error_module
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt)        , intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)        , intent(inout) :: csml(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3))

    real(rt)        , parameter :: small = 1.e-8_rt

    integer          :: i, j, k, g
    integer          :: n, iq, ipassive
    real(rt)         :: kineng, rhoinv
    real(rt)         :: vel(3)
    real(rt) :: a_half, a_dot
    real(rt) :: dtdxaold, dtdyaold, dtdzaold, small_pres_over_dens

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)

#ifndef AMREX_USE_CUDA
          do i = lo(1), hi(1)
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::ctoprim ',i, j, k
                print *,'>>> ... negative density ', uin(i,j,k,URHO)
                call amrex_error("Error:: advection_util_nd.f90 :: ctoprim")
             else if (uin(i,j,k,URHO) .lt. small_dens) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::ctoprim ',i, j, k
                print *,'>>> ... small density ', uin(i,j,k,URHO)
                call amrex_error("Error:: advection_util_nd.f90 :: ctoprim")
             endif
          end do
#endif
          do i = lo(1), hi(1)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = uin(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Note the dual energy formulation is enforced elsewhere
             q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv

             ! Note we assume QTEMP doesn't exists
!             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

         enddo
       enddo
    enddo

    ! Load passively advected quatities into q
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       iq = qpass_map(ipassive)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                q(i,j,k,iq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    small_pres_over_dens = small_pres / small_dens
        
    ! get gamc, p, T, c, csml using q state
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Define the soundspeed from the EOS
             !             call nyx_eos_soundspeed(c(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT))
             call nyx_eos_soundspeed(qaux(i,j,k,QC), q(i,j,k,QRHO), q(i,j,k,QREINT))

             ! Set csmal based on small_pres and small_dens
             csml(i,j,k) = sqrt(gamma_const * small_pres_over_dens)
             
             ! Convert "e" back to "rho e"
             q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)
             
             ! Pressure = (gamma - 1) * rho * e
             q(i,j,k,QPRES) = gamma_minus_1 * q(i,j,k,QREINT)

             ! Note we are not using eos_type, and don't require an update call to eos

             ! Note we are not storing dpdr and dpde since these are computable based on q
          enddo
       enddo
    enddo
    ! Compute flattening coef for slope calculations in driver
      
  end subroutine ca_ctoprim



  subroutine ca_srctoprim(lo, hi, &
       q,     q_lo,   q_hi, &
       qaux, qa_lo,  qa_hi, &
       grav,  g_lo,   g_hi, &
       src, src_lo, src_hi, &
       srcQ,srQ_lo, srQ_hi, a_old, a_new, dt) bind(c,name='ca_srctoprim')

    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, NQAUX, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     npassive, upass_map, qpass_map, &
                                     nadv, small_dens, small_pres, &
                                     gamma_const, gamma_minus_1, use_flattening, QFA, QFS, UFS, UFA
    use amrex_constants_module, only: ZERO, HALF, ONE, THREE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3),   qa_hi(3)
    integer, intent(in) ::   g_lo(3),   g_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)
    real(rt), intent(in) :: a_new, a_old, dt
!grav explicitly includes spacedim
    real(rt)        , intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt)        , intent(in   ) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)        , intent(in   ) :: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)
    real(rt)        , intent(in   ) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt)        , intent(inout) :: srcQ(srQ_lo(1):srQ_hi(1),srQ_lo(2):srQ_hi(2),srQ_lo(3):srQ_hi(3),QVAR)

    integer          :: i, j, k
    integer          :: n, iq, ipassive, iadv, ispec
    real(rt)         :: rhoinv, a_half, a_dot, dpde, dpdr

    !$gpu
    a_half = HALF * (a_old + a_new)
    a_dot   = (a_new - a_old) / dt
      
    srcQ(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

    ! compute srcQ terms
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoinv = ONE / q(i,j,k,QRHO)

             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
                  grav(i,j,k,1)
             srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
                  grav(i,j,k,2)
             srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
                  grav(i,j,k,3)
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
                  q(i,j,k,QV)*src(i,j,k,UMY) - &
                  q(i,j,k,QW)*src(i,j,k,UMZ) - &
                  a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)


             dpde = gamma_minus_1 * q(i,j,k,QRHO)
             dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO)
             srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoInv &
                  + dpdr * srcQ(i,j,k,QRHO)

             !!!!!!!!!!!UFS and UFA mapping should be identical to ipassive
!             if (UFS .gt. 0) then
!                  do ispec = 1,nspec+naux
!                     srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
!                  enddo
!               end if ! UFS > 0

!               do iadv = 1,nadv
!                  srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
!               enddo


          enddo
       enddo
    enddo

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       iq = qpass_map(ipassive)

       ! we already accounted for velocities above
       if (iq == QU .or. iq == QV .or. iq == QW) cycle

       ! we may not be including the ability to have species sources,
       ! so check to make sure that we are < NQSRC
       if (iq > QVAR) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                srcQ(i,j,k,iq) = ( src(i,j,k,n) - q(i,j,k,iq) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo
          enddo
       enddo

    enddo

  end subroutine ca_srctoprim

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  !> @brief this computes the *node-centered* divergence
  !!
  subroutine divu(lo, hi, &
       q, q_lo, q_hi, &
       dx, div, div_lo, div_hi) bind(C, name='divu')

    use meth_params_module, only : QU, QV, QW, QVAR
    use amrex_constants_module, only : HALF, FOURTH, ONE, ZERO
    use prob_params_module, only : dg, coord_type
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv
    real(rt) :: rl, rr, rc, ul, ur, vt, vb

    !$gpu

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ux = FOURTH*( &
                  + q(i        ,j        ,k        ,QU) - q(i-1*dg(1),j        ,k        ,QU) &
                  + q(i        ,j        ,k-1*dg(3),QU) - q(i-1*dg(1),j        ,k-1*dg(3),QU) &
                  + q(i        ,j-1*dg(2),k        ,QU) - q(i-1*dg(1),j-1*dg(2),k        ,QU) &
                  + q(i        ,j-1*dg(2),k-1*dg(3),QU) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QU) ) * dxinv

             vy = FOURTH*( &
                  + q(i        ,j        ,k        ,QV) - q(i        ,j-1*dg(2),k        ,QV) &
                  + q(i        ,j        ,k-1*dg(3),QV) - q(i        ,j-1*dg(2),k-1*dg(3),QV) &
                  + q(i-1*dg(1),j        ,k        ,QV) - q(i-1*dg(1),j-1*dg(2),k        ,QV) &
                  + q(i-1*dg(1),j        ,k-1*dg(3),QV) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QV) ) * dyinv

             wz = FOURTH*( &
                  + q(i        ,j        ,k        ,QW) - q(i        ,j        ,k-1*dg(3),QW) &
                  + q(i        ,j-1*dg(2),k        ,QW) - q(i        ,j-1*dg(2),k-1*dg(3),QW) &
                  + q(i-1*dg(1),j        ,k        ,QW) - q(i-1*dg(1),j        ,k-1*dg(3),QW) &
                  + q(i-1*dg(1),j-1*dg(2),k        ,QW) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu

  
    subroutine ca_apply_av(lo, hi, idir, dx, &
       div, div_lo, div_hi, &
       uin, uin_lo, uin_hi, &
       flux, f_lo, f_hi, dt) bind(c, name="apply_av")

    use amrex_constants_module, only: ZERO, FOURTH, ONE
    use meth_params_module, only: NVAR, difmag
    use prob_params_module, only: dg

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3), dt
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1
    real(rt) :: area(3)
    real(rt) :: vol, volinv

!!!    !$gpu

    vol     = dx(1) * dx(2) * dx(3)
    volinv  = ONE / vol

    area(1)   =         dx(2) * dx(3)
    area(2)   = dx(1)         * dx(3)
    area(3)   = dx(1) * dx(2)

    do n = 1, NVAR

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i,j+1*dg(2),k        ) + &
                        div(i,j,k+1*dg(3)) + div(i,j+1*dg(2),k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i-1*dg(1),j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i+1*dg(1),j,k        ) + &
                        div(i,j,k+1*dg(3)) + div(i+1*dg(1),j,k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j-1*dg(2),k,n))

                else

                   div1 = FOURTH * (div(i,j        ,k) + div(i+1*dg(1),j        ,k) + &
                        div(i,j+1*dg(2),k) + div(i+1*dg(1),j+1*dg(2),k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-1*dg(3),n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1
                flux(i,j,k,n) = flux(i,j,k,n) * area(idir) * dt

             end do
          end do
       end do

    end do

  end subroutine ca_apply_av

    !> @brief Normalize the fluxes of the mass fractions so that
  !! they sum to 0.  This is essentially the CMA procedure that is
  !! defined in Plewa & Muller, 1999, A&A, 342, 179.
  !!
  subroutine ca_normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")

    use network, only: nspec
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UFS

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n
    real(rt) :: sum, fac

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             sum = ZERO

             do n = UFS, UFS+nspec-1
                sum = sum + flux(i,j,k,n)
             end do

             if (sum .ne. ZERO) then
                fac = flux(i,j,k,URHO) / sum
             else
                fac = ONE
             end if

             do n = UFS, UFS+nspec-1
                flux(i,j,k,n) = flux(i,j,k,n) * fac
             end do

          end do
       end do
    end do

  end subroutine ca_normalize_species_fluxes
! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine ca_consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                      hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
                      flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                      flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                      flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                      divu_nd,divu_cc,d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                      lo,hi,dx,dy,dz,dt,a_old,a_new)

      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : difmag, NVAR, URHO, UMX, UMZ, &
           UEDEN, UEINT, UFS, normalize_species, gamma_minus_1

      implicit none

      integer lo(3), hi(3)
      integer   uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer   hsrc_l1,  hsrc_l2,  hsrc_l3,  hsrc_h1,  hsrc_h2,  hsrc_h3
      integer d_l1,   d_l2,   d_l3,   d_h1,   d_h2,   d_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3

      real(rt)  :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      real(rt)  :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt)  :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      real(rt)  :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      real(rt)  :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      real(rt)  :: divu_nd(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(rt)  :: divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
      real(rt)  :: dx, dy, dz, dt, a_old, a_new

      real(rt) :: div1, a_half, a_oldsq, a_newsq
      real(rt) :: area1, area2, area3
      real(rt) :: vol, volinv, a_newsq_inv
      real(rt) :: a_half_inv, a_new_inv, dt_a_new
      integer  :: i, j, k, n

      a_half  = HALF * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      vol     = dx * dy * dz
      volinv  = ONE / vol

      a_half_inv  = ONE / a_half
      a_new_inv   = ONE / a_new
      dt_a_new    = dt / a_new
      a_newsq_inv = ONE / a_newsq

      do n = 1, NVAR
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  ! Density
                  if (n .eq. URHO) then
                      hydro_src(i,j,k,n) = &
                          ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv  ) * a_half_inv

                  ! Momentum
                  else if (n .ge. UMX .and. n .le. UMZ) then
                     hydro_src(i,j,k,n) =  &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv

                  ! (rho E)
                  else if (n .eq. UEDEN) then
                     hydro_src(i,j,k,n) =  &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * a_half * volinv &
                          +   a_half * (a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT)

                  ! (rho e)
                  else if (n .eq. UEINT) then

                     hydro_src(i,j,k,n) =  &
                          ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                           +flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                           +flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * a_half * volinv &
                           +a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                           -a_half*dt*(gamma_minus_1 * uin(i,j,k,n)) * divu_cc(i,j,k)

                  ! (rho X_i) and (rho adv_i) and (rho aux_i)
                  else
                     hydro_src(i,j,k,n) = &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv 
                  endif

                  hydro_src(i,j,k,n) = hydro_src(i,j,k,n) / dt

               enddo
            enddo
         enddo
      enddo

      do n = 1, NVAR
         if (n .eq. URHO) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_half_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_half_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_half_inv
         else if (n.ge.UMX .and. n.le.UMZ) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_new_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_new_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_new_inv
         else if (n.eq.UEINT .or. n.eq.UEDEN) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * (a_half * a_newsq_inv)
            flux2(:,:,:,n) = flux2(:,:,:,n) * (a_half * a_newsq_inv)
            flux3(:,:,:,n) = flux3(:,:,:,n) * (a_half * a_newsq_inv)
         else
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_half_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_half_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_half_inv
         end if
      end do

      end subroutine ca_consup
  
end module advection_module
