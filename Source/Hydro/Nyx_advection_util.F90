module advection_module
      use amrex_fort_module, only : rt => amrex_real
contains

     AMREX_CUDA_FORT_DEVICE subroutine update_state(lo,hi, &
                             uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                             uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                             src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                             hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
                             divu_cc,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                             sum_state, s_lo, s_hi, &
                             dt,a_old,a_new,print_fortran_warnings) &
                             bind(C, name="ca_fort_update_state")
 
      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, &
                                     gamma_minus_1, normalize_species
      use ca_enforce_module, only : ca_enforce_nonnegative_species
      use enforce_density_module, only : ca_enforce_minimum_density, ca_enforce_minimum_density_1cell
      use amrex_error_module, only : amrex_error

      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) ::print_fortran_warnings
      integer, intent(in) ::  uin_l1,   uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer, intent(in) ::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer, intent(in) ::   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
      integer, intent(in) ::  hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
      integer, intent(in) ::  d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer, intent(in) :: s_lo(3), s_hi(3)

      real(rt), intent(in)  ::       uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt), intent(out) ::      uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      real(rt), intent(in)  ::       src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NVAR)
      real(rt), intent(in)  :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt), intent(in)  ::   divu_cc(   d_l1:   d_h1,   d_l2:   d_h2,   d_l3:   d_h3)
      real(rt), intent(inout) :: sum_state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3),NVAR)
      real(rt), intent(in)  ::  dt, a_old, a_new

      real(rt) :: a_half, a_oldsq, a_newsq
      real(rt) :: a_new_inv, a_newsq_inv, a_half_inv, dt_a_new
      integer  :: i, j, k, n
      integer :: noprint

      noprint = 0

      a_half     = HALF * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new

      a_half_inv  = ONE / a_half
      a_new_inv   = ONE / a_new
      a_newsq_inv = ONE / a_newsq

      do n = 1, NVAR

         ! Actually do the update here
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  ! Density
                  if (n .eq. URHO) then
                     uout(i,j,k,n) = uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + dt *  src(i,j,k,n) * a_half_inv

                  ! Momentum
                  else if (n .ge. UMX .and. n .le. UMZ) then
                     uout(i,j,k,n) = a_old * uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + dt   * src(i,j,k,n)
                     uout(i,j,k,n) = uout(i,j,k,n) * a_new_inv

                  ! (rho E)
                  else if (n .eq. UEDEN) then
                     uout(i,j,k,n) =  a_oldsq * uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + a_half  * dt * src(i,j,k,n)  
                     uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv

                  ! (rho e)
                  else if (n .eq. UEINT) then

                     uout(i,j,k,n) =  a_oldsq*uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + a_half * dt * src(i,j,k,n) 

                     uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv 

!                    We don't do this here because we are adding all of this term explicitly in hydro_src
!                    uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv / &
!                        ( ONE + a_half * dt * (HALF * gamma_minus_1 * divu_cc(i,j,k)) * a_newsq_inv )

                  ! (rho X_i) and (rho adv_i) and (rho aux_i)
                  else
                     uout(i,j,k,n) = uin(i,j,k,n) +  dt * hydro_src(i,j,k,n) &
                                    + dt * src(i,j,k,n) * a_half_inv

                  endif

               enddo
            enddo
         enddo
      enddo

      ! Enforce the density >= small_dens.  Make sure we do this immediately after consup.
      call ca_enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                   uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                   sum_state,s_lo,s_hi, &
                                   lo,hi,print_fortran_warnings)

!      call ca_enforce_minimum_density_1cell(lo, hi, &
!                                   uout,s_lo, s_hi, &
!                                   print_fortran_warnings)

      ! Enforce species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                       uout_h1,uout_h2,uout_h3,lo,hi,noprint)

      ! Re-normalize the species
      if (normalize_species .eq. 1) then
         call ca_normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                    lo,hi)
      end if

    end subroutine update_state

! :::
! ::: ------------------------------------------------------------------
! :::

    AMREX_CUDA_FORT_DEVICE subroutine fort_add_grav_source(lo,hi,&
                               uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               dt,a_old,a_new) &
                               bind(C, name="ca_fort_add_grav_source")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
           UEDEN, grav_source_type

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer  gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3

      real(rt)  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt) uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      real(rt) grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,  gv_l3:  gv_h3,3)
      real(rt) dt
      real(rt) a_old, a_new

      real(rt) :: a_half, a_oldsq, a_newsq, a_newsq_inv
      real(rt) :: rho
      real(rt) :: SrU, SrV, SrW, SrE
      real(rt) :: rhoInv, dt_a_new
      real(rt) :: old_rhoeint, new_rhoeint, old_ke, new_ke
      integer          :: i, j, k

      a_half  = 0.5d0 * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      a_newsq_inv = 1.d0 / a_newsq

      dt_a_new    = dt / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      ! Add gravitational source terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 
               old_rhoeint = uout(i,j,k,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rho    = uin(i,j,k,URHO)
               rhoInv = 1.0d0 / rho

               SrU = rho * grav(i,j,k,1)
               SrV = rho * grav(i,j,k,2)
               SrW = rho * grav(i,j,k,3)

               ! We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
               uout(i,j,k,UMX)   = uout(i,j,k,UMX) + SrU * dt_a_new
               uout(i,j,k,UMY)   = uout(i,j,k,UMY) + SrV * dt_a_new
               uout(i,j,k,UMZ)   = uout(i,j,k,UMZ) + SrW * dt_a_new

               if (grav_source_type .eq. 1) then

                   ! This does work (in 1-d)
                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   SrE = uin(i,j,k,UMX) * grav(i,j,k,1) + &
                         uin(i,j,k,UMY) * grav(i,j,k,2) + &
                         uin(i,j,k,UMZ) * grav(i,j,k,3)
                   uout(i,j,k,UEDEN) = (a_newsq*uout(i,j,k,UEDEN) + SrE * (dt*a_half)) * a_newsq_inv

               else if (grav_source_type .eq. 3) then

                   new_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                     uout(i,j,k,URHO) 
                   uout(i,j,k,UEDEN) = old_rhoeint + new_ke
               else 
#ifndef AMREX_USE_CUDA
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
#else
                  uout(i,j,k,UEDEN)=-3.d200
#endif
               end if

            enddo
         enddo
      enddo

      end subroutine fort_add_grav_source

  
AMREX_CUDA_FORT_DEVICE subroutine ca_ctoprim(lo, hi, &
       uin, uin_lo, uin_hi, &
       q,     q_lo,   q_hi, &
       qaux, qa_lo,  qa_hi) bind(c,name='ca_ctoprim')

    use eos_module, only : nyx_eos_soundspeed
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, &
                                     UEDEN, UEINT, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, QFA, QFS, &!                                   UEDEN, UEINT, UTEMP, &
!                                   QRHO, QU, QV, QW, &
!                                   QREINT, QPRES, QTEMP, QGAME, QFS, QFX, &
                                     !                                   QVAR, QC, QGAMC, QGC, QDPDR, QDPDE,
                                   QC, NQAUX, &
                                   npassive, upass_map, qpass_map, &
                                   small_dens, small_pres, &
                                   gamma_minus_1, use_flattening

    use amrex_constants_module, only: ZERO, HALF, ONE
    use amrex_error_module
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
!    integer, intent(in) :: ca_lo(3), ca_hi(3)

    real(rt)        , intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
!    real(rt)        , intent(inout) :: csml(ca_lo(1):ca_hi(1),ca_lo(2):ca_hi(2),ca_lo(3):ca_hi(3),1)

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
!             call nyx_eos_soundspeed(qaux(i,j,k,QC), q(i,j,k,QRHO), q(i,j,k,QREINT))

             qaux(i,j,k,QC) = sqrt((gamma_minus_1+ONE) * q(i,j,k,QREINT)*gamma_minus_1)

!             This is hard-coded into riemann since it doesn't vary per-cell
!             ! Set csmal based on small_pres and small_dens
!             csml(i,j,k,1) = sqrt((gamma_minus_1+ONE) * small_pres_over_dens)
             
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



  AMREX_CUDA_FORT_DEVICE subroutine ca_srctoprim(lo, hi, &
       q,     q_lo,   q_hi, &
       qaux, qa_lo,  qa_hi, &
       grav,  g_lo,   g_hi, &
       src, src_lo, src_hi, &
       srcQ,srQ_lo, srQ_hi, a_old, a_new, dt) bind(c,name='ca_srctoprim')

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, NQAUX, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, &
                                     npassive, upass_map, qpass_map, &
                                     nadv, small_dens, small_pres, &
                                     gamma_minus_1, use_flattening, use_srcQ_in_trace
    use amrex_constants_module, only: ZERO, HALF, ONE, THREE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3),   qa_hi(3)
    integer, intent(in) ::   g_lo(3),   g_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)
    real(rt), intent(in), value :: a_new, a_old, dt
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
!    srcQ = ZERO

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
                if(use_srcQ_in_trace.eq.1) then
                srcQ(i,j,k,iq) = ( src(i,j,k,n) - q(i,j,k,iq) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             else
                srcQ(i,j,k,iq) = ( src(i,j,k,n) ) / &
                    q(i,j,k,QRHO)
             endif

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
  AMREX_CUDA_FORT_DEVICE subroutine divu(lo, hi, &
       q, q_lo, q_hi, &
       dx, div, div_lo, div_hi) bind(C, name='divu')

    use meth_params_module, only : QU, QV, QW, QVAR
    use amrex_constants_module, only : HALF, FOURTH, ONE, ZERO
!    use prob_params_module, only : dg
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv, dg(3)
    real(rt) :: rl, rr, rc, ul, ur, vt, vb

    !$gpu

    dg = 1

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ux = FOURTH*( &
                  + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                  + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                  + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                  + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                  + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                  + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                  + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                  + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                  + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                  + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                  + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                  + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu

    AMREX_CUDA_FORT_DEVICE subroutine ca_apply_av(lo, hi, &
       div, div_lo, div_hi, &
       uin, uin_lo, uin_hi, &
       flux, f_lo, f_hi, idir, dx, dt) bind(c, name="apply_av")

    use amrex_constants_module, only: ZERO, FOURTH, ONE
    use meth_params_module, only: NVAR, difmag, use_area_dt_scale_apply

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir
    real(rt), intent(in   ), value :: dt

    real(rt) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1
    real(rt) :: area(3)
    real(rt) :: dg(3)
    real(rt) :: vol, volinv

    dg = 1
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

                   div1 = FOURTH * (div(i,j,k        ) + div(i,j+1      ,k        ) + &
                        div(i,j,k+1      ) + div(i,j+1      ,k+1      ))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i-1      ,j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i+1      ,j,k        ) + &
                        div(i,j,k+1      ) + div(i+1      ,j,k+1      ))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j-1      ,k,n))

                else

                   div1 = FOURTH * (div(i,j        ,k) + div(i+1      ,j        ,k) + &
                        div(i,j+1      ,k) + div(i+1      ,j+1      ,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-1      ,n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1
                if(use_area_dt_scale_apply.eq.1) then
                   flux(i,j,k,n) = flux(i,j,k,n) * area(idir) * dt
                endif

             end do
          end do
       end do

    end do

  end subroutine ca_apply_av

  !> @brief Normalize the fluxes of the mass fractions so that
  !! they sum to 0.  This is essentially the CMA procedure that is
  !! defined in Plewa & Muller, 1999, A&A, 342, 179.
  !!
  AMREX_CUDA_FORT_DEVICE subroutine ca_normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UFS, use_const_species, npassive, upass_map

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n, ipassive
    real(rt) :: sum, fac

    !$gpu

    if(UFS .gt. 0) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                sum = ZERO

                do ipassive = 1, npassive
                   n = upass_map(ipassive)
                   sum = sum + flux(i,j,k,n)
                end do
                
                if (sum .ne. ZERO) then
                   fac = flux(i,j,k,URHO) / sum
                else
                   fac = ONE
                end if

                do ipassive = 1, npassive
                   n = upass_map(ipassive)               
                   flux(i,j,k,n) = flux(i,j,k,n) * fac
                end do
                
             end do
          end do
       end do
    endif
    
  end subroutine ca_normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

      AMREX_CUDA_FORT_DEVICE subroutine ca_normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR, URHO, UFS, upass_map, npassive

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
      real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n,ipassive
      real(rt) :: fac,sum

      if (UFS .gt. 0) then

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            sum = 0.d0
            do ipassive = 1, npassive
               n = upass_map(ipassive)
!            do n = UFS, UFS+nspec-1
               sum = sum + u(i,j,k,n)
            end do
            if (sum .ne. 0.d0) then
               fac = u(i,j,k,URHO) / sum
            else
               fac = 1.d0
            end if
            do ipassive = 1, npassive
               n = upass_map(ipassive)
!            do n = UFS, UFS+nspec-1
               u(i,j,k,n) = u(i,j,k,n) * fac
            end do
         end do
      end do
      end do

      end if ! UFS > 0

      end subroutine ca_normalize_new_species

! :::
! ::: ------------------------------------------------------------------
! :::

    AMREX_CUDA_FORT_DEVICE subroutine ca_consup(lo, hi, uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                      hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
                      flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                      flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                      flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                      qx, qx_lo, qx_hi, &
#if AMREX_SPACEDIM >= 2
                        qy, qy_lo, qy_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        qz, qz_lo, qz_hi, &
#endif
                      pdivu, pdivu_lo, pdivu_hi, &
!                      divu_cc,d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                      dx,dt,a_old,a_new) bind(C, name="ca_consup")

      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : NVAR, URHO, UMX, UMZ, &
           UEDEN, UEINT, gamma_minus_1, NGDNV, &
           use_pressure_law_pdivu, use_area_dt_scale_apply
      
      !use advection_util_module, only : calc_pdivu
      
      implicit none

      integer lo(3), hi(3)
      integer   uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer   hsrc_l1,  hsrc_l2,  hsrc_l3,  hsrc_h1,  hsrc_h2,  hsrc_h3
      integer pdivu_lo(3), pdivu_hi(3)
!      integer d_l1,   d_l2,   d_l3,   d_h1,   d_h2,   d_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    
      real(rt)  :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      real(rt)  :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt)  :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      real(rt)  :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      real(rt)  :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
      real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
      real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
      real(rt), intent(inout) :: pdivu(pdivu_lo(1):pdivu_hi(1),pdivu_lo(2):pdivu_hi(2),pdivu_lo(3):pdivu_hi(3))
!      real(rt)  :: divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
      real(rt), intent(in)  :: dx(3)
      real(rt), intent(in), value  :: dt, a_old, a_new

      real(rt) :: div1, a_half, a_oldsq, a_newsq
      real(rt) :: area1, area2, area3
      real(rt) :: vol, volinv, a_newsq_inv
      real(rt) :: a_half_inv, a_new_inv, dt_a_new
      real(rt) :: src_eint
      integer  :: i, j, k, n

      a_half  = HALF * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      vol     = dx(1) * dx(2) * dx(3)
      volinv  = ONE / vol

      a_half_inv  = ONE / a_half
      a_new_inv   = ONE / a_new
      dt_a_new    = dt / a_new
      a_newsq_inv = ONE / a_newsq

      area1   =      dx(2) * dx(3)
      area2   = dx(1)      * dx(3)
      area3   = dx(1) * dx(2)

      call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, &
#if AMREX_SPACEDIM >= 2
                    qy, qy_lo, qy_hi, &
                    area2, &
#endif
#if AMREX_SPACEDIM == 3
                    qz, qz_lo, qz_hi, &
                    area3, &
#endif
                    vol, &
                    dx, pdivu, pdivu_lo, pdivu_hi)
    

      if (use_area_dt_scale_apply .eq. 1) then

         do n = 1, NVAR
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     
                     ! Density
                     if (n .eq. URHO) then
                        hydro_src(i,j,k,n) = &
                             ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                             + flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                             + flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv * a_half_inv
                        
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
                        
                        if (use_pressure_law_pdivu .eq. 0) then
                           src_eint = &
                                +(a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                                -a_half*dt*pdivu(i,j,k))
                        else
                           src_eint = &
                                +(a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                                -a_half*dt*(gamma_minus_1 * uin(i,j,k,n)) * pdivu(i,j,k))
                        end if
                        
                        hydro_src(i,j,k,n) =  &
                             ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                             +flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                             +flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * a_half * volinv & !hydro_src matches old with factor of dt here
                             +src_eint

                        ! (rho X_i) and (rho adv_i) and (rho aux_i)
                     else
                        hydro_src(i,j,k,n) = &
                             ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                             + flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                             + flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv * a_half_inv
                     end if

                     hydro_src(i,j,k,n) = hydro_src(i,j,k,n) / dt

                  enddo
               enddo
            enddo
         enddo

      else
         do n = 1, NVAR
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     
                     ! Density
                     if (n .eq. URHO) then
                        hydro_src(i,j,k,n) = &
                             ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) ) * area1 &
                             + ( flux2(i,j,k,n) - flux2(i,j+1,k,n) ) * area2 &
                             + ( flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * area3) * volinv * a_half_inv
                        
                     ! Momentum
                     else if (n .ge. UMX .and. n .le. UMZ) then
                        hydro_src(i,j,k,n) =  &
                             ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) ) * area1 &
                             + ( flux2(i,j,k,n) - flux2(i,j+1,k,n) ) * area2 &
                             + ( flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * area3) * volinv  

                     ! (rho E)
                     else if (n .eq. UEDEN) then
                        hydro_src(i,j,k,n) =  &
                             ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) ) * area1 &
                             + ( flux2(i,j,k,n) - flux2(i,j+1,k,n) ) * area2 &
                             + ( flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * area3) * a_half * volinv  &
                             +   a_half * (a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT)
                        

                     ! (rho e)
                     else if (n .eq. UEINT) then
                        
                        if (use_pressure_law_pdivu .eq. 0) then
                           src_eint = &
                                +(a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                                -a_half*dt*pdivu(i,j,k)) / dt
                        else
                           src_eint = &
                                +(a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                                -a_half*dt*(gamma_minus_1 * uin(i,j,k,n)) * pdivu(i,j,k)) / dt
                        end if
                        
                        hydro_src(i,j,k,n) =  &
                             ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) ) * area1 &
                             + (  flux2(i,j,k,n) - flux2(i,j+1,k,n) ) * area2 &
                             + (  flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * area3) * a_half * volinv  &
                             +src_eint

                     ! (rho X_i) and (rho adv_i) and (rho aux_i)
                     else
                        hydro_src(i,j,k,n) = &
                             ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) ) * area1 &
                             + ( flux2(i,j,k,n) - flux2(i,j+1,k,n) ) * area2 &
                             + ( flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * area3) * volinv * a_half_inv
                     endif

                  enddo
               enddo
            enddo
         enddo
      endif
      
    end subroutine ca_consup

  !> @brief This is a basic multi-dimensional shock detection algorithm.
  !! This implementation follows Flash, which in turn follows
  !! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
  !!
  !! The spirit of this follows the shock detection in Colella &
  !! Woodward (1984)
  !!
  subroutine ca_shock(lo, hi, &
       q, qd_lo, qd_hi, &
       shk, s_lo, s_hi, &
       dx) bind(C, name="ca_shock")

    use meth_params_module, only : QPRES, QU, QV, QW, QVAR
    use prob_params_module, only : coord_type
    use amrex_constants_module, only: ZERO, HALF, ONE
    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    real(rt), intent(inout) :: shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: i, j, k

    real(rt) :: dxinv, dyinv, dzinv
    real(rt) :: div_u
    real(rt) :: px_pre, px_post, py_pre, py_post, pz_pre, pz_post
    real(rt) :: e_x, e_y, e_z, d
    real(rt) :: p_pre, p_post, pjump

    real(rt), parameter :: small = 1.e-10_rt
    real(rt), parameter :: eps = 0.33e0_rt

    real(rt) :: rm, rc, rp

    !$gpu

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

#ifndef AMREX_USE_CUDA
    if (coord_type /= 0) then
       call amrex_error("ERROR: invalid geometry in shock()")
    endif
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! construct div{U}
             if (coord_type == 0) then
                ! Cartesian
                div_u = HALF*(q(i+1,j,k,QU) - q(i-1,j,k,QU))*dxinv
#if (AMREX_SPACEDIM >= 2)
                div_u = div_u + HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))*dyinv
#endif
#if (AMREX_SPACEDIM == 3)
                div_u = div_u + HALF*(q(i,j,k+1,QW) - q(i,j,k-1,QW))*dzinv
#endif
             elseif (coord_type == 1) then
                ! r-z
                rc = dble(i + HALF)*dx(1)
                rm = dble(i - 1 + HALF)*dx(1)
                rp = dble(i + 1 + HALF)*dx(1)

#if (AMREX_SPACEDIM == 1)
                div_u = HALF*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx(1))
#endif
#if (AMREX_SPACEDIM == 2)
                div_u = HALF*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx(1)) + &
                     HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))/dx(2)
#endif

             elseif (coord_type == 2) then
                ! 1-d spherical
                rc = dble(i + HALF)*dx(1)
                rm = dble(i - 1 + HALF)*dx(1)
                rp = dble(i + 1 + HALF)*dx(1)

                div_u = HALF*(rp**2*q(i+1,j,k,QU) - rm**2*q(i-1,j,k,QU))/(rc**2*dx(1))

#ifndef AMREX_USE_CUDA
             else
                call amrex_error("ERROR: invalid coord_type in shock")
#endif
             endif


             ! find the pre- and post-shock pressures in each direction
             if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < ZERO) then
                px_pre  = q(i+1,j,k,QPRES)
                px_post = q(i-1,j,k,QPRES)
             else
                px_pre  = q(i-1,j,k,QPRES)
                px_post = q(i+1,j,k,QPRES)
             endif

#if (AMREX_SPACEDIM >= 2)
             if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < ZERO) then
                py_pre  = q(i,j+1,k,QPRES)
                py_post = q(i,j-1,k,QPRES)
             else
                py_pre  = q(i,j-1,k,QPRES)
                py_post = q(i,j+1,k,QPRES)
             endif
#else
             py_pre = ZERO
             py_post = ZERO
#endif

#if (AMREX_SPACEDIM == 3)
             if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < ZERO) then
                pz_pre  = q(i,j,k+1,QPRES)
                pz_post = q(i,j,k-1,QPRES)
             else
                pz_pre  = q(i,j,k-1,QPRES)
                pz_post = q(i,j,k+1,QPRES)
             endif
#else
             pz_pre = ZERO
             pz_post = ZERO
#endif

             ! use compression to create unit vectors for the shock direction
             e_x = (q(i+1,j,k,QU) - q(i-1,j,k,QU))**2
#if (AMREX_SPACEDIM >= 2)
             e_y = (q(i,j+1,k,QV) - q(i,j-1,k,QV))**2
#else
             e_y = ZERO
#endif
#if (AMREX_SPACEDIM == 3)
             e_z = (q(i,j,k+1,QW) - q(i,j,k-1,QW))**2
#else
             e_z = ZERO
#endif
             d = ONE/(e_x + e_y + e_z + small)

             e_x = e_x*d
             e_y = e_y*d
             e_z = e_z*d

             ! project the pressures onto the shock direction
             p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre
             p_post = e_x*px_post + e_y*py_post + e_z*pz_post

             ! test for compression + pressure jump to flag a shock
             if (p_pre == ZERO) then
                ! this can happen if U = 0, so e_x, ... = 0
                pjump = ZERO
             else
                pjump = eps - (p_post - p_pre)/p_pre
             endif

             if (pjump < ZERO .and. div_u < ZERO) then
                shk(i,j,k) = ONE
             else
                shk(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine ca_shock

  
  !> @brief this computes the *node-centered* divergence
  !!
  AMREX_CUDA_FORT_DEVICE subroutine calc_pdivu(lo, hi, &
       q1, q1_lo, q1_hi, &
       !       area1, a1_lo, a1_hi, &
       area1, &
#if AMREX_SPACEDIM >= 2
       q2, q2_lo, q2_hi, &
       !       area2, a2_lo, a2_hi, &
              area2, &
#endif
#if AMREX_SPACEDIM == 3
       q3, q3_lo, q3_hi, &
!       area3, a3_lo, a3_hi, &
       area3,  &
#endif
       vol, &
       !vol, v_lo, v_hi, &
       dx, pdivu, div_lo, div_hi)

    use meth_params_module, only : GDPRES, GDU, GDV, GDW, NGDNV, &
         use_pressure_law_pdivu
    use amrex_constants_module, only : HALF
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: pdivu(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))

    integer, intent(in) :: q1_lo(3), q1_hi(3)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt), intent(in) :: area1
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(in) :: area2
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt), intent(in) :: area3
#endif
    real(rt), intent(in) :: vol
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if(use_pressure_law_pdivu.eq.1) then
#if AMREX_SPACEDIM == 1
             pdivu(i,j,k) = (q1(i+1,j,k,GDU)*area1 - q1(i,j,k,GDU)*area1) / vol
#endif

#if AMREX_SPACEDIM == 2
             pdivu(i,j,k) = ((q1(i+1,j,k,GDU)*area1 - q1(i,j,k,GDU)*area1) + &
                  (q2(i,j+1,k,GDV)*area2 - q2(i,j,k,GDV)*area2) ) / vol
#endif

#if AMREX_SPACEDIM == 3
             pdivu(i,j,k) = &
                  (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                  (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)
#endif
          else
             
!!!!!!!!!!!!!!!!!!!!!!This gives different pdivu contribution to hydro_src than gamma
#if AMREX_SPACEDIM == 1
             pdivu(i,j,k) = HALF * &
                  (q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES))* &
                  (q1(i+1,j,k,GDU)*area1 - q1(i,j,k,GDU)*area1) / vol
#endif


#if AMREX_SPACEDIM == 2
             pdivu(i,j,k) = HALF*( &
                  (q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                  (q1(i+1,j,k,GDU)*area1 - q1(i,j,k,GDU)*area1) + &
                  (q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                  (q2(i,j+1,k,GDV)*area2 - q2(i,j,k,GDV)*area2) ) / vol
#endif

#if AMREX_SPACEDIM == 3
             pdivu(i,j,k) = &
                  HALF*(q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                 (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                 HALF*(q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                 (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                 HALF*(q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)
#endif
          endif
          enddo
       enddo
    enddo

  end subroutine calc_pdivu

    AMREX_CUDA_FORT_DEVICE  subroutine scale_flux(lo, hi, &
#if AMREX_SPACEDIM == 1
                          qint, qi_lo, qi_hi, &
#endif
                          flux, f_lo, f_hi, &
                          area, dt, a_old, a_new) bind(c, name="scale_flux")

    use meth_params_module, only: NVAR, UMX, GDPRES, NGDNV, &
           use_area_dt_scale_apply, URHO, UMX, UMZ, &
           UEDEN, UEINT
    use amrex_constants_module, only: HALF, ONE

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
#if AMREX_SPACEDIM == 1
    integer,  intent(in   ) :: qi_lo(3), qi_hi(3)
#endif
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
!    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: dt, area
    real(rt), intent(in   ), value :: a_old, a_new

    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
#if AMREX_SPACEDIM == 1
    real(rt), intent(in   ) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NGDNV)
#endif
    integer :: i, j, k, n

    real(rt) :: a_half, a_newsq
    real(rt) :: a_newsq_inv
    real(rt) :: a_half_inv, a_new_inv

    a_half  = HALF * (a_old + a_new)
    a_newsq = a_new * a_new
    
    a_half_inv  = ONE / a_half
    a_new_inv   = ONE / a_new
    a_newsq_inv = ONE / a_newsq

    !$gpu
    if(use_area_dt_scale_apply.ne.1) then
       do n = 1, NVAR
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   flux(i,j,k,n) = dt * flux(i,j,k,n) * area
                enddo
             enddo
          enddo
       enddo
    endif

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (n .eq. URHO) then
                   flux(i,j,k,n) = flux(i,j,k,n) * (1.d0/a_half)
                else if (n.ge.UMX .and. n.le.UMZ) then
                   flux(i,j,k,n) = flux(i,j,k,n) * a_new_inv
                else if (n.eq.UEINT .or. n.eq.UEDEN) then
                   flux(i,j,k,n) = flux(i,j,k,n) * (a_half * a_new_inv * a_new_inv)
                else
                   flux(i,j,k,n) = flux(i,j,k,n) * (1.d0/a_half)
                end if
             enddo
          enddo
       enddo
    enddo

  end subroutine scale_flux
  
end module advection_module
