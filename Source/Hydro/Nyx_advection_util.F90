module advection_module
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
    use amrex_fort_module, only : rt => amrex_real
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
                                     gamma_const, gamma_minus_1, use_flattening
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
    integer          :: n, iq, ipassive
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

end module advection_module
