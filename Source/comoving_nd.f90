module comoving_nd_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  contains

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_integrate_comoving_a(old_a,new_a,dt) &
         bind(C, name="fort_integrate_comoving_a")

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, &
                                                comoving_OmR, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a, dt
        real(rt), intent(  out) :: new_a

        real(rt), parameter :: xacc = 1.0d-8
        real(rt) :: H_0, OmL
        real(rt) :: Delta_t, prev_soln
        real(rt) :: start_a, end_a, start_slope, end_slope
        integer  :: iter, j, nsteps

        if (comoving_h .eq. 0.0d0) then
          new_a = old_a
          return
        endif

        H_0 = comoving_h * Hubble_const
        OmL = 1.d0 - comoving_OmM - comoving_OmR 

        prev_soln = 2.0d0 ! 0<a<1 so a=2 will do as "wrong" solution
        do iter = 1, 11  ! max allowed iterations
          nsteps = 2**(iter-1)
          Delta_t = dt/nsteps
          end_a = old_a

          do j = 1, nsteps
            ! This uses RK2 to integrate the ODE:
            !   da / dt = H_0 * sqrt(OmM/a + OmR/a^2 + OmL*a^2)
            start_a = end_a

            ! Compute the slope at the old time
            if (comoving_type > 0) then
                start_slope = H_0*dsqrt(comoving_OmM/start_a + comoving_OmR/(start_a*start_a) + OmL*start_a**2)
            else
                start_slope = comoving_h
            end if

            ! Compute a provisional value of ln(a) at the new time 
            end_a = start_a + start_slope * Delta_t
            
            ! Compute the slope at the new time
            if (comoving_type > 0) then
                end_slope = H_0*dsqrt(comoving_OmM/end_a + comoving_OmR/(end_a*end_a) + OmL*end_a**2)
            else
                end_slope = comoving_h 
            end if
       
            ! Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5d0 * (start_slope + end_slope) * Delta_t
          enddo

          new_a  = end_a
          if (abs(1.0d0-new_a/prev_soln) .le. xacc) return
          prev_soln = new_a

        enddo

      end subroutine fort_integrate_comoving_a

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_integrate_comoving_a_to_z(old_a,z_value,dt) &
         bind(C, name="fort_integrate_comoving_a_to_z")

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, &
                                                comoving_OmR, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a, z_value
        real(rt), intent(inout) :: dt

        real(rt), parameter :: xacc = 1.0d-8
        real(rt) :: H_0, OmL
        real(rt) :: Delta_t
        real(rt) :: start_a, end_a, start_slope, end_slope
        real(rt) :: a_value
        integer  :: j, nsteps

        if (comoving_h .eq. 0.0d0) &
          call amrex_error("fort_integrate_comoving_a_to_z: Shouldn't be setting plot_z_values if not evolving a")

        H_0 = comoving_h * Hubble_const
        OmL = 1.d0 - comoving_OmM - comoving_OmR
      
        ! Translate the target "z" into a target "a"
        a_value = 1.d0 / (1.d0 + z_value)

        ! Use lots of steps if we want to nail the z_value
        nsteps = 1024

        ! We integrate a, but stop when a = a_value (or close enough)
        Delta_t = dt/nsteps
        end_a = old_a
        do j = 1, nsteps
            ! This uses RK2 to integrate the ODE:
            !   da / dt = H_0 * sqrt(OmM/a + OmR/a^2 + OmL*a^2)
            start_a = end_a

            ! Compute the slope at the old time
            if (comoving_type > 0) then
                start_slope = H_0*dsqrt(comoving_OmM/start_a + comoving_OmR/(start_a*start_a) + OmL*start_a**2)
            else
                start_slope = comoving_h
            end if

            ! Compute a provisional value of ln(a) at the new time 
            end_a = start_a + start_slope * Delta_t

            ! Compute the slope at the new time
            if (comoving_type > 0) then
                end_slope = H_0*dsqrt(comoving_OmM/end_a + comoving_OmR/(end_a*end_a) + OmL*end_a**2)
            else
                end_slope = comoving_h 
            end if
       
            ! Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5d0 * (start_slope + end_slope) * Delta_t

            ! We have crossed from a too small to a too big in this step
            if ( (end_a - a_value) * (start_a - a_value) < 0) then
                dt = ( (  end_a - a_value) * dble(j  ) + &
                       (a_value - start_a) * dble(j+1) ) / (end_a - start_a) * Delta_t
                exit
            end if
        end do

      end subroutine fort_integrate_comoving_a_to_z

! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

!       subroutine fort_get_omb(frac) &
!          bind(C, name="fort_get_omb")

!         use comoving_module, only: comoving_OmB, comoving_OmM

!         real(rt) :: frac

!         frac = comoving_OmB / comoving_OmM

!       end subroutine fort_get_omb


! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

!       subroutine fort_get_omm(omm) &
!          bind(C, name="fort_get_omm")

!         use comoving_module, only: comoving_OmM

!         real(rt) :: omm

!         omm = comoving_OmM

!       end subroutine fort_get_omm

! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

        subroutine fort_get_hubble(hubble) &
           bind(C, name="fort_get_hubble")

          use comoving_module, only: comoving_h

          real(rt) :: hubble

          hubble = comoving_h

        end subroutine fort_get_hubble


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_omb(omb) &
         bind(C, name="fort_set_omb")

        use comoving_module, only: comoving_OmB

        real(rt), intent(in) :: omb

        comoving_OmB = omb

      end subroutine fort_set_omb


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_omm(omm) &
         bind(C, name="fort_set_omm")

        use comoving_module, only: comoving_OmM

        real(rt), intent(in) :: omm

        comoving_OmM = omm

      end subroutine fort_set_omm


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_omr(omr) &
         bind(C, name="fort_set_omr")

        use comoving_module, only: comoving_OmR

        real(rt), intent(in) :: omr

        comoving_OmR = omr

      end subroutine fort_set_omr

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_hubble(hubble) &
         bind(C, name="fort_set_hubble")

        use comoving_module, only: comoving_h

        real(rt), intent(in) :: hubble

        comoving_h = hubble

      end subroutine fort_set_hubble

end module comoving_nd_module
