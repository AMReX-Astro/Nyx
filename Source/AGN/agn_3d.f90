  subroutine nyx_compute_overlap(np, particles, ng, ghosts, delta_x) &
       bind(c,name='nyx_compute_overlap')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: agn_particle_t

    integer             , intent(in   ) :: np, ng
    type(agn_particle_t), intent(inout) :: particles(np)
    type(agn_particle_t), intent(in   ) :: ghosts(ng)
    real(amrex_real)    , intent(in   ) :: delta_x(3)

    real(amrex_real) dx, dy, dz, r2
    real(amrex_real) cutoff
    integer i, j

    cutoff = delta_x(1)
    
    do i = 1, np
       do j = i+1, np

          dx = particles(i)%pos(1) - particles(j)%pos(1)
          dy = particles(i)%pos(2) - particles(j)%pos(2)
          dz = particles(i)%pos(3) - particles(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else
             ! We only remove the newer (aka of lower mass) particle
             if (particles(i)%mass .lt. particles(j)%mass)  then
                particles(i)%id = -1
             else
                particles(j)%id = -1
             end if
          end if

       end do
    end do

    do i = 1, np
       do j = 1, ng

          dx = particles(i)%pos(1) - ghosts(j)%pos(1)
          dy = particles(i)%pos(2) - ghosts(j)%pos(2)
          dz = particles(i)%pos(3) - ghosts(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else

             ! We only remove a particle if it is both 1) valid 2) newer (aka of lower mass)
             if (particles(i)%mass .lt. ghosts(j)%mass)  &
                particles(i)%id = -1
             
          end if

      end do
    end do

  end subroutine nyx_compute_overlap

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_merge_particles(np, particles, ng, ghosts, delta_x) &
       bind(c,name='agn_merge_particles')

    use iso_c_binding

    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst
    use particle_mod      , only: agn_particle_t

    integer             , intent(in   ) :: np, ng
    type(agn_particle_t), intent(inout) :: particles(np)
    type(agn_particle_t), intent(in   ) :: ghosts(ng)
    real(amrex_real)    , intent(in   ) :: delta_x(3)

    real(amrex_real) dx, dy, dz, r2
    real(amrex_real) du, dv, dw, vrelsq
    real(amrex_real) xmom, ymom, zmom
    real(amrex_real) cutoff, larger_mass
    integer i, j

    cutoff = delta_x(1)
    
    do i = 1, np
       do j = i+1, np

          ! Distance between particles
          dx = particles(i)%pos(1) - particles(j)%pos(1)
          dy = particles(i)%pos(2) - particles(j)%pos(2)
          dz = particles(i)%pos(3) - particles(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else

             ! Relative velocity
             du = particles(i)%vel(1) - particles(j)%vel(1)
             dv = particles(i)%vel(2) - particles(j)%vel(2)
             dw = particles(i)%vel(3) - particles(j)%vel(3)

             vrelsq = du * du + dv * dv + dw * dw

             larger_mass = max(particles(i)%mass, particles(j)%mass)

             if ( (vrelsq * r2) .lt. (Gconst * larger_mass)**2 ) then

                ! Total momentum of both particles before merging
                xmom = particles(i)%vel(1) * particles(i)%mass + particles(j)%vel(1) * particles(j)%mass
                ymom = particles(i)%vel(2) * particles(i)%mass + particles(j)%vel(2) * particles(j)%mass
                zmom = particles(i)%vel(3) * particles(i)%mass + particles(j)%vel(3) * particles(j)%mass
 
                ! Combine masses and put onto particle "i" if that is the heavier particle
                if (particles(i)%mass .ge. particles(j)%mass)  then
                   particles(i)%mass = particles(i)%mass + particles(j)%mass

                   ! Put total momentum onto the one particle 
                   particles(i)%vel(1) = xmom / particles(i)%mass
                   particles(i)%vel(2) = ymom / particles(i)%mass
                   particles(i)%vel(3) = zmom / particles(i)%mass
   
                   ! Set particle ID of particle j to -1
                   particles(j)%id = -1

                else
                   particles(j)%mass = particles(i)%mass + particles(j)%mass

                   ! Put total momentum onto the one particle 
                   particles(j)%vel(1) = xmom / particles(j)%mass
                   particles(j)%vel(2) = ymom / particles(j)%mass
                   particles(j)%vel(3) = zmom / particles(j)%mass
   
                   ! Set particle ID of particle i to -1
                   particles(i)%id = -1
                end if

             end if
          end if

       end do
    end do

    do i = 1, np
       do j = 1, ng

          dx = particles(i)%pos(1) - ghosts(j)%pos(1)
          dy = particles(i)%pos(2) - ghosts(j)%pos(2)
          dz = particles(i)%pos(3) - ghosts(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle

          else

             ! Relative velocity
             du = particles(i)%vel(1) - ghosts(j)%vel(1)
             dv = particles(i)%vel(2) - ghosts(j)%vel(2)
             dw = particles(i)%vel(3) - ghosts(j)%vel(3)
             vrelsq = du * du + dv * dv + dw * dw

             larger_mass = max(particles(i)%mass, ghosts(j)%mass)

             if ( (vrelsq * r2) .lt. (Gconst * larger_mass)**2) then

                ! The bigger particle is in the valid region so we put all the mass onto it
                if ( particles(i)%mass .gt. ghosts(j)%mass ) then

                   ! Total momentum of both particles before merging
                   xmom = particles(i)%vel(1) * particles(i)%mass + ghosts(j)%vel(1) * ghosts(j)%mass
                   ymom = particles(i)%vel(2) * particles(i)%mass + ghosts(j)%vel(2) * ghosts(j)%mass
                   zmom = particles(i)%vel(3) * particles(i)%mass + ghosts(j)%vel(3) * ghosts(j)%mass
   
                   ! Combine masses and put onto particle "i"
                   particles(i)%mass = particles(i)%mass + ghosts(j)%mass
   
                   ! Put total momentum onto the one particle
                   particles(i)%vel(1) = xmom / particles(i)%mass
                   particles(i)%vel(2) = ymom / particles(i)%mass
                   particles(i)%vel(3) = zmom / particles(i)%mass
   
                ! The bigger particle "j" is in the ghost region  -- we will let its grid update that particle
                !  so here we just invalidate particle "i"
                else

                   ! Set particle ID of particle i to -1
                   particles(i)%id = -1

                end if

             end if
          end if


      end do
    end do

  end subroutine agn_merge_particles

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_particle_velocity(np, particles, state_old, sold_lo, sold_hi, state_new, snew_lo, snew_hi, &
                                   dx, add_energy) &
       bind(c,name='agn_particle_velocity')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN
    use particle_mod      , only: agn_particle_t

    integer,              intent(in   )        :: np
    integer,              intent(in   )        :: sold_lo(3), sold_hi(3)
    integer,              intent(in   )        :: snew_lo(3), snew_hi(3)
    integer,              intent(in   )        :: add_energy
    type(agn_particle_t), intent(inout)        :: particles(np)
    real(amrex_real),     intent(in   )        :: state_old &
         (sold_lo(1):sold_hi(1),sold_lo(2):sold_hi(2),sold_lo(3):sold_hi(3),NVAR)
    real(amrex_real),     intent(in   )        :: state_new &
         (snew_lo(1):snew_hi(1),snew_lo(2):snew_hi(2),snew_lo(3):snew_hi(3),NVAR)
    real(amrex_real),     intent(in   )        :: dx(3)

    integer          :: i, j, k, n
    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: mass, momx, momy, momz, E
    
    vol = dx(1) * dx(2) * dx(3)

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       ! momx, momy, momz, E: momentum and total energy.

       momx = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMX) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMX)) * weight) * vol

       momy = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMY) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMY)) * weight) * vol

       momz = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMZ) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMZ)) * weight) * vol

       E = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) - &
                state_old(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN)) * weight) * vol

       mass = particles(n)%mass

       ! Update velocity of particle so as to reduce momentum in the amount
       ! of the difference between old and new state.
       particles(n)%vel(1) = particles(n)%vel(1) - momx / mass
       particles(n)%vel(2) = particles(n)%vel(2) - momy / mass
       particles(n)%vel(3) = particles(n)%vel(3) - momz / mass
       
       ! Update energy if particle isn't brand new
       if (add_energy .gt. 0) then
          particles(n)%energy = particles(n)%energy - E / mass
       endif

    end do

  end subroutine agn_particle_velocity

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_accrete_mass(np, particles, state, density_lost, slo, shi, eps_rad, eps_coupling, dt, dx) &
       bind(c,name='agn_accrete_mass')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst, pi, eddington_const, c_light
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT
    use particle_mod      , only: agn_particle_t

    integer,              intent(in   )        :: np, slo(3), shi(3)
    type(agn_particle_t), intent(inout)        :: particles(np)
    real(amrex_real),     intent(in   )        :: state &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(amrex_real),     intent(inout)        :: density_lost &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real),     intent(in   )        :: eps_rad, eps_coupling, dt, dx(3)

    integer          :: i, j, k, n
    real(amrex_real) :: avg_rho, avg_csq, avg_speedsq
    real(amrex_real) :: c, denom, mass, mdot, m_edd

    real(amrex_real), parameter :: alpha = 10.d0
    real(amrex_real), parameter :: bondi_const = alpha * 4.0d0*pi * Gconst*Gconst
    real(amrex_real), parameter :: max_frac_removed = 0.5
    real(amrex_real) :: speedsq &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real) :: csq &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: mass_change

    ! Remember state holds primitive variables:
    ! velocity, not momentum.
    speedsq = state(:,:,:,UMX)**2 + state(:,:,:,UMY)**2 + state(:,:,:,UMZ)**2

    do i = slo(1), shi(1)
       do j = slo(2), shi(2)
          do k = slo(3), shi(3)
             ! state(i,j,k,UEINT) is internal energy divided by density
             if (state(i,j,k,URHO) .lt. 0.) then
                print *, 'density = ', state(i,j,k,URHO)
                continue
             endif
             if (state(i,j,k,UEINT) .lt. 0.) then
                print *, 'energy = ', state(i,j,k,UEINT)
                continue
             endif
             call nyx_eos_soundspeed(c, state(i,j,k,URHO), state(i,j,k,UEINT))
             csq(i, j, k) = c*c
          end do
       end do
    end do

    vol = dx(1) * dx(2) * dx(3)

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       avg_rho =     sum(  state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * weight)
       avg_speedsq = sum(speedsq(i-1:i+1, j-1:j+1, k-1:k+1) * weight)
       avg_csq =     sum(    csq(i-1:i+1, j-1:j+1, k-1:k+1) * weight)

       denom = (avg_speedsq + avg_csq)**1.5

       mass = particles(n)%mass
 
       ! Bondi accretion 
       mdot = bondi_const * mass * mass * avg_rho / denom

       ! Eddington limit
       m_edd = eddington_const * mass / eps_rad
       print *,'MDOT, MEDD: ', mdot, m_edd
       mdot = min(mdot, m_edd)

       mass_change = mdot * dt

       ! From each stencil cell, we'll reduce density by
       ! mass_change * weight[cell] / vol.
       ! But the most density that will be removed from cell is
       ! max_frac_removed * density[cell].
       ! This corresponds to maximum mass removal of
       ! Hence the most mass that will be removed from cell is
       ! max_frac_removed * density[cell] * vol.
       ! That is the maximum allowable mass_change * weight[cell].
       mass_change = MIN(mass_change, &
            max_frac_removed * &
            MAXVAL(state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * vol / weight))
       
       ! Increase the mass of the particle by mass_change
       particles(n)%mass = particles(n)%mass + mass_change * ( 1.d0 - eps_rad)

       density_lost(i-1:i+1, j-1:j+1, k-1:k+1) = &
            density_lost(i-1:i+1, j-1:j+1, k-1:k+1) + &
            mass_change * weight / vol

       particles(n)%energy = particles(n)%energy + mass_change * &
            c_light**2 * eps_coupling * eps_rad
    end do

  end subroutine agn_accrete_mass

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_release_energy(np, particles, state, slo, shi, T_min, dx) &
       bind(c,name='agn_release_energy')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: k_B, m_proton
    use eos_module
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT
    use particle_mod      , only: agn_particle_t

    integer,              intent(in   )        :: np, slo(3), shi(3)
    type(agn_particle_t), intent(inout)        :: particles(np)
    real(amrex_real),     intent(inout)        :: state &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(amrex_real),     intent(in   )        :: T_min
    real(amrex_real),     intent(in   )        :: dx(3)

    integer          :: i, j, k, n

    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: E_crit_over_rho, avg_rho

    vol = dx(1) * dx(2) * dx(3)

    E_crit_over_rho = 1.5 * k_B * T_min * vol / m_proton

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       avg_rho = sum(state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * weight)

       if (particles(n)%energy > E_crit_over_rho * avg_rho) then

          print *, 'RELEASING ENERGY of particle at ', particles(n)%pos
          print *, 'particle energy: ', particles(n)%energy
          print *, 'threshold: ', E_crit_over_rho * avg_rho
          print *, 'avg_rho = ', avg_rho
          print *, 'k_B = ', k_B
          print *, 'T_min = ', T_min

          state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) = &
          state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) + &
          particles(n)%energy * weight

          state(i-1:i+1, j-1:j+1, k-1:k+1, UEINT) = &
          state(i-1:i+1, j-1:j+1, k-1:k+1, UEINT) + &
          particles(n)%energy * weight

          particles(n)%energy = 0.
               
       endif

    end do

  end subroutine agn_release_energy

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_length_frac(frac, x, dx)

    use amrex_fort_module, only : amrex_real
    use bl_constants_module, only : ZERO, ONE, HALF

    real(amrex_real), intent(out  )  :: frac(-1:1)
    real(amrex_real), intent(in   )  :: x
    real(amrex_real), intent(in   )  :: dx

    integer :: i
    real(amrex_real) :: offset

    i = x / dx
    offset = x / dx - i
    if (offset < HALF) then ! offset in range 0 : 0.5
       frac(-1) = HALF - offset
       frac(0) = HALF + offset
       frac(1) = ZERO
    else ! offset in range 0.5 : 1
       frac(-1) = ZERO
       frac(0) = ONE + HALF - offset
       frac(1) = offset - HALF
    endif
    
  end subroutine get_length_frac

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_weights(weight, pos, dx)

    use amrex_fort_module, only : amrex_real
    use bl_constants_module, only : ZERO, ONE, HALF

    real(amrex_real), intent(out  )  :: weight(-1:1, -1:1, -1:1)
    real(amrex_real), intent(in   )  :: pos(3)
    real(amrex_real), intent(in   )  :: dx(3)

    integer :: d, i, j, k
    real(amrex_real) :: frac(-1:1, 3)

    do d = 1, 3
       call get_length_frac(frac(:, d), pos(d), dx(d))
    end do

    do k = -1, 1
       do j = -1, 1
          do i = -1, 1
             weight(i, j, k) = frac(i, 1) * frac(j, 2) * frac(k, 3)
          end do
       end do
    end do

  end subroutine get_weights
