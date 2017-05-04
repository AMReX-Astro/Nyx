  subroutine nyx_compute_overlap(particles, ns, np, my_id, ghosts, ng, delta_x) &
       bind(c,name='nyx_compute_overlap')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np, ng
    integer,          intent(inout)        :: my_id(np)
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: ghosts(3, ng)
    real(amrex_real), intent(in   )        :: delta_x(3)

    real(amrex_real) dx, dy, dz, r2
    real(amrex_real) cutoff
    integer i, j

    cutoff = delta_x(1)
    
    do i = 1, np
       do j = i+1, np

          dx = particles(1, i) - particles(1, j)
          dy = particles(2, i) - particles(2, j)
          dz = particles(3, i) - particles(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else
             ! set particle id of particle j to -1
             my_id(j) = -1
          end if

       end do
    end do

    do i = 1, np
       do j = 1, ng

          dx = particles(1, i) - ghosts(1, j)
          dy = particles(2, i) - ghosts(2, j)
          dz = particles(3, i) - ghosts(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else

             ! We only remove a particle if it is both 1) valid 2) newer (aka of lower mass)
             if (particles(4,i) .lt. particles(4,j))  &
                my_id(i) = -1
             
          end if

      end do
    end do

  end subroutine nyx_compute_overlap

  subroutine agn_merge_particles(particles, ns, np, my_id, ghosts, ng, delta_x) &
       bind(c,name='agn_merge_particles')

    use iso_c_binding

    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst

    integer,          intent(in   ), value :: ns, np, ng
    integer,          intent(inout)        :: my_id(np)
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: ghosts(3, ng)
    real(amrex_real), intent(in   )        :: delta_x(3)

    real(amrex_real) dx, dy, dz, r2, r, coef
    real(amrex_real) du, dv, dw, vrelsq
    real(amrex_real) xmom, ymom, zmom
    real(amrex_real) cutoff, larger_mass
    integer i, j

    cutoff = delta_x(1)
    
    do i = 1, np
       do j = i+1, np

          ! Distance between particles
          dx = particles(1, i) - particles(1, j)
          dy = particles(2, i) - particles(2, j)
          dz = particles(3, i) - particles(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          else

             ! Relative velocity
             du = particles(5, i) - particles(5, j)
             dv = particles(6, i) - particles(6, j)
             dw = particles(7, i) - particles(7, j)
             vrelsq = du * du + dv * dv + dw * dw

             larger_mass = max(particles(4,i), particles(4,j))

             if ( (vrelsq * r2) .lt. (Gconst * larger_mass)**2 ) then

                ! Total momentum of both particles before merging
                xmom = particles(5,i) * particles(4,i) + particles(5,j) * particles(4,j)
                ymom = particles(6,i) * particles(6,i) + particles(6,j) * particles(6,j)
                zmom = particles(7,i) * particles(7,i) + particles(7,j) * particles(7,j)
 
                ! Combine masses and put onto particle "i"
                particles(4,i) = particles(4,i) + particles(4,j)

                ! Put total momentum onto the one particle 
                particles(5,i) = xmom / particles(4,i)
                particles(6,i) = ymom / particles(4,i)
                particles(7,i) = zmom / particles(4,i)

                ! Set particle ID of particle j to -1
                my_id(j) = -1

             end if
          end if

       end do
    end do

    do i = 1, np
       do j = 1, ng

          dx = particles(1, i) - ghosts(1, j)
          dy = particles(2, i) - ghosts(2, j)
          dz = particles(3, i) - ghosts(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle

          else

             ! Relative velocity
             du = particles(5, i) - particles(5, j)
             dv = particles(6, i) - particles(6, j)
             dw = particles(7, i) - particles(7, j)
             vrelsq = du * du + dv * dv + dw * dw

             larger_mass = max(particles(4,i), particles(4,j))

             if ( (vrelsq * r2) .lt. (Gconst * larger_mass)**2) then

                ! The bigger particle is in the valid region so we put all the mass onto it
                if ( particles(i,4) .gt. particles(j,4) ) then

                   ! Total momentum of both particles before merging
                   xmom = particles(5,i) * particles(4,i) + particles(5,j) * particles(4,j)
                   ymom = particles(6,i) * particles(6,i) + particles(6,j) * particles(6,j)
                   zmom = particles(7,i) * particles(7,i) + particles(7,j) * particles(7,j)
   
                   ! Combine masses and put onto particle "i"
                   particles(4,i) = particles(4,i) + particles(4,j)
   
                   ! Put total momentum onto the one particle
                   particles(5,i) = xmom / particles(4,i)
                   particles(6,i) = ymom / particles(4,i)
                   particles(7,i) = zmom / particles(4,i)
   
                   ! Set particle ID of particle j to -1 (we don't really need to since this is a ghost particle we will throw away)
                   my_id(j) = -1
   
                ! The bigger particle "j" is in the ghost region  -- we will let its grid update that particle
                !  so here we just invalidate particle "i"
                else

                   ! Set particle ID of particle i to -1
                   my_id(i) = -1

                end if

             end if
          end if


      end do
    end do

  end subroutine agn_merge_particles

  subroutine agn_particle_velocity(particles, ns, np, state_old, sold_lo, sold_hi, state_new, snew_lo, snew_hi, &
                                   dx, add_energy) &
       bind(c,name='agn_particle_velocity')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN

    integer,          intent(in   ), value :: ns, np
    integer,          intent(in   )        :: sold_lo(3), sold_hi(3)
    integer,          intent(in   )        :: snew_lo(3), snew_hi(3)
    integer,          intent(in   )        :: add_energy
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: state_old &
         (sold_lo(1):sold_hi(1),sold_lo(2):sold_hi(2),sold_lo(3):sold_hi(3),NVAR)
    real(amrex_real), intent(in   )        :: state_new &
         (snew_lo(1):snew_hi(1),snew_lo(2):snew_hi(2),snew_lo(3):snew_hi(3),NVAR)
    real(amrex_real), intent(in   )        :: dx(3)

    integer          :: i, j, k, n
    integer          :: ii, jj, kk
    real(amrex_real) :: sum_ux, sum_uy, sum_uz, sum_rE, dv
    
    dv = dx(1) * dx(2) * dx(3)

    do n = 1, np
 
          ! Components 1,2,3 are (x,y,z)
          i = particles(1,n) / dx(1)
          j = particles(2,n) / dx(2)
          k = particles(3,n) / dx(3)

          ! Set to zero
          sum_ux = 0.d0
          sum_uy = 0.d0
          sum_uz = 0.d0
          sum_rE = 0.d0

          ! Sum over 3^3 cube of cells with (i,j,k) at the center
          do kk = k-1,k+1
          do jj = j-1,j+1
          do ii = i-1,i+1
             sum_ux = sum_ux + (state_new(ii,jj,kk,UMX) - state_old(ii,jj,kk,UMX))
             sum_uy = sum_uy + (state_new(ii,jj,kk,UMY) - state_old(ii,jj,kk,UMY))
             sum_uz = sum_uz + (state_new(ii,jj,kk,UMZ) - state_old(ii,jj,kk,UMZ))
             sum_rE = sum_rE + (state_new(ii,jj,kk,UEDEN) - state_old(ii,jj,kk,UEDEN))
          end do
          end do
          end do

          ! sum_ux, sum_uy, sum_uz contain components of momentum density (rho u) and total energy density (rho E)
          ! added up over cells.  Multiply them by cell volume to get momentum (m u) and total energy (m E).
          sum_ux = sum_ux * dv
          sum_uy = sum_uy * dv
          sum_uz = sum_uz * dv
          sum_rE = sum_rE * dv

          ! Components 4 is mass; components 5,6,7 are (u,v,w)
          particles(5,n) = particles(5,n) - sum_ux / particles(4,n)
          particles(6,n) = particles(6,n) - sum_uy / particles(4,n)
          particles(7,n) = particles(7,n) - sum_uz / particles(4,n)

          ! Component 8 is energy
          if (add_energy .gt. 0) &
             particles(8,n) = particles(8,n) - sum_rE / particles(4,n)

    end do

  end subroutine agn_particle_velocity

  subroutine agn_accrete_mass(particles, ns, np, state, slo, shi, eps_rad, dt, dx) &
       bind(c,name='agn_accrete_mass')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst, pi, m_proton, c_light
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT

    integer,          intent(in   ), value :: ns, np
    integer,          intent(in   )        :: slo(3), shi(3)
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: state &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(amrex_real), intent(in   )        :: eps_rad, dt, dx(3)

    integer          :: i, j, k, n
    integer          :: ii, jj, kk
    real(amrex_real) :: avg_rho, avg_usq, avg_vsq, avg_wsq, avg_csq
    real(amrex_real) :: c, e, denom, mass, mdot, m_edd

    real(amrex_real) :: fourpi, alpha, sigma

    fourpi = 4.d0 * pi

    alpha = 10.d0

    ! Thomson cross-section -- this is totally made up here
    sigma = 1.d0
    
    do n = 1, np

          mass = particles(4,n)
 
          ! Components 1,2,3 are (x,y,z)
          i = particles(1,n) / dx(1)
          j = particles(2,n) / dx(2)
          k = particles(3,n) / dx(3)

          ! Set to zero for each (i,j,k)
          avg_rho = 0.d0
          avg_usq = 0.d0
          avg_vsq = 0.d0
          avg_wsq = 0.d0
          avg_csq = 0.d0

          ! Sum over 3^3 cube of cells with (i,j,k) at the center
          do kk = k-1,k+1
          do jj = j-1,j+1
          do ii = i-1,i+1
             avg_rho = avg_rho + state(ii,jj,kk,URHO)
             avg_usq = avg_usq + (state(ii,jj,kk,UMX) / state(ii,jj,kk,URHO))**2
             avg_vsq = avg_vsq + (state(ii,jj,kk,UMY) / state(ii,jj,kk,URHO))**2
             avg_wsq = avg_wsq + (state(ii,jj,kk,UMZ) / state(ii,jj,kk,URHO))**2

             e = state(ii,jj,kk,UEINT) / state(ii,jj,kk,URHO)
             call nyx_eos_soundspeed(c,state(ii,jj,kk,URHO),e)
             avg_csq  = avg_csq + c*c

          end do
          end do
          end do 
          avg_rho = avg_rho / 27.d0
          avg_usq = avg_usq / 27.d0
          avg_vsq = avg_vsq / 27.d0
          avg_wsq = avg_wsq / 27.d0
          avg_csq = avg_csq / 27.d0

          denom = (avg_usq + avg_vsq + avg_wsq + avg_csq)**1.5

          ! Bondi accretion 
          mdot = alpha * fourpi * Gconst * Gconst * avg_rho * mass * mass * avg_rho / denom

          ! Eddington limit
          m_edd = fourpi * Gconst * mass * m_proton / (eps_rad * sigma * c_light)

          mdot = min(mdot, m_edd)

          ! Increase the mass of the particle by Mdot * dt
          particles(4,n) = particles(4,n) + mdot * dt * ( 1.d0 - eps_rad)

    end do

  end subroutine agn_accrete_mass
