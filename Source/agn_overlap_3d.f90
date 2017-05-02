! subroutine nyx_compute_overlap(particles, ns, np, my_id, ghosts, ng, delta_x, density, dlo, dhi) &
  subroutine nyx_compute_overlap(particles, ns, np, my_id, ghosts, ng, delta_x) &
       bind(c,name='nyx_compute_overlap')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np, ng
!   integer,          intent(in   )        :: dlo(3), dhi(3)
    integer,          intent(inout)        :: my_id(np)
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: ghosts(3, ng)
    real(amrex_real), intent(in   )        :: delta_x(3)
!   real(amrex_real), intent(inout)        :: density &
!        (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    real(amrex_real) dx, dy, dz, r2, r, coef
    real(amrex_real) cutoff, min_r, mass
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

!   do i = 1, np
!      do j = 1, ng

!         dx = particles(1, i) - ghosts(1, j)
!         dy = particles(2, i) - ghosts(2, j)
!         dz = particles(3, i) - ghosts(3, j)

!         r2 = dx * dx + dy * dy + dz * dz

!         if (r2 .gt. cutoff*cutoff) then
!            cycle
!         else
!            ! set particle id of particle j to -1
!         end if

!     end do
!   end do

  end subroutine nyx_compute_overlap

  subroutine agn_particle_velocity(particles, ns, np, state_old, sold_lo, sold_hi, state_new, snew_lo, snew_hi, start_comp, dx) &
       bind(c,name='agn_particle_velocity')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np
    integer,          intent(in   )        :: sold_lo(3), sold_hi(3)
    integer,          intent(in   )        :: snew_lo(3), snew_hi(3)
    integer,          intent(in   )        :: start_comp
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: state_old &
         (sold_lo(1):sold_hi(1),sold_lo(2):sold_hi(2),sold_lo(3):sold_hi(3),0:3)
    real(amrex_real), intent(in   )        :: state_new &
         (snew_lo(1):snew_hi(1),snew_lo(2):snew_hi(2),snew_lo(3):snew_hi(3),0:3)
    real(amrex_real), intent(in   )        :: dx(3)

    integer          :: i, j, k, n
    integer          :: ii, jj, kk
    real(amrex_real) :: sum_ux, sum_uy, sum_uz
    
    do n = 1, np
 
          ! Components 1,2,3 are (x,y,z)
          i = particles(1,n) / dx(1)
          j = particles(2,n) / dx(2)
          k = particles(3,n) / dx(3)

          ! Set to zero
          sum_ux = 0.d0
          sum_uy = 0.d0
          sum_uz = 0.d0

          ! Sum over 3^3 cube of cells with (i,j,k) at the center
          do kk = k-1,k+1
          do jj = j-1,j+1
          do ii = i-1,i+1
             sum_ux = sum_ux + (state_new(ii,jj,kk,start_comp  ) - state_old(ii,jj,kk,start_comp  ))
             sum_uy = sum_uy + (state_new(ii,jj,kk,start_comp+1) - state_old(ii,jj,kk,start_comp+1))
             sum_uz = sum_uz + (state_new(ii,jj,kk,start_comp+2) - state_old(ii,jj,kk,start_comp+2))
          end do
          end do
          end do

          ! Components 4 is mass; components 5,6,7 are (u,v,w)
          particles(5,n) = particles(5,n) - sum_ux / particles(4,n)
          particles(6,n) = particles(6,n) - sum_uy / particles(4,n)
          particles(7,n) = particles(7,n) - sum_uz / particles(4,n)

          print *,'PARTICLE XVEL ', particles(5,n), state_new(i,j,k,start_comp), state_old(i,j,k,start_comp)

    end do

  end subroutine agn_particle_velocity
