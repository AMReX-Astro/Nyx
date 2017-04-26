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
