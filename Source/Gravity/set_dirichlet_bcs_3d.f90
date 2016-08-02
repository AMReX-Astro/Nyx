
     subroutine fort_set_homog_bcs(lo,hi,domlo,domhi, &
                                   phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,dx);
 
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3),domlo(3),domhi(3)
     integer         ,intent(in   ) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     double precision,intent(  out) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
     double precision,intent(in   ) :: dx(3)

     phi = 0.d0
 
     end subroutine fort_set_homog_bcs

     ! ************************************************************************************
     ! Loops over the particles and adds their monopole contribution to phi outside the boundary
     ! ************************************************************************************

     subroutine fort_add_monopole_bcs(lo,hi,domlo,domhi, &
                                      np,part_locs,part_mass, &
                                      phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,dx);

     use fundamental_constants_module, only : Gconst
 
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3),domlo(3),domhi(3),np
     integer         ,intent(in   ) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     double precision,intent(in   ) :: part_locs(0:3*np-1)
     double precision,intent(in   ) :: part_mass(0:  np-1)
     double precision,intent(  out) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
     double precision,intent(in   ) :: dx(3)

     ! Local variables
     integer          :: i,j,k,n
     double precision :: x,y,z,r
     double precision :: x0,y0,z0

     phi = 0.d0

     ! Define phi = -G M / r where r is distance from particle and M is mass of particle

     do k = lo(3)-1,hi(3)+1
         do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1

                   if (i.lt.domlo(1) .or. i.gt.domhi(1) .or. &
                       j.lt.domlo(2) .or. j.gt.domhi(2) .or. &
                       k.lt.domlo(3) .or. k.gt.domhi(3)) then

                       if (i.lt.domlo(1)) then
                           x0 = (dble(i+1))*dx(1)
                       else if (i.gt.domhi(1)) then
                           x0 = (dble(i))*dx(1)
                       else 
                           x0 = (dble(i)+0.5d0)*dx(1)
                       end if
    
                       if (j.lt.domlo(2)) then
                           y0 = (dble(j+1))*dx(2)
                       else if (j.gt.domhi(2)) then
                           y0 = (dble(j))*dx(2)
                       else 
                           y0 = (dble(j)+0.5d0)*dx(2)
                       end if
    
                       if (k.lt.domlo(3)) then
                           z0 = (dble(k+1))*dx(3)
                       else if (k.gt.domhi(3)) then
                           z0 = (dble(k))*dx(3)
                       else 
                           z0 = (dble(k)+0.5d0)*dx(3)
                       end if

                       do n = 0, np-1
                           x = x0 - part_locs(3*n  )
                           y = y0 - part_locs(3*n+1)
                           z = z0 - part_locs(3*n+2)
                           r = sqrt(x*x + y*y + z*z)
                           phi(i,j,k) = phi(i,j,k) - Gconst * part_mass(n) / r 
                       end do
 
                   end if
              end do
         end do
     end do
 
     end subroutine fort_add_monopole_bcs
