
     subroutine fort_set_dirichlet_bcs(lo,hi,domlo,domhi, &
                           phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,dx);

     use amrex_fort_module, only : rt => amrex_real
     use bl_constants_module         , only : M_PI
     use probdata_module             , only : a1,a3,center
     use meth_params_module          , only : NVAR
     use fundamental_constants_module, only : Gconst
 
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3),domlo(3),domhi(3)
     integer         ,intent(in   ) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     real(rt),intent(  out) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
     real(rt),intent(in   ) :: dx(3)

     integer          :: i,j,k
     real(rt) :: x,y,z,rsq,zsq
     real(rt) :: a1sq,a3sq
     real(rt) :: e,esq,b,lambda,h,ah

     a1sq = a1**2
     a3sq = a3**2

     esq = 1.d0 - a3sq/a1sq
      e   = dsqrt(esq)

     if (lo(1).eq.domlo(1)) then
        i = lo(1)-1
        x = 0.d0 - center(1)
        do k = lo(3), hi(3)
           z = (dble(k)+0.5d0)*dx(3) - center(3)
           zsq = z**2
              do j = lo(2), hi(2)
                 y = (dble(j)+0.5d0)*dx(2) - center(2)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )

           end do
        end do
     end if

     if (hi(1).eq.domhi(1)) then
        i = hi(1)+1
        x = 1.d0 - center(1)
        do k = lo(3), hi(3)
           z = (dble(k)+0.5d0)*dx(3) - center(3)
           zsq = z**2
              do j = lo(2), hi(2)
                 y = (dble(j)+0.5d0)*dx(2) - center(2)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
           end do
        end do
     end if

     if (lo(2).eq.domlo(2)) then
        j = lo(2)-1
        y = 0.d0 - center(2)
        do k = lo(3), hi(3)
           z = (dble(k)+0.5d0)*dx(3) - center(3)
           zsq = z**2
              do i = lo(1), hi(1)
                 x = (dble(i)+0.5d0)*dx(1) - center(1)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
           end do
        end do
     end if

     if (hi(2).eq.domhi(2)) then
        j = hi(2)+1
        y = 1.d0 - center(2)
        do k = lo(3), hi(3) 
           z = (dble(k)+0.5d0)*dx(3) - center(3)
           zsq = z**2
              do i = lo(1), hi(1)
                 x = (dble(i)+0.5d0)*dx(1) - center(1)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
           end do
        end do
     end if

     if (lo(3).eq.domlo(3)) then
        k = lo(3)-1
        z = 0.d0 - center(3)
        zsq = z**2
        do j = lo(2), hi(2)
           y = (dble(j)+0.5d0)*dx(2) - center(2)
              do i = lo(1), hi(1)
                 x = (dble(i)+0.5d0)*dx(1) - center(1)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
           end do
        end do
     end if

     if (hi(3).eq.domhi(3)) then
        k = hi(3)+1
        z = 1.d0 - center(3)
        zsq = z**2
        do j = lo(2), hi(2)
           y = (dble(j)+0.5d0)*dx(2) - center(2)
              do i = lo(1), hi(1)
                 x = (dble(i)+0.5d0)*dx(1) - center(1)
                    rsq = x**2 + y**2
                    b = a1sq + a3sq - rsq - zsq
                    lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                     (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                    h = a1 * e / (dsqrt(a3sq + lambda))
                    ah = datan(h)
                    phi(i,j,k) = &
                       2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                       ah - 0.5d0 / (a1sq * esq) * &
                       ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
           end do
        end do
     end if
 
     end subroutine fort_set_dirichlet_bcs

