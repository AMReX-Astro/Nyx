
! ::: -----------------------------------------------------------

      subroutine hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                         domlo,domhi,delta,xlo,time,bc)
 
      use meth_params_module
 
      implicit none
      include 'AMReX_bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
 
      integer n
 
      do n = 1,NVAR
          call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
               adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
               domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end subroutine hypfill

! ::: -----------------------------------------------------------

      subroutine denfill(den,den_l1,den_l2,den_l3,den_h1,den_h2,den_h3, &
                         domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'AMReX_bc_types.fi'
      integer den_l1,den_l2,den_l3,den_h1,den_h2,den_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision den(den_l1:den_h1,den_l2:den_h2,den_l3:den_h3)

      call filcc(den,den_l1,den_l2,den_l3,den_h1,den_h2,den_h3,domlo,domhi,delta,xlo,bc)

      end subroutine denfill

! ::: -----------------------------------------------------------

      subroutine generic_fill(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'AMReX_bc_types.fi'
      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)

      call filcc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc)

      end subroutine generic_fill
