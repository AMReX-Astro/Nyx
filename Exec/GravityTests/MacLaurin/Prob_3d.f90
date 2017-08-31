
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      use network   , only : network_init
      implicit none

      integer init, namlen
      integer name(namlen)
      real(rt) problo(3), probhi(3)

      integer untin,i

      namelist /fortin/ max_num_part, a1, a3

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      center(1:3) = 0.5d0 * (problo(1:3) + probhi(1:3))

      end

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
      subroutine fort_initdata(level,time,lo,hi, &
                               ns, state   ,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                               nd, diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                               delta,xlo,xhi)  &
                               bind(C, name="fort_initdata")

      use amrex_fort_module, only : rt => amrex_real
      use bl_constants_module         , only : M_PI
      use probdata_module             , only : a1,a3,center
      use fundamental_constants_module, only : Gconst
      use meth_params_module          , only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFA, &
                                               TEMP_COMP

      implicit none
      integer level, ns, nd
      integer lo(3), hi(3)
      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt) xlo(3), xhi(3), time, delta(3)
      real(rt)    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      real(rt) diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)

      ! Local variables
      integer          :: i,j,k
      integer          :: ii,jj,kk
      real(rt) :: x,y,z,dist,rsq,zsq
      real(rt) :: xxlo,yylo,zzlo
      real(rt) :: b,lambda,h,ah,e,esq
      real(rt) :: a1sq,a3sq
      real(rt) :: capA1, capA3
      real(rt) :: dx_fine
      real(rt) :: contrib_rho, contrib_ufa

      dx_fine = 0.5d0 * delta(1)

      a1sq = a1**2
      a3sq = a3**2

      esq = 1.d0 - a3sq/a1sq
      e   = dsqrt(esq)

      capA1 = dsqrt(1.d0 - esq) / e**3 * dasin(e) - (1.d0 - esq) / esq
      capA3 = 2.d0/esq  - (2.d0/e**3) * dsqrt(1.d0-esq) * dasin(e)

      do k = lo(3), hi(3)
         zzlo = (dble(k))*delta(3) - center(3)

         do j = lo(2), hi(2)
            yylo = (dble(j))*delta(2) - center(2)

            do i = lo(1), hi(1)

               xxlo = (dble(i))*delta(1) - center(1)

               contrib_rho = 0.d0
               contrib_ufa = 0.d0

               do kk = 0,1
               do jj = 0,1
               do ii = 0,1

                  x = xxlo + (dble(ii)+0.5d0) * dx_fine
                  y = yylo + (dble(jj)+0.5d0) * dx_fine
                  z = zzlo + (dble(kk)+0.5d0) * dx_fine

                  rsq = x**2 + y**2
                  zsq = z**2

                  dist = dsqrt( rsq/a1sq + zsq/a3sq )

                  b = a1sq + a3sq - rsq - zsq
                  lambda = 0.5d0 * ( -b + dsqrt( b**2 - 4.d0 * &
                                   (a1sq * a3sq - rsq * a3sq - zsq * a1sq) ) )
                  h = a1 * e / (dsqrt(a3sq + lambda))
                  ah = datan(h)
   
                  if (dist .lt. 1.d0) then 
                     contrib_rho = contrib_rho + 1.d0
                     contrib_ufa = contrib_ufa + &
                         M_PI * Gconst * ( &
                         2.d0 * capA1 * a1sq - capA1 * rsq + capA3 * (a3sq - zsq) ) 
                  else
                     contrib_rho = contrib_rho + 1.d-16
                     contrib_ufa = contrib_ufa + &
                         2.d0 * M_PI * Gconst * a1 * a3 / e * ( &
                         ah - 0.5d0 / (a1sq * esq) * &
                        ( rsq * (ah - h / (1.d0+h**2) ) + 2.d0 * zsq * (h - ah) ) )
                  endif

               end do
               end do
               end do

               state(i,j,k,URHO) = contrib_rho / 8.d0
               state(i,j,k,UFA ) = contrib_ufa / 8.d0

               state(i,j,k,UMX:UMZ) = 0.0d0

               state(i,j,k,UEINT) = state(i,j,k,URHO)
               state(i,j,k,UEDEN) = state(i,j,k,URHO)

               if (UFS .gt. 0) then
                   state(i,j,k,UFS  ) = state(i,j,k,URHO)
                   state(i,j,k,UFS+1) = 0.d0
               end if

               diag_eos(i,j,k,TEMP_COMP) = 1000.d0

            enddo
         enddo
      enddo

      end subroutine fort_initdata
