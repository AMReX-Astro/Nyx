
      subroutine derforcex(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the x-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use forcing_spect_module
      use amrex_constants_module, only : TWO, ONE, HALF, ZERO, M_PI, M_SQRT_2
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer :: i, j, k
      integer :: m, mi, mj, mk
      integer :: alloc
      real(rt) :: accel

      integer :: num_phases(3)
      real(rt) :: delta_phase(3), phase_lo(3)
      real(rt) :: buf(num_modes) 
      real(rt) :: phasefct_init_even(num_modes), phasefct_init_odd(num_modes)
      real(rt) :: phasefct_mult_even(num_modes,3), phasefct_mult_odd(num_modes,3)
      real(rt), allocatable :: phasefct_even_x(:), phasefct_even_y(:), phasefct_even_z(:) 
      real(rt), allocatable :: phasefct_odd_x(:), phasefct_odd_y(:), phasefct_odd_z(:) 

      delta_phase(:) = TWO*M_PI * delta(:) / (prob_hi(:) - prob_lo(:)) ! phase increment per cell
      phase_lo(:) = (dble(lo(:)) + HALF) * delta_phase(:)              ! phase of low corner
 
      ! compute initial phase factors and multiplying factors
      do m = 1, num_modes
         i = wavevectors(1,m)
         j = wavevectors(2,m)
         k = wavevectors(3,m)

         phasefct_init_even(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * cos(k*phase_lo(3)) - &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * sin(k*phase_lo(3))

         phasefct_init_odd(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * sin(k*phase_lo(3)) + &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * cos(k*phase_lo(3))

         phasefct_mult_even(m,1) = cos(i*delta_phase(1));
         phasefct_mult_odd (m,1) = sin(i*delta_phase(1));

         phasefct_mult_even(m,2) = cos(j*delta_phase(2));
         phasefct_mult_odd (m,2) = sin(j*delta_phase(2));

         phasefct_mult_even(m,3) = cos(k*delta_phase(3));
         phasefct_mult_odd (m,3) = sin(k*delta_phase(3));
      end do

      num_phases(:) = (hi(:)-lo(:)+1)*num_modes

      allocate(phasefct_even_x(num_phases(1)), phasefct_even_y(num_phases(2)), phasefct_even_z(num_phases(3)), &
               phasefct_odd_x(num_phases(1)),  phasefct_odd_y(num_phases(2)),  phasefct_odd_z(num_phases(3)), &
               STAT=alloc)

      if (alloc > 0) call bl_abort('failed to allocate arrays for phase factors')      
 
      ! initialize phase factors for each coordinate axis
      do m = 1, num_modes
         phasefct_even_x(m) = ONE 
         phasefct_odd_x(m)  = ZERO
         phasefct_even_y(m) = ONE 
         phasefct_odd_y(m)  = ZERO 
         phasefct_even_z(m) = phasefct_init_even(m)  
         phasefct_odd_z(m)  = phasefct_init_odd(m)
      end do

      do i = lo(1)+1,hi(1)
         mi = (i-lo(1))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_x(mi-num_modes);
              phasefct_even_x(mi) = phasefct_mult_even(m,1) * phasefct_even_x(mi-num_modes) - &
                                    phasefct_mult_odd (m,1) * phasefct_odd_x(mi-num_modes)
              phasefct_odd_x(mi)  = phasefct_mult_even(m,1) * phasefct_odd_x(mi-num_modes) + &
                                    phasefct_mult_odd (m,1) * buf(m)
              mi = mi + 1
         end do
      end do         

      do j = lo(2)+1,hi(2)
         mj = (j-lo(2))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_y(mj-num_modes);
              phasefct_even_y(mj) = phasefct_mult_even(m,2) * phasefct_even_y(mj-num_modes) - &
                                    phasefct_mult_odd (m,2) * phasefct_odd_y(mj-num_modes)
              phasefct_odd_y(mj)  = phasefct_mult_even(m,2) * phasefct_odd_y(mj-num_modes) + &
                                    phasefct_mult_odd (m,2) * buf(m)
              mj = mj + 1
         end do
      end do         

      do k = lo(3)+1, hi(3)
         mk = (k-lo(3))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_z(mk-num_modes);
              phasefct_even_z(mk) = phasefct_mult_even(m,3) * phasefct_even_z(mk-num_modes) - &
                                    phasefct_mult_odd (m,3) * phasefct_odd_z(mk-num_modes)
              phasefct_odd_z(mk)  = phasefct_mult_even(m,3) * phasefct_odd_z(mk-num_modes) + &
                                    phasefct_mult_odd (m,3) * buf(m)
              mk = mk + 1
         end do
      end do

      ! compute acceleration component in physical space
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               accel = ZERO

               mi = (i-lo(1))*num_modes + 1 ! offset in x-direction
               mj = (j-lo(2))*num_modes + 1 ! offset in y-direction
               mk = (k-lo(3))*num_modes + 1 ! offset in z-direction
  
               do m = 1, num_modes
                  ! sum up even modes
                  accel = accel + &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_even_z(mk) - &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_odd_z(mk))  * modes_even(m,1)
                  ! sum up odd modes
                  accel = accel - &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_odd_z(mk)  + &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_even_z(mk)) * modes_odd(m,1)

                  mi = mi + 1
                  mj = mj + 1
                  mk = mk + 1
               end do

               force(i,j,k,1) = M_SQRT_2 * dat(i,j,k,1) * accel

            end do
         end do
      end do

      end subroutine derforcex

!-----------------------------------------------------------------------

      subroutine derforcey(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the y-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use forcing_spect_module
      use amrex_constants_module, only : TWO, ONE, HALF, ZERO, M_PI, M_SQRT_2
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer :: i, j, k
      integer :: m, mi, mj, mk
      integer :: alloc
      real(rt) :: accel

      integer :: num_phases(3)
      real(rt) :: delta_phase(3), phase_lo(3)
      real(rt) :: buf(num_modes) 
      real(rt) :: phasefct_init_even(num_modes), phasefct_init_odd(num_modes)
      real(rt) :: phasefct_mult_even(num_modes,3), phasefct_mult_odd(num_modes,3)
      real(rt), allocatable :: phasefct_even_x(:), phasefct_even_y(:), phasefct_even_z(:) 
      real(rt), allocatable :: phasefct_odd_x(:), phasefct_odd_y(:), phasefct_odd_z(:) 

      delta_phase(:) = TWO*M_PI * delta(:) / (prob_hi(:) - prob_lo(:)) ! phase increment per cell
      phase_lo(:) = (dble(lo(:)) + HALF) * delta_phase(:)              ! phase of low corner
 
      ! compute initial phase factors and multiplying factors
      do m = 1, num_modes
         i = wavevectors(1,m)
         j = wavevectors(2,m)
         k = wavevectors(3,m)

         phasefct_init_even(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * cos(k*phase_lo(3)) - &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * sin(k*phase_lo(3))

         phasefct_init_odd(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * sin(k*phase_lo(3)) + &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * cos(k*phase_lo(3))

         phasefct_mult_even(m,1) = cos(i*delta_phase(1));
         phasefct_mult_odd (m,1) = sin(i*delta_phase(1));

         phasefct_mult_even(m,2) = cos(j*delta_phase(2));
         phasefct_mult_odd (m,2) = sin(j*delta_phase(2));

         phasefct_mult_even(m,3) = cos(k*delta_phase(3));
         phasefct_mult_odd (m,3) = sin(k*delta_phase(3));
      end do

      num_phases(:) = (hi(:)-lo(:)+1)*num_modes

      allocate(phasefct_even_x(num_phases(1)), phasefct_even_y(num_phases(2)), phasefct_even_z(num_phases(3)), &
               phasefct_odd_x(num_phases(1)),  phasefct_odd_y(num_phases(2)),  phasefct_odd_z(num_phases(3)), &
               STAT=alloc)

      if (alloc > 0) call bl_abort('failed to allocate arrays for phase factors')      
 
      ! initialize phase factors for each coordinate axis
      do m = 1, num_modes
         phasefct_even_x(m) = ONE 
         phasefct_odd_x(m)  = ZERO
         phasefct_even_y(m) = ONE 
         phasefct_odd_y(m)  = ZERO 
         phasefct_even_z(m) = phasefct_init_even(m)  
         phasefct_odd_z(m)  = phasefct_init_odd(m)
      end do

      do i = lo(1)+1,hi(1)
         mi = (i-lo(1))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_x(mi-num_modes);
              phasefct_even_x(mi) = phasefct_mult_even(m,1) * phasefct_even_x(mi-num_modes) - &
                                    phasefct_mult_odd (m,1) * phasefct_odd_x(mi-num_modes)
              phasefct_odd_x(mi)  = phasefct_mult_even(m,1) * phasefct_odd_x(mi-num_modes) + &
                                    phasefct_mult_odd (m,1) * buf(m)
              mi = mi + 1
         end do
      end do         

      do j = lo(2)+1,hi(2)
         mj = (j-lo(2))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_y(mj-num_modes);
              phasefct_even_y(mj) = phasefct_mult_even(m,2) * phasefct_even_y(mj-num_modes) - &
                                    phasefct_mult_odd (m,2) * phasefct_odd_y(mj-num_modes)
              phasefct_odd_y(mj)  = phasefct_mult_even(m,2) * phasefct_odd_y(mj-num_modes) + &
                                    phasefct_mult_odd (m,2) * buf(m)
              mj = mj + 1
         end do
      end do         

      do k = lo(3)+1, hi(3)
         mk = (k-lo(3))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_z(mk-num_modes);
              phasefct_even_z(mk) = phasefct_mult_even(m,3) * phasefct_even_z(mk-num_modes) - &
                                    phasefct_mult_odd (m,3) * phasefct_odd_z(mk-num_modes)
              phasefct_odd_z(mk)  = phasefct_mult_even(m,3) * phasefct_odd_z(mk-num_modes) + &
                                    phasefct_mult_odd (m,3) * buf(m)
              mk = mk + 1
         end do
      end do

      ! compute acceleration component in physical space
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               accel = ZERO

               mi = (i-lo(1))*num_modes + 1 ! offset in x-direction
               mj = (j-lo(2))*num_modes + 1 ! offset in y-direction
               mk = (k-lo(3))*num_modes + 1 ! offset in z-direction
  
               do m = 1, num_modes
                  ! sum up even modes
                  accel = accel + &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_even_z(mk) - &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_odd_z(mk))  * modes_even(m,2)
                  ! sum up odd modes
                  accel = accel - &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_odd_z(mk)  + &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_even_z(mk)) * modes_odd(m,2)

                  mi = mi + 1
                  mj = mj + 1
                  mk = mk + 1
               end do

               force(i,j,k,1) = M_SQRT_2 * dat(i,j,k,1) * accel

            end do
         end do
      end do

      end subroutine derforcey

!-----------------------------------------------------------------------

      subroutine derforcez(force,force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine computes the z-component of the forcing term
!
      use amrex_fort_module, only : rt => amrex_real
      use forcing_spect_module
      use amrex_constants_module, only : TWO, ONE, HALF, ZERO, M_PI, M_SQRT_2
      use probdata_module    , only : prob_lo, prob_hi

      implicit none

      integer          lo(3), hi(3)
      integer          force_l1,force_l2,force_l3,force_h1,force_h2,force_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) force(force_l1:force_h1,force_l2:force_h2,force_l3:force_h3,nv)
      real(rt) dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer :: i, j, k
      integer :: m, mi, mj, mk
      integer :: alloc
      real(rt) :: accel

      integer :: num_phases(3)
      real(rt) :: delta_phase(3), phase_lo(3)
      real(rt) :: buf(num_modes) 
      real(rt) :: phasefct_init_even(num_modes), phasefct_init_odd(num_modes)
      real(rt) :: phasefct_mult_even(num_modes,3), phasefct_mult_odd(num_modes,3)
      real(rt), allocatable :: phasefct_even_x(:), phasefct_even_y(:), phasefct_even_z(:) 
      real(rt), allocatable :: phasefct_odd_x(:), phasefct_odd_y(:), phasefct_odd_z(:) 

      delta_phase(:) = TWO*M_PI * delta(:) / (prob_hi(:) - prob_lo(:)) ! phase increment per cell
      phase_lo(:) = (dble(lo(:)) + HALF) * delta_phase(:)              ! phase of low corner
 
      ! compute initial phase factors and multiplying factors
      do m = 1, num_modes
         i = wavevectors(1,m)
         j = wavevectors(2,m)
         k = wavevectors(3,m)

         phasefct_init_even(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * cos(k*phase_lo(3)) - &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * sin(k*phase_lo(3))

         phasefct_init_odd(m) = &
            (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
             sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * sin(k*phase_lo(3)) + &
            (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
             sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * cos(k*phase_lo(3))

         phasefct_mult_even(m,1) = cos(i*delta_phase(1));
         phasefct_mult_odd (m,1) = sin(i*delta_phase(1));

         phasefct_mult_even(m,2) = cos(j*delta_phase(2));
         phasefct_mult_odd (m,2) = sin(j*delta_phase(2));

         phasefct_mult_even(m,3) = cos(k*delta_phase(3));
         phasefct_mult_odd (m,3) = sin(k*delta_phase(3));
      end do

      num_phases(:) = (hi(:)-lo(:)+1)*num_modes

      allocate(phasefct_even_x(num_phases(1)), phasefct_even_y(num_phases(2)), phasefct_even_z(num_phases(3)), &
               phasefct_odd_x(num_phases(1)),  phasefct_odd_y(num_phases(2)),  phasefct_odd_z(num_phases(3)), &
               STAT=alloc)

      if (alloc > 0) call bl_abort('failed to allocate arrays for phase factors')      
 
      ! initialize phase factors for each coordinate axis
      do m = 1, num_modes
         phasefct_even_x(m) = ONE 
         phasefct_odd_x(m)  = ZERO
         phasefct_even_y(m) = ONE 
         phasefct_odd_y(m)  = ZERO 
         phasefct_even_z(m) = phasefct_init_even(m)  
         phasefct_odd_z(m)  = phasefct_init_odd(m)
      end do

      do i = lo(1)+1,hi(1)
         mi = (i-lo(1))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_x(mi-num_modes);
              phasefct_even_x(mi) = phasefct_mult_even(m,1) * phasefct_even_x(mi-num_modes) - &
                                    phasefct_mult_odd (m,1) * phasefct_odd_x(mi-num_modes)
              phasefct_odd_x(mi)  = phasefct_mult_even(m,1) * phasefct_odd_x(mi-num_modes) + &
                                    phasefct_mult_odd (m,1) * buf(m)
              mi = mi + 1
         end do
      end do         

      do j = lo(2)+1,hi(2)
         mj = (j-lo(2))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_y(mj-num_modes);
              phasefct_even_y(mj) = phasefct_mult_even(m,2) * phasefct_even_y(mj-num_modes) - &
                                    phasefct_mult_odd (m,2) * phasefct_odd_y(mj-num_modes)
              phasefct_odd_y(mj)  = phasefct_mult_even(m,2) * phasefct_odd_y(mj-num_modes) + &
                                    phasefct_mult_odd (m,2) * buf(m)
              mj = mj + 1
         end do
      end do         

      do k = lo(3)+1, hi(3)
         mk = (k-lo(3))*num_modes + 1
         do m = 1, num_modes
              buf(m) = phasefct_even_z(mk-num_modes);
              phasefct_even_z(mk) = phasefct_mult_even(m,3) * phasefct_even_z(mk-num_modes) - &
                                    phasefct_mult_odd (m,3) * phasefct_odd_z(mk-num_modes)
              phasefct_odd_z(mk)  = phasefct_mult_even(m,3) * phasefct_odd_z(mk-num_modes) + &
                                    phasefct_mult_odd (m,3) * buf(m)
              mk = mk + 1
         end do
      end do

      ! compute acceleration component in physical space
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               accel = ZERO

               mi = (i-lo(1))*num_modes + 1 ! offset in x-direction
               mj = (j-lo(2))*num_modes + 1 ! offset in y-direction
               mk = (k-lo(3))*num_modes + 1 ! offset in z-direction
  
               do m = 1, num_modes
                  ! sum up even modes
                  accel = accel + &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_even_z(mk) - &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_odd_z(mk))  * modes_even(m,3)
                  ! sum up odd modes
                  accel = accel - &
                          ((phasefct_even_x(mi) * phasefct_even_y(mj) - &
                            phasefct_odd_x(mi)  * phasefct_odd_y(mj))  * phasefct_odd_z(mk)  + &
                           (phasefct_even_x(mi) * phasefct_odd_y(mj)  + &
                            phasefct_odd_x(mi)  * phasefct_even_y(mj)) * phasefct_even_z(mk)) * modes_odd(m,3)

                  mi = mi + 1
                  mj = mj + 1
                  mk = mk + 1
               end do

               force(i,j,k,1) = M_SQRT_2 * dat(i,j,k,1)* accel

            end do
         end do
      end do

      end subroutine derforcez

