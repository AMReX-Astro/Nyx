

      subroutine average3dgen (c, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3,f, 
     &f_l1, f_l2, f_l3, f_h1, f_h2, f_h3,lo, hi, nc)
      implicit none
      integer nc
      integer c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
      integer f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
      integer lo(3)
      integer hi(3)
      DOUBLE PRECISION f(f_l1:f_h1, f_l2:f_h2, f_l3:f_h3,nc)
      DOUBLE PRECISION c(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3,nc)

      integer i, i2, i2p1, j, j2, j2p1, k, k2, k2p1, n

      do n = 1, nc
         do k = lo(3), hi(3)
            k2 = 2*k
            k2p1 = k2 + 1
            do j = lo(2), hi(2)
               j2 = 2*j
               j2p1 = j2 + 1
               do i = lo(1), hi(1)
                  i2 = 2*i
                  i2p1 = i2 + 1
                  c(i,j,k,n) = (+ f(i2p1,j2p1,k2 ,n) + f(i2,j2p1,k2 ,n)+
     & f(i2p1,j2 ,k2 ,n) + f(i2,j2 ,k2 ,n)+ f(i2p1,j2p1,k2
     &p1,n) + f(i2,j2p1,k2p1,n)+ f(i2p1,j2 ,k2p1,n) + f(i2
     &,j2 ,k2p1,n))*0.125D0
               end do
            end do
         end do
      end do

      end

      subroutine interp3dgen (f, f_l1, f_l2, f_l3, f_h1, f_h2, f_h3,c, c
     &_l1, c_l2, c_l3, c_h1, c_h2, c_h3,lo, hi, nc)
      implicit none
      integer nc
      integer f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
      integer c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
      integer lo(3)
      integer hi(3)
      DOUBLE PRECISION f(f_l1:f_h1, f_l2:f_h2, f_l3:f_h3,nc)
      DOUBLE PRECISION c(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3,nc)

      integer i, i2, i2p1, j, j2, j2p1, k, k2, k2p1, n

!     MultiGrid::relax(...) does only V-cycles (not F-cycles), and for V-cycles, 
!     piecewise-constant interpolation performs better than linear interpolation,
!     as measured both by run-time and number of V-cycles for convergence.
      do n = 1, nc
         do k = lo(3), hi(3)
            k2 = 2*k
            k2p1 = k2 + 1
            do j = lo(2), hi(2)
               j2 = 2*j
               j2p1 = j2 + 1

               do i = lo(1), hi(1)
                  i2 = 2*i
                  i2p1 = i2 + 1

                  f(i2p1,j2p1,k2  ,n) = c(i,j,k,n) + f(i2p1,j2p1,k2  ,n)
                  f(i2  ,j2p1,k2  ,n) = c(i,j,k,n) + f(i2  ,j2p1,k2  ,n)
                  f(i2p1,j2  ,k2  ,n) = c(i,j,k,n) + f(i2p1,j2  ,k2  ,n)
                  f(i2  ,j2  ,k2  ,n) = c(i,j,k,n) + f(i2  ,j2  ,k2  ,n)
                  f(i2p1,j2p1,k2p1,n) = c(i,j,k,n) + f(i2p1,j2p1,k2p1,n)
                  f(i2  ,j2p1,k2p1,n) = c(i,j,k,n) + f(i2  ,j2p1,k2p1,n)
                  f(i2p1,j2  ,k2p1,n) = c(i,j,k,n) + f(i2p1,j2  ,k2p1,n)
                  f(i2  ,j2  ,k2p1,n) = c(i,j,k,n) + f(i2  ,j2  ,k2p1,n)

               end do
            end do
         end do
      end do

      end
