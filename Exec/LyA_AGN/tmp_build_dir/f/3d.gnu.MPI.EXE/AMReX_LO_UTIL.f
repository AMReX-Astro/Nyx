

c     polyInterpCoeff:
c  
c     This routine returns the Lagrange interpolating coefficients for a
c     polynomial through N points, evaluated at xInt (see Numerical Recipes,
c     v2, p102, e.g.):
c
c            (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
c    P(x) = ----------------------- y1  + ... + ------------------------  yN
c           (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
c
c     P(xInt) = sum_(i=1)^(N) y[i]*c[i]
c
      subroutine polyInterpCoeff(xInt, x, N, c)
      implicit none
      integer N, i, j
      DOUBLE PRECISION xInt, x(N), c(N), num, den
      do j=1,N
         num = 1.0D0
         den = 1.0D0
         do i = 1,j-1
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         c(j) = num/den
      end do
      return
      end

