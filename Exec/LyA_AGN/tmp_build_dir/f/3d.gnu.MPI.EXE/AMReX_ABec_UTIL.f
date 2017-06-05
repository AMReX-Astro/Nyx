

c-----------------------------------------------------------------------
c
c     Tridiagonal solve
c
      subroutine tridiag(a,b,c,r,u,n)

      integer n
      integer nmax

      DOUBLE PRECISION a(n)
      DOUBLE PRECISION b(n)
      DOUBLE PRECISION c(n)
      DOUBLE PRECISION r(n)
      DOUBLE PRECISION u(n)

      parameter (nmax = 4098)

      integer j
      DOUBLE PRECISION bet
      DOUBLE PRECISION gam(nmax)
      if (n .gt. nmax ) call bl_error('tridiag: size exceeded')
      if (b(1) .eq. 0) call bl_error('tridiag: CANT HAVE B(1) = ZERO')

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) call bl_error('tridiag: TRIDIAG FAILED')
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

      return
      end

