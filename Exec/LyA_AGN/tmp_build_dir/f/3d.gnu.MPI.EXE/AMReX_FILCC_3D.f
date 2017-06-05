

c ::: -----------------------------------------------------------
c ::: This routine is intended to be a generic fill function
c ::: for cell centered data.  It knows how to exrapolate,
c ::: and reflect data and can be used to suppliment problem
c ::: specific fill functions (ie. 3).
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: q        <=  array to fill
c ::: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3   => index extent of q array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of q array
c ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
c ::: 
c ::: NOTE: corner data not used in computing soln but must have
c :::       reasonable values for arithmetic to live
c ::: -----------------------------------------------------------
      subroutine filcc(q,q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,domlo,domhi,
     &dx,xlo,bc)

      implicit none

      integer    q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer    domlo(3), domhi(3)
      DOUBLE PRECISION     xlo(3), dx(3)
      DOUBLE PRECISION     q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3)
      integer    bc(3,2)

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(q_l1,domlo(1))
      ie = min(q_h1,domhi(1))
      js = max(q_l2,domlo(2))
      je = min(q_h2,domhi(2))
      ks = max(q_l3,domlo(3))
      ke = min(q_h3,domhi(3))

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-domhi(1))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-domhi(2))
      ndwn = max(0,domlo(3)-q_l3)
      nup  = max(0,q_h3-domhi(3))

c
c     ::::: first fill sides
c
      if (nlft .gt. 0) then
         ilo = domlo(1)

         if (bc(1,1) .eq. 3) then
c     set all ghost cell values to a prescribed dirichlet
c     value; in this example, we have chosen 1
c	    do i = 1, nlft
c           do k = q_l3,q_h3
c           do j = q_l2,q_h2
c              q(ilo-i,j,k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(1,1) .eq. 2) then
            do i = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-i,j,k) = q(ilo,j,k)
            end do
            end do
            end do
         else if (bc(1,1) .eq. 4) then
            do i = 2, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-i,j,k) = q(ilo,j,k)
            end do
            end do
            end do
            if (ilo+2 .le. ie) then
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-1,j,k) = (15*q(ilo,j,k) - 10*q(ilo+1,j,k) + 3*q(
     &ilo+2,j,k)) * 0.125D0
               end do
               end do
            else  
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-1,j,k) = 0.5D0*(3*q(ilo,j,k) - q(ilo+1,j,k))
               end do
               end do
            end if
         else if (bc(1,1) .eq. 1) then
            do i = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-i,j,k) = q(ilo+i-1,j,k)
            end do
            end do
            end do
         else if (bc(1,1) .eq. -1) then
            do i = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-i,j,k) = -q(ilo+i-1,j,k)
            end do
            end do
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         if (bc(1,2) .eq. 3) then
c	    do i = 1, nrgt
c           do k = q_l3,q_h3
c           do j = q_l2,q_h2
c              q(ihi+i,j,k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(1,2) .eq. 2) then
            do i = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+i,j,k) = q(ihi,j,k)
            end do
            end do
            end do
         else if (bc(1,2) .eq. 4) then
            do i = 2, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+i,j,k) = q(ihi,j,k)
            end do
            end do
            end do
            if (ihi-2 .ge. is) then
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+1,j,k) = (15*q(ihi,j,k) - 10*q(ihi-1,j,k) + 3*q(
     &ihi-2,j,k)) * 0.125D0
               end do
               end do
            else
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+1,j,k) = 0.5D0*(3*q(ihi,j,k) - q(ihi-1,j,k))
               end do
               end do
            end if
         else if (bc(1,2) .eq. 1) then
            do i = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+i,j,k) = q(ihi-i+1,j,k)
            end do
            end do
            end do
         else if (bc(1,2) .eq. -1) then
            do i = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+i,j,k) = -q(ihi-i+1,j,k)
            end do
            end do
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         if (bc(2,1) .eq. 3) then
c	    do j = 1, nbot
c           do k = q_l3,q_h3
c           do i = q_l1,q_h1
c              q(i,jlo-j,k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(2,1) .eq. 2) then
            do j = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-j,k) = q(i,jlo,k)
            end do
            end do
            end do
         else if (bc(2,1) .eq. 4) then
            do j = 2, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-j,k) = q(i,jlo,k)
            end do
            end do
            end do
            if (jlo+2 .le. je) then
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-1,k) = (15*q(i,jlo,k) - 10*q(i,jlo+1,k) + 3*q(
     &i,jlo+2,k)) * 0.125D0
               end do
               end do
            else
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-1,k) = 0.5D0*(3*q(i,jlo,k) - q(i,jlo+1,k))
               end do
               end do
            end if
         else if (bc(2,1) .eq. 1) then
            do j = 1, nbot 
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-j,k) = q(i,jlo+j-1,k)
            end do
            end do
            end do
         else if (bc(2,1) .eq. -1) then
            do j = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-j,k) = -q(i,jlo+j-1,k)
            end do
            end do
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         if (bc(2,2) .eq. 3) then
c	    do j = 1, ntop
c           do k = q_l3,q_h3
c           do i = q_l1,q_h1
c              q(i,jhi+j,k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(2,2) .eq. 2) then
            do j = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+j,k) = q(i,jhi,k)
            end do
            end do
            end do
         else if (bc(2,2) .eq. 4) then
            do j = 2, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+j,k) = q(i,jhi,k)
            end do
            end do
            end do
            if (jhi-2 .ge. js) then
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+1,k) = (15*q(i,jhi,k) - 10*q(i,jhi-1,k) + 3*q(
     &i,jhi-2,k)) * 0.125D0
               end do
               end do
            else
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+1,k) = 0.5D0*(3*q(i,jhi,k) - q(i,jhi-1,k))
               end do
               end do
            end if
         else if (bc(2,2) .eq. 1) then
            do j = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+j,k) = q(i,jhi-j+1,k)
            end do
            end do
            end do
         else if (bc(2,2) .eq. -1) then
            do j = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+j,k) = -q(i,jhi-j+1,k)
            end do
            end do
            end do
         end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

         if (bc(3,1) .eq. 3) then
c	    do k = 1, ndwn
c           do j = q_l2,q_h2
c           do i = q_l1,q_h1
c              q(i,j,klo-k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(3,1) .eq. 2) then
            do k = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-k) = q(i,j,klo)
            end do
            end do
            end do
         else if (bc(3,1) .eq. 4) then
            do k = 2, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-k) = q(i,j,klo)
            end do
            end do
            end do
            if (klo+2 .le. ke) then
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-1) = (15*q(i,j,klo) - 10*q(i,j,klo+1) + 3*q(
     &i,j,klo+2)) * 0.125D0
               end do
               end do
            else
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-1) = 0.5D0*(3*q(i,j,klo) - q(i,j,klo+1))
               end do
               end do
            end if
         else if (bc(3,1) .eq. 1) then
            do k = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-k) = q(i,j,klo+k-1)
            end do
            end do
            end do
         else if (bc(3,1) .eq. -1) then
            do k = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-k) = -q(i,j,klo+k-1)
            end do
            end do
            end do
         end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

         if (bc(3,2) .eq. 3) then
c	    do k = 1, nup
c           do j = q_l2,q_h2
c           do i = q_l1,q_h1
c              q(i,j,khi+k) = 1.d0
c           end do
c           end do
c	    end do
         else if (bc(3,2) .eq. 2) then
            do k = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+k) = q(i,j,khi)
            end do
            end do
            end do
         else if (bc(3,2) .eq. 4) then
            do k = 2, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+k) = q(i,j,khi)
            end do
            end do
            end do
            if (khi-2 .ge. ks) then
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+1) = (15*q(i,j,khi) - 10*q(i,j,khi-1) + 3*q(
     &i,j,khi-2)) * 0.125D0
               end do
               end do
            else
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+1) = 0.5D0*(3*q(i,j,khi) - q(i,j,khi-1))
               end do
               end do
            end if
         else if (bc(3,2) .eq. 1) then
            do k = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+k) = q(i,j,khi-k+1)
            end do
            end do
            end do
         else if (bc(3,2) .eq. -1) then
            do k = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+k) = -q(i,j,khi-k+1)
            end do
            end do
            end do
         end if
      end if
c
c    First correct the i-j edges and all corners
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(nbot .gt. 0 .and. bc(
     &2,1) .eq. 4) ) then
         if (jlo+2 .le. je) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.5D0 * 0.125D0 * (15*q(ilo-1,jlo,k) -
     & 10*q(ilo-1,jlo+1,k) + 3*q(ilo-1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.5D0 * 0.5D0 * (3*q(ilo-1,jlo,k) - q(
     &ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.5D0 * 0.125D0 * (
     &15*q(ilo,jlo-1,k) - 10*q(ilo+1,jlo-1,k) + 3*q(ilo+2,jlo
     &-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.5D0 * 0.5D0 * (3*
     &q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. 4) then
            if (klo+2 .le. ke) then
               q(ilo-1,jlo-1,klo-1) = 0.125D0 * ((15*q(ilo-1,jlo-1,klo) 
     &- 10*q(ilo-1,jlo-1,klo+1) +3*q(ilo-1,jlo-1,klo+2)) )
            else
               q(ilo-1,jlo-1,klo-1) = 0.5D0 * (3*q(ilo-1,jlo-1,klo) - q(
     &ilo-1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. 4) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jlo-1,khi+1) = 0.125D0 * ((15*q(ilo-1,jlo-1,khi) 
     &- 10*q(ilo-1,jlo-1,khi-1) +3*q(ilo-1,jlo-1,khi-2)) )
            else
               q(ilo-1,jlo-1,khi+1) = 0.5D0 * (3*q(ilo-1,jlo-1,khi) - q(
     &ilo-1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(ntop .gt. 0 .and. bc(
     &2,2) .eq. 4) ) then
         if (jhi-2 .ge. js) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5D0 * 0.125D0 * (15*q(ilo-1,jhi,k) -
     & 10*q(ilo-1,jhi-1,k) + 3*q(ilo-1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5D0 * 0.5D0 * (3*q(ilo-1,jhi,k) - q(
     &ilo-1,jhi-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5D0 * 0.125D0 * (
     &15*q(ilo,jhi+1,k) - 10*q(ilo+1,jhi+1,k) + 3*q(ilo+2,jhi
     &+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5D0 * 0.5D0 * (3*
     &q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. 4) then
            if (klo+2 .le. ke) then
               q(ilo-1,jhi+1,klo-1) = 0.125D0 * ((15*q(ilo-1,jhi+1,klo) 
     &- 10*q(ilo-1,jhi+1,klo+1) +3*q(ilo-1,jhi+1,klo+2)) )
            else
               q(ilo-1,jhi+1,klo-1) = 0.5D0 * (3*q(ilo-1,jhi+1,klo) - q(
     &ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. 4) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jhi+1,khi+1) = 0.125D0 * ((15*q(ilo-1,jhi+1,khi) 
     &- 10*q(ilo-1,jhi+1,khi-1) +3*q(ilo-1,jhi+1,khi-2)) )
            else
               q(ilo-1,jhi+1,khi+1) = 0.5D0 * (3*q(ilo-1,jhi+1,khi) - q(
     &ilo-1,jhi+1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(nbot .gt. 0 .and. bc(
     &2,1) .eq. 4) ) then
         if (jlo+2 .le. je) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5D0 * 0.125D0 * (15*q(ihi+1,jlo,k) -
     & 10*q(ihi+1,jlo+1,k) + 3*q(ihi+1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5D0 * 0.5D0 * (3*q(ihi+1,jlo,k) - q(
     &ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5D0 * 0.125D0 * (
     &15*q(ihi,jlo-1,k) - 10*q(ihi-1,jlo-1,k) + 3*q(ihi-2,jlo
     &-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5D0 * 0.5D0 * (3*
     &q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. 4) then
            if (klo+2 .le. ke) then
               q(ihi+1,jlo-1,klo-1) = 0.125D0 * (15*q(ihi+1,jlo-1,klo) -
     & 10*q(ihi+1,jlo-1,klo+1) +3*q(ihi+1,jlo-1,klo+2))
            else
               q(ihi+1,jlo-1,klo-1) = 0.5D0 * (3*q(ihi+1,jlo-1,klo) - q(
     &ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. 4) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jlo-1,khi+1) = 0.125D0 * (15*q(ihi+1,jlo-1,khi) -
     & 10*q(ihi+1,jlo-1,khi-1) +3*q(ihi+1,jlo-1,khi-2))
            else
               q(ihi+1,jlo-1,khi+1) = 0.5D0 * (3*q(ihi+1,jlo-1,khi) - q(
     &ihi+1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(ntop .gt. 0 .and. bc(
     &2,2) .eq. 4) ) then
         if (jhi-2 .ge. js) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5D0 * 0.125D0 * (15*q(ihi+1,jhi,k) -
     & 10*q(ihi+1,jhi-1,k) + 3*q(ihi+1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5D0 * 0.5D0 * (3*q(ihi+1,jhi,k) - q(
     &ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5D0 * 0.125D0 * (
     &15*q(ihi,jhi+1,k) - 10*q(ihi-1,jhi+1,k) + 3*q(ihi-2,jhi
     &+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5D0 * 0.5D0 * (3*
     &q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. 4) then
            if (klo+2 .le. ke) then
               q(ihi+1,jhi+1,klo-1) = 0.125D0 *(15*q(ihi+1,jhi+1,klo) - 
     &10*q(ihi+1,jhi+1,klo+1) +3*q(ihi+1,jhi+1,klo+2))
            else
               q(ihi+1,jhi+1,klo-1) = 0.5D0 * (3*q(ihi+1,jhi+1,klo) - q(
     &ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. 4) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jhi+1,khi+1) = 0.125D0 *(15*q(ihi+1,jhi+1,khi) - 
     &10*q(ihi+1,jhi+1,khi-1) +3*q(ihi+1,jhi+1,khi-2))
            else
               q(ihi+1,jhi+1,khi+1) = 0.5D0 * (3*q(ihi+1,jhi+1,khi) - q(
     &ihi+1,jhi+1,khi-1))
            end if
         end if

      end if
c
c    Next correct the i-k edges
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(ndwn .gt. 0 .and. bc(
     &3,1) .eq. 4) ) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5D0 * 0.125D0 * (15*q(ilo-1,j,klo) -
     & 10*q(ilo-1,j,klo+1) + 3*q(ilo-1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5D0 * 0.5D0 * (3*q(ilo-1,j,klo) - q(
     &ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5D0 * 0.125D0 * (
     &15*q(ilo,j,klo-1) - 10*q(ilo+1,j,klo-1) + 3*q(ilo+2,j,k
     &lo-1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5D0 * 0.5D0 * (3*
     &q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(nup .gt. 0 .and. bc(3
     &,2) .eq. 4) ) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5D0 * 0.125D0 * (15*q(ilo-1,j,khi) -
     & 10*q(ilo-1,j,khi-1) + 3*q(ilo-1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5D0 * 0.5D0 * (3*q(ilo-1,j,khi) - q(
     &ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5D0 * 0.125D0 * (
     &15*q(ilo,j,khi+1) - 10*q(ilo+1,j,khi+1) + 3*q(ilo+2,j,k
     &hi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5D0 * 0.5D0 * (3*
     &q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(ndwn .gt. 0 .and. bc(
     &3,1) .eq. 4) ) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5D0 * 0.125D0 *(15*q(ihi+1,j,klo) - 
     &10*q(ihi+1,j,klo+1) + 3*q(ihi+1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5D0 * 0.5D0 * (3*q(ihi+1,j,klo) - q(
     &ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5D0 * 0.125D0 *(1
     &5*q(ihi,j,klo-1) - 10*q(ihi-1,j,klo-1) + 3*q(ihi-2,j,kl
     &o-1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5D0 * 0.5D0 * (3*
     &q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(nup .gt. 0 .and. bc(3
     &,2) .eq. 4) ) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5D0 * 0.125D0 * (15*q(ihi+1,j,khi) -
     & 10*q(ihi+1,j,khi-1) + 3*q(ihi+1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5D0 * 0.5D0 * (3*q(ihi+1,j,khi) - q(
     &ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5D0 * 0.125D0 * (
     &15*q(ihi,j,khi+1) - 10*q(ihi-1,j,khi+1) + 3*q(ihi-2,j,k
     &hi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5D0 * 0.5D0 * (3*
     &q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
            end do
         end if
      end if
c
c    Next correct the j-k edges
c
      if ((nbot .gt. 0 .and. bc(2,1) .eq. 4) .and.(ndwn .gt. 0 .and. bc(
     &3,1) .eq. 4) ) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5D0 * 0.125D0 *(15*q(i,jlo-1,klo) - 
     &10*q(i,jlo-1,klo+1) + 3*q(i,jlo-1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5D0 * 0.5D0 * (3*q(i,jlo-1,klo) - q(
     &i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5D0 * 0.125D0 * (
     &15*q(i,jlo,klo-1) - 10*q(i,jlo+1,klo-1) + 3*q(i,jlo+2,k
     &lo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5D0 * 0.5D0 * (3*
     &q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nbot .gt. 0 .and. bc(2,1) .eq. 4) .and.(nup .gt. 0 .and. bc(3
     &,2) .eq. 4) ) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5D0 * 0.125D0 * (15*q(i,jlo-1,khi) -
     & 10*q(i,jlo-1,khi-1) + 3*q(i,jlo-1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5D0 * 0.5D0 *(3*q(i,jlo-1,khi) - q(i
     &,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5D0 * 0.125D0 * (
     &15*q(i,jlo,khi+1) - 10*q(i,jlo+1,khi+1) + 3*q(i,jlo+2,k
     &hi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5D0 * 0.5D0 *(3*q
     &(i,jlo,khi+1) - q(i,jlo+1,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0 .and. bc(2,2) .eq. 4) .and.(ndwn .gt. 0 .and. bc(
     &3,1) .eq. 4) ) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5D0 * 0.125D0 * (15*q(i,jhi+1,klo) -
     & 10*q(i,jhi+1,klo+1) + 3*q(i,jhi+1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5D0 * 0.5D0 * (3*q(i,jhi+1,klo) - q(
     &i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5D0 * 0.125D0 * (
     &15*q(i,jhi,klo-1) - 10*q(i,jhi-1,klo-1) + 3*q(i,jhi-2,k
     &lo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5D0 * 0.5D0 * (3*
     &q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0 .and. bc(2,2) .eq. 4) .and.(nup .gt. 0 .and. bc(3
     &,2) .eq. 4) ) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5D0 * 0.125D0 * (15*q(i,jhi+1,khi) -
     & 10*q(i,jhi+1,khi-1) + 3*q(i,jhi+1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5D0 * 0.5D0 * (3*q(i,jhi+1,khi) - q(
     &i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5D0 * 0.125D0 * (
     &15*q(i,jhi,khi+1) - 10*q(i,jhi-1,khi+1) + 3*q(i,jhi-2,k
     &hi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5D0 * 0.5D0 * (3*
     &q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
            end do
         end if
      end if

      end

      subroutine hoextraptocc(q,q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,domlo
     &,domhi,dx,xlo)

      implicit none

      integer    q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer    domlo(3), domhi(3)
      DOUBLE PRECISION     xlo(3), dx(3)
      DOUBLE PRECISION     q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3)

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(q_l1,domlo(1))
      ie = min(q_h1,domhi(1))
      js = max(q_l2,domlo(2))
      je = min(q_h2,domhi(2))
      ks = max(q_l3,domlo(3))
      ke = min(q_h3,domhi(3))

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-domhi(1))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-domhi(2))
      ndwn = max(0,domlo(3)-q_l3)
      nup  = max(0,q_h3-domhi(3))
c
c     First fill sides.
c
      if (nlft .gt. 0) then
         ilo = domlo(1)

         do i = 2, nlft
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-i,j,k) = q(ilo,j,k)
               end do
            end do
         end do
         if (ilo+2 .le. ie) then
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-1,j,k) = 3*q(ilo,j,k) - 3*q(ilo+1,j,k) +q(ilo+2,
     &j,k)
               end do
            end do
         else  
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-1,j,k) = 2*q(ilo,j,k) - q(ilo+1,j,k)
               end do
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         do i = 2, nrgt
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+i,j,k) = q(ihi,j,k)
               end do
            end do
         end do
         if (ihi-2 .ge. is) then
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+1,j,k) = 3*q(ihi,j,k) - 3*q(ihi-1,j,k) +q(ihi-2,
     &j,k)
               end do
            end do
         else
            do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+1,j,k) = 2*q(ihi,j,k) - q(ihi-1,j,k)
               end do
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         do j = 2, nbot
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-j,k) = q(i,jlo,k)
               end do
            end do
         end do
         if (jlo+2 .le. je) then
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-1,k) = 3*q(i,jlo,k) - 3*q(i,jlo+1,k) +q(i,jlo+
     &2,k)
               end do
            end do
         else
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-1,k) = 2*q(i,jlo,k) - q(i,jlo+1,k)
               end do
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         do j = 2, ntop
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+j,k) = q(i,jhi,k)
               end do
            end do
         end do
         if (jhi-2 .ge. js) then
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+1,k) = 3*q(i,jhi,k) - 3*q(i,jhi-1,k) +q(i,jhi-
     &2,k)
               end do
            end do
         else
            do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+1,k) = 2*q(i,jhi,k) - q(i,jhi-1,k)
               end do
            end do
         end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

         do k = 2, ndwn
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-k) = q(i,j,klo)
               end do
            end do
         end do
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-1) = 3*q(i,j,klo) - 3*q(i,j,klo+1) +q(i,j,kl
     &o+2)
               end do
            end do
         else
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-1) = 2*q(i,j,klo) - q(i,j,klo+1)
               end do
            end do
         end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

         do k = 2, nup
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+k) = q(i,j,khi)
               end do
            end do
         end do
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+1) = 3*q(i,j,khi) - 3*q(i,j,khi-1) +q(i,j,kh
     &i-2)
               end do
            end do
         else
            do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+1) = 2*q(i,j,khi) - q(i,j,khi-1)
               end do
            end do
         end if
      end if
c
c    First correct the i-j edges and all corners
c
      if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.5D0 *(3*q(ilo-1,jlo,k) - 3*q(ilo-1,j
     &lo+1,k) +q(ilo-1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.5D0 *(2*q(ilo-1,jlo,k) - q(ilo-1,jlo
     &+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.5D0 *(3*q(ilo,jlo
     &-1,k) - 3*q(ilo+1,jlo-1,k) + q(ilo+2,jlo-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.5D0 *(2*q(ilo,jlo
     &-1,k) - q(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ilo-1,jlo-1,klo-1) = (3*q(ilo-1,jlo-1,klo) - 3*q(ilo-1,
     &jlo-1,klo+1) +q(ilo-1,jlo-1,klo+2))
            else
               q(ilo-1,jlo-1,klo-1) = (2*q(ilo-1,jlo-1,klo) - q(ilo-1,jl
     &o-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jlo-1,khi+1) = (3*q(ilo-1,jlo-1,khi) - 3*q(ilo-1,
     &jlo-1,khi-1) +q(ilo-1,jlo-1,khi-2))
            else
               q(ilo-1,jlo-1,khi+1) =(2*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo
     &-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5D0 *(3*q(ilo-1,jhi,k) - 3*q(ilo-1,j
     &hi-1,k) + q(ilo-1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5D0 *(2*q(ilo-1,jhi,k) - q(ilo-1,jhi
     &-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5D0 *(3*q(ilo,jhi
     &+1,k) - 3*q(ilo+1,jhi+1,k) + q(ilo+2,jhi+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5D0 *(2*q(ilo,jhi
     &+1,k) - q(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ilo-1,jhi+1,klo-1) = (3*q(ilo-1,jhi+1,klo) - 3*q(ilo-1,
     &jhi+1,klo+1) +q(ilo-1,jhi+1,klo+2))
            else
               q(ilo-1,jhi+1,klo-1) =(2*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi
     &+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jhi+1,khi+1) =(3*q(ilo-1,jhi+1,khi) - 3*q(ilo-1,j
     &hi+1,khi-1) +q(ilo-1,jhi+1,khi-2))
            else
               q(ilo-1,jhi+1,khi+1) =(2*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi
     &+1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5D0 *(3*q(ihi+1,jlo,k) - 3*q(ihi+1,j
     &lo+1,k) + q(ihi+1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5D0 *(2*q(ihi+1,jlo,k) - q(ihi+1,jlo
     &+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5D0 *(3*q(ihi,jlo
     &-1,k) - 3*q(ihi-1,jlo-1,k) + q(ihi-2,jlo-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5D0 *(2*q(ihi,jlo
     &-1,k) - q(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ihi+1,jlo-1,klo-1) =(3*q(ihi+1,jlo-1,klo) - 3*q(ihi+1,j
     &lo-1,klo+1) +q(ihi+1,jlo-1,klo+2))
            else
               q(ihi+1,jlo-1,klo-1) =(2*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo
     &-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jlo-1,khi+1) =(3*q(ihi+1,jlo-1,khi) - 3*q(ihi+1,j
     &lo-1,khi-1) +q(ihi+1,jlo-1,khi-2))
            else
               q(ihi+1,jlo-1,khi+1) =(2*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo
     &-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5D0 *(3*q(ihi+1,jhi,k) - 3*q(ihi+1,j
     &hi-1,k) + q(ihi+1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5D0 *(2*q(ihi+1,jhi,k) - q(ihi+1,jhi
     &-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5D0 *(3*q(ihi,jhi
     &+1,k) - 3*q(ihi-1,jhi+1,k) + q(ihi-2,jhi+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5D0 *(2*q(ihi,jhi
     &+1,k) - q(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ihi+1,jhi+1,klo-1) =(3*q(ihi+1,jhi+1,klo) - 3*q(ihi+1,j
     &hi+1,klo+1) +q(ihi+1,jhi+1,klo+2))
            else
               q(ihi+1,jhi+1,klo-1) =(2*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi
     &+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jhi+1,khi+1) =(3*q(ihi+1,jhi+1,khi) - 3*q(ihi+1,j
     &hi+1,khi-1) +q(ihi+1,jhi+1,khi-2))
            else
               q(ihi+1,jhi+1,khi+1) =(2*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi
     &+1,khi-1))
            end if
         end if

      end if
c
c    Next correct the i-k edges
c
      if ((nlft .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5D0 *(3*q(ilo-1,j,klo) - 3*q(ilo-1,j
     &,klo+1) + q(ilo-1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5D0 *(2*q(ilo-1,j,klo) - q(ilo-1,j,k
     &lo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5D0 *(3*q(ilo,j,k
     &lo-1) - 3*q(ilo+1,j,klo-1) + q(ilo+2,j,klo-1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5D0 *(2*q(ilo,j,k
     &lo-1) - q(ilo+1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5D0 *(3*q(ilo-1,j,khi) - 3*q(ilo-1,j
     &,khi-1) + q(ilo-1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5D0 *(2*q(ilo-1,j,khi) - q(ilo-1,j,k
     &hi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5D0 *(3*q(ilo,j,k
     &hi+1) - 3*q(ilo+1,j,khi+1) + q(ilo+2,j,khi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5D0 *(2*q(ilo,j,k
     &hi+1) - q(ilo+1,j,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5D0 *(3*q(ihi+1,j,klo) - 3*q(ihi+1,j
     &,klo+1) + q(ihi+1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5D0 *(2*q(ihi+1,j,klo) - q(ihi+1,j,k
     &lo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5D0 *(3*q(ihi,j,k
     &lo-1) - 3*q(ihi-1,j,klo-1) + q(ihi-2,j,klo-1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5D0 *(2*q(ihi,j,k
     &lo-1) - q(ihi-1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5D0 *(3*q(ihi+1,j,khi) - 3*q(ihi+1,j
     &,khi-1) + q(ihi+1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5D0 *(2*q(ihi+1,j,khi) - q(ihi+1,j,k
     &hi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5D0 *(3*q(ihi,j,k
     &hi+1) - 3*q(ihi-1,j,khi+1) + q(ihi-2,j,khi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5D0 *(2*q(ihi,j,k
     &hi+1) - q(ihi-1,j,khi+1))
            end do
         end if
      end if
c
c    Next correct the j-k edges
c
      if ((nbot .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5D0 *(3*q(i,jlo-1,klo) - 3*q(i,jlo-1
     &,klo+1) + q(i,jlo-1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5D0 *(2*q(i,jlo-1,klo) - q(i,jlo-1,k
     &lo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5D0 *(3*q(i,jlo,k
     &lo-1) - 3*q(i,jlo+1,klo-1) + q(i,jlo+2,klo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5D0 *(2*q(i,jlo,k
     &lo-1) - q(i,jlo+1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nbot .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5D0 *(3*q(i,jlo-1,khi) - 3*q(i,jlo-1
     &,khi-1) + q(i,jlo-1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5D0 *(2*q(i,jlo-1,khi) - q(i,jlo-1,k
     &hi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5D0 *(3*q(i,jlo,k
     &hi+1) - 3*q(i,jlo+1,khi+1) + q(i,jlo+2,khi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5D0 *(2*q(i,jlo,k
     &hi+1) - q(i,jlo+1,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5D0 *(3*q(i,jhi+1,klo) - 3*q(i,jhi+1
     &,klo+1) + q(i,jhi+1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5D0 *(2*q(i,jhi+1,klo) - q(i,jhi+1,k
     &lo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5D0 *(3*q(i,jhi,k
     &lo-1) - 3*q(i,jhi-1,klo-1) + q(i,jhi-2,klo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5D0 *(2*q(i,jhi,k
     &lo-1) - q(i,jhi-1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5D0 *(3*q(i,jhi+1,khi) - 3*q(i,jhi+1
     &,khi-1) + q(i,jhi+1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5D0 *(2*q(i,jhi+1,khi) - q(i,jhi+1,k
     &hi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5D0 *(3*q(i,jhi,k
     &hi+1) - 3*q(i,jhi-1,khi+1) + q(i,jhi-2,khi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5D0 *(2*q(i,jhi,k
     &hi+1) - q(i,jhi-1,khi+1))
            end do
         end if
      end if

      end
