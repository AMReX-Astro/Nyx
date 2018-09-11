
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on overdensity
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nc        => number of components in density array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine tag_overdensity(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                 set,clear, &
                                 den,denl1,denl2,denl3,denh1,denh2,denh3, &
                                 lo,hi,nc,domlo,domhi,delta,level,avg_den)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      implicit none

      integer set, clear, nc, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      real(rt) delta(3), avg_den

      real(rt) :: over_den

      integer i, j, k

      over_den = 1.1d0 * avg_den

!     Tag on regions of overdensity
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( den(i,j,k,1) .gt. over_den ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_overdensity

! ::: -----------------------------------------------------------
! ::: This routine will tag a specific region.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_region(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            set,lo,hi,domlo,domhi, &
                            delta,xlo,problo,level)
      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      implicit none

      integer          :: set,level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: delta(3), xlo(3), problo(3)

      integer          :: ilo,jlo,klo,ihi,jhi,khi
      integer          :: i,j,k

      if (level .eq. 0) then
         ilo = (domhi(1)+1)*1/4
         ihi = (domhi(1)+1)*3/4 - 1
         jlo = (domhi(2)+1)*1/4
         jhi = (domhi(2)+1)*3/4 - 1
         klo = (domhi(3)+1)*1/4
         khi = (domhi(3)+1)*3/4 - 1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            tag(i,j,k) = set
         end do
         end do
         end do
      else
         ilo = (domhi(1)+1)*3/8
         ihi = (domhi(1)+1)*5/8 - 1
         jlo = (domhi(2)+1)*3/8
         jhi = (domhi(2)+1)*5/8 - 1
         klo = (domhi(3)+1)*3/8
         khi = (domhi(3)+1)*5/8 - 1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (i.ge.ilo .and. i.le.ihi .and. &
                j.ge.jlo .and. j.le.jhi .and. &
                k.ge.klo .and. k.le.khi) then
               tag(i,j,k) = set
            end if
         end do
         end do
         end do
      endif

      end subroutine tag_region
