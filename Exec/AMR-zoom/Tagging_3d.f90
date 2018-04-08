! ::: -----------------------------------------------------------
! ::: This routine will tag cells based on position
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

      integer  :: set, clear, nc, level
      integer  :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer  :: denl1,denl2,denl3,denh1,denh2,denh3
      integer  :: lo(3), hi(3), domlo(3), domhi(3)
      integer  :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      real(rt) :: delta(3), avg_den

      integer  :: i,j,k
      integer  :: ref_size(3), center(3), ilo(3), ihi(3)

      ref_size = domhi / (2*2**(level+1))
      center   = (domhi-domlo+1) / 2
      ilo      = max(center-ref_size+1, lo)
      ihi      = min(center+ref_size,   hi)

      ! Tag the region
      tag(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3)) = set

      end subroutine tag_overdensity
