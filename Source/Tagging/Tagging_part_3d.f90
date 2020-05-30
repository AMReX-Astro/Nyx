
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the number of particles.
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
      subroutine tag_part_cnt_err(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                  set,clear, &
                                  var,varl1,varl2,varl3,varh1,varh2,varh3, &
                                  lo,hi,nd,domlo,domhi, &
                                  delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: var(varl1:varh1,varl2:varh2,varl3:varh3,nd)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: ilo,jlo,klo,ihi,jhi,khi
      integer          :: i,j,k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         if (var(i,j,k,1) .gt. max_num_part) then
            tag(i,j,k) = set
         end if
      end do
      end do
      end do

      end subroutine tag_part_cnt_err

