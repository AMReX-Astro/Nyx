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

      integer  :: set, clear, nc, level
      integer  :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer  :: denl1,denl2,denl3,denh1,denh2,denh3
      integer  :: lo(3), hi(3), domlo(3), domhi(3)
      integer  :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      real(rt) :: delta(3), avg_den

      integer  :: i,j,k
      real(rt) :: m_cell, V_cell 

      V_cell = delta(1)*delta(2)*delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               m_cell = den(i,j,k,1)*V_cell
               if ( m_cell .gt. 3.5d9 ) then 
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_overdensity

