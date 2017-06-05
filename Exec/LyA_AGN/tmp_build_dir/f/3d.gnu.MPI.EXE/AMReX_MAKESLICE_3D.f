

      subroutine makemaskslice(data1,data1_l1, data1_l2, data1_l3, data1
     &_h1, data1_h2, data1_h3,data2,data2_l1, data2_l2, data2_l3, data
     &2_h1, data2_h2, data2_h3,data_min,data_max,slice_val,slice,dx)

      implicit none

      integer data1_l1, data1_l2, data1_l3, data1_h1, data1_h2, data1_h3
      integer data2_l1, data2_l2, data2_l3, data2_h1, data2_h2, data2_h3
      integer data1(data1_l1:data1_h1, data1_l2:data1_h2, data1_l3:data1
     &_h3)
      integer data2(data2_l1:data2_h1, data2_l2:data2_h2, data2_l3:data2
     &_h3)
      DOUBLE PRECISION dx(3)
      DOUBLE PRECISION data_min,data_max,slice_val
      integer slice

      integer dlo1(3), dhi1(3)
      integer dlo2(3), dhi2(3)
      integer index,idx1,idx2,found,i,j,k,kval,kval1,kval2
      integer leftval,rightval
      DOUBLE PRECISION zval,tval,floatdata

      kval     = 0
      leftval  = 0
      rightval = 0

      dlo1(1) = data1_l1
      dlo1(2) = data1_l2
      dlo1(3) = data1_l3
      dhi1(1) = data1_h1
      dhi1(2) = data1_h2
      dhi1(3) = data1_h3
      dlo2(1) = data2_l1
      dlo2(2) = data2_l2
      dlo2(3) = data2_l3
      dhi2(1) = data2_h1
      dhi2(2) = data2_h2
      dhi2(3) = data2_h3

      if ((dlo2(3).ne.0).or.(dhi2(3).ne.0)) then
         write (6,*) "second data must be thin!"
      end if

      if (slice.eq.0) then
         idx1=1
         idx2=2
         index=3
      else if (slice.eq.1) then
         idx1=1
         idx2=3
         index=2
      else
         idx1=2
         idx2=3
         index=1
      end if

      if ( (dlo1(idx1).ne.dlo2(1)).or.(dlo1(idx2).ne.dlo2(2)).or.(dhi1(i
     &dx1).ne.dhi2(1)).or.(dhi1(idx2).ne.dhi2(2)) ) then
         write (6,*) "thin data must jive with 3d data!"
      end if

      found=0
      do k=dlo1(index)-1,dhi1(index)
         zval=(k+0.5D0)*dx(index)
         if ((slice_val.ge.zval).and.(slice_val.le.zval+dx(index)).and.(
     &found.eq.0)) then
            found=1
            kval=k 
         end if
      end do
      if (found.eq.0) then
         write (6,*) "slice data out of range!"
      end if
      kval1=kval
      kval2=kval+1
      if (kval1.lt.dlo1(index)) then
         kval1=kval2
      end if
      if (kval2.gt.dhi1(index)) then
         kval2=kval1
      end if

      do i=dlo1(idx1),dhi1(idx1)
         do j=dlo1(idx2),dhi1(idx2)
            zval=(kval+0.5D0)*dx(index)
            if (slice.eq.0) then
               leftval=data1(i,j,kval1)
               rightval=data1(i,j,kval2)
            else if (slice.eq.1) then
               leftval=data1(i,kval1,j)
               rightval=data1(i,kval2,j)
            else if (slice.eq.2) then
               leftval=data1(kval1,i,j)
               rightval=data1(kval2,i,j)
            end if
            tval=(slice_val-zval)/dx(index)
            floatdata=leftval*(1.0D0-tval)+rightval*tval

            if (floatdata.gt.0.0D0) then
               data2(i,j,0)=1
            else
               data2(i,j,0)=0
            end if
         end do
      end do

      end

      subroutine makeslice(data1,data1_l1, data1_l2, data1_l3, data1_h1,
     & data1_h2, data1_h3,data2,data2_l1, data2_l2, data2_l3, data2_h1
     &, data2_h2, data2_h3,data_min,data_max,slice_val,slice,dx)

      implicit none

      integer data1_l1, data1_l2, data1_l3, data1_h1, data1_h2, data1_h3
      integer data2_l1, data2_l2, data2_l3, data2_h1, data2_h2, data2_h3
      DOUBLE PRECISION data1(data1_l1:data1_h1, data1_l2:data1_h2, data1
     &_l3:data1_h3)
      DOUBLE PRECISION data2(data2_l1:data2_h1, data2_l2:data2_h2, data2
     &_l3:data2_h3)
      DOUBLE PRECISION dx(3)
      DOUBLE PRECISION data_min,data_max,slice_val
      integer slice

      integer index,idx1,idx2,found,i,j,k,kval,kval1,kval2
      DOUBLE PRECISION zval,leftval,rightval,tval

      integer dlo1(3),dhi1(3)
      integer dlo2(3),dhi2(3)

      kval     = 0
      leftval  = 0
      rightval = 0

      dlo1(1) = data1_l1
      dlo1(2) = data1_l2
      dlo1(3) = data1_l3
      dhi1(1) = data1_h1
      dhi1(2) = data1_h2
      dhi1(3) = data1_h3
      dlo2(1) = data2_l1
      dlo2(2) = data2_l2
      dlo2(3) = data2_l3
      dhi2(1) = data2_h1
      dhi2(2) = data2_h2
      dhi2(3) = data2_h3

      if ((dlo2(3).ne.0).or.(dhi2(3).ne.0)) then
         write (6,*) "second data must be thin!"
      end if

      if (slice.eq.0) then
         idx1=1
         idx2=2
         index=3
      else if (slice.eq.1) then
         idx1=1
         idx2=3
         index=2
      else
         idx1=2
         idx2=3
         index=1
      end if

      if ( (dlo1(idx1).ne.dlo2(1)).or.(dlo1(idx2).ne.dlo2(2)).or.(dhi1(i
     &dx1).ne.dhi2(1)).or.(dhi1(idx2).ne.dhi2(2)) ) then
         write (6,*) "thin data must jive with 3d data!"
      end if

      found=0
      do k=dlo1(index)-1,dhi1(index)
         zval=(k+0.5D0)*dx(index)
         if ((slice_val.ge.zval).and.(slice_val.le.zval+dx(index)).and.(
     &found.eq.0)) then
            found=1
            kval=k 
         end if
      end do
      if (found.eq.0) then
         write (6,*) "slice data out of range!"
      end if
      kval1=kval
      kval2=kval+1
      if (kval1.lt.dlo1(index)) then
         kval1=kval2
      end if
      if (kval2.gt.dhi1(index)) then
         kval2=kval1
      end if

      do i=dlo1(idx1),dhi1(idx1)
         do j=dlo1(idx2),dhi1(idx2)
            zval=(kval+0.5D0)*dx(index)
            if (slice.eq.0) then
               leftval=data1(i,j,kval1)
               rightval=data1(i,j,kval2)
            else if (slice.eq.1) then
               leftval=data1(i,kval1,j)
               rightval=data1(i,kval2,j)
            else if (slice.eq.2) then
               leftval=data1(kval1,i,j)
               rightval=data1(kval2,i,j)
            end if
            tval=(slice_val-zval)/dx(index)
            data2(i,j,0)=leftval*(1.0D0-tval)+rightval*tval
         end do
      end do

      end

