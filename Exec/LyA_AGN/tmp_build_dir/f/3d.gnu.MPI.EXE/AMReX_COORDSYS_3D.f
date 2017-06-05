

c :: ----------------------------------------------------------
c :: SETVOL
c ::             Compute the volume of each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  vol         <=  volume array
c ::  vlo,vhi      => index limits of vol array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine setvol(reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3,
     &vol,vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3,offset,dx,co
     &ord)
       implicit none
       integer    reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3
       integer    vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3
       integer    coord
       DOUBLE PRECISION     dx(3), offset(3)
       DOUBLE PRECISION vol(vol_l1:vol_h1, vol_l2:vol_h2, vol_l3:vol_h3)

       integer    i, j, k
       DOUBLE PRECISION     v

       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          v = dx(1)*dx(2)*dx(3)
          do k = reg_l3, reg_h3
             do j = reg_l2, reg_h2
                do i = reg_l1, reg_h1
                   vol(i,j,k) = v
                end do
             end do
          end do
       else
          write(6,*) "FORT_SETVOLUME not define for coord = ",coord
          call bl_abort(" ")
       end if

       end

c :: ----------------------------------------------------------
c :: SETAREA
c ::             Compute the area of given cell face
c ::
c :: INPUTS / OUTPUTS:
c ::  area        <=  area array
c ::  alo,ahi      => index limits of area array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 =cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine setarea(reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3
     &,area,area_l1, area_l2, area_l3, area_h1, area_h2, area_h3,offs
     &et,dx,dir,coord)
       implicit none
       integer    reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3
       integer    area_l1, area_l2, area_l3, area_h1, area_h2, area_h3
       integer    coord, dir
       DOUBLE PRECISION     dx(3), offset(3)
       DOUBLE PRECISION area(area_l1:area_h1, area_l2:area_h2, area_l3:a
     &rea_h3)

       integer    i, j, k
       DOUBLE PRECISION     fa

       fa = 1.d0

       if (coord .eq. 0) then
          if (dir .eq. 0) then
             fa = dx(2)*dx(3)
          else if (dir .eq. 1) then
             fa = dx(1)*dx(3)
          else if (dir .eq. 2) then
             fa = dx(1)*dx(2)
          else
             write(6,*) "FORT_SETAREA: invalid dir = ",dir
             call bl_abort(" ")
          end if
          do k = reg_l3, reg_h3
             do j = reg_l2, reg_h2
                do i = reg_l1, reg_h1
                   area(i,j,k) = fa
                end do
             end do
          end do
       else
          write(6,*) "FORT_SETAREA not define for coord = ",coord
          call bl_abort(" ")
       end if

       end

c :: SETDLOGA
c ::             Compute  d(log(A))/dr in each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  dloga        <=  dloga array
c ::  dlo,dhi      => index limits of dloga array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine setdloga(dloga,dloga_l1, dloga_l2, dloga_l3, dloga_h1,
     & dloga_h2, dloga_h3,offset,dx,dir,coord)
       implicit none
       integer dloga_l1, dloga_l2, dloga_l3, dloga_h1, dloga_h2, dloga_h
     &3
       integer    coord
       DOUBLE PRECISION     dx(3), offset(3)
       DOUBLE PRECISION dloga(dloga_l1:dloga_h1, dloga_l2:dloga_h2, dlog
     &a_l3:dloga_h3)
       integer dir

       integer    i, j, k

       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          do k = dloga_l3, dloga_h3
             do j = dloga_l2, dloga_h2
                do i = dloga_l1, dloga_h1
                   dloga(i,j,k) = 0.0D0
                end do
             end do
          enddo
       else 
         write(6,*)' non-cartesian not allowed in 3D yet'
         call bl_abort(" ")

      endif

      end
