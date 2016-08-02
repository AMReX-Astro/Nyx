! ::
! :: ----------------------------------------------------------
! ::

      !===========================================================================
      ! This version is called assuming no threaded loops, so the threading happens here... 
      !===========================================================================
      subroutine fort_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                                  uout_h1,uout_h2,uout_h3,lo,hi, &
                                                  print_fortran_warnings)

        use meth_params_module, only : NVAR
        use threadbox_module, only : get_lo_hi
        use bl_constants_module
  
        implicit none

        integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
        double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
        integer          :: print_fortran_warnings

        ! Local variables
        integer              :: lo(3), hi(3)
        integer              :: ib,jb,kb
        integer              :: nb(3), boxsize(3), tlo(3), thi(3), iblock, iblockxy, nblocks, nblocksxy
        integer, parameter   :: xblksize = 2048, yblksize = 8, zblksize = 8
        integer, allocatable :: bxlo(:), bxhi(:), bylo(:), byhi(:), bzlo(:), bzhi(:)

        boxsize = hi-lo+1
        nb(1) = max(boxsize(1)/xblksize, 1)
        nb(2) = max(boxsize(2)/yblksize, 1)
        nb(3) = max(boxsize(3)/zblksize, 1)

        allocate(bxlo(0:nb(1)-1))
        allocate(bxhi(0:nb(1)-1))
        allocate(bylo(0:nb(2)-1))
        allocate(byhi(0:nb(2)-1))
        allocate(bzlo(0:nb(3)-1))
        allocate(bzhi(0:nb(3)-1))

        call get_lo_hi(boxsize(1), nb(1), bxlo, bxhi)
        call get_lo_hi(boxsize(2), nb(2), bylo, byhi)
        call get_lo_hi(boxsize(3), nb(3), bzlo, bzhi)

        nblocksxy = nb(1)*nb(2)
        nblocks   = nb(1)*nb(2)*nb(3)

        !$omp parallel private(ib,jb,kb,tlo,thi,iblock,iblockxy)
        !$omp do
        do iblock = 0, nblocks-1

           kb = iblock / nblocksxy
           iblockxy = iblock - kb*nblocksxy
           jb = iblockxy / nb(1)
           ib = iblockxy - jb*nb(1)
           
           tlo(1) = lo(1) + bxlo(ib)
           thi(1) = lo(1) + bxhi(ib)
     
           tlo(2) = lo(2) + bylo(jb)
           thi(2) = lo(2) + byhi(jb)
     
           tlo(3) = lo(3) + bzlo(kb)
           thi(3) = lo(3) + bzhi(kb)
     
           call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,tlo,thi,&
                                            print_fortran_warnings)

        end do
        !$omp end do
        !$omp end parallel
      
      end subroutine fort_enforce_nonnegative_species
! ::
! :: ----------------------------------------------------------
! ::

      !===========================================================================
      ! This version is called from within threaded loops so *no* OMP here ...
      !===========================================================================
      subroutine enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                             uout_h1,uout_h2,uout_h3, &
                                             lo,hi,print_fortran_warnings)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer          :: print_fortran_warnings
      double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n
      integer          :: int_dom_spec
      logical          :: any_negative
      double precision :: dom_spec,x

      double precision, parameter :: eps = -1.0d-16

      if (UFS .gt. 0) then

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         any_negative = .false.
         !
         ! First deal with tiny undershoots by just setting them to zero.
         !
         do n = UFS, UFS+nspec-1
           if (uout(i,j,k,n) .lt. 0.d0) then
              x = uout(i,j,k,n)/uout(i,j,k,URHO)
              if (x .gt. eps) then
                 uout(i,j,k,n) = 0.d0
              else
                 any_negative = .true.
              end if
           end if
         end do
         !
         ! We know there are one or more undershoots needing correction.
         !
         if (any_negative) then
            !
            ! Find the dominant species.
            !
            int_dom_spec = UFS
            dom_spec     = uout(i,j,k,int_dom_spec)

            do n = UFS,UFS+nspec-1
              if (uout(i,j,k,n) .gt. dom_spec) then
                dom_spec     = uout(i,j,k,n)
                int_dom_spec = n
              end if
           end do
           !
           ! Now take care of undershoots greater in magnitude than 1e-16.
           !
           do n = UFS, UFS+nspec-1

              if (uout(i,j,k,n) .lt. 0.d0) then

                 x = uout(i,j,k,n)/uout(i,j,k,URHO)
                 !
                 ! Here we only print the bigger negative values.
                 !
                 if (print_fortran_warnings .gt. 0 .and. x .lt. -1.d-2) then
                    !
                    ! A critical region since we usually can't write from threads.
                    !
                    print *,'Correcting nth negative species ',n-UFS+1
                    print *,'   at cell (i,j,k)              ',i,j,k
                    print *,'Negative (rho*X) is             ',uout(i,j,k,n)
                    print *,'Negative      X  is             ',x
                    print *,'Filling from dominant species   ',int_dom_spec-UFS+1
                    print *,'  which had X =                 ',&
                             uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                 end if
                 !
                 ! Take enough from the dominant species to fill the negative one.
                 !
                 uout(i,j,k,int_dom_spec) = uout(i,j,k,int_dom_spec) + uout(i,j,k,n)
                 !
                 ! Test that we didn't make the dominant species negative.
                 !
                 if (uout(i,j,k,int_dom_spec) .lt. 0.d0) then
                    print *,' Just made nth dominant species negative ',int_dom_spec-UFS+1,' at ',i,j,k
                    print *,'We were fixing species ',n-UFS+1,' which had value ',x
                    print *,'Dominant species became ',uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                    call bl_error("Error:: Nyx_3d.f90 :: fort_enforce_nonnegative_species")
                 end if
                 !
                 ! Now set the negative species to zero.
                 !
                 uout(i,j,k,n) = 0.d0

              end if

           enddo
         end if
      enddo
      enddo
      enddo

      end if ! UFS > 0

      end subroutine enforce_nonnegative_species
