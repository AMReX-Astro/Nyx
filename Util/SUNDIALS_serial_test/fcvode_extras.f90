module fcvode_extras

  implicit none

  contains

    integer(c_int) function RhsFnReal(tn, yvec, fvec, rpar, neq) &
           result(ierr) bind(C,name='RhsFnReal')

      use, intrinsic :: iso_c_binding
      implicit none


      real(c_double), value :: tn
      integer(c_int), value :: neq
!      type(c_ptr), value    :: sunvec_y
!      type(c_ptr), value    :: sunvec_f
!      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double) :: yvec(neq)
      real(c_double) :: fvec(neq)
      real(c_double), intent(inout) :: rpar(neq*4)
      real(c_double) :: energy(neq)

!      print*, "r1", rpar(1)
!      print*, "r2", rpar(2)
!      print*, rpar(3)
!      print*, rpar(4)
      call f_rhs(neq, tn, yvec, fvec, rpar, 0)
!      print*, "after r1", rpar(1)
!      print*, "after r2", rpar(2)
!      print*, "after r3", rpar(3)
!      print*, "after r4", rpar(4)

      ierr = 0
    end function RhsFnReal
end module fcvode_extras
