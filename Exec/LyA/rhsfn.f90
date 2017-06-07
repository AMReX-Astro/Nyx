module rhs
  implicit none

  contains

    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(C,name='RhsFn')

      use, intrinsic :: iso_c_binding
      use fnvector_serial
      use cvode_interface
      implicit none

      real(c_double), value :: tn
      type(c_ptr), value    :: sunvec_y
      type(c_ptr), value    :: sunvec_f
      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), pointer :: yvec(:)
      real(c_double), pointer :: fvec(:)

      real(c_double) :: energy

      integer(c_long), parameter :: neq = 1

      ! get data arrays from SUNDIALS vectors
      call N_VGetData_Serial(sunvec_y, neq, yvec)
      call N_VGetData_Serial(sunvec_f, neq, fvec)

      call f_rhs(1, tn, yvec(1), energy, 0.0, 0)

      fvec(1) = energy

      ierr = 0
    end function RhsFn
end module rhs
