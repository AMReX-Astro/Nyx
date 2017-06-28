module fcvode_extras

  implicit none

  contains

    subroutine fcvode_wrapper(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, T_out, ne_out, e_out)

        use amrex_fort_module, only : rt => amrex_real
        use vode_aux_module, only: rho_vode, T_vode, ne_vode
        use cvode_interface
        use fnvector_serial
        use, intrinsic :: iso_c_binding

        implicit none

        real(rt), intent(in   ) :: dt
        real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
        type(c_ptr), value :: cvmem
        type(c_ptr), value :: sunvec_y
        real(rt), intent(  out) ::         T_out,ne_out,e_out

        real(c_double) :: atol, rtol
        real(c_double) :: time, tout
        integer(c_long), intent(in) :: neq
        real(c_double), pointer, intent(in) :: yvec(:)

        integer(c_int) :: ierr

        real(c_double) :: t_soln

        T_vode   = T_in
        ne_vode  = ne_in
        rho_vode = rho_in

        ! Initialize the integration time
        time = 0.d0

        ! We will integrate "e" in time. 
        yvec(1) = e_in

        ! Set the tolerances.  
        atol = 1.d-4 * e_in
        rtol = 1.d-4

        ierr = FCVodeReInit(cvmem, time, sunvec_y)
        ierr = FCVodeSStolerances(CVmem, rtol, atol)

        ierr = FCVode(CVmem, dt, sunvec_y, time, CV_NORMAL)

        e_out  = yvec(1)
        T_out  = T_vode
        ne_out = ne_vode

    end subroutine fcvode_wrapper

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

end module fcvode_extras
