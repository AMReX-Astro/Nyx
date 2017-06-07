module fcvode_wrapper_mod
  implicit none
  contains

    subroutine fcvode_wrapper(dt, rho_in, T_in, ne_in, e_in, cvmem, sunvec_y, T_out, ne_out, e_out)

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

        ! Set the number of independent variables -- this should be just "e"
        integer(c_long), parameter :: NEQ = 1

        real(c_double), pointer :: y(:)

        real(c_double) :: atol, rtol
        real(c_double) :: time, tout

        integer(c_int) :: ierr

        real(c_double) :: t_soln

        T_vode   = T_in
        ne_vode  = ne_in
        rho_vode = rho_in

        ! Initialize the integration time
        time = 0.d0

        call N_VGetData_Serial(sunvec_y, neq, y)

        ! We will integrate "e" in time. 
        y(1) = e_in

        ! Set the tolerances.  
        atol = 1.d-4 * e_in
        rtol = 1.d-4

        ierr = FCVodeReInit(cvmem, time, sunvec_y)
        ierr = FCVodeSStolerances(CVmem, rtol, atol)

        ierr = FCVode(CVmem, dt, sunvec_y, time, CV_NORMAL)

        e_out  = y(1)
        T_out  = T_vode
        ne_out = ne_vode

    end subroutine fcvode_wrapper

end module fcvode_wrapper_mod
