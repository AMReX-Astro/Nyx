subroutine ext_src_hc(lo, hi, old_state, os_l1, os_l2, os_l3, os_h1, os_h2, os_h3, &
                              new_state, ns_l1, ns_l2, ns_l3, ns_h1, ns_h2, ns_h3, &
                              old_diag , od_l1, od_l2, od_l3, od_h1, od_h2, od_h3, &
                              new_diag , nd_l1, nd_l2, nd_l3, nd_h1, nd_h2, nd_h3, &
                      src, src_l1, &
                      src_l2, src_l3, src_h1, src_h2, src_h3, problo, dx, time, z, dt)
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   old_state_* : double arrays
!       The state vars
!   new_state_* : double arrays
!       The state vars
!   src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   problo : double array (3)
!       The low corner of the entire domain
!   dx : double array (3)
!       The cell size of this level.
!   time : double
!       The current time, in Mpc km^-1 s ~ 10^12 yr.
!      z : double
!       The current z = 1 / a - 1
!   dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   src : double array (dims) @todo
!       @todo
!
    use meth_params_module, only : NVAR, UEDEN, UEINT, heat_cool_type
    use fundamental_constants_module
    use atomic_rates_module, only: interp_to_this_z

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: os_l1, os_l2, os_l3
    integer, intent(in) :: os_h1, os_h2, os_h3
    integer, intent(in) :: ns_l1, ns_l2, ns_l3
    integer, intent(in) :: ns_h1, ns_h2, ns_h3
    integer, intent(in) :: od_l1, od_l2, od_l3
    integer, intent(in) :: od_h1, od_h2, od_h3
    integer, intent(in) :: nd_l1, nd_l2, nd_l3
    integer, intent(in) :: nd_h1, nd_h2, nd_h3
    integer, intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    double precision, intent(in) :: old_state(os_l1:os_h1, &
                                              os_l2:os_h2, &
                                              os_l3:os_h3, NVAR)
    double precision, intent(in) :: new_state(ns_l1:ns_h1, &
                                              ns_l2:ns_h2, &
                                              ns_l3:ns_h3, NVAR)
    double precision, intent(in) :: old_diag (od_l1:od_h1, &
                                              od_l2:od_h2, &
                                              od_l3:od_h3, 2)
    double precision, intent(in) :: new_diag (nd_l1:nd_h1, &
                                              nd_l2:nd_h2, &
                                              nd_l3:nd_h3, 2)
    double precision, intent(in) :: problo(3), dx(3), z, dt, time

    double precision, intent(out) :: src(src_l1:src_h1, src_l2:src_h2, &
                                         src_l3:src_h3, NVAR)

    double precision, allocatable :: tmp_state(:,:,:,:)

    integer          :: i, j, k
    integer          :: src_lo(3),src_hi(3)
    integer          :: max_iter, min_iter
    double precision :: a, half_dt

    ! Make a copy of the state so we can evolve it then throw it away
    allocate(tmp_state(ns_l1:ns_h1,ns_l2:ns_h2,ns_l3:ns_h3,NVAR))
    tmp_state(:,:,:,:) = new_state(:,:,:,:)

     a = 1.d0 / (1.d0+z)

    ! Note that when we call this routine to compute the "old" source,
    !      both "old_state" and "new_state" are acutally the "old" state.
    ! When we call this routine to compute the "new" source,
    !      both "old_state" is in fact the "old" state and
    !           "new_state" is in fact the "new" state

    src_lo(1) = src_l1
    src_lo(2) = src_l2
    src_lo(3) = src_l3
    src_hi(1) = src_h1
    src_hi(2) = src_h2
    src_hi(3) = src_h3

    call interp_to_this_z(z)

    half_dt = 0.5d0 * dt
    if (heat_cool_type .eq. 1) then
        call integrate_state_hc(src_lo,src_hi,tmp_state,ns_l1,ns_l2,ns_l3, ns_h1,ns_h2,ns_h3, &
                                              new_diag ,nd_l1,nd_l2,nd_l3, nd_h1,nd_h2,nd_h3, &
                                a,half_dt,min_iter,max_iter)
    else if (heat_cool_type .eq. 3) then
        call integrate_state_vode(src_lo,src_hi,tmp_state,ns_l1,ns_l2,ns_l3, ns_h1,ns_h2,ns_h3, &
                                                new_diag ,nd_l1,nd_l2,nd_l3, nd_h1,nd_h2,nd_h3, &
                                  a,half_dt,min_iter,max_iter)
    endif
    do k = src_l3, src_h3
        do j = src_l2, src_h2
            do i = src_l1, src_h1
                  src(i,j,k,UEINT) = (tmp_state(i,j,k,UEINT) - new_state(i,j,k,UEINT)) * a / half_dt
                  src(i,j,k,UEDEN) = src(i,j,k,UEINT)
            end do
        end do
    end do

end subroutine ext_src_hc
