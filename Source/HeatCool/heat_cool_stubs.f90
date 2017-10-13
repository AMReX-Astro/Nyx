subroutine ext_src_hc(lo, hi, old_state, old_state_l1, old_state_l2, &
                      old_state_l3, old_state_h1, old_state_h2, old_state_h3, &
                      new_state, new_state_l1, new_state_l2, new_state_l3, &
                      new_state_h1, new_state_h2, new_state_h3, src, src_l1, &
                      src_l2, src_l3, src_h1, src_h2, src_h3, problo, dx, time, z, dt)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: old_state_l1, old_state_l2, old_state_l3
    integer, intent(in) :: old_state_h1, old_state_h2, old_state_h3
    integer, intent(in) :: new_state_l1, new_state_l2, new_state_l3
    integer, intent(in) :: new_state_h1, new_state_h2, new_state_h3
    integer, intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    real(rt), intent(in) :: old_state(old_state_l1:old_state_h1, &
                                              old_state_l2:old_state_h2, &
                                              old_state_l3:old_state_h3, NVAR)
    real(rt), intent(in) :: new_state(new_state_l1:new_state_h1, &
                                              new_state_l2:new_state_h2, &
                                              new_state_l3:new_state_h3, NVAR)
    real(rt), intent(in) :: problo(3), dx(3), z, dt, time
 
    real(rt), intent(out) :: src(src_l1:src_h1, src_l2:src_h2, &
                                         src_l3:src_h3, NVAR)

    src = 0.d0
end subroutine ext_src_hc

subroutine integrate_state(lo, hi, state, state_l1, state_l2, &
                           state_l3, state_h1, state_h2, state_h3, &
                           dx, time, a, half_dt, min_iter, max_iter) bind(C, name="integrate_state")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: state_l1, state_l2, state_l3
    integer, intent(in) :: state_h1, state_h2, state_h3
    real(rt), intent(inout) :: state(state_l1:state_h1, state_l2:state_h2, &
                                             state_l3:state_h3, NVAR)
    real(rt), intent(in   ) :: dx(3), time, a, half_dt
    integer,  intent(inout) :: min_iter, max_iter

end subroutine integrate_state
 
module adjust_heat_cool_module

  implicit none
 
  contains

    subroutine adjust_heat_cool(lo,hi, &
                                u_old,uo_l1,uo_l2,uo_l3,uo_h1,uo_h2,uo_h3, &
                                u_new,un_l1,un_l2,un_l3,un_h1,un_h2,un_h3, &
                                src_old, so_l1,so_l2,so_l3,so_h1,so_h2,so_h3, &
                                src_new, sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                a_old, a_new, dt)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: uo_l1, uo_l2, uo_l3, uo_h1, uo_h2, uo_h3
      integer          :: un_l1, un_l2, un_l3, un_h1, un_h2, un_h3
      integer          :: so_l1,so_l2,so_l3,so_h1,so_h2,so_h3
      integer          :: sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3
      real(rt) ::   u_old(uo_l1:uo_h1,uo_l2:uo_h2,uo_l3:uo_h3,NVAR)
      real(rt) ::   u_new(un_l1:un_h1,un_l2:un_h2,un_l3:un_h3,NVAR)
      real(rt) :: src_old(so_l1:so_h1,so_l2:so_h2,so_l3:so_h3,NVAR)
      real(rt) :: src_new(sn_l1:sn_h1,sn_l2:sn_h2,sn_l3:sn_h3,NVAR)
      real(rt) :: a_old, a_new, dt

    end subroutine adjust_heat_cool

end module adjust_heat_cool_module

! unused VODE stubs if we are not doing heating/cooling
module vode_aux_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt) :: z_vode
end module vode_aux_module
