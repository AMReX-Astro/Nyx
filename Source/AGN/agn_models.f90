! Module providing routines implementing AGN feedback subgrid models.
!
! Paul Ricker (2/2013)

module agn_models

!==============================================================================

    use amrex_fort_module, only : rt => amrex_real

use agn_geometry
use agn_mesh
use agn_random
use agn_params

implicit none

public

contains

!==============================================================================

! agn_compute_acc_rate: Given current properties of a black hole and the
!                       surrounding gas, determine the mass accretion rate
!                       (Mdot) onto the black hole using a chosen subgrid model.
!
! arguments: agn          AGN data structure for black hole
!            Mdot         Variable to return mass accretion rate
!            t            Current time
!            m            Mesh data structure
!            rho          Gas density on mesh
!            P            Gas pressure on mesh
!            E            Gas specific total energy on mesh
!            accmodel     Choice of accretion model (see constant definitions)

subroutine agn_compute_acc_rate(agn, Mdot, time, &
                                state, s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                lo, hi, dx, accmodel)

use meth_params_module, only: gamma_const, gamma_minus_1, NVAR, URHO, UEINT

type(agn_particle), intent(in)   :: agn
real(rt), intent(in   )  :: time, dx(:)
integer         , intent(in   )  :: accmodel
real(rt), intent(  out)  :: Mdot

integer         ,intent(in   ) :: lo(:),hi(:)
integer         ,intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
real(rt),intent(in   ) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)

real(rt)              :: rho, P_over_rho
real(rt)              :: rhoinf, csinf, nh
real(rt), allocatable :: cs(:,:,:)

integer                       :: i,j,k

allocate(cs(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3))

do k = lo(3),hi(3)
do j = lo(2),hi(2)
do i = lo(1),hi(1)
    P_over_rho = gamma_minus_1*state(i,j,k,UEINT)/state(i,j,k,URHO)
    cs(i,j,k)  = sqrt(gamma_const*P_over_rho)
end do
end do
end do

select case (accmodel)

  case (AGN_ACC_ALPHABONDI)
     rhoinf = average_over_sphere(state(:,:,:,URHO), &
                                  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  agn%x, agn%y, agn%z, agn_ac_accradius, dx)
      csinf = average_over_sphere(cs, &
                                  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  agn%x, agn%y, agn%z, agn_ac_accradius, dx)

     Mdot   = agn_ac_ab_alpha * 4*pi*Gconst**2 * agn%m**2 * rhoinf/csinf**3

  case (AGN_ACC_BETABONDI)

     rhoinf = average_over_sphere(state(:,:,:,URHO), &
                                  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  agn%x, agn%y, agn%z, agn_ac_accradius, dx)
      csinf = average_over_sphere(cs, &
                                  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  agn%x, agn%y, agn%z, agn_ac_accradius, dx)

     Mdot   = agn_ac_ab_alpha * 4.d0*pi*Gconst**2 * agn%m**2 * rhoinf/csinf**3
     nh     = rhoinf/m_p
     if (nh > agn_ac_beta_nhstar) &
       Mdot = (nh/agn_ac_beta_nhstar)**agn_ac_ab_alpha * Mdot

  case (AGN_ACC_POPE)

    stop "agn_compute_acc_rate: Pope accretion case not yet implemented"

  case default

    stop "agn_compute_acc_rate: invalid accretion model requested"

end select

deallocate(cs)

Mdot = min(Mdot, agn_eddington_rate(agn%m))

return
end subroutine agn_compute_acc_rate

!------------------------------------------------------------------------------

! agn_compute_feedback: Given current properties of a black hole and the
!                       surrounding gas, along with an accretion rate, compute
!                       rates of mass, momentum, and energy deposition on the
!                       grid according to a chosen subgrid model.
!
! arguments: agn          AGN data structure for black hole
!            Mdot         Mass accretion rate
!            t, dt        Current time and (forward) timestep
!            m            Mesh data structure
!            rho          Gas density on mesh
!            P            Gas pressure on mesh
!            vx, vy, vz   Gas velocity components on mesh
!            E            Gas specific total energy on mesh
!            src          Array to receive rates of change
!            fbmodel      Choice of feedback model (see constant definitions)

subroutine agn_compute_feedback(agn, Mdot, t, dt, &
                                state, s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                src, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                dx, fbmodel)

use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN

type(agn_particle), intent(inout) :: agn
real(rt), intent(in   )  :: Mdot, t, dt
integer, intent(in)              :: fbmodel
real(rt), intent(in   )  :: dx(:)

integer         , intent(in   )  :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
real(rt), intent(in   )  :: state(s_l1:,s_l2:,s_l3:,:)

integer, intent(in)              :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
real(rt), intent(inout)  :: src(src_l1:,src_l2:,src_l3:,:)

real(rt)                 :: rholocal, Rdep, empty(1) = 0.
real(rt)                 :: rdisp, xbub, ybub, zbub, mudisp, &
                                    sintdisp, phidisp
real(rt)                 :: Rbub, Ebub, Mjet, Pjet, Ejet
integer                          :: ibh, jbh, kbh
real(rt)                 :: wparams(6), rjet, hjet
real(rt)                 :: d_zone

! Determine size of depletion region
d_zone = (dx(1)+dx(2)+dx(3))/3.d0

! TODO: replace with proper averaging
ibh = floor(agn%x/dx(1))
jbh = floor(agn%y/dx(2))
kbh = floor(agn%z/dx(3))

rholocal = state(ibh,jbh,kbh,URHO)

Rdep = (Mdot*dt * rholocal / (4./3.*pi * agn_fb_gasdepfac))**(1./3.)
Rdep = max(Rdep, agn_fb_depradius*d_zone)

! Remove gas that will be consumed by/processed by AGN

call map_sphere_onto_mesh(agn%x, agn%y, agn%z, Rdep, -Mdot, &
                          src(:,:,:,URHO), &
                          src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                          dx)

! Compute the spatially dependent feedback rates depending on the model chosen

! Various forms of "radio mode" feedback...

if (Mdot < agn_fb_twomode_threshold*agn_eddington_rate(agn%m)) then

  select case (fbmodel)

    case (AGN_FB_BUBBLE)

      if (agn%m-agn%mlast >= agn_fb_bub_delta*agn%mlast) then
        Ebub     = agn_fb_bub_mecheff*agn_fb_fbeff*c_light**2 * (agn%m-agn%mlast)
        Rbub     = agn_fb_bub_R0 * ( (Ebub/agn_fb_bub_E0) * &
                                     (agn_fb_bub_rho0/rholocal) )**0.2
        rdisp    = Rbub * random_unif()
        mudisp   = 2*random_unif() - 1
        phidisp  = 2*pi*random_unif()
        sintdisp = sqrt(1. - mudisp**2)
        xbub     = agn%x + rdisp*sintdisp*cos(phidisp)
        ybub     = agn%y + rdisp*sintdisp*sin(phidisp)
        zbub     = agn%z + rdisp*mudisp
        call map_sphere_onto_mesh(xbub, ybub, zbub, Rbub, Ebub/dt, &
                                  src(:,:,:,UEDEN), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  dx)
        agn%mlast = agn%m
      endif

    case (AGN_FB_JET)

      Mjet = agn_fb_jet_massload * Mdot * dt
      Pjet = sqrt(2.d0*agn_fb_fbeff) * Mdot * c_light * dt
      Ejet = agn_fb_fbeff * Mdot * c_light**2 * (1.d0-1.d0/agn_fb_jet_massload) * dt
      rjet = agn_fb_jet_cylradius
      hjet = agn_fb_jet_cylheight
      wparams = (/ rjet, hjet, agn%nx, agn%ny, agn%nz, rjet /)
      call map_cylinder_onto_mesh(agn%x, agn%y, agn%z, rjet, hjet, &
                                  agn%nx, agn%ny, agn%nz, Mjet/dt, &
                                  src(:,:,:,URHO), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  weight_cattaneo, wparams,dx)
      call map_cylinder_onto_mesh(agn%x, agn%y, agn%z, rjet, hjet, &
                                  agn%nx, agn%ny, agn%nz, &
                                  (Pjet*agn%nx + Mjet*agn%vx)/dt, &
                                  src(:,:,:,UMX), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  weight_cattaneo, wparams,dx)
      call map_cylinder_onto_mesh(agn%x, agn%y, agn%z, rjet, hjet, &
                                  agn%nx, agn%ny, agn%nz, &
                                  (Pjet*agn%ny + Mjet*agn%vy)/dt, &
                                  src(:,:,:,UMY), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  weight_cattaneo, wparams,dx)
      call map_cylinder_onto_mesh(agn%x, agn%y, agn%z, rjet, hjet, &
                                  agn%nx, agn%ny, agn%nz, &
                                  (Pjet*agn%nz + Mjet*agn%vz)/dt, &
                                  src(:,:,:,UMZ), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  weight_cattaneo, wparams,dx)
      call map_cylinder_onto_mesh(agn%x, agn%y, agn%z, rjet, hjet, &
                                  agn%nx, agn%ny, agn%nz, Ejet/dt, &
                                  src(:,:,:,UEDEN), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  weight_cattaneo, wparams,dx)

    case (AGN_FB_FIXBUB)

      if (agn%m-agn%mlast >= agn_fb_bub_delta*agn%mlast) then
        Ebub = agn_fb_bub_mecheff*agn_fb_fbeff*c_light**2 * (agn%m-agn%mlast)
        Rbub = agn_fb_fixbub_radius
        call map_sphere_onto_mesh(agn%x, agn%y, agn%z, Rbub, Ebub/dt, &
                                  src(:,:,:,UEDEN), &
                                  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                  dx)
        agn%mlast = agn%m
      endif

    case default

      stop "agn_compute_feedback: invalid feedback model requested"

  end select

! ... or "quasar mode" feedback if accretion rate is high enough

else

  Ebub = agn_fb_qsoeff*agn_fb_fbeff*c_light**2 * Mdot*dt
  Rbub = 8*d_zone
  call map_sphere_onto_mesh(agn%x, agn%y, agn%z, Rbub, Ebub/dt, &
                            src(:,:,:,UEDEN), &
                            src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                            dx)

endif

return
end subroutine agn_compute_feedback
!==============================================================================

end module agn_models
