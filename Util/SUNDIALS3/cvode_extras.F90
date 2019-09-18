module cvode_extras
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  contains

    subroutine ode_eos_setup(a,half_dt) &
      bind(C,name="fort_ode_eos_setup")
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, &
                                   gamma_minus_1
!    use amrex_constants_module, only: M_PI
    use eos_params_module
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: JH_vode, JHe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   T_zhi, T_zheii, inhomogeneous_on
    real(rt), intent(inout) :: a
    real(rt), intent(inout) :: half_dt
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
     real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)
    real(rt) :: rho_vode, T_vode, ne_vode
    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, half_dt)
    z_end = 1.0d0/a_end - 1.0d0

    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    if (inhomogeneous_on) then
       STOP "Do not currently support inhomogenous_on with box"
       !H_reion_z = diag_eos(i,j,k,ZHI_COMP)
       if (z .gt. H_reion_z) then
          JH_vode = 0
       else
          JH_vode = 1
       endif
    endif

    end subroutine ode_eos_setup

    AMREX_CUDA_FORT_DEVICE subroutine ode_eos_finalize(e_out, rpar, num_eq) &
      bind(C,name="fort_ode_eos_finalize")
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, &
                                   gamma_minus_1
    use amrex_constants_module, only: M_PI
    use eos_params_module
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT, nyx_eos_T_given_Re_device
    use fundamental_constants_module
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: JH_vode, JHe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   T_zhi, T_zheii, inhomogeneous_on
    integer, value :: num_eq
   real(rt), intent(inout) :: e_out(num_eq)
   real(rt), intent(inout) :: rpar(4*num_eq)

    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
     real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)
    real(rt) :: rho_vode, T_vode, ne_vode
    real(rt) :: a
!     attributes(managed) :: T_vode,Ne_vode

      T_vode=rpar(1)
      ne_vode=rpar(2)
      rho_vode=rpar(3)
      a=1/(rpar(4)+1)


      if (e_out(1) .lt. 0.d0) then
#ifndef AMREX_USE_CUDA
         !$OMP CRITICAL
         print *,'negative e exiting strang integration ',z, rho/mean_rhob, e_out
         !call flush(6)
         !$OMP END CRITICAL
#endif
         T_vode  = 10.0
         ne_vode = 0.0
         mu     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_vode)
         e_out  = T_vode / (gamma_minus_1 * mp_over_kB * mu)
         !                    call amrex_abort('bad e out of strang')
      end if
      
      ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
      call nyx_eos_T_given_Re_device(JH_vode, JHe_vode, T_vode, ne_vode, rho_vode, e_out(1), a, species)
      
      ! Instanteneous heating from reionization:
      T_H = 0.0d0
      if (inhomogeneous_on .or. flash_h) then
         if ((H_reion_z  .lt. z) .and. (H_reion_z  .ge. z_end)) T_H  = (1.0d0 - species(2))*max((T_zhi-T_vode), 0.0d0)
      endif
      
      T_He = 0.0d0
      if (flash_he) then
         if ((He_reion_z .lt. z) .and. (He_reion_z .ge. z_end)) T_He = (1.0d0 - species(5))*max((T_zheii-T_vode), 0.0d0)
      endif
      
      if ((T_H .gt. 0.0d0) .or. (T_He .gt. 0.0d0)) then
         T_vode = T_vode + T_H + T_He                            ! For simplicity, we assume
         ne_vode = 1.0d0 + YHELIUM                              !    completely ionized medium at
         if (T_He .gt. 0.0d0) ne_vode = ne_vode + YHELIUM        !    this point.  It's a very minor
         mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_vode)   !    detail compared to the overall approximation.
         e_out  = T_vode / (gamma_minus_1 * mp_over_kB * mu)
         call nyx_eos_T_given_Re_device(JH_vode, JHe_vode, T_vode, ne_vode, rho_vode, e_out(1), a, species)
      endif
      rpar(1)=T_vode
      rpar(2)=ne_vode
      rpar(3)=rho_vode
      rpar(4)=z_vode      

    end subroutine ode_eos_finalize

    
    AMREX_CUDA_FORT_DEVICE integer(c_int) function RhsFnReal(tn, yvec, fvec, rpar, neq) &
           result(ierr) bind(C,name='RhsFnReal')

      use, intrinsic :: iso_c_binding
      use f_kernel_rhs_dev

      implicit none


      real(rt), value :: tn
      integer(c_int), value :: neq
!      type(c_ptr), value    :: sunvec_y                                                                                                                                                                            
!      type(c_ptr), value    :: sunvec_f                                                                                                                                                                            
!      type(c_ptr), value    :: user_data                                                                                                                                                                           

      ! pointers to data in SUNDAILS vectors                                                                                                                                                                        
      real(rt) :: yvec(neq)
      real(rt) :: fvec(neq)
      real(rt), intent(inout) :: rpar(neq*4)
      real(rt) :: energy(neq)

!      print*, "r1", rpar(1)                                                                                                                                                                                        
!      print*, "r2", rpar(2)                                                                                                                                                                                        
 !     print*, rpar(3)                                                                                                                                                                                              
  !    print*, rpar(4)                                                                                                                                                                                              
      call f_rhs_rpar(tn, yvec, fvec, rpar)
!      print*, fvec
   !   print*, "after r1", rpar(1)                                                                                                                                                                                  
    !  print*, "after r2", rpar(2)                                                                                                                                                                                 
     ! print*, "after r3", rpar(3)                                                                                                                                                                                  
!      print*, "after r4", rpar(4)                                                                                                                                                                                  

      ierr = 0
    end function RhsFnReal

end module cvode_extras
