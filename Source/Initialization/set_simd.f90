subroutine set_simd (simd_width_in) bind(C, name='set_simd')

   use misc_params, only: simd_width
   implicit none

   integer, intent(in) :: simd_width_in

   simd_width = simd_width_in

end subroutine set_simd
