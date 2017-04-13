
      subroutine get_grav_const(Gconst_out)

         use amrex_fort_module, only : rt => amrex_real
         use fundamental_constants_module, only: Gconst

         real(rt) :: Gconst_out

         Gconst_out = Gconst

      end subroutine get_grav_const

