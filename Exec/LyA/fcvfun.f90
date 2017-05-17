subroutine fcvfun(t, y, ydot, ipar, rpar, ier)
  implicit none

  double precision :: t, y(*), ydot(*), rpar(*), energy
  integer*8 :: ipar(*)
  integer :: ier

  call f_rhs(1, t, y(1), energy, rpar(1), ipar(1))

  ydot(1) = energy

  ier = 0
end subroutine fcvfun
