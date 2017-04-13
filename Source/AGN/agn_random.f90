! Random number generator wrapper module
!
! Paul Ricker (2/2013)

module agn_random

!==============================================================================

    use amrex_fort_module, only : rt => amrex_real
    implicit none

contains

!==============================================================================

! init_random_seed: Example random number seed initialization routine from 
!                   gfortran docs

subroutine init_random_seed()

integer                            :: i, n, clock
integer, dimension(:), allocatable :: seed

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed)

deallocate(seed)
end subroutine init_random_seed

!------------------------------------------------------------------------------

! random_unif: Return uniform random deviate in [0, 1).

real(rt) function random_unif()

real(rt) :: x

call random_number(x)

random_unif = x

return
end function random_unif

!==============================================================================

end module agn_random
