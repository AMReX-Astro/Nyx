! This module initializes and sets the Fourier modes for the computation of the physical force field

module forcing_spect_module

  use amrex_fort_module, only : bl_spacedim, rt => amrex_real
  use bl_types

  implicit none

  integer :: num_modes
  integer, allocatable :: wavevectors(:,:)
  real(rt), allocatable :: modes_even(:,:), modes_odd(:,:)

  public fort_alloc_spect, fort_set_wavevector, fort_set_modes

contains

  ! allocate arrays
  subroutine fort_alloc_spect(length) &
             bind(C, name="fort_alloc_spect")

    integer,  intent(in) :: length

    integer alloc

    if (length > 1) then
       num_modes = length

       allocate (modes_even(num_modes,bl_spacedim), modes_odd(num_modes,bl_spacedim), &
                 wavevectors(bl_spacedim,num_modes), STAT=alloc)
       if (alloc > 0) call bl_abort('failed to allocate arrays for forcing modes')
    else
       call bl_abort('number of forcing modes must be positive')
    end if

  end subroutine fort_alloc_spect

  ! set m-th nonzero wavevector of spectrum
  subroutine fort_set_wavevector(kvect, m) &
             bind(C, name="fort_set_wavevector")

    integer, intent(in) :: m
    integer, intent(in) :: kvect(bl_spacedim)

    if ((m.lt.0).or.(m.ge.num_modes)) then
	call bl_abort('invalid index of wavevector')
    end if

    wavevectors(:,m+1) = kvect(:) ! index conversion from C to Fortran
    print *, m+1, "k = ", wavevectors(1,m+1), wavevectors(2,m+1), wavevectors(3,m+1)
  end subroutine fort_set_wavevector

  ! set modes
  subroutine fort_set_modes(even, odd, length, comp) &
             bind(C, name="fort_set_modes")

    integer,  intent(in) :: length, comp
    real(rt), intent(in) :: even(length), odd(length)

    integer m

    if ((length.ne.num_modes).or.(comp.ge.bl_spacedim)) then
	call bl_abort('dimensions of input arrays do not match')
    end if

    modes_even(:,comp+1) = even(:)
    modes_odd(:,comp+1) = odd(:)

    do m = 1, num_modes
        print *, comp, m, modes_even(m,comp+1), modes_odd(m,comp+1)
    end do

  end subroutine fort_set_modes

end module forcing_spect_module

