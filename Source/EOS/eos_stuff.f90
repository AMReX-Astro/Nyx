!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
!
!   NOTE: in the helmholtz EOS, we use Abar, which is the weighted
!   ion mass.  Abar does not include any free electron contribution.
!
! There are 2 ways to compute mu, depending on whether the gas is
! completely ionized or not.
!
! For a neutral gas, the mean molecular weight is:
!
!   1/mu = sum_k { X_k / A_k }
!
! For a completely ionized gas, the mean molecular weight is:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! At the moment, we will select between these 2 ionization extremes
! (completely neutral vs. completely ionized) by a hard-coded
! parameter, eos_assume_neutral
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module eos_module

  use bl_types
  use bl_space
  use bl_constants_module, only: M_PI, ONE
  use network, only: nspec, aion, zion
  use atomic_rates_module, only: XHYDROGEN

  implicit none

  private

  integer, parameter, private :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter, private :: eos_input_re = 5   ! density, internal energy are inputs

  logical, parameter, private :: eos_assume_neutral = .false.

  private nspec, aion, zion

  public eos_init_small_pres, nyx_eos_T_given_Re, nyx_eos_S_given_Re, &
         nyx_eos_soundspeed, nyx_eos_given_RT, eos

contains

  !---------------------------------------------------------------------------
  ! Nyx interfaces 
  !---------------------------------------------------------------------------

  subroutine eos_init_small_pres(R, T, Ne, P, comoving_a)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: P
     real(kind=dp_t), intent(in   ) :: R, T, Ne
     real(kind=dp_t), intent(in   ) :: comoving_a

     ! Local variables
     logical :: do_diag

     real(kind=dp_t) :: xn_eos(nspec)
     real(kind=dp_t) :: temp_eos
     real(kind=dp_t) :: den_eos
     real(kind=dp_t) :: e_eos
     real(kind=dp_t) :: p_eos
     real(kind=dp_t) :: cv_eos
     real(kind=dp_t) :: dpdt_eos
     real(kind=dp_t) :: dpdr_eos
     real(kind=dp_t) :: dedt_eos
     real(kind=dp_t) ::    s_eos
     real(kind=dp_t) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = (comoving_a*comoving_a*comoving_a)

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_rt, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Pressure must be converted from proper to comoving coordinates
    P  = p_eos * comoving_a_cubed

  end subroutine eos_init_small_pres

  subroutine nyx_eos_soundspeed(c, R, e)

     use meth_params_module, only: gamma_const, gamma_minus_1

     ! In/out variables
     real(kind=dp_t), intent(in   ) :: R, e
     real(kind=dp_t), intent(  out) :: c

     ! Pressure
     real(kind=dp_t) :: P

     P = R * e * gamma_minus_1

     ! sound speed
     c = sqrt(gamma_const * P / R)

  end subroutine nyx_eos_soundspeed

  subroutine nyx_eos_T_given_Re(T, Ne, R, e, comoving_a, pt_index)

     ! In/out variables
     real(kind=dp_t),           intent(inout) :: T, Ne
     real(kind=dp_t),           intent(in   ) :: R, e
     real(kind=dp_t),           intent(in   ) :: comoving_a
     integer        , optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     real(kind=dp_t) :: xn_eos(nspec)
     real(kind=dp_t) :: temp_eos
     real(kind=dp_t) :: den_eos
     real(kind=dp_t) :: e_eos
     real(kind=dp_t) :: p_eos
     real(kind=dp_t) :: cv_eos
     real(kind=dp_t) :: dpdt_eos
     real(kind=dp_t) :: dpdr_eos
     real(kind=dp_t) :: dedt_eos
     real(kind=dp_t) ::    s_eos
     real(kind=dp_t) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = comoving_a**3

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T 
     e_eos = e 

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    T  = temp_eos

    Ne = 1.d0

  end subroutine nyx_eos_T_given_Re

  subroutine nyx_eos_S_given_Re(S, R, e, T, Ne, comoving_a, pt_index)

     implicit none

     ! In/out variables
     real(kind=dp_t),           intent(  out) :: S
     real(kind=dp_t),           intent(in   ) :: R, e, T, Ne
     real(kind=dp_t),           intent(in   ) :: comoving_a
     integer, optional,         intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     real(kind=dp_t) :: xn_eos(nspec)
     real(kind=dp_t) :: temp_eos
     real(kind=dp_t) :: den_eos
     real(kind=dp_t) :: e_eos
     real(kind=dp_t) :: p_eos
     real(kind=dp_t) :: cv_eos
     real(kind=dp_t) :: dpdt_eos
     real(kind=dp_t) :: dpdr_eos
     real(kind=dp_t) :: dedt_eos
     real(kind=dp_t) ::    s_eos

     do_diag = .false.

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / (comoving_a*comoving_a*comoving_a)

     temp_eos = T
     e_eos = e

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos,  &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Should this be weighted by comoving a??
    S  = s_eos 

  end subroutine nyx_eos_S_given_Re

  subroutine nyx_eos_given_RT(e, P, R, T, Ne, comoving_a, pt_index)

     ! In/out variables
     real(kind=dp_t),           intent(  out) :: e, P
     real(kind=dp_t),           intent(in   ) :: R, T, Ne
     real(kind=dp_t),           intent(in   ) :: comoving_a
     integer,         optional, intent(in   ) :: pt_index(:)


     ! Local variables
     logical :: do_diag
     
     real(kind=dp_t) :: xn_eos(nspec)
     real(kind=dp_t) :: temp_eos
     real(kind=dp_t) :: den_eos
     real(kind=dp_t) :: e_eos
     real(kind=dp_t) :: p_eos
     real(kind=dp_t) :: cv_eos
     real(kind=dp_t) :: dpdt_eos
     real(kind=dp_t) :: dpdr_eos
     real(kind=dp_t) :: dedt_eos
     real(kind=dp_t) ::    s_eos
     real(kind=dp_t) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = (comoving_a*comoving_a*comoving_a)

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T 

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_rt, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Pressure must be converted from proper to comoving coordinates
    P  = p_eos * comoving_a_cubed

    e  = e_eos

  end subroutine nyx_eos_given_RT

  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(input, dens, temp, &
                 xmass, &
                 pres, eint, c_v, &
                 dPdT, dPdR, dEdT, &
                 entropy, &
                 do_eos_diag, &
                 pt_index)

    use bl_error_module
    use fundamental_constants_module, only: k_B, n_A, hbar
    use meth_params_module, only: gamma_minus_1

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! xmass    -- the mass fractions of the individual isotopes
! pres     -- the pressure (dyn/cm**2)
! eint     -- the internal energy (erg/g)
! c_v      -- specific heat at constant volume
! ne       -- number density of electrons + positrons
! dPdT     -- d pressure/ d temperature
! dPdR     -- d pressure/ d density
! dEdT     -- d energy/ d temperature
! entropy  -- entropy (erg/g/K)  NOTE: presently the entropy expression is 
!             valid only for an ideal MONATOMIC gas (gamma = 5/3).
!
! input = 1 means temp, dens    , and xmass are inputs, return eint    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input

    real(kind=dp_t) :: dens, temp
    real(kind=dp_t) :: xmass(nspec)
    real(kind=dp_t) :: pres, eint
    real(kind=dp_t) :: c_v
    real(kind=dp_t) :: dPdT, dPdR, dedT
    real(kind=dp_t) :: entropy

    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    real(kind=dp_t) :: ymass(nspec)    
    real(kind=dp_t) :: mu
    real(kind=dp_t) :: sum_y
    real(kind=dp_t) :: m_nucleon_over_kB
    real(kind=dp_t) :: t1,t2,t3

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: n

    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------

    sum_y  = 0.d0

    if (eos_assume_neutral) then
       ! assume completely neutral atoms
       do n = 1, nspec
          ymass(n) = xmass(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    else
       ! assume completely ionized species
       do n = 1, nspec
          ymass(n) = xmass(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    endif

    mu = 1.d0/sum_y

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------

    ! These are both very large numbers so we pre-divide
    m_nucleon_over_kB = m_nucleon / k_B

    if (input .NE. eos_input_rt) then

       ! dens, energy, and xmass are inputs
       
       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       temp = eint * mu * m_nucleon_over_kB * gamma_minus_1

    endif

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    ! compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation
    
    pres = dens * temp / mu / m_nucleon_over_kB

    eint = pres/gamma_minus_1/dens

    ! entropy (per gram) of an ideal monoatomic gas (the Sactur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.

    t1 = (mu*m_nucleon);     t1 = t1*t1*sqrt(t1)
    t2 = (k_B*temp);         t2 = t2*sqrt(t2)
    t3 = (2*M_PI*hbar*hbar); t3 = t3*sqrt(t3)

    entropy = (1.d0 / (mu*m_nucleon_over_kB)) * (2.5d0 + log(t1/dens*t2/t3))

    ! compute the thermodynamic derivatives and specific heats
    dPdT = pres/temp
    dPdR = pres/dens
    dedT = eint/temp

    c_v = dedT

  end subroutine eos

end module eos_module
