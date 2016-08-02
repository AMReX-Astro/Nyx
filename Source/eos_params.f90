
! This module stores the runtime EOS species IF they are defined to be constants.  
! These parameter are initialized in set_eos_params().

module eos_params_module

  implicit none

  double precision, save ::  h_species
  double precision, save :: he_species

end module eos_params_module
