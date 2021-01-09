set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files inputs inputs.no_hydro_64 inputs.rt ascent_actions.yaml )
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)


if (NOT Nyx_NEUTRINO_PARTICLES)
   message( WARNING "\nNyx_NEUTRINOS_PARTICLES=OFF:"
      " LyA_Neutrinos may not work as expected")
endif ()

if (NOT Nyx_NEUTRINO_DARK_PARTICLES)
   message( WARNING "\nNyx_NEUTRINOS_DARK_PARTICLES=OFF:"
      " LyA_Neutrinos may not work as expected")
endif ()

nyx_setup_executable(_sources _input_files )

unset(_sources)
unset(_input_files)
