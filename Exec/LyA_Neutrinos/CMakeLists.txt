set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files inputs inputs.no_hydro_64 inputs.rt ascent_actions.yaml )
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)

#set(_props AGN ) #must set Nyx_AGN to true for this case, but AGN is broken
set(_props)
if (Nyx_NEUTRINO_PARTICLES)
   list(APPEND _props NEUTRINO_PARTICLES)
endif ()
if (Nyx_NEUTRINO_DARK_PARTICLES)
   list(APPEND _props NEUTRINO_DARK_PARTICLES)
endif ()

nyx_setup_executable(_sources _input_files
   PROPERTIES ${_props} )

unset(_sources)
unset(_input_files)
