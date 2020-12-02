set(_sources     Prob.H)
set(_input_files inputs inputs.no_hydro_64 inputs.rt ascent_actions.yaml )
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)

# must set Nyx_AGN to true for this case
option( Nyx_NEUTRINO_PARTICLES "" YES)
option( Nyx_NEUTRINO_DARK_PARTICLES "" YES)

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
