set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files inputs inputs.summit inputs_gimlet_in_transit.dsc)
list(APPEND _input_files inputs_nohydro.rt inputs.rt inputs.rt.garuda)
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)
list(APPEND _input_files ascent_actions_nohydro.yaml ascent_actions.yaml ascent_actions_slicefile.yaml)

# setup_nyx_executable(_sources _input_files)
nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
