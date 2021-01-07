set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files inputs inputs.cuda inputs_nohydro.rt inputs.rt)
list(APPEND _input_files  64sssss_20mpc.nyx TREECOOL_middle)

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
