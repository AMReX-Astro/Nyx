set(_sources     Prob.H)
set(_input_files inputs inputs.rt )
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)

nyx_setup_executable(_sources _input_files PROPERTIES AGN)

unset(_sources)
unset(_input_files)
