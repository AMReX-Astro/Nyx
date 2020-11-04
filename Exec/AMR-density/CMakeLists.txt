set(_sources     Prob.H)
set(_input_files inputs inputs.cuda inputs_nohydro.rt inputs.rt)

# setup_nyx_executable(_sources _input_files)
nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
