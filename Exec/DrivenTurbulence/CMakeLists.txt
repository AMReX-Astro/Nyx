set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files inputs inputs.regtest)

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
