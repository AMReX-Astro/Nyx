set(_sources     Prob.H)
set(_input_files inputs inputs.regtest)

# Internal option for this test
set(Nyx_FORCING  YES CACHE INTERNAL "")

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
