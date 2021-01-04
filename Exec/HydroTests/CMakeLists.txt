set(_sources     Prob.cpp Prob_error.cpp)
set(_input_files
   inputs.3d.sph.testsuite.sedov
   inputs.regtest.sod
   inputs.regtest.test3
   inputs-sod-x
   inputs-sod-z
   inputs-test2-y
   inputs-test3-x
   inputs-test3-z
   inputs.3d.sph.sedov
   inputs.regtest.sedov
   inputs.regtest.test2
   inputs.restarttest.test2
   inputs-sod-y
   inputs-test2-x
   inputs-test2-z
   inputs-test3-y
   )

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
