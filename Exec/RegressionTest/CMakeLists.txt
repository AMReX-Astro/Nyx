set(_sources     Prob.H)
set(_input_files
   fixed_grids.mass
   inputs.2line
   inputs.2line.restarttest
   inputs.mass.nosub
   inputs.mass.sub
   Make.package
   particle_file.2line
   particle_file.mass
   )

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
