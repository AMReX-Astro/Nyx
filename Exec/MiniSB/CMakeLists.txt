set(_sources     Prob.H)
set(_input_files
   inputs.32
   inputs.32.plot_z
   inputs.32.ref
   inputs.analysis
   particle_file.small
   fixed_grids
   )

nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
