
# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 400

nyx.ppm_type         = 1
nyx.ppm_reference    = 1
nyx.use_colglaz      = 0
nyx.corner_coupling  = 1

nyx.add_ext_src      = 0
nyx.heat_cool_type   = 3
nyx.strang_split     = 1

#This is 1e-8 times the lowest density in plt00000
nyx.small_dens = 5.162470e1

#This is 1e-5 times the constant temparature in plt00000
nyx.small_temp = 1.e-2

nyx.do_santa_barbara = 1
nyx.init_sb_vels     = 1
gravity.ml_tol = 1.e-10
gravity.sl_tol = 1.e-10

nyx.initial_z = 159.0
nyx.final_z = 0.0

#File written during the run: nstep | time | dt | redshift | a
amr.data_log = runlog
#amr.grid_log = grdlog

#This is how we restart from a checkpoint and write an ascii particle file
#Leave this commented out in cvs version
#amr.restart = chk00600
#max_step = 4
#particles.particle_output_file = particle_output

gravity.gravity_type = PoissonGrav
gravity.no_sync      = 1
gravity.no_composite = 1
mg.bottom_solver = 4

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic =  1     1     1
geometry.coord_sys   =  0

geometry.prob_lo     =  0     0     0

#Domain size in Mpc
geometry.prob_hi     =  28.49002849  28.49002849  28.49002849

amr.n_cell           = 1024 1024 1024
amr.max_grid_size    = 32
amr.blocking_factor  = 32

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
nyx.lo_bc       =  0   0   0
nyx.hi_bc       =  0   0   0

# WHICH PHYSICS
nyx.do_hydro = 1
nyx.do_grav  = 1

# COMOVING
nyx.comoving_OmM = 0.27
nyx.comoving_OmB = 0.045
nyx.comoving_h   = 0.71e0

# PARTICLES
nyx.do_dm_particles = 1

# >>>>>>>>>>>>>  PARTICLE INIT OPTIONS <<<<<<<<<<<<<<<<
#  "AsciiFile"        "Random"	    "Cosmological"
# >>>>>>>>>>>>>  PARTICLE INIT OPTIONS <<<<<<<<<<<<<<<<
nyx.particle_init_type = BinaryFile
nyx.binary_particle_file = 1024s_20mpc.nyx
particles.nparts_per_read = 2097152

# >>>>>>>>>>>>>  PARTICLE MOVE OPTIONS <<<<<<<<<<<<<<<<
#  "Gravitational"    "Random"
# >>>>>>>>>>>>>  PARTICLE MOVE OPTIONS <<<<<<<<<<<<<<<<
nyx.particle_move_type = Gravitational


# TIME STEP CONTROL
nyx.relative_max_change_a = 0.02    # max change in scale factor
particles.cfl             = 0.5     # 'cfl' for particles 
nyx.cfl                   = 0.5     # cfl number for hyperbolic system
nyx.init_shrink           = 1.0     # scale back initial timestep
nyx.change_max            = 1.1     # factor by which timestep can change
nyx.dt_cutoff             = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
nyx.print_fortran_warnings = 0
nyx.sum_interval      = -1      # timesteps between computing mass
nyx.v                 = 1       # verbosity in Nyx.cpp
gravity.v             = 1       # verbosity in Gravity.cpp
amr.v                 = 1       # verbosity in Amr.cpp
mg.v                  = 1       # verbosity in Amr.cpp
particles.v           = 1       # verbosity in Particle class

# REFINEMENT / REGRIDDING
amr.max_level          = 0        # maximum level number allowed
#amr.ref_ratio          = 2 2 2 2
#amr.regrid_int         = 4 4 4 4
#amr.n_error_buf        = 0 0 0 8
#amr.refine_grid_layout = 1
#amr.regrid_on_restart  = 1
#amr.blocking_factor    = 32
#amr.nosub              = 1

# CHECKPOINT FILES
amr.checkpoint_files_output = 1
amr.check_file      = chk
amr.check_int       = 20
amr.checkpoint_nfiles = 86

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt
amr.plot_int        = 20
amr.plot_nfiles     = 86
#nyx.plot_z_values   = 4.0 3.5 3.0 2.5 2.25 2.0

amr.plot_vars        = density
#amr.derive_plot_vars = particle_mass_density

#PROBIN FILENAME
amr.probin_file = probin

