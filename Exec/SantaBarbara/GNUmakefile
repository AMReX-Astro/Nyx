# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../amrex

HPGMG_DIR ?= ../../Util/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = intel  # gnu
USE_OMP = FALSE
USE_MPI = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = FALSE
DEBUG     = FALSE

USE_HPGMG = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
