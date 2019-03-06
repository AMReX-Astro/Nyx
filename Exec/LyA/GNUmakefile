# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

CVODE_LIB_DIR ?= $(CVODE_LIB)
CVODE_LIB_DIR = /ccs/home/jmsexton/jan_restart/sundials/instdir/lib64

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = pgi  # gnu
USE_MPI = TRUE
USE_OMP = FALSE

USE_HENSON = FALSE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE
USE_SUNDIALS_3x4x = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

