# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

SUNDIALS_ROOT ?= $(CVODE_LIB)

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = intel  # gnu
USE_MPI = TRUE
USE_OMP = TRUE
USE_CUDA = FALSE
USE_CXX_PROB = TRUE

ifeq ($(USE_CXX_PROB),TRUE)
  DEFINES += -DCXX_PROB
endif

PROFILE       = TRUE
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
USE_SUNDIALS = FALSE
USE_CONST_SPECIES = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

