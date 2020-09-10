# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

CVODE_LIB_DIR ?= $(CVODE_LIB)
BOOST_INLUDE_DIR := ${OLCF_BOOST_ROOT}/include/boost

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = TRUE
USE_OMP = FALSE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = FALSE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = TRUE
USE_CVODE = FALSE
USE_CONST_SPECIES = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
