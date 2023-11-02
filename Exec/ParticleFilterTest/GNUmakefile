# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../subprojects/amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = TRUE
# Use with Async IO
MPI_THREAD_MULTIPLE = FALSE
USE_OMP = TRUE
USE_CUDA = FALSE

USE_SUNDIALS      = TRUE
USE_FORT_ODE = FALSE
SUNDIALS_ROOT ?=  ${CURDIR}/$(TOP)/../sundials/instdir

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

# physics
DIM      = 3
USE_HEATCOOL = TRUE
USE_SAVE_REACT = FALSE
USE_AGN = FALSE
ifeq ($(NO_HYDRO),TRUE)
USE_SDC = FALSE
USE_SUNDIALS = FALSE
USE_FUSED = FALSE
else
ifeq ($(USE_HEATCOOL),TRUE)
USE_SDC = TRUE
USE_SUNDIALS = TRUE
ifeq ($(USE_HIP),TRUE)
USE_FUSED ?= $(USE_HIP)
endif
USE_FUSED ?= $(USE_CUDA)
else
USE_SDC = FALSE
USE_SUNDIALS = FALSE
USE_FUSED = FALSE
endif
endif
USE_CONST_SPECIES = TRUE
NEUTRINO_PARTICLES = FALSE
NEUTRINO_DARK_PARTICLES = FALSE

USE_OWN_BCS = FALSE

# Halo finder
BOOST_INLUDE_DIR := ${OLCF_BOOST_ROOT}/include/boost
REEBER = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
