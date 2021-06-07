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
USE_AGN = FALSE
ifeq ($(NO_HYDRO),TRUE)
USE_SDC = FALSE
USE_SUNDIALS = FALSE
USE_FUSED = FALSE
else
ifeq ($(USE_HEATCOOL),TRUE)
USE_SDC = TRUE
USE_SUNDIALS = TRUE
USE_FUSED ?= $(USE_CUDA)
else
USE_SDC = FALSE
USE_SUNDIALS = FALSE
USE_FUSED = FALSE
endif
endif
USE_CONST_SPECIES = TRUE
NEUTRINO_PARTICLES = TRUE
NEUTRINO_DARK_PARTICLES = TRUE

USE_OWN_BCS = FALSE

Bpack := ./Make.package
Blocs := .

ifeq ($(USE_FUSED),TRUE)
ifeq ($(USE_CUDA),TRUE)
     LIBRARIES += -L$(SUNDIALS_ROOT)/lib -lsundials_cvode_fused_cuda
else
     LIBRARIES += -L$(SUNDIALS_ROOT)/lib -lsundials_cvode_fused_stubs
endif
endif

include $(TOP)/Exec/Make.Nyx

