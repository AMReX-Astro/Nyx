# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

HPGMG_DIR ?= ../../Util/hpgmg/finite-volume
CVODE_LIB_DIR ?= $(CVODE_LIB)

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = intel  # gnu
USE_MPI = TRUE
USE_OMP = TRUE

PROFILE       = TRUE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

USE_HPGMG = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

HDF5_HOME   = /global/homes/s/sbyna/software/hdf5develop/install
USE_HDF5 = TRUE
ifeq ($(USE_HDF5), TRUE)
DEFINES += -DBL_HDF5 
INCLUDE_LOCATIONS += $(HDF5_HOME)/include
LIBRARIES         += -lhdf5 -L$(HDF5_HOME)/lib
# Edison
#  XTRALIBS += -lhdf5 -L$(HDF5_HOME)/lib -L/usr/common/software/zlib/1.2.7/intel/lib -lz
endif
