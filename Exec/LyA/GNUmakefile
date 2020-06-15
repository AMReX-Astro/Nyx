# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME  ?= ../../../amrex
ASCENT_HOME ?=/project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/

#This install location assumes modules: pgi/19.10 cuda/10.1.243
USE_FUSED = TRUE
ifeq ($(USE_FUSED),TRUE)
INSTALL_PREFIX ?= /gpfs/alpine/proj-shared/ast160/software/sundials-5.3.0-fused-lasterror/instdir
else
INSTALL_PREFIX ?= /gpfs/alpine/proj-shared/ast160/software/sundials/instdir
endif
#CVODE_LIB_DIR ?= $(CVODE_LIB)
CVODE_LIB_DIR ?= $(INSTALL_PREFIX)/lib64/
BOOST_INLUDE_DIR := ${OLCF_BOOST_ROOT}/include/boost

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = TRUE
USE_OMP = TRUE
USE_CUDA = FALSE

USE_SENSEI_INSITU = FALSE
USE_HENSON = FALSE
REEBER = FALSE

USE_ASCENT_INSITU = FALSE

ifeq ($(USE_ASCENT_INSITU),TRUE)
  ASCENT_HOME ?= NOT_SET
  ifneq ($(ASCENT_HOME),NOT_SET)
    include $(ASCENT_HOME)/share/ascent/ascent_config.mk
  endif
  USE_CONDUIT = TRUE
  USE_ASCENT = TRUE
endif

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE
NEUTRINO_PARTICLES = FALSE
NEUTRINO_DARK_PARTICLES = FALSE

ifeq ($(USE_CUDA),TRUE)
	USE_SUNDIALS_3x4x = TRUE
	USE_SUNDIALS = TRUE
	COMP = pgi
endif

USE_OWN_BCS = FALSE

Bpack := ./Make.package
Blocs := .

ifeq ($(USE_FUSED),TRUE)
     LIBRARIES += -L$(CVODE_LIB_DIR) -lsundials_cvode_fused_cuda
endif

include $(TOP)/Exec/Make.Nyx

