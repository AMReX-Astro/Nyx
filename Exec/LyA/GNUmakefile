# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= /project/projectdirs/nyx/src/BoxLib
#BOXLIB_HOME = /home/vince/Development/BLMaster/BoxLib

HPGMG_DIR ?= /global/homes/f/friesen/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..


# compilation options
COMP    = gcc
FCOMP   = gfortran
USE_MPI = TRUE
USE_OMP = FALSE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE

USE_CXX11     = TRUE

DEFINES += -DBL_RANK_IN_PLT
DEFINES += -DDEBUG_AFAP

GIMLET = TRUE
REEBER = FALSE

GIMLET_DIR ?= /home/vince/Development/gimlet
# Gimlet needs FFTW MPI.
FFTW_INC ?= /usr/include
FFTW_DIR ?= /usr/lib/x86_64-linux-gnu

REEBER_HOME ?= $(HOME)/devel/Reeber/dev
# Reeber needs Boost (both headers and libraries) and diy2.
BOOST_INCLUDE_DIR ?= /usr/local/include
BOOST_LIB_DIR ?= /usr/local/lib
DIY2_INCLUDE_DIR ?= /usr/local/include

# Vince's fake Reeber (just usleeps) since he didn't have access to Reeber
# source while developing sidecars.
FAKE_REEBER = FALSE

USE_IN_TRANSIT = TRUE

USE_HPGMG = FALSE
HPGMG_FCYCLES = FALSE
HPGMG_POST_F_CYCLE_TYPE = V
HPGMG_HELMHOLTZ = FALSE
HPGMG_STENCIL_VARIABLE_COEFFICIENT = FALSE
HPGMG_USE_SUBCOMM = TRUE
HPGMG_BOTTOM_SOLVER= BICGSTAB
HPGMG_SMOOTHER = GSRB

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE

# Analysis
USE_IN_SITU = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

vpath %.a   $(LIBRARY_LOCATIONS)
