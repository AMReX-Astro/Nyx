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
DEBUG     = TRUE

BACKTRACE     = TRUE

USE_CXX11     = TRUE

DEFINES += -DBL_RANK_IN_PLT
DEFINES += -DDEBUG_AFAP

FFTW_INC ?= /usr/include
FFTW_DIR ?= /usr/lib/x86_64-linux-gnu

GIMLET_DIR ?= /home/vince/Development/gimlet
GIMLET = TRUE
REEBER = FALSE
ifeq ($(REEBER), TRUE)
  DEFINES += -DREEBER
  BOOST_DIR ?= /opt/local
  BOOST_INCLUDE_DIR ?= $(HOME)/Downloads/boost_1_60_0
  BOOST_LIB_DIR ?= $(HOME)/Downloads/boost_1_60_0/stage/lib
  REEBER_HOME ?= $(HOME)/devel/Reeber/dev
endif


FAKE_REEBER = FALSE
ifeq ($(FAKE_REEBER), TRUE)
  DEFINES += -DFAKE_REEBER
  BOOST_DIR =
endif

USE_IN_TRANSIT = TRUE
ifeq ($(USE_IN_TRANSIT), TRUE)
DEFINES += -DIN_TRANSIT
endif


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


