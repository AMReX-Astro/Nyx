# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../amrex
#AMREX_HOME = /home/vince/Development/BLMaster/amrex

HPGMG_DIR ?= ../../../hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = TRUE
USE_OMP = FALSE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

GIMLET_DIR ?= /home/vince/Development/gimlet
# Gimlet needs FFTW MPI.
FFTW_INC ?= /usr/include
FFTW_DIR ?= /usr/lib/x86_64-linux-gnu

REEBER_HOME ?= /project/projectdirs/nyx/ghweber/reeber2
# Reeber needs Boost (both headers and libraries) and diy2.
BOOST_INCLUDE_DIR ?= /project/projectdirs/nyx/ghweber/boost-1.61.0-noarch/include
DIY_INCLUDE_DIR ?= /project/projectdirs/nyx/ghweber/diy/include

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

USE_AGN = TRUE

# units
USE_CGS = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

