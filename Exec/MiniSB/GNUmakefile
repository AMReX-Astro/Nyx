# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= /project/projectdirs/nyx/src/BoxLib

HPGMG_DIR ?= $(HOME)/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# Reeber (Edison)
BOOST_DIR ?= /project/projectdirs/nyx/ghweber/boost-1.58.0-edison-gcc-4.9.2
DIY2_INCLUDE_DIR ?= /project/projectdirs/nyx/ghweber/diy2/include
REEBER_HOME ?= /project/projectdirs/nyx/ghweber/reeber2

# compilation options
COMP    = gcc
FCOMP   = gfortran

USE_MPI = TRUE
USE_OMP = FALSE

# Analysis
REEBER = FALSE
#DEFINES += -DREEBER_PERSISTENT_INTEGRAL_TRACE_VTCS

USE_HPGMG = FALSE
HPGMG_FCYCLES = TRUE
HPGMG_POST_F_CYCLE_TYPE = V
HPGMG_HELMHOLTZ = FALSE
HPGMG_STENCIL_VARIABLE_COEFFICIENT = FALSE
HPGMG_USE_SUBCOMM = TRUE
HPGMG_BOTTOM_SOLVER= BICGSTAB
HPGMG_SMOOTHER = GSRB

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
