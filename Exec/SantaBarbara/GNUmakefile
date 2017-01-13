# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= /project/projectdirs/nyx/src/amrex

HPGMG_DIR ?= /global/homes/f/friesen/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# Reeber
BOOST_DIR ?= /opt/local
REEBER_HOME ?= $(HOME)/devel/Reeber/dev

# compilation options
COMP    = gcc
FCOMP   = gfortran
USE_OMP = FALSE
USE_MPI = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE

USE_HPGMG = FALSE
HPGMG_FCYCLES = TRUE
HPGMG_POST_F_CYCLE_TYPE = V
HPGMG_HELMHOLTZ = FALSE
HPGMG_STENCIL_VARIABLE_COEFFICIENT = FALSE
HPGMG_USE_SUBCOMM = TRUE
HPGMG_BOTTOM_SOLVER= BICGSTAB
HPGMG_SMOOTHER = GSRB

# Use single precision for particles?
#
# If you set this to be TRUE be sure to do a "make realclean"
# to make sure all your code is consistent.
#
USE_SINGLE_PRECISION_PARTICLES = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
