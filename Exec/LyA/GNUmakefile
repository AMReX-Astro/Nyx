# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME = $(HOME)/BoxLib

HPGMG_DIR = $(HOME)/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc                                                                                             
TOP = ../..

# compilation options
COMP    = gcc
FCOMP   = gfortran
USE_MPI = TRUE
USE_OMP = TRUE

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

REEBER_HOME ?= $(HOME)/devel/Reeber/dev
# Reeber needs Boost (both headers and libraries) and diy2.
BOOST_INCLUDE_DIR ?= /usr/local/include
BOOST_LIB_DIR ?= /usr/local/lib
DIY2_INCLUDE_DIR ?= /usr/local/include

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

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

