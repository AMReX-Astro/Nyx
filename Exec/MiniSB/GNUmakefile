# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= /project/projectdirs/nyx/src/amrex

HPGMG_DIR ?= $(HOME)/hpgmg/finite-volume

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = FALSE

# Analysis
REEBER = FALSE
#DEFINES += -DREEBER_PERSISTENT_INTEGRAL_TRACE_VTCS

USE_HPGMG = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
