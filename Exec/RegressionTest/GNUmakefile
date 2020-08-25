# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = FALSE

USE_CXX_PROB = TRUE

ifeq ($(USE_CXX_PROB),TRUE)
  DEFINES += -DCXX_PROB
endif

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
