# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
