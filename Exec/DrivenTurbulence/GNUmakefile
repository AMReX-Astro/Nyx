# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../subprojects/amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gcc

# NOTE THAT THIS CODE IS NOT THREAD-SAFE RIGHT NOW -- DO NOT RUN WITH USE_OMP
USE_MPI = TRUE
USE_OMP = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3

USE_HEATCOOL = FALSE

USE_CONST_SPECIES = TRUE

USE_SDC = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

