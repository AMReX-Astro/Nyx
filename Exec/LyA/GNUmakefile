# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

CVODE_LIB_DIR ?= $(CVODE_LIB)

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = FALSE
USE_OMP = FALSE

USE_SENSEI_INSITU = FALSE
USE_HENSON = FALSE

USE_ASCENT_INSITU=TRUE 
ASCENT_HOME=/project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/
include $(ASCENT_HOME)/share/ascent/ascent_config.mk
ASCENT_DIR = /project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/
#CONDUIT_DIR = /project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/
USE_CONDUIT = TRUE
USE_ASCENT = TRUE

PROFILE       = TRUE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = FALSE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE
USE_SUNDIALS_3x4x = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

