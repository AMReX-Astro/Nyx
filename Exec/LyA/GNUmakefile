# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME  ?= ../../../amrex
#ASCENT_HOME ?=/project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/
ASCENT_HOME ?=/gpfs/alpine/world-shared/csc340/software/ascent/0.5.3-pre/summit/cuda/gnu/ascent-install/

#This install location assumes modules: pgi/19.10 cuda/10.1.243
USE_FUSED = TRUE
ifeq ($(USE_FUSED),TRUE)
INSTALL_PREFIX ?= /gpfs/alpine/world-shared/ast160/software/sundials-5.3.0-fused-lasterror/instdir
else
INSTALL_PREFIX ?= /gpfs/alpine/world-shared/ast160/software/sundials/instdir
endif
#CVODE_LIB_DIR ?= $(CVODE_LIB)
CVODE_LIB_DIR ?= $(INSTALL_PREFIX)/lib64/
BOOST_INLUDE_DIR := ${OLCF_BOOST_ROOT}/include/boost

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_MPI = TRUE
# Use with Async IO
MPI_THREAD_MULTIPLE = FALSE
USE_OMP = TRUE
USE_CUDA ?= FALSE

USE_SENSEI_INSITU = FALSE
USE_HENSON = FALSE
REEBER = FALSE

USE_ASCENT_INSITU = FALSE

ifeq ($(USE_ASCENT_INSITU),TRUE)
  ASCENT_HOME ?= NOT_SET
  ifneq ($(ASCENT_HOME),NOT_SET)
    include $(ASCENT_HOME)/share/ascent/ascent_config.mk
  endif
  USE_CONDUIT = TRUE
  USE_ASCENT = TRUE
endif

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE
NEUTRINO_PARTICLES = FALSE
NEUTRINO_DARK_PARTICLES = FALSE

ifeq ($(USE_CUDA),TRUE)
	USE_SUNDIALS_3x4x = TRUE
	USE_SUNDIALS = TRUE
	COMP = pgi
	USE_OMP = FALSE
ifeq ($(REEBER), TRUE)
#6.4.0
#	NVCC_CCBIN=/sw/summit/gcc/6.4.0/bin/g++
#        BOOST_INCLUDE_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l6qu5nrynyuksm4ex7gswb6rsvvhsj42/include
#	LIBRARY_LOCATIONS += /sw/summit/gcc/6.4.0/lib64/
#        LIBRARIES += -Wl,-rpath, /sw/summit/gcc/6.4.0/lib64/
#7.4.0
	NVCC_CCBIN=/sw/summit/gcc/7.4.0/bin/g++
	BOOST_INCLUDE_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-7.4.0/boost-1.66.0-lmjl4oabmluitwt3kdwc6vks6riuuqqm/include
	LIBRARY_LOCATIONS += /sw/summit/gcc/7.4.0/lib64/
        LIBRARIES += -Wl,-rpath, /sw/summit/gcc/7.4.0/lib64/
endif
endif

USE_OWN_BCS = FALSE

Bpack := ./Make.package
Blocs := .

ifeq ($(USE_FUSED),TRUE)
ifeq ($(USE_CUDA),TRUE)
     LIBRARIES += -L$(CVODE_LIB_DIR) -lsundials_cvode_fused_cuda
endif
endif

include $(TOP)/Exec/Make.Nyx

