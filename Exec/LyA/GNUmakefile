# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME  ?= ../../../amrex
#ASCENT_HOME ?=/project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/
ASCENT_HOME ?=/gpfs/alpine/world-shared/csc340/software/ascent/0.5.2-pre/summit/cuda/gnu/ascent-install

#This install location assumes modules: pgi/19.10 cuda/10.1.243
USE_FUSED = TRUE
ifeq ($(USE_FUSED),TRUE)
INSTALL_PREFIX ?= /gpfs/alpine/world-shared/ast160/software/sundials-5.3.0-fused-lasterror/instdir
else
INSTALL_PREFIX ?= /gpfs/alpine/world-shared/ast160/software/sundials/instdir
endif
#CVODE_LIB_DIR ?= $(CVODE_LIB)
CVODE_LIB_DIR ?= $(INSTALL_PREFIX)/lib64/
BOOST_INLUDE_DIR ?= ${OLCF_BOOST_ROOT}/include/boost
HDF5_HOME ?= ${OLCF_HDF5_ROOT}

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

USE_HDF5 ?= FALSE
USE_HDF5_ASYNC ?= FALSE

USE_ASCENT_INSITU = FALSE

ifeq ($(USE_ASCENT_INSITU),TRUE)
  ASCENT_HOME ?= NOT_SET
  ifneq ($(ASCENT_HOME),NOT_SET)
    include $(ASCENT_HOME)/share/ascent/ascent_config.mk
  endif
  USE_CONDUIT = TRUE
  USE_ASCENT = TRUE
endif

ifeq ($(USE_HDF5_ASYNC), TRUE)
HDF5_HOME= /gpfs/alpine/world-shared/csc300/gnu_build/hdf5/async_threadsafe/hdf5/
ABT_HOME = /gpfs/alpine/csc300/world-shared/gnu_build/argobots-1.0/install/
ASYNC_HOME = /gpfs/alpine/world-shared/csc300/gnu_build/async/src 
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
ifeq ($(NO_HYDRO),TRUE)
	USE_SUNDIALS_3x4x = FALSE
	USE_SUNDIALS = FALSE
	USE_FUSED = FALSE
else
	USE_SUNDIALS_3x4x = TRUE
	USE_SUNDIALS = TRUE
	COMP = pgi
endif
	USE_OMP = FALSE
ifeq ($(REEBER), TRUE)
#6.4.0
#	NVCC_CCBIN=/sw/summit/gcc/6.4.0/bin/g++
#        BOOST_INCLUDE_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l6qu5nrynyuksm4ex7gswb6rsvvhsj42/include
#	LIBRARY_LOCATIONS += /sw/summit/gcc/6.4.0/lib64/
#        LIBRARIES += -Wl,-rpath,/sw/summit/gcc/6.4.0/lib64/
#7.4.0
	NVCC_CCBIN=/sw/summit/gcc/7.4.0/bin/g++
	BOOST_INCLUDE_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-7.4.0/boost-1.66.0-lmjl4oabmluitwt3kdwc6vks6riuuqqm/include
	LIBRARY_LOCATIONS += /sw/summit/gcc/7.4.0/lib64/
        LIBRARIES += -Wl,-rpath,/sw/summit/gcc/7.4.0/lib64/
endif
ifeq ($(USE_ASCENT_INSITU), TRUE)
	NVCC_CCBIN=/sw/summit/gcc/6.4.0/bin/g++
#        BOOST_INCLUDE_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l6qu5nrynyuksm4ex7gswb6rsvvhsj42/include
	LIBRARY_LOCATIONS += /sw/summit/gcc/6.4.0/lib64/
        LIBRARIES += -Wl,-rpath,/sw/summit/gcc/6.4.0/lib64/
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

ifeq ($(USE_ASCENT_INSITU),TRUE)
#ifneq ($(USE_CUDA),TRUE)
     XTRAOBJS += $(ASCENT_DEVICE_OBJECT_MPI)
#     LIBRARIES += $(ASCENT_MPI_LIB_FLAGS)
      LIBRARIES += -L $(ASCENT_DIR)/lib \
                       -lascent_mpi \
                       -lascent_flow \
                       -lascent_lodepng \
                       -lrover_mpi $(DRAY_LINK_RPATH) $(DRAY_LIB_FLAGS) $(ASCENT_VTKH_MPI_LIB_FLAGS) $(ASCENT_VTKM_LIB_FLAGS) $(ASCENT_CONDUIT_MPI_LIB_FLAGS) $(ASCENT_MFEM_LIB_FLAGS) $(ASCENT_PYTHON_LIBS) $(ASCENT_OPENMP_LINK_FLAGS) -L $(ASCENT_CUDA_LIB_FLAGS)
#endif
endif

ifeq ($(USE_HDF5), TRUE)
DEFINES += -DAMREX_USE_HDF5
INCLUDE_LOCATIONS += $(HDF5_HOME)/include
LIBRARIES         += -L$(HDF5_HOME)/lib -lhdf5 -lz -ldl -Wl,-rpath,$(HDF5_HOME)/lib
endif

# To use HDF5 asynchronous I/O VOL connector, follow the instructions at https://bitbucket.hdfgroup.org/projects/HDF5VOL/repos/async/browse
ifeq ($(USE_HDF5_ASYNC), TRUE)
DEFINES   	  += -DAMREX_USE_HDF5_ASYNC -DAMREX_MPI_THREAD_MULTIPLE
INCLUDE_LOCATIONS += $(ABT_HOME)/include $(ASYNC_HOME)
LIBRARIES 	  += -L$(ABT_HOME)/lib -L$(ASYNC_HOME) -lh5async -labt -Wl,-rpath,$(ABT_HOME)/lib -Wl,-rpath,$(ASYNC_HOME)
endif

include $(TOP)/Exec/Make.Nyx

