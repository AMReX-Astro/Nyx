# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

HPGMG_DIR ?= ../../Util/hpgmg/finite-volume
CVODE_LIB_DIR ?= /project/projectdirs/ecpcosmo/Code/sundials/sundials-intel/lib


# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

DEFINES += -DAMREX_USE_CATALYST
# If we want to output image data when AMR refinement is disabled and number of blocks/grids equals number of MPI ranks
#DEFINES += -DCATALYST_IMAGEDATA

# these are for full Catalyst edition
# result from ./paraview-config --include vtkPVPythonCatalyst
INCLUDE_LOCATIONS += /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/include/paraview-5.4

XTRALIBS += -Wl,-rpath,/global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4:/global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVPythonCatalyst-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVCatalyst-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVServerManagerApplication-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVServerManagerCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVServerImplementationCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVClientServerCoreCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersParallel-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersExtraction-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersModeling-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersGeometry-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersProgrammable-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVVTKExtensionsCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPVCommon-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkClientServer-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkIOXMLParser-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersSources-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersGeneral-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkFiltersCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkIOImage-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libprotobuf.so /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkjsoncpp-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkPythonInterpreter-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/libpython2.7.so /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkWrappingPython27Core-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkParallelMPI-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkParallelCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkIOLegacy-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkIOCore-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonExecutionModel-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonDataModel-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonSystem-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtksys-pv5.4.so.1 -ldl /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonMisc-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonTransforms-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonMath-pv5.4.so.1 /global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4/libvtkCommonCore-pv5.4.so.1 -Wl,-rpath-link,/global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib/paraview-5.4:/global/homes/a/acbauer/CODE/ParaView/sbbuild/install/lib


# these are for non-rendering Catalyst edition
# result from ./paraview-config --include vtkPVPythonCatalyst
#INCLUDE_LOCATIONS += /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/include/paraview-5.4

#XTRALIBS +=   -Wl,-rpath,/global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4:/opt/python/17.06.1/lib /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVPythonCatalyst-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVCatalyst-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVServerManagerApplication-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVServerManager-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVServerImplementationCore-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVClientServer-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVExtensions-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkPVCore-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libprotobuf.so /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkjsoncpp-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkWrapping-pv5.4.so.1 /opt/python/17.06.1/lib/libpython2.7.so /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkWrappingPython27Core-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkParallel-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkFilters-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkIO-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtkCommon-pv5.4.so.1 /global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4/libvtksys-pv5.4.so.1 -ldl -Wl,-rpath-link,/global/cscratch1/sd/acbauer/ParaViewCatalyst/install2/lib/paraview-5.4

# compilation options
COMP    = intel  # gnu
USE_MPI = TRUE
USE_OMP = TRUE

PROFILE       = TRUE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

USE_HPGMG = TRUE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

