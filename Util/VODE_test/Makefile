all: fmain.exe fmain_vode.exe


# Set SUNDIALS_INSTALL as an environment variable for a sundials installation 
# of a version older than 3.0 which has built cvode

prefix       = ${SUNDIALS_INSTALL}
exec_prefix  = ${SUNDIALS_INSTALL}
includedir   = ${SUNDIALS_INSTALL}/include
libdir       = ${SUNDIALS_INSTALL}/lib

F90         = /usr/bin/gfortran
F90FLAGS    = -O3
F90_LDFLAGS = -O3
F90_LIBS    =  -lm /usr/lib/x86_64-linux-gnu/librt.so

LINKFLAGS = -Wl,-rpath,${SUNDIALS_INSTALL}/lib

EXAMPLES = atomic_rates comoving_nd comoving_params constants_cosmo cvode_interface eos_hc eos_params fcvode_extras fnvector_serial f_rhs meth_params misc_params reion_aux_module vode_aux integrate_state_vode_3d ../BLAS/daxpy ../BLAS/ddot  ../BLAS/dscal  ../BLAS/idamax  ../BLAS/dcopy ../VODE/dewset ../VODE/dgesl ../VODE/dvhin  ../VODE/dvjac   ../VODE/dvnlsd ../VODE/dvode ../VODE/dvsol ../VODE/dvstep ../VODE/ixsav ../VODE/xsetf ../VODE/dacopy ../VODE/dgbfa ../VODE/dgefa ../VODE/dumach ../VODE/dvindy ../VODE/dvjust ../VODE/dvnorm ../VODE/dvset ../VODE/dvsrco ../VODE/iumach ../VODE/xerrwd ../VODE/xsetun ../VODE/dgbsl

OBJECTS  = ${EXAMPLES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .f90 .f

 %.o: %.mod

.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $< -o $(*F).o

.f.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $< -o $(@D)/$(*F).o

# -----------------------------------------------------------------------------------------

all: fmain_vode.exe fmain.exe

fmain_vode.exe: ${OBJECTS} fmain_vode.o
	${F90} ${F90_LDFLAGS} -o fmain_vode.exe -I${includedir} ${OBJECTS} fmain_vode.o ${F90_LIBS} -L${libdir} -lsundials_fcvode -lsundials_cvode -lsundials_fcvode -lsundials_fnvecserial -lsundials_nvecserial
	echo "${F90} ${F90_LDFLAGS} -o fmain_vode.exe -I${includedir} ${OBJECTS} fmain_vode.o ${F90_LIBS} -L${libdir} -lsundials_fcvode -lsundials_cvode -lsundials_fcvode -lsundials_fnvecserial -lsundials_nvecserial" 

fmain.exe: ${OBJECTS} fmain.o
	${F90} ${F90_LDFLAGS} -o fmain.exe -I${includedir} ${OBJECTS} fmain.o ${F90_LIBS} -L${libdir} -lsundials_fcvode -lsundials_cvode -lsundials_fcvode -lsundials_fnvecserial -lsundials_nvecserial
	echo "${F90} ${F90_LDFLAGS} -o fmain_vode.exe -I${includedir} ${OBJECTS} fmain_vode.o ${F90_LIBS} -L${libdir} -lsundials_fcvode -lsundials_cvode -lsundials_fcvode -lsundials_fnvecserial -lsundials_nvecserial" 

clean:
	rm -f ${OBJECTS}

