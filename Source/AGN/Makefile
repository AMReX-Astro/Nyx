OBJECT = agnug.o agn_models.o agn_mesh.o agn_geometry.o agn_random.o

.SUFFIXES:

.SUFFIXES: .f90 .o

%.o : %.f90
	gfortran -fdefault-real-8 -c $<

agnug : $(OBJECT)
	gfortran -o agnug $(OBJECT)

clean :
	rm -f agnug *.o

# module dependencies

agn_models.o : agn_geometry.o agn_mesh.o agn_random.o

agn_mesh.o   : agn_geometry.o

agnug.o      : agn_models.o
