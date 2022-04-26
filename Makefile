FC=mpif90
FFLAGS=-fdefault-real-8

MACRO=

HDF5compile=-I/usr/include/hdf5/openmpi
HDF5link=-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5

INCLUDE+=$(HDF5compile)
LIB+=$(HDF5link)

# the name of the executable
TARGET=./spec_field_reader.x

# the source files
OBJS=poly.o spec_geometry.o spec_field.o spec_state.o spec_io.o main.o


all :	$(OBJS)
	$(FC) ${FFLAGS} ${MACRO} -o $(TARGET) $(MODULES) $(OBJS) $(LIB)

clean :
	rm *.o
	rm *.mod
	rm $(TARGET)

%.o : %.f90
	$(FC) ${FFLAGS} ${MACRO} ${INCLUDE} -c $<
