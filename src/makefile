
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#-------------------------------------------------------------------------

CFLAGS = -std=c99

OBJECTS = roaring.o utility.o material.o init.o MPF-mc.o

MPF-mc: $(OBJECTS) chkopts 
	-${CLINKER} -o MPF-mc.exe $(OBJECTS) ${PETSC_TS_LIB}

clean::
	rm -f *.o *.exe

run: 
	-@${MPIEXEC} -n $(N) ./$(FILE).exe $(ARGS)
