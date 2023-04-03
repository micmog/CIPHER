
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#-------------------------------------------------------------------------

CFLAGS = -std=c99

OBJECTS = typedef.o utility.o material.o init.o cipher.o

cipher: $(OBJECTS) 
	-${CLINKER} -o cipher.exe $(OBJECTS) ${PETSC_TS_LIB}

clean::
	rm -f *.o *.exe

