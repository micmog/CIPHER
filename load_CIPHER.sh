#!/bin/bash

module purge
PATH=$BASE_PATH
LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
module load tools/env/proxy
module load mpi/intel-19.1/openmpi/4.1.1 
export PETSC_DIR=$CIPHER_DIR/petsc
export PETSC_ARCH=cipher
PATH=$PETSC_DIR/$PETSC_ARCH/bin:$CIPHER_DIR/bin:$PATH
LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
