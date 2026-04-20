# Before running this script,
# ensure that you have already cloned the CIPHER repo using
# git clone --recurse-submodules https://github.com/micmog/CIPHER.git
# This should create a `CIPHER` subdirectory in the current directory
# i.e. CIPHER_DIR (defined below).

export CIPHER_DIR=/mnt/eps01-rds/jf01-home01/shared/software/apps/cipher/CIPHER
export PETSC_DIR=$CIPHER_DIR/petsc
export PETSC_ARCH=cipher

module load mpi/gcc/openmpi/4.1.8-gcc-14.2.0

cd $PETSC_DIR

# Update PETSc to latest version
git pull origin main

./configure --download-metis --download-parmetis --download-chaco \
  --download-triangle --download-ctetgen --download-pragmatic \
  --download-eigen --download-hypre --download-ml --download-hdf5 \
  --download-zlib --download-yaml --download-p4est --with-pthread \
  --with-debugging=0  --download-fblaslapack=1
make all 2>&1 | tee make-all.log
make check 2>&1 | tee make-check.log

cd $CIPHER_DIR
make install 2>&1 | tee make-install.log
