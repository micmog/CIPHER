# A Phase-Field Simulation Tool for Microstructure Evolution in Multi-Component Alloy Systems

## What is this?

A phase-field research code for microstructure evolution in multi-component alloy systems. Features include:

- Local truccation error estimates and adaptive time stepping.
- Designed for large number of phases (10-10000), with constant computational complexity, i.e. computational cost is independent of the number of phases.
- Efficient grand-canonical-based phase-field implementation with direct use of CALPHAD thermodynamic databases for multi-component systems [1].
- MPI parallelization and scalable.

## Requirements

This software requires PETSc v3.11 and assumes PETSC_DIR and PETSC_ARCH environment variables have been set [2]. Also included is the the Roaring BitMap implementation from [3].

## Usage

To compile:
```bash
cd src
make clean
make MPF-mc
```

Running an example:
```bash
cd example/GrainBoundaryPrecipitate
./MPF-mc.exe -finaltime 100000.0 -mintimestep 1e-9 -maxtimestep 1e+3 -l 5 -ptol 1e-2 -ctol 1e-3 -kI 0.2 -geomfile GrainBoundaryPrecipitate.geom -interfacefile N3.interface -outfreq 10000 -outfile GrainBoundaryPrecipitate
```

## References

[1] Grand-canonical phase-field implementation: https://arxiv.org/abs/1906.10503  
[2] PETSc: https://www.mcs.anl.gov/petsc/  
[3] Roaring BitMap: https://github.com/RoaringBitmap/CRoaring

## License

You are free to use and redistribute this is a free software under the terms of the GNU General Public License v3.0.
 
