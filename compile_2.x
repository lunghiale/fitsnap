### LAMMPS must have been compiled as a library 
### The files LAMMPS.F90 LAMMS-wrapper.h LAMMPS-wrapper.cpp need to be copied from the lammps repository into 
### the fitsnap folder.
### The compilation scripts has to be run in the same folder as the fitsnap source code.

export LAMMPS_DIR=<path to the lammps distribution>

mpic++ -I ${LAMMPS_DIR}/src -c LAMMPS-wrapper.cpp 
mpif90 -c LAMMPS.F90

mpif90 -c  *.f90 -g
mpif90 -O3 -L ${LAMMPS_DIR}/src/  *.o -llammps_mpi -lmpi_cxx -lstdc++ -lm -g -llapack -o fitsnap.x
