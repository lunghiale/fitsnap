export LAMMPS_DIR=/home/Alessandro/Documents/Sources/lammps-16Mar18

mpic++ -I ${LAMMPS_DIR}/src -c LAMMPS-wrapper.cpp 
mpif90 -c LAMMPS.F90

mpif90 -c  *.f90 -g
mpif90 -O3 -L ${LAMMPS_DIR}/src/  *.o -llammps_mpi -lmpi_cxx -lstdc++ -lm -g -llapack -o fitsnap.x
###mpif90 -O3 -L ${LAMMPS_DIR}/lib/reax -L ${LAMMPS_DIR}/src/  *.o -llammps_mpi -lmpi_cxx -lstdc++ -lm -g -llapack -o fitsnap.x
