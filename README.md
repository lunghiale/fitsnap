# fitsnap

The code generates a SNAP potential in LAMMPS format from the fitting of reference energies and/or atomic forces.

# Usage

The execution of the code with no command-line arguments prints out the available keywords. 
fitsnap.x 

     -inp <lammps_input_file_heading>
     -datas <dataset_info_file>
     -pot_fix <LAMMPS_potential_file>
     -pot_fit <SNAP_potential_info_file>
     -ener perform the fitting of reference energies
     -force perform the fitting of reference energies
     -skip_fit does not perform a fitting but simply test the potential
     -refine <max new configurations> <MD temperature> <threshold for Gaussian metric>


