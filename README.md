# fitsnap

The code generates a SNAP potential in LAMMPS format from the fitting of reference energies and/or atomic forces.

# Citation

A. Lunghi, S. Sanvito, A unified picture of the covalent bond within quantumaccurate force fields: From organic molecules to metallic complexes’ reactivity. Sci. Adv. 5, eaaw2210 (2019)

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
     -compress <value of L2 norm constraint>
     -E0 <it can be either 1 or 0. Use 1 if all the training geometries contain the same number of atoms/atomic kind and 0 in any other case.>
     -refine <max new configurations> <MD temperature> <threshold for Gaussian metric>
     
     <lammps_input_file_heading>: Input file necessary to initialize a LAMMPS simulation.

     <dataset_info_file>: Contain the information about the training set used to fit SNAP.
     
     Line1: $(N: number of lammps' data files)
     Line2: $(name of data file 1) $(number of configurations) $(name of .xyz file containing all the geometries) $(if keyword -ener, name of .xyz file containing all the energies) $(if keyword -force, name of the file containing all the forces) $(weight)
     Line3: $(name of data file 2) $(number of configurations) $(name of .xyz file containing all the geometries) $(if keyword -ener, name of .xyz file containing all the energies) $(if keyword -force, name of the file containing all the forces) $(weight)
     ...For each LAMMPS data file
     
     The energies file should contain 1 column of numbers
     The structures file should contain a list of geometries in the xyz format
     The forces file should contain a list of forces in the xyz format but without the atomic labels column.
     

     <LAMMPS_potential_file>: Input file to set up a potantial in LAMMPS, including the SNAP potential that needs to be determined.

     <SNAP_potential_info_file>: File containing all the information on the SNAP potential to be determined.

     Line1: $(Radial Cutoff), $(Order of the bispectrum component 2J), $(Number of bispectrum components for the selected order), 0
     Line2: $(Number of atomic kinds)
     Line3: $(Label atomic kind 1) 1 $(element radius) $(weight) $(sigma gaussian metric)
     Line4: $(Label atomic kind 3) 2 $(element radius) $(weight) $(sigma gaussian metric)
     Line5: $(Label atomic kind 3) 3 $(element radius) $(weight) $(sigma gaussian metric)
     ...For each atomic kind
     
     Example: 
     
     Generation of SNAP potential for Water using 100 configurations and energies. 
     The SNAP potential is constructed for the Oxygen and Hydrogen atomic kinds using a bispectrum order of 2J=8 (56 parameters per kind).
     
     Command: "fitsnap.x -inp lammps_inp -datas data_file -pot_fix pot.fix -pot_fit pot.fit -ener"
     
     lammps_inp: 
          units real
          atom_style atomic
          pair_style snap
          
     data_file:
          1
          data.Water 100 Water.xyz Water.ener 1.0
     
     pot.fix:
          pair_coeff * * snapcoeff O H snapparam O H
          
     pot.fit:
          3.1000 8 56 0
          2
          O 1 0.5 1.0 1.0
          H 2 0.5 1.0 1.0
          
   








     
     


