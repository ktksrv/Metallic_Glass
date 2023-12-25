# Utilities
This directory contains the functions that are required for basic operations supporting the complete workflow

#### **_merge.py 
These functions merge the output from different CPUs to produce the required final output, when running the code on parallel settings on multi-CPU nodes

#### read_lammps_out.py
These functions can be used to read the atom-style data for a cluster from LAMMPS data and dump files

#### random_cluster_generator.py
It contains the OVITO-based code to extract the specified size of the STZ cluster from a random location inside the 32,800-atom metallic glass generated from LAMMPS. One can modify the bounds of the simulation cell to use it for their purpose.

#### misc.py
Contains miscellaneous functions used by other sections of the workflow
