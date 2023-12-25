# Bonding Calculations and Yield Estimation
This directory contains the code used for bond calculations on Vitreloy-105 using a statistical mechanics-based approach to estimate the status of bonding between any given pair of atoms at a given applied stress state.


## Method for Bond Cut-offs
1. Nearest -Neighbour distance for each of the atoms (with the same or another species, such as Pd-Ni, Cu-Cu, etc) for a large cluster (~1000 atoms) is calculated across all the deformation states applied using the deformation gradient approach
2. A probability distribution is then generated for the obtained nearest-neighbor distances for each atom with the same or different species
3. The minima of this distribution for various atoms are again plotted as a single distribution, and the maximum of this distribution is taken as the bond cutoff for the pair of species being investigated.

The code contains, the bond calculations done for 15 possible bond lengths in Vitreloy-105 glass, it's a general script that can used to analyze any multi-component system, by adding or removing required interactions based on the type of glass under investigation. Only the section defining the atom types and potential has to be changed to analyze a different system. The function for this functionality is explained below.

### Normal_Stress_Matrix_bond_len_multi(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc)
* **no_atoms** - defines the STZ cluster size
* **max_alpha, max_beta** - defines the upper bound of the tunable parameters of the deformation gradient
* **total_disp_steps** - defines the density of mesh in the displacement space of α and β ( higher means finer displacement grid on which stresses will be computed)
* **iteration** - defines the STZ cluster number of a given size
* **total_proc,proc_per_task,current_proc** - defines the parameters for leveraging parallel CPU architecture as all the calculations are embarrassingly parallel. 

After using the merging script (in the ../Utilities directory) to merge the output from various CPUs, the output will be a (no_atoms x (total_disp_steps * total_disp_steps)) matrix _M_, whose _M<sub>ij</sub>_ element represents _i<sup>th</sup>_ atom's nearest neighbor distance with same species or different species at the j<sup>th</sup> step of deformation out of total (total_disp_steps * total_disp_steps) steps

## Method for Yield estimation of the STZ cluster
Once the bond cutoff for all the possible pairing of speicies in the multi-component glass is specified, one can define a paramter call bond status changes at a given strain step, in comparision to the last one. This parameter basically tracks the summation of bond breakage or new bond formation inside the cluster at a given strain step. Then using the cummulative value of all such changes uptill the current strain step, one can define the cumulative bond status changes
A sudden change in the slope of cummulative bond order change with increasing strain correlates directly with the yield point of the cluster. This functionality is included for Vitreloy-105 in the function,

### yield_strain_multi_bs(max_alpha,max_beta,total_disp_steps,no_atoms,iteration,no_segs)
* **no_atoms** - defines the STZ cluster size
* **max_alpha, max_beta** - defines the upper bound of the tunable parameters of the deformation gradient
* **total_disp_steps** - defines the density of mesh in the displacement space of α and β ( higher means finer displacement grid on which stresses will be computed)
* **iteration** - defines the STZ cluster number of a given size
* **no_segs** - to track the sudden change in slop, the cummulative bond frequency plot is fitted with a piece wise linear function, the number of segments is defined by this paramter
* 
This function will give yield stress and strain at zero normal stress as output, that would be used by post-processing codes to compute the Mohr-Coloumb parameter

**NOTE: Some support functions that are used in this section of code, will require functions from the misc.py in the utilities section**








