# Post-Processing
This directory contains code for finding the local Mohr-Coloumb parameter for an STZ cluster

## Method
1. The iso-normal contour is extracted from the normal stress matrix generated from Athermal Quasi-Static Simulations for a given normal stress
2. For the strains on the iso-normal contours, shear stress is evaluated and plotted against the shear strain to estimate the yield shear strain using the sudden change in slope of cumulative bond status change parameter
3. Repeating this process for different normal stresses, we can plot normalized normal stress vs normalized yield shear stress, whose slope yields the local Mohr-Coloumb parameter for the cluster 

The procedure is implemented using,
##### MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc,shift)
* **no_atoms** - defines the STZ cluster size
* **max_alpha, max_beta** - defines the upper bound of the tunable parameters of the deformation gradient
* **total_disp_steps** - defines the density of mesh in the displacement space of α and β ( higher means finer displacement grid on which stresses will be computed)
* **lower, upper, steps** - defines the upper and lower limit of normal stresses and the number of steps in between these limits, where the sampling of specific normal stress will take place
* **iteration** - defines the STZ cluster number of a given size
* **total_proc,proc_per_task,current_pro**c - defines the parameters for leveraging parallel CPU architecture as all the calculations are embarrassingly parallel.
* **shift** - defines the gap between two sucessive  large drops in shear stress ( if such a situation arises in a stress-strain plot toggle this parameter to match with bond data)  


