# Athermal Quasi-Static Simulations
This directory contains the code to run Athermal Quasi-Static Simulations on the STZ cluster extracted from the Bulk Metallic Glass generated using LAMMPS

## Method
1. The required combination of shear and normal strain is imposed on the STZ cluster by applying the deformation gradient matrix onto its initial coordinates.
2. The α and β parameters of the deformation gradient are varied to produce a variety of combinations of stress states ( in the current workflow α - (-0.2,0.2) and β 
-(0,0.4) )
3. At each value of α and β, after relaxing the interior atoms, while keeping the surface atoms fixed, the stress tensor is computed using a virial formulation of the stress on an atomistic scale.

The following procedure is implemented using two main functions 
#### Normal_Stress_Matrix_p(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc)
no_atoms - defines the STZ cluster size
max_alpha, max_beta - defines the upper bound of the tunable parameters of the deformation gradient
total_disp_steps - defines the density of mesh in the displacement space of α and β ( higher means finer displacement grid on which stresses will be computed)
iteration - defines the STZ cluster number
total_proc,proc_per_task,current_proc - defines the parameters for leveraging parallel CPU architecture as all the calculations are embarrassingly parallel. 

The output will be a (total_disp_steps x total_disp_steps) matrix, whose each element represents the normal stress at the given stress state ( defined for that element at a given value of α and β )

#### surface_atoms(initial,no_atoms)
initial - define the initial coordinates of the relaxed STZ cluster
no_atoms - defines the STZ cluster size
This OVITO ported code section identifies the IDs of atoms at the surface of the clusters and gives it as the output to be used by the above function

 Markup : * Bullet list
              * Nested bullet
                  * Sub-nested bullet etc
          * Bullet list item 2
