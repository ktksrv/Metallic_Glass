import numpy as np 
no_task = 10
no_atoms = 1013
total_disp_steps = 169
batch_ns_size = int((total_disp_steps+1)/no_task)
# ns = np.zeros([int(total_disp_steps+1), int(total_disp_steps+1)])
ns = np.zeros([int(total_disp_steps+1), 350])
for iteration in [1,2,3,4,5]:
    for i in range(no_task):
        path = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(i,no_atoms,iteration,i)
        ns_batch = np.loadtxt(path)
        ns[0+batch_ns_size*i:batch_ns_size+batch_ns_size*i,:] = ns_batch
        np.savetxt('/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1//Latest_full_shear_trial/Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration),ns)











