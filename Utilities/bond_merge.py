import numpy as np 
no_task = 10
no_atoms = 1013
bl_1_1_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_1-1.txt")
bl_1_2_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_1-2.txt")
bl_1_3_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_1-3.txt")
bl_1_4_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_1-4.txt")
bl_1_5_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_1-5.txt")
bl_2_2_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_2-2.txt")
bl_2_3_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_2-3.txt")
bl_2_4_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_2-4.txt")
bl_2_5_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_2-5.txt")
bl_3_3_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_3-3.txt")
bl_3_4_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_3-4.txt")
bl_3_5_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_3-5.txt")
bl_4_4_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_4-4.txt")
bl_4_5_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_4-5.txt")
bl_5_5_intial = np.loadtxt("/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_0/bond_len_5-5.txt")


for i in range(1,10):
    path_1_1 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_1-1.txt'.format(i)
    path_1_2 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_1-2.txt'.format(i)
    path_1_3 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_1-3.txt'.format(i)
    path_1_4 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_1-4.txt'.format(i)
    path_1_5 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_1-5.txt'.format(i)
    path_2_2 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_2-2.txt'.format(i)
    path_2_3 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_2-3.txt'.format(i)
    path_2_4 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_2-4.txt'.format(i)
    path_2_5 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_2-5.txt'.format(i)
    path_3_3 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_3-3.txt'.format(i)
    path_3_4 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_3-4.txt'.format(i)
    path_3_5 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_3-5.txt'.format(i)
    path_4_4 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_4-4.txt'.format(i)
    path_4_5 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_4-5.txt'.format(i)
    path_5_5 = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/bond_len_5-5.txt'.format(i)

    bl_batch_1_1 = np.loadtxt(path_1_1)
    bl_batch_1_2 = np.loadtxt(path_1_2)
    bl_batch_1_3 = np.loadtxt(path_1_3)
    bl_batch_1_4 = np.loadtxt(path_1_4)
    bl_batch_1_5 = np.loadtxt(path_1_5)
    bl_batch_2_2 = np.loadtxt(path_2_2)
    bl_batch_2_3 = np.loadtxt(path_2_3)
    bl_batch_2_4 = np.loadtxt(path_2_4)
    bl_batch_2_5 = np.loadtxt(path_2_5)
    bl_batch_3_3 = np.loadtxt(path_3_3)
    bl_batch_3_4 = np.loadtxt(path_3_4)
    bl_batch_3_5 = np.loadtxt(path_3_5)
    bl_batch_4_4 = np.loadtxt(path_4_4)
    bl_batch_4_5 = np.loadtxt(path_4_5)
    bl_batch_5_5 = np.loadtxt(path_5_5)

    bl_1_1_intial  = np.column_stack((bl_1_1_intial,bl_batch_1_1))
    bl_1_2_intial  = np.column_stack((bl_1_2_intial,bl_batch_1_2))
    bl_1_3_intial  = np.column_stack((bl_1_3_intial,bl_batch_1_3))
    bl_1_4_intial  = np.column_stack((bl_1_4_intial,bl_batch_1_4))
    bl_1_5_intial  = np.column_stack((bl_1_5_intial,bl_batch_1_5))
    bl_2_2_intial  = np.column_stack((bl_2_2_intial,bl_batch_2_2))
    bl_2_3_intial  = np.column_stack((bl_2_3_intial,bl_batch_2_3))
    bl_2_4_intial  = np.column_stack((bl_2_4_intial,bl_batch_2_4))
    bl_2_5_intial  = np.column_stack((bl_2_5_intial,bl_batch_2_5))
    bl_3_3_intial  = np.column_stack((bl_3_3_intial,bl_batch_3_3))
    bl_3_4_intial  = np.column_stack((bl_3_4_intial,bl_batch_3_4))
    bl_3_5_intial  = np.column_stack((bl_3_5_intial,bl_batch_3_5))
    bl_4_4_intial  = np.column_stack((bl_4_4_intial,bl_batch_4_4))
    bl_4_5_intial  = np.column_stack((bl_4_5_intial,bl_batch_4_5))
    bl_5_5_intial  = np.column_stack((bl_5_5_intial,bl_batch_5_5))


np.savetxt("bl_merged_1_1.txt",bl_1_1_intial)
np.savetxt("bl_merged_1_2.txt",bl_1_2_intial)
np.savetxt("bl_merged_1_3.txt",bl_1_3_intial)
np.savetxt("bl_merged_1_4.txt",bl_1_4_intial)
np.savetxt("bl_merged_1_5.txt",bl_1_5_intial)
np.savetxt("bl_merged_2_2.txt",bl_2_2_intial)
np.savetxt("bl_merged_2_3.txt",bl_2_3_intial)
np.savetxt("bl_merged_2_4.txt",bl_2_4_intial)
np.savetxt("bl_merged_2_5.txt",bl_2_5_intial)
np.savetxt("bl_merged_3_3.txt",bl_3_3_intial)
np.savetxt("bl_merged_3_4.txt",bl_3_4_intial)
np.savetxt("bl_merged_3_5.txt",bl_3_5_intial)
np.savetxt("bl_merged_4_4.txt",bl_4_4_intial)
np.savetxt("bl_merged_4_5.txt",bl_4_5_intial)
np.savetxt("bl_merged_5_5.txt",bl_5_5_intial)