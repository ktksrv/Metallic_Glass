def dist_NN(ID1,ID2,no_atoms):
    import STZ
    import math
    coordinates,_=STZ.extract_fdump('dump.force',no_atoms)
    atom_1_coord = coordinates[ID1-1,:]
    atom_2_coord = coordinates [ID2-1,:]
    interatomic_dist_12= math.dist(atom_1_coord,atom_2_coord)
    return interatomic_dist_12

def Normal_Stress_Matrix_bond_len_multi(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import numpy as np
    import STZ

    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(d_alpha,no_of_task)
    no_pass = len(batch_alpha[current_proc])*len(d_beta)
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu
    type_3_count = np.count_nonzero(particle_type==3) # Ni
    type_4_count = np.count_nonzero(particle_type==4) # Al
    type_5_count = np.count_nonzero(particle_type==5) # Ti


    nn_dist_mat_1_1=np.zeros([type_1_count, no_pass])
    nn_dist_mat_1_2=np.zeros([type_1_count, no_pass])
    nn_dist_mat_1_3=np.zeros([type_1_count, no_pass])
    nn_dist_mat_1_4=np.zeros([type_1_count, no_pass])
    nn_dist_mat_1_5=np.zeros([type_1_count, no_pass])
    nn_dist_mat_2_2=np.zeros([type_2_count, no_pass])
    nn_dist_mat_2_3=np.zeros([type_2_count, no_pass])
    nn_dist_mat_2_4=np.zeros([type_2_count, no_pass])
    nn_dist_mat_2_5=np.zeros([type_2_count, no_pass])
    nn_dist_mat_3_3=np.zeros([type_3_count, no_pass])
    nn_dist_mat_3_4=np.zeros([type_3_count, no_pass])
    nn_dist_mat_3_5=np.zeros([type_3_count, no_pass])
    nn_dist_mat_4_4=np.zeros([type_4_count, no_pass])
    nn_dist_mat_4_5=np.zeros([type_4_count, no_pass])
    nn_dist_mat_5_5=np.zeros([type_5_count, no_pass])

    n_NN_1_1,n_NN_1_2,n_NN_1_3,n_NN_1_4,n_NN_1_5,n_NN_2_2,n_NN_2_3,n_NN_2_4,n_NN_2_5,n_NN_3_3,n_NN_3_4,n_NN_3_5,n_NN_4_4,n_NN_4_5,n_NN_5_5 = STZ.initial_NN_list_multi(initial,particle_type)
    y=-1
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    for i in batch_alpha[current_proc]:
        alpha=i
        for j in range(0,total_disp_steps+1):
            y=y+1
            beta=(max_beta/total_disp_steps)*j
            deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
            d_initial=np.zeros([no_atoms,3])
            for k in range(no_atoms):
                d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
            x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
            initialization_block='''
            dimension 3
            units metal
            boundary s s s
            atom_style atomic
            timestep 0.001
            region myregion block {} {} {} {} {} {}  units box
            create_box 5 myregion
            mass 1 91.224
            mass 2 63.546
            mass 3 58.693
            mass 4 26.982
            mass 5 47.867 
            pair_style eam/alloy
            pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)
            
            create_atoms=['create_atoms {} single {} {} {}'.format(particle_type[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)

            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            dump_modify fcal sort id
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            #Zr-Zr
            for atom_ID in range(type_1_count):
                ID1 = n_NN_1_1[atom_ID,0]
                ID2 = n_NN_1_1[atom_ID,1]
                nn_dist_mat_1_1[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Zr-Cu
            for atom_ID in range(type_1_count):
                ID1 = n_NN_1_2[atom_ID,0]
                ID2 = n_NN_1_2[atom_ID,1]
                nn_dist_mat_1_2[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Zr-Ni
            for atom_ID in range(type_1_count):
                ID1 = n_NN_1_3[atom_ID,0]
                ID2 = n_NN_1_3[atom_ID,1]
                nn_dist_mat_1_3[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)


            #Zr-Al
            for atom_ID in range(type_1_count):
                ID1 = n_NN_1_4[atom_ID,0]
                ID2 = n_NN_1_4[atom_ID,1]
                nn_dist_mat_1_4[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)


            #Zr-Ti
            for atom_ID in range(type_1_count):
                ID1 = n_NN_1_5[atom_ID,0]
                ID2 = n_NN_1_5[atom_ID,1]
                nn_dist_mat_1_5[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Cu-Cu 
            for atom_ID in range(type_2_count):
                ID1 = n_NN_2_2[atom_ID,0]
                ID2 = n_NN_2_2[atom_ID,1]
                nn_dist_mat_2_2[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Cu-Ni
            for atom_ID in range(type_2_count):
                ID1 = n_NN_2_3[atom_ID,0]
                ID2 = n_NN_2_3[atom_ID,1]
                nn_dist_mat_2_3[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)  

            #Cu-Al
            for atom_ID in range(type_2_count):
                ID1 = n_NN_2_4[atom_ID,0]
                ID2 = n_NN_2_4[atom_ID,1]
                nn_dist_mat_2_4[atom_ID,y] = dist_NN(ID1,ID2,no_atoms) 

            #Cu-Ti
            for atom_ID in range(type_2_count):
                ID1 = n_NN_2_5[atom_ID,0]
                ID2 = n_NN_2_5[atom_ID,1]
                nn_dist_mat_2_5[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)   

            #Ni-Ni
            for atom_ID in range(type_3_count):
                ID1 = n_NN_3_3[atom_ID,0]
                ID2 = n_NN_3_3[atom_ID,1]
                nn_dist_mat_3_3[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Ni-Al
            for atom_ID in range(type_3_count):
                ID1 = n_NN_3_4[atom_ID,0]
                ID2 = n_NN_3_4[atom_ID,1]
                nn_dist_mat_3_4[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)
    
            #Ni-Ti
            for atom_ID in range(type_3_count):
                ID1 = n_NN_3_5[atom_ID,0]
                ID2 = n_NN_3_5[atom_ID,1]
                nn_dist_mat_3_5[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)
    
            #Al-Al
            for atom_ID in range(type_4_count):
                ID1 = n_NN_4_4[atom_ID,0]
                ID2 = n_NN_4_4[atom_ID,1]
                nn_dist_mat_4_4[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)
    
            #Al-Ti
            for atom_ID in range(type_4_count):
                ID1 = n_NN_4_5[atom_ID,0]
                ID2 = n_NN_4_5[atom_ID,1]
                nn_dist_mat_4_5[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)

            #Ti-Ti
            for atom_ID in range(type_5_count):
                ID1 = n_NN_5_5[atom_ID,0]
                ID2 = n_NN_5_5[atom_ID,1]
                nn_dist_mat_5_5[atom_ID,y] = dist_NN(ID1,ID2,no_atoms)
    
        
        np.savetxt('bond_len_1-1.txt',nn_dist_mat_1_1)
        np.savetxt('bond_len_1-2.txt',nn_dist_mat_1_2)
        np.savetxt('bond_len_1-3.txt',nn_dist_mat_1_3)
        np.savetxt('bond_len_1-4.txt',nn_dist_mat_1_4)
        np.savetxt('bond_len_1-5.txt',nn_dist_mat_1_5)
        np.savetxt('bond_len_2-2.txt',nn_dist_mat_2_2)
        np.savetxt('bond_len_2-3.txt',nn_dist_mat_2_3)
        np.savetxt('bond_len_2-4.txt',nn_dist_mat_2_4)
        np.savetxt('bond_len_2-5.txt',nn_dist_mat_2_5)
        np.savetxt('bond_len_3-3.txt',nn_dist_mat_3_3)
        np.savetxt('bond_len_3-4.txt',nn_dist_mat_3_4)
        np.savetxt('bond_len_3-5.txt',nn_dist_mat_3_5)
        np.savetxt('bond_len_4-4.txt',nn_dist_mat_4_4)
        np.savetxt('bond_len_4-5.txt',nn_dist_mat_4_5)
        np.savetxt('bond_len_5-5.txt',nn_dist_mat_5_5)

def initial_NN_list_multi(initial,particle_type):
    import STZ
    import numpy as np
    no_atoms = initial.shape[0]
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu 
    type_3_count = np.count_nonzero(particle_type==3) # Ni 
    type_4_count = np.count_nonzero(particle_type==4) # Al   
    type_5_count = np.count_nonzero(particle_type==5) # Ti   
    type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,sorted_initial_ID = STZ.init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms)

    #Zr-Zr
    n_NN_1_1 = np.zeros([type_1_count,2])
    n_NN_1_1[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_1[i,:].argsort()
        n_NN_1_1[i,1] = sorted_initial_ID[per_atom_dist_sort_arg[1]]
    n_NN_1_1 = n_NN_1_1.astype(int)

    #Zr-Cu
    n_NN_1_2 = np.zeros([type_1_count,2])
    n_NN_1_2[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_2[i,:].argsort()
        n_NN_1_2[i,1] = sorted_initial_ID[type_1_count+per_atom_dist_sort_arg[1]]
    n_NN_1_2 = n_NN_1_2.astype(int)

    #Zr-Ni
    n_NN_1_3 = np.zeros([type_1_count,2])
    n_NN_1_3[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_3[i,:].argsort()
        n_NN_1_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_1_3 = n_NN_1_3.astype(int)

    #Zr-Al
    n_NN_1_4 = np.zeros([type_1_count,2])
    n_NN_1_4[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_4[i,:].argsort()
        n_NN_1_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_1_4 = n_NN_1_4.astype(int)

    #Zr-Ti
    n_NN_1_5 = np.zeros([type_1_count,2])
    n_NN_1_5[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_5[i,:].argsort()
        n_NN_1_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_1_5 = n_NN_1_5.astype(int)

    #Cu-Cu
    n_NN_2_2 = np.zeros([type_2_count,2])
    n_NN_2_2[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_2[i,:].argsort()
        n_NN_2_2[i,1] = sorted_initial_ID[type_1_count+per_atom_dist_sort_arg[1]]
    n_NN_2_2 = n_NN_2_2.astype(int)

    #Cu-Ni
    n_NN_2_3 = np.zeros([type_2_count,2])
    n_NN_2_3[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_3[i,:].argsort()
        n_NN_2_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_2_3 = n_NN_2_3.astype(int)

    #Cu-Al
    n_NN_2_4 = np.zeros([type_2_count,2])
    n_NN_2_4[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_4[i,:].argsort()
        n_NN_2_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_2_4 = n_NN_2_4.astype(int)

    #Cu-Ti
    n_NN_2_5 = np.zeros([type_2_count,2])
    n_NN_2_5[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_5[i,:].argsort()
        n_NN_2_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_2_5 = n_NN_2_5.astype(int)
       
    #Ni-Ni
    n_NN_3_3 = np.zeros([type_3_count,2])
    n_NN_3_3[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_3[i,:].argsort()
        n_NN_3_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_3_3 = n_NN_3_3.astype(int)

    #Ni-Al
    n_NN_3_4 = np.zeros([type_3_count,2])
    n_NN_3_4[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_4[i,:].argsort()
        n_NN_3_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_3_4 = n_NN_3_4.astype(int)

    #Ni-Ti
    n_NN_3_5 = np.zeros([type_3_count,2])
    n_NN_3_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_5[i,:].argsort()
        n_NN_3_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_3_5 = n_NN_3_5.astype(int)

    #Al-Al
    n_NN_4_4 = np.zeros([type_4_count,2])
    n_NN_4_4[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count:type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_4_count):
        per_atom_dist_sort_arg = type_4_4[i,:].argsort()
        n_NN_4_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_4_4 = n_NN_4_4.astype(int)

    #Al-Ti
    n_NN_4_5 = np.zeros([type_4_count,2])
    n_NN_4_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count:type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_4_count):
        per_atom_dist_sort_arg = type_4_5[i,:].argsort()
        n_NN_4_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_4_5 = n_NN_4_5.astype(int)


    #Ti-Ti
    n_NN_5_5 = np.zeros([type_5_count,2])
    n_NN_5_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count:type_5_count+type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_5_count):
        per_atom_dist_sort_arg = type_5_5[i,:].argsort()
        n_NN_5_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_5_5 = n_NN_5_5.astype(int)


    return n_NN_1_1,n_NN_1_2,n_NN_1_3,n_NN_1_4,n_NN_1_5,n_NN_2_2,n_NN_2_3,n_NN_2_4,n_NN_2_5,n_NN_3_3,n_NN_3_4,n_NN_3_5,n_NN_4_4,n_NN_4_5,n_NN_5_5

def init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms):
    import STZ
    import numpy as np
    import math
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu 
    type_3_count = np.count_nonzero(particle_type==3) # Ni 
    type_4_count = np.count_nonzero(particle_type==4) # Al   
    type_5_count = np.count_nonzero(particle_type==5) # Ti   
    type_1_1 = np.zeros([type_1_count,type_1_count])
    type_1_2 = np.zeros([type_1_count,type_2_count])
    type_1_3 = np.zeros([type_1_count,type_3_count])
    type_1_4 = np.zeros([type_1_count,type_4_count])
    type_1_5 = np.zeros([type_1_count,type_5_count])
    type_2_2 = np.zeros([type_2_count,type_2_count])
    type_2_3 = np.zeros([type_2_count,type_3_count])
    type_2_4 = np.zeros([type_2_count,type_4_count])
    type_2_5 = np.zeros([type_2_count,type_5_count])
    type_3_3 = np.zeros([type_3_count,type_3_count])
    type_3_4 = np.zeros([type_3_count,type_4_count])
    type_3_5 = np.zeros([type_3_count,type_5_count])
    type_4_4 = np.zeros([type_4_count,type_4_count])
    type_4_5 = np.zeros([type_4_count,type_5_count])
    type_5_5 = np.zeros([type_5_count,type_5_count])

    sort = np.argsort(particle_type)
    type_sorted_initial = initial[sort]
    type_sorted_initial_coord = type_sorted_initial[:,1:]
    for i in range(type_1_count):
        for j in range(type_1_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j,:]
            type_1_1[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_1_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_1_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_1_3[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_1_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_1_5[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_2_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_2_3[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_2_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_2_4[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_2_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_2_5[i,j]= math.dist(atom_i_coord,atom_j_coord)    
    
    for i in range(type_3_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_3_3[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_3_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_3_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_3_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_3_5[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_4_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_4_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_4_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_4_5[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_5_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count+type_4_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_5_5[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    return type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,sort

# Yield Estimation
def yield_strain_multi_bs(max_alpha,max_beta,total_disp_steps,no_atoms,iteration,no_segs):
    import STZ
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_alpha = d_alpha[130:300]  # Use range
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_1_4,initial_interatomic_mat_1_5,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_2_4,initial_interatomic_mat_2_5,initial_interatomic_mat_3_3,initial_interatomic_mat_3_4,initial_interatomic_mat_3_5,initial_interatomic_mat_4_4,initial_interatomic_mat_4_5,initial_interatomic_mat_5_5, sort_order = STZ.init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms)
    cs = plt.contour(Beta,Alpha,ns,0)
    dat0=cs.allsegs[1][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0))
    initial_bonded_matrix_1_1,initial_bonded_matrix_1_2,initial_bonded_matrix_1_3,initial_bonded_matrix_1_4,initial_bonded_matrix_1_5,initial_bonded_matrix_2_2,initial_bonded_matrix_2_3,initial_bonded_matrix_2_4,initial_bonded_matrix_2_5,initial_bonded_matrix_3_3,initial_bonded_matrix_3_4,initial_bonded_matrix_3_5,initial_bonded_matrix_4_4,initial_bonded_matrix_4_5,initial_bonded_matrix_5_5 = STZ.bond_mat_from_interatomic_dist_mat_multi(initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_1_4,initial_interatomic_mat_1_5,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_2_4,initial_interatomic_mat_2_5,initial_interatomic_mat_3_3,initial_interatomic_mat_3_4,initial_interatomic_mat_3_5,initial_interatomic_mat_4_4,initial_interatomic_mat_4_5,initial_interatomic_mat_5_5,3.35036517,3.55296326,3.29013609,3.20843578,3.38126814,3.21116134,3.09788097,3.74233114,3.50770858,3.29733712,3.84690821,3.60998013,3.64419147,3.70396619,5)
    bond_change_frequency = np.zeros(len(dat0))
    dat_transfer = np.zeros(2)
    for j in range(len(dat0)):
        current_bonded_matrix_1_1,current_bonded_matrix_1_2,current_bonded_matrix_1_3,current_bonded_matrix_1_4,current_bonded_matrix_1_5,current_bonded_matrix_2_2,current_bonded_matrix_2_3,current_bonded_matrix_2_4,current_bonded_matrix_2_5,current_bonded_matrix_3_3,current_bonded_matrix_3_4,current_bonded_matrix_3_5,current_bonded_matrix_4_4,current_bonded_matrix_4_5,current_bonded_matrix_5_5,_,ss_lvl[j] = STZ.bond_matrix_stress_config_multi(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,sort_order,iteration)
        s_strain[j]=dat0[j][0]
        bond_status_change_matrix_1_1 = current_bonded_matrix_1_1-initial_bonded_matrix_1_1
        bond_status_change_matrix_1_2 = current_bonded_matrix_1_2-initial_bonded_matrix_1_2
        bond_status_change_matrix_1_3 = current_bonded_matrix_1_3-initial_bonded_matrix_1_3
        bond_status_change_matrix_1_4 = current_bonded_matrix_1_4-initial_bonded_matrix_1_4
        bond_status_change_matrix_1_5 = current_bonded_matrix_1_5-initial_bonded_matrix_1_5
        bond_status_change_matrix_2_2 = current_bonded_matrix_2_2-initial_bonded_matrix_2_2
        bond_status_change_matrix_2_3 = current_bonded_matrix_2_3-initial_bonded_matrix_2_3
        bond_status_change_matrix_2_4 = current_bonded_matrix_2_4-initial_bonded_matrix_2_4
        bond_status_change_matrix_2_5 = current_bonded_matrix_2_5-initial_bonded_matrix_2_5
        bond_status_change_matrix_3_3 = current_bonded_matrix_3_3-initial_bonded_matrix_3_3
        bond_status_change_matrix_3_4 = current_bonded_matrix_3_4-initial_bonded_matrix_3_4
        bond_status_change_matrix_3_5 = current_bonded_matrix_3_5-initial_bonded_matrix_3_5
        bond_status_change_matrix_4_4 = current_bonded_matrix_4_4-initial_bonded_matrix_4_4
        bond_status_change_matrix_4_5 = current_bonded_matrix_4_5-initial_bonded_matrix_4_5
        # bond_status_change_matrix_5_5 = current_bonded_matrix_5_5-initial_bonded_matrix_5_5

        initial_bonded_matrix_1_1 = current_bonded_matrix_1_1
        initial_bonded_matrix_1_2 = current_bonded_matrix_1_2
        initial_bonded_matrix_1_3 = current_bonded_matrix_1_3
        initial_bonded_matrix_1_4 = current_bonded_matrix_1_4
        initial_bonded_matrix_1_5 = current_bonded_matrix_1_5
        initial_bonded_matrix_2_2 = current_bonded_matrix_2_2
        initial_bonded_matrix_2_3 = current_bonded_matrix_2_3
        initial_bonded_matrix_2_4 = current_bonded_matrix_2_4
        initial_bonded_matrix_2_5 = current_bonded_matrix_2_5
        initial_bonded_matrix_3_3 = current_bonded_matrix_3_3
        initial_bonded_matrix_3_4 = current_bonded_matrix_3_4
        initial_bonded_matrix_3_5 = current_bonded_matrix_3_5
        initial_bonded_matrix_4_4 = current_bonded_matrix_4_4
        initial_bonded_matrix_4_5 = current_bonded_matrix_4_5
        # initial_bonded_matrix_3_3 = current_bonded_matrix_5_5
        bond_change_frequency[j] = np.count_nonzero(bond_status_change_matrix_1_1)+np.count_nonzero(bond_status_change_matrix_1_2)+np.count_nonzero(bond_status_change_matrix_1_3)+np.count_nonzero(bond_status_change_matrix_1_4)+np.count_nonzero(bond_status_change_matrix_1_5)+np.count_nonzero(bond_status_change_matrix_2_2)+np.count_nonzero(bond_status_change_matrix_2_3)+np.count_nonzero(bond_status_change_matrix_2_4)+np.count_nonzero(bond_status_change_matrix_2_5)+np.count_nonzero(bond_status_change_matrix_3_3)+np.count_nonzero(bond_status_change_matrix_3_4)+np.count_nonzero(bond_status_change_matrix_3_5)+np.count_nonzero(bond_status_change_matrix_4_4)+np.count_nonzero(bond_status_change_matrix_4_5)

    bond_change_frequency_trunc = np.delete(bond_change_frequency,0)
    cummulative_sum = 0
    bond_change_cumulative = np.zeros(len(bond_change_frequency_trunc))
    for i in range(len(bond_change_frequency_trunc)):
        cummulative_sum = cummulative_sum+bond_change_frequency[i+1]
        bond_change_cumulative[i] = cummulative_sum

    bond_change_cumulative_averaged= STZ.moving_average(bond_change_cumulative,4)
    bond_status_data =np.column_stack((dat0[0:-4,0],bond_change_cumulative_averaged))                                                               
    X, Y = bond_status_data[0:360,0], bond_status_data[0:360,1]
    px, py = STZ.segments_fit(X, Y, no_segs)
    thresh_zero = px[1]+0.025
    stress_diff = np.zeros(len(ss_lvl))
    for i in range(len(ss_lvl)-1):
        stress_diff[i] = ss_lvl[i+1]-ss_lvl[i]
    drop_strains = s_strain[stress_diff<0]
    drop_strains_useful = drop_strains[drop_strains< thresh_zero]
    if(len(drop_strains_useful)!=0):
        y_interp = interpolate.interp1d(s_strain,ss_lvl)
        tau_o_zero = y_interp(drop_strains_useful).max()
        strain_index = y_interp(drop_strains_useful).argmax()
        yield_strain = drop_strains_useful[strain_index]
    else:
        strains_useful = s_strain[s_strain<thresh_zero]
        y_interp = interpolate.interp1d(s_strain,ss_lvl)
        tau_o_zero = y_interp(strains_useful).max()
        strain_index = y_interp(strains_useful).argmax()
        yield_strain = strains_useful[strain_index]

    dat_transfer[0] = yield_strain
    dat_transfer[1] = tau_o_zero

    plt.figure(figsize=[6,4], dpi=300)
    plt.title('{} atoms STZ at normalised normal stress = 0.0'.format(no_atoms))
    ax1 = plt.subplot()
    color = 'tab:red'
    ax1.set_xlabel('Shear Strain')
    ax1.set_ylabel('Shear Stress (GPa)', color=color) 
    l1, = ax1.plot(s_strain,ss_lvl/10**9, color=color)
    ax1.plot(yield_strain, tau_o_zero/10**9, marker='x', markersize=10, color='r' )
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('Cummulative Bond Status Changes', color=color)
    l2, = ax2.plot(bond_status_data[0:360,0], bond_status_data[0:360,1], ".",color=color,)
    ax2.plot(px, py, "-or")
    ax2.tick_params(axis='y', labelcolor=color)
    plt.legend([l2, l1], ["Bond status changes", "Stress vs Strain"])
    plt.savefig('Clusters/{}_cluster/cluster_{}/bond_freq_yield_0.0.png'.format(no_atoms,iteration), facecolor='white', transparent=False)
    np.savetxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration), dat_transfer)
    plt.clf()


def bond_mat_from_interatomic_dist_mat_multi(type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,thresh_1_1,thresh_1_2,thresh_1_3,thresh_1_4,thresh_1_5,thresh_2_2,thresh_2_3,thresh_2_4,thresh_2_5,thresh_3_3,thresh_3_4,thresh_3_5,thresh_4_4,thresh_4_5,thresh_5_5):
    import numpy as np
    bond_mat_1_1 = np.where(type_1_1<thresh_1_1,1,0)
    bond_mat_1_2 = np.where(type_1_2<thresh_1_2,1,0)
    bond_mat_1_3 = np.where(type_1_3<thresh_1_3,1,0)
    bond_mat_1_4 = np.where(type_1_4<thresh_1_4,1,0)
    bond_mat_1_5 = np.where(type_1_5<thresh_1_5,1,0)
    bond_mat_2_2 = np.where(type_2_2<thresh_2_2,1,0)
    bond_mat_2_3 = np.where(type_2_3<thresh_2_3,1,0)
    bond_mat_2_4 = np.where(type_2_4<thresh_2_4,1,0)
    bond_mat_2_5 = np.where(type_2_5<thresh_2_5,1,0)
    bond_mat_3_3 = np.where(type_3_3<thresh_3_3,1,0)
    bond_mat_3_4 = np.where(type_3_4<thresh_3_4,1,0)
    bond_mat_3_5 = np.where(type_3_5<thresh_3_5,1,0)
    bond_mat_4_4 = np.where(type_4_4<thresh_4_4,1,0)
    bond_mat_4_5 = np.where(type_4_5<thresh_4_5,1,0)
    bond_mat_5_5 = np.where(type_5_5<thresh_5_5,1,0)

    return bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5

def bond_matrix_stress_config_multi(alpha,beta,initial,type_particle,no_atoms,sorted_ID,iteration):
    import numpy as np
    import STZ
    from lammps import lammps
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block {} {} {} {} {} {}  units box
    create_box 5 myregion
    mass 1 91.224
    mass 2 63.546
    mass 3 58.693
    mass 4 26.982
    mass 5 47.867 
    pair_style eam/alloy
    pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
    '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

    create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    create_atoms_str = '\n'.join(str(e) for e in create_atoms)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    compute force all property/atom fx fy fz
    dump fcal all custom 1 dump.force id type x y z fx fy fz
    dump_modify fcal sort id
    run 1 
    '''

    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,force,particle_type=STZ.extract_fdump_type('dump.force',no_atoms)
    coordinates_sorted = coordinates[sorted_ID]
    type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,_ = STZ.init_interatomic_dist_matrix_multi(coordinates_sorted,particle_type,no_atoms)
    bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5 = STZ.bond_mat_from_interatomic_dist_mat_multi(type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,3.35036517,3.55296326,3.29013609,3.20843578,3.38126814,3.21116134,3.09788097,3.74233114,3.50770858,3.29733712,3.84690821,3.60998013,3.64419147,3.70396619,5)
    ns,ss=STZ.stress_calc(coordinates,force)
    return bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5, ns ,ss
