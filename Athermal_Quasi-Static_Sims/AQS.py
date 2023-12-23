
def stress_calc(coordinates,force):
    import numpy as np
    from scipy.spatial import ConvexHull
    sigma=np.zeros([3,3])
    for i in range(len(force)):
        F = force[i,:]
        F=F.reshape(3,1)
        R = coordinates[i,:]
        sigma = sigma + np.kron(F,R)
    hull= ConvexHull(coordinates) 
    volume=hull.volume
    sigma=-sigma/volume
    sigma= sigma*1.6021766208e+11   #eV/ang^3  to SI
    return sigma[2,2], sigma[0,2]


# Vitreloy - 105
def Normal_Stress_Matrix_p(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import numpy as np
    import STZ

    ## NORMAL STRESS MATRIX GENERATION ##
    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    reduced_alpha = d_alpha[130:300]  ## actual use alpha range
    # reduced_alpha = d_alpha ## actual use alpha range
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(reduced_alpha,no_of_task)
    ns=np.zeros([int(len(reduced_alpha)/no_of_task),len(d_beta)])
    initial,type_particle=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    y=-1
    for i in batch_alpha[current_proc]:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
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

            create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)


            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(no_atoms,iteration,current_proc),ns)


# Pd-Ni-P
def Normal_Stress_Matrix_p(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import numpy as np
    import STZ

    ## NORMAL STRESS MATRIX GENERATION ##
    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    reduced_alpha = d_alpha[130:300]  ## actual use aplha rnage
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(reduced_alpha,no_of_task)
    ns=np.zeros([int(len(reduced_alpha)/no_of_task),len(d_beta)])
    initial,type_particle=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    y=-1
    for i in batch_alpha[current_proc]:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
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
            create_box  4 myregion
            mass 1 106.42 
            mass 2 63.546
            mass 3 58.693
            mass 4 30.974
            pair_style nep nep.population160.generation1013200.txt
            pair_coeff * * 
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

            create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)


            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(no_atoms,iteration,current_proc),ns)

# Cu-Cu
def Normal_Stress_Matrix_p(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import matplotlib.pyplot as plt
    import numpy as np
    import STZ

    ## NORMAL STRESS MATRIX GENERATION ##
    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    reduced_alpha = d_alpha[125:265]  ## actual use alpha range
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(reduced_alpha,no_of_task)
    ns=np.zeros([int(len(reduced_alpha)/no_of_task),len(d_beta)])
    initial=STZ.extract_dat('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    y=-1
    for i in batch_alpha[current_proc]:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
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
            create_box  1 myregion
            mass 1 63.546

            pair_style eam
            pair_coeff 1 1 Cu_u6.eam
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

            create_atoms=['create_atoms 1 single {} {} {}'.format(d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)

            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(no_atoms,iteration,current_proc),ns)