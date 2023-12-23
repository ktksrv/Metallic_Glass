def MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc,shift):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import signal, interpolate, stats
    from scipy.optimize import curve_fit
    import STZ
    no_of_task = int(total_proc/proc_per_task)
    dat_trans = np.loadtxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration))
    tau_o_zero = dat_trans[1]
    bs_change = dat_trans[0]
    upper_l = bs_change+0.01
    upper_l_neg = bs_change+shift        # value added depends upon difference in strain in potential yield points (smaller sizes less value -> not much shift in yield point with normal stress )
    upper_l_pos = bs_change+shift

    factors_ns=np.linspace(lower,upper,steps)
    factors_ns = np.delete(factors_ns,int((steps-1)/2))
    factor_ns_batch = np.split(factors_ns,no_of_task)
    factor = factor_ns_batch[current_proc]
    normal_s = factor*tau_o_zero

    string='Normal_stress_'
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    d_alpha = np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    # d_alpha = d_alpha[130:300]
    d_beta = np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    cs = plt.contour(Beta,Alpha,ns,normal_s)
    dat0=cs.allsegs[0][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0)) 
    for j in range(len(dat0)):
        _,ss_lvl[j]=STZ.stress_config_grad(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,iteration)
        s_strain[j]=dat0[j][0]
    string=string+str(round(normal_s[0]/tau_o_zero,4)) 
    ss_data = np.column_stack((s_strain, ss_lvl))

    if(normal_s>=0):
        upper_l_corr = upper_l_pos
    else:
        upper_l_corr = upper_l_neg

    stress_diff = np.zeros(len(ss_lvl))
    for p in range(len(ss_lvl)-1):
        stress_diff[p] = ss_lvl[p+1]-ss_lvl[p]
    drop_strains = s_strain[stress_diff<0]
    drop_strains_useful = drop_strains[drop_strains<upper_l_corr]
    if(len(drop_strains_useful)==0):
        stress_useful = ss_lvl[s_strain<upper_l_corr]
        max_shear_lvl = stress_useful.max()
        strain_arg = stress_useful.argmax()
        max_yeild_strain = s_strain[strain_arg]
    else:
        y_interp = interpolate.interp1d(s_strain,ss_lvl)
        max_shear_lvl = y_interp(drop_strains_useful).max()
        index = y_interp(drop_strains_useful).argmax()
        max_yeild_strain = drop_strains[index]

    plt.figure(figsize=[6,4], dpi=300)
    plt.title('{} atoms STZ at normalised normal stress = {}'.format(no_atoms,round(factor[0],4)))
    plt.plot(s_strain,ss_lvl/10**9)
    plt.plot(max_yeild_strain,max_shear_lvl/10**9, marker='x', markersize=10, color='r' )
    plt.xlabel('Shear strain')
    plt.ylabel('Shear Stress (GPa)')
    np.savetxt('Clusters/{}_cluster/cluster_{}/{}.txt'.format(no_atoms,iteration,string),ss_data)
    plt.savefig('Clusters/{}_cluster/cluster_{}/{}.png'.format(no_atoms,iteration,string), facecolor='white', transparent=False)
    plt.clf()
    string='Normal_stress_'

    MC_data=np.column_stack((factor[0],max_shear_lvl/tau_o_zero))
    np.savetxt('Clusters/{}_cluster/cluster_{}/MC_data_{}'.format(no_atoms,iteration,iteration),MC_data)


# Vitreloy-105
def stress_config_grad(alpha,beta,initial,type_particle,no_atoms,iteration):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()     
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
    ns,ss=STZ.stress_calc(coordinates,force)
    return ns, ss

Pd-Ni-P

def stress_config_grad(alpha,beta,initial,type_particle,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
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
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "    

    surface_group='\ngroup surface id {}'.format(surface_str)
    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
    ns,ss=STZ.stress_calc(coordinates,force)
    return ns, ss

Cu-Cu
def stress_config_grad(alpha,beta,initial,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()  
    lmp=lammps()
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
    ns,ss=STZ.stress_calc(coordinates,force)

    return ns, ss