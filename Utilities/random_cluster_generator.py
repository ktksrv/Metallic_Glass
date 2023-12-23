def random_cluster_centroid(x_ul,x_ll,y_ul,y_ll,z_ul,z_ll):
    import numpy as np
    x= x_ll+(x_ul-x_ll)*np.random.random()
    y= y_ll+(y_ul-x_ll)*np.random.random()
    z= z_ll+(z_ul-z_ll)*np.random.random()
    return x, y, z

def random_cluster_generator(size,iteration):
    from ovito.data import Particles, DataCollection, SimulationCell,ParticleType
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier
    from lammps import lammps
    from ovito.io import export_file
    import numpy as np
    import STZ

    # Random Cluster of given size
    current_count_diff=10 #flag
    #initial BMG size and parameters which define the usable region for cluster extraction
    no_atoms=32800                                  
    #radius for midpoint distance of the ranges below
    # r1,r2,r3 = 80,170,240   #Cu-Cu    
    # r1,r2,r3 = 85,180,250   #Pd-Ni-P 
    r1,r2,r3 = 105,220,310   #vit-105

    if size>0 and size<500 :
        r_initial=r1
    elif size>500 and size<1000 :
        r_initial=r2
    elif size>1000 and size<1500 :
        r_initial=r3
    r_updated=r_initial
    # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(54,20,54,20,90,54) #Cu-Cu
    # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(45,30,45,30,100,40) #Pd-Ni-P
    x_rand,y_rand,z_rand= STZ.random_cluster_centroid(69,21,69,21,106,74) #vit-105
    initial,type_particle=STZ.extract_dat_type('Initial_BMG.dat',no_atoms)
    iterations=0
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    #use cell vectors for given congifuration from ovito
    #Cu-Cu
    # cell[...] = [[72.56,0,0,0],                               
    #             [0,72.56,0,0],
    #             [0,0,145.12,0]]
    # #Pd-Ni-P
    # cell[...] = [[60.2781,0,0,6.14095],                               
    #             [0,60.2781,0,6.14095],
    #             [0,0,120.556,12.2819]]
    #vit
    cell[...] = [[90.51,0,0,0],                               
                [0,90.51,0,0],
                [0,0,181.02,0]]
    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    type_prop = particles.create_property('Particle Type')
    type_prop.types.append(ParticleType(id = 1, name = 'Zr', color = (0.0,1.0,0.0)))
    type_prop.types.append(ParticleType(id = 2, name = 'Cu', color = (1.0,0.0,0.0)))
    type_prop.types.append(ParticleType(id = 3, name = 'Ni', color = (0.0,0.0,1.0)))
    type_prop.types.append(ParticleType(id = 4, name = 'Al', color = (0.0,0.5,0.5)))
    type_prop.types.append(ParticleType(id = 5, name = 'Ti', color = (0.5,0.5,0.5)))

    for ind in range(no_atoms):
        type_prop[ind] = type_particle[ind]
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    while (current_count_diff!=0):
        if iterations==1000:
            # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(54,20,54,20,90,54) #Cu-Cu
            # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(45,30,45,30,100,40) #Pd-Ni-P
            x_rand,y_rand,z_rand= STZ.random_cluster_centroid(69,21,69,21,106,74) #vit-105
            iterations=0
        pipeline.modifiers.append(ExpressionSelectionModifier(expression = '(Position.X-{})^2+(Position.Y-{})^2+(Position.Z-{})^2 <{}'.format(x_rand,y_rand,z_rand,r_updated)))
        data= pipeline.compute()
        current_count=data.attributes['ExpressionSelection.count']
        current_count_diff=size-current_count
        del pipeline.modifiers[0]
        iterations+=1
        r_updated= r_updated+0.01*current_count_diff

    pipeline.modifiers.append(ExpressionSelectionModifier(expression = '(Position.X-{})^2+(Position.Y-{})^2+(Position.Z-{})^2 > {}'.format(x_rand,y_rand,z_rand,r_updated)))
    data = pipeline.compute()
    count=data.attributes['ExpressionSelection.count']
    pipeline.modifiers.append(DeleteSelectedModifier())
    export_file(pipeline, "Clusters/{}_cluster/cluster_{}/cluster_{}.dat".format(size,iteration,size), "lammps/data")
    lmp=lammps()
    # Minimize the cluster

    # Cu-Cu
    # minimization_block='''
    # dimension 3
    # units metal
    # boundary s s s
    # atom_style atomic
    # timestep 0.001

    # read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    # mass 1 63.546
    # pair_style eam
    # pair_coeff 1 1 Cu_u6.eam
    # minimize 0 1e-4 100000 100000


    # dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z
    # run 1

    # '''.format(size,iteration,size,size,iteration,size)
    # lmp.commands_string(minimization_block)


    # #Pd-Ni-P
    # minimization_block='''
    # dimension 3
    # units metal
    # boundary s s s
    # atom_style atomic
    # timestep 0.001

    # read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    # mass 1 106.42 
    # mass 2 63.546
    # mass 3 58.693
    # mass 4 30.974
    # pair_style nep nep.population160.generation1013200.txt
    # pair_coeff * * 

    # minimize 0 1e-4 100000 100000


    # dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z

    # run 1

    # '''.format(size,iteration,size,size,iteration,size)
    # lmp.commands_string(minimization_block)


    #vit-105
    minimization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001

    read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    mass 1 91.224
    mass 2 63.546
    mass 3 58.693
    mass 4 26.982
    mass 5 47.867 
    pair_style eam/alloy
    pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti 

    minimize 0 1e-4 100000 100000


    dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z

    run 1

    '''.format(size,iteration,size,size,iteration,size)
    lmp.commands_string(minimization_block)


    # Affine Transformations
    initial,type_particle=STZ.extract_dump_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(size,iteration,size),size)
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    #Transformed Simulation box, which should contain all the sheared atoms
    # Cu-Cu
    # cell[...] = [[200,0,0,0],                               
    #             [0,72.56,0,0],
    #             [0,0,145.12,0]]
    #Pd-Ni-P
    # cell[...] = [[180,0,0,6.14095],                               
    #             [0,60.2781,0,6.14095],
    #             [0,0,120.556,12.2819]]
    #vit-105
    cell[...] = [[250,0,0,0],                               
                [0,90.51,0,0],
                [0,0,181.02,0]]

    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    type_prop = particles.create_property('Particle Type')
    type_prop.types.append(ParticleType(id = 1, name = 'Zr', color = (0.0,1.0,0.0)))
    type_prop.types.append(ParticleType(id = 2, name = 'Cu', color = (1.0,0.0,0.0)))
    type_prop.types.append(ParticleType(id = 3, name = 'Ni', color = (0.0,0.0,1.0)))
    type_prop.types.append(ParticleType(id = 4, name = 'Al', color = (0.0,0.5,0.5)))
    type_prop.types.append(ParticleType(id = 4, name = 'Ti', color = (0.0,0.5,0.5)))

    for ind in range(size):
        type_prop[ind] = type_particle[ind]
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    export_file(pipeline, "Clusters/{}_cluster/cluster_{}/cluster_{}.dat".format(size,iteration,size),"lammps/data")

    