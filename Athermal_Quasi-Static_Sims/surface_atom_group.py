def surface_atoms(initial,no_atoms):
    from ovito.data import Particles, DataCollection, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import ConstructSurfaceModifier
    import numpy as np
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    cell[...] = [[200,0,0,0],                               #use cell vectors for given congifuration from ovito
                [0,72.56,0,0],
                [0,0,145.12,0]]
    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial)
    id=range(1,no_atoms+1)
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    pipeline.modifiers.append(ConstructSurfaceModifier(method = ConstructSurfaceModifier.Method.AlphaShape, radius=2.55, select_surface_particles=True))
    data= pipeline.compute()
    selection=np.array(data.particles["Selection"])
    surface_ID_all=selection*id
    surface_ID=surface_ID_all[surface_ID_all !=0]
    
    return  surface_ID

def init_surface_atom_str(no_atoms,iteration):
    import STZ
    initial = STZ.extract_dat('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    d_initial = initial[:,1:]
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "    

    surface_group='\ngroup surface id {}'.format(surface_str)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "w")
    text_file.write(surface_group)
    text_file.close()