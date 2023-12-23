import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate, stats
no_task = 8
no_atoms = 270
lower = -1
upper =  1
for iteration in [3]:
    ns = np.zeros([no_task,2])
    for i in range(no_task):
        path = '/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_{}/Clusters/{}_cluster/cluster_{}/MC_data_{}'.format(i,no_atoms,iteration,iteration)
        ns_batch = np.loadtxt(path)
        ns[i,:] = ns_batch
    ns = np.insert(ns,4,[0,1],axis=0)
    np.savetxt('/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_trial/Clusters/{}_cluster/cluster_{}/MC_combined.txt'.format(no_atoms,iteration),ns)
    ####                                                              MC plot

    plt.figure(figsize=[6,4], dpi=300)
    plt.scatter(ns[:,0],ns[:,1], marker='.',s=150)
    plt.title(' MC Test (Linear fit)')
    plt.xlabel("Normalised normal stress")
    plt.ylabel("Normalised shear stress")
    slope, intercept, r, _, _ = stats.linregress(ns[:,0],ns[:,1])
    x=np.linspace(lower,upper,10)
    y=slope*x+intercept
    plt.plot(x,y,'r')
    plt.text(lower,1,'y={}x+{}'.format(round(slope,5),round(intercept,5)))
    plt.text(lower,0.9,'R\u00b2={}'.format(r**2))
    plt.savefig('/shared/LAMMPS/lammps/build_nep_BMG_generation/workflow/vit-1/Latest_full_shear_trial/Clusters/{}_cluster/cluster_{}/MC_test.png'.format(no_atoms,iteration), facecolor='white', transparent=False)
    plt.clf()


