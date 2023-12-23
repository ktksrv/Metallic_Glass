def extract_dat(path,no_atoms):                                                     
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=8, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    out=np.column_stack([ID,x,y,z])

    return out

def extract_dat_type(path,no_atoms):
                                                              
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=8, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    typ = np.array(df.type)
    out=np.column_stack([ID,x,y,z])
    return out,typ

def extract_dat(path,no_atoms):
                                                              
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=8, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    out=np.column_stack([ID,x,y,z])

    return out
    
def extract_dump(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    out=np.column_stack([ID,x,y,z])
    return out

def extract_dump_type(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    typ = np.array(df.type)
    out=np.column_stack([ID,x,y,z])
    return out,typ

def extract_fdump(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z', 'F_x', 'F_y', 'F_z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    fx=np.array(df.F_x)
    fy=np.array(df.F_y)
    fz=np.array(df.F_z)

    out_c=np.column_stack([x,y,z])                                 
    out_f=np.column_stack([fx,fy,fz]) 

    return out_c, out_f