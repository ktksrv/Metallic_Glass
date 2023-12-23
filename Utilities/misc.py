def box_coordinates (def_initial):
    margin = 2 # based on type of atom and its unit cell
    x_min = def_initial[:,0].min()-margin
    x_max = def_initial[:,0].max()+margin
    y_min = def_initial[:,1].min()-margin
    y_max = def_initial[:,1].max()+margin
    z_min = def_initial[:,2].min()-margin
    z_max = def_initial[:,2].max()+margin
    return x_min,x_max,y_min,y_max,z_min,z_max
def moving_average(a, n):
    import numpy as np
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def segments_fit(X, Y, count):
    from scipy import optimize
    import numpy as np
    import pylab as pl
    xmin = X.min()
    xmax = X.max()

    seg = np.full(count - 1, (xmax - xmin) / count)

    px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
    py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.01].mean() for x in px_init])

    def func(p):
        seg = p[:count - 1]
        py = p[count - 1:]
        px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        return px, py

    def err(p):
        px, py = func(p)
        Y2 = np.interp(X, px, py)
        return np.mean((Y - Y2)**2)

    r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    return func(r.x)

