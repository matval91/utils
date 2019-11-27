import utils.plot_utils as pu

col, col2, styles, my_cmap, dpi = pu.define_colors()

def plot_spaceE(tfbm, ax, color, ls, label):
    """
    plot 1D (pitch, int_space (int_E (fdist)))
    """
    x = tfbm.dict_dim['pitch']; y = tfbm.f_spaceE_int
    xlabel = r'$xi=(v_parallel/v)$'; ylabel = r'f'
    pu._plot_1d(x, y, xlabel, ylabel, ax=ax, color=color, ls=ls, label=label)
    
def plot_spacep(tfbm, ax, color, ls, label):
    """
    plot 1D (energy, int_space (int_pitch (fdist)))
    """  
    x = tfbm.dict_dim['E']
    if x[0]<1:
        x/=1.602e-19*1e-3
    y = tfbm.f_spacep_int
    xlabel = r'$E [keV]$'; ylabel = r'f [1/keV]'
    pu._plot_1d(x, y, xlabel, ylabel, ax=ax, color=color, ls=ls, label=label)

def plot_Epitch(tfbm, ax):
    """
    plot 2D (pitch, energy, int_space(fdist))
    """
    x,y = tfbm.dict_dim['pitch'], tfbm.dict_dim['E']*1e-3
    title = tfbm.runid + ' ' + str(tfbm.time)
    pu._plot_2d(x,y,  r'$\xi$', r'E [keV]', tfbm.f_space_int,title)

def plot_space_rz(tfbm, ax):
    """
    """
    x,y = tfbm.dict_dim['R'], tfbm.dict_dim['z']
    title = tfbm.runid + ' ' + str(tfbm.time)
    pu._plot_2d(x,y, r'R [m]', r'z [m]', dist=tfbm.f_Ep_int.T, title=title, \
                wallrz=[tfbm.R_w, tfbm.z_w], surf=[tfbm.rsurf, tfbm.zsurf])
