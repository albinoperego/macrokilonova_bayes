import numpy as np
import units

class ExpansionModelSingleSpherical(object):

    def __init__(self,expansion_model):
        if (expansion_model=='GK'):
            self.expansion_model_single_spherical = GK_expansion_model
        else:
            print('No other expansion model presently implemented\n')
            exit(-1)
        
    def __call__(self,omega,m_ej,v_rms,v_min,n_v,vscale,kappa, **kwargs):
        return self.expansion_model_single_spherical(omega,m_ej,v_rms,v_min,n_v,vscale,kappa, **kwargs)

def mass_gt_v(v,mej,v_exp):
  return mej*units.Msun*(1.0+func_vel(v/v_exp))  #[g]

def func_vel(x):
  return 35.*x**7/112.-105.*x**5/80.+35.*x**3/16.-35.*x/16.

def t_diff_v(kappa,v,m_v,omega):
  return np.sqrt(kappa*m_v/(omega*v*units.c*units.c))  #[s]

def t_fs_v(kappa,v,m_v,omega):
  return np.sqrt(1.5*kappa*m_v/(omega*v*v*units.c*units.c))  #[s]

def GK_expansion_model(Omega,m_ej,v_rms,v_min,n_v,vscale,kappa,**kwargs):

    v_max  = 3.*v_rms
    if (vscale == 'linear'):
        vel    = np.linspace(v_min,v_max,n_v)
    elif (vscale == 'log'):
        vel    = np.logspace(np.log10(v_min),np.log10(v_max),n_v)
    else:
        print('Wrong scale for expansion velocity')
        exit(0)
    m_vel  = mass_gt_v(vel,m_ej,v_max)
    t_diff = t_diff_v(kappa,vel,m_vel,Omega)
    t_fs   = t_fs_v(kappa,vel,m_vel,Omega)
    return vel,m_vel,t_diff,t_fs


################

if __name__=="__main__":
    M = ExpansionModelSingleSpherical("GK")
    Omega = 0.2
    m_ej = 0.1
    v_rms = 0.2
    v_min = 1.e-6
    n_v = 100
    kappa = 10.
    vscale = 'linear'
    x1,x2,x3,x4 = M(Omega,m_ej,v_rms,v_min,n_v,vscale,kappa)
    from pylab import *
    plot(x1,x2)
    show()
    print(x1)
    print(x2)
    print(x3)
    print(x4)
