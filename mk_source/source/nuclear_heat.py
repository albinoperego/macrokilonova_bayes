import numpy as np
import units

class NuclearHeat(object):

    def __init__(self,eps_ye_dep):
        if eps_ye_dep:
            self.heat_rate = heat_rate_w_ye_dependence
        else:
            self.heat_rate = heat_rate_wo_ye_dependence

    def __call__(self,alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff,**kwargs):
        return self.heat_rate(alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff,**kwargs)

def calc_eps_nuc(kappa,time,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc):
  if (kappa < 1.):
    tmp=eps0 * yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time)
  elif (kappa > 10.):
    tmp=eps0
  else:
    tmp=0.5*eps0*(1 + yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time))
  return tmp 

def yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time):
  time_day=time*units.sec2day
  tmp = min(max(4*time_day-4.,-20),20)   # t_eps_nuc still missing!
  return a_eps_nuc + b_eps_nuc/(1.+np.exp(tmp))

def heat_rate_w_ye_dependence(alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff, **kwargs):
    a_eps_nuc = kwargs['cnst_a_eps_nuc']
    b_eps_nuc = kwargs['cnst_b_eps_nuc']
    t_eps_nuc = kwargs['cnst_t_eps_nuc']
    kappa     = kwargs['opacity']
    eps_nuc = calc_eps_nuc(kappa,t,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc)
    eps_th = ET(time_sec=t,mass_ej=m_ej,omega=Omega,vel=v_rms,cnst_eff=0.333)
    return eps_nuc*(0.5 - 1./np.pi * np.arctan((t-t0)/sigma0))**alpha * (eps_th/0.5)

def heat_rate_wo_ye_dependence(alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff):
    eps_nuc = eps_cnst_fac * eps0
    eps_th = ET(time_sec=t,mass_ej=m_ej,omega=Omega,vel=v_rms,cnst_eff=0.333)
    return eps_nuc*(0.5 - 1./np.pi * np.arctan((t-t0)/sigma0))**alpha * (eps_th/0.5)
