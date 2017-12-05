import numpy as np
import units

#class NuclearHeat(object):

#    def __init__(self)
#        if

#    def __call__(self,)
#        return 

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

def heat_rate(alpha,t,t0,sigma0,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc,ET,m_ej,Omega,v_rms,kappa,cnst_eff):
    eps_nuc = calc_eps_nuc(kappa,t,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc)
    eps_th = ET(time_sec=t,mass_ej=m_ej,omega=Omega,vel=v_rms,cnst_eff=0.333)
    return eps_nuc*(0.5 - 1./np.pi * np.arctan((t-t0)/sigma0))**alpha * (eps_th/0.5)
