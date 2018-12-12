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
        return eps0 * yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time)
    elif (kappa > 10.):
        return eps0
    else:
        return 0.5*eps0*(1 + yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time))

def yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time):
    time_day=time*units.sec2day
    if (np.isscalar(time_day)):
        tmp = min(max(4*time_day-4.,-20),20)   # t_eps_nuc still missing!
    else:
        T1=np.zeros(len(time_day))
        T1[(4.*time_day-4.)>(-20)]=4.*time_day[(4.*time_day-4.)>(-20)]-4.
        T1[(4.*time_day-4.)<(-20)]=-20.
        tmp=np.zeros(len(time_day))
        tmp[T1<20.]=T1[T1<20.]
        tmp[T1>20.]=20.
    return a_eps_nuc + b_eps_nuc/(1.+np.exp(tmp))

def heat_rate_w_ye_dependence(alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff, **kwargs):
    a_eps_nuc = kwargs['cnst_a_eps_nuc']
    b_eps_nuc = kwargs['cnst_b_eps_nuc']
    t_eps_nuc = kwargs['cnst_t_eps_nuc']
    kappa     = kwargs['opacity']
    eps_nuc = calc_eps_nuc(kappa,t,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc)
    if(kwargs['normalize']):
        m_ej=m_ej/units.Msun
        v_rms=v_rms/units.c
    eps_th = ET(time_sec=t,mass_ej=m_ej,omega=Omega,vel=v_rms,cnst_eff=cnst_eff)
    return eps_nuc*(0.5 - units.oneoverpi * np.arctan((t-t0)/sigma0))**alpha * (2.*eps_th)  #here I am dividing eps_th by 0.5

def heat_rate_wo_ye_dependence(alpha,t,t0,sigma0,eps0,ET,m_ej,Omega,v_rms,cnst_eff, **kwargs):
    eps_nuc = eps0
    if(kwargs['normalize']):
        m_ej=m_ej/units.Msun
        v_rms=v_rms/units.c    
    eps_th = ET(time_sec=t,mass_ej=m_ej,omega=Omega,vel=v_rms,cnst_eff=cnst_eff)
    return eps_nuc*(0.5 - units.oneoverpi * np.arctan((t-t0)/sigma0))**alpha * (2.*eps_th)  #here I am dividing eps_th by 0.5

