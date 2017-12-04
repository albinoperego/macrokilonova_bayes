import numpy as np
import units
from scipy import interpolate

class Thermalization(object):

    def __init__(self,therm_eff_model):
        if (therm_eff_model=='BKWM'):
            self.therm_efficiency = BKWM_therm_efficiency
        elif (therm_eff_model=='cnst'):
            self.therm_efficiency = cnst_therm_efficiency
        else:
            print('Unknown thermal efficiency model\n')
            exit(-1)
        
    def __call__(self, **kwargs):
        return self.therm_efficiency( **kwargs)

def therm_efficiency_params(m,omega,v):
# define the interpolation functions
    m_iso = 4.*np.pi/omega * m
    x = [np.log10(1.e-3),np.log10(5e-3),np.log10(1e-2),np.log10(5e-2)]
    y = [0.1,0.2,0.3]
    a = [[2.01,0.81,0.56,0.27],[4.52,1.90,1.31,0.55],[8.16,3.20,2.19,0.95]]
    b = [[0.28,0.19,0.17,0.10],[0.62,0.28,0.21,0.13],[1.19,0.45,0.31,0.15]]
    d = [[1.12,0.86,0.74,0.60],[1.39,1.21,1.13,0.90],[1.52,1.39,1.32,1.13]]
    fa = interpolate.interp2d(x, y, a, kind='linear')
    fb = interpolate.interp2d(x, y, b, kind='linear')
    fd = interpolate.interp2d(x, y, d, kind='linear')
# assign the values of the mass and velocity
    xnew=np.log10(m_iso)   #mass     [Msun]
    ynew=v                 #velocity [c]
# compute the parameters by linear interpolation in the table
    return [fa(xnew,ynew),fb(xnew,ynew),fd(xnew,ynew)]

def BKWM_therm_efficiency(**kwargs):
    t = kwargs['time_sec']   # s
    m = kwargs['mass_ej']    # m_sun
    omega = kwargs['omega']  # rad
    v = kwargs['vel']        # c
    if any( [t is None, m is None, omega is None, v is None] ):
        print('Error. For the heating efficiency, user must specify time, mass, angle and velocity\n')
        exit(-1)
    coeff=therm_efficiency_params(m,omega,v)
    time_days=t*units.sec2day
    tmp = float(2.*coeff[1]*time_days**coeff[2])
    return float(0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp ))

def cnst_therm_efficiency(**kwargs):
    thef = kwargs['cnst_eff']
    if (thef is None):
        print('Error, missing constant for the thermal efficiency\n')
        exit(-1)   
    if (thef < 0. or thef > 1.):
            print('Error, fixed heating efficiency must be between 0 and 1\n')
            exit(-1)
    return thef

if __name__=="__main__":
    M = Thermalization("cnst")
    x1 = M(cnst_eff=0.333)
    print(x1)
    M = Thermalization("BKWM")
    x1 = M(time_sec=2e5,mass_ej=0.1,omega=0.3,vel=0.02,cnst_eff=0.333)
    print(x1)
