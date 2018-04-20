import numpy as np
import units
from scipy import interpolate

class Thermalization(object):

    def __init__(self,therm_eff_model):
        if (therm_eff_model=='BKWM'):
            self.therm_efficiency = BKWM_therm_efficiency
        # 2-D
            x = [np.log10(1.e-3),np.log10(5e-3),np.log10(1e-2),np.log10(5e-2)]
            y = [0.1,0.2,0.3]
            a = [[2.01,0.81,0.56,0.27],[4.52,1.90,1.31,0.55],[8.16,3.20,2.19,0.95]]
            b = [[0.28,0.19,0.17,0.10],[0.62,0.28,0.21,0.13],[1.19,0.45,0.31,0.15]]
            d = [[1.12,0.86,0.74,0.60],[1.39,1.21,1.13,0.90],[1.52,1.39,1.32,1.13]]
        # define the interpolation functions
            self.fa = interpolate.interp2d(x, y, a, kind='linear')
            self.fb = interpolate.interp2d(x, y, b, kind='linear')
            self.fd = interpolate.interp2d(x, y, d, kind='linear')
        elif (therm_eff_model=='BKWM_1d'):
            self.therm_efficiency = BKWM_1d_therm_efficiency
        # 1-D
        # Uses a single variable [m_ej/(v_ej^2)] instead of using m_ej and v_ej separately.
            x_barnes = [0.011,0.025,0.0556,0.1,0.111,0.125,0.25,0.5,0.5556,1.,1.25,5.]
            a_barnes = [8.16,4.52,3.20,2.01,2.19,1.90,1.31,0.81,0.95,0.56,0.55,0.27]
            b_barnes = [1.19,0.62,0.45,0.28,0.31,0.28,0.21,0.19,0.15,0.17,0.13,0.10]
            d_barnes = [1.52,1.39,1.39,1.12,1.32,1.21,1.13,0.86,1.13,0.74,0.90,0.60]
        # define the interpolation functions
            self.fa_1d=interpolate.interp1d(x_barnes,a_barnes,fill_value='extrapolate')
            self.fb_1d=interpolate.interp1d(x_barnes,b_barnes,fill_value='extrapolate')
            self.fd_1d=interpolate.interp1d(x_barnes,d_barnes,fill_value='extrapolate') 
        elif (therm_eff_model=='cnst'):
            self.therm_efficiency = cnst_therm_efficiency
        else:
            print('Unknown thermal efficiency model\n')
            exit(-1)
        
    def __call__(self, **kwargs):
        return self.therm_efficiency(self, **kwargs)

    def therm_efficiency_params(self, m,omega,v):
        m_iso = units.fourpi/omega * m
        # assign the values of the mass and velocity
        xnew=np.log10(m_iso)   #mass     [Msun]
        ynew=v                 #velocity [c]
        # compute the parameters by linear interpolation in the table
        return (self.fa(xnew,ynew),self.fb(xnew,ynew),self.fd(xnew,ynew))

    def therm_efficiency_params_1d(self, m,omega,v):
        m_iso = units.fourpi/omega * m
        # assign the value of x=m/v^2
        xnew = m_iso/(v*v)
        # compute the parameters by 1-d interpolation
        return (self.fa_1d(xnew),self.fb_1d(xnew),self.fd_1d(xnew))


def BKWM_therm_efficiency(cls, **kwargs):
    t = kwargs['time_sec']   # s
    m = kwargs['mass_ej']    # m_sun
    omega = kwargs['omega']  # rad
    v = kwargs['vel']        # c
    if any( [t is None, m is None, omega is None, v is None] ):
        print('Error. For the heating efficiency, user must specify time, mass, angle and velocity\n')
        exit(-1)
    coeff=cls.therm_efficiency_params(m,omega,v)
    time_days=t*units.sec2day
    if(np.isscalar(time_days)):
        tmp = float(2.*coeff[1]*time_days**coeff[2])
        return float(0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp ))
    else:
        tmp = (2.*coeff[1]*time_days**coeff[2])
        return (0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp ))

def BKWM_1d_therm_efficiency(cls, **kwargs):
    t = kwargs['time_sec']   # s
    m = kwargs['mass_ej']    # m_sun
    omega = kwargs['omega']  # rad
    v = kwargs['vel']        # c
    if any( [t is None, m is None, omega is None, v is None] ):
        print('Error. For the heating efficiency, user must specify time, mass, angle and velocity\n')
        exit(-1)
    coeff=cls.therm_efficiency_params_1d(m,omega,v)
    time_days=t*units.sec2day
    if(np.isscalar(time_days)):
        tmp = float(2.*coeff[1]*time_days**coeff[2])
        return float(0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp ))
    else:
        tmp = (2.*coeff[1]*time_days**coeff[2])
        return (0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp ))        
        

def cnst_therm_efficiency(cls,  **kwargs):
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
