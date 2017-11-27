import numpy as np

class ObserverProjection(object):

    def __init__(self,n_rays):

        if n_rays==12:
            self.read_flux_factors = read_flux_factors_12rays
        elif n_rays==24:
            self.read_flux_factors = read_flux_factors_24rays
        else:
            print("Error in the number of rays for the flux factor")
            exit(-1)

    def __call__(self,angle,delta):
        return self.read_flux_factors(angle,delta)


def read_flux_factors_12rays(angles,delta_angle):
    if (any(angle < 0. or angle > 90. for angle in angles)):
        print('Error, the observer angles must be between 0 and 90')
        exit(0)
    flux_factor1=np.array([np.loadtxt('./output_12_rays/flux_obs_ang_'+str(int(np.floor(angle)))+'.0_nrays_12.dat',unpack=True,usecols=([4])) for angle in angles])
    flux_factor2=np.array([np.loadtxt('./output_12_rays/flux_obs_ang_'+str(int(np.floor(angle+delta_angle)))+'.0_nrays_12.dat',unpack=True,usecols=([4])) for angle in angles])
    alpha = np.array([(angle - float(np.floor(angle)))/delta_angle for angle in angles])
    return np.array([(1.-a)*f1 + a*f2 for f1,f2,a in zip(flux_factor1,flux_factor2,alpha)])

if __name__=="__main__":
    M = ObserverProjection(12)
    ff = M([15.,26.],1.)
    print(ff)
