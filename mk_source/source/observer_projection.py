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


def read_flux_factors_12rays(angle,delta_angle):
    if ( angle < 0. or angle > 90.):
        print('Error, the observer angles must be between 0 and 90')
        exit(0)
    flux_factor1=np.loadtxt('./output_12_rays/flux_obs_ang_'+str(int(np.floor(angle)))+'.0_nrays_12.dat',unpack=True,usecols=([4]))
    flux_factor2=np.loadtxt('./output_12_rays/flux_obs_ang_'+str(int(np.floor(angle+delta_angle)))+'.0_nrays_12.dat',unpack=True,usecols=([4]))
    alpha = (angle - float(np.floor(angle)))/delta_angle
    return (1.-alpha)*flux_factor1 + alpha*flux_factor2 

if __name__=="__main__":
    M = ObserverProjection(12)
    ff = M(26.,1.)
    print(ff)
