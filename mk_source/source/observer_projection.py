import numpy as np
from scipy import  interpolate

class ObserverProjection(object):

    def __init__(self, n_slices, dist_slices):
        
        self.n_rays = n_slices
        if (dist_slices == "uniform"):
            self.data_location = "output_"+str(self.n_rays)+"_rays"
        else:
            self.data_location = "output_cos_"+str(self.n_rays)+"_rays"

        self.flux_interpolant = self.read_flux_factors(self.n_rays)

    def __call__(self,angle):
        return np.array([f(angle) for f in self.flux_interpolant])

    def read_flux_factors(self, nrays):
        view_angles = np.linspace(0.0,90.0,91)
        flux_factors = np.array([np.loadtxt(self.data_location+'/ff_ray_%d.dat'%i,unpack=True,usecols=([1])) for i in range(nrays)])
        flux_interpolant = [interpolate.interp1d(view_angles,f) for f in flux_factors]
        return flux_interpolant

if __name__=="__main__":
    M = ObserverProjection(12)
    ff = M(31.0)
    print(ff)
