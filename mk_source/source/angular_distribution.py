import numpy as np

class AngularDistribution(object):
    """
    angular distribution class
    accepts uniform in the angle
    and uniform in the cosine of the angle
    """
    def __init__(self, angular_law, n_slices):

        if n_slices % 2 == 1:
            print('Error: n_slices must be an even number!')
            exit(-1)
        if angular_law=="uniform":
            self.angular_distribution = self.uniform_ang
        elif angular_law=="cos_uniform":
            self.angular_distribution = self.cos_uniform_ang
        else:
            print("Unknown angular distribution")
            exit(0)
        
    def __call__(self,n,omega_frac):
        return self.angular_distribution(n,omega_frac)

    def uniform_ang(self,n,omega_fraction):
        delta = np.pi/2./float(n)
        a = np.array([ [delta*i,delta*(i+1)] for i in range(int(n))])
        o = np.array([ 2.*np.pi*(np.cos(x[0]) - np.cos(x[1])) for x in a])
        return a,o*omega_fraction

    def cos_uniform_ang(self,n,omega_fraction):
        delta = 1./float(n)
        a = np.array([ [np.arccos(delta*float(i)),np.arccos(delta*float(i-1))] for i in range(int(n),0,-1)])
        o = np.array([ 2.*np.pi*(np.cos(x[0]) - np.cos(x[1])) for x in a])
        return a,o*omega_fraction
                          
if __name__=="__main__":
    N = 12
    M = AngularDistribution("uniform", 12)
    a,o = M(10)
    M = AngularDistribution("cos_uniform", 12)
    a,o = M(10)
