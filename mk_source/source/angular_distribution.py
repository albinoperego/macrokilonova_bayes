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
            self.angular_distribution = uniform_ang
        elif angular_law=="cos_uniform":
            self.angular_distribution = cos_uniform_ang
        else:
            print("Unknown angular distribution")
            exit(0)
        
    def __call__(self,n):
        return self.angular_distribution(n)

def uniform_ang(n):
    delta = np.pi/2./float(n)
    a = np.array([ [delta*i,delta*(i+1)] for i in range(n)])
    o = np.array([ 2.*np.pi*(np.cos(x[0]) - np.cos(x[1])) for x in a])
    return a,o

def cos_uniform_ang(n):
    delta = 1./float(n)
    a = np.array([ [np.arccos(delta*i),np.arccos(delta*(i-1))] for i in range(n,0,-1)])
    o = np.array([ 2.*np.pi*(np.cos(x[0]) - np.cos(x[1])) for x in a])
    return a,o
                          
if __name__=="__main__":
    N = 12
    M = AngularDistribution("uniform", 12)
    a,o = M(10)
    print a
    print o,np.sum(o)
    M = AngularDistribution("cos_uniform", 12)
    a,o = M(10)
    print a
    print o,np.sum(o)
