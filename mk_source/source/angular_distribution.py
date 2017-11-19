import numpy as np

class AngularDistribution(object):

    def __init__(self,angular_law):

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
    M = AngularDistribution("uniform")
    a,o = M(10)
    print a
    print o,np.sum(o)
    M = AngularDistribution("cos_uniform")
    a,o = M(10)
    print a
    print o,np.sum(o)
