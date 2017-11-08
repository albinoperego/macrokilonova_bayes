import numpy as np

class Filter(object):

    def __init__(self, name, nu_central, width, filename=None):

        self.name = name
        self.nu_central = nu_central
        self.width = width
        if filename is not None:
            self.t, self.magnitude, self.magnitude_error = self.read_measurements(filename)
    
    def read_measurements(self, filename):
        self.t, self.magnitude, self.magnitude_error = np.loadtxt(filename, unpack = True)
        return self.t, self.magnitude, self.magnitude_error

    def residuals(self, model):
        self.res = (model-self.magnitude)/self.magnitude_error
        return self.res

if __name__=="__main__":
    names = ['K','M']
    freqs = [500,700]
    widths = [10,20]

    misure = [Filter(n,l,w) for n,l,w in zip(names,freqs,widths)]
    logL = -0.5*np.sum([m.residuals(modello)**2 for m in misure])
