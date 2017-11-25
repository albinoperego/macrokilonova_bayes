import numpy as np
import filter_dictionary
import copy

class Filter(object):

#    def __init__(self, name, nu_central, width, filename=None):

#        self.name = name
#        self.nu_central = nu_central
##        self.width = width
#        if filename is not None:
#            self.t, self.magnitude, self.magnitude_error = self.read_measurements(filename)
    
    def read_measurements(self, filename):
        self.t, self.magnitude, self.magnitude_error = np.loadtxt(filename, unpack = True)
        return self.t, self.magnitude, self.magnitude_error

    def residuals(self, model):
        self.res = (model-self.magnitude)/self.magnitude_error
        return self.res

if __name__=="__main__":

    dic_filt = filter_dictionary.filters
    lambda_vec = copy.deepcopy(dic_filt.keys())
    lambda_vec.sort()

    f = Filter()


    measures={}
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==1):
           measures[ilambda] = {}
           t_tot = np.asarray([])
           m_tot = np.asarray([])
           sm_tot = np.asarray([])
           for fname in dic_filt[ilambda]["filename"]:
               print(fname)
               t,m,sm = f.read_measurements("filter_data/"+fname)
               print(t)
               t_tot = np.append(t_tot,t)
               m_tot = np.append(m_tot,m)
               sm_tot = np.append(sm_tot,sm)

           measures[ilambda]['time']  = t_tot
           measures[ilambda]['mag']   = m_tot
           measures[ilambda]['sigma'] = sm_tot
           measures[ilambda]['name']  = dic_filt[ilambda]["name"]
               
    print(measures)           

#    misure = [Filter(n,l,w) for n,l,w in zip(names,freqs,widths)]
#    logL = -0.5*np.sum([m.residuals(modello)**2 for m in misure])
