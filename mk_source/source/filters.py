import numpy as np
import filter_dictionary
import copy
import units
from scipy import interpolate

class Filters(object):

    def __init__(self,usage):
        if (usage=="measures"):
            self.read_filters = read_filter_measures
        elif(usage=="properties"):
            self.read_filters = read_filter_properties
        else:
            print("Wrong usage for filters.")
            exit(0)

    def __call__(self):
        return self.read_filters()

def read_filter_measures():

    dic_filt,lambda_vec = read_filter_properties()

# load the measured magnitudes
    measures={}
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==1):
            measures[ilambda] = {}
            t_tot = np.asarray([])
            m_tot = np.asarray([])
            sm_tot = np.asarray([])
            for fname in dic_filt[ilambda]["filename"]:
                t,m,sm = read_measurements("filter_data/"+fname)
                t_tot = np.append(t_tot,t)
                m_tot = np.append(m_tot,m)
                sm_tot = np.append(sm_tot,sm)

            measures[ilambda]['time']  = t_tot
            measures[ilambda]['mag']   = m_tot
            measures[ilambda]['sigma'] = sm_tot
            measures[ilambda]['name']  = dic_filt[ilambda]["name"]

    return dic_filt,lambda_vec,measures


def read_filter_properties():

# load the filter informations 
    dic_filt = filter_dictionary.filters

# sort by wavelength
    lambda_vec = copy.deepcopy(dic_filt.keys())
    lambda_vec.sort()

    return dic_filt,lambda_vec

def planckian(nu,T_plk):
    tmp = (units.h*nu)/(units.kB*T_plk)
    return (2.*units.h*nu**3)/(units.c**2)/(np.exp(tmp)-1.)
 
def read_measurements(filename):
    t, magnitude, magnitude_error = np.loadtxt(filename, unpack = True)
    return t, magnitude, magnitude_error

def m_filter(lam,T,rad,dist,ff):
    fnu = calc_fnu(lam,T,rad,dist,ff)
    return -2.5*np.log10(fnu)-48.6

def calc_fnu(lam,temp,rad,dist,ff):
  ff1 = ff[:len(ff)/2]
  ff2 = ff[len(ff)/2:]
  tmp1 = np.array([r**2 * f * planckian(units.c/(100.*lam),T) for r,f,T in zip(rad,ff1,temp)])
  tmp2 = np.array([r**2 * f * planckian(units.c/(100.*lam),T) for r,f,T in zip(rad[::-1],ff2,temp[::-1])])
  return np.sum(tmp1+tmp2)/dist**2

def calc_magnitudes(ff,time,rad_ray,T_ray,lambda_vec,dic_filt,D):

    ordered_T = np.asarray([list(x) for x in zip(*T_ray)])
    ordered_R = np.asarray([list(x) for x in zip(*rad_ray)])

    mag_model={}
    mag_model[0]=time
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==1):
            mag_model[ilambda] = np.array([m_filter(dic_filt[ilambda]["lambda"],T,R,D,ff) for T,R in zip(ordered_T,ordered_R)])

    return mag_model 

def calc_residuals(data,model):
    res = {}
    for ilambda in data.keys():
        if ilambda == 0: 
            continue
        fmag = interpolate.interp1d(model[0],model[ilambda])
        res[ilambda] = np.array([(fmag(t)-m)/sm   for t,m,sm in zip(data[ilambda]['time'],data[ilambda]['mag'],data[ilambda]['sigma'])])
      
    return res


if __name__=="__main__":

    FT = Filters("properties") 
    dic_filt,lambda_vec = FT()

    FT2 = Filters("measures") 
    dic_filt,lambda_vec,mag = FT2()

