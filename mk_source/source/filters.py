import numpy as np
import filter_dictionary
import copy
from scipy import interpolate
import matplotlib.pyplot as plt
import dered_cardelli as drd
import units_mkn

# CORRECTION FOR REDDENING: if 'True' uses 'dered_cardelli' module to correct data for reddening.
dered_correction = True

# UPPER LIMITS: if 'True' it includes every data in the files inside 'filter_data' 
#               folder (also upper limits and uncertain data).
# Note: Upper limits are flagged with '-1' errors, uncertain data 
#       (that are presented in Villar table without errors but are not upper limits) with '123456789' errors.
upper_limits = False

class Filters(object):

    def __init__(self,usage):
        if (usage=="measures"):
            self.read_filters = read_filter_measures
        elif(usage=="properties"):
            self.read_filters = read_filter_properties
        else:
            print("Wrong usage for filters.")
            exit(0)
        self.dic_filt = None
        self.lambda_vec = None
        self.measures = None

    def __call__(self,folder):
        if (self.dic_filt is None) and (self.lambda_vec is None) and (self.measures is None):
            self.dic_filt, self.lambda_vec, self.measures = self.read_filters(folder)
        return self.dic_filt, self.lambda_vec, self.measures


def read_filter_properties(folder):
    # load the filter informations
    dic_filt = filter_dictionary.filters
    # sort by wavelength
    lambda_vec = copy.deepcopy(list(dic_filt.keys()))
    lambda_vec.sort()
    measures = {}
    return dic_filt,lambda_vec,measures


def read_measurements(filename):
    # print("reading: {}".format(filename))
    t, magnitude, magnitude_error = np.loadtxt(filename, unpack = True)
    return t, magnitude, magnitude_error


def read_filter_measures(folder_name):
    dic_filt,lambda_vec,dummy = read_filter_properties(folder_name)
# load the measured magnitudes
    measures={}
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==1):
            measures[ilambda] = {}
            t_tot = np.asarray([])
            m_tot = np.asarray([])
            sm_tot = np.asarray([])
            for fname in dic_filt[ilambda]["filename"]:

                t,m,sm = read_measurements(folder_name+"/"+fname)
                if(upper_limits):
                    t_tot = np.append(t_tot,t)
                    m_tot = np.append(m_tot,m)
                    sm_tot = np.append(sm_tot,sm)
                else:
                    if (type(sm) is np.float64):                  # load t,mag,err for a single data point
                        if(sm>0. and sm!=123456789):
                            t_tot = np.append(t_tot,t)
                            m_tot = np.append(m_tot,m)
                            sm_tot = np.append(sm_tot,sm)
                    else:
                        for i in range(0,len(sm)):                # load t,mag,err for data arrays (more than 1 data point)
                            if(sm[i]>0. and sm[i]!=123456789):
                                t_tot=np.append(t_tot,t[i])
                                m_tot=np.append(m_tot,m[i])
                                sm_tot=np.append(sm_tot,sm[i])
                if (dered_correction):                            # if dered_correction is 'True' it corrects the magnitudes [M_dered = M - correction]
                    m_tot -= drd.dered_CCM(np.asarray([ilambda]))

            measures[ilambda]['time']  = t_tot
            measures[ilambda]['mag']   = m_tot
            measures[ilambda]['sigma'] = sm_tot
            measures[ilambda]['name']  = dic_filt[ilambda]["name"]

    return dic_filt,lambda_vec,measures





def planckian(nu,T_plk):
    tmp = (units_mkn.h * nu) / (units_mkn.kB * T_plk)
    return (2. * units_mkn.h * nu ** 3) / (units_mkn.c2) / (np.exp(tmp) - 1.)

def m_filter(lam,T,rad,dist,ff):
    fnu = calc_fnu(lam,T,rad,dist,ff)
    return -2.5*np.log10(fnu)-48.6

def calc_fnu(lam,temp,rad,dist,ff):
  ff1 = ff[:len(ff)//2]
  ff2 = ff[len(ff)//2:]
  tmp1 = np.array([r * r * f * planckian(units_mkn.c / (100. * lam), T) for r, f, T in zip(rad, ff1, temp)])
  tmp2 = np.array([r * r * f * planckian(units_mkn.c / (100. * lam), T) for r, f, T in zip(rad[::-1], ff2, temp[::-1])])
  return np.sum(tmp1+tmp2)/(dist*dist)

def calc_lum_iso(llum,ff):
  tmp1 = np.array([lum * f for lum,f in zip(llum,ff[:len(ff)//2])])
  tmp2 = np.array([lum * f for lum,f in zip(llum[::-1],ff[len(ff)//2:])])
  return np.sum(tmp1+tmp2,axis=0)

def calc_magnitudes(ff,time,rad_ray,T_ray,lambda_vec,dic_filt,D,t0):

    ordered_T = np.asarray([list(x) for x in zip(*T_ray)])
    ordered_R = np.asarray([list(x) for x in zip(*rad_ray)])

    mag_model={}
    mag_model[0]=time
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==1):
            mag_model[ilambda] = np.array([m_filter(dic_filt[ilambda]["lambda"],T,R,D,ff) for T,R in zip(ordered_T,ordered_R)])
    return mag_model 

def calc_residuals(data,model,t0):
    res = {}
    for ilambda in data.keys():
        if ilambda == 0: 
            continue
#        fmag = interpolate.interp1d(model[0],model[ilambda], copy=False, bounds_error=None, fill_value=np.nan, assume_sorted=True)
        fmag = np.interp((data[ilambda]['time']-t0) * units_mkn.day2sec, model[0], model[ilambda])
#        res[ilambda] = np.array([(fmag-m)/sm   for t,m,sm in zip(data[ilambda]['time'],data[ilambda]['mag'],data[ilambda]['sigma'])])
        res[ilambda] = (fmag-data[ilambda]['mag'])  /data[ilambda]['sigma']
    return res



def calc_all_residuals(ff,time,rad_ray,T_ray,lambda_vec,dic_filt,D,t0,data):

    res = {}
    for ilambda in lambda_vec:
        if (dic_filt[ilambda]["active"]==0):
            continue
        if (len(data[ilambda]['time'])==0):
            continue

#  here time delay must be implemented
        R_intrp = np.array([np.interp((data[ilambda]['time']-t0) * units_mkn.day2sec, time, x) for x in rad_ray])
        T_intrp = np.array([np.interp((data[ilambda]['time']-t0) * units_mkn.day2sec, time, x) for x in T_ray])

        ordered_T = np.asarray([list(x) for x in zip(*T_intrp)])
        ordered_R = np.asarray([list(x) for x in zip(*R_intrp)])

        mag_model = np.array([m_filter(dic_filt[ilambda]["lambda"],x,y,D,ff) for x,y in zip(ordered_T,ordered_R) ])

        res[ilambda] = (mag_model-data[ilambda]['mag'])  /data[ilambda]['sigma']

    return res
   
if __name__=="__main__":

    FT = Filters("properties") 
    dic_filt,lambda_vec = FT()

    FT2 = Filters("measures") 
    dic_filt,lambda_vec,mag = FT2()

