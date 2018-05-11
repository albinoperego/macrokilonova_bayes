import numpy as np
import filter_dictionary
import copy
from scipy import interpolate
import matplotlib.pyplot as plt
import dered_cardelli as drd
import units


dered_correction = True
upper_limits = True

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
    def __call__(self):
        if (self.dic_filt is None) and (self.lambda_vec is None) and (self.measures is None):
            self.dic_filt, self.lambda_vec, self.measures = self.read_filters()
        return self.dic_filt, self.lambda_vec, self.measures


def read_filter_properties():
    # load the filter informations
    dic_filt = filter_dictionary.filters
    # sort by wavelength
    lambda_vec = copy.deepcopy(list(dic_filt.keys()))
    lambda_vec.sort()
    return dic_filt,lambda_vec


def read_measurements(filename):
    t, magnitude, magnitude_error = np.loadtxt(filename, unpack = True)
    return t, magnitude, magnitude_error


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
                print (fname)
#
                t,m,sm = read_measurements("filter_data_AT2017gfo/"+fname)
                if(upper_limits):
                    t_tot = np.append(t_tot,t)
                    m_tot = np.append(m_tot,m)
                    sm_tot = np.append(sm_tot,sm)
                else:
                    if (type(sm) is np.float64):
                        if(sm>0. and sm!=123456789):
                            t_tot = np.append(t_tot,t)
                            m_tot = np.append(m_tot,m)
                            sm_tot = np.append(sm_tot,sm)
                    else:
                        for i in range(0,len(sm)):
                            if(sm[i]>0. and sm[i]!=123456789):
                                t_tot=np.append(t_tot,t[i])
                                m_tot=np.append(m_tot,m[i])
                                sm_tot=np.append(sm_tot,sm[i])
#
#                print(ilambda)
                if (dered_correction):
                    m_tot -= drd.dered_CCM(np.asarray([ilambda]))
                for i in range(0,len(m_tot)):
                    if(sm_tot[i]==123456789):
                        plt.scatter(t_tot[i]-57982.529,m_tot[i],marker='x',c='black',s=15)
                    elif(sm_tot[i]<0.):
                        plt.scatter(t_tot[i]-57982.529,m_tot[i],marker='v',c='red',s=15)
                    else:
                        plt.errorbar(t_tot[i]-57982.529,m_tot[i],yerr=sm_tot[i],marker='o',c='blue',markersize=3)
                plt.ylim([15,28])
                plt.xlim([0,20])
                plt.gca().invert_yaxis()
                plt.title('%s' %fname)
                plt.show()
#
            measures[ilambda]['time']  = t_tot
            measures[ilambda]['mag']   = m_tot
            measures[ilambda]['sigma'] = sm_tot
            measures[ilambda]['name']  = dic_filt[ilambda]["name"]
    return dic_filt,lambda_vec,measures

read_filter_measures()

