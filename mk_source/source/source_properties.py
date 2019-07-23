import numpy as np
import units_mkn

class SourceProperties(object):

    def __init__(self,source_label):
        '''
        class to initialize the global properties of the source:
            D:                   distance [cm]
            filter_data_folder:  name of the folder with the filter data
            t0:                  universal time of the observation [days]
            view_angle:          viewing angle between the source and the observer [deg]
        '''

        if (source_label == 'default'):
            self.D = 40.e+6 * units_mkn.pc2cm
            self.view_angle = 90.
            self.filter_data_folder=None
            self.t0 = 0.e0

        elif (source_label == 'AT2017gfo'):
            self.D = 40.e+6 * units_mkn.pc2cm
            self.view_angle = 180/12.
            self.filter_data_folder='filter_data_AT2017gfo'
            self.t0 = 57982.529

        elif (source_label.__contains__('AT2017gfo') and source_label.__contains__('view_angle')):
            self.D = 40.e+6 * units_mkn.pc2cm
            self.view_angle = float(source_label.split('view_angle=')[-1])
            self.filter_data_folder='filter_data_AT2017gfo'
            self.t0 = 57982.529
            print("Source properies changed")
            # exit(1)
        else:
            print('unknown source!')
            print('please specify a source or default')
            exit()

    def init_time(self,tscale,time_min,time_max,n_time,mag):  
        if (tscale == 'linear'):
            time = np.linspace(time_min,time_max,num=n_time)
        elif (tscale == 'log'):
            time = np.logspace(np.log10(time_min),np.log10(time_max),num=n_time)
        elif (tscale == 'measures'):
            toll = 0.05
            all_time = []
            for ilambda in mag.keys():
                if (len(mag[ilambda]['time']>0)):
                    for i in range(len(mag[ilambda]['time'])):
                        all_time.append(mag[ilambda]['time'][i]-self.t0)
            all_time = sorted(np.array(all_time))
            time = []
            i = 0
            while (i < len(all_time)):
                delta = (1.+2.*toll) * all_time[i]
                i_start = i
                while (all_time[i] < delta):
                    i = i + 1
                time.append(0.5*(all_time[i]+all_time[i_start]))
                if (i == len(all_time)-1):
                    break
            time = sorted(np.array(time) * units_mkn.day2sec)
        else:
            print('Error! Wrong option for the time scale')
            exit(-1)
        return np.array(time)

    def time_measures(self,mag,toll):
        
        all_time = []
        for ilambda in mag.keys():
            if (len(mag[ilambda]['time']>0)):
                for i in range(len(mag[ilambda]['time'])):
                    all_time.append(mag[ilambda]['time'][i]-self.t0)
        all_time = sorted(np.array(all_time))
        time = []
        i = 0
        while (i < len(all_time)):
            delta = (1.+2.*toll) * all_time[i]
            i_start = i
            while (all_time[i] < delta):
                i = i + 1
            time.append(0.5*(all_time[i]+all_time[i_start]))
            if (i == len(all_time)-1):
                break
        time = sorted(np.array(time) * units_mkn.day2sec)
        return np.array(time)

