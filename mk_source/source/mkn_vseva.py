import angular_distribution as ad
import ejecta as ej
import filters as ft
import math
import matplotlib.pyplot as plt
import numpy as np
import observer_projection as op
import source_properties as sp
import units


class MKN(object):
    """
    Driver for the kilonova model
    """

    def __init__(self,
                 # number of different components of the ejecta
                 Nshell,
                 # dictionary of global parameters defining basic properties of the ejecta
                 glob_params,
                 # dictionary of global parameters defining basic properties of the ejecta to be sampled
                 glob_vars,
                 # dictionary of ejecta parameters defining its composition and geometry
                 ejecta_params,
                 # dictionary of shell parameters defining basic properties of the shell
                 shell_params,
                 source_name,
                 **kwargs):

        super(MKN,self).__init__(**kwargs)

        # initializing the global properties of the source
        print('I am initializing the global properties of the source')
        SP = sp.SourceProperties(source_name)
        self.D = SP.D
        self.view_angle = SP.view_angle

        # initialize the angular distribution
        print('I am initializing the angles')
        self.n_slices = glob_params['n slices']
        self.dist_slices = glob_params['dist slices']
        self.AD = ad.AngularDistribution(self.dist_slices,self.n_slices)
        self.angular_distribution, self.omega_distribution = self.AD(self.n_slices/2)   # due to the symmetry abount the equatorial plane, the number of independent slices is half

        # initialize the filters
        print('I am initializing the filters')
        if (source_name == 'default'):
            self.FT = ft.Filters("properties")
        else:
            self.FT = ft.Filters("measures")
        self.dic_filt,self.lambda_vec,self.mag = self.FT(SP.filter_data_folder)


        #initialize the time
        print('I am initializing the global time')
        self.time_min = glob_params['time min']     
        self.time_max = glob_params['time max']     
        self.n_time   = glob_params['n time']       
        self.tscale   = glob_params['scale for t']
        self.t0 = SP.t0
        if (self.tscale == "measures" and source_name=='default'):
            print('')
            print("no measures available to set the time (default option)")
            print("please use linear or log scale")
            exit() 

        self.time = SP.init_time(self.tscale,self.time_min,self.time_max,self.n_time,self.mag)

        #initialize the observer location
        print('I am initializing the observer orientation')
        self.FF = op.ObserverProjection(self.n_slices,self.dist_slices)
#        self.flux_factor = self.FF(SP.view_angle)

        # register the geometry of the ejecta and create an Ejecta object
        self.glob_params   = glob_params
        self.glob_vars     = glob_vars
        self.ejecta_params = ejecta_params
        self.shell_params  = shell_params
        self.ejecta        = ej.Ejecta(len(self.ejecta_params), self.ejecta_params.keys(), self.ejecta_params)
       
        print('I am initializing the components')
        self.E = ej.Ejecta(Nshell, self.ejecta_params.keys(), self.ejecta_params)

#############################
#  LIGHT CURVE CALCULATION  #
#############################

    def lightcurve(self,shell_vars):

        self.shell_vars = shell_vars 

        # compute the lightcurve
        r_ph, L_bol, T_eff = self.E.lightcurve(self.angular_distribution,
                                               self.omega_distribution,
                                               self.time,
                                               self.shell_vars,
                                               self.glob_vars,
                                               self.glob_params)

        return r_ph,L_bol,T_eff

    def compute_log_likelihood(self,residuals):
        # compute the likelihood
        logL = 0.
        for ilambda in residuals.keys():
            logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))
        
        return logL

    def compute_log_likelihood(self,residuals):
        # compute the likelihood
        logL = 0.
        for ilambda in residuals.keys():
            logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))
        
        return logL

    def log_likelihood(self,r_ph,T_eff):

        self.flux_factor = self.FF(self.view_angle)

        # compute the residuals
        if (source_name != 'default'):
            print('i.e., I am computing residuals')
            residuals = ft.calc_all_residuals(self.flux_factor,self.time,r_ph,T_eff,self.lambda_vec,self.dic_filt,self.D,self.t0,self.mag)

        # compute the likelihood
            print('and then I am computing the likelihood')
            logL = self.compute_log_likelihood(residuals)

#            logL = 0.
#            for ilambda in residuals.keys():
#                logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

            print('logL')
            print(logL)

        return logL

#########################
# write out the output  #
#########################

    def write_output(self,r_ph,T_eff):

        self.model_mag = ft.calc_magnitudes(self.flux_factor,self.time,r_ph,T_eff,self.lambda_vec,self.dic_filt,self.D,self.t0)

        file_output = 'mkn_model.txt'
        g = open(file_output,'w')

        print('file name:',file_output)

        g.write('%20s' %('time'))
        for ilambda in self.dic_filt.keys():
            g.write('%20s' %(self.dic_filt[ilambda]['name']))
        g.write('\n')
    
        for i in range(len(self.time)):
            g.write('%20s' %(self.time[i]))
            for ilambda in self.model_mag.keys():
                if (ilambda == 0):
                    continue
                g.write('%20s' %(self.model_mag[ilambda][i]))
            g.write('\n')

        g.close()


################################
# plot some of the lightcurves #
################################

    def plot_output_sep(self):

        fig1 = plt.figure()
        for ilambda in self.dic_filt.keys():
            if (self.dic_filt[ilambda]['plot'] !=1):
                continue
    
            plt.plot(self.time*units.sec2day,self.model_mag[ilambda])
    
            if (source_name != 'default'):
                if(len(self.mag[ilambda]['sigma'])!=0):
                    plt.errorbar(self.mag[ilambda]['time']-self.t0,self.mag[ilambda]['mag'],yerr=self.mag[ilambda]['sigma'],fmt='o')
            plt.title('%s' %self.dic_filt[ilambda]['name'])
            plt.xlim(0.1,10)
            plt.ylim(27,15)
            plt.show()
        

    def plot_output_tog(self):

        fig1 = plt.figure()
        for ilambda in self.dic_filt.keys():
            if (self.dic_filt[ilambda]['plot'] !=1):
                continue
            plt.plot(self.time*units.sec2day,self.model_mag[ilambda])
            if (source_name != 'default'):
                if(len(self.mag[ilambda]['sigma'])!=0):
                    plt.errorbar(self.mag[ilambda]['time']-self.t0,self.mag[ilambda]['mag'],yerr=self.mag[ilambda]['sigma'],fmt='o')
        plt.title('Lightcurves')
        plt.xlim(0.1,10)
        plt.ylim(27,15)
        plt.show()


if __name__=='__main__':

#dictionary with the global parameters of the model
    glob_params = {'lc model'   :'grossman',    # model for the lightcurve (grossman or villar)    !fixed Grossman et al 2014, Villar et al 2017
                   'v_min'      :1.e-7,         # minimal velocity for the Grossman model          !fixed
                   'n_v'        :400,           # number of points for the Grossman model          !fixed
                   'vscale'     :'linear',      # scale for the velocity in the Grossman model     !fixed
                   'sigma0'     :0.11,          # parameter for the nuclear heating rate           !fixed Korobkin et al 2012
                   'alpha'      :1.3,           # parameter for the nuclear heating rate           !fixed Korobkin et al 2012
                   't0eps'      :1.3,           # parameter for the nuclear heating rate           !fixed Korobkin et al 2012
                   'cnst_eff'   :0.3333,        # parameter for the constant heating efficiency    !fixed
                   'n slices'   :30,            # number for the number of slices along the polar angle [12,18,24,30]   !fixed
                   'dist slices':'cos_uniform', # discretization law for the polar angle [uniform or cos_uniform]       !fixed
                   'time min'   :3600.,         # minimum time [s]                                  !fixed
                   'time max'   :2000000.,      # maximum time [s]                                  !fixed
                   'n time'     :200,           # integer number of bins in time                    !fixed
                   'scale for t':'linear'       # kind of spacing in time [log - linear - measures] !fixed
                   }
    
    source_name = 'AT2017gfo'   # name of the source or "default"   !it could be changed! 
    # source_name = 'default'   # name of the source or "default"
    
    # dictionary for the parameters of each shell
    ejecta_params = {}
    ejecta_params['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM', 'eps_ye_dep':True}  #fixed
    ejecta_params['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM', 'eps_ye_dep':True}  #fixed
    ejecta_params['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}  #fixed
    # BKWM_1d

#   BKWM         Barnes et al ApJ 2016
#   eps_ye_dep   Perego et al ApJL 2017 
 
    
    # dictionary for the variables of each shell
    shell_vars={}
    
    shell_vars['dynamics'] = {'xi_disk':None,
                              'm_ej':0.04,           # dynamical ejecta mass [Msun], to be changed!
                              'central_vel':0.33,    # dynamical ejecta mean velocity [c], to be changed!
                              'low_lat_vel':None,
                              'high_lat_vel':None,
                              'step_angle_vel':None,
                              'low_lat_op':30.,      # equatorial opacity, usually for Lanthanide rich matter, > 10. [cm^2/g]
                              'high_lat_op':1.0,     # polar opacity, usually Lanthanide-free, ~ 1. [cm^2/g]
                              'step_angle_op':math.radians(45.)}  # fixed
    
    shell_vars['wind'] = {'xi_disk':0.05,            # nu-driven wind simulations says 0.01-0.05, but B field can increase it (<~ 0.1)
                          'm_ej':None,
                          'step_angle_mass':math.radians(60.),  #fixed
                          'high_lat_flag':True,      # fixed
                          'central_vel':0.08,        # 0.04 c < v < 0.1
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':5.0,         # 1-5
                          'high_lat_op':0.5,        # 0.1-1
                          'step_angle_op':math.radians(30.)}   # fixed
    
    shell_vars['secular'] = {'xi_disk':0.2,         # 0.1-0.4
                             'm_ej':None,
                             'central_vel':0.06,    # 0.04 < v < 0.1
                             'low_lat_vel':None,
                             'high_lat_vel':None,
                             'step_angle_vel':None,
                             'central_op':5.0,      # 1-10
                             'low_lat_op':None,
                             'high_lat_op':None,
                             'step_angle_op':None}
    
    # dictionary for the global variables
    glob_vars = {'m_disk':0.12,    # to be changed (1.e-3,0.2) Msun
                 'eps0':1.5e19,    # it could be changed (4e18 - 2e19)
                 'T_floor_LA':1000.,   # fixed (Villar et al 2017)
                 'T_floor_Ni':3000.,   # fixed (Villar et al 2017)
                 'a_eps_nuc':0.5,      # fixxed, Ye correction to heating rate (Perego et al ApJL 2017)
                 'b_eps_nuc':2.5,      # ''
                 't_eps_nuc':1.0}      # ''
    
    print('I am initializing the model')
    model = MKN(3,glob_params, glob_vars, ejecta_params, shell_vars,source_name)

    print('I am computing the light curves')
    r_ph,L_bol,T_eff = model.lightcurve(shell_vars) 

    print('I am computing the likelihood')
    logL =  model.log_likelihood(r_ph,T_eff)

    write_output = True
    if (write_output):
        print('I am printing out the output')
        model.write_output(r_ph,T_eff)

    plot_separately = True            # Choose to plot all lightcurves in different bands on the same plot
    plot_together = False             # or to plot lightcurve and data in each band on different plots

    if (plot_separately):
        model.plot_output_sep()
    elif (plot_together):
        model.plot_output_tog()
