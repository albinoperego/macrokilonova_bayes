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
                 glob_params,
                 # dictionary of global parameters defining basic properties of the ejecta to be sampled
                 glob_vars,
                 # dictionary of ejecta parameters defining its composition and geometry
                 ejecta_params,
                 # dictionary of shell parameters defining basic properties of the shell
                 ejecta_vars,
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
        self.angular_distribution, self.omega_distribution = self.AD(self.n_slices/2,glob_params['omega frac'])   # due to the symmetry abount the equatorial plane, the number of independent slices is half

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
        self.ejecta_vars    = ejecta_vars
#        self.ejecta        = ej.Ejecta(len(self.ejecta_params),self.ejecta_params.keys(),self.ejecta_params)
       
        print('I am initializing the components')
        self.E = ej.Ejecta(len(self.ejecta_params.keys()), self.ejecta_params.keys(), self.ejecta_params)

#############################
#  LIGHT CURVE CALCULATION  #
#############################

#    def lightcurve(self,ejecta_vars,NR_data,NR_filename):

#        self.ejecta_vars = ejecta_vars 

        # compute the lightcurve
#        r_ph, L_bol, T_eff = self.E.lightcurve(self.angular_distribution,
#                                               self.omega_distribution,
#                                               self.time,
#                                               self.ejecta_vars,
#                                               self.glob_vars,
#                                               self.glob_params,
#                                               NR_data,
#                                               NR_filename)

#        return r_ph,L_bol,T_eff

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

        else:
            logL = 0.

        return logL

#########################
# write out the output  #
#########################

    def write_output(self,r_ph,T_eff,L_bol):

        self.model_lum = ft.calc_lum_iso(L_bol,self.flux_factor)

        self.model_mag = ft.calc_magnitudes(self.flux_factor,self.time,r_ph,T_eff,self.lambda_vec,self.dic_filt,self.D,self.t0)

        file_output = 'mkn_model.txt'
        g = open(file_output,'w')

        print('file name:',file_output)

        g.write('%20s %20s' %('time','luminosity'))
        for ilambda in self.dic_filt.keys():
            g.write('%20s' %(self.dic_filt[ilambda]['name']))
        g.write('\n')
        g.write('%20s %20s' %('[s]','[erg/s]'))
        g.write('\n')
    
        for i in range(len(self.time)):
            g.write('%20s %20s' %(self.time[i],self.model_lum[i]))
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
    glob_params = {'lc model'   :'grossman',    # model for the lightcurve (grossman or villar)  
                   'mkn model'  :'iso1comp',    # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                   'omega frac' :0.5,           #
                   'rad shell'  :False,         #
                   'v_min'      :1.e-7,         # minimal velocity for the Grossman model
                   'n_v'        :400,           # number of points for the Grossman model
                   'vscale'     :'linear',      # scale for the velocity in the Grossman model
                   'sigma0'     :0.11,          # parameter for the nuclear heating rate
                   'alpha'      :1.3,           # parameter for the nuclear heating rate
                   't0eps'      :1.3,           # parameter for the nuclear heating rate
                   'cnst_eff'   :0.3333,        # parameter for the constant heating efficiency
                   'n slices'   :30,            # number for the number of slices along the polar angle [12,18,24,30]
                   'dist slices':'cos_uniform', # discretization law for the polar angle [uniform or cos_uniform]
                   'time min'   :3600.,         # minimum time [s]
                   'time max'   :2000000.,      # maximum time [s]
                   'n time'     :200,           # integer number of bins in time
                   'scale for t':'log',    # kind of spacing in time [log - linear - measures]
                   'NR_data'    :False,         # use (True) or not use (False) NR profiles
                   'NR_filename':'../example_NR_data/DD2_M125125_LK/outflow_1/ejecta_profile.dat'           # path of the NR profiles, necessary if NR_data is True
                   }
    
    #source_name = 'AT2017gfo'   # name of the source or "default"
    source_name = 'default'   # name of the source or "default"
    
    # dictionary for the global variables
    glob_vars = {'m_disk':0.12,
                 'eps0':1.5e19, 
                 'T_floor_LA':1000., 
                 'T_floor_Ni':3500., 
                 'a_eps_nuc':0.5,
                 'b_eps_nuc':2.5,
                 't_eps_nuc':1.0}

###############################
# Template for isotropic case # 
###############################
    
    # hardcoded ejecta geometric and thermal parameters for the spherical case
    ejecta_params_iso = {}
    ejecta_params_iso['dynamics'] = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'uniform'}
    ejecta_params_iso['wind']     = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
    ejecta_params_iso['secular']  = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}

    # set of shell parameters to be sampled on
    ejecta_vars_iso={}

    ejecta_vars_iso['dynamics'] = {'xi_disk'        :None,
                                  'm_ej'           :0.04,
                                  'step_angle_mass':None,
                                  'high_lat_flag'  :None,
                                  'central_vel'    :0.24,
                                  'high_lat_vel'   :None,
                                  'low_lat_vel'    :None,
                                  'step_angle_vel' :None,
                                  'central_op'     :30.,
                                  'high_lat_op'    :None,
                                  'low_lat_op'     :None,
                                  'step_angle_op'  :None}

    ejecta_vars_iso['secular'] = {'xi_disk'        :None,
                                 'm_ej'           :0.05,
                                 'step_angle_mass':None,
                                 'high_lat_flag'  :None,
                                 'central_vel'    :0.06,
                                 'high_lat_vel'   :None,
                                 'low_lat_vel'    :None,
                                 'step_angle_vel' :None,
                                 'central_op'     :5.0,
                                 'low_lat_op'     :None,
                                 'high_lat_op'    :None,
                                 'step_angle_op'  :None}

    ejecta_vars_iso['wind'] = {'xi_disk'        :None,
                              'm_ej'           :0.02,
                              'step_angle_mass':None,
                              'high_lat_flag'  :True,
                              'central_vel'    :0.08,
                              'high_lat_vel'   :None,
                              'low_lat_vel'    :None,
                              'step_angle_vel' :None,
                              'central_op'     :1.0,
                              'high_lat_op'    :None,
                              'low_lat_op'     :None,
                              'step_angle_op'  :None}

#################################
# Template for anisotropic case # 
#################################
    
    # hardcoded ejecta geometric and thermal parameters for the aspherical case
    ejecta_params_aniso = {}
    ejecta_params_aniso['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
    ejecta_params_aniso['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
    ejecta_params_aniso['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}

    # set of shell parameters to be sampled on
    ejecta_vars_aniso={}

    ejecta_vars_aniso['dynamics'] = {'xi_disk'        :None,
                                    'm_ej'           :0.04,
                                    'step_angle_mass':None,
                                    'high_lat_flag'  :None,
                                    'central_vel'    :0.33,
                                    'high_lat_vel'   :None,
                                    'low_lat_vel'    :None,
                                    'step_angle_vel' :None,
                                    'central_op'     :None,
                                    'high_lat_op'    :10.,
                                    'low_lat_op'     :30.,
                                    'step_angle_op'  :math.radians(45.)}

    ejecta_vars_aniso['secular'] = {'xi_disk'        :0.2,
                                   'm_ej'           :None,
                                   'step_angle_mass':None,
                                   'high_lat_flag'  :None,
                                   'central_vel'    :0.06,
                                   'high_lat_vel'   :None,
                                   'low_lat_vel'    :None,
                                   'step_angle_vel' :None,
                                   'central_op'     :5.0,
                                   'low_lat_op'     :None,
                                   'high_lat_op'    :None,
                                   'step_angle_op'  :None}
 
    ejecta_vars_aniso['wind'] = {'xi_disk'        :0.1,
                                'm_ej'           :None,
                                'step_angle_mass':math.radians(60.),
                                'high_lat_flag'  :True,
                                'central_vel'    :0.08,
                                'high_lat_vel'   :None,
                                'low_lat_vel'    :None,
                                'step_angle_vel' :None,
                                'central_op'     :None,
                                'high_lat_op'    :0.1,
                                'low_lat_op'     :5.0,
                                'step_angle_op'  :math.radians(30.)}


##########################################################
# choose the appropriate set of parameters and variables #
##########################################################

    if (glob_params['mkn model'] == 'iso1comp'):
        ejecta_params = {}
        ejecta_vars = {}
        ejecta_params['dynamics'] = ejecta_params_iso['dynamics']
        ejecta_vars['dynamics']    = ejecta_vars_iso['dynamics']
    elif (glob_params['mkn model'] == 'iso2comp'):
        ejecta_params = {}
        ejecta_vars = {}
        ejecta_params['dynamics'] = ejecta_params_iso['dynamics']
        ejecta_vars['dynamics']    = ejecta_vars_iso['dynamics']
        ejecta_params['secular'] = ejecta_params_iso['secular']
        ejecta_vars['secular']    = ejecta_vars_iso['secular']
    elif (glob_params['mkn model'] == 'iso3comp'):
        ejecta_params = ejecta_params_iso
        ejecta_vars    = ejecta_vars_iso
    elif (glob_params['mkn model'] == 'aniso1comp'):
        ejecta_params = {}
        ejecta_vars = {}
        ejecta_params['dynamics'] = ejecta_params_aniso['dynamics']
        ejecta_vars['dynamics']    = ejecta_vars_aniso['dynamics']
    elif (glob_params['mkn model'] == 'aniso2comp'):
        ejecta_params = {}
        ejecta_vars = {}
        ejecta_params['dynamics'] = ejecta_params_aniso['dynamics']
        ejecta_vars['dynamics']    = ejecta_vars_aniso['dynamics']
        ejecta_params['secular'] = ejecta_params_aniso['secular']
        ejecta_vars['secular']    = ejecta_vars_aniso['secular']
    elif (glob_params['mkn model'] == 'aniso3comp'):
        ejecta_params = ejecta_params_aniso
        ejecta_vars    = ejecta_vars_aniso

    
    print('I am initializing the model')
    model = MKN(glob_params,glob_vars,ejecta_params,ejecta_vars,source_name)

    print('I am computing the light curves')
#    r_ph,L_bol,T_eff = model.lightcurve(ejecta_vars,glob_params['NR_data'],glob_params['NR_filename']) 
    r_ph,L_bol,T_eff = model.E.lightcurve(model.angular_distribution,
                                          model.omega_distribution,
                                          model.time,
                                          model.ejecta_vars,
                                          model.ejecta_params,
                                          model.glob_vars,
                                          model.glob_params)

    print('I am computing the likelihood')
    logL =  model.log_likelihood(r_ph,T_eff)

    write_output = True
    if (write_output):
        print('I am printing out the output')
        model.write_output(r_ph,T_eff,L_bol)

    plot_separately = False            # Choose to plot all lightcurves in different bands on the same plot
    plot_together = False             # or to plot lightcurve and data in each band on different plots

    if (plot_separately):
        model.plot_output_sep()
    elif (plot_together):
        model.plot_output_tog()


