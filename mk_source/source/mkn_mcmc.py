import angular_distribution as ad
import corner
import emcee
#from emcee.utils import MPIPool
import ejecta as ej
import filters as ft
import math
import matplotlib.pyplot as plt
import numpy as np
import observer_projection as op
from optparse import OptionParser
import os
import sys
import source_properties as sp
import units


class MKN(object):
    """
    Driver for the kilonova model
    """

    def __init__(self,
                 glob_params,
                 ejecta_params,
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
        self.dic_filt,self.lambda_vec,self.mag = self.FT(SP.filter_data_folder,glob_params["time min"]*units.sec2day+SP.t0,glob_params["time max"]*units.sec2day+SP.t0)


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

        # register the geometry of the ejecta and create an Ejecta object
        print('I am initializing the components')
        self.E = ej.Ejecta(len(ejecta_params.keys()), ejecta_params.keys(), ejecta_params)

    def compute_log_likelihood(self,residuals):
        # compute the likelihood
        logL = 0.
        for ilambda in residuals.keys():
            logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))
        return logL

    def log_likelihood(self,r_ph,T_eff):

        # check the viewing angle and overwrite, if necessary
        if (glob_vars['view_angle'] != None):
            self.view_angle=glob_vars['view_angle']

        # check the source distance and overwrite, if necessary
        if (glob_vars['source_distance'] != None):
            self.D = glob_vars['source_distance']*units.pc2cm

        self.flux_factor = self.FF(self.view_angle)

        # compute the residuals
        if (source_name != 'default'):
            residuals = ft.calc_all_residuals(self.flux_factor,self.time,r_ph,T_eff,self.lambda_vec,self.dic_filt,self.D,self.t0,self.mag)

        # compute the likelihood
            logL = self.compute_log_likelihood(residuals)

        else:
            logL = 0.

        return logL

    def log_prior(self,theta,vmin,vmax):
        for i in range(len(theta)):
            if ( (theta[i] < vmin[i]) or (theta[i] > vmax[i])):
                return -np.inf
        return 0.0

    def log_like(self,theta,list_var,ejecta_vars,ejecta_params,glob_vars,glob_params):

        for val,item in zip(theta,list_var):
            if (item == 'gv_m_disk'):
                glob_vars['m_disk'] = val
            elif (item == 'gv_eps0'):
                glob_vars['eps0'] = val 
            elif (item == 'gv_view_angle'):
                glob_vars['view_angle'] = val 
            elif (item == 'gv_source_distance'):
                glob_vars['source_distance'] = val 
            elif (item == 'gv_T_floor_LA'):
                glob_vars['T_floor_LA'] = val
            elif (item == 'gv_T_floor_Ni'):
                glob_vars['T_floor_Ni'] = val
            elif (item == 'gv_a_eps_nuc'):
                glob_vars['a_eps_nuc'] = val
            elif (item == 'gv_b_eps_nuc'):
                glob_vars['b_eps_nuc'] = val
            elif (item == 'gv_t_eps_nuc'):
                glob_vars['t_eps_nuc'] = val
            elif (item == 'ev_dyn_xi_disk'):
                ejecta_vars['dynamics']['xi_disk'] = val
            elif (item == 'ev_dyn_m_ej'):
                ejecta_vars['dynamics']['m_ej'] = val
            elif (item == 'ev_dyn_step_angle_mass'):
                ejecta_vars['dynamics']['m_step_angle_mass'] = val
            elif (item == 'ev_dyn_central_vel'):
                ejecta_vars['dynamics']['central_vel'] = val
            elif (item == 'ev_dyn_high_lat_vel'):
                ejecta_vars['dynamics']['high_lat_vel'] = val
            elif (item == 'ev_dyn_low_lat_vel'):
                ejecta_vars['dynamics']['low_lat_vel'] = val
            elif (item == 'ev_dyn_step_angle_vel'):
                ejecta_vars['dynamics']['step_angle_vel'] = val
            elif (item == 'ev_dyn_central_op'):
                ejecta_vars['dynamics']['central_op'] = val
            elif (item == 'ev_dyn_high_lat_op'):
                ejecta_vars['dynamics']['high_lat_op'] = val
            elif (item == 'ev_dyn_low_lat_op'):
                ejecta_vars['dynamics']['low_lat_op'] = val
            elif (item == 'ev_dyn_step_angle_op'):
                ejecta_vars['dynamics']['step_angle_op'] = val
            elif (item == 'ev_sec_xi_disk'):
                ejecta_vars['secular']['xi_disk'] = val
            elif (item == 'ev_sec_m_ej'):
                ejecta_vars['secular']['m_ej'] = val
            elif (item == 'ev_sec_step_angle_mass'):
                ejecta_vars['secular']['m_step_angle_mass'] = val
            elif (item == 'ev_sec_central_vel'):
                ejecta_vars['secular']['central_vel'] = val
            elif (item == 'ev_sec_high_lat_vel'):
                ejecta_vars['secular']['high_lat_vel'] = val
            elif (item == 'ev_sec_low_lat_vel'):
                ejecta_vars['secular']['low_lat_vel'] = val
            elif (item == 'ev_sec_step_angle_vel'):
                ejecta_vars['secular']['step_angle_vel'] = val
            elif (item == 'ev_sec_central_op'):
                ejecta_vars['secular']['central_op'] = val
            elif (item == 'ev_sec_high_lat_op'):
                ejecta_vars['secular']['high_lat_op'] = val
            elif (item == 'ev_sec_low_lat_op'):
                ejecta_vars['secular']['low_lat_op'] = val
            elif (item == 'ev_sec_step_angle_op'):
                ejecta_vars['secular']['step_angle_op'] = val
            elif (item == 'ev_wind_xi_disk'):
                ejecta_vars['wind']['xi_disk'] = val
            elif (item == 'ev_wind_m_ej'):
                ejecta_vars['wind']['m_ej'] = val
            elif (item == 'ev_wind_step_angle_mass'):
                ejecta_vars['wind']['m_step_angle_mass'] = val
            elif (item == 'ev_wind_central_vel'):
                ejecta_vars['wind']['central_vel'] = val
            elif (item == 'ev_wind_high_lat_vel'):
                ejecta_vars['wind']['high_lat_vel'] = val
            elif (item == 'ev_wind_low_lat_vel'):
                ejecta_vars['wind']['low_lat_vel'] = val
            elif (item == 'ev_wind_step_angle_vel'):
                ejecta_vars['wind']['step_angle_vel'] = val
            elif (item == 'ev_wind_central_op'):
                ejecta_vars['wind']['central_op'] = val
            elif (item == 'ev_wind_high_lat_op'):
                ejecta_vars['wind']['high_lat_op'] = val
            elif (item == 'ev_wind_low_lat_op'):
                ejecta_vars['wind']['low_lat_op'] = val
            elif (item == 'ev_wind_step_angle_op'):
                ejecta_vars['wind']['step_angle_op'] = val
            else:
                print('Unknown variable for the mcmc! please check...')
                exit()
 

        #print('I am computing the light curves')
        r_ph,L_bol,T_eff = self.E.lightcurve(self.angular_distribution,
                                              self.omega_distribution,
                                              self.time,
                                              ejecta_vars,
                                              ejecta_params,
                                              glob_vars,
                                              glob_params)

        #print('I am computing the likelihood')
        logL =  model.log_likelihood(r_ph,T_eff)

        return logL

# ===================================================================================

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('--outdir',default=None,type='string',metavar='DIR',help='Directory for output')
    parser.add_option('--nthreads',default=None,type='int',metavar='N',help='Number of threads (default = 1)')
    parser.add_option('--nwalkers',default=20,type='int',metavar='n',help='Number of walkers (default = 100)')
    parser.add_option('--nsteps',default=1000,type='int',metavar='m',help='Number of steps (default = 1000)')
    parser.add_option('--nburn',default=100,type='int',metavar='k',help='Number of burn-in steps (default = 100)')
    (opts,args)=parser.parse_args()

#dictionary with the global parameters of the model
    glob_params = {'lc model'   :'grossman',    # model for the lightcurve (grossman or villar)  
                   'mkn model'  :'iso1comp',    # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                   'omega frac' :1.0,           #
                   'rad shell'  :True,          #
                   'v_min'      :1.e-7,         # minimal velocity for the Grossman model
                   'n_v'        :400,           # number of points for the Grossman model
                   'vscale'     :'linear',      # scale for the velocity in the Grossman model
                   'sigma0'     :0.11,          # parameter for the nuclear heating rate
                   'alpha'      :1.3,           # parameter for the nuclear heating rate
                   't0eps'      :1.3,           # parameter for the nuclear heating rate
                   'cnst_eff'   :0.3333,        # parameter for the constant heating efficiency
                   'n slices'   :12,            # number for the number of slices along the polar angle [12,18,24,30]
                   'dist slices':'cos_uniform', # discretization law for the polar angle [uniform or cos_uniform]
                   'time min'   :3600.,         # minimum time [s]
                   'time max'   :432000.,       # maximum time [s]
                   'n time'     :200,           # integer number of bins in time
                   'scale for t':'measures',    # kind of spacing in time [log - linear - measures]
                   'NR_data'    :False,         # use (True) or not use (False) NR profiles
                   'NR_filename':'../example_NR_data/DD2_M125125_LK/outflow_1/ejecta_profile.dat'           # path of the NR profiles, necessary if NR_data is True
                   }
    
    source_name = 'AT2017gfo'   # name of the source or "default"
#    source_name = 'default'   # name of the source or "default"
    
    # dictionary for the global variables
    glob_vars = {'m_disk':0.12,
                 'eps0':1.5e19,
                 'view_angle':None,      # [deg]; it uses the one in source properties
                 'source_distance':None, # [pc] ; it uses the one in source properties
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
    ejecta_params_iso['dynamics'] = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}
    ejecta_params_iso['wind']     = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}
    ejecta_params_iso['secular']  = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}

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
    ejecta_params_aniso['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}
    ejecta_params_aniso['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}
    ejecta_params_aniso['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True,'v_law':'poly'}

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
    model = MKN(glob_params,ejecta_params,source_name)

    def log_prob(theta,var_list,vmin,vmax,ejecta_vars,ejecta_params,glob_vars,glob_params):
        lp = model.log_prior(theta,vmin,vmax)
        if not np.isfinite(lp):
            return -np.inf
        return lp + model.log_like(theta,var_list,ejecta_vars,ejecta_params,glob_vars,glob_params)

    ndim     = 6
    nwalkers = opts.nwalkers

    mcmc_var_list=[
                   #'gv_m_disk',
                   'gv_eps0',
                   #'gv_view_angle',
                   #'gv_source_distance',
                   'gv_T_floor_LA',
                   'gv_T_floor_Ni',
                   #'gv_a_eps_nuc',
                   #'gv_b_eps_nuc',
                   #'gv_t_eps_nuc',
                   #'ev_dyn_xi_disk', 
                   'ev_dyn_m_ej', 
                   #'ev_dyn_step_angle_mass', 
                   'ev_dyn_central_vel', 
                   #'ev_dyn_high_lat_vel', 
                   #'ev_dyn_low_lat_vel', 
                   #'ev_dyn_step_angle_vel', 
                   'ev_dyn_central_op', 
                   #'ev_dyn_high_lat_op', 
                   #'ev_dyn_low_lat_op', 
                   #'ev_dyn_step_angle_op', 
                   #'ev_sec_xi_disk', 
                   #'ev_sec_m_ej', 
                   #'ev_sec_step_angle_mass', 
                   #'ev_sec_central_vel', 
                   #'ev_sec_high_lat_vel', 
                   #'ev_sec_low_lat_vel', 
                   #'ev_sec_step_angle_vel', 
                   #'ev_sec_central_op', 
                   #'ev_sec_high_lat_op', 
                   #'ev_sec_low_lat_op', 
                   #'ev_sec_step_angle_op', 
                   #'ev_wind_xi_disk', 
                   #'ev_wind_m_ej', 
                   #'ev_wind_step_angle_mass', 
                   #'ev_wind_central_vel', 
                   #'ev_wind_high_lat_vel', 
                   #'ev_wind_low_lat_vel', 
                   #'ev_wind_step_angle_vel', 
                   #'ev_wind_central_op', 
                   #'ev_wind_high_lat_op', 
                   #'ev_wind_low_lat_op', 
                   #'ev_wind_step_angle_op', 
                   ]

    mcmc_var_minmax=[
                     #[],   #'gv_m_disk'
                     [1.e18,1.e20],   #'gv_eps0'
                     #[]    # gv_view_angle',
                     #[]    # gv_source_distance',
                     [500.,5000.],   #'gv_T_floor_LA'
                     [500.,5000.],   #'gv_T_floor_Ni'
                     #[],   #'gv_a_eps_nuc'
                     #[],   #'gv_b_eps_nuc'
                     #[],   #'gv_t_eps_nuc'
                     #[],   #'ev_dyn_xi_disk', 
                     [0.001,0.2],   #'ev_dyn_m_ej', 
                     #[],   #'ev_dyn_step_angle_mass', 
                     [0.01,0.33],   #'ev_dyn_central_vel', 
                     #[],   #'ev_dyn_high_lat_vel', 
                     #[],   #'ev_dyn_low_lat_vel', 
                     #[],   #'ev_dyn_step_angle_vel', 
                     [0.1,50.],   #'ev_dyn_central_op', 
                     #[],   #'ev_dyn_high_lat_op', 
                     #[],   #'ev_dyn_low_lat_op', 
                     #[],   #'ev_dyn_step_angle_op', 
                     #[],   #'ev_sec_xi_disk', 
                     #[],   #'ev_sec_m_ej', 
                     #[],   #'ev_sec_step_angle_mass', 
                     #[],   #'ev_sec_central_vel', 
                     #[],   #'ev_sec_high_lat_vel', 
                     #[],   #'ev_sec_low_lat_vel', 
                     #[],   #'ev_sec_step_angle_vel', 
                     #[],   #'ev_sec_central_op', 
                     #[],   #'ev_sec_high_lat_op', 
                     #[],   #'ev_sec_low_lat_op', 
                     #[],   #'ev_sec_step_angle_op', 
                     #[],   #'ev_wind_xi_disk', 
                     #[],   #'ev_wind_m_ej', 
                     #[],   #'ev_wind_step_angle_mass', 
                     #[],   #'ev_wind_central_vel', 
                     #[],   #'ev_wind_high_lat_vel', 
                     #[],   #'ev_wind_low_lat_vel', 
                     #[],   #'ev_wind_step_angle_vel', 
                     #[],   #'ev_wind_central_op', 
                     #[],   #'ev_wind_high_lat_op', 
                     #[],   #'ev_wind_low_lat_op', 
                     #[],   #'ev_wind_step_angle_op', 
                   ]

    vmin = np.asarray([item[0] for item in mcmc_var_minmax])
    vmax = np.asarray([item[1] for item in mcmc_var_minmax])

    # flat initialization of the walkers
    p0 = [vmin + np.random.uniform(0.,1.,ndim)*(vmax-vmin) for i in range(nwalkers)]

#    pool = MPIPool()
#    if not pool.is_master():
#        pool.wait()
#        sys.exit(0)

    # initialize the sampler
    sampler = emcee.EnsembleSampler(opts.nwalkers, 
                                    ndim, 
                                    log_prob, 
                                    args=(mcmc_var_list,
                                          vmin,
                                          vmax,
                                          ejecta_vars,
                                          ejecta_params,
                                          glob_vars,
                                          glob_params), 
                                    threads=opts.nthreads)

    # burn-in to better initialize the walkers
    print('Starting the burn-in phase')
    pos,prob,state = sampler.run_mcmc(p0,opts.nburn)
    print('Finished the burn-in phase')

    if (opts.outdir is not 'None'):
        try:
            os.stat(opts.outdir)
        except:
            os.mkdir(opts.outdir)

        f = open(opts.outdir+"/chain.dat", "w")
        f.close()

    for i, result in enumerate(sampler.sample(pos, iterations=opts.nsteps)):
        if (i+1) % 100 == 0:

            # save an intermediate step
            print('Step %i/%i completed!' %(i+1,opts.nsteps))
            np.save(opts.outdir+"/last_intermediate_step.npy", sampler.chain)

            if (opts.outdir is not 'None'):
                f = open(opts.outdir+"/chain.dat", "a")
                position = result[0]
                for k in range(position.shape[0]):
                #f.write("{0:4d} {1:s}\n".format(k, " ".join(str(position[k]))))
                    f.write("%i %30s \n" %(k,position[k]))
                f.close()

            # position = result[0]
            #f = open("chain.dat", "a")
            #for k in range(position.shape[0]):
            #    #f.write("{0:4d} {1:s}\n".format(k, " ".join(str(position[k]))))
            #    f.write("{0:4d} {1:s}\n".format(k,position[k]))
            #f.close()

    # save final results
    print('All steps completed!')
    np.save(opts.outdir+"/last_step.npy", sampler.chain)

    exit()

    samples = sampler.chain[:, 10:, :].reshape((-1, ndim))

    eps0_mcmc, T_ni_mcmc, T_la_mcmc, m_mcmc, v_mcmc, k_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],axis=0)))

    print('T_ni mcmc',T_ni_mcmc)
    print('T_la mcmc',T_la_mcmc)
    print('eps0 mcmc',eps0_mcmc)
    print('m mcmc',m_mcmc)
    print('v mcmc',v_mcmc)
    print('k mcmc',k_mcmc)
   
#    '''
#    fig2 = plt.figure(2)
#    for i in range(nwalkers):
#        plt.plot(sampler.chain[i,:,0])

#    fig3 = plt.figure(3)
#    for i in range(nwalkers):
#        plt.plot(sampler.chain[i,:,1])
#
#    fig4 = plt.figure(4)
#    for i in range(nwalkers):
#        plt.plot(sampler.chain[i,:,2])
#    '''

    if (opts.outdir is not 'None'):
        fig = corner.corner(samples, labels=["$m$", "$v$", "$k$"])
        fig.savefig(opts.outdir+"/triangle.png")

#    pool.close()
