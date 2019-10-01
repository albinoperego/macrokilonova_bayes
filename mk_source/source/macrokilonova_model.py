#!/usr/bin/env python3.5
## usage: sbatch <script>
#SBATCH --account=INF18_teongrav
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=68
#SBATCH --partition=knl_usr_prod
###SBATCH --partition=knl_usr_dbg
#SBATCH --job-name=myjob
#SBATCH --mail-user=walter.delpozzo@unipi.it
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00

import cpnest.model
import itertools as it
import math
import numpy as np
from optparse import OptionParser
import os
import sys

import angular_distribution as ad
import ejecta as ej
import filters as ft
import mkn
import numpy as np
import observer_projection as op
import source_properties as sp
import units

# Necessary to add cwd to path when script run
# by SLURM (since it executes a copy)
sys.path.append(os.getcwd())

class MacroKilonovaModel(cpnest.model.Model):
    """
    Kilonova model for cpnest
    """

    names = []
    bounds = []

    def __init__(self,
                 # number of different components in the ejecta
                 glob_params,
                 # dictionary of global parameters defining basic properties of the ejecta to be sampled
                 glob_vars,
                 # dictionary of ejecta parameters defining its composition and geometry
                 ejecta_params,
                 # dictionary of shell parameters defining basic properties of the shell
                 ejecta_vars,
                 source_name,
                 **kwargs):

        self.MKN = mkn.MKN(glob_params,
                           glob_vars,
                           ejecta_params,
                           ejecta_vars,
                           source_name)

        self.glob_params   = glob_params
        self.glob_vars     = glob_vars
        self.ejecta_params = ejecta_params
        self.ejecta_vars    = ejecta_vars

        # set of global variables
        model_parameters=['distance', 'view_angle']
        model_bounds = {'view_angle':[0.0,90.0],
                        'distance':[39.9,40.1]}
        
        for item in self.glob_vars.keys():
            model_parameters.append(item)
            model_bounds[item] = self.glob_vars[item]

        # now add the single shell parameters
        for s in self.ejecta_vars:
            for item in self.ejecta_vars[s]:
                v = self.ejecta_vars[s][item]
                if v is not None and type(v) is not bool:
                    model_parameters.append(item+'_%s'%s)
                    model_bounds[item+'_%s'%s] = v

        for n in model_parameters:
            self.names.append(n)
            self.bounds.append(model_bounds[n])

    def fill_control_structures(self, x):
        # populate the relevant parameters
        # set of global variables
        for item in self.glob_vars.keys():
            self.glob_vars[item] = x[item]

        # now add the single shell parameters
        for s in self.ejecta_vars:
            for item in self.ejecta_vars[s]:
                v = self.ejecta_vars[s][item]
                if v is not None and type(v) is not bool:
                    self.ejecta_vars[s][item] = x[item+'_%s'%s]
        # impose the mass constraint
#        self.ejecta_vars['secular']['xi_disk'] = min(self.ejecta_vars['secular']['xi_disk'], 1.0 - self.ejecta_vars['wind']['xi_disk'])
#        x['xi_disk_secular'] = self.ejecta_vars['secular']['xi_disk']

    # define force for hamiltonian montecarlo 
    def force(self, x):
        return np.zeros(1, dtype = {'names':x.names, 'formats':['f8' for _ in x.names]})

    def log_likelihood(self,x):

        self.flux_factor = self.MKN.FF(x['view_angle'])

        r_ph, L_bol, T_eff = self.MKN.E.lightcurve(self.MKN.angular_distribution, 
                                                   self.MKN.omega_distribution,
                                                   self.MKN.time,
                                                   self.MKN.ejecta_vars,
                                                   self.MKN.ejecta_params,
                                                   self.MKN.glob_vars,
                                                   self.MKN.glob_params)

        # compute the residuals
        D = x['distance']*1e6*units.pc2cm
        residuals = ft.calc_all_residuals(self.flux_factor,self.MKN.time,r_ph,T_eff,self.MKN.lambda_vec,self.MKN.dic_filt,D,self.MKN.t0,self.MKN.mag)
        logL = self.MKN.compute_log_likelihood(residuals)

        return logL
            
    def log_prior(self, x):
        logP = 0.0
        if np.isfinite(super(MacroKilonovaModel,self).log_prior(x)):
            self.fill_control_structures(x)
            logP += -0.5*(x['distance']-40.0)**2/16.0
            return logP
        else:
            return -np.inf

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to test/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=1,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--nlive',default=1000,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=1000,type='int',metavar='m',help='max MCMC points')
    parser.add_option('--poolsize',default=1000,type='int',metavar='k',help='number of points in the ensemble sampler pool')
    parser.add_option('--seed',default=1234,type='int',metavar='s',help='seed for the pseudo-random chain')
    (opts,args)=parser.parse_args()
    np.seterr(all='ignore')

    source_name = 'AT2017gfo'

    #dictionary with the global parameters of the model not to be fit
    glob_params = {'lc model'   :'grossman',    # model for the lightcurve (grossman or villar)
                   'mkn model'  :'iso1comp',    # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                   'omega frac' :1.0,           #
                   'rad shell'  :True,          #
                   'v_min'      :1.e-7,         # minimal velocity for the Grossman model
                   'n_v'        :200,           # number of points for the Grossman model
                   'vscale'     :'linear',      # scale for the velocity in the Grossman model
                   'sigma0'     :0.11,          # parameter for the nuclear heating rate
                   'alpha'      :1.3,           # parameter for the nuclear heating rate
                   't0eps'      :1.3,           # parameter for the nuclear heating rate
                   'cnst_eff'   :0.3333,        # parameter for the constant heating efficiency
                   'n slices'   :30,            # number for the number of slices along the polar angle [12,18,24,30]
                   'dist slices':'cos_uniform', # discretization law for the polar angle [uniform or cos_uniform]
                   'time min'   :3600.,         # minimum time [s]
                   'time max'   :2000000.,      # maximum time [s]
                   'n time'     :20,            # integer number of bins in time
                   'scale for t':'measures',    # kind of spacing in time [log - linear - measures]
                   'NR_data'    :False,         # use (True) or not use (False) NR profiles
                   'NR_filename':'../example_NR_data/DD2_M125125_LK/outflow_1/ejecta_profile.dat'
                       # path of the NR profiles, necessary if NR_data is True
                  }

    # set of global parameters to be fit
    glob_vars = {'m_disk':[0.0001, 0.5],
                 'eps0':  [2.e17, 2.5e20],
                 'T_floor_LA':[100.,3000.],
                 'T_floor_Ni':[3000.,8000.],
                 'a_eps_nuc':[0.499, 0.501], # to fix
                 'b_eps_nuc':[2.49, 2.51],   # to fix
                 't_eps_nuc':[0.99, 1.01]}   # to fix

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
                                  'm_ej'           :[5.e-4,1.e-1],
                                  'step_angle_mass':None,
                                  'high_lat_flag'  :None,
                                  'central_vel'    :[0.001, 0.333],
                                  'high_lat_vel'   :None,
                                  'low_lat_vel'    :None,
                                  'step_angle_vel' :None,
                                  'central_op'     :[0.1,40.],
                                  'high_lat_op'    :None,
                                  'low_lat_op'     :None,
                                  'step_angle_op'  :None}

    ejecta_vars_iso['secular'] = {'xi_disk'        :None,
                                 'm_ej'           :[5.e-4,2.e-2],
                                 'step_angle_mass':None,
                                 'high_lat_flag'  :None,
                                 'central_vel'    :[0.001, 0.333],
                                 'high_lat_vel'   :None,
                                 'low_lat_vel'    :None,
                                 'step_angle_vel' :None,
                                 'central_op'     :[0.1,20.0],
                                 'low_lat_op'     :None,
                                 'high_lat_op'    :None,
                                 'step_angle_op'  :None}

    ejecta_vars_iso['wind'] = {'xi_disk'        :None,
                              'm_ej'           :[5.e-4,2.e-2],
                              'step_angle_mass':None,
                              'high_lat_flag'  :True,
                              'central_vel'    :[0.001, 0.333],
                              'high_lat_vel'   :None,
                              'low_lat_vel'    :None,
                              'step_angle_vel' :None,
                              'central_op'     :[0.1,30.],
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
                                    'm_ej'           :[5.e-4,2.e-2],
                                    'step_angle_mass':None,
                                    'high_lat_flag'  :None,
                                    'central_vel'    :[0.001, 0.333],
                                    'high_lat_vel'   :None,
                                    'low_lat_vel'    :None,
                                    'step_angle_vel' :None,
                                    'central_op'     :None,
                                    'high_lat_op'    :[0.1,2.0],
                                    'low_lat_op'     :[1.0,20.0],
                                    'step_angle_op'  :[0.0,np.pi/2.0]}

    ejecta_vars_aniso['secular'] = {'xi_disk'        :[0.0,1.0],
                                   'm_ej'           :None,
                                   'step_angle_mass':None,
                                   'high_lat_flag'  :None,
                                   'central_vel'    :[0.001, 0.333],
                                   'high_lat_vel'   :None,
                                   'low_lat_vel'    :None,
                                   'step_angle_vel' :None,
                                   'central_op'     :[0.01,20.0],
                                   'low_lat_op'     :None,
                                   'high_lat_op'    :None,
                                   'step_angle_op'  :None}
 
    ejecta_vars_aniso['wind'] = {'xi_disk'        :[0.0,1.0],
                                'm_ej'           :None,
                                'step_angle_mass':[0.0,np.pi/2.0],
                                'high_lat_flag'  :True,
                                'central_vel'    :[0.001, 0.333],
                                'high_lat_vel'   :None,
                                'low_lat_vel'    :None,
                                'step_angle_vel' :None,
                                'central_op'     :None,
                                'high_lat_op'    :[0.01,1.0],
                                'low_lat_op'     :[1.0,20.0],
                                'step_angle_op'  :[0.0,np.pi/2.0]}


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



    if opts.out_dir is None:
        opts.out_dir='./test/'

    if opts.full_run:
        model = MacroKilonovaModel(glob_params, glob_vars, ejecta_params, ejecta_vars,source_name)
        work  = cpnest.CPNest(model,
                              verbose=2,
                              poolsize=opts.poolsize,
                              nthreads=opts.threads,
                              nlive=opts.nlive,
                              maxmcmc=opts.maxmcmc,
                              output=opts.out_dir,
                              seed=opts.seed)
        work.run()
        work.get_posterior_samples(filename = 'posterior.dat')
"""
    fig1 = plt.figure()
    band_list = ['U','B','g','V','R']
    for band in band_list:
    for ilambda in model_mag.keys():
    try:
    if(dic_filt[ilambda]['name']==band):
    plt.plot(time/24./60./60.,model_mag[ilambda],color=dic_filt[ilambda]['color'],label=dic_filt[ilambda]['name'])
    plt.scatter(mag[ilambda]['time']-57982.529,mag[ilambda]['mag'],color=dic_filt[ilambda]['color'])
    except KeyError:
    continue
    plt.title('UV and visible bands')
    plt.xlabel('Time [day]')
    plt.ylabel('AB magnitude')
    plt.legend(loc='0')
    plt.xlim(0.1,10)
    plt.ylim(27,15)
"""

