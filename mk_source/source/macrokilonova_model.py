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

        self.MKN = mkn.MKN(glob_params,
                  glob_vars,
                  ejecta_params,
                  shell_params,
                  source_name)

        self.glob_params   = glob_params
        self.glob_vars     = glob_vars
        self.ejecta_params = ejecta_params
        self.shell_params  = shell_params

        # set of global variables
        model_parameters=['distance', 'view_angle']
        model_bounds = {'view_angle':[0.0,90.0],
                        'distance':[39.9,40.1]}
        
        for item in self.glob_vars.keys():
            model_parameters.append(item)
            model_bounds[item] = self.glob_vars[item]

        # now add the single shell parameters
        for s in self.shell_params:
            for item in self.shell_params[s]:
                v = self.shell_params[s][item]
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
        for s in self.shell_params:
            for item in self.shell_params[s]:
                v = self.shell_params[s][item]
                if v is not None and type(v) is not bool:
                    self.shell_params[s][item] = x[item+'_%s'%s]
        # impose the mass constraint
        self.shell_params['secular']['xi_disk'] = min(self.shell_params['secular']['xi_disk'], 1.0 - self.shell_params['wind']['xi_disk'])
        x['xi_disk_secular'] = self.shell_params['secular']['xi_disk']

    def log_likelihood(self,x):

        self.flux_factor = self.MKN.FF(x['view_angle'])

        r_ph, L_bol, T_eff = self.MKN.ejecta.lightcurve(self.MKN.angular_distribution, self.MKN.omega_distribution,
                                                  self.MKN.time,
                                                  self.MKN.shell_params,
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
    glob_params = {'lc model'   :'villar',      # model for the lightcurve (grossman or villar)  
                   'v_min'      :1.e-7,         # minimal velocity for the Grossman model
                   'n_v'        :10,            # number of points for the Grossman model
                   'vscale'     :'linear',      # scale for the velocity in the Grossman model
                   'sigma0'     :0.11,          # parameter for the nuclear heating rate
                   'alpha'      :1.3,           # parameter for the nuclear heating rate
                   't0eps'      :1.3,           # parameter for the nuclear heating rate
                   'cnst_eff'   :0.3333,        # parameter for the constant heating efficiency
                   'n slices'   :12,            # number for the number of slices along the polar angle [12,18,24,30]
                   'dist slices':'uniform',     # discretization law for the polar angle [uniform or cos_uniform]
                   'time min'   :3600.,         # minimum time [s]
                   'time max'   :2000000.,      # maximum time [s]
                   'n time'     :20,            # integer number of bins in time
                   'scale for t':'linear'       # kind of spacing in time [log - linear - measures]
                  }

    # set of global parameters to be fit
    glob_vars = {'m_disk':[0.0001, 0.5],
                 'eps0':  [2.e17, 2.5e20],
                 'T_floor_LA':[100.,3000.],
                 'T_floor_Ni':[3000.,8000.],
                 'a_eps_nuc':[0.499, 0.501],
                 'b_eps_nuc':[2.49, 2.51],
                 't_eps_nuc':[0.99, 1.01]}
    
    # hardcoded ejecta geometric and thermal parameters
    ejecta_params = {}
    ejecta_params['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True}
    ejecta_params['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   ,'therm_model':'BKWM','eps_ye_dep':True}
    ejecta_params['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform','therm_model':'BKWM','eps_ye_dep':True}

    # set of shell parameters to be sampled on
    shell_vars={}

    shell_vars['dynamics'] = {'xi_disk':None,
                              'm_ej':[5.e-4,2.e-2],
                              'central_vel':[0.001, 0.333],
                              'low_lat_vel':None,
                              'high_lat_vel':None,
                              'step_angle_vel':None,
                              'low_lat_op':[1.0,20.0],
                              'high_lat_op':[0.1,2.0],
                              'step_angle_op':[0.0,np.pi/2.0]}

    shell_vars['wind'] = {'xi_disk':[0.0,1.0],
                          'm_ej':None,
                          'step_angle_mass':[0.0,np.pi/2.0],
                          'high_lat_flag':True,
                          'central_vel':[0.001, 0.333],
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':[1.0,20.0],
                          'high_lat_op':[0.01,1.0],
                          'step_angle_op':[0.0,np.pi/2.0]}

    shell_vars['secular'] = {'xi_disk':[0.0,1.0],
                             'm_ej':None,
                             'central_vel':[0.001, 0.333],
                             'low_lat_vel':None,
                             'high_lat_vel':None,
                             'step_angle_vel':None,
                             'central_op':[0.01,20.0],
                             'low_lat_op':None,
                             'high_lat_op':None,
                             'step_angle_op':None}

    if opts.out_dir is None:
        opts.out_dir='./test/'

    if opts.full_run:
        model = MacroKilonovaModel(glob_params, glob_vars, ejecta_params, shell_vars,source_name)
        work  = cpnest.CPNest(model,
                            verbose=2,
                            Poolsize=opts.poolsize,
                            Nthreads=opts.threads,
                            Nlive=opts.nlive,
                            maxmcmc=opts.maxmcmc,
                            output=opts.out_dir,
                            seed=opts.seed)
        work.run()

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

