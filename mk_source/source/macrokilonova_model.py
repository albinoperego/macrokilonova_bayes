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

import numpy as np
import math
import cpnest.model
import sys
import os
from optparse import OptionParser
import itertools as it

import angular_distribution as ad
import filters as ft
import lightcurve as lc
import numpy as np
import observer_projection as op
import units
import ejecta as ej

# Necessary to add cwd to path when script run
# by SLURM (since it executes a copy)
sys.path.append(os.getcwd())

class MacroKilonovaModel(cpnest.model.Model):
    """
    Kilonova model
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
                 n_slices = 12,
                 dist_slices = "uniform",
                 **kwargs):

        super(MacroKilonovaModel,self).__init__(**kwargs)
        #initialize the time
        self.time_min = 3600.      #
        self.time_max = 2000000.   #
        self.n_time = 20
        self.tscale   = 'log'
        # initialize global time
        if (self.tscale == 'linear'):
            self.time = np.linspace(self.time_min,self.time_max,num=self.n_time)
        elif (self.tscale == 'log'):
            self.time = np.logspace(np.log10(self.time_min),np.log10(self.time_max),num=self.n_time)
        else:
            print('Error! Wrong option for the time scale')
            exit(-1)
        
        
        # number and distribution of slices for the angular integrals
        self.n_slices = n_slices
        self.dist_slices = dist_slices
        # initialize the angular distribution
        print("Initialising angular distribution")
        self.AD = ad.AngularDistribution(self.dist_slices,self.n_slices)
        self.ang_dist, self.omega_dist = self.AD(self.n_slices/2)

        
        # initialize the filters
        print("Initialising filters")
        self.FT = ft.Filters("measures")
        self.dic_filt,self.lambda_vec,self.mag = self.FT()
        
        # initialise the view angle
        print("Initialising observer projection")
        self.FF = op.ObserverProjection(self.n_slices,self.dist_slices)
        
        # register the geometry of the ejecta and create an Ejecta object
        self.glob_params  = glob_params
        self.glob_vars  = glob_vars
        self.ejecta_params = ejecta_params
        self.shell_params = shell_params
        self.ejecta = ej.Ejecta(len(self.ejecta_params), self.ejecta_params.keys(), self.ejecta_params)
        
        # set of global variables

        model_parameters=['distance', 'view_angle']
        model_bounds = {'view_angle':[0.0,90.0],
                        'distance':[35,45]}
        
        for item in self.glob_vars.keys():
            model_parameters.append(item)
            model_bounds[item] = self.glob_vars[item]

#        for shell, p in self.ejecta_params.iteritems():
#            print "initialising",shell,"component"
#            # check for mass distribution choice and add relevant parameter
#
#            if p['mass_dist'] =="step":
#                model_parameters.append('step_angle_mass_%s'%shell)
#                model_bounds['step_angle_mass_%s'%shell] = [20.0,80.0]
#
#            # check for velocity distribution choice and add relevant parameters
#            if p['vel_dist'] =="step":
#                model_parameters.append('step_angle_vel_%s'%shell)
#                model_bounds['step_angle_vel_%s'%shell] = [20.0,80.0]
#                model_parameters.append('high_lat_vel_%s'%shell)
#                model_bounds['high_lat_vel_%s'%shell] = [1e-1,7e-1]
#                model_parameters.append('low_lat_vel_%s'%shell)
#                model_bounds['low_lat_vel_%s'%shell] = [1e-1,7e-1]
#            elif p['vel_dist'] =="uniform":
#                model_parameters.append('central_vel_%s'%shell)
#                model_bounds['central_vel_%s'%shell] = [1e-1,7e-1]
#            else:
#                model_parameters.append('min_vel_%s'%shell)
#                model_bounds['min_vel_%s'%shell] = [1e-1,7e-1]
#                model_parameters.append('max_vel_%s'%shell)
#                model_bounds['max_vel_%s'%shell] = [1e-1,7e-1] # max_vel_dyn > min_vel_dyn
#
#            # check for opacity distribution choice and add relevant parameters
#            if p['op_dist']=="step":
#                model_parameters.append('step_angle_op_%s'%shell)
#                model_bounds['step_angle_op_%s'%shell] = [20.0,80.0]
#                model_parameters.append('high_lat_op_%s'%shell)
#                model_bounds['high_lat_op_%s'%shell] = [1e-1,1.0]
#                model_parameters.append('low_lat_op_%s'%shell)
#                model_bounds['low_lat_op_%s'%shell] = [1.0,50.0]
#            elif p['op_dist']=="uniform":
#                model_parameters.append('central_op_%s'%shell)
#                model_bounds['central_op_%s'%shell] = [1e-1,50.0]
#            else:
#                model_parameters.append('min_op_%s'%shell)
#                model_bounds['min_op_%s'%shell] = [1e-1,1.0]
#                model_parameters.append('max_op_%s'%shell)
#                model_bounds['max_op_%s'%shell] = [1.0,50.0] # max_op_dyn > min_op_dyn
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
#        print self.names, self.bounds
#        exit()

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
        self.shell_params['dynamics']['xi_disk'] = 1.0 - self.shell_params['wind']['xi_disk'] - self.shell_params['secular']['xi_disk']
        x['xi_disk_dynamics'] = self.shell_params['dynamics']['xi_disk']
        
    def log_likelihood(self,x):

        self.flux_factor = self.FF(x['view_angle'])

        r_ph, L_bol, T_eff = self.ejecta.lightcurve(self.ang_dist, self.omega_dist,
                                                  self.time,
                                                  self.shell_params,
                                                  self.glob_vars,
                                                  self.glob_params)

        # compute the magnitudes from a certain distance
        D = x['distance']*1e6*units.pc2cm
        model_mag = ft.calc_magnitudes(self.flux_factor,self.time,r_ph,T_eff,self.lambda_vec,self.dic_filt,D)

        # compute the residuals
        residuals = ft.calc_residuals(self.mag,model_mag)

        # compute the likelihood
        logL = 0.
        for ilambda in residuals.keys():
            logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

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
    (opts,args)=parser.parse_args()
    np.seterr(all='ignore')
    # set of global parameters not to be fit
    glob_params = {'v_min':1.e-7,
                   'n_v':10,
                   'vscale':'linear',
                   'sigma0':0.11,
                   'alpha':1.3,
                   't0eps':1.3,
                   'cnst_eff':0.3333}
    
    # set of global parameters to be fit
    glob_vars = {'m_disk':[0.0001, 0.5],
                 'eps0':  [2.e17, 2.5e20],
                 'a_eps_nuc':[0.495, 0.505],
                 'b_eps_nuc':[2.45, 2.55],
                 't_eps_nuc':[0.95, 1.05]}
    
    # hardcoded ejecta geometric and thermal parameters
    ejecta_params = {}
    ejecta_params['dynamics']   = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    ejecta_params['wind']       = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    ejecta_params['secular']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}

    # set of shell parameters to be sampled on
    shell_vars={}

    shell_vars['dynamics'] = {'xi_disk':[0.0, 1.0],
                              'm_ej':None,
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
        model = MacroKilonovaModel(glob_params, glob_vars, ejecta_params, shell_vars)
        work  = cpnest.CPNest(model,
                            verbose=2,
                            Poolsize=opts.poolsize,
                            Nthreads=opts.threads,
                            Nlive=opts.nlive,
                            maxmcmc=opts.maxmcmc,
                            output=opts.out_dir)
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

