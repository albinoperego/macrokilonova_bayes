#!/usr/bin/env python
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
        self.n_time = 50
        self.tscale   = 'linear'
        # initialize global time
        if (self.tscale == 'linear'):
            self.time = np.linspace(self.time_min,self.time_max,num=self.n_time)
        elif (tscale == 'log'):
            self.time = np.logspace(np.log10(time_min),np.log10(time_max),num=n_time)
        else:
            print('Error! Wrong option for the time scale')
            exit(-1)
        
        
        # number and distribution of slices for the angular integrals
        self.n_slices = n_slices
        self.dist_slices = dist_slices
        # initialize the angular distribution
        print "Initialising angular distribution"
        self.AD = ad.AngularDistribution(self.dist_slices,self.n_slices)
        self.ang_dist, self.omega_dist = self.AD(self.n_slices/2)

        
        # initialize the filters
        print "Initialising filters"
        self.FT = ft.Filters("measures")
        self.dic_filt,self.lambda_vec,self.mag = self.FT()
        
        # initialise the view angle
        print "Initialising observer projection"
        self.FF = op.ObserverProjection(self.n_slices,self.dist_slices)
        
        # register the geometry of the ejecta and create an Ejecta object
        self.glob_params  = glob_params
        self.glob_vars  = glob_vars
        self.ejecta_params = ejecta_params
        self.shell_params = shell_params
        self.ejecta = ej.Ejecta(len(self.ejecta_params), self.ejecta_params.keys(), self.ejecta_params)
        
        # set of global variables

        model_parameters = self.glob_vars.keys()
        model_parameters.append('distance')
        model_parameters.append('view_angle')
        model_bounds = {'view_angle':[0.0,90.0],
                        'm_disk':[1e-4,1e-2],
                        'eps0':[2e18,2e19],
                        'a_eps_nuc':[0.4,0.8],
                        'b_eps_nuc':[1.0,3.0],
                        't_eps_nuc':[0.5,1.5],
                        'distance':[30,50]}

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

                if type(v) is float or type(v) is np.float64:
                    model_parameters.append(item+'_%s'%s)
                    model_bounds[item+'_%s'%s] = [v-0.75*v,v+0.75*v]

        for n in model_parameters:
            self.names.append(n)
            self.bounds.append(model_bounds[n])
#        print self.names
#        exit()
    def log_likelihood(self,x):
        # populate the relevant parameters
        # ================================
        # set of global variables
        for item in self.glob_vars.keys():
            self.glob_vars[item] = x[item]

#        for shell, p in self.ejecta_params.iteritems():
#            # check for mass distribution choice and add relevant parameter
#            if p['mass_dist'] =="step":
#                self.ejecta_params[shell]['step_angle_mass'] = x['step_angle_mass_%s'%shell]
#
#            # check for velocity distribution choice and add relevant parameters
#            if p['vel_dist'] =="step":
#                self.ejecta_params[shell]['step_angle_vel'] = x['step_angle_vel_%s'%shell]
#                self.ejecta_params[shell]['high_lat_vel'] = x['high_lat_vel_%s'%shell]
#                self.ejecta_params[shell]['low_lat_vel'] = x['low_lat_vel_%s'%shell]
#            elif p['vel_dist'] =="uniform":
#                self.ejecta_params[shell]['central_vel'] = x['central_vel_%s'%shell]
#            else:
#                self.ejecta_params[shell]['min_vel'] = x['min_vel_%s'%shell]
#                self.ejecta_params[shell]['max_vel'] = x['max_vel_%s'%shell]
#
#            # check for opacity distribution choice and add relevant parameters
#            if p['op_dist']=="step":
#                self.ejecta_params[shell]['step_angle_op'] = x['step_angle_op_%s'%shell]
#                self.ejecta_params[shell]['high_lat_op'] = x['high_lat_op_%s'%shell]
#                self.ejecta_params[shell]['low_lat_op'] = x['low_lat_op_%s'%shell]
#
#            elif p['op_dist']=="uniform":
#                self.ejecta_params[shell]['central_op'] = x['central_op_%s'%shell]
#            else:
#                self.ejecta_params[shell]['min_op'] = x['min_op_%s'%shell]
#                self.ejecta_params[shell]['max_op'] = x['max_op_%s'%shell]

        # now add the single shell parameters
        for s in self.shell_params:
            for item in self.shell_params[s]:
                v = self.shell_params[s][item]
                if type(v) is float:
                    self.shell_params[s][item] = x[item+'_%s'%s]
        # ================================

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
        logP = 0.
        if np.isfinite(super(MacroKilonovaModel,self).log_prior(x)):
            logP += -0.5*(x['distance']-40.0)**2/16.0
            return logP
        else:
            return -np.inf

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to test/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=1,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--nlive',default=100,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=10,type='int',metavar='m',help='max MCMC points')
    parser.add_option('--poolsize',default=16,type='int',metavar='k',help='numer of points in the ensemble sampler pool')
    (opts,args)=parser.parse_args()

    # set of global parameters not to be fit
    glob_params = {'v_min':1.e-7,
                   'n_v':100,
                   'vscale':'linear',
                   'sigma0':0.11,
                   'alpha':1.3,
                   't0eps':1.3,
                   'cnst_eff':0.3333}
    
    # set of global parameters to be fit
    glob_vars = {'m_disk':0.3,
                 'eps0':1.2e19,
                 'a_eps_nuc':0.5,
                 'b_eps_nuc':2.5,
                 't_eps_nuc':1.0}
    
    # hardcoded ejecta geometric and thermal parameters
    ejecta_params = {}
    ejecta_params['dynamics']   = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    ejecta_params['wind']       = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    ejecta_params['secular']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}

    # set of shell parameters to be sampled on
    shell_vars={}

    shell_vars['dynamics'] = {'xi_disk':None,
                              'm_ej':0.05,
                              'central_vel':0.2,
                              'low_lat_vel':None,
                              'high_lat_vel':None,
                              'step_angle_vel':None,
                              'low_lat_op':10.,
                              'high_lat_op':0.1,
                              'step_angle_op':np.radians(45.)}

    shell_vars['wind'] = {'xi_disk':0.1,
                          'm_ej':None,
                          'step_angle_mass':np.radians(60.),
                          'high_lat_flag':True,
                          'central_vel':0.05,
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':1.0,
                          'high_lat_op':0.1,
                          'step_angle_op':np.radians(45.)}

    shell_vars['secular'] = {'xi_disk':0.2,
                             'm_ej':None,
                             'central_vel':0.02,
                             'low_lat_vel':None,
                             'high_lat_vel':None,
                             'step_angle_vel':None,
                             'central_op':5.0,
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


