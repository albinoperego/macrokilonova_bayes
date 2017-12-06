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

model_parameters = ['m_ej_dyn',
                    'central_vel_dyn',
                    'central_op_dyn',
                    'm_ej_wind',
                    'central_vel_wind',
                    'central_op_wind',
                    'm_ej_sec',
                    'central_vel_sec',
                    'central_op_sec']

model_bounds = {'m_ej_dyn':[0.0,0.1],
                'central_op_dyn':[0.0,1.0],
                'central_vel_dyn':[0.0,1.0],
                'm_ej_wind':[0.0,0.1],
                'central_vel_wind':[0.0,1.0],
                'central_op_wind':[0.0,0.1],
                'm_ej_sec':[0.0,0.1],
                'central_vel_sec':[0.0,1.0],
                'central_op_sec':[0.0,1.0]}


class MacroKilonovaModel(cpnest.model.Model):
    """
    Kilonova model
    """

    names = []
    bounds = []
    
    for n in model_parameters:
        names.append(n)
        bounds.append(model_bounds[n])

    def __init__(self,
                 # angular distribution
                 n_slices=12,
                 dist_slices = "uniform",
                 # observer location
                 view_angle = 15.,
                 view_angle_delta = 1.,
                 # set which components must be activated
                 dyn_flag  = True,
                 wind_flag = True,
                 sec_flag  = True,
                 # values to initialize the global time
                 time_min = 3600.,      # one hour
                 time_max = 2000000.,   # 30 days
                 n_time   = 200,
                 tscale   = 'linear',
                 # values to initialize the ray velocity
                 v_min  = 1.e-6,
                 n_v    = 200,
                 vscale = 'linear',
                 # nuclear heating values
                 eps_ye_dep = True,
                 eps0 = 1.2e+19,
                 sigma0 = 0.11,
                 alpha = 1.3,
                 t0eps = 1.3,
                 cnst_eff = 0.3333,
                 a_eps_nuc = 0.5,
                 b_eps_nuc = 2.5,
                 t_eps_nuc = 1.,
                 # mass dynamic ejecta
                 mass_dist_law_dyn   ='sin2',
                 step_angle_mass_dyn =None,
                 high_lat_flag_dyn   =None,
                 # velocity dynamic ejecta
                 vel_dist_law_dyn    ='uniform',
                 min_vel_dyn         =None,
                 max_vel_dyn         =None,
                 step_angle_vel_dyn  =None,
                 high_lat_vel_dyn    =None,
                 low_lat_vel_dyn     =None,
                 # opacity dynamic ejecta
                 kappa_dist_law_dyn  ='step',
                 central_op_dyn      =None,
                 min_op_dyn          =None,
                 max_op_dyn          =None,
                 step_angle_op_dyn   =math.radians(45.),
                 high_lat_op_dyn     =0.1,
                 low_lat_op_dyn      =10.,
                 # mass wind ejecta
                 mass_dist_law_wind  ='step',
                 m_ej_wind           =0.01,
                 step_angle_mass_wind=math.radians(60.),
                 high_lat_flag_wind  =1,
                 # velocity wind ejecta
                 vel_dist_law_wind   ='uniform',
                 central_vel_wind    =0.05,
                 min_vel_wind        =None,
                 max_vel_wind        =None,
                 step_angle_vel_wind =None,
                 high_lat_vel_wind   =None,
                 low_lat_vel_wind    =None,
                 # opacity wind ejecta
                 kappa_dist_law_wind ='step',
                 central_op_wind     =None,
                 min_op_wind         =None,
                 max_op_wind         =None,
                 step_angle_op_wind  =math.radians(45.),
                 high_lat_op_wind    =0.1,
                 low_lat_op_wind     =1.0,
                 # mass secular ejecta
                 mass_dist_law_sec   ='sin2',
                 m_ej_sec            =0.02,
                 step_angle_mass_sec =None,
                 high_lat_flag_sec   =None,
                 # velocity secular ejecta
                 vel_dist_law_sec    ='uniform',
                 central_vel_sec     =0.02,
                 min_vel_sec         =None,
                 max_vel_sec         =None,
                 step_angle_vel_sec  =None,
                 high_lat_vel_sec    =None,
                 low_lat_vel_sec     =None,
                 # opacity secular ejecta
                 kappa_dist_law_sec  ='uniform',
                 central_op_sec      =5.0,
                 min_op_sec          =None,
                 max_op_sec          =None,
                 step_angle_op_sec   =None,
                 high_lat_op_sec     =None,
                 low_lat_op_sec      =None,
                 **kwargs):

        super(MacroKilonovaModel,self).__init__(**kwargs)
    
        self.n_slices = n_slices
        self.dist_slices = dist_slices
        # initialize the angular distribution
        self.AD = ad.AngularDistribution(dist_slices,n_slices)
        self.ang_dist, self.omega_dist = self.AD(self.n_slices/2)

        print('')
        print('I have initialized the angles')

        # initialize the filters
        self.FT = ft.Filters("measures")
        self.dic_filt,self.lambda_vec,self.mag = self.FT()

        print('')
        print('I have initialized the filters')

        #initialize the observer location
        self.view_angle = view_angle
        self.view_angle_delta = view_angle_delta
        self.FF = op.ObserverProjection(self.n_slices)
        self.flux_factor = self.FF(self.view_angle,self.view_angle_delta)

        print('')
        print('I have initialized the observer location')

        # set which components must be activated
        self.dyn_flag  = dyn_flag
        self.wind_flag = wind_flag
        self.sec_flag  = sec_flag

        # values to initialize the global time
        self.time_min = time_min      # one hour
        self.time_max = time_max   # 30 days
        self.n_time   = n_time
        self.tscale   = 'linear'

        # values to initialize the ray velocity
        self.v_min  = v_min
        self.n_v    = n_v
        self.vscale = 'linear'

        # nuclear heating values
        self.eps_ye_dep = eps_ye_dep
        self.eps0 = eps0
        self.sigma0 = sigma0
        self.alpha = alpha
        self.t0eps = t0eps
        self.cnst_eff = cnst_eff
        self.a_eps_nuc = a_eps_nuc
        self.b_eps_nuc = b_eps_nuc
        self.t_eps_nuc = t_eps_nuc
        
        # mass dynamic ejecta
        self.mass_dist_law_dyn   = mass_dist_law_dyn
        self.step_angle_mass_dyn = step_angle_mass_dyn
        self.high_lat_flag_dyn   = high_lat_flag_dyn
        
        # velocity dynamic ejecta
        self.vel_dist_law_dyn    = vel_dist_law_dyn
        self.min_vel_dyn         = min_vel_dyn
        self.max_vel_dyn         = max_vel_dyn
        self.step_angle_vel_dyn  = step_angle_vel_dyn
        self.high_lat_vel_dyn    = high_lat_vel_dyn
        self.low_lat_vel_dyn     = low_lat_vel_dyn
        
        # opacity dynamic ejecta
        self.kappa_dist_law_dyn  = kappa_dist_law_dyn
        self.central_op_dyn      = central_op_dyn
        self.min_op_dyn          = min_op_dyn
        self.max_op_dyn          = max_op_dyn
        self.step_angle_op_dyn   = step_angle_op_dyn
        self.high_lat_op_dyn     = high_lat_op_dyn
        self.low_lat_op_dyn      = low_lat_op_dyn
        
        # mass wind ejecta
        self.mass_dist_law_wind  = mass_dist_law_wind
        self.step_angle_mass_wind= step_angle_mass_wind
        self.high_lat_flag_wind  = high_lat_flag_wind
        # velocity wind ejecta
        self.vel_dist_law_wind   = vel_dist_law_wind
        self.min_vel_wind        = min_vel_wind
        self.max_vel_wind        = max_vel_wind
        self.step_angle_vel_wind = step_angle_vel_wind
        self.high_lat_vel_wind   = high_lat_vel_wind
        self.low_lat_vel_wind    = low_lat_vel_wind
        # opacity wind ejecta
        self.kappa_dist_law_wind = kappa_dist_law_wind
        self.central_op_wind     = central_op_wind
        self.min_op_wind         = min_op_wind
        self.max_op_wind         = max_op_wind
        self.step_angle_op_wind  = step_angle_op_wind
        self.high_lat_op_wind    = high_lat_op_wind
        self.low_lat_op_wind     = low_lat_op_wind
        # mass secular ejecta
        self.mass_dist_law_sec   = mass_dist_law_sec
        self.step_angle_mass_sec = step_angle_mass_sec
        self.high_lat_flag_sec   = high_lat_flag_sec
        # velocity secular ejecta
        self.vel_dist_law_sec    = vel_dist_law_sec
        self.min_vel_sec         = min_vel_sec
        self.max_vel_sec         = max_vel_sec
        self.step_angle_vel_sec  = step_angle_vel_sec
        self.high_lat_vel_sec    = high_lat_vel_sec
        self.low_lat_vel_sec     = low_lat_vel_sec
        # opacity secular ejecta
        self.kappa_dist_law_sec  = kappa_dist_law_sec
        self.central_op_sec      = central_op_sec
        self.min_op_sec          = min_op_sec
        self.max_op_sec          = max_op_sec
        self.step_angle_op_sec   = step_angle_op_sec
        self.high_lat_op_sec     = high_lat_op_sec
        self.low_lat_op_sec      = low_lat_op_sec
    
    def log_likelihood(self,x):
        # compute the lightcurve
        time,r_ph_tot,L_bol_tot,T_eff_tot = lc.lightcurve(self.dyn_flag,
                                                          self.wind_flag,
                                                          self.sec_flag,
                                                          self.ang_dist,
                                                          self.omega_dist,
                                                          'GK',
                                                          'BKWM',
                                                          self.time_min,
                                                          self.time_max,
                                                          self.n_time,
                                                          self.tscale,
                                                          self.v_min,
                                                          self.n_v,
                                                          self.vscale,
                                                          self.eps_ye_dep,
                                                          self.eps0,
                                                          self.sigma0,
                                                          self.alpha,
                                                          self.t0eps,
                                                          self.cnst_eff,
                                                          self.a_eps_nuc,
                                                          self.b_eps_nuc,
                                                          self.t_eps_nuc,
                                                          # mass dynamic ejecta
                                                          mass_dist_law_dyn   = self.mass_dist_law_dyn,
                                                          m_ej_dyn            = x['m_ej_dyn'],
                                                          step_angle_mass_dyn = self.step_angle_mass_dyn,
                                                          high_lat_flag_dyn   = self.high_lat_flag_dyn,
                                                          # velocity dynamic ejecta
                                                          vel_dist_law_dyn    = self.vel_dist_law_dyn,
                                                          central_vel_dyn     = x['central_vel_dyn'],
                                                          min_vel_dyn         = self.min_vel_dyn,
                                                          max_vel_dyn         = self.max_vel_dyn,
                                                          step_angle_vel_dyn  = self.step_angle_vel_dyn,
                                                          high_lat_vel_dyn    = self.high_lat_vel_dyn,
                                                          low_lat_vel_dyn     = self.low_lat_vel_dyn,
                                                          # opacity dynamic ejecta
                                                          kappa_dist_law_dyn  = self.kappa_dist_law_dyn,
                                                          central_op_dyn      = x['central_op_dyn'],
                                                          min_op_dyn          = self.min_op_dyn,
                                                          max_op_dyn          = self.max_op_dyn,
                                                          step_angle_op_dyn   = self.step_angle_op_dyn,
                                                          high_lat_op_dyn     = self.high_lat_op_dyn,
                                                          low_lat_op_dyn      = self.low_lat_op_dyn,
                                                          # mass wind ejecta
                                                          mass_dist_law_wind  = self.mass_dist_law_wind,
                                                          m_ej_wind           = x['m_ej_wind'],
                                                          step_angle_mass_wind= self.step_angle_mass_wind,
                                                          high_lat_flag_wind  = self.high_lat_flag_wind,
                                                          # velocity wind ejecta
                                                          vel_dist_law_wind   = self.vel_dist_law_wind,
                                                          central_vel_wind    = x['central_vel_wind'],
                                                          min_vel_wind        = self.min_vel_wind,
                                                          max_vel_wind        = self.max_vel_wind,
                                                          step_angle_vel_wind = self.step_angle_vel_wind,
                                                          high_lat_vel_wind   = self.high_lat_vel_wind,
                                                          low_lat_vel_wind    = self.low_lat_vel_wind,
                                                          # opacity wind ejecta
                                                          kappa_dist_law_wind = self.kappa_dist_law_wind,
                                                          central_op_wind     = x['central_op_wind'],
                                                          min_op_wind         = self.min_op_wind,
                                                          max_op_wind         = self.max_op_wind,
                                                          step_angle_op_wind  = self.step_angle_op_wind,
                                                          high_lat_op_wind    = self.high_lat_op_wind,
                                                          low_lat_op_wind     = self.low_lat_op_wind,
                                                          # mass secular ejecta
                                                          mass_dist_law_sec   = self.mass_dist_law_sec,
                                                          m_ej_sec            = x['m_ej_sec'],
                                                          step_angle_mass_sec = self.step_angle_mass_sec,
                                                          high_lat_flag_sec   = self.high_lat_flag_sec,
                                                          # velocity secular ejecta
                                                          vel_dist_law_sec    = self.vel_dist_law_sec,
                                                          central_vel_sec     = x['central_vel_sec'],
                                                          min_vel_sec         = self.min_vel_sec,
                                                          max_vel_sec         = self.max_vel_sec,
                                                          step_angle_vel_sec  = self.step_angle_vel_sec,
                                                          high_lat_vel_sec    = self.high_lat_vel_sec,
                                                          low_lat_vel_sec     = self.low_lat_vel_sec,
                                                          # opacity secular ejecta
                                                          kappa_dist_law_sec  = self.kappa_dist_law_sec,
                                                          central_op_sec      = x['central_op_sec'],
                                                          min_op_sec          = self.min_op_sec,
                                                          max_op_sec          = self.max_op_sec,
                                                          step_angle_op_sec   = self.step_angle_op_sec,
                                                          high_lat_op_sec     = self.high_lat_op_sec,
                                                          low_lat_op_sec      = self.low_lat_op_sec)
        # compute the magnitudes from a certain distance
        D = 40.e+6*units.pc2cm
        model_mag = ft.calc_magnitudes(self.flux_factor,time,r_ph_tot,T_eff_tot,self.lambda_vec,self.dic_filt,D)

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
            return logP
        else:
            return -np.inf

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to gw150914/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=1,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--nlive',default=1024,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=10,type='int',metavar='m',help='max MCMC points')
    parser.add_option('--poolsize',default=16,type='int',metavar='k',help='numer of points in the ensemble sampler pool')
    (opts,args)=parser.parse_args()

    if opts.out_dir is None:
        opts.out_dir='./test/'

    if opts.full_run:
        model = MacroKilonovaModel()
        work=cpnest.CPNest(model,
                           verbose=3,
                           Poolsize=opts.poolsize,
                           Nthreads=opts.threads,
                           Nlive=opts.nlive,
                           maxmcmc=opts.maxmcmc,
                           output=opts.out_dir)
        work.run()

