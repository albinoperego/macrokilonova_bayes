import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import angular_distribution as ad
import filters as ft
import lightcurve as lc
import observer_projection as op
import units
import ejecta as ej

# initialize the filters
print("Initialising filters")
FT = ft.Filters("measures")
dic_filt,lambda_vec,mag = FT()
time_min = 3600.      #
time_max = 2000000.   #
n_time = 20
# initialize global time
time = np.logspace(np.log10(time_min),np.log10(time_max),num=n_time)

# read the chain and plot the maximum likelihood point, for now
chain = np.genfromtxt('quick_test/chain_1024_1234.txt',names=True)[-1]

# set of global parameters not to be fit
glob_params = {'v_min':1.e-7,
               'n_v':10,
               'vscale':'linear',
               'sigma0':0.11,
               'alpha':1.3,
               't0eps':1.3,
               'cnst_eff':0.3333}

# set of global parameters to be fit
glob_vars = {'m_disk':chain['m_disk'],
             'eps0':  chain['eps0'],
             'a_eps_nuc':chain['a_eps_nuc'],
             'b_eps_nuc':chain['b_eps_nuc'],
             't_eps_nuc':chain['t_eps_nuc']}

# hardcoded ejecta geometric and thermal parameters
ejecta_params = {}
ejecta_params['dynamics']   = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
ejecta_params['wind']       = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
ejecta_params['secular']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}

# set of shell parameters to be sampled on
shell_params={}

shell_params['dynamics'] = {'xi_disk':chain['xi_disk_dynamics'],
                            'm_ej':None,#chain['m_ej_dynamics'],
                          'central_vel':chain['central_vel_dynamics'],
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':chain['low_lat_op_dynamics'],
                          'high_lat_op':chain['high_lat_op_dynamics'],
                          'step_angle_op':chain['step_angle_op_dynamics']}

shell_params['wind'] = {'xi_disk':chain['xi_disk_wind'],
                      'm_ej':None,
                      'step_angle_mass':chain['step_angle_mass_wind'],
                      'high_lat_flag':True,
                      'central_vel':chain['central_vel_wind'],
                      'low_lat_vel':None,
                      'high_lat_vel':None,
                      'step_angle_vel':None,
                      'low_lat_op':chain['low_lat_op_wind'],
                      'high_lat_op':chain['high_lat_op_wind'],
                      'step_angle_op':chain['step_angle_op_wind']}

shell_params['secular'] = {'xi_disk':chain['xi_disk_secular'],
                         'm_ej':None,
                         'central_vel':chain['central_vel_secular'],
                         'low_lat_vel':None,
                         'high_lat_vel':None,
                         'step_angle_vel':None,
                         'central_op':chain['central_op_secular'],
                         'low_lat_op':None,
                         'high_lat_op':None,
                         'step_angle_op':None}
# m_disk    a_eps_nuc    t_eps_nuc    eps0    b_eps_nuc    low_lat_op_dynamics    xi_disk_dynamics    step_angle_op_dynamichigh_lat_op_dynamics    central_vel_dynamics    xi_disk_secular    central_vel_secular    central_op_secular    low_lat_op_wind    xi_disk_wind    step_angle_mass_wind    central_vel_wind    step_angle_op_wind    high_lat_op_wind
#                                                                                        -4.350513e+04

#posteriors = np.genfromtxt('test/posterior.dat',names=True)
# number and distribution of slices for the angular integrals
n_slices = 12
dist_slices = 'uniform'
# initialize the angular distribution
print("Initialising angular distribution")
AD = ad.AngularDistribution(dist_slices,n_slices)
ang_dist, omega_dist = AD(n_slices/2)
FF = op.ObserverProjection(n_slices,dist_slices)
flux_factor = FF(chain['view_angle'])

ejecta = ej.Ejecta(len(ejecta_params), ejecta_params.keys(), ejecta_params)
r_ph, L_bol, T_eff = ejecta.lightcurve(ang_dist, omega_dist,
                                        time,
                                        shell_params,
                                        glob_vars,
                                        glob_params)
    
# compute the magnitudes from a certain distance
D = chain['distance']*1e6*units.pc2cm
model_mag = ft.calc_magnitudes(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,D)
print model_mag.keys()
fig = plt.figure()
ax = fig.add_subplot(111)
band_list = ['U','V','B','R','K','J']
colors = iter(cm.rainbow_r(np.linspace(0, 1, len(band_list))))

for ilambda in model_mag.keys():
    if ilambda in dic_filt.keys():
        if(dic_filt[ilambda]['name'] in band_list):
            c = next(colors)
            ax.errorbar(mag[ilambda]['time']-57982.529,mag[ilambda]['mag'],yerr=mag[ilambda]['sigma'],color=c, label=dic_filt[ilambda]['name'])
            ax.plot(time/24./60./60.,model_mag[ilambda],color=c,linestyle='dashed')
    else: print "key",ilambda,"not found"
plt.xscale('log')
plt.title('UV and visible bands')
plt.xlabel('Time [day]')
plt.ylabel('AB magnitude')
plt.legend(loc='0')
plt.xlim(0.1,30)
plt.ylim(27,15)
plt.show()
