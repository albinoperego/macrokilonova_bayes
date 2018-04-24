import angular_distribution as ad
import filters as ft
import numpy as np
import observer_projection as op
import units
import math
import ejecta as ej
import matplotlib.pyplot as plt

# global properties of the source
D = 40.e+6*units.pc2cm    # source distance [cm]
t0 = 57982.529            # universal time of the observation [days]

# initializ the angular distribution
n_slices = 30                   # number of slices in which the polar angle is discretized [12-18-24-30]
dist_slices = "cos_uniform"     # discretization law for the polar angle [uniform and cos_uniform]
AD = ad.AngularDistribution(dist_slices,n_slices)
angular_distribution, omega_distribution = AD(n_slices/2)   # due to the symmetry abount the equatorial plane, the number of independent slices is half

print('')
print('I have initialized the angles')

# initialize the filters
FT = ft.Filters("measures")
dic_filt,lambda_vec,mag = FT()

print('')
print('I have initialized the filters')

#initialize the observer location
view_angle = 30.                # viewing angle between the source and the observer
FF = op.ObserverProjection(n_slices,dist_slices)
flux_factor = FF(view_angle)#,view_angle_delta)

print('')
print('I have initialized the observer location')

    #initialize the time
time_min = 3600.      # minimum time [s]
time_max = 2000000.   # maximum time [s]
n_time   = 200        # number of bins in time
#tscale   = 'linear'      # kind of spacing in time [log - linear]
tscale   = 'measures'      # kind of spacing in time [log - linear - measures]
# initialize global time
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
                all_time.append(mag[ilambda]['time'][i]-t0)
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
    time = sorted(np.array(time)*units.day2sec)
    time = np.array(time)
else:
    print('Error! Wrong option for the time scale')
    exit(-1)

print('')
print('I have initialized the global time')

params = {}
params['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM_1d', 'eps_ye_dep':True}
params['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM_1d', 'eps_ye_dep':True}
params['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM_1d', 'eps_ye_dep':True}
# BKWM_1d

E = ej.Ejecta(3, params.keys(), params)

shell_vars={}

shell_vars['dynamics'] = {'xi_disk':None,
                          'm_ej':0.04,
                          'central_vel':0.33,
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':30.,
                          'high_lat_op':10.0,
                          'step_angle_op':math.radians(45.)}

shell_vars['wind'] = {'xi_disk':0.02,
                      'm_ej':None,
                      'step_angle_mass':math.radians(60.),
                      'high_lat_flag':True,
                      'central_vel':0.08,
                      'low_lat_vel':None,
                      'high_lat_vel':None,
                      'step_angle_vel':None,
                      'low_lat_op':0.5,
                      'high_lat_op':5.0,
                      'step_angle_op':math.radians(30.)}

shell_vars['secular'] = {'xi_disk':0.2,
                         'm_ej':None,
                         'central_vel':0.06,
                         'low_lat_vel':None,
                         'high_lat_vel':None,
                         'step_angle_vel':None,
                         'central_op':5.0,
                         'low_lat_op':None,
                         'high_lat_op':None,
                         'step_angle_op':None}

glob_params = {'v_min':1.e-7, 
               'n_v':400, 
               'vscale':'linear',
               'sigma0':0.11,
               'alpha':1.3,
               't0eps':1.3,
               'cnst_eff':0.3333,
               'lc model':'grossman'}  #villar or grossman

glob_vars = {'m_disk':0.09,
             'eps0':1.5e19, 
             'T_floor_LA':1000., 
             'T_floor_Ni':3000., 
             'a_eps_nuc':0.5,
             'b_eps_nuc':2.5,
             't_eps_nuc':1.0}

for j in range(1000):

    r_ph, L_bol, T_eff = E.lightcurve(angular_distribution,
                           omega_distribution,
                           time,
                           shell_vars,
                           glob_vars,
                           glob_params)

# compute the residuals
residuals = ft.calc_all_residuals(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,D,t0,mag)


# compute the likelihood
logL = 0.
for ilambda in residuals.keys():
    logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

print('')
print('logL')
print(logL)

##############################
# write out the output
############################## 

write_output = False
if (write_output):

    model_mag = ft.calc_magnitudes(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,D,t0)

    g = open('mkn_model.txt','w')

    g.write('%20s' %('time'))
    for ilambda in mag.keys():
        g.write('%20s' %(mag[ilambda]['name']))
    g.write('\n')
    
    for i in range(len(time)):
        g.write('%20s' %(time[i]))
        for ilambda in model_mag.keys():
            if (ilambda == 0):
                continue
            g.write('%20s' %(model_mag[ilambda][i]))
        g.write('\n')

    g.close()


###############################
# plot some of the lightcurves 
############################### 

plot_separately = False           # Choose to plot all lightcurves in different bands on the same plot
plot_together = False            # or to plot lightcurve and data in each band on different plots

if (plot_separately):            # plot lightcurve and data in each band separately
    dic,lambdas,misure = ft.read_filter_measures()
    fig1 = plt.figure()
    for ilambda in mag.keys():
        if (dic[ilambda]['plot'] !=1):
            continue
        if(len(misure[ilambda]['sigma'])!=0):
            plt.plot(time*units.sec2day,model_mag[ilambda])
            plt.errorbar(misure[ilambda]['time']-t0,misure[ilambda]['mag'],yerr=misure[ilambda]['sigma'],fmt='o')
            plt.title('%s' %misure[ilambda]['name'])
            plt.xlim(0.1,10)
            plt.ylim(27,15)
            plt.show()
elif (plot_together):        # plot lightcurves for every band in the same plot
    fig1 = plt.figure()
    dic,lambdas,misure = ft.read_filter_measures()
    for ilambda in mag.keys():
        if (dic[ilambda]['plot'] !=1):
            continue
        plt.plot(time*units.sec2day,model_mag[ilambda])
    plt.title('Lightcurves')
    plt.xlim(0.1,10)
    plt.ylim(27,15)
    plt.show()


print_output = False
if (print_output):
    plt.plot(time, model_mag[7],'-')
    plt.show()
    
