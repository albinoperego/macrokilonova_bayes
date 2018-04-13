import angular_distribution as ad
import filters as ft
import numpy as np
import observer_projection as op
import units
import math
import ejecta as ej
import matplotlib.pyplot as plt

# global properties of the source
D = 40.e+6*units.pc2cm
t0 = 57982.529

# initialize the angular distribution
n_slices = 12
dist_slices = "cos_uniform"
AD = ad.AngularDistribution(dist_slices,n_slices)
angular_distribution, omega_distribution = AD(n_slices/2)

print('')
print('I have initialized the angles')

# initialize the filters
FT = ft.Filters("measures")
dic_filt,lambda_vec,mag = FT()

print('')
print('I have initialized the filters')

#initialize the observer location
view_angle = 90.
FF = op.ObserverProjection(n_slices,dist_slices)
flux_factor = FF(view_angle)#,view_angle_delta)

print('')
print('I have initialized the observer location')

    #initialize the time
time_min = 3600.      #
time_max = 2000000.   #
n_time = 200
tscale   = 'log'
# initialize global time
if (tscale == 'linear'):
    time = np.linspace(time_min,time_max,num=n_time)
elif (tscale == 'log'):
    time = np.logspace(np.log10(time_min),np.log10(time_max),num=n_time)
else:
    print('Error! Wrong option for the time scale')
    exit(-1)

print('')
print('I have initialized the global time')

params = {}
params['dynamics']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM_1d', 'eps_ye_dep':True}
params['wind']    = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM_1d', 'eps_ye_dep':True}
params['secular'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM_1d', 'eps_ye_dep':True}

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
               'lc model':'villar'}  #villar or grossman

glob_vars = {'m_disk':0.09,
             'eps0':1.5e19, 
             'T_floor_LA':1000., 
             'T_floor_Ni':3000., 
             'a_eps_nuc':0.5,
             'b_eps_nuc':2.5,
             't_eps_nuc':1.0}

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

write_output = True
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

plot_separately = True           # Choose to plot all lightcurves in different bands on the same plot 
plot_together = False	         # or to plot lightcurve and data in each band on different plots

if (plot_separately):			# plot lightcurve and data in each band separately
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
elif (plot_together) :							# plot lightcurves for every band in the same plot
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
    
