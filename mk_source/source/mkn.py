import angular_distribution as ad
import ejecta as ej
import filters as ft
import math
import matplotlib.pyplot as plt
import numpy as np
import observer_projection as op
import source_properties as sp
import units

#dictionary with the global parameters of the model
glob_params = {'lc model'   :'grossman',    # model for the lightcurve (grossman or villar)  
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
               'scale for t':'measures'     # kind of spacing in time [log - linear - measures]
               }

source_name = 'AT2017gfo'   # name of the source or "default"
#source_name = 'default'   # name of the source or "default"

# dictionary for the parameters of each shell
params = {}
params['dynamics'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM', 'eps_ye_dep':True}
params['wind']     = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step'   , 'therm_model':'BKWM', 'eps_ye_dep':True}
params['secular']  = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}
# BKWM_1d


# dictionary for the variables of each shell
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

# dictionary for the global variables
glob_vars = {'m_disk':0.09,
             'eps0':1.5e19, 
             'T_floor_LA':1000., 
             'T_floor_Ni':3000., 
             'a_eps_nuc':0.5,
             'b_eps_nuc':2.5,
             't_eps_nuc':1.0}


#############################
#  LIGHT CURVE CALCULATION  #
#############################


# initializing the global properties of the source
SP = sp.SourceProperties(source_name)
print('')
print('I have initialized the global properties of the source')

# initialize the angular distribution
n_slices = glob_params['n slices']
dist_slices = glob_params['dist slices']
AD = ad.AngularDistribution(dist_slices,n_slices)
angular_distribution, omega_distribution = AD(n_slices/2)   # due to the symmetry abount the equatorial plane, the number of independent slices is half
print('')
print('I have initialized the angles')

# initialize the filters
if (source_name == 'default'):
    FT = ft.Filters("properties")
else:
    FT = ft.Filters("measures")

dic_filt,lambda_vec,mag = FT(SP.filter_data_folder)
print('')
print('I have initialized the filters')

#initialize the observer location
FF = op.ObserverProjection(n_slices,dist_slices)
flux_factor = FF(SP.view_angle)
print('')
print('I have initialized the observer orientation')

#initialize the time
time_min = glob_params['time min']     
time_max = glob_params['time max']     
n_time   = glob_params['n time']       
tscale   = glob_params['scale for t']
if (tscale == "measures" and source_name=='default'):
    print('')
    print("no measures available to set the time (default option)")
    print("please use linear or log scale")
    exit() 

time = SP.init_time(tscale,time_min,time_max,n_time,mag)
print('')
print('I have initialized the global time')

E = ej.Ejecta(3, params.keys(), params)

# compute the lightcurve
r_ph, L_bol, T_eff = E.lightcurve(angular_distribution,
                           omega_distribution,
                           time,
                           shell_vars,
                           glob_vars,
                           glob_params)
print('')
print('I have computed the lightcurve')


# compute the residuals
if (source_name != 'default'):
    residuals = ft.calc_all_residuals(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,SP.D,SP.t0,mag)

    print('')
    print('I have computed the residuals')

    # compute the likelihood
    logL = 0.
    for ilambda in residuals.keys():
        logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

    print('')
    print('I have computed the likehood')
    print('logL')
    print(logL)

#########################
# write out the output  #
#########################

write_output = True
if (write_output):

    print('')
    print('I am printing out the output')

    model_mag = ft.calc_magnitudes(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,SP.D,SP.t0)

    file_output = 'mkn_model.txt'
    g = open(file_output,'w')

    print('file name:',file_output)

    g.write('%20s' %('time'))
    for ilambda in dic_filt.keys():
        g.write('%20s' %(dic_filt[ilambda]['name']))
    g.write('\n')
    
    for i in range(len(time)):
        g.write('%20s' %(time[i]))
        for ilambda in model_mag.keys():
            if (ilambda == 0):
                continue
            g.write('%20s' %(model_mag[ilambda][i]))
        g.write('\n')

    g.close()


################################
# plot some of the lightcurves #
################################

plot_separately = False           # Choose to plot all lightcurves in different bands on the same plot
plot_together = False             # or to plot lightcurve and data in each band on different plots

if (plot_separately):            # plot lightcurve and data in each band separately
    fig1 = plt.figure()
    for ilambda in dic_filt.keys():
        if (dic_filt[ilambda]['plot'] !=1):
            continue

        plt.plot(time*units.sec2day,model_mag[ilambda])

        if (source_name != 'default'):
            if(len(misure[ilambda]['sigma'])!=0):
                plt.errorbar(misure[ilambda]['time']-SP.t0,misure[ilambda]['mag'],yerr=misure[ilambda]['sigma'],fmt='o')
        plt.title('%s' %dic_filt[ilambda]['name'])
        plt.xlim(0.1,10)
        plt.ylim(27,15)
        plt.show()

elif (plot_together):        # plot lightcurves for every band in the same plot
    fig1 = plt.figure()
    for ilambda in dic_filt.keys():
        if (dic_filt[ilambda]['plot'] !=1):
            continue
        plt.plot(time*units.sec2day,model_mag[ilambda])
    plt.title('Lightcurves')
    plt.xlim(0.1,10)
    plt.ylim(27,15)
    plt.show()

