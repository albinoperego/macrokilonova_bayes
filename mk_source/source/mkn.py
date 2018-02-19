import angular_distribution as ad
import filters as ft
import lightcurve as lc
import numpy as np
import observer_projection as op
import units
import math
import ejecta as ej
import matplotlib.pyplot as plt

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
params['dynamics']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
params['wind']    = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
params['secular'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform', 'therm_model':'BKWM', 'eps_ye_dep':True}

E = ej.Ejecta(3, params.keys(), params)

shell_vars={}

shell_vars['dynamics'] = {'xi_disk':None,
                      'm_ej':0.0004,
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
               'cnst_eff':0.3333}

glob_vars = {'m_disk':0.0009,
             'eps0':1.5e19, 
             'a_eps_nuc':0.5,
             'b_eps_nuc':2.5,
             't_eps_nuc':1.0}

r_ph, L_bol, T_eff = E.lightcurve(angular_distribution,
                           omega_distribution,
                           time,
                           shell_vars,
                           glob_vars,
                           glob_params)

# compute the magnitudes from a certain distance
D = 40.e+6*units.pc2cm
model_mag = ft.calc_magnitudes(flux_factor,time,r_ph,T_eff,lambda_vec,dic_filt,D)

# compute the residuals
residuals = ft.calc_residuals(mag,model_mag)

# compute the likelihood
logL = 0.
for ilambda in residuals.keys():
    logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

print('')
print('logL')
print(logL)

write_output = True
if (write_output):
    g = open('mkn_SFHo_M144139_LK_90.txt','w')

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

fig1 = plt.figure()
#band_list = ['U']
band_list = ['U','B','g','V','R','r']
for band in band_list:
    for ilambda in mag.keys():
        if(mag[ilambda]['name']==band):
            plt.plot(time/24./60./60.,model_mag[ilambda])
#            plt.semilogx(time/24./60./60.,model_mag[ilambda])
#            print('time')
#            print(mag[ilambda]['time']/24./60./60.)
#            print('mag')
#            print(mag[ilambda]['mag'])
#            plt.scatter(mag[ilambda]['time']-57982.529,mag[ilambda]['mag'])

plt.xlim(0.1,10)
plt.ylim(27,15)


fig2 = plt.figure()
#band_list = ['U']
band_list = ['i','z','J','H','Ks']
for band in band_list:
    for ilambda in mag.keys():
        if(mag[ilambda]['name']==band):
            plt.plot(time/24./60./60.,model_mag[ilambda])
#            print('time')
#            print(mag[ilambda]['time']/24./60./60.)
#            print('mag')
#            print(mag[ilambda]['mag'])
#            plt.scatter(mag[ilambda]['time']-57982.529,mag[ilambda]['mag'])

plt.xlim(0.1,10)
plt.ylim(27,15)
plt.show()


print_output = False
if (print_output):
    plt.plot(time, model_mag[7],'-')
    plt.show()
    


#import matplotlib.pyplot as plt
#print np.shape(r_ph), np.shape(L_bol)
#for j in range(E.ncomponents):
#    for k in range(L_bol.shape[1]): plt.plot(time, L_bol[j, 0, :],'.')
#plt.show()

## set which components must be activated
#dyn_flag  = True
#wind_flag = True
#sec_flag  = True
#
## values to initialize the global time
#time_min = 3600.      # one hour
#time_max = 2000000.   # 30 days
#n_time   = 200
#tscale   = 'linear'
#
## values to initialize the ray velocity
#v_min  = 1.e-6
#n_v    = 200
#vscale = 'linear'
#
## nuclear heating values
#eps_ye_dep = True
#eps0 = 1.2e+19
#sigma0 = 0.11
#alpha = 1.3
#t0eps = 1.3
#cnst_eff = 0.3333
#a_eps_nuc = 0.5
#b_eps_nuc = 2.5
#t_eps_nuc = 1.
#
## compute the lightcurve
#time,r_ph_tot,L_bol_tot,T_eff_tot = lc.lightcurve(dyn_flag,wind_flag,sec_flag,ang_dist,omega_dist,'GK','BKWM',time_min,time_max,n_time,tscale,v_min,n_v,vscale,eps_ye_dep,eps0,sigma0,alpha,t0eps,cnst_eff,a_eps_nuc,b_eps_nuc,t_eps_nuc,
## mass dynamic ejecta
#        mass_dist_law_dyn   ='sin2',
#        m_ej_dyn            =0.005,
#        step_angle_mass_dyn =None,
#        high_lat_flag_dyn   =None,
## velocity dynamic ejecta
#        vel_dist_law_dyn    ='uniform',
#        central_vel_dyn     =0.2,
#        min_vel_dyn         =None,
#        max_vel_dyn         =None,
#        step_angle_vel_dyn  =None,
#        high_lat_vel_dyn    =None,
#        low_lat_vel_dyn     =None,
## opacity dynamic ejecta
#        kappa_dist_law_dyn  ='step',
#        central_op_dyn      =None,
#        min_op_dyn          =None,
#        max_op_dyn          =None,
#        step_angle_op_dyn   =math.radians(45.),
#        high_lat_op_dyn     =0.1,
#        low_lat_op_dyn      =10.,
## mass wind ejecta
#        mass_dist_law_wind  ='step',
#        m_ej_wind           =0.01,
#        step_angle_mass_wind=math.radians(60.),
#        high_lat_flag_wind  =1,
## velocity wind ejecta
#        vel_dist_law_wind   ='uniform',
#        central_vel_wind    =0.05,
#        min_vel_wind        =None,
#        max_vel_wind        =None,
#        step_angle_vel_wind =None,
#        high_lat_vel_wind   =None,
#        low_lat_vel_wind    =None,
## opacity wind ejecta
#        kappa_dist_law_wind ='step',
#        central_op_wind     =None,
#        min_op_wind         =None,
#        max_op_wind         =None,
#        step_angle_op_wind  =math.radians(45.),
#        high_lat_op_wind    =0.1,
#        low_lat_op_wind     =1.0,
## mass secular ejecta
#        mass_dist_law_sec   ='sin2',
#        m_ej_sec            =0.02,
#        step_angle_mass_sec =None,
#        high_lat_flag_sec   =None,
## velocity secular ejecta
#        vel_dist_law_sec    ='uniform',
#        central_vel_sec     =0.02,
#        min_vel_sec         =None,
#        max_vel_sec         =None,
#        step_angle_vel_sec  =None,
#        high_lat_vel_sec    =None,
#        low_lat_vel_sec     =None,
## opacity secular ejecta
#        kappa_dist_law_sec  ='uniform',
#        central_op_sec      =5.0,
#        min_op_sec          =None,
#        max_op_sec          =None,
#        step_angle_op_sec   =None,
#        high_lat_op_sec     =None,
#        low_lat_op_sec      =None)
#
## compute the magnitudes from a certain distance
#D = 40.e+6*units.pc2cm
#model_mag = ft.calc_magnitudes(flux_factor,time,r_ph_tot,T_eff_tot,lambda_vec,dic_filt,D)
#
## compute the residuals
#residuals = ft.calc_residuals(mag,model_mag)
#
## compute the likelihood
#logL = 0.
#for ilambda in residuals.keys():
#    logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))
#
#print('')
#print('logL')
#print(logL)
#
#write_output = True
#if (write_output):
#    g = open('mkn_output.txt','w')
#
#    g.write('%20s' %('time'))
#    for ilambda in mag.keys():
#        g.write('%20s' %(mag[ilambda]['name']))
#    g.write('\n')
#    
#    for i in range(len(time)):
#        g.write('%20s' %(time[i]))
#        for ilambda in model_mag.keys():
#            if (ilambda == 0):
#                continue
#            g.write('%20s' %(model_mag[ilambda][i]))
#        g.write('\n')
#
#    g.close()
