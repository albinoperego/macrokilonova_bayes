import angular_distribution as ad
import filters as ft
import lightcurve as lc
import numpy as np
import observer_projection as op
import units
import math

n_slices = 12
dist_slices = "uniform"
# initialize the angular distribution
AD = ad.AngularDistribution(dist_slices,n_slices)
ang_dist, omega_dist = AD(n_slices/2)

print('')
print('I have initialized the angles')

# initialize the filters
FT = ft.Filters("measures")
dic_filt,lambda_vec,mag = FT()

print('')
print('I have initialized the filters')

#initialize the observer location
view_angle = 15.
view_angle_delta = 1.
FF = op.ObserverProjection(n_slices)
flux_factor = FF(view_angle,view_angle_delta)

print('')
print('I have initialized the observer location')

# set which components must be activated
dyn_flag  = True
wind_flag = True
sec_flag  = True

# values to initialize the global time
time_min = 3600.      # one hour
time_max = 2000000.   # 30 days
n_time   = 200
tscale   = 'linear'

# values to initialize the ray velocity
v_min  = 1.e-6
n_v    = 200
vscale = 'linear'

# compute the lightcurve
time,r_ph_tot,L_bol_tot,T_eff_tot = lc.lightcurve(dyn_flag,wind_flag,sec_flag,ang_dist,omega_dist,'GK','BKWM',time_min,time_max,n_time,tscale,v_min,n_v,vscale,
# mass dynamic ejecta
        mass_dist_law_dyn   ='sin2',
        m_ej_dyn            =0.005,
        step_angle_mass_dyn =None,
        high_lat_flag_dyn   =None,
# velocity dynamic ejecta
        vel_dist_law_dyn    ='uniform',
        central_vel_dyn     =0.2,
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
        low_lat_op_sec      =None)

# compute the magnitudes from a certain distance
D = 40.e+6*units.pc2cm
model_mag = ft.calc_magnitudes(flux_factor,time,r_ph_tot,T_eff_tot,lambda_vec,dic_filt,D)

# compute the residuals
residuals = ft.calc_residuals(mag,model_mag)

# compute the likelihood
logL = 0.
for ilambda in residuals.keys():
    logL += -0.5*np.sum(np.array([res*res for res in residuals[ilambda]]))

print('')
print('logL')
print(logL)

for ilambda in mag.keys():
    print(ilambda,mag[ilambda]['name'])

write_output = True
if (write_output):
    g = open('mkn_output.txt','w')

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
