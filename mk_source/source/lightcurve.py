import numpy as np
import angular_distribution
import mass_angular_distribution
import thermalization
import velocity_angular_distribution
import opacity_angular_distribution
import initialize_components
import units
import initialize_components

def T_eff_calc(Lum,dOmega,r_ph):
  return (Lum/(dOmega*r_ph**2*units.sigma_SB))**(1./4.)

def lightcurve(dyn_flag,wind_flag,sec_flag,ang_dist,omega_dist,exp_model,therm_model,**kwargs):

# initialize global time
    time_min = 36000.      #
    time_max = 100000.   #
    n_time = 200
    time = np.linspace(time_min,time_max,n_time)

    n_ang = len(ang_dist)

    if(dyn_flag == True):

        mass_dist_law_d  = kwargs['mass_dist_law_dyn']
        m_ej_d           = kwargs['m_ej_dyn']
        step_angle_mass_d  = kwargs['step_angle_mass_dyn']
        high_lat_flag_d  = kwargs['high_lat_flag_dyn']
        vel_dist_law_d   = kwargs['vel_dist_law_dyn']
        central_vel_d    = kwargs['central_vel_dyn'] 
        min_vel_d        = kwargs['min_vel_dyn']      
        max_vel_d        = kwargs['max_vel_dyn']      
        step_angle_vel_d   = kwargs['step_angle_vel_dyn'] 
        high_lat_vel_d   = kwargs['high_lat_vel_dyn'] 
        low_lat_vel_d    = kwargs['low_lat_vel_dyn']  
        kappa_dist_law_d = kwargs['kappa_dist_law_dyn']
        central_op_d     = kwargs['central_op_dyn']    
        min_op_d         = kwargs['min_op_dyn']        
        max_op_d         = kwargs['max_op_dyn']        
        step_angle_op_d    = kwargs['step_angle_op_dyn']   
        high_lat_op_d    = kwargs['high_lat_op_dyn']  
        low_lat_op_d     = kwargs['low_lat_op_dyn']    


        M,V,O = initialize_components.InitializeComponent(mass_dist_law_d,vel_dist_law_d,kappa_dist_law_d)
        ET = thermalization.Thermalization(therm_model)
        r_ph_dyn, L_bol_dyn = initialize_components.expansion_angular_distribution(M,V,O,ET,exp_model,ang_dist,omega_dist,m_ej_d,time,
          step_angle_mass = step_angle_mass_d,
          high_lat_flag = high_lat_flag_d,
          vel_dist_law  = vel_dist_law_d,
          central_vel   = central_vel_d,
          min_vel       = min_vel_d,
          max_vel       = max_vel_d,
          step_angle_vel  = step_angle_vel_d,
          high_lat_vel  = high_lat_vel_d,
          low_lat_vel   = low_lat_vel_d,
          kappa_dist_law= kappa_dist_law_d,
          central_op    = central_op_d,
          min_op        = min_op_d,
          max_op        = max_op_d,
          step_angle_op   = step_angle_op_d, 
          high_lat_op   = high_lat_op_d,
          low_lat_op    = low_lat_op_d)

    else:
        r_ph_dyn = np.full((n_ang,n_time),units.small)
        L_bol_dyn = np.full((n_ang,n_time),units.small)


    if(wind_flag == True):

        mass_dist_law_w  = kwargs['mass_dist_law_wind']
        m_ej_w           = kwargs['m_ej_wind']
        step_angle_mass_w  = kwargs['step_angle_mass_wind']
        high_lat_flag_w  = kwargs['high_lat_flag_wind']
        vel_dist_law_w   = kwargs['vel_dist_law_wind']
        central_vel_w    = kwargs['central_vel_wind'] 
        min_vel_w        = kwargs['min_vel_wind']      
        max_vel_w        = kwargs['max_vel_wind']      
        step_angle_vel_w   = kwargs['step_angle_vel_wind'] 
        high_lat_vel_w   = kwargs['high_lat_vel_wind'] 
        low_lat_vel_w    = kwargs['low_lat_vel_wind']  
        kappa_dist_law_w = kwargs['kappa_dist_law_wind']
        central_op_w     = kwargs['central_op_wind']    
        min_op_w         = kwargs['min_op_wind']        
        max_op_w         = kwargs['max_op_wind']        
        step_angle_op_w    = kwargs['step_angle_op_wind']   
        high_lat_op_w    = kwargs['high_lat_op_wind']  
        low_lat_op_w     = kwargs['low_lat_op_wind']    

        M,V,O = initialize_components.InitializeComponent(mass_dist_law_w,vel_dist_law_w,kappa_dist_law_w)


        ET = thermalization.Thermalization(therm_model)

        r_ph_wind, L_bol_wind = initialize_components.expansion_angular_distribution(M,V,O,ET,exp_model,ang_dist,omega_dist,m_ej_w,time,step_angle_mass=step_angle_mass_w,high_lat_flag = high_lat_flag_w,vel_dist_law=vel_dist_law_w,central_vel = central_vel_w,min_vel= min_vel_w,max_vel=max_vel_w,step_angle_vel=step_angle_vel_w,high_lat_vel=high_lat_vel_w,low_lat_vel=low_lat_vel_w,kappa_dist_law= kappa_dist_law_w,central_op=central_op_w,min_op=min_op_w,max_op=max_op_w,step_angle_op=step_angle_op_w,high_lat_op=high_lat_op_w,low_lat_op=low_lat_op_w)

    else:
        r_ph_wind = np.full((n_ang,n_time),units.small)
        L_bol_wind = np.full((n_ang,n_time),units.small)

    if(sec_flag == True):
        mass_dist_law_s  = kwargs['mass_dist_law_sec']
        m_ej_s           = kwargs['m_ej_sec']
        step_angle_mass_s  = kwargs['step_angle_mass_sec']
        high_lat_flag_s  = kwargs['high_lat_flag_sec']
        vel_dist_law_s   = kwargs['vel_dist_law_sec']
        central_vel_s    = kwargs['central_vel_sec'] 
        min_vel_s        = kwargs['min_vel_sec']      
        max_vel_s        = kwargs['max_vel_sec']      
        step_angle_vel_s   = kwargs['step_angle_vel_sec'] 
        high_lat_vel_s   = kwargs['high_lat_vel_sec'] 
        low_lat_vel_s    = kwargs['low_lat_vel_sec']  
        kappa_dist_law_s = kwargs['kappa_dist_law_sec']
        central_op_s     = kwargs['central_op_sec']    
        min_op_s         = kwargs['min_op_sec']        
        max_op_s         = kwargs['max_op_sec']        
        step_angle_op_s    = kwargs['step_angle_op_sec']   
        high_lat_op_s    = kwargs['high_lat_op_sec']  
        low_lat_op_s     = kwargs['low_lat_op_sec']    

        M,V,O = initialize_components.InitializeComponent(mass_dist_law_s,vel_dist_law_s,kappa_dist_law_s)
        ET = thermalization.Thermalization(therm_model)
        r_ph_sec, L_bol_sec = initialize_components.expansion_angular_distribution(M,V,O,ET,exp_model,ang_dist,omega_dist,m_ej_s,time,
          step_angle_mass = step_angle_mass_s,
          high_lat_flag = high_lat_flag_s,
          vel_dist_law  = vel_dist_law_s,
          central_vel   = central_vel_s,
          min_vel       = min_vel_s, 
          max_vel       = max_vel_s,     
          step_angle_vel  = step_angle_vel_s,
          high_lat_vel  = high_lat_vel_s,
          low_lat_vel   = low_lat_vel_s,
          kappa_dist_law= kappa_dist_law_s,
          central_op    = central_op_s,
          min_op        = min_op_s,
          max_op        = max_op_s,
          step_angle_op   = step_angle_op_s, 
          high_lat_op   = high_lat_op_s,
          low_lat_op    = low_lat_op_s)

    else:
        r_ph_sec = np.full((n_ang,n_time),units.small)
        L_bol_sec = np.full((n_ang,n_time),units.small)

    r_ph_tot = np.maximum(np.maximum(r_ph_dyn,r_ph_wind),r_ph_sec)
    L_bol_tot = L_bol_dyn + L_bol_wind + L_bol_sec
  
    tmp = []
    for k in range(n_ang):
        tmp.append(np.array([T_eff_calc(L,omega_dist[k],R) for L,R in zip(L_bol_tot[k,:],r_ph_tot[k,:])]))
        T_eff_tot = np.asarray(tmp)
    
    return r_ph_tot,L_bol_tot,T_eff_tot 


if __name__=="__main__":
    AD = angular_distribution.AngularDistribution("uniform")
    ang_dist, omega_dist = AD(12)

    dyn_flag = False
    wind_flag = True
    sec_flag = False

    r_ph_tot,L_bol_tot,T_eff_tot = lightcurve(dyn_flag,wind_flag,sec_flag,ang_dist,omega_dist,'GK','BKWM',
# mass dynamic ejecta
        mass_dist_law_dyn='sin2',
        m_ej_dyn=0.005,
        step_angle_mass_dyn=None,
        high_lat_flag_dyn=None,
# velocity dynamic ejecta
        vel_dist_law_dyn ='uniform',
        central_vel_dyn  =0.2,
        min_vel_dyn      =None,
        max_vel_dyn      =None,
        step_angle_vel_dyn =None,
        high_lat_vel_dyn =None,
        low_lat_vel_dyn  =None,
# opacity dynamic ejecta
        kappa_dist_law_dyn='step',
        central_op_dyn    =None,
        min_op_dyn        =None,
        max_op_dyn        =None,
        step_angle_op_dyn   =0.33,
        high_lat_op_dyn   =0.1,
        low_lat_op_dyn    =10.,
# mass wind ejecta
        mass_dist_law_wind='step',
        m_ej_wind         =0.01,
        step_angle_mass_wind=1.1,
        high_lat_flag_wind=1,
# velocity wind ejecta
        vel_dist_law_wind ='uniform',
        central_vel_wind  =0.1,
        min_vel_wind      =None,
        max_vel_wind      =None,
        step_angle_vel_wind =None,
        high_lat_vel_wind =None,
        low_lat_vel_wind  =None,
# opacity wind ejecta
        kappa_dist_law_wind='step',
        central_op_wind    =None,
        min_op_wind        =None,
        max_op_wind        =None,
        step_angle_op_wind   =0.9,
        high_lat_op_wind   =0.1,
        low_lat_op_wind    =1.0,
# mass secular ejecta
        mass_dist_law_sec='sin2',
        m_ej_sec         =0.02,
        step_angle_mass_sec=None,
        high_lat_flag_sec=None,
# velocity secular ejecta
        vel_dist_law_sec ='uniform',
        central_vel_sec  =0.06,
        min_vel_sec      =None,
        max_vel_sec      =None,
        step_angle_vel_sec =None,
        high_lat_vel_sec =None,
        low_lat_vel_sec  =None,
# opacity secular ejecta
        kappa_dist_law_sec='uniform',
        central_op_sec    =5.0,
        min_op_sec        =None,
        max_op_sec        =None,
        step_angle_op_sec =None,
        high_lat_op_sec   =None,
        low_lat_op_sec    =None)

    print r_ph_tot

