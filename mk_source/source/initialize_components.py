import numpy as np
import mass_angular_distribution
import thermalization
import velocity_angular_distribution
import opacity_angular_distribution
import expansion_model_single_spherical
from scipy import interpolate
import units

def InitializeComponent(mass_dist,vel_dist,op_dist):
    MAD = mass_angular_distribution.MassAngularDistribution(mass_dist)
    VAD = velocity_angular_distribution.VelocityAngularDistribution(vel_dist)
    OAD = opacity_angular_distribution.OpacityAngularDistribution(op_dist)
    return MAD,VAD,OAD

def FillComponent(MM,VV,OO,m_tot,angular_distribution,**kwargs):
    print ('I am filling')
    x = MM(m_tot,angular_distribution,**kwargs)
    print ('x done')
    y = VV(angular_distribution,**kwargs)
    z = OO(angular_distribution,**kwargs)
    return x,y,z

def r_ph_calc(vel,t):
    return (vel*units.c)*t

def calc_eps_nuc(kappa,time,val_min):
  if (kappa > 10.):
    tmp=val_min * ft_eps_nuc(time)
  elif (kappa < 1.):
    tmp=val_min
  else:
    tmp=0.5*(val_min + val_min*ft_eps_nuc(time))
  return tmp 

def ft_eps_nuc(time):
  time_day=time*units.sec2day
  tmp = min(max(4*time_day-4.,-20),20)
  return 0.5 + 2.5/(1.+np.exp(tmp))

def heat_rate(eps_nuc,alpha,eps_th,t,t0,sigma0):
  return eps_nuc*(0.5 - 1./np.pi * np.arctan((t-t0)/sigma0))**alpha * (eps_th/0.5)

def bol_lum(alpha,eps_nuc,eps_th,t,t0,sigma0,m_em):
  eps=heat_rate(eps_nuc,alpha,eps_th,t,t0,sigma0)
  return m_em*eps

def expansion_angular_distribution(MM,VV,OO,ET,model_name,angular_distribution,omega_distribution,m_tot,time,**kwargs):

    M = expansion_model_single_spherical.ExpansionModelSingleSpherical(model_name)
    m_ej_dist,v_rms_dist,kappa_dist = FillComponent(MM,VV,OO,m_tot,angular_distribution,**kwargs)   
    v_min = 0.0001


    r_ph = []
    L_bol = []

    for Omega,m_ej,v_rms,kappa in zip(omega_distribution,m_ej_dist,v_rms_dist,kappa_dist):

        vel,m_vel,t_diff,t_fs = M(Omega,m_ej,v_rms,v_min,kappa)

        print('I am before t_diff extremes')
        print 'Min t_diff',min(t_diff)
        print 'Max t_diff',max(t_diff)


        f_vel_t_diff = interpolate.interp1d(t_diff[::-1],vel[::-1])
        f_vel_t_fs   = interpolate.interp1d(t_fs[::-1],vel[::-1])

        f_m_vel_t_diff = interpolate.interp1d(t_diff[::-1],m_vel[::-1])
        f_m_vel_t_fs   = interpolate.interp1d(t_fs[::-1],m_vel[::-1])
        
        print('Interpolation done')
#        for t in time:
#            print t,f_vel_t_diff(t)

        v_diff  = np.array([f_vel_t_diff(t) for t in time])
        v_fs    = np.array([f_vel_t_fs(t) for t in time])
        mv_diff = np.array([f_m_vel_t_diff(t) for t in time])
        mv_fs   = np.array([f_m_vel_t_fs(t) for t in time])

        print list(v_diff)

#        m_rad = np.array([(mv_diff(t)-mv_fs(t)) for t in time])
        m_rad = mv_diff-mv_fs
    
        r_ph.append(np.array([r_ph_calc(v,t) for v,t in zip(v_fs,time)]))

        e_th = np.array([ET(time_sec=t,mass_ej=m,omega=Omega,vel=v,cnst_eff=0.333) for t,m,omega,v in zip(time,m_ej_dist,omega_distribution,v_rms_dist)])

        eps0 = 2.e+18
        eps_nuc = np.array([calc_eps_nuc(k,t,eps0) for k,t in zip(kappa_dist,time)])

        sigma0 = 0.11
        alpha = 1.3
        t0eps = 1.3
        L_bol.append(np.array([bol_lum(alpha,en,et,t,t0eps,sigma0,m) for en,et,t,m in zip(eps_nuc,e_th,time,m_rad)]))

    return np.asarray(r_ph),np.asarray(L_bol)


if __name__=="__main__":
    MM,VV,OO = InitializeComponent("sin2","step","step")
    ET = thermalization.Thermalization("BKWM")
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    omega_distribution = [0.01,0.2,0.5]
    time_min = 36000.      #
    time_max = 172800.   #
    n_time = 200
    time = np.linspace(time_min,time_max,n_time)
    r_ph, L_bol = expansion_angular_distribution(MM,VV,OO,ET,"GK",angular_distribution,omega_distribution,0.1,time,low_lat_vel=0.2,high_lat_vel=0.001,step_angle_vel=1.0,low_lat_op=0.2,high_lat_op=0.001,step_angle_op=1.0)

    print r_ph

#    x,y,z = FillComponent(1.0,angular_distribution,low_lat_vel=0.2,high_lat_vel=0.001,step_angle_vel=1.0,low_lat_op=0.2,high_lat_op=0.001,step_angle_op=1.0)

#    x = M(1.0,angular_distribution)
#    print(x)
#    y = V(angular_distribution,low_lat_velop=0.2,high_lat_velop=0.001,step_angle=1.0)
#    print(y)
#    z = V(angular_distribution,low_lat_velop=1.,high_lat_velop=0.001,step_angle=1.2)
#    print(z)

