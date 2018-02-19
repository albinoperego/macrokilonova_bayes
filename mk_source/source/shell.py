import numpy as np
import angular_distribution
import mass_angular_distribution
import thermalization
import velocity_angular_distribution
import opacity_angular_distribution
import initialize_components
import units
import initialize_components
import observer_projection as op
import nuclear_heat
import scipy.integrate as integrate
from expansion_model_single_spherical import ExpansionModelSingleSpherical

def T_eff_calc(Lum,dOmega,r_ph):
    return (Lum/(dOmega*r_ph**2*units.sigma_SB))**(1./4.)

class Shell(object):
    """
    Ejecta shell class
    Required arguments are:
    name,          # name of the shell     (dynamical, wind, secular)
    mass_dist,     # mass distribution     (step, uniform, continuos)
    vel_dist,      # velocity distribution (step, uniform, continuos)
    op_dist,       # opacity distribution  (step, uniform, continuos)
    therm_model,   # thermal model         (BKWM, cnst)
    eps_ye_dep,    # nuclear heating model (Ye depedence True, False)
    """
    def __init__(self, name, params, **kwargs):
        assert len(params) >= 5
        self.name            = name
        self.ET              = thermalization.Thermalization(params['therm_model'])
        self.EPSNUC          = nuclear_heat.NuclearHeat(params['eps_ye_dep'])
        self.mass_dist       = mass_angular_distribution.MassAngularDistribution(params['mass_dist'])
        self.vel_dist        = velocity_angular_distribution.VelocityAngularDistribution(params['vel_dist'])
        self.op_dist         = opacity_angular_distribution.OpacityAngularDistribution(params['op_dist'])
        self.expansion_model = ExpansionModelSingleSpherical('GK')
    
    def update(self, m_tot, angular_distribution,**kwargs):
        x = self.mass_dist(m_tot,angular_distribution,**kwargs)
        y = self.vel_dist(angular_distribution,**kwargs)
        z = self.op_dist(angular_distribution,**kwargs)
        return x, y, z

    def expansion_angular_distribution(self,
                                       angular_distribution,
                                       omega_distribution,
                                       time,
                                       shell_vars,
                                       glob_vars,
                                       glob_params,
                                       **kwargs):
                                       
# assign the global model variables
        v_min    = glob_params['v_min']
        n_v      = glob_params['n_v']
        vscale   = glob_params['vscale']
        sigma0   = glob_params['sigma0']
        alpha    = glob_params['alpha']
        t0eps    = glob_params['t0eps']
        cnst_eff = glob_params['cnst_eff']

# assign the global variables
        eps0      = glob_vars['eps0']
        a_eps_nuc = glob_vars['a_eps_nuc']
        b_eps_nuc = glob_vars['b_eps_nuc']
        t_eps_nuc = glob_vars['t_eps_nuc']

        if (shell_vars['m_ej'] == None):
            m_tot = np.float(glob_vars['m_disk']) * np.float(shell_vars['xi_disk'])
        elif (shell_vars['xi_disk'] == None):
            m_tot = np.float(shell_vars['m_ej'])

        self.ejected_mass,self.velocity_rms,self.opacity = self.update(m_tot,angular_distribution,**shell_vars)#,**kwargs)
        self.physical_radius = []
        self.Lbol = []
        
        for omega,m_ej,v_rms,kappa in zip(omega_distribution,self.ejected_mass,self.velocity_rms,self.opacity):

            vel,m_vel,t_diff,t_fs = self.expansion_model(omega,m_ej,v_rms,v_min,n_v,vscale,kappa)
            v_diff  = np.interp(time, t_diff[::-1], vel[::-1])
            v_fs    = np.interp(time, t_fs[::-1], vel[::-1])
            mv_diff = np.interp(time, t_diff[::-1], m_vel[::-1])
            mv_fs   = np.interp(time, t_fs[::-1], m_vel[::-1])
            m_rad = mv_diff-mv_fs

            self.physical_radius.append(self.r_ph_calc(v_fs, time))

            self.Lbol.append([self.bolometric_luminosity(m,
                                                    t,
                                                    t0eps,
                                                    sigma0,
                                                    eps0,
                                                    a_eps_nuc,
                                                    b_eps_nuc,
                                                    t_eps_nuc,
                                                    m_ej,
                                                    omega,
                                                    v_rms,
                                                    kappa,
                                                    cnst_eff,
                                                    alpha) for t,m in zip(time,m_rad)])
        self.physical_radius = np.array(self.physical_radius)
        self.Lbol = np.array(self.Lbol)
        tmp = []

        for k in range(len(angular_distribution)):
            tmp.append(np.array([T_eff_calc(L,omega_distribution[k],R) for L,R in zip(self.Lbol[k,:],self.physical_radius[k,:])]))
            self.Teff = np.asarray(tmp)
        return self.physical_radius, self.Lbol, self.Teff

#=======================================================================

# TIME-SCALE
    def t_d(self,omega,k,m_ej,v_ej):
        beta=13.4
        m = (4.*np.pi)*m_ej/omega
        return(np.sqrt((2.*k*m)/(beta*v_ej*units.c)))

# THERMALISATION EFFICIENCY
#    a=ipl.interpolation_a(m_ej/Msun,v_ej/c)
#    b=ipl.interpolation_b(m_ej/Msun,v_ej/c)
#    d=ipl.interpolation_d(m_ej/Msun,v_ej/c)
    def epsilon(self,t,a,b,d):
#                    return(0.36*(math.exp(-a*t/86400.)+((math.log(1.+2.*b*((t/86400.)**d)))/(2.*b*((t/86400.)**d)))))
        return(0.5)

# RADIOACTIVE HEATING RATE
    def L_in(self,t,m_ej):
        t_0=1.3
        sigma=0.11
        return(4.e18*m_ej*((0.5-((np.arctan((t-t_0)/(sigma)))/np.pi))**1.3))

# BOLOMETRIC LUMINOSITY
    def L(self,omega,t,k,m_ej,v_ej,a,b,d):
        x = np.linspace(0.001,t,500)
        y = []
        for i in range(0,len(x)):
            y.append((self.L_in(x[i],m_ej)*self.epsilon(x[i],a,b,d)*(np.exp((x[i]**2.)/((self.t_d(omega,k,m_ej,v_ej))**2.)))*(x[i]/self.t_d(omega,k,m_ej,v_ej))))
        x=np.asarray(x)
        y=np.asarray(y)
        val1=integrate.trapz(y,x)
        return(val1*np.exp(-(t**2.)/((self.t_d(omega,k,m_ej,v_ej))**2.))/self.t_d(omega,k,m_ej,v_ej))

    def villar(self,time,omega,m_ej,v_ej,k):

        m_ej = m_ej * units.Msun
        v_ej = v_ej * units.c

        a = 0.1
        b = 0.1
        d = 0.1

        L_bol,R_phot,T_BB=[],[],[]
        for i in range(0,len(time)):
            L_bol.append(self.L(omega,time[i],k,m_ej,v_ej,a,b,d))

# BLACKBODY TEMPERATURE AND PHOTOSPHERE RADIUS
            R_phot.append(v_ej*time[i])
            T_BB.append((L_bol[i]/(omega*units.sigma_SB*(R_phot[i]**2.)))**0.25)

        return(L_bol,R_phot,T_BB)


    def expansion_angular_distribution_villar(self,
                                              angular_distribution,
                                              omega_distribution,
                                              time,
                                              shell_vars,
                                              glob_vars,
                                              glob_params,
                                              **kwargs):
                                       
# assign the global model variables
        v_min    = glob_params['v_min']
        n_v      = glob_params['n_v']
        vscale   = glob_params['vscale']
        sigma0   = glob_params['sigma0']
        alpha    = glob_params['alpha']
        t0eps    = glob_params['t0eps']
        cnst_eff = glob_params['cnst_eff']

# assign the global variables
        eps0      = glob_vars['eps0']
        a_eps_nuc = glob_vars['a_eps_nuc']
        b_eps_nuc = glob_vars['b_eps_nuc']
        t_eps_nuc = glob_vars['t_eps_nuc']

        if (shell_vars['m_ej'] == None):
            m_tot = np.float(glob_vars['m_disk']) * np.float(shell_vars['xi_disk'])
        elif (shell_vars['xi_disk'] == None):
            m_tot = np.float(shell_vars['m_ej'])

        self.ejected_mass,self.velocity_rms,self.opacity = self.update(m_tot,angular_distribution,**shell_vars)#,**kwargs)
        self.physical_radius = []
        self.Lbol = []
        self.Teff = []
        
        for omega,m_ej,v_rms,kappa in zip(omega_distribution,self.ejected_mass,self.velocity_rms,self.opacity):

            print(omega)
            print(m_ej)
            print(v_rms)
            print(kappa)
            print(time)

            tmp1,tmp2,tmp3 = self.villar(time,omega,m_ej,v_rms,kappa)

            self.Lbol.append(tmp1)
            self.physical_radius.append(tmp2)
            self.Teff.append(tmp3)

        return self.physical_radius, self.Lbol, self.Teff

#=======================================================================

    def bolometric_luminosity(self, m_rad, time,
                              t0eps,sigma0,eps0,
                              a_eps_nuc,b_eps_nuc,t_eps_nuc,
                              m_ej,
                              omega,v_rms,kappa,
                              cnst_eff,alpha):
        
        eps = self.EPSNUC(alpha,time,t0eps,sigma0,
                          eps0,self.ET,m_ej,omega,v_rms,
                          cnst_eff,
                          cnst_a_eps_nuc=a_eps_nuc,
                          cnst_b_eps_nuc=b_eps_nuc,
                          cnst_t_eps_nuc=t_eps_nuc,
                          opacity=kappa)
        return m_rad*eps

    def r_ph_calc(self, vel, t):
        return (vel*units.c)*t




#####################

if __name__=="__main__":
    glob_params = {'v_min':1.e-7,
               'n_v':100,
               'vscale':'linear',
               'sigma0':0.11,
               'alpha':1.3,
               't0eps':1.3,
               'cnst_eff':0.3333}

    glob_vars = {'m_disk':0.1,
             'eps0':1.2e19,
             'a_eps_nuc':0.5,
             'b_eps_nuc':2.5,
             't_eps_nuc':1.0}
    params = {}
    params['dynamics']    = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    shell_vars={}

    shell_vars['dynamics'] = {'xi_disk':None,
                          'm_ej':0.005,
                          'central_vel':0.2,
                          'low_lat_vel':None,
                          'high_lat_vel':None,
                          'step_angle_vel':None,
                          'low_lat_op':10.,
                          'high_lat_op':0.1,
                          'step_angle_op':np.radians(45.)}
    s = Shell('dynamics',params['dynamics'])
    import matplotlib.pyplot as plt
    for k in range(r_ph.shape[0]): plt.plot(r_ph[k, :], L_bol[k, :],'.')
    plt.show()
