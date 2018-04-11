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
import interpolation_barnes as itp

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

    def calc_Tfloor(self,kappa,T_floor_LA,T_floor_Ni):
        if (kappa < 1.):
            return T_floor_Ni
        elif (kappa > 10.):
            return T_floor_LA
        else:
            return 0.5*(T_floor_Ni + T_floor_LA)

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
        T_floor_Ni = glob_vars['T_floor_Ni']
        T_floor_LA = glob_vars['T_floor_LA']

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


            Ltmp = np.array([self.bolometric_luminosity(m,
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

            rtmp = self.r_ph_calc(v_fs, time)

            Tf = self.calc_Tfloor(kappa,T_floor_LA,T_floor_Ni)

            self.physical_radius.append(np.array([min(r,np.sqrt(L/(4.*np.pi*units.sigma_SB*Tf**4))) for r,L in zip(rtmp,Ltmp)]))

            self.Lbol.append(Ltmp)

        return self.physical_radius, self.Lbol

#======================================================================#

#				VILLAR

# TIME-SCALE
    def t_d(self,omega,k,m_ej,v_ej):
        beta=13.4
        m = (4.*np.pi)*m_ej/omega
        return(np.sqrt((2.*k*m)/(beta*v_ej*units.c)))

# THERMALISATION EFFICIENCY
# Note: If input parameter is outside the interpolation range (at least one of a,b or d is 0) returns thermalisation efficiency=0.5 
    def epsilon(self,t,a,b,d):
        if (a*b*d==0):
            print('\n Mass and velocity outside interpolation range for Barnes a,b,d parameters ... Returning Thermalisation Efficiency = 0.5\n')
            return(0.5)
        else:
            return(0.36*(np.exp(-a*t/86400.)+((np.log(1.+2.*b*((t/86400.)**d)))/(2.*b*((t/86400.)**d)))))
#       return(0.5)

# RADIOACTIVE HEATING RATE
    def L_in(self,t,m_ej):
        t_0=1.3; sigma=0.11
        return(4.e18*m_ej*((0.5-((np.arctan((t-t_0)/(sigma)))/np.pi))**1.3))

# Epsilon_Y
# Factor that introduces heating rate dependence on opacity from Perego et al. (2017)
    e_min=0.5; e_max=2.5; t_e=86400.
    def e_y(self,t,k):
	    if(k>2.):
		    return(1.)
	    else:
		    return(self.e_min+self.e_max/(1.+np.exp(4.*(t/self.t_e-1.))))                

# BOLOMETRIC LUMINOSITY
# Note: to avoid overflow in e^((t/td)^2) a 'cut' is introduced. If (t/td)^2 > 500, it is fixed to 500. Now L is computed also at later times.
    def L(self,omega,t,k,m_ej,v_ej,a,b,d):
        td=self.t_d(omega,k,m_ej,v_ej)
        init_x = np.logspace(-3,np.log10(t[0]),100.)
        init_y = (self.L_in(init_x,m_ej)*self.epsilon(init_x,a,b,d)*self.e_y(init_x,k)*np.exp((init_x**2)/(td**2))*init_x/td)
        init_int=integrate.trapz(init_y,init_x)
        cut = t
        cut = (cut**2)/((td**2)*500)
        cut[t<td*np.sqrt(500)]=1
        integral = integrate.cumtrapz(self.L_in(t,m_ej)*self.epsilon(t,a,b,d)*self.e_y(t,k)*np.exp((t/(td*cut))**2)*(t/td),t,initial=0)
        integral+=init_int
        Lum=integral*np.exp(-(t**2.)/(td**2.))/td
        return(Lum)

    def villar(self,time,omega,m_ej,v_ej,k,T_floor_LA,T_floor_Ni):
		# Interpolation from Barnes et al. (2016) using module 'interpolation_barnes'.
		# Now it uses a single variable [m_ej/(v_ej^2)] instead of using m_ej and v_ej separately.
        a = itp.int_a(m_ej/(v_ej**2))
        b = itp.int_b(m_ej/(v_ej**2))
        d = itp.int_d(m_ej/(v_ej**2))
        m_ej = m_ej * units.Msun
        v_ej = v_ej * units.c
        L_bol=self.L(omega,time,k,m_ej,v_ej,a,b,d)
        T_floor = self.calc_Tfloor(k,T_floor_LA,T_floor_Ni)
        rtmp = v_ej*time
        R_phot=np.array([min(r,np.sqrt(L/(4.*np.pi*units.sigma_SB*T_floor**4))) for r,L in zip(rtmp,L_bol)])
#        T_BB=(L_bol/(4.*np.pi*units.sigma_SB*(R_phot**2.)))**0.25
        return(L_bol,R_phot)#,T_BB)


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
        T_floor_Ni = glob_vars['T_floor_Ni']
        T_floor_LA = glob_vars['T_floor_LA']

        if (shell_vars['m_ej'] == None):
            m_tot = np.float(glob_vars['m_disk']) * np.float(shell_vars['xi_disk'])
        elif (shell_vars['xi_disk'] == None):
            m_tot = np.float(shell_vars['m_ej'])

        self.ejected_mass,self.velocity_rms,self.opacity = self.update(m_tot,angular_distribution,**shell_vars)#,**kwargs)
        self.physical_radius = []
        self.Lbol = []
#        self.Teff = []
        
        for omega,m_ej,v_rms,kappa in zip(omega_distribution,self.ejected_mass,self.velocity_rms,self.opacity):

            #print(omega)
            #print(m_ej)
            #print(v_rms)
            #print(kappa)
            #print(time)

            tmp1,tmp2 = self.villar(time,omega,m_ej,v_rms,kappa,T_floor_LA,T_floor_Ni)

            self.Lbol.append(tmp1)
            self.physical_radius.append(tmp2)

        return self.physical_radius, self.Lbol

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
