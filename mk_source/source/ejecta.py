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
import filters as ft
import observer_projection as op

from expansion_model_single_spherical import ExpansionModelSingleSpherical
from shell import Shell

def T_eff_calc(Lum,dOmega,r_ph):
    return (Lum/(dOmega*r_ph**2*units.sigma_SB))**(1./4.)

class Ejecta(object):

    def __init__(self, n_shells, names, params, *args, **kwargs):
        assert len(names) == n_shells
        self.ncomponents = n_shells
        self.components = [Shell(n, params[n], **kwargs) for n in names]
    
    def lightcurve(self,
                   angular_distribution,
                   omega_distribution,
                   time,
                   shell_vars,
                   glob_vars,
                   glob_params,
                   **kwargs):

        model_grossman = True
        model_villar = False

        photospheric_radii = []
        bolometric_luminosities = []
        for c in self.components:
          
            if (model_grossman):
                r, Lb = c.expansion_angular_distribution(angular_distribution,
                                                    omega_distribution,
                                                     time,
                                                     shell_vars[c.name],
                                                     glob_vars,
                                                     glob_params,
                                                     **kwargs)
            elif (model_villar):    
                r, Lb = c.expansion_angular_distribution_villar(angular_distribution,
                                                     omega_distribution,
                                                     time,
                                                     shell_vars[c.name],
                                                     glob_vars,
                                                     glob_params,
                                                     **kwargs)
#                print('\n Using VILLAR model\n')
            else:
                print('Please choose an available method')
                exit()

            photospheric_radii.append(r)
            bolometric_luminosities.append(Lb)
        
        # select the photospheric radius as the maximum between the different single
        # photospheric radii
        self.photospheric_radius = np.zeros(np.shape(r))
        for i in range(len(self.components)):
            self.photospheric_radius = np.maximum(self.photospheric_radius,photospheric_radii[i])

        # define the total bolometric luminosity as the sum of the different single luminosities
        self.bolometric_luminosity = np.sum(bolometric_luminosities, axis = 0)

        # compute the effective BB temperature based on the photospheric radius and luminosity
        self.T_eff = []
        for k in range(len(angular_distribution)):
            self.T_eff.append(np.array([T_eff_calc(L,omega_distribution[k],R) for L,R in zip(self.bolometric_luminosity[k,:],self.photospheric_radius[k,:])]))

        return np.array(self.photospheric_radius), np.array(self.bolometric_luminosity), np.asarray(self.T_eff)

if __name__=="__main__":
    params = {}
    params['wind'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    params['secular'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    params['dynamical'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}

    E = Ejecta(3, params.keys(), params)
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    omega_distribution = [0.01,0.2,0.5]
    time_min = 36000.      #
    time_max = 172800.   #
    n_time = 2000
    m_tot = 0.1
    v_min = 1e-7
    n_v = 100
    vscale = 'linear'
    eps0 = 1e19
    sigma0 = 1.0
    alpha = 0.1
    t0eps = 1.0
    cnst_eff = 1.0
    a_eps_nuc = 1.0
    b_eps_nuc = 1.0
    t_eps_nuc = 1.0
    time = np.linspace(time_min,time_max,n_time)
    r_ph, L_bol, Teff = E.lightcurve(angular_distribution,
                               omega_distribution,
                               m_tot,
                               time,
                               v_min,
                               n_v,
                               vscale,
                               eps0,
                               sigma0,
                               alpha,
                               t0eps,
                               cnst_eff,
                               a_eps_nuc,
                               b_eps_nuc,
                               t_eps_nuc,
                               low_lat_vel=0.2,
                               high_lat_vel=0.001,
                               step_angle_vel=1.0,
                               low_lat_op=0.2,
                               high_lat_op=0.001,
                               step_angle_op=1.0)
    import matplotlib.pyplot as plt
    for j in range(E.ncomponents):
        for k in range(L_bol.shape[1]): plt.plot(time, L_bol[j, 0, :],'.', label=E.components[j].name)
    plt.legend()
    plt.show()
