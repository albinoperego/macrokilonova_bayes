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
from expansion_model_single_spherical import ExpansionModelSingleSpherical
from shell import Shell

def T_eff_calc(Lum,dOmega,r_ph):
    return (Lum/(dOmega*r_ph**2*units.sigma_SB))**(1./4.)

class Ejecta(object):

    def __init__(self, n_shells, names, params, *args, **kwargs):
        assert len(names) == n_shells
        self.components = [Shell(n, params[n], **kwargs) for n in names]
    
    def lightcurve(self,
                                       angular_distribution,
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
                                       **kwargs):
        physical_radius = []
        self.bolometric_luminosity = []
        for c in self.components:
            r, Lb = c.expansion_angular_distribution(angular_distribution,
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
                                                     **kwargs)
            physical_radius.append(r)
            self.bolometric_luminosity.append(Lb)
        r_ph_tot = np.maximum(np.maximum(r_ph_dyn,r_ph_wind),r_ph_sec)
        
        self.physical_radius = physical_radius[0]
        for k in np.arange(1,len(physical_radius)): self.physical_radius = np.maximum(self.physical_radius,r[k])

        self.total_bolometric_luminosity = None
        for b in self.bolometric_luminosity:
            if self.total_bolometric_luminosity is None: self.total_bolometric_luminosity = b
            else: self.total_bolometric_luminosity += b

        tmp = []
        for k in range(len(angular_distribution)):
            tmp.append(np.array([T_eff_calc(L,omega_distribution[k],R) for L,R in zip(self.bolometric_luminosity[k,:],self.physical_radius[k,:])]))
            self.T_eff_tot = np.asarray(tmp)
        return time, self.physical_radius, self.bolometric_luminosity, self.T_eff_tot

if __name__=="__main__":
    params = {}
    params['wind'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    params['secular'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    E = Ejecta(2, params.keys(), params)
    print E.components
