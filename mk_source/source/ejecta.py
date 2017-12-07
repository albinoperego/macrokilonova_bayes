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

class Ejecta(object):

    def __init__(self, n_shells, names, params, *args, **kwargs):
        assert len(names) == n_shells
        self.components = [Shell(n, params[n], **kwargs) for n in names]

if __name__=="__main__":
    params = {}
    params['wind'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    params['secular'] = {'mass_dist':'uniform', 'vel_dist':'step', 'op_dist':'step', 'therm_model':'BKWM', 'eps_ye_dep':True}
    E = Ejecta(2, params.keys(), params)
    print E.components
