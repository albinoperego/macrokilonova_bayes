import numpy as np

class MassAngularDistribution(object):

    def __init__(self,angular_law_distribution):

        if angular_law_distribution=="uniform":
            self.mass_angular_distribution = uniform_mass_distribution
        elif angular_law_distribution=="sin":
            self.mass_angular_distribution = sin_mass_distribution
        elif angular_law_distribution=="cos2":
            self.mass_angular_distribution = cos2_mass_distribution
        elif angular_law_distribution=="sin2":
            self.mass_angular_distribution = sin2_mass_distribution
        elif angular_law_distribution=="step":
            self.mass_angular_distribution = step_mass_distribution
        else:
            print("Unknown mass angular distribution")
            exit(0)
        
    def __call__(self, m_tot, angles_array, **kwargs):
        return self.mass_angular_distribution(m_tot, angles_array, **kwargs)

def uniform_mass_distribution(m_tot,angles_array, **kwargs):
    return np.array([m_tot * 0.5 * (np.cos(a[0])-np.cos(a[1])) for a in angles_array])

def sin_mass_distribution(m_tot,angles_array, **kwargs):
    return np.array([(m_tot / np.pi) * (a[1]-a[0] - ( np.sin(a[1])*np.cos(a[1]) - np.sin(a[0])*np.cos(a[0]))) for a in angles_array])
                          
def cos2_mass_distribution(m_tot,angles_array, **kwargs):
    return np.array([m_tot * 0.5 * (np.cos(a[0])**3 - np.cos(a[1])**3) for a in angles_array])

def sin2_mass_distribution(m_tot,angles_array, **kwargs):
    return np.array([m_tot * 0.0625 * (np.cos(3.*a[1]) - 9.*np.cos(a[1]) - np.cos(3.*a[0]) + 9.*np.cos(a[0]) ) for a in angles_array])

def step_mass_distribution(m_tot, angles_array, **kwargs):

    step_angle = kwargs['step_angle_mass']
    high_lat_flag = kwargs['high_lat_flag']
    if step_angle is None:
        print("Error! Must specify a step angle in radians! exiting\n")
        exit(-1)

    if high_lat_flag is None:
        print("Error! User must specify high or low latitude side! exiting\n")
        exit(-1)
    elif (high_lat_flag != 1 and high_lat_flag != 0):
        print("Error! User must specify 1 for high latitude and 0 for low latitude! exiting\n")
        exit(-1)    

    o = []
    if (high_lat_flag == 1):
        for a in angles_array:
            if np.sin(0.5*(a[1]+a[0])) < np.sin(step_angle):
                o.append(m_tot * 0.5 / (1.-np.cos(step_angle)) * (np.cos(a[0])- np.cos(a[1])))
            else:
                o.append(np.maximum(m_tot * 1.e-4,1.e-5))
    else:
         for a in angles_array:
            if np.sin(0.5*(a[1]+a[0])) > np.sin(step_angle):
                o.append(m_tot * 0.5 / (1.-np.cos(step_angle)) * (np.cos(a[0])- np.cos(a[1])))
            else:
                o.append(np.maximum(m_tot * 1.e-4,1.e-5))

    return np.array(o)

if __name__=="__main__":
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    M = MassAngularDistribution("step")
    x = M(1.0,angular_distribution, step_angle_mass=1.1,high_lat_flag=0)
