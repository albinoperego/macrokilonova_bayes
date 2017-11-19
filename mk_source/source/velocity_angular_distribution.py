import numpy as np

class VelocityAngularDistribution(object):

    def __init__(self,angular_law_distribution):

        if angular_law_distribution=="uniform":
            self.velocity_angular_distribution = uniform_velocity_distribution
        elif angular_law_distribution=="sin":
            self.velocity_angular_distribution = sin_velocity_distribution
        elif angular_law_distribution=="cos2":
            self.velocity_angular_distribution = cos2_velocity_distribution
        elif angular_law_distribution=="sin2":
            self.velocity_angular_distribution = sin2_velocity_distribution
        elif angular_law_distribution=="abscos":
            self.velocity_angular_distribution = abscos_velocity_distribution
        elif angular_law_distribution=="step":
            self.velocity_angular_distribution = step_velocity_distribution
        else:
            print("Unknown velocity angular distribution")
            exit(0)
        
    def __call__(self, angles_array, **kwargs):
        return self.velocity_angular_distribution(angles_array, **kwargs)

def uniform_velocity_distribution(angles_array, **kwargs):
    central_vel = kwargs['central_vel']
    if central_vel is None:
        print("Error! user must specify a velocity ! exiting\n")
        exit(-1)
    return np.array([central_vel for a in angles_array])

def sin_velocity_distribution(angles_array, **kwargs):
    min_vel = kwargs['min_vel']
    max_vel = kwargs['max_vel']
    if min_vel is None:
        print("Error! user must specify a minimum velocity! exiting\n")
        exit(-1)
    if max_vel is None:
        print("Error! user must specify a maximum velocity! exiting\n")
        exit(-1)
    delta_vel = max_vel-min_vel
    return np.array([min_vel + delta_vel*np.sin(0.5*(a[1]+a[0])) for a in angles_array])

def sin2_velocity_distribution(angles_array, **kwargs):
    min_vel = kwargs['min_vel']
    max_vel = kwargs['max_vel']
    if min_vel is None:
        print("Error! user must specify a minimum velocity! exiting\n")
        exit(-1)
    if max_vel is None:
        print("Error! user must specify a maximum velocity! exiting\n")
        exit(-1)
    delta_vel = max_vel-min_vel
    return np.array([min_vel + delta_vel*(np.sin(0.5*(a[1]+a[0]))**2) for a in angles_array])

def cos2_velocity_distribution(angles_array, **kwargs):
    min_vel = kwargs['min_vel']
    max_vel = kwargs['max_vel']
    if min_vel is None:
        print("Error! user must specify a minimum velocity! exiting\n")
        exit(-1)
    if max_vel is None:
        print("Error! user must specify a maximum velocity! exiting\n")
        exit(-1)
    delta_vel = max_vel-min_vel
    return np.array([min_vel + delta_vel*(np.cos(0.5*(a[1]+a[0]))**2) for a in angles_array])

def abscos_velocity_distribution(angles_array, **kwargs):
    min_vel = kwargs['min_vel']
    max_vel = kwargs['max_vel']
    if min_vel is None:
        print("Error! user must specify a minimum velocity! exiting\n")
        exit(-1)
    if max_vel is None:
        print("Error! user must specify a maximum velocity! exiting\n")
        exit(-1)
    delta_vel = max_vel-min_vel
    return np.array([min_vel + delta_vel*(abs(np.cos(0.5*(a[1]+a[0])))) for a in angles_array])

def step_velocity_distribution(angles_array, **kwargs):
    step_angle_vel   = kwargs['step_angle_vel']
    high_lat_vel = kwargs['high_lat_vel']
    low_lat_vel  = kwargs['low_lat_vel']
    if step_angle_vel is None:
        print("Error! Must specify a step angle in radians! exiting\n")
        exit(-1)
    if high_lat_vel is None:
        print("Error! user must specify a high latitude velocity! exiting\n")
        exit(-1)
    if low_lat_vel is None:
        print("Error! user must specify a low latitude velocity! exiting\n")
        exit(-1)
    o = []
    for a in angles_array:
        if np.sin(0.5*(a[1]+a[0])) < np.sin(step_angle_vel):
            o.append(high_lat_vel)
        else:
            o.append(low_lat_vel)
    return np.array(o)

if __name__=="__main__":
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    V = VelocityAngularDistribution("step")
    x = V(angular_distribution,low_lat_vel=1.,high_lat_vel=0.001,step_angle_vel=1.0)
    print x
