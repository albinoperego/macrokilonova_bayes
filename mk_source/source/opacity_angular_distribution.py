import numpy as np

class OpacityAngularDistribution(object):

    def __init__(self,angular_law_distribution):

        if angular_law_distribution=="uniform":
            self.opacity_angular_distribution = uniform_opacity_distribution
        elif angular_law_distribution=="sin":
            self.opacity_angular_distribution = sin_opacity_distribution
        elif angular_law_distribution=="cos2":
            self.opacity_angular_distribution = cos2_opacity_distribution
        elif angular_law_distribution=="sin2":
            self.opacity_angular_distribution = sin2_opacity_distribution
        elif angular_law_distribution=="abscos":
            self.opacity_angular_distribution = abscos_opacity_distribution
        elif angular_law_distribution=="step":
            self.opacity_angular_distribution = step_opacity_distribution
        else:
            print("Unknown opacity angular distribution")
            exit(0)
        
    def __call__(self, angles_array, **kwargs):
        return self.opacity_angular_distribution(angles_array, **kwargs)

def uniform_opacity_distribution(angles_array, **kwargs):
    central_op = kwargs['central_op']
    if central_op is None:
        print("Error! user must specify a opacity! exiting\n")
        exit(-1)
    return np.array([central_op for a in angles_array])

def sin_opacity_distribution(angles_array, **kwargs):
    min_op = kwargs['min_op']
    max_op = kwargs['max_op']
    if min_op is None:
        print("Error! user must specify a minimum opacity! exiting\n")
        exit(-1)
    if max_op is None:
        print("Error! user must specify a maximum opacity! exiting\n")
        exit(-1)
    delta_op = max_op-min_op
    return np.array([min_op + delta_op*np.sin(0.5*(a[1]+a[0])) for a in angles_array])

def sin2_opacity_distribution(angles_array, **kwargs):
    min_op = kwargs['min_op']
    max_op = kwargs['max_op']
    if min_op is None:
        print("Error! user must specify a minimum opacity! exiting\n")
        exit(-1)
    if max_op is None:
        print("Error! user must specify a maximum opacity! exiting\n")
        exit(-1)
    delta_op = max_op-min_op
    return np.array([min_op + delta_op*(np.sin(0.5*(a[1]+a[0]))**2) for a in angles_array])

def cos2_opacity_distribution(angles_array, **kwargs):
    min_op = kwargs['min_op']
    max_op = kwargs['max_op']
    if min_op is None:
        print("Error! user must specify a minimum opacity! exiting\n")
        exit(-1)
    if max_op is None:
        print("Error! user must specify a maximum opacity! exiting\n")
        exit(-1)
    delta_op = max_op-min_op
    return np.array([min_op + delta_op*(np.cos(0.5*(a[1]+a[0]))**2) for a in angles_array])

def abscos_opacity_distribution(angles_array, **kwargs):
    min_op = kwargs['min_op']
    max_op = kwargs['max_op']
    if min_op is None:
        print("Error! user must specify a minimum opacity! exiting\n")
        exit(-1)
    if max_op is None:
        print("Error! user must specify a maximum opacity! exiting\n")
        exit(-1)
    delta_op = max_op-min_op
    return np.array([min_op + delta_op*(abs(np.cos(0.5*(a[1]+a[0])))) for a in angles_array])

def step_opacity_distribution(angles_array, **kwargs):
    step_angle_op   = kwargs['step_angle_op']
    high_lat_op = kwargs['high_lat_op']
    low_lat_op  = kwargs['low_lat_op']
    if step_angle_op is None:
        print("Error! Must specify a step angle in radians! exiting\n")
        exit(-1)
    if high_lat_op is None:
        print("Error! user must specify a high latitude opacity! exiting\n")
        exit(-1)
    if low_lat_op is None:
        print("Error! user must specify a low latitude opacity! exiting\n")
        exit(-1)
    o = []
    for a in angles_array:
        if np.sin(0.5*(a[1]+a[0])) < np.sin(step_angle_op):
            o.append(high_lat_op)
        else:
            o.append(low_lat_op)
    return np.array(o)

if __name__=="__main__":
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    V = OpacityAngularDistribution("step")
    x = V(angular_distribution,low_lat_op=1.,high_lat_op=0.001,step_angle_op=1.0)
