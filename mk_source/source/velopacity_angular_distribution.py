import numpy as np

class VelopacityAngularDistribution(object):

    def __init__(self,angular_law_distribution):

        if angular_law_distribution=="uniform":
            self.velopacity_angular_distribution = uniform_velopacity_distribution
        elif angular_law_distribution=="sin":
            self.velopacity_angular_distribution = sin_velopacity_distribution
        elif angular_law_distribution=="cos2":
            self.velopacity_angular_distribution = cos2_velopacity_distribution
        elif angular_law_distribution=="sin2":
            self.velopacity_angular_distribution = sin2_velopacity_distribution
        elif angular_law_distribution=="abscos":
            self.velopacity_angular_distribution = abscos_velopacity_distribution
        elif angular_law_distribution=="step":
            self.velopacity_angular_distribution = step_velopacity_distribution
        else:
            print("Unknown velocity or opacity angular distribution")
            exit(0)
        
    def __call__(self, angles_array, **kwargs):
        return self.velopacity_angular_distribution(angles_array, **kwargs)

def uniform_velopacity_distribution(angles_array, **kwargs):
    central_velop = kwargs['central_velop']
    if central_velop is None:
        print("Error! user must specify a velocity or opacity! exiting\n")
        exit(-1)
    return np.array([central_velop for a in angles_array])

def sin_velopacity_distribution(angles_array, **kwargs):
    min_velop = kwargs['min_velop']
    max_velop = kwargs['max_velop']
    if min_velop is None:
        print("Error! user must specify a minimum velocity or opacity! exiting\n")
        exit(-1)
    if max_velop is None:
        print("Error! user must specify a maximum velocity or opacity! exiting\n")
        exit(-1)
    delta_velop = max_velop-min_velop
    return np.array([min_velop + delta_velop*np.sin(0.5*(a[1]+a[0])) for a in angles_array])

def sin2_velopacity_distribution(angles_array, **kwargs):
    min_velop = kwargs['min_velop']
    max_velop = kwargs['max_velop']
    if min_velop is None:
        print("Error! user must specify a minimum velocity or opacity! exiting\n")
        exit(-1)
    if max_velop is None:
        print("Error! user must specify a maximum velocity or opacity! exiting\n")
        exit(-1)
    delta_velop = max_velop-min_velop
    return np.array([min_velop + delta_velop*(np.sin(0.5*(a[1]+a[0]))**2) for a in angles_array])

def cos2_velopacity_distribution(angles_array, **kwargs):
    min_velop = kwargs['min_velop']
    max_velop = kwargs['max_velop']
    if min_velop is None:
        print("Error! user must specify a minimum velocity or opacity! exiting\n")
        exit(-1)
    if max_velop is None:
        print("Error! user must specify a maximum velocity or opacity! exiting\n")
        exit(-1)
    delta_velop = max_velop-min_velop
    return np.array([min_velop + delta_velop*(np.cos(0.5*(a[1]+a[0]))**2) for a in angles_array])

def abscos_velopacity_distribution(angles_array, **kwargs):
    min_velop = kwargs['min_velop']
    max_velop = kwargs['max_velop']
    if min_velop is None:
        print("Error! user must specify a minimum velocity or opacity! exiting\n")
        exit(-1)
    if max_velop is None:
        print("Error! user must specify a maximum velocity or opacity! exiting\n")
        exit(-1)
    delta_velop = max_velop-min_velop
    return np.array([min_velop + delta_velop*(abs(np.cos(0.5*(a[1]+a[0])))) for a in angles_array])

def step_velopacity_distribution(angles_array, **kwargs):
    step_angle   = kwargs['step_angle']
    high_lat_velop = kwargs['high_lat_velop']
    low_lat_velop  = kwargs['low_lat_velop']
    if step_angle is None:
        print("Error! Must specify a step angle in radians! exiting\n")
        exit(-1)
    if high_lat_velop is None:
        print("Error! user must specify a high latitude velocity or opacity! exiting\n")
        exit(-1)
    if low_lat_velop is None:
        print("Error! user must specify a low latitude velocity or opacity! exiting\n")
        exit(-1)
    o = []
    for a in angles_array:
        if np.sin(0.5*(a[1]+a[0])) < np.sin(step_angle):
            o.append(high_lat_velop)
        else:
            o.append(low_lat_velop)
    return np.array(o)

if __name__=="__main__":
    angular_distribution = [(0,1),(1,2),(2,3.1415)]
    V = VelopacityAngularDistribution("step")
    x = V(angular_distribution,low_lat_velop=1.,high_lat_velop=0.001,step_angle=1.0)
    print x
