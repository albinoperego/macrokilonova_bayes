import angular_distribution as ad
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

def smooth_profile(lval):
    val_smooth = []
    for i in range(len(lval)):
        if (i == 0):
            val_smooth.append(10. ** ((lval[0])))
        elif (i == 1):
            val_smooth.append(10. ** ((lval[i - 1] + lval[i] + lval[i + 1]) / 3.))
        elif (i == 2):
            val_smooth.append(10. ** ((lval[i - 2] + lval[i - 1] + lval[i] + lval[i + 1] + lval[i + 2]) / 5.))
        elif (i == len(lval) - 3):
            val_smooth.append(10. ** ((lval[i - 2] + lval[i - 1] + lval[i] + lval[i + 1] + lval[i + 2]) / 5.))
        elif (i == len(lval) - 2):
            val_smooth.append(10. ** ((lval[i - 1] + lval[i] + lval[i + 1]) / 3.))
        elif (i == len(lval) - 1):
            val_smooth.append(10. ** ((lval[i])))
        else:
            val_smooth.append(10. ** ((lval[i - 3] + lval[i - 2] + lval[i - 1] + lval[i] + lval[i + 1] + lval[
                i + 2] + lval[i + 3]) / 7.))
    val_smooth = np.asarray(val_smooth)
    return val_smooth

def importNRprofiles(fname, ang_dist, kappa_low=1., kappa_high=30.):

    th, mass, vel, ye = np.loadtxt(fname, unpack=True, usecols=(0,1,2,3))

    for i in range(len(th)):
        if (mass[i] < 1.e-9):
            mass[i] = 1.e-9
            vel[i] = 0.1
            ye[i] = 0.1

# find the central angles
    th_central = []
    for a in ang_dist:
        th_central.append(0.5*(a[1]+a[0]))

    lmass = np.log10(mass)
    mass_smooth = smooth_profile(lmass)
    tmp1 = np.sum(mass)
    tmp2 = np.sum(mass_smooth)
    mass_smooth = tmp1/tmp2 * mass_smooth
#    print(sum(mass_smooth))

# smooth the ye
    lye = np.log10(ye)
    ye_smooth = smooth_profile(lye)
    tmp1 = np.sum(ye*mass)
    tmp2 = np.sum(ye_smooth*mass)
    ye_smooth = tmp1/tmp2 * ye_smooth

# smooth the velocity
    lmom = np.log10(mass*vel)
    vel_smooth = smooth_profile(lmom)
    tmp1 = np.sum(vel*mass)
    #tmp1 = 1.
    tmp2 = np.sum(vel_smooth)
    #tmp2 = 1.
    vel_smooth = (tmp1/tmp2 * vel_smooth) / mass_smooth


    func_m = interpolate.interp1d(th, mass_smooth)
    func_vel = interpolate.interp1d(th, vel_smooth)
    func_ye = interpolate.interp1d(th, ye_smooth)

    profile_m = func_m(th_central)
    profile_ye = func_ye(th_central)
    profile_vel = func_vel(th_central)

#    print(2.*sum(profile_m))

    profile_m = profile_m * np.sum(mass_smooth)/np.sum(2.*profile_m)

    profile_kappa = np.zeros(len(th_central))

    for i in range(len(th_central)):
        if (profile_ye[i] < 0.25):
            profile_kappa[i] = kappa_high
        else: # (profile_ye[i] > 0.25):
            profile_kappa[i] = kappa_low

#    fig1 = plt.figure()
#    plt.scatter(th,np.log10(mass))
#    plt.scatter(th,np.log10(mass_smooth),color='g')
#    plt.scatter(th_central,np.log10(profile_m),color='r')

#    fig2 = plt.figure()
#    plt.scatter(th,ye)
#    plt.scatter(th,ye_smooth,color='g')
#    plt.scatter(th_central,profile_ye,color='r')

#    fig3 = plt.figure()
#    plt.scatter(th,vel,color='g')
#    plt.scatter(th_central,profile_vel,color='r')

#    plt.show()
#     exit('done')
    return profile_m, profile_vel, profile_kappa
