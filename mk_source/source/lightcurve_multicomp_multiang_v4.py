import numpy as np
from scipy import interpolate
import sys
#import matplotlib.pyplot as plt

##############################################################
##############################################################

# constants

##############################################################
##############################################################
c        = 2.99792458e10     #[cm/s]
Msun     = 1.98855e33        #[g]
sec2day  = 1./(60.*60.*24.)    #[day/s]
sigma_SB = 5.6704e-5         #[erg/cm^2/s/K^4]
h        = 6.6260755e-27     #[erg*s]
kB       = 1.380658e-16      #[erg/K]
pc2cm    = 3.085678e+18      #[cm/pc]
sec2hour = 1./(60.*60.)
day2hour = 24.

##############################################################
##############################################################

# parameters

##############################################################
##############################################################

m_ej_dyn_tot_sample = [5.e-4,5.e-3,0.01,0.02,0.05]  # PARAMETER [solar masses]
m_disk_sample = [0.01,0.08,0.10,0.12,0.15,0.20]   # PARAMETER
xi_w_sample = [0.001,0.05,0.10,0.15,0.20]    # PARAMETER 
xi_vis_sample = [0.001,0.10,0.20,0.30,0.40]    # PARAMETER 


ye_dist_dyn   = str(sys.argv[1])            #'step45'  # PARAMETER
kappa_dyn_min = float(sys.argv[2])  #1.    
ye_dist_w     = str(sys.argv[3])     #'step45'  # PARAMETER
kappa_min_w   = float(sys.argv[4])   #1.    # PARAMETER
kappa_max_w   = float(sys.argv[5])  #5.    # PARAMETER
kappa_min_vis = float(sys.argv[6])  #5.  # PARAMETER

v_dyn_max     = 0.3  # PARAMETER
v_wind_max    = 0.08  # PARAMETER
v_vis_max     = 0.05  # PARAMETER
eps0   = 2.0e18  # PARAMETER [erg/g/s]


v_dyn_dist    ='unif_mean'  # PARAMETER
kappa_dyn_max = 30.  #10.   # PARAMETER

s_dist_pc = 40e+6    # PARAMETER [pc]

ang_deg_sample=[0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0]     # PARAMETER 
#ang_deg_sample=[50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0]     # PARAMETER 

n_v    = 1000
n_time = 100


##############################################################
##############################################################

# filters

##############################################################
##############################################################

# photometric filters
n_filters=13
filter_name=[" " for x in range(n_filters)]
filter_center   =np.zeros(n_filters)
filter_center_nu=np.zeros(n_filters)

filter_name[0]  ='B'
filter_center[0]=445.   #[nm]

filter_name[1]  ='V'
filter_center[1]=551.   #[nm]

filter_name[2]  ='R'
filter_center[2]=658.   #[nm]

filter_name[3]  ='J'
filter_center[3]=1220.  #[nm]

filter_name[4]  ='Ks'
filter_center[4]=2150.   #[nm]

filter_name[5]  ='g'
filter_center[5]=520.   #[nm]

filter_name[6]  ='r'
filter_center[6]=670.   #[nm]

filter_name[7]  ='i'
filter_center[7]=790.   #[nm]

filter_name[8]  ='z'
filter_center[8]=910.   #[nm]

filter_name[9]  ='H'
filter_center[9]=1630.   #[nm]

filter_name[10]  ='U'
filter_center[10]=365.   #[nm]

filter_name[11]  ='I'
filter_center[11]=806.   #[nm]

filter_name[12]  ='Y'
filter_center[12]=1020.   #[nm]

# convert the filter to frequences
filter_center_nu = c/(1.e-7*filter_center)

filter2use = np.zeros(n_filters)
filter2chi = np.zeros(n_filters)

filter2use[0] = 0.  #1  #B
filter2use[1] = 0.  #1  #V
filter2use[2] = 0.      #R
filter2use[3] = 0.      #J
filter2use[4] = 1.      #K
filter2use[5] = 1.      #g
filter2use[6] = 1.      #r
filter2use[7] = 1.      #i
filter2use[8] = 1.      #z
filter2use[9] = 0.      #H
filter2use[10] = 0.     #U
filter2use[11] = 0.     #I
filter2use[12] = 0.     #Y

filter2chi[0] = 1.  #1  #B
filter2chi[1] = 1.  #1  #V
filter2chi[2] = 1.      #R
filter2chi[3] = 1.      #J
filter2chi[4] = 1.      #K
filter2chi[5] = 1.      #g
filter2chi[6] = 1.      #r
filter2chi[7] = 1.      #i
filter2chi[8] = 1.      #z
filter2chi[9] = 0.      #H
filter2chi[10] = 0.     #U
filter2chi[11] = 1.     #I
filter2chi[12] = 1.     #Y

num_filter2use = sum(1 for i in range(len(filter2use)) if filter2use[i] > 0.1)

t_band = []
m_band = []

# B band
t_band.append([4.5,6.5,8.5])
m_band.append([22.729,23.807,24.581])

# V band
t_band.append([4.5,6.5,8.5,10.5])
m_band.append([21.082,22.361,23.152,23.761])

# R band
t_band.append([4.5   ,6.5   ,8.5   ,10.5  ])
m_band.append([20.238,21.268,22.501,23.379])

# J band
t_band.append([0.48,0.51,1.46,2.46,3.46,4.46,6.47,7.46,8.45,9.45,10.46,11.46])
m_band.append([17.88,17.82,17.45,17.66,17.86,18.08,18.74,19.07,19.69,20.06,20.94,21.16])

# Ks band
t_band.append([0.47,0.50,1.32,1.46,2.45,3.45,4.45,6.46,7.45,8.45,9.45,10.45,12.46,14.46,17.45])
m_band.append([18.62,18.64,17.86,17.77,17.67,17.54,17.6,17.84,17.95,18.25,18.49,18.74,19.34,20.02,20.77 ])

# g band
t_band.append([0.55,1.5,2.5,3.5,5.5])
m_band.append([17.32,18.553,20.314,20.940,22.4])

# r band
t_band.append([0.5  ,1.5   ,2.5   ,3.5   ,5.5  ,6.5   ])
m_band.append([17.14,17.813,19.178,19.740,20.74,21.310])

# i band
t_band.append([0.55,1.5,2.5,6.5])
m_band.append([16.984,17.664,18.354,20.329])

# z band
t_band.append([0.55,1.5,3.5 ,4.5,5.5,6.5,8.5,10.5,11.5,14.5])
m_band.append([16.854,17.62,18.3,18.927,19.160,19.627,20.606,22.014 ,22.815,23.342])

# H band
t_band.append([])
m_band.append([])

# U band
t_band.append([])
m_band.append([])

# Y band
t_band.append([0.6 ,1.5 ,2.56,3.59,7.56])
m_band.append([17.5,17.2,17.5,17.8,19.2])

# I band
t_band.append([4.5,6.5,8.5,10.5,11.5])
m_band.append([19.284,20.142,21.133,22.052,23.])

mpeak             = np.zeros(n_filters)
mpeak_upper_limit = np.zeros(n_filters)
mpeak_lower_limit = np.zeros(n_filters)

#print(m_band[0][:])

for j in range(n_filters):
  try: 
    mpeak[j] = min(m_band[j][:])
  except ValueError:
    mpeak[j] = -1.

mpeak_lower_limit = mpeak - 0.5
mpeak_upper_limit = mpeak + 0.5

##############################################################
##############################################################

# source-observer geometry

##############################################################
##############################################################

# distance of the source
s_dist_cm = s_dist_pc*pc2cm   #[cm]

# define the observer location
#ang_rad=ang_deg*np.pi/180.

##############################################################
##############################################################

# model geometry

##############################################################
##############################################################

# define the rays
n_k = 12
th_low      = np.zeros(n_k)
th          = np.zeros(n_k)
th_high     = np.zeros(n_k)
solid_angle = np.zeros(n_k)
sth         = np.zeros(n_k)
cth         = np.zeros(n_k)
sth2        = np.zeros(n_k)
cth2        = np.zeros(n_k)

th_low[0]=0.
delta=np.pi/float(n_k)
th[0]=0.5*delta
th_high[0]=delta
solid_angle[0] = 2.*np.pi*(np.cos(th_low[0])-np.cos(th_high[0]))
for k in range(1,n_k):
  th_low[k]      = th_low[k-1] +delta
  th[k]          = th[k-1]     +delta
  th_high[k]     = th_high[k-1]+delta
  solid_angle[k] = 2.*np.pi*(np.cos(th_low[k])-np.cos(th_high[k]))
for k in range(n_k):
  sth[k]=np.sin(th[k])
  cth[k]=np.cos(th[k])
  sth2[k]=sth[k]*sth[k]
  cth2[k]=cth[k]*cth[k]


##############################################################
##############################################################

# ejecta properties

##############################################################
##############################################################

def mass_ang_dist(n_rays,m_tot,law,m_rays):
  if (law=='unif'):
    for k in range(n_rays):
      m_rays[k]=m_tot * 0.5 * (np.cos(th_low[k])- np.cos(th_high[k]))
  elif (law=='sin'):
    for k in range(n_rays):
      m_rays[k]=m_tot / np.pi * (th_high[k]-th_low[k] - ( np.sin(th_high[k])*np.cos(th_high[k]) - np.sin(th_low[k])*np.cos(th_low[k])))
  elif (law=='cos2'):
    for k in range(n_rays):
      m_rays[k]=m_tot * 0.5 * (np.cos(th_low[k])**3 - np.cos(th_high[k])**3 )
  elif (law=='sin2'):
    for k in range(n_rays):
      m_rays[k]=m_tot * 0.0625 * (np.cos(3.*th_high[k]) - 9.*np.cos(th_high[k]) - np.cos(3.*th_low[k]) + 9.*np.cos(th_low[k]) )
  elif (law=='step60'):
    for k in range(n_rays):
      if (th[k] < np.pi/3. or th[k] > 2.*np.pi/3.):
        m_rays[k]=m_tot * 0.5 / (1.-np.cos(np.pi/3.)) * (np.cos(th_low[k])- np.cos(th_high[k]))
      else:
        m_rays[k]=max(m_tot * 1.e-4,1.e-5)
  elif (law=='step45'):
    for k in range(n_rays):
      if (th[k] < np.pi/4. or th[k] > 3.*np.pi/4.):
        m_rays[k]=m_tot * 0.5 / (1.-np.cos(np.pi/4.)) * (np.cos(th_low[k])- np.cos(th_high[k]))
      else:
        m_rays[k]=max(m_tot * 1.e-4,1.e-5)
  elif (law=='step30'):
    for k in range(n_rays):
      if (th[k] < np.pi/6. or th[k] > 5.*np.pi/6.):
        m_rays[k]=m_tot * 0.5 / (1.-np.cos(np.pi/6.)) * (np.cos(th_low[k])- np.cos(th_high[k]))
      else:
        m_rays[k]=max(m_tot * 1.e-4,1.e-5)
  else:
    print('Wrong law selection!')
    sys.exit()
  return 0


def vel_dist(n_rays,scale,law,g_min,g_max,g_rays):
  f_th=np.zeros(n_k)
  if (law=='sin'):
    for k in range(n_k):
      f_th[k]=sth[k]
  elif (law=='abscos'):
    for k in range(n_k):
      f_th[k]=abs(cth[k])
  elif (law=='sin2'):
    for k in range(n_k):
      f_th[k]=sth2[k]
  elif (law=='cos2'):
    for k in range(n_k):
      f_th[k]=cth2[k] 
  elif (law=='unif_max'):
    f_th[:]=1.
  elif (law=='unif_min'):
    f_th[:]=0.
  elif (law=='unif_mean'):
    f_th[:]=0.5
  elif (law=='step60'):
    for k in range(n_rays):
      if (th[k] < np.pi/3. or th[k] > 2.*np.pi/3.):
        f_th[k]=1.
      else:
        f_th[k]=0.
  elif (law=='step45'):
    for k in range(n_rays):
      if (th[k] < np.pi/4. or th[k] > 3.*np.pi/4.):
        f_th[k]=1.
      else:
        f_th[k]=0.
  elif (law=='step30'):
    for k in range(n_rays):
      if (th[k] < np.pi/6. or th[k] > 5.*np.pi/6.):
        f_th[k]=1.
      else:
        f_th[k]=0.
  elif (law=='invstep60'):
    for k in range(n_rays):
      if (th[k] < np.pi/3. or th[k] > 2.*np.pi/3.):
        f_th[k]=0.
      else:
        f_th[k]=1.
  elif (law=='invstep45'):
    for k in range(n_rays):
      if (th[k] < np.pi/4. or th[k] > 3.*np.pi/4.):
        f_th[k]=0.
      else:
        f_th[k]=1.
  elif (law=='invstep30'):
    for k in range(n_rays):
      if (th[k] < np.pi/6. or th[k] > 5.*np.pi/6.):
        f_th[k]=0.
      else:
        f_th[k]=1.
  else:
    print('Wrong law selection!')
    sys.exit()

  if (scale=='log'):
    h_min = np.log10(g_min)
    h_max = np.log10(g_max)
    for k in range(n_rays):
      g_rays[k] = 10.**(h_min + (h_max-h_min)*f_th[k])
  elif (scale=='lin'):
    for k in range(n_rays):
      g_rays[k] = g_min + (g_max-g_min)*f_th[k]
  else:
    print('Wrong scale selection!')
    sys.exit()

  return 0

def ye_dist(n_rays,law,unif_value,s_ray):
  if (law=='step60'):
    for k in range(n_rays):
      if (th[k] < np.pi/3. or th[k] > 2.*np.pi/3.):
        s_ray[k]=1.
      else:
        s_ray[k]=0.
  elif (law=='step45'):
    for k in range(n_rays):
      if (th[k] < np.pi/4. or th[k] > 3.*np.pi/4.):
        s_ray[k]=1.
      else:
        s_ray[k]=0.
  elif (law=='step30'):
    for k in range(n_rays):
      if (th[k] < np.pi/6. or th[k] > 5.*np.pi/6.):
        s_ray[k]=1.
      else:
        s_ray[k]=0.
  elif (law=='unif'):
    for k in range(n_rays):
      s_ray[k]=unif_value
  else:
    print('Wrong law selection!')
    sys.exit()
  return 0

def kappa_dist(n_rays,ye,val_min,val_max,k_rays):
  for k in range(n_rays):
    if (ye[k]==1.):
      k_rays[k]=val_min
    elif (ye[k]==0.):
      k_rays[k]=val_max
    elif (ye[k]==0.5):
      k_rays[k]=0.5*(val_max+val_min)

#############################################################
#############################################################

# nuclear heating parameters

##############################################################
##############################################################

t0eps  = 1.3   #[s]
sigma0 = 0.11  #[s]
alpha  = 1.3

#############################################################
#############################################################

# define useful functions

#############################################################
#############################################################

def mass_gt_v(v,mej,v_exp,omega):
  return mej * Msun * ( 1.0 + func_vel(v/v_exp) )  #[g]

def func_vel(x):
  return 35.*x**7/112. - 105.*x**5/80. + 35.*x**3/16.- 35.*x/16.

def t_diff_v(kappa,v,m_v,omega):
  return np.sqrt(kappa*m_v/(omega*v*c*c))  #[s]

def t_fs_v(kappa,v,m_v,omega):
  return np.sqrt(1.5*kappa*m_v/(omega*v*v*c*c))  #[s]

def bol_lum(alpha,eps_nuc,eps_th,t,t0,m_em):
  eps=heat_rate(eps_nuc,alpha,eps_th,t,t0)
  return m_em*eps

def heat_rate(eps_nuc,alpha,eps_th,t,t0):
  return eps_nuc*(0.5 - 1./np.pi * np.arctan((t-t0)/sigma0))**alpha * (eps_th/0.5)

def r_ph_calc(vel,t):
  return (vel*c)*t

def T_eff_calc(Lum,dOmega,r_ph):
  return (Lum/(dOmega*r_ph**2*sigma_SB))**(1./4.)

def planckian(nu,T_plk):
  return (2.*h*nu**3)/(c*c)/(np.exp((h*nu)/(kB*T_plk))-1.)

def m_filter(name,T_ray,rad_ray,dist,ff):
  for i in range(n_filters):
    if (name == filter_name[i]):
      break
    if (i == n_filters):
      sys.exit("Problem with the filter name!")
  fnu = calc_fnu(filter_center_nu[i],T_ray,rad_ray,dist,ff)
  return -2.5*np.log10(fnu)-48.6

def calc_fnu(nu,T,rad,dist,ff):
  tmp = 0.
  for k in range (n_k):
    tmp1 = ff[k]*planckian(nu,T[k]) # [erg/cm^2/s/Hz]
    tmp = tmp + tmp1 * (rad[k]/dist)**2               # flux [erg/s/cm^2]
  return tmp

def model_single_spherical(n_v,Omega,m_ej,v_exp,v_min,kappa,eps_th,vel,t_diff,t_fs,m_vel):  #,L_b,r_ph,T_eff):
  for i in range(n_v):
    vel[i]    = v_min + float(i)*(v_exp-v_min)/float(n_v)
    m_vel[i]  = mass_gt_v(vel[i],m_ej,v_exp,Omega)
    t_diff[i] = t_diff_v(kappa,vel[i],m_vel[i],Omega)
    t_fs[i]   = t_fs_v(kappa,vel[i],m_vel[i],Omega)
    eps_th[i] = therm_efficiency(t_diff[i],m_ej,Omega,v_exp)  #to be checked!
  return 0

def therm_efficiency_params(m,omega,v):
# define the interpolation functions
  m_iso = 4.*np.pi/omega * m
  x = [np.log10(1.e-3),np.log10(5e-3),np.log10(1e-2),np.log10(5e-2)]
  y = [0.1,0.2,0.3]
  a = [[2.01,0.81,0.56,0.27],[4.52,1.90,1.31,0.55],[8.16,3.20,2.19,0.95]]
  b = [[0.28,0.19,0.17,0.10],[0.62,0.28,0.21,0.13],[1.19,0.45,0.31,0.15]]
  d = [[1.12,0.86,0.74,0.60],[1.39,1.21,1.13,0.90],[1.52,1.39,1.32,1.13]]
  fa = interpolate.interp2d(x, y, a, kind='linear')
  fb = interpolate.interp2d(x, y, b, kind='linear')
  fd = interpolate.interp2d(x, y, d, kind='linear')
# assign the values of the mass and velocity
  xnew=np.log10(m_iso)   #mass     [Msun]
  ynew=v             #velocity [c]
# compute the parameters by linear interpolation in the table
  return [fa(xnew,ynew),fb(xnew,ynew),fd(xnew,ynew)]

def therm_efficiency(t,m,omega,v):
  coeff=np.zeros(3)
  coeff=therm_efficiency_params(m,omega,v)
# time "t" is in days
  time_days=t*sec2day
  tmp = 2.*coeff[1]*time_days**coeff[2]
  return 0.36*(np.exp(-coeff[0]*time_days) + np.log(1.+tmp)/tmp )
#  return 0.5

def calc_eps_nuc(ye,time,val_min):
  if (ye==1.):
    tmp=val_min * ft_eps_nuc(time)
  elif (ye==0.):
    tmp=val_min
  elif (ye==0.5):
    tmp=0.5*(val_min + val_min*ft_eps_nuc(time))
  else:
    print('Wrong Ye value!')
    sys.exit()
  return tmp 

def ft_eps_nuc(time):
  time_day=time*sec2day
  tmp = min(max(4*time_day-4.,-20),20)
  return 0.5 + 2.5/(1.+np.exp(tmp))

m_ej_dyn  = np.zeros(n_k)
v_min_dyn = np.zeros(n_k)
v_exp_dyn = np.zeros(n_k)
kappa_dyn = np.zeros(n_k)
ye_dyn    = np.zeros(n_k)

m_ej_w  = np.zeros(n_k)
v_min_w = np.zeros(n_k)
v_exp_w = np.zeros(n_k)
kappa_w = np.zeros(n_k)
ye_w = np.zeros(n_k)

m_ej_vis  = np.zeros(n_k)
v_min_vis = np.zeros(n_k)
v_exp_vis = np.zeros(n_k)
kappa_vis = np.zeros(n_k)
ye_vis = np.zeros(n_k)

eps_th_dyn = np.zeros((n_v,n_k))
eps_th_w = np.zeros((n_v,n_k))
eps_th_vis = np.zeros((n_v,n_k))
vel_dyn = np.zeros((n_v,n_k))
vel_w = np.zeros((n_v,n_k))
vel_vis = np.zeros((n_v,n_k))
t_diff_dyn = np.zeros((n_v,n_k))
t_diff_w = np.zeros((n_v,n_k))
t_diff_vis = np.zeros((n_v,n_k))
t_fs_dyn = np.zeros((n_v,n_k))
t_fs_w = np.zeros((n_v,n_k))
t_fs_vis = np.zeros((n_v,n_k))
m_vel_dyn = np.zeros((n_v,n_k))
m_vel_w = np.zeros((n_v,n_k))
m_vel_vis = np.zeros((n_v,n_k))

time_r_min = np.zeros(n_k)
time_r_max = np.zeros(n_k)

time=np.zeros(n_time)

L_bol_tot =np.zeros((n_time,n_k))
r_ph_tot  =np.zeros((n_time,n_k))
T_eff_tot =np.zeros((n_time,n_k))

L_obs =np.zeros((n_time))
Lum = np.zeros(n_time)
mag_tot=np.zeros((n_time,n_filters))

  
chi2_file='output/chi2.txt'
output_file='output/output.txt'

f = open(output_file,'w')
f.close()

f = open(chi2_file,'w')
f.close()

n_ang = len(ang_deg_sample)
flux_fact=np.zeros((n_ang,n_k))
for mm in range(n_ang):
  ang_deg=ang_deg_sample[mm]
# load the projection factor
  flux_file='../../mk_source/source/output_12_rays/flux_obs_ang_'+str(ang_deg)+'_nrays_'+str(n_k)+'.dat'
  tmp_col1,tmp_col2=np.loadtxt(flux_file,unpack=True,usecols=(0,4))
  for k in range(n_k):
    flux_fact[mm,k]=float(tmp_col2[k])

for ii in range(len(m_ej_dyn_tot_sample)):
  m_ej_dyn_tot = m_ej_dyn_tot_sample[ii]
  print('I am doing m_ej_dyn_tot',m_ej_dyn_tot)
  for jj in range(len(m_disk_sample)):
    m_disk=m_disk_sample[jj]
    print('I am doing m_disk',m_disk)
    for kk in range(len(xi_w_sample)):
      xi_w = xi_w_sample[kk]
      print('I am doing xi_w',xi_w)
      for ll in range(len(xi_vis_sample)):
        xi_vis = xi_vis_sample[ll]
        print('I am doing xi_vis',xi_vis)

        # dynamic ejecta
        istat = mass_ang_dist(n_k,m_ej_dyn_tot,'sin2',m_ej_dyn)
        istat = vel_dist(n_k,'lin','unif_min',1.e-6,1.e-6,v_min_dyn)
        istat = vel_dist(n_k,'lin',v_dyn_dist,v_dyn_max,v_dyn_max,v_exp_dyn)
        istat = ye_dist(n_k,ye_dist_dyn,0.5,ye_dyn)
        istat = kappa_dist(n_k,ye_dyn,kappa_dyn_min,kappa_dyn_max,kappa_dyn)
        
        #for k in range(n_k):
        #  print(k,m_ej_dyn[k],v_exp_dyn[k],ye_dyn[k],kappa_dyn[k])
        #print()
          
        # wind ejecta
        m_ej_w_tot = m_disk*xi_w  #solar masses
        istat = mass_ang_dist(n_k,m_ej_w_tot,'step60',m_ej_w)
        istat = vel_dist(n_k,'lin','unif_min',1.e-6,1.e-6,v_min_w)
        istat = vel_dist(n_k,'lin','unif_max',0.06,v_wind_max,v_exp_w)
        istat = ye_dist(n_k,ye_dist_w,0.5,ye_w)
        istat = kappa_dist(n_k,ye_w,kappa_min_w,kappa_max_w,kappa_w)
        
        #for k in range(n_k):
        #  print(k,m_ej_w[k],v_exp_w[k],ye_w[k],kappa_w[k])
        #print()
        
        # viscous ejecta
        m_ej_vis_tot = m_disk*xi_vis  #solar masses
        istat = mass_ang_dist(n_k,m_ej_vis_tot,'sin2',m_ej_vis)
        istat = vel_dist(n_k,'lin','unif_min',1.e-6,1.e-6,v_min_vis)
        istat = vel_dist(n_k,'lin','unif_max',0.05,v_vis_max,v_exp_vis)
        istat = ye_dist(n_k,'unif',0.5,ye_vis)
        kappa_max_vis = kappa_min_vis   # here we assume uniform, so we take min=max
        istat = kappa_dist(n_k,ye_vis,kappa_min_vis,kappa_max_vis,kappa_vis)
          
        #for k in range(n_k):
        #  print(k,m_ej_vis[k],v_exp_vis[k],ye_vis[k],kappa_vis[k])
          
        for k in range(n_k):
          
          istat = model_single_spherical(n_v,solid_angle[k],m_ej_dyn[k],v_exp_dyn[k],v_min_dyn[k],kappa_dyn[k],eps_th_dyn[:,k],vel_dyn[:,k],t_diff_dyn[:,k],t_fs_dyn[:,k],m_vel_dyn[:,k])
          
          istat = model_single_spherical(n_v,solid_angle[k],m_ej_w[k],v_exp_w[k],v_min_w[k],kappa_w[k],eps_th_w[:,k],vel_w[:,k],t_diff_w[:,k],t_fs_w[:,k],m_vel_w[:,k])
          
          istat = model_single_spherical(n_v,solid_angle[k],m_ej_vis[k],v_exp_vis[k],v_min_vis[k],kappa_vis[k],eps_th_vis[:,k],vel_vis[:,k],t_diff_vis[:,k],t_fs_vis[:,k],m_vel_vis[:,k])
          
          time_r_min[k]=1.1*max(t_fs_w[n_v-1,k],max(t_fs_dyn[n_v-1,k],t_fs_vis[n_v-1,k]))
          time_r_max[k]=0.9*min(t_diff_w[0,k],min(t_diff_dyn[0,k],t_diff_vis[0,k]))
          
        # global time extremes
        time_min=max(max(time_r_min),3600.)
        time_max=min(min(time_r_max),1728000.)
          
        for i in range(n_time):
          time[i] = time_min + (time_max-time_min)*float(i)/float(n_time)
          
        for k in range(n_k):
          
          fun_v_diff_dyn_t = interpolate.interp1d(t_diff_dyn[::-1,k],vel_dyn[::-1,k])
          fun_v_fs_dyn_t   = interpolate.interp1d(t_fs_dyn[::-1,k],vel_dyn[::-1,k])
          fun_mv_dyn_v     = interpolate.interp1d(vel_dyn[:,k],m_vel_dyn[:,k])
          fun_v_diff_w_t = interpolate.interp1d(t_diff_w[::-1,k],vel_w[::-1,k])
          fun_v_fs_w_t   = interpolate.interp1d(t_fs_w[::-1,k],vel_w[::-1,k])
          fun_mv_w_v     = interpolate.interp1d(vel_w[:,k],m_vel_w[:,k])
          fun_v_diff_vis_t = interpolate.interp1d(t_diff_vis[::-1,k],vel_vis[::-1,k])
          fun_v_fs_vis_t   = interpolate.interp1d(t_fs_vis[::-1,k],vel_vis[::-1,k])
          fun_mv_vis_v     = interpolate.interp1d(vel_vis[:,k],m_vel_vis[:,k])
          
          fun_eps_th_dyn_t = interpolate.interp1d(t_diff_dyn[::-1,k],eps_th_dyn[::-1,k])
          fun_eps_th_w_t   = interpolate.interp1d(t_diff_w[::-1,k],eps_th_w[::-1,k])
          fun_eps_th_vis_t = interpolate.interp1d(t_diff_vis[::-1,k],eps_th_vis[::-1,k])
          
          for i in range(n_time):
          
            v_diff_dyn = fun_v_diff_dyn_t(time[i])
            v_diff_w   = fun_v_diff_w_t(time[i])
            v_diff_vis = fun_v_diff_vis_t(time[i])
          
            v_fs_dyn = fun_v_fs_dyn_t(time[i])
            v_fs_w   = fun_v_fs_w_t(time[i])
            v_fs_vis = fun_v_fs_vis_t(time[i])
          
            m_v_diff_dyn = fun_mv_dyn_v(v_diff_dyn)
            m_v_diff_w   = fun_mv_w_v(v_diff_w)
            m_v_diff_vis = fun_mv_vis_v(v_diff_vis)
          
            m_v_fs_dyn = fun_mv_dyn_v(v_fs_dyn)
            m_v_fs_w   = fun_mv_w_v(v_fs_w)
            m_v_fs_vis = fun_mv_vis_v(v_fs_vis)
          
            m_emis_dyn = m_v_diff_dyn - m_v_fs_dyn
            m_emis_w   = m_v_diff_w - m_v_fs_w
            m_emis_vis = m_v_diff_vis - m_v_fs_vis
          
            e_th_dyn = fun_eps_th_dyn_t(time[i])
            e_th_w   = fun_eps_th_w_t(time[i])
            e_th_vis = fun_eps_th_vis_t(time[i])
          
            eps_nuc_dyn = calc_eps_nuc(ye_dyn[k],time[i],eps0)
            eps_nuc_w   = calc_eps_nuc(ye_w[k],time[i],eps0)
            eps_nuc_vis = calc_eps_nuc(ye_vis[k],time[i],eps0)
          
            L_bol_dyn  = bol_lum(alpha,eps_nuc_dyn,e_th_dyn,time[i],t0eps,m_emis_dyn)
            L_bol_w    = bol_lum(alpha,eps_nuc_w,e_th_w,time[i],t0eps,m_emis_w)
            L_bol_vis  = bol_lum(alpha,eps_nuc_vis,e_th_vis,time[i],t0eps,m_emis_vis)
          
            r_ph_dyn   = r_ph_calc(v_fs_dyn,time[i])
            r_ph_w     = r_ph_calc(v_fs_w,time[i])
            r_ph_vis   = r_ph_calc(v_fs_vis,time[i])
          
            r_ph_tot[i,k]  = max(r_ph_dyn,r_ph_w,r_ph_vis)
            L_bol_tot[i,k] = L_bol_dyn + L_bol_w + L_bol_vis
            T_eff_tot[i,k] = T_eff_calc(L_bol_tot[i,k],solid_angle[k],r_ph_tot[i,k])

        for mm in range(n_ang):

          ang_deg=ang_deg_sample[mm]

          # string with all the parameters
          parameter_string=[m_ej_dyn_tot,v_dyn_max,ye_dist_dyn,kappa_dyn_min,kappa_dyn_max,m_disk,xi_w,v_wind_max,ye_dist_w,kappa_min_w,kappa_max_w,xi_vis,v_vis_max,kappa_min_vis,eps0,s_dist_pc,ang_deg]
          
          for i in range(n_time):
            for k in range(n_k):
              L_obs[i] = L_obs[i] + 4.*np.pi/solid_angle[k] * flux_fact[mm,k]/np.pi * L_bol_tot[i,k]
          
          for i in range(n_time):
            for j in range(n_filters):
              mag_tot[i,j] = m_filter(filter_name[j],T_eff_tot[i,:],r_ph_tot[i,:],float(s_dist_cm),flux_fact[mm,:])

#####
# check on the filter peaks
#####
          min_mag = np.zeros(n_filters)
          count = 0
          for j in range(n_filters):
            ntmp = np.argmin(mag_tot[:,j])
            min_mag[j] = mag_tot[ntmp,j]
            if (filter2use[j] > 0.1):
              if (min_mag[j] > mpeak_lower_limit[j] and min_mag[j] < mpeak_upper_limit[j]):
                count = count + 1

          if (count < num_filter2use ):
            continue

########
#  chi2
########

          jump = 0
          chi2 = np.zeros(n_filters)
          for kk in range(n_filters):

            if (filter2chi[kk] < 0.1):
              continue
  
            counter = 0
            y = np.zeros(n_time)
            for l in range(n_time):
              y[l] = float(mag_tot[l,kk])   #apparent magnitudes
            fun = interpolate.interp1d(time,y,kind='linear')
            nb = len(t_band[:][kk])
            for jj in range(nb):
              if (t_band[jj][kk] > time[0]  and t_band[jj][kk] < time[n_time-1]):
                counter = counter + 1
                mag_interp = fun(t_band[jj][kk])
                chi2[kk] = chi2[kk] + (mag_interp - m_band[jj][kk])**2
            if (counter > 0):
              chi2[kk] = chi2[kk]/counter
            else:
              chi2[kk] = -1.
              jump = 1
              break

            chi2_tot = chi2_tot + chi2[kk]

          if (jump == 1):
            continue


          f = open(output_file,'a')
          
          f.write('%s\n' %parameter_string)
          for i in range(n_time):
            lum_phys = sum(L_bol_tot[i,:])
            f.write('%16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e \n' %(time[i],L_obs[i],lum_phys,mag_tot[i,0],mag_tot[i,1],mag_tot[i,2],mag_tot[i,3],mag_tot[i,4],mag_tot[i,5],mag_tot[i,6],mag_tot[i,7],mag_tot[i,8],mag_tot[i,9],mag_tot[i,10],mag_tot[i,11],mag_tot[i,12]))
        
          f.write("\n") 
          f.close
          
          g = open(chi2_file,'a')
          g.write('%s %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e %16e \n' %(parameter_string,chi2[0:12],chi2_tot))
  
          g.write("\n") 
          g.close


