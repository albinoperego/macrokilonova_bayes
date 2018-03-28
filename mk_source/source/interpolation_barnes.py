import scipy.interpolate as itp

# INTERPOLATION of A,B,D PARAMETERS (from Barnes et al. 2016) for THERMALISATION EFFICIENCY 
# Now it uses a single variable [m_ej/(v_ej^2)] instead of using m_ej and v_ej separately.
# Input: (M_ej/Msun) / (v_ej/c)**2
# Note: if input is outside the interpolation range, the output is 0. Then the thermalisation efficiency will be fixed at 0.5

x_barnes = [0.011,0.025,0.0556,0.1,0.111,0.125,0.25,0.5,0.5556,1.,1.25,5.]
a_barnes = [8.16,4.52,3.20,2.01,2.19,1.90,1.31,0.81,0.95,0.56,0.55,0.27]
b_barnes = [1.19,0.62,0.45,0.28,0.31,0.28,0.21,0.19,0.15,0.17,0.13,0.10]
d_barnes = [1.52,1.39,1.39,1.12,1.32,1.21,1.13,0.86,1.13,0.74,0.90,0.60]
interp_a=itp.interp1d(x_barnes,a_barnes,fill_value='extrapolate')
interp_b=itp.interp1d(x_barnes,b_barnes,fill_value='extrapolate')
interp_d=itp.interp1d(x_barnes,d_barnes,fill_value='extrapolate')
def int_a(x):
	aa = interp_a(x)
	if(aa>0):
		return(aa)
	else:
		return(0)	
def int_b(x):
	bb = interp_b(x)
	if(bb>0):
		return(bb)
	else:
		return(0)
def int_d(x):
	dd = interp_d(x)
	if(dd>0):
		return(dd)
	else:
		return(0)
