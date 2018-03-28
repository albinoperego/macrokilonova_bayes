import numpy as np

# MAGNITUDE CORRECTION FOR DEREDDENING
# Input wavelength in nanometers!!!

def dered_CCM(wave, EBV=0.105, R_V=3.1):
    '''
    - Input: 
		wave 1D array (Nanometers)
		EBV: E(B-V) (default 0.105)
		R_V: Reddening coefficient to use (default 3.1)    
    - Output: 
		A_lambda correction according to the CCM89 Law.
    '''
    x = 1000./ wave  #Convert to inverse microns
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ## Infrared ##
    mask = (x > 0.3) & (x < 1.1)
    if np.any(mask):
        a[mask] =  0.574 * x[mask]**(1.61)
        b[mask] = -0.527 * x[mask]**(1.61)

    ## Optical/NIR ##
    mask = (x >= 1.1) & (x < 3.3)
    if np.any(mask):
        xxx = x[mask] - 1.82
        # c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085, #Original
        #        0.01979, -0.77530,  0.32999 ]               #coefficients
        # c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434, #from CCM89
        #       -0.62251,  5.30260, -2.09002 ]
        c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,     #New coefficients
              -1.718,   -0.827,    1.647, -0.505 ]         #from O'Donnell
        c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,     #(1994)
               11.102,    5.491,  -10.805,  3.347 ]
        a[mask] = np.poly1d(c1[::-1])(xxx)
        b[mask] = np.poly1d(c2[::-1])(xxx)

    ## Mid-UV ##
    mask = (x >= 3.3) & (x < 8.0)
    if np.any(mask):
        F_a = np.zeros_like(x[mask])
        F_b = np.zeros_like(x[mask])
        mask1 = x[mask] > 5.9
        if np.any(mask1):
            xxx = x[mask][mask1] - 5.9
            F_a[mask1] = -0.04473 * xxx**2 - 0.009779 * xxx**3
        a[mask] = 1.752 - 0.316*x[mask] - (0.104 / ( (x[mask]-4.67)**2 + 0.341 )) + F_a
        b[mask] = -3.090 + 1.825*x[mask] + (1.206 / ( (x[mask]-4.62)**2 + 0.263 )) + F_b

    ## Far-UV ##
    mask = (x >= 8.0) & (x < 11.0)
    if np.any(mask):
        xxx = x[mask] - 8.0
        c1 = [ -1.073, -0.628,  0.137, -0.070 ]
        c2 = [ 13.670,  4.257, -0.420,  0.374 ]
        a[mask] = np.poly1d(c1[::-1])(xxx)
        b[mask] = np.poly1d(c2[::-1])(xxx)

    #Now compute extinction correction
    A_V = R_V * EBV
    A_lambda = A_V * (a + b/R_V)
    return A_lambda
