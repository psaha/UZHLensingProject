"""ellip"""

"""This program imports the lens data and calculates the parameterised form of
the mass distribution for an isothermal ellipsoid from the functional form of 
the gravitational potential. It also defines N, R and maximgpos such that 
parameters are given in the correct units. For use with mass_ellip or readmass"""


import numpy as np
import pickle
mname = 'ASW000102p/WM4H5RZXQZ_hires'
#mname = 'ASW0000h2m/IHRULOMX6D'
#mname = 'WM4H5RZXQZ'
mname = 'ASW0002b6m/003250'
fil = open(mname+'.pkl')
chutney = pickle.load(fil)
ensem = chutney['grids']                                    #ensem = the ensemble of 200 free-form mass distributions for the lens
pixrad = chutney['pixrad']                                  #pixrad = radius in number of pixels
N = 2*pixrad+1
R = chutney['maprad']                                       #maprad = radius from central point to central point of outer tile (in arcseconds)
maximgpos = chutney['maximgpos']


# See eqns (33-35) from Keeton astro-ph/0102341
def poten_SIE(x,y,reinst,ell,ell_pa):                       #parameterised function for isothermal ellipsoid
    pa = (ell_pa+90)*np.pi/180
    q = np.sqrt((1-ell)/(1+ell))
    reinst *= np.sqrt((1+q*q)/(2*q*q))
    cs,sn = np.cos(pa),np.sin(pa)
    x,y = cs*x + sn*y, -sn*x + cs*y
    A = reinst*q/np.sqrt(1-q*q)
    B = np.sqrt((1-q*q)/(q*q*x*x + y*y + 1e-12))
    phix = A*np.arctan(B*x)
    phiy = A*np.arctanh(B*y)
    return x*phix + y*phiy


# kappa(zl=0.5,zs=1) = 0.4312*kappa_inf(zl=0.5)             #red shift fudge factor


def grid(params):
    reinst,ell,ell_pa = params[0],params[1],params[2]       #defining grav ptl as describing points on a 2D surface (not actually used unless you want graph of grav ptl)
    x = np.linspace(-R,R,N)
    X,Y = np.meshgrid(x,x)
    F = poten_SIE(X,Y,reinst,ell,ell_pa)
    return F


#calculate parameterised functional form of mass distribution
def profile(params):
    reinst,ell,ell_pa = params[0],params[1],params[2]       #parameters: Einstein radius, Ellipticity and Position angle of ellipticity
    S = 3
    r = R + 0.5*(1-1./S)*R/pixrad                           #oversampling, taking care about size of grid
    x = np.linspace(-r,r,N*S)
    X,Y = np.meshgrid(x,x)
    F = poten_SIE(X,Y,reinst,ell,ell_pa)
    M = 0*F
    M[1:-1,1:-1] = F[2:,1:-1] + F[:-2,1:-1] + F[1:-1,2:] + F[1:-1,:-2] \
                 - 4*F[1:-1,1:-1]                           #taking the second difference of the gravitational potential
    M = M*(pixrad/R)**2 / 2                                     #dividing by "delta(x)^2" to obtain grad^2 of potential i.e. mass distrib.
    K = np.ndarray(shape=(N,N))                             #undersampling again to go back to original grid size
    for i in range(N):
        for j in range(N):
            K[i,j] = np.sum(M[i*S:(i+1)*S,j*S:(j+1)*S])
    return K
    
