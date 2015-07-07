import numpy as np
#import matplotlib.pyplot as pl
import pickle
mname = 'ASW000102p/WM4H5RZXQZ_hires' #hi res makes a difference
fil = open(mname+'.pkl')
chutney = pickle.load(fil)
ensem = chutney['grids']
pixrad = chutney['pixrad']
N = 2*pixrad+1
R = chutney['maprad']

# See eqns (33-35) from Keeton astro-ph/0102341
def poten_SIE(x,y,reinst,ell,ell_pa): #ell_pa?
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
#s=0 but why? s is the scale radius but what does that mean?
#"in limit of singular (s=0) and spherical (q=1) model, b = Einst radius"
#but we are assuming b=reinst and s=0 and not q=1
    #all sims have density go to infinity at origin so s=0 (makes eqns simpler) and this is a reasonable approximation to make for real galaxies (often black hole at centre)
    
# kappa(zl=0.5,zs=1) = 0.4312*kappa_inf(zl=0.5)



def grid(params):
    reinst,ell,ell_pa = params[0],params[1],params[2]    
    x = np.linspace(-R,R,N)
    X,Y = np.meshgrid(x,x)
    F = poten_SIE(X,Y,reinst,ell,ell_pa)
    return F

def profile(params):
    reinst,ell,ell_pa = params[0],params[1],params[2]
    S = 3
    r = R + 0.5*(1-1./S)*R/pixrad
    x = np.linspace(-r,r,N*S)
    X,Y = np.meshgrid(x,x)
    F = poten_SIE(X,Y,reinst,ell,ell_pa)
    M = 0*F
    M[1:-1,1:-1] = F[2:,1:-1] + F[:-2,1:-1] + F[1:-1,2:] + F[1:-1,:-2] \
                 - 4*F[1:-1,1:-1]
    K = np.ndarray(shape=(N,N))
    for i in range(N):
        for j in range(N):
            K[i,j] = np.sum(M[i*S:(i+1)*S,j*S:(j+1)*S])
    return K
    
#params = [10,0.5,45]
#K = profile(params)

#R = (N-1)/2
#x = np.linspace(-R,R,N)
#X,Y = np.meshgrid(x,x)

# M = np.log(abs(M)+1e-12)
#lev = np.linspace(np.amin(K),np.amax(K),21)
#pl.contour(X,Y,K, levels=[0,1,2,3,4,5])
#pl.axes().set_aspect('equal')
#pl.show()        

