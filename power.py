#new model...
import numpy as np
import pickle
mname = 'ASW000102p/WM4H5RZXQZ_hires'
#mname = 'ASW0000h2m/IHRULOMX6D'
fil = open(mname+'.pkl')
chutney = pickle.load(fil)
ensem = chutney['grids']
pixrad = chutney['pixrad']
N = 2*pixrad+1
R = chutney['maprad']

def power(X,Y,a,b,h,t):    
#    x,y = X*np.cos(phi) - Y*np.sin(phi),X*np.sin(phi) + Y*np.cos(phi)
#    R = np.sqrt(q*q*x*x + y*y + 1e-12)
    R = (a*X*X + b*Y*Y + 2*h*X*Y + 1e-12)
    t = 1
    kappa = R**-t
    return kappa
    
def profile(params):
    phi,q,b,t = params[0],params[1],params[2],params[3]
    x = np.linspace(-R,R,N)
    X,Y = np.meshgrid(x,x)
    P = power(X,Y,phi,q,b,t)
    return P

