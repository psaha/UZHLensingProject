#new model... seems to plot lots of contours in the centre and then a single contour further out - generally a bad fit
import numpy as np
import pickle

"""
mname = 'ASW000102p/WM4H5RZXQZ_hires'
#mname = 'ASW0000h2m/IHRULOMX6D'
fil = open(mname+'.pkl')
chutney = pickle.load(fil)
ensem = chutney['grids']
pixrad = chutney['pixrad']
N = 2*pixrad+1
R = chutney['maprad']
"""

#open data file
mname = 'ASW0007k4r/012771'
#mname = 'ASW0000h2m/007022'
#mname = 'ASW0000h2m/IHRULOMX6D'
#mname = 'ASW000102p/WM4H5RZXQZ_hires' #hi res makes a difference
#mname = 'gribbles'
fil = open(mname+'.pkl')
ensem = pickle.load(fil)
N = ensem[0].shape[0]                      #N and R now defined within ellip
R = (N-1)/2




def power(X,Y,phi,q,b,t):    
    q = 0.8 #value that makes ellipse closest to mean contours (but not change)
    x,y = X*np.cos(phi) - Y*np.sin(phi),X*np.sin(phi) + Y*np.cos(phi)
    R = np.sqrt(q*q*x*x + y*y + 1e-12)
#    R = (a*X*X + b*Y*Y + 2*h*X*Y + 1e-12)
    kappa = (1-t/2)*(abs(b)/R)**t
    return kappa
    
def profile(params):
    phi,q,b,t = params[0],params[1],params[2],params[3]
    x = np.linspace(-R,R,N)
    X,Y = np.meshgrid(x,x)
    P = power(X,Y,phi,q,b,t)
    return P

