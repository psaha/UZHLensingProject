import numpy as np
import matplotlib.pyplot as pl

# See eqns (33-35) from Keeton astro-ph/0102341
def poten_SIE(x,y,reinst,ell,ell_pa):
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

# kappa(zl=0.5,zs=1) = 0.4312*kappa_inf(zl=0.5)

N = 25
S = 5
R = (N-1)/2
x = np.linspace(-R,R,N*S)
X,Y = np.meshgrid(x,x)

F = poten_SIE(X,Y,10,0.5,45)

lev = np.linspace(np.amin(F),np.amax(F),21)
pl.contour(X,Y,F, levels=lev)
pl.axes().set_aspect('equal')
pl.show()        

