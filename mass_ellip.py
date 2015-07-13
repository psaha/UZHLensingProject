"""mass_ellip"""

"""This program is a version of readmass specifically for modelling
mass distributions of lenses using the isothermal ellipsoid parameterised 
functional form found in ellip.py
The lenses themselves are also imported here from ellip.py"""

#****************
"""Set up data"""
#****************

#import libraries
import numpy as np                          #for calculations 
import matplotlib.pyplot as pl              #for plotting graphs
import scipy.optimize as opt                #for the least squares fit
from metrop import samp


#import lens data (ensemble of 200 free-form mass distributions), N, R, parameterised functional form, maximgpos
from ellip import ensem, N, R, profile, maximgpos, mname

print 'Lens: '
print mname


#defining the data as describing N*N points on a 2D surface in the range -R to R
x = np.linspace(-R,R,N)
X,Y = np.meshgrid(x,x)


#define red shift fudge factor
fudge = 0.4312

#*****************************
"""Calculate inertia tensor"""
#*****************************

#calculate mean of ensemble
for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))   #reshape 2D array as 1D array 
    ensem1d = fudge*ensem1d                 #multiply by the fudge factor for red shift
    if m==0:    
        sum = ensem1d
    else:    
        sum = sum+ensem1d                   #sum elements in the ensemble
mean = sum/len(ensem)                       #calculate mean


#calculate moment of inertia tensor (a.k.a covariance matrix) and its eigenvectors and eigenvalues
for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))   #reshape 2D array as 1D array
    ensem1d = fudge * ensem1d               #multiply by the fudge factor for red shift
    diff = ensem1d - mean                   #delta(k) = datum k - mean <k>
    out = np.outer(diff,diff)               #outer products of pairs of delta(k)
    if m==0:                                #create MoI tensor (outer products of pairs of values)
        outsum = out
    else:
        outsum = outsum + out
outsum = outsum/len(ensem)                  #scale MoI by tensor size of ensemble (appears to make no difference to output)
vals,vecs = np.linalg.eigh(outsum)          #find eigenvecs/vals of MoI tensor


#*****************************************
"""Define parameterised functional form"""
#*****************************************

#define masking function for which selects only the data points that are inside the image of the lens
mask = (1-np.sign(X*X+Y*Y-maximgpos*maximgpos))/2
mask = np.reshape(mask,(N**2))


#define residuals function: parameterised form - mean - projections along principle axes multiplied by the masking function
def residuals(params):
    f = profile(params)                         #f = k(param)
    f = np.reshape(f,(N**2))                    #reshape f into 1D array
    f -= mean                                   #take away the mean
    for m in range(1,6):
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]  #removing projections along principle axes
    return mask*f                               

#**************************
"""Linear regression fit"""
#**************************

#set initial parameter values and perform linear regression fit of f to the ensemble
ini = [1,0.1,0.1]                             #initial values for parameters
lsq = opt.leastsq(residuals,ini)[0]             #perform parameter optimisation on residuals


#Print out parameters
param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
#print 'Einstein radius = {0:.3}, Ellipticity = {1:.3}, Position angle of ellipticity = {2:.3}'.format(param1,param2,param3) #prints values of optimised parameters
print 'Least squares parameters'
print '%.2f' %param1, '\t', '%.2f' %param2, '\t', '%.2f' %param3

#*************
"""MCMC fit"""
#*************

def lnprob(params):
    r = residuals(params)
    return np.sum(-r*r)

lo = [1,0,-90]
hi = [5,0.5,90]
S = 100
lnp,pars = samp(lnprob,lo,hi,S)
print lsq
lsq = np.mean(pars,axis=0)
sig = np.var(pars,axis=0)**.5
print 'mcmc estimate and errors:'
print lsq
print sig

#Print out parameters
param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
#print 'Einstein radius = {0:.3}, Ellipticity = {1:.3}, Position angle of ellipticity = {2:.3}'.format(param1,param2,param3) #prints values of optimised parameters
print 'MCMC parameters'
print '%.2f' %param1, '\t', '%.2f' %param2, '\t', '%.2f' %param3


#*********************************
"""Output graphs"""
#*********************************


"""Plot parameterised model"""
F = profile(lsq)                                  #F = k(param) with optimised parameters
lev = np.linspace(0,5,11)
pl.contour(X,Y,F, levels=lev)             #contour plot of parameterised model
"""plot colour-filled contours"""
#bar = pl.contourf(X,Y,F,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)
#pl.title('Param')
#pl.show()


"""Plot mean"""
#meanplot = np.reshape(mean,(N,N))
#pl.contour(X,Y,meanplot, levels=[0,1,2,3,4])     #plot graph of mean
"""plot colour-filled contours"""
#lev = np.linspace(np.amin(meanplot),np.amax(meanplot),10)
#bar = pl.contourf(X,Y,meanplot,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)
#pl.title('Param and Mean')
#pl.show()


"""Plot change: points on the moment of inertia ellipse that are closest to the parameterised functional form"""
F1d = np.reshape(F,N**2)
change = mean
for m in range(1,6):
    change += np.inner(F1d,vecs[:,-m])*vecs[:,-m] #adding projections (of the parameterised form along the eigenvectors) to the mean
H = np.reshape(change,(N,N))
lev = np.linspace(0,5,11)
#pl.contour(X,Y,H, levels=lev)                   #plot graph of 'change' on same graph - these are the points on the MoI ellipse that are closest to the parameterised form
"""plot colour-filled contours"""
bar = pl.contourf(X,Y,H,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)
pl.title('Parameterised form and "change"')
pl.show()


"""Plot residuals function with optimised parameters, i.e. the difference between parameterised model and 'change'"""
K = residuals(lsq)
K = np.reshape(K,(N,N))
#pl.contour(X,Y,K, levels=[0,1,2,3,4])            #contour plot of residuals
"""plot colour-filled contours"""
lmax = np.amax(abs(K))
lev = np.linspace(-lmax,lmax,50)
bar = pl.contourf(X,Y,K,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)
pl.axes().set_aspect('equal')
pl.title('Residuals: param - change')
pl.show()        


