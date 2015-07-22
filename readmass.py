"""readmass""" #currently set up to work with power.py but can be changed

#import libraries
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize as opt
#import pickle

"""
#open data file
mname = 'ASW0007k4r/012771'
#mname = 'ASW0000h2m/007022'
#mname = 'ASW0000h2m/IHRULOMX6D'
#mname = 'ASW000102p/WM4H5RZXQZ_hires' #hi res makes a difference
#mname = 'gribbles'
fil = open(mname+'.pkl')
ensem = pickle.load(fil)
"""

from ellip import ensem, N, R, profile, maximgpos #, grid      #ellip is set up for lenses 102p and h2m/IHRU... (hires)

#from power import ensem, N, R, profile                           #an experiment with a new model...

#defining the data as describing points on a 2D surface
#N = ensem[0].shape[0]                      #N and R now defined within ellip
#R = (N-1)/2
x = np.linspace(-R,R,N)
X,Y = np.meshgrid(x,x)


fudge = 0.4312                              #red shift fudge factor


for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))   #reshape 2D array as 1D array 
    ensem1d = fudge*ensem1d                 #multiply by the fudge factor for red shift
    if m==0:    
        sum = ensem1d
    else:    
        sum = sum+ensem1d                   #sum elements in the ensemble
mean = sum/len(ensem)                       #calculate mean


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


"""Alternative parameterised form not currently used here"""
def Xprof(params):
    a,n,b,c,h=params[0],params[1],params[2],params[3],params[4]
    #return a*(1+X*X+Y*Y)**-n               #define test parametrized functional form k(param)
    return a*(h+c*X*X+Y*Y+b*X*Y)**-n

#from ellip import profile, grid             #Import isothermal ellipsoid functional form
#profile=Xprof
#mask = (1-np.sign(X*X+Y*Y-maximgpos*maximgpos))/2
#mask = np.reshape(mask,(N**2))

def residuals(params):
    f = profile(params)                         #f = k(param)
    f = np.reshape(f,(N**2))                    #reshape f into 1D array
    f -= mean                                   #changed 'change' (which was an experiment) to 'mean'#f = f-change
    for m in range(1,6):
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]  #removing projections along principle axes
    return f #mask*f


#ini = [1.8,0.8,-0.1,0.3]                             #initial values for parameters
ini = [1.8,0.8,-0.1]
lsq = opt.leastsq(residuals,ini)[0]             #perform least squares fit on f


"""Print out parameters"""
param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
#param4 = lsq[3]
#maybe user can input required no of sf?
#print 'Param1 (Einstein radius) = {0:.3e}, Param2 (Ellipticity) = {1:.3e}, Param3 (Position angle of ellipticity) = {2:.3e}'.format(param1,param2,param3) #prints values of optimised parameters
print 'phi = {0:.3e}, q = {1:.3e}, b = {2:.3e}'.format(param1,param2,param3)#, t = {3:.3e}'.format(param1,param2,param3,param4)

#grav = grid(lsq)                                 #For plotting grav ptl
#lev = np.linspace(np.amin(grav),np.amax(grav),21)
#pl.contour(X,Y,grav, levels=lev)                 #plot of gravitational potential


"""Plot parameterised model"""
F = profile(lsq)                                  #F = k(param) with optimised parameters
##F1d = np.reshape(F,N**2)
##pl.plot(F1d)
lev = np.linspace(0,5,11)
pl.contour(X,Y,F, levels=lev)             #plot graph of parameterised model
"""plot colour-filled contours"""
#bar = pl.contourf(X,Y,F,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)
#pl.title('Param')
#pl.show()


"""Plot mean"""
meanplot = np.reshape(mean,(N,N))
#pl.contour(X,Y,meanplot, levels=[0,1,2,3,4])     #plot graph of mean
"""plot colour-filled contours"""
lev = np.linspace(np.amin(meanplot),np.amax(meanplot),10)
bar = pl.contourf(X,Y,meanplot,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)
pl.title('Param and Mean')
pl.show()


"""Plot change"""
F1d = np.reshape(F,N**2)
change = mean
for m in range(1,6):
    change += np.inner(F1d,vecs[:,-m])*vecs[:,-m] #adding projections (of the parameterised form along the eigenvectors) to the mean
##pl.plot(change)
H = np.reshape(change,(N,N))
lev = np.linspace(0,5,11)
#pl.contour(X,Y,H, levels=lev)            #plot graph of 'change' on same graph - these are the points on the MoI ellipse that are closest to the parameterised form
"""plot colour-filled contours"""
#bar = pl.contourf(X,Y,H,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)
#pl.title('Param and Change')
#pl.show()


"""Plot difference between parameterised model and 'change'"""
K = residuals(lsq)
K = np.reshape(K,(N,N))
#pl.contour(X,Y,K, levels=[0,1,2,3,4])            #plot graph of change - k(param)
"""plot colour-filled contours"""
lmax = np.amax(abs(K))
lev = np.linspace(-lmax,lmax,50)
bar = pl.contourf(X,Y,K,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)
pl.axes().set_aspect('equal')
pl.title('Param - Change')
pl.show()        



