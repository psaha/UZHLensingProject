#testfile: parameterisation with good fit but runtime warning

#import libraries
import pickle
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize as opt

#open data file
mname = 'ASW0007k4r/012771'
mname = 'ASW0000h2m/007022'
mname = 'ASW0000h2m/IHRULOMX6D'
mname = 'ASW000102p/WM4H5RZXQZ'
#mname = 'gribbles'
fil = open(mname+'.pkl')
ensem = pickle.load(fil)

#defining the data as describing points on a 2D surface
N = ensem[0].shape[0]
R = (N-1)/2
x = np.linspace(-R,R,N)
X,Y = np.meshgrid(x,x)


for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))   #reshape 2D array as 1D array 
    if m==0:    
        sum = ensem1d
    else:    
        sum = sum+ensem1d                   #sum elements in the ensemble
mean = sum/len(ensem)                       #calculate mean


for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))   #reshape 2D array as 1D array
    diff = ensem1d - mean                   #delta(k) = datum k - mean <k>
    out = np.outer(diff,diff)               #outer products of pairs of delta(k)
    if m==0:                                #create MoI tensor
        outsum = out
    else:
        outsum = outsum + out
outsum = outsum/len(ensem)                  #scale MoI by tensor size of ensemble
vals,vecs = np.linalg.eigh(outsum)          #find eigenvecs/vals of MoI tensor


def profile(params):
    a,b,n,h,c=params[0],params[1],params[2],params[3],params[4]
    #return a*(h+X*X+Y*Y+b*X*Y)**-n                #define test parameterised functional form k(param)
    return a*(h+c*X*X+Y*Y+b*X*Y)**-n
    #return a*(1+X*X+Y*Y)**-n

def residuals(params):
    f = profile(params)                     #f = k(param)
    f = np.reshape(f,(N**2))                #reshape f into 1D array
    f -= mean                               #f = f-mean (used to be f-change)
    for m in range(1,6):
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]         #removing projections along principle axes
    return f


ini = [1,-1.3,3.7,1.6,1]                               #initial values for parameters
lsq = opt.leastsq(residuals,ini)[0]                    #perform least squares fit on f
#print(lsq)

param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
param4 = lsq[3]
param5 = lsq[4]

print 'param1 = {0:.3e}, param2 = {1:.3e}, param3 = {2:.3e}, param4 = {3:.3e}, param5 = {4:.3e}'.format(param1,param2,param3,param4,param5) #prints values of optimised parameters



F = profile(lsq)                                    #F = k(param) with optimised parameters
pl.contour(X,Y,F, levels=[0,1,2,3,4])               #plot graph of parametrized model
"""plot colour-filled contours"""
lev = np.linspace(np.amin(F),np.amax(F),10)
#bar = pl.contourf(X,Y,F,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)



meanplot = np.reshape(mean,(N,N))                   #reshape mean as 2D array
#pl.contour(X,Y,meanplot, levels=[0,1,2,3,4])       #plot graph of mean
"""plot colour-filled contours"""
lev = np.linspace(np.amin(meanplot),np.amax(meanplot),10)
#bar = pl.contourf(X,Y,meanplot,levels=lev,cmap=pl.cm.seismic)
#pl.colorbar(bar)


F1d = np.reshape(F,N**2)
change = mean
for m in range(1,6):
    change += np.inner(F1d,vecs[:,-m])*vecs[:,-m]   #adding projections (of the parameterised form along the eigenvectors) to the mean
H = np.reshape(change,(N,N))
#pl.contour(X,Y,H, levels=[0,1,2,3,4])              #plot graph of 'change' on same graph - these are the points on the MoI ellipse that are closest to the parameterised form
"""plot colour-filled contours"""
lev = np.linspace(np.amin(H),np.amax(H),10)
bar = pl.contourf(X,Y,H,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)


pl.axes().set_aspect('equal')
pl.show()        

K = H-F
#pl.contour(X,Y,H, levels=[0,1,2,3,4])            #plot graph of 'change' on same graph - these are the points on the MoI ellipse that are closest to the parameterised form
"""plot colour-filled contours"""
lmax = np.amax(abs(K))
lev = np.linspace(-lmax,lmax,10)
bar = pl.contourf(X,Y,K,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)
pl.axes().set_aspect('equal')
pl.show()        

