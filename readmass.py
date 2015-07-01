import pickle
import numpy as np
import matplotlib.pyplot as pl

mname = 'ASW0007k4r/012771'
fil = open(mname+'.pkl')
ensem = pickle.load(fil)

N = ensem[0].shape[0]
R = (N-1)/2
x = np.linspace(-R,R,N)
X,Y = np.meshgrid(x,x)



for m in range(len(ensem)):
    model1d = np.reshape(ensem[m],(N**2))    
    if m==0:    
        sum = model1d
    else:    
        sum = sum+model1d
mean = sum/len(ensem)
plotmean = np.reshape(mean,(N,N))

#lev = np.linspace(0,10,21)
#pl.contour(X,Y,plotmean, levels=[0,1,2,3])
#pl.axes().set_aspect('equal')
#pl.show()

for m in range(len(ensem)):
    model1d = np.reshape(ensem[m],(N**2))
    diff = model1d - mean
    out = np.outer(diff,diff)
    if m==0:
        outsum = out
    else:
        outsum = outsum + out
outsum = outsum/len(ensem)
vals,vecs = np.linalg.eigh(outsum)

change = (vals[-1]**0.5)*vecs[:,-1] + mean
plotchange = np.reshape(change,(N,N))

# print vals

def profile(params):
    a,n=params[0],params[1]#,params[2],params[3]
    return a*(1+X*X+Y*Y)**-n

def residuals(params):
    f = profile(params)
    f = np.reshape(f,(N**2))
    f -= change
#    df = 5*[0]
    for m in range(1,6):
#        df[m-1] = np.inner(f,vecs[:,-m])/vals[-m]
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]
    return f

import scipy.optimize as opt
ini = [3,0.5]
lsq = opt.leastsq(residuals,ini)[0]
print(lsq)



F = profile(lsq)
# F = np.reshape(mean,(N,N))
lev = np.linspace(np.amin(F),np.amax(F),21)
pl.contour(X,Y,F, levels=[0,1,2,3,4])
F = np.reshape(change,(N,N))
pl.contour(X,Y,F, levels=[0,1,2,3,4])
pl.axes().set_aspect('equal')
pl.show()

