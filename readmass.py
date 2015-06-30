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
F = ensem[0]

G = np.reshape(F,(N**2))

H = np.reshape(G,(N,N))

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
##outsummod = outsum/len(ensem)
vals,vecs = np.linalg.eigh(outsum)

print vals

change = (vals[-1]**0.5)*vecs[:,-1] + mean
plotchange = np.reshape(change,(N,N))

lev = np.linspace(0,10,21)
pl.contour(X,Y,plotchange, levels=[0,1,2,3])
pl.axes().set_aspect('equal')
pl.show()

