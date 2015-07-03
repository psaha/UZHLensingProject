#testfile: parameterisation with good fit but runtime warning

#import libraries
import pickle
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize as opt

#open data file
mname = 'ASW0007k4r/012771'
mname = 'ASW0000h2m/007022'
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
plotmean = np.reshape(mean,(N,N))           #reshape mean as 2D array

#lev = np.linspace(0,10,21)
#pl.contour(X,Y,plotmean, levels=[0,1,2,3])
#pl.axes().set_aspect('equal')
#pl.show()                                  #plot graph of mean

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

change = (vals[-1]**0.5)*vecs[:,-1] + mean  #<k> + sqrt(val)*vec for largest val (explain physics again?)
plotchange = np.reshape(change,(N,N))       #reshape as 2D array for plotting

# print vals

def profile(params):
    a,b,n,h,c=params[0],params[1],params[2],params[3],params[4]
    #return a*(h+X*X+Y*Y+b*X*Y)**-n                #define test parametrized functional form k(param)
    return a*(h+c*X*X+Y*Y+b*X*Y)**-n

def residuals(params):
    f = profile(params)                     #f = k(param)
    f = np.reshape(f,(N**2))                #reshape f into 1D array
    f -= change                             #f = f-change
#    df = 5*[0]
    for m in range(1,6):
#        df[m-1] = np.inner(f,vecs[:,-m])/vals[-m]      #chi-squared attempt
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]          #removing projections along principle axes
    return f


ini = [1,-1.3,3.7,1.6,1]                               #initial values for parameters
lsq = opt.leastsq(residuals,ini)[0]         #perform least squares fit on f
#print(lsq)

param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
param4 = lsq[3]
param5 = lsq[4]

print 'param1 = {0:.3e}, param2 = {1:.3e}, param3 = {2:.3e}, param4 = {3:.3e}, param5 = {4:.3e}'.format(param1,param2,param3,param4,param5) #prints values of optimised parameters



F = profile(lsq)                            #profile?
# F = np.reshape(mean,(N,N))
lev = np.linspace(np.amin(F),np.amax(F),11) 
pl.contour(X,Y,F, levels=lev)       #plot graph of parametrized model
F = np.reshape(change,(N,N))
pl.contour(X,Y,F, levels=lev)       #plot graph of <k> + sqrt(val)*vec for largest val on same graph
pl.axes().set_aspect('equal')
pl.show()                                   

#pl.savefig() for graph
