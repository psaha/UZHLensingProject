#readmass
#small edit to test branching
#second branching test

#import libraries
import pickle
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize as opt

#open data file
mname = 'ASW0007k4r/012771'
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
    if m==0:                                #create MoI tensor (outer products of pairs of values)
        outsum = out
    else:
        outsum = outsum + out
outsum = outsum/len(ensem)                  #scale MoI by tensor size of ensemble
vals,vecs = np.linalg.eigh(outsum)          #find eigenvecs/vals of MoI tensor

change = (vals[-1]**0.5)*vecs[:,-1] + mean  #<k> + sqrt(val)*vec for largest val (explain physics again?)
plotchange = np.reshape(change,(N,N))       #reshape as 2D array for plotting

# print vals

def profile(params):
    a,n=params[0],params[1]#,params[2],params[3]
    return a*(1+X*X+Y*Y)**-n                #define test parametrized functional form k(param)

def residuals(params):
    f = profile(params)                     #f = k(param)
    f = np.reshape(f,(N**2))                #reshape f into 1D array
    f -= mean                             #changed 'change' (which was an experiment) to 'mean', consider changing this in all scripts or at least think about it#f = f-change
#    df = 5*[0]
    for m in range(1,6):
#        df[m-1] = np.inner(f,vecs[:,-m])/vals[-m]      #chi-squared attempt
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]          #?? sthg to do with delta(k(param))
    return f


ini = [3,0.5]                               #initial values for parameters
lsq = opt.leastsq(residuals,ini)[0]         #perform least squares fit on f
print(lsq)



F = profile(lsq)
# F = np.reshape(mean,(N,N))
lev = np.linspace(np.amin(F),np.amax(F),21) 
pl.contour(X,Y,F, levels=[0,1,2,3,4])       #plot graph of parametrized model
F = np.reshape(change,(N,N))
pl.contour(X,Y,F, levels=[0,1,2,3,4])       #plot graph of <k> + sqrt(val)*vec for largest val on same graph
pl.axes().set_aspect('equal')
pl.show()                                   

#pl.savefig() for graph
