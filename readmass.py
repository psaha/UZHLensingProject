#readmass

#import libraries
import pickle
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize as opt

#open data file
mname = 'ASW0007k4r/012771'
mname = 'ASW0000h2m/007022'
mname = 'gribbles'
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


"""Alternative method for calculating mean"""

for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))
    mean_m = np.mean(ensem1d[m])
    if m==0:
        mean1 = np.array([mean_m])
    else:
        mean1 = np.append([mean1], [mean_m])

"""Why doesn't it work?"""

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
outsum = outsum/len(ensem)                  #scale MoI by tensor size of ensemble - makes no difference...
vals,vecs = np.linalg.eigh(outsum)          #find eigenvecs/vals of MoI tensor

#An experiment, no longer used:
#change = (vals[-1]**0.5)*vecs[:,-1] + mean  #<k> + sqrt(val)*vec for largest val (explain physics again?)
#plotchange = np.reshape(change,(N,N))       #reshape as 2D array for plotting


def profile(params):
    a,n=params[0],params[1]#,params[2],params[3]
    return a*(1+X*X+Y*Y)**-n                #define test parametrized functional form k(param)

from ellip import profile, grid

#print 'The functional form used is %s with parameters a and n' % model #prints k(param)

def residuals(params):
    f = profile(params)                     #f = k(param)
    f = np.reshape(f,(N**2))                #reshape f into 1D array
    f -= mean                             #changed 'change' (which was an experiment) to 'mean'#f = f-change
#    df = 5*[0]
    for m in range(1,6):
#        df[m-1] = np.inner(f,vecs[:,-m])/vals[-m]      #chi-squared attempt
        f -= np.inner(f,vecs[:,-m])*vecs[:,-m]          #removing projections along principle axes
    return f


ini = [3,0.5,0]                               #initial values for parameters
lsq = opt.leastsq(residuals,ini)[0]         #perform least squares fit on f
#print(lsq)

param1 = lsq[0]
param2 = lsq[1]
param3 = lsq[2]
#maybe user can input required no of sf?
print 'Einstein radius = {0:.3e}, Ellipticity = {1:.3e}, ell_pa (?projection angle) = {2:.3e}'.format(param1,param2,param3) #prints values of optimised parameters

G = grid(lsq)                               #For plotting grav ptl

F = profile(lsq)                            #F = k(param) with optimised parameters

F1d = np.reshape(F,N**2)

# F = np.reshape(mean,(N,N))
lev = np.linspace(np.amin(G),np.amax(G),21) #alternative graph plotting, currently unused
#pl.contour(X,Y,G, levels=lev)              #plot of gravitational potential
pl.contour(X,Y,F, levels=[0,1,2,3,4])       #plot graph of parametrized model

change = mean
for m in range(1,6):
    change += np.inner(F1d,vecs[:,-m])*vecs[:,-m] #adding projections (of the parameterised form along the eigenvectors) to the mean

H = np.reshape(change,(N,N))
pl.contour(X,Y,H, levels=[0,1,2,3,4])       #plot graph of 'change' on same graph - these are the points on the MoI ellipse that are closest to the parameterised form
pl.axes().set_aspect('equal')
#pl.xlabel('x axis')                        #some experiments with graph formatting
#pl.ylabel('y axis')
#pl.title('Title')
#pl.text(0,0,r'Text')
#pl.xlim(-20,20)
#pl.ylim(-20,20)
#pl.grid(True)
pl.show()        

#pl.savefig() for graph
