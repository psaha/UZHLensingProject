"""Extra code from readmass moved here for tidying-up purposes"""

"""Alternative method for calculating mean"""

for m in range(len(ensem)):
    ensem1d = np.reshape(ensem[m],(N**2))
    mean_m = np.mean(ensem1d[m])
    if m==0:
        mean1 = np.array([mean_m])
    else:
        mean1 = np.append([mean1], [mean_m])

"""Calculates the wrong mean"""


"""Plot graph of mean"""
#lev = np.linspace(0,10,21)
#pl.contour(X,Y,plotmean, levels=[0,1,2,3])
#pl.axes().set_aspect('equal')
#pl.show()                                  #plot graph of mean


"""An experiment, no longer used:"""
#change = (vals[-1]**0.5)*vecs[:,-1] + mean  #<k> + sqrt(val)*vec for largest val (explain physics again?)
#plotchange = np.reshape(change,(N,N))       #reshape as 2D array for plotting

"""Chi-squared remnants"""
def residuals(params):
    f = profile(params)                     #f = k(param)
    f = np.reshape(f,(N**2))                #reshape f into 1D array
    f -= mean                             #changed 'change' (which was an experiment) to 'mean'#f = f-change
    df = 5*[0]
    for m in range(1,6):
        df[m-1] = np.inner(f,vecs[:,-m])/vals[-m]      #chi-squared attempt
    return df


"""Graph axis labelling etc"""
#pl.xlabel('x axis')                        #some experiments with graph formatting
#pl.ylabel('y axis')
#pl.title('Title')
#pl.text(0,0,r'Text')
#pl.xlim(-20,20)
#pl.ylim(-20,20)
#pl.grid(True)
#pl.savefig() for graph


"""Stuff from ellip for graph plotting"""
#params = [10,0.5,45]
#K = profile(params)

#R = (N-1)/2
#x = np.linspace(-R,R,N)
#X,Y = np.meshgrid(x,x)

# M = np.log(abs(M)+1e-12)
#lev = np.linspace(np.amin(K),np.amax(K),21)
#pl.contour(X,Y,K, levels=[0,1,2,3,4,5])
#pl.axes().set_aspect('equal')
#pl.show()        






