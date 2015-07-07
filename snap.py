import numpy as np
import matplotlib.pyplot as pl

x = np.linspace(-5,5,20)
X,Y = np.meshgrid(x,x)

F = (X-2)**2 - (Y-3)**2
lev = np.linspace(np.amin(F),np.amax(F),10)

bar = pl.contourf(X,Y,F,levels=lev,cmap=pl.cm.seismic)
pl.colorbar(bar)

pl.axes().set_aspect('equal')
pl.show()
