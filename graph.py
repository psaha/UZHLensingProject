"""Plots graphs comparing input and output parameters of simulated lenses"""

import numpy as np
import matplotlib.pyplot as pl

fil = open('compare_params.txt')
lines = fil.readlines()                             #Read data from file

#Define x and y coordinates for einstein radius, ellipticity and angle of ellipticity
xeinst = []
xellip = []
xangle = []
yeinst = []
yellip = []
yangle = []

xeinstd = []                                        #Separate out the data for doubles
xellipd = []
xangled = []
yeinstd = []
yellipd = []
yangled = []

for i in lines:
    if 'Original ' in i:
        orig = i[28:].split()
        xeinst.append(float(orig[0]))
        xellip.append(float(orig[1]))
        xangle.append(float(orig[2]))
    if 'pixrad 12 lsq' in i:
        param = i[28:].split()
        yeinst.append(float(param[0]))
        yellip.append(float(param[1]))
        yangle.append(float(param[2]))
    if 'Originald' in i:                            #Separate out the data for doubles
        origd = i[29:].split()
        xeinstd.append(float(origd[0]))
        xellipd.append(float(origd[1]))
        xangled.append(float(origd[2]))
    if 'pixradd' in i:
        paramd = i[29:].split()
        yeinstd.append(float(paramd[0]))
        yellipd.append(float(paramd[1]))
        yangled.append(float(paramd[2]))

    
"""Plot einstein radii"""
#pl.plot(xeinst,yeinst,'.')
pl.xlim(0,2)
pl.ylim(0,2)
pl.plot(xeinst,yeinst,'b.',label="not doubles")
pl.plot(xeinstd,yeinstd,'r.',label="doubles")
#pl.legend()
pl.plot([0, 2], [0, 2], 'b--')
pl.axes().set_aspect('equal')
pl.title('Einstein radius comparison')
pl.xlabel('Original parameters')
pl.ylabel('Output parameters')
pl.grid(True)
pl.show()

"""Plot ellipticity and angles"""
for n in range(len(xellip)):
    R = xellip[n]
    r = yellip[n]
    T = (xangle[n]+90)*np.pi/180
    t = (yangle[n]+90)*np.pi/180
    x = [R*np.cos(T),r*np.cos(t)]
    y = [R*np.sin(T),r*np.sin(t)]
    pl.plot(x[:-1],y[:-1],'.')                      #True parameters
    #pl.plot(x[1],y[1],'*')
    pl.plot(x,y,'b')
"""doubles"""
for n in range(len(xellipd)):
    R = xellipd[n]
    r = yellipd[n]
    T = (xangled[n]+90)*np.pi/180
    t = (yangled[n]+90)*np.pi/180
    x = [R*np.cos(T),r*np.cos(t)]
    y = [R*np.sin(T),r*np.sin(t)]
    pl.plot(x[:-1],y[:-1],'.')                      #True parameters
    #pl.plot(x[1],y[1],'*')
    pl.plot(x,y,'r')
#pl.plot(xellip,yellip,'.')
pl.xlim(-0.5,0.5)
pl.ylim(-0.5,0.5)
pl.xlabel('ellip*cos(angle)')
pl.ylabel('ellip*sin(angle)')
pl.axes().set_aspect('equal')
pl.show()


    
xangle = np.array(xangle)
yangle = np.array(yangle)
diff = yangle - xangle
    
xangled = np.array(xangled)
yangled = np.array(yangled)
diffd = yangled - xangled

pl.plot(np.cos(diff),np.sin(diff),'b.')
pl.plot(np.cos(diffd),np.sin(diffd),'r.')
#pl.xlim(0,2)
#pl.ylim(0,2)
pl.grid(True)
pl.axes().set_aspect('equal')
pl.show()

