"""Plots graphs comparing input and output parameters of simulated lenses"""

import numpy as np
import matplotlib.pyplot as pl

fil = open('Com_params_w_oldsims.txt')
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
yeinsto = []
yellipo = []
yangleo = []
yeinstod = []
yellipod = []
yangleod = []

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
        origd = i[28:].split()
        xeinstd.append(float(origd[0]))
        xellipd.append(float(origd[1]))
        xangled.append(float(origd[2]))
    if 'pixradd' in i:
        paramd = i[29:].split()
        yeinstd.append(float(paramd[0]))
        yellipd.append(float(paramd[1]))
        yangled.append(float(paramd[2]))
    if 'Old ' in i:
        old = i[28:].split()
        yeinsto.append(float(old[0]))
        yellipo.append(float(old[1]))
        yangleo.append(float(old[2]))      
    if 'Oldd' in i:
        oldd = i[29:].split()
        yeinstod.append(float(oldd[0]))
        yellipod.append(float(oldd[1]))
        yangleod.append(float(oldd[2]))    



"""Plot einstein radii"""
pl.xlim(0,2)
pl.ylim(0,2)
pl.plot(xeinst,yeinst,'b.',label="quads")
pl.plot(xeinstd,yeinstd,'r.',label="doubles")

#old
pl.plot(xeinst,yeinsto,'g*',label="quads, others' models")
pl.plot(xeinstd,yeinstod,'y*',label="doubles, others'models")

#pl.legend()
pl.plot([0, 2], [0, 2], 'b--')
pl.axes().set_aspect('equal')
pl.title('Einstein radius comparison')
pl.xlabel('Original parameters')
pl.ylabel('Output parameters')
pl.grid(True)
pl.show()


"""Plot ellipticity"""
pl.xlim(0,0.5)
pl.ylim(0,0.5)
pl.plot(xellip,yellip,'b.',label="not doubles")
pl.plot(xellipd,yellipd,'r.',label="doubles")

#old
pl.plot(xellip,yellipo,'g*')
pl.plot(xellipd,yellipod,'y*')

#pl.legend()
pl.plot([0, 0.5], [0, 0.5], 'b--')
pl.axes().set_aspect('equal')
pl.title('Ellipticity comparison')
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

#old
for n in range(len(xellip)):
    R = xellip[n]
    r = yellip[n]
    T = (xangle[n]+90)*np.pi/180
    t = (yangleo[n]+90)*np.pi/180
    x = [R*np.cos(T),r*np.cos(t)]
    y = [R*np.sin(T),r*np.sin(t)]
    pl.plot(x[:-1],y[:-1],'.')                      #True parameters
    #pl.plot(x[1],y[1],'*')
    pl.plot(x,y,'g')
#old doubles
for n in range(len(xellipd)):
    R = xellipd[n]
    r = yellipd[n]
    T = (xangled[n]+90)*np.pi/180
    t = (yangleod[n]+90)*np.pi/180
    x = [R*np.cos(T),r*np.cos(t)]
    y = [R*np.sin(T),r*np.sin(t)]
    pl.plot(x[:-1],y[:-1],'.')                      #True parameters
    #pl.plot(x[1],y[1],'*')
    pl.plot(x,y,'y')
   
#pl.plot(xellip,yellip,'.')
pl.xlim(-0.5,0.5)
pl.ylim(0,0.5)
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

"""
"""
#histograms (no doubles)
yeinst.sort()
pl.hist(yeinst)
pl.xlabel('Einstein radius/arcseconds')
pl.title('Distribution of Einstein radii')
pl.xlim(0.8,1.4)
pl.show()

yellip.sort()
pl.hist(yellip)
pl.xlabel('Ellipticity/units')
pl.title('Distribution of Ellipticities')
#pl.xlim(0,4)
pl.show()

yangle.sort()
pl.hist(yangle)
pl.xlabel('Angle of ellipticity/degrees')
pl.title('Distribution of angles of ellipticity')
pl.xlim(-90,90)
pl.show()

#histograms of old sims (no doubles)
yeinsto.sort()
pl.hist(yeinsto)
pl.xlabel('Einstein radius/arcseconds')
pl.title('Distribution of Einstein radii')
#pl.xlim(0.8,1.4)
pl.show()

yellipo.sort()
pl.hist(yellipo)
pl.xlabel('Ellipticity/units')
pl.title('Distribution of Ellipticities')
#pl.xlim(0,4)
pl.show()

yangleo.sort()
pl.hist(yangleo)
pl.xlabel('Angle of ellipticity/degrees')
pl.title('Distribution of angles of ellipticity')
pl.xlim(-90,90)
pl.show()
