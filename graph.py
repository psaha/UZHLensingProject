import numpy as np
import matplotlib.pyplot as pl

fil = open('compare_params.txt')
lines = fil.readlines()

xeinst = []
xellip = []
xangle = []
yeinst = []
yellip = []
yangle = []
for i in lines:
    if 'Original' in i:
        orig = i[28:].split()
        xeinst.append(float(orig[0]))
        xellip.append(float(orig[1]))
        xangle.append(float(orig[2]))
    if 'pixrad 12 lsq' in i:
        print i
        param = i[28:].split()
        yeinst.append(float(param[0]))
        yellip.append(float(param[1]))
        yangle.append(float(param[2]))
    
xangle = np.array(xangle)
yangle = np.array(yangle)
diff = yangle - xangle
    
pl.plot(xeinst,yeinst,'.')
pl.xlim(0,2)
pl.ylim(0,2)
pl.axes().set_aspect('equal')
pl.show()

for n in range(len(xellip)):
    R = xellip[n]
    r = yellip[n]
    T = (xangle[n]+90)*np.pi/180
    t = (yangle[n]+90)*np.pi/180
    x = [R*np.cos(T),r*np.cos(t)]
    y = [R*np.sin(T),r*np.sin(t)]
    pl.plot(x[:-1],y[:-1],'.')
    pl.plot(x,y)
#pl.plot(xellip,yellip,'.')
#pl.xlim(0,0.5)
#pl.ylim(0,0.5)
pl.axes().set_aspect('equal')
pl.show()

pl.plot(np.cos(diff),np.sin(diff),'.')
#pl.xlim(0,2)
#pl.ylim(0,2)
pl.grid(True)
pl.axes().set_aspect('equal')
pl.show()
    
#pl.xlabel('x axis')                        #some experiments with graph formatting
#pl.ylabel('y axis')
#pl.title('Title')
#pl.text(0,0,r'Text')
#pl.xlim(-20,20)
#pl.ylim(-20,20)
#pl.grid(True)
#pl.savefig() for graph