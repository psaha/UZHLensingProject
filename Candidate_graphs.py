import numpy as np
import matplotlib.pyplot as pl

fil = open('candidate_params_angles.txt')
lines = fil.readlines() 

"""Plot histograms"""
einst = []

for i in lines:
    if 'ASW' in i:
        results = i[23:].split()
        einst.append(float(results[0]))
        
einst.sort()
pl.hist(einst)
pl.xlabel('Einstein radius/arcseconds')
pl.title('Distribution of Einstein radii')
pl.xlim(0,4)
pl.show()

ellip = []

for i in lines:
    if 'ASW' in i:
        results = i[23:].split()
        ellip.append(float(results[1]))
        
ellip.sort()
pl.hist(ellip)
pl.xlabel('Ellipticity/units')
pl.title('Distribution of Ellipticities')
#pl.xlim(0,4)
pl.show()

angle = []

for i in lines:
    if 'ASW' in i:
        results = i[23:].split()
        angle.append(float(results[2]))
        
angle.sort()
pl.hist(angle)
pl.xlabel('Angle of ellipticity/degrees')
pl.title('Distribution of angles of ellipticity')
pl.xlim(-90,90)
pl.show()


"""Plot quads and doubles as separate points"""
einst = []
ellip = []
angle = []
einstd = []
ellipd = []
angled = []

for i in lines:
    if 'DOUBLE' in i:
        resultsd = i[23:].split()
        einstd.append(float(resultsd[0]))
        ellipd.append(float(resultsd[1]))
        angled.append(float(resultsd[2]))
    elif 'ASW' in i:
        results = i[23:].split()
        einst.append(float(results[0]))
        ellip.append(float(results[1]))
        angle.append(float(results[2]))

einstd.sort()
ellipd.sort()
angled.sort()
einst.sort()
ellip.sort()
angle.sort()


for i in einst:
    pl.plot(i,0,'b.')
for i in einstd:
    pl.plot(i,0,'r.')
pl.ylim(-0.01,0.01)
#ax.set_yticklabels([])
#pl.show()

