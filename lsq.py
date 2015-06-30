# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 14:29:34 2015

@author: psaha
"""

import scipy.optimize as opt

def residuals(params):
    return (params[0]-25,params[1]+3,params[2]-1)

ini = [0,0,0]    
ans = opt.leastsq(residuals,ini)[0]
print(ans)


    