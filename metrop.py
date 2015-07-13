from numpy.random import random, random_sample
from numpy import array, matrix, exp

def samp(fun,lo,hi,N):
    lo,hi = array(lo),array(hi)
    w = (lo+hi)/2.
    wlyst = [w+0]
    lnplyst = [fun(w)]
    n = 0
    while True:
        dw = 2*random_sample(len(w)) - 1
        w += (hi-lo) * dw/10
        for k in range(len(w)):
            if w[k] < lo[k]:
                w[k] += hi[k] - lo[k]
            if w[k] > hi[k]:
                w[k] -= hi[k] - lo[k]
        lnp = fun(w)
        if random() > exp(lnp-lnplyst[n%N]):
            w = wlyst[n%N] + 0
            lnp = lnplyst[n%N]
        n += 1
        if n%50 == 0:
            print n
        if n < N:
            wlyst += [w+0]
            lnplyst += [lnp]
        else:
            wlyst[n%N] = w + 0
            lnplyst[n%N] = lnp
            if lnplyst[n%N] < lnplyst[(n+1)%N]:
                return (array(lnplyst),array(wlyst))


