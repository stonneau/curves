from __future__ import division
from sympy import *

from sympy.tensor.array import derive_by_array

t, T = symbols('t T')
t1, T2 = symbols('t T')

#~ f = -log(t*t - T)
#~ f = -log(1- t*t - T)
f = -log(1+ t*t - T)
#~ a  = 1.
#~ f = a * exp(a*(t*t-T))
#~ f = -1/(t*t - T)
#~ f = log(t*t - T)


hess = list(derive_by_array(derive_by_array(f, (t, T)), (t, T)))

m = Matrix(( [hess[0] , hess[1]], [hess[2], hess[3]]  ))

print "hessian"
print simplify(m)

a = simplify(m.det())

print 'hessian neg', simplify(a)

import numpy as np

def fu(t,T):
    assert (t*t > T)
    #~ return -np.log(t)
    #~ return np.exp(t*t-T)
    return -np.log(t*t-T)
    #~ print -1/(t*t-T)
    #~ return -1/(t*t-T)

from random import uniform

def randomval():
    t = uniform(0,100)
    while True:
        T = uniform(np.sqrt(2)+0.001 ,100)
        if(t*t > T):
            return (t,T)
            
#~ for i in range(1000):
    #~ (t1,T1) = randomval()
    #~ (t2,T2) = randomval()
    #~ for j in range(100):
        #~ alpha = uniform(0,1)
        #~ nt = alpha * t1 + (1-alpha) * t2
        #~ nT = alpha * T1 + (1-alpha) * T2
        #~ if(nt*nt > nT):
            #~ if (fu(nt,nT) > max(fu(t1,T1),  fu(t2,T2) ) + 0.0001):
                #~ print "failed"
                #~ print "int", fu(nt,nT)
                #~ print "t1", fu(t1,T1)
                #~ print "t2", fu(t2,T2) 
                #~ print "t", nt
                #~ print "T", nT 
                #~ assert(False)

#~ print "ok"

#this is negative since T < ti**2
from sympy.plotting import plot3d
#~ plot3d(f)
#~ plot3d(f, (t,0,2),(T,4,10))
plot3d(f, (t,0,11),(T,np.sqrt(2),8))
#~ plot3d(f, (t,0,1000),(T,0,1000))
#~ plot3d(-1 / (0.001 * (t*t - T)), (t,0.001,10),(T,0.001,10))
#~ plot3d(exp(t*t + T)exp(t*t + T), (t,0.001,10),(T,0.001,10))

#~ plot3d(-f)
