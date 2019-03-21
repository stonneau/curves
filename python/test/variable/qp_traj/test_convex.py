from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np



colors=['b', 'g', 'r', 'c', 'm', 'y', 'k']
labels=['Distance', 'Velocity', 'Acceleration', 'Jerk']

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

rangeval = 40
rg2 = rangeval / 2
delta = 0.5

#~ xx = []; tt = []
#~ for i in range(rangeval):
    #~ for j in range(rangeval):
        #~ xx = xx + [(i-rg2) *delta ]
        #~ tt = tt + [j *delta] # t always positive

#~ zz = [ xx[i]*xx[i] + 0.2* xx[i] + tt[i] * tt[i] for i in range(len(xx)) ]

t=[]; T=[]; z = []; zp = []; zzp =[]
for i in range(rangeval):
    for j in range(rangeval):
        ti = i *delta
        Ti = j*delta
        if(ti*ti >= Ti):
            t = t + [ti]
            T = T + [Ti]
            z = z + [- np.log(ti*ti - Ti)]

#~ alpha=[]; beta=[]; hess = []
#~ ti = 1
#~ T = 1.
#~ for i in range(rangeval):
    #~ for j in range(rangeval):
        #~ ai = (i) *delta
        #~ bi = (j) *delta
        #~ hess = hess + [2*ai*ai*ti*ti + ai*ai * T + bi*bi - 4 * ti * bi* ai]
        #~ alpha = alpha + [ai]
        #~ beta  = beta  + [bi]
#~ zz_1 = [el for el in zz if el <=1.]

ax.set_xlabel('t')
ax.set_ylabel('T')
ax.set_zlabel('hess')
#~ ax.set_xlabel('a')
#~ ax.set_ylabel('b')
#~ ax.set_zlabel('z')

ax.scatter(t,T,z,color='r', linewidth = 0.01)
#~ ax.scatter(alpha,beta,hess,color='r', linewidth = 0.01)
#~ ax.scatter(xx,tt,zz_1,color='b', linewidth = 0.01)
#~ plt.show()

print "min z", min(z)
print "min hess", min(hess)

from __future__ import division
from sympy import *

from sympy.tensor.array import derive_by_array

t, T = symbols('t T')

f = -log(t*t - T)

list(derive_by_array(f, (t, T)))
hess = list(derive_by_array(derive_by_array(f, (t, T)), (t, T)))

m = Matrix(( [hess[0] , hess[1]], [hess[2], hess[3]]  ))

a = simplify(m.det())

