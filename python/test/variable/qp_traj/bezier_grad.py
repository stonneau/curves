from __future__ import division
from sympy import *

from sympy.tensor.array import derive_by_array

xs0, x0, xg0, t = symbols('xs0 x0 xg0 t')
xs1, x1, xg1    = symbols('xs1 x1 xg1')
#~ y0, y1, t = symbols('y0 y1 t')

init_printing(use_unicode=True)

from scipy.special import binom as bino

npts = 3
pis = [xs0,x0,xg0]
deg = npts -1

factors = [bino(deg,i) * t**i * (1-t)**(deg-i)*pis[i] for i in range(len(pis))]
eq = sum(factors);

f1 = (eq*eq - 1)**2


grad = derive_by_array(f1, (x0, x1))
#~ grad = derive_by_array(f2, (x0, y0,x1,y1))
#~ grad = derive_by_array(f2, (x0, y0,x1,y1))
    

#~ grad = derive_by_array(F, (x0, x1)).reshape(2,2)
#~ hess = derive_by_array(grad, (x0, x1))
#~ hess = derive_by_array(grad, (x0, y0,x1,y1))

#~ m = Matrix(( [hess[0] , hess[1]], [hess[2], hess[3]]  ))


print "gradient f1"
print grad

#~ print "hessian f1"
#~ print simplify(hess)


from sympy.plotting import plot3d

