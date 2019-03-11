#~ % Boyd & Vandenberghe, "Convex Optimization"
#~ % Joëlle Skaf - 08/23/05
#~ %
#~ % Solved a QCQP with 3 inequalities:
#~ %           minimize    1/2 x'*P0*x + q0'*r + r0
#~ %               s.t.    1/2 x'*Pi*x + qi'*r + ri <= 0   for i=1,2,3
#~ % and verifies that strong duality holds.

import cvxpy as cp
from numpy.random import randn
from numpy import eye
eps =0.000001


degree = 3
n = degree
nvars = degree+1
dim = 3
rdim = 2

def one_step(Pi, Pip1, V, zero_start = False, zero_end = False, x_start = None): #TODO: add V and A    
    x  = [cp.Variable(rdim) for i in range(nvars)] 
    ti = cp.Variable(1)
    
    P    = Pi  [0][:,:rdim]       ; p  = Pi  [1]
    P1   = Pip1[0][:,:rdim]       ; p1 = Pip1[1]
    Vn   = V[0][:,:rdim]*n        ; v = V[1]
    
    constraints = [
    #~ ti >= 0.0,
    ti == 1.,
    P*(x[0]) <= p,
    P*(x[3]) <= p,
    P*(x[0] + x[1]) <= p,
    P*(x[3] - x[2]) <= p,
    #end point constraint    
    P1*x[3] <= p1,
    #~ ]
    #velocities
    Vn*x[1] <= v * ti,
    Vn*x[2] <= v * ti,
    Vn*(x[3]-x[2]-x[1]-x[0]) <= v * ti,
    #accelerations
    #~ cp.quad_form(ti,ones((1,1)))*a + Ann1*(x[3]-x[2]-2*x[1]-x[0]) <= 0.
    ]
    
    #zero init vel
    if zero_start:
        constraints = constraints + [x[1] == zeros(rdim)]
    if zero_end:
        constraints = constraints + [x[2] == zeros(rdim)]
        
    if x_start is not None:
        constraints = constraints + [x[0] == x_start[:rdim]]
    #~ if x_end is not None:
        #~ constraints = constraints + [x[0] == x_end
    
    obj = cp.Minimize(ti)
    prob = cp.Problem(obj, constraints)
    prob.solve()
    return prob, x, ti[0]
    
    


from hpp_spline import bezier

def bezierFromVal(xis, ti):
    wps = zeros((len(xis),dim))
    wps[0,:rdim] = array(xis[0].value)
    wps[1,:rdim] = wps[0,:rdim] + array(xis[1].value)
    wps[3,:rdim] = array(xis[3].value)
    wps[2,:rdim] = wps[3,:rdim] - array(xis[2].value )
    return bezier(wps.transpose(), ti)

if __name__ == '__main__':
    from  cord_methods import *
    from  plot_cord import *
    pDef, inequalities_per_phase = genProblemDef(nvars,2)
    
    def oneVelIneq():
        v = ones(4); #v[:-3] = -v[3:]
        V = zeros((4,2)); 
        V[:2,:] = identity(2); 
        V[2:,:] =-identity(2); 
        return (V, v)
        
    def one():
        V = oneVelIneq()
        Pi   = inequalities_per_phase[0]
        Pip1 = inequalities_per_phase[1]
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.start.reshape((-1,))

        prob, xis, ti = one_step(Pi, Pip1, V, zero_start = True, zero_end = True, x_start=x_start)
        b = bezierFromVal(xis, abs(ti.value))
        plotBezier(b, "b", label = None, linewidth = 3.0)
        plotControlPoints(b, "b",linewidth=2)
        
        print("status:", prob.status)
        print("optimal value", ti.value)
        
        
        plt.plot(x_start[0],x_start[1],"g",linewidth=12)
        plt.show()
        return b
        
    b = one()
    
