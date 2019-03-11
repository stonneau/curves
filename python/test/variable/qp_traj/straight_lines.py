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

#find the minimum time sequence of straight lines connecting from start to end position along union of convex sets

def one_step(Pis, V, x_start = None, x_end = None): #TODO: add V and A    
    num_phases = len(Pis)
    num_vars = nvars * num_phases
    x   = [cp.Variable(rdim) for i in range(num_vars)] 
    tis = [cp.Variable(1) for i in range(num_phases)]
    
    constraints = []
    
    for j in range(len(Pis)):
        idx = j*nvars
        Pi = Pis[j]
        print "PI ", Pi  [1]
        P    = Pi[0][:,:rdim]       ; p  = Pi  [1]
        Vn   = V[0][:,:rdim]*n        ; v = V[1]
        ti = tis[j]
        constraints = constraints + [
        ti >= 0.0001,
        #positions
        P*(x[idx]) <= p,
        P*(x[idx+3]) <= p,
        P*(x[idx] + x[idx+1]) <= p,
        P*(x[idx+3] - x[idx+2]) <= p,
        #zero start/end velocities
        x[idx+1] == zeros(rdim),
        x[idx+2] == zeros(rdim),
        #velocities
        Vn*x[idx+1] <= v * ti,
        Vn*x[idx+2] <= v * ti,
        Vn*(x[idx+3]-x[idx+2]-x[idx+1]-x[idx]) <= v * ti,
        ]
        
        print "idx3", idx+3
        
    #continuity constraints
    for j in range(1,len(Pis)):
        idx = j*nvars
        constraints = constraints + [
        x[idx] == x[idx-1] #just positions because velocities at 0
        ]
            
    if x_start is not None:
        constraints = constraints + [x[0]  == x_start[:rdim]]
    #~ if x_end is not None:
        #~ constraints = constraints + [x[-1] == x_end[:rdim]]
    #~ if x_end is not None:
        #~ constraints = constraints + [x[0] == x_end
    
    obj = cp.Minimize(ti)
    prob = cp.Problem(obj, constraints)
    prob.solve(verbose=False)
    return prob, x, tis
    
    


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
    
    def oneVelIneq():
        v = ones(4); #v[:-3] = -v[3:]
        V = zeros((4,2)); 
        V[:2,:] = identity(2); 
        V[:2,:] =-identity(2); 
        return (V, v)
        
    def one(nphase=7):
        pDef, inequalities_per_phase = genProblemDef(nvars,nphase)
        V = oneVelIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.start.reshape((-1,))

        prob, xis, tis = one_step(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end)
        for i in range(len(tis)):
            print "tis", tis[i]
            ti = abs(tis[i].value)[0]
            b = bezierFromVal(xis[i*4:i*4+4], abs(ti))
            plotBezier(b, colors[i], label = None, linewidth = 3.0)
            plotControlPoints(b, colors[i],linewidth=2)
            print("optimal value", ti)
        
        print("status:", prob.status)
        
        
        plt.plot(x_start[0],x_start[1],"g",linewidth=12)
        plt.plot(x_end[0],x_end[1],"g",linewidth=12)
        plt.show()
        return b
        
    b = one()
    
