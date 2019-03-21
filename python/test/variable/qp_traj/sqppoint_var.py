
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array
from numpy.linalg import norm
eps =0.000001


degree = 4
n = degree
nvars = degree+1
dim = 3
rdim = 2

#find the minimum time sequence of straight lines connecting from start to end position along union of convex sets

def solve_straight_lines(Pis, V, x_start = None, x_end = None, rafacc = True): #TODO: add V and A    
    num_phases = len(Pis)
    num_vars = nvars * num_phases
    x   = [cp.Variable(rdim) for i in range(num_vars)] 
    tis = [cp.Variable(1) for i in range(num_phases)]
    
    constraints = [ti >= 0.00 for ti in tis]
    
    Vn   = V[0][:,:rdim]*n ; v = V[1]
    A    = Vn*(n-1)        ; a = V[1]
    
    for j in range(len(Pis)):
        idx = j*nvars
        pos = x[idx:idx+nvars]
        velt = [pos[k+1]-pos[k] for k in range(nvars-1)] 
        acctt = [velt[k+1]-velt[k] for k in range(nvars-2)] 
        P    = Pis[j][0][:,:rdim]       ; p  = Pis[j][1]
        ti = tis[j]
        
        #positions
        constraints = constraints + [P *el <= p    for el in pos]
        #velocities
        constraints = constraints + [Vn*el <= v*ti for el in velt]
        #velocities
        #~ constraints = constraints + [A*el <= a*ti for el in acctt]
        
        if (j>0):            
            P_prev    = Pis[j-1][0][:,:rdim]; p_prev  = Pis[j-1]  [1] 
            #position constraints           
            constraints = constraints + [
            pos[0] == x[idx-1] ]
            #assuming same time for crossing previous polytope 
            #soft continuity velocity: previous waypoints is in previous polytope
            P_n1prev = pos[0] - velt[0]
            P_n2prev = acctt[0] + 2 * P_n1prev - pos[0]
            constraints = constraints + [
            P_prev*(P_n1prev) <= p_prev,
            P_prev*(P_n2prev) <= p_prev,
            ]
        
        if j < len(Pis)-1:
            P_next    = Pis[j+1][0][:,:rdim]; p_next  = Pis[j+1]  [1]
            #assuming same time for crossing previous polytope 
            #soft continuity velocity: previous waypoints is in previous polytope
            P1next = pos[-1] + velt[-1]
            P2next = acctt[-1] + 2 * P1next - pos[0]
            constraints = constraints + [
            P_next*(P1next) <= p_next,
            P_next*(P2next) <= p_next,
            ]
            
    if x_start is not None:
        constraints = constraints + [x[0]  == x_start[:rdim]]
    if x_end is not None:
        constraints = constraints + [x[-1]  == x_end[:rdim]]
        
    #zero vel/ end constraints:
    constraints = constraints + [x[0]  == x[1], x[0]  == x[2],x[-1]  == x[-2],x[-3]  == x[-1]]
        
    obj = cp.Minimize(sum(tis))
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS",verbose=True)
    #~ prob.solve(verbose=True)
    
    tis, xs = tovals(tis), tovals(x)
    
    for j in range(1,len(tis)):
        idx     = j*nvars
        pos     = xs[idx:idx+nvars]
        vel     = [n     / tis[j] * (pos[k+1] - pos[k]) for k in range(nvars-1) ]
        accs    = [(n-1) / tis[j] * (vel[k+1] - vel[k]) for k in range(nvars-2) ]
        idprev  = (j-1)*nvars
        posprev = xs[idprev:idprev+nvars]
        velprev = [n     / tis[j-1] * (posprev[k+1] - posprev[k]) for k in range(nvars-1) ]
        accprev = [(n-1) / tis[j-1] * (velprev[k+1] - velprev[k]) for k in range(nvars-2) ]
        #take whatever velocity is the fastest
        #~ if(norm(velprev[-1]) > norm(vel[0])):
        if(tis[j] < tis[j-1]):
            pos[1] = pos[0] + tis[j] * velprev[-1] / n
            #~ pos[2] = accprev[-1] / (n*(n-1)) * tis[j] * tis[j] +  2 * pos[1] - pos[0]
            xs[idx:idx+nvars] = pos
        else:
            posprev[-2] = posprev[-1] - tis[j-1] * vel[0] / n
        #~ if rafacc:
            #~ if(tis[j]**2 < tis[j-1]**2):
                #~ pos[2] = accprev[-1] / (n*(n-1)) * tis[j] * tis[j] +  2 * pos[1] - pos[0]
            #~ else:
                #~ posprev[-3] = accs[0] / (n*(n-1)) * tis[j-1] * tis[j-1] + 2 * posprev[-2] - posprev[-1]
        xs[idx:idx+nvars] = pos
        xs[idprev:idprev+nvars] = posprev
    
    return prob, xs, tis
    
    

def bezierFromVal(xis, ti):
    wps = zeros((len(xis),dim))
    for i in range(nvars):
        wps[i,:rdim] = array(xis[i])
    
    return bezier(wps.transpose(), ti)
    
from hpp_spline import bezier


def tovals(variables):
    return [el.value for el in variables]

    

from  cord_methods import *
from  plot_cord import *

def boundIneq():
    v = ones(4)*1.; #v[:-3] = -v[3:]
    V = zeros((4,2)); 
    V[:2,:] = identity(2); 
    V[2:,:] =-identity(2); 
    return (V, v)

if __name__ == '__main__':
    
        
    def one(nphase=2):
        plt.close()
        pDef, inequalities_per_phase = genProblemDef(6,nphase)
        V = boundIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))

        prob, xis, tis = solve_straight_lines(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end)
        
        #~ print "times", tis
        for i in range(len(tis)):
            ti = abs(tis[i][0])
            b = bezierFromVal(xis[i*nvars:i*nvars+nvars], abs(ti))
            plotBezier(b, colors[i], label = None, linewidth = 3.0)
            #~ plotBezier(b, "r", label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[i],linewidth=2)
            
        
        #~ prob, xis, tis = solve_straight_lines(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, rafacc=False)
        
        #~ print "times", tis
        #~ for i in range(len(tis)):
            #~ ti = abs(tis[i][0])
            #~ b = bezierFromVal(xis[i*nvars:i*nvars+nvars], abs(ti))
            #~ plotBezier(b, colors[i], label = None, linewidth = 3.0)
            #~ plotBezier(b, "b", label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[i],linewidth=2)
            
            
        
        
        plt.show()
        
        #~ plt.show()
        return b, tis
        
    b, tis = one(2)
    
