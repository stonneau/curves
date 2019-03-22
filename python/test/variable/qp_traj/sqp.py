
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array
eps =0.000001


degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2

#find the minimum time sequence of straight lines connecting from start to end position along union of convex sets

def solve_straight_lines(Pis, V, x_start = None, x_end = None): #TODO: add V and A    
    num_phases = len(Pis)
    num_vars = nvars * num_phases
    x   = [cp.Variable(rdim) for i in range(num_vars)] 
    tis = [cp.Variable(1) for i in range(num_phases)]
    Tis = [cp.Variable(1) for i in range(num_phases)]
    
    constraints = [Ti >= 0. for Ti in Tis]
    constraints = constraints + [ti >= 0.00 for ti in tis]
    
    for j in range(len(Pis)):
        idx = j*nvars
        Pi = Pis[j]
        P    = Pi[0][:,:rdim]       ; p  = Pi  [1]
        Vn   = V[0][:,:rdim]*n        ; v = V[1]
        A    = Vn*(n-1)        ; a = V[1]
        ti = tis[j]
        Ti = Tis[j]
        constraints = constraints + [
        #~ ti >= 0.0001,
        #positions
        P*(x[idx]) <= p,
        P*(x[idx+5]) <= p,
        P*(x[idx  ] + x[idx+1]) <= p,
        P*(x[idx+5] - x[idx+4]) <= p,
        P*(x[idx+2] +2*x[idx+1] + x[idx]) <= p,
        P*(x[idx+5] -2*x[idx+4] + x[idx+3] ) <= p,
        #zero start/end velocities
        #~ x[idx+1] == zeros(rdim),
        #~ x[idx+2] == zeros(rdim),
        #velocities
        Vn*x[idx+1] <= v * ti,
        Vn*x[idx+4] <= v * ti,
        Vn*(x[idx+2]+x[idx+1]) <= v * ti,
        Vn*(x[idx+4]-x[idx+3]) <= v * ti,
        Vn*(x[idx+5]-2*x[idx+4]+x[idx+3]-x[idx+2]-2*x[idx+1]-x[idx]) <= v * ti,
        #~ Vn*(x[idx+3]-x[idx+2]-x[idx+1]-x[idx]) <= v * ti,
        ]
        
        dc = [x[idx+1],x[idx+2]+x[idx+1],x[idx+5]-2*x[idx+4]+x[idx+3]-x[idx+2]-2*x[idx+1]-x[idx],x[idx+4]-x[idx+3],x[idx+1] ]
        ddc = [dc[k+1]-dc[k] for k in range(len(dc)-1)]
        constraints = constraints + [A * ddci <= a * Ti for ddci in ddc]
        
        
    #continuity constraints
    for j in range(0,len(Pis)):
        
        idx = j*nvars
        P    = Pis[j][0][:,:rdim]       ; p  = Pis[j][1]
        
        
        #velocity constraints
        P0 = x[idx]
        P1  = x[idx] + x[idx+1]
        P2  = x[idx+2] + 2* P1 - P0
        PN = x[idx+5]
        P_n1 = PN - x[idx+4]
        P_n2 =  x[idx+3] + 2 * P_n1  - PN
        c0t = x[idx+1]
        ddc0t = x[idx+2]
        ddc3t = x[idx+3]
        c4t = x[idx+4]
        
        if (j>0):            
            P_prev    = Pis[j-1][0][:,:rdim]; p_prev  = Pis[j-1]  [1] 
            #position constraints           
            constraints = constraints + [
            x[idx] == x[idx-1] ]
            #assuming same time for crossing previous polytope 
            #soft continuity velocity: previous waypoints is in previous polytope
            P_n1prev = P0 - c0t
            P_n2prev = ddc0t + 2 * P_n1prev - P0
            constraints = constraints + [
            P_prev*(P_n1prev) <= p_prev,
            P_prev*(P_n2prev) <= p_prev,
            ]
        
        if j < len(Pis)-1:
            P_next    = Pis[j+1][0][:,:rdim]; p_next  = Pis[j+1]  [1]
            #assuming same time for crossing previous polytope 
            #soft continuity velocity: previous waypoints is in previous polytope
            P1next = PN + c4t
            P2next = ddc3t + 2 * P1next - P0
            constraints = constraints + [
            P_next*(P1next) <= p_next,
            P_next*(P2next) <= p_next,
            ]
        
        
        
        
        
        
        
        
        
        #~ constraints = constraints + [
        #~ P*(x[idx] + x[idx-2]) <= p,      #velocity in next phase
        #~ Pm*(x[idx] - x[idx+1]) <= pm,      #velocity in previous phase
        #~ P*([])
        #~ P*( x[idx-3]+ 2* (x[idx] + x[idx+1]) -x[idx]) <= p,      #previous acceleration in next phase
        #~ Pm*(x[idx+2]+ 2* (x[idx-2] + x[idx]) -x[idx]) <= pm      #acceleration in previous phase
        #~ x[idx+1] == x[idx-2], #just positions because velocities at 0
        #~ x[idx+2] == x[idx-3] #just positions because velocities at 0
        #~ ,x[idx+3] == x[idx-4] #just positions because velocities at 0 ?????
        #~ ,x[idx+4] == x[idx-5] #just positions because velocities at 0
        #~ ]
            
    if x_start is not None:
        constraints = constraints + [x[0]  == x_start[:rdim]]
        constraints = constraints + [x[1]  == zeros(rdim)]
        constraints = constraints + [x[1]  == x[2]]
    if x_end is not None:
        constraints = constraints + [x[-1]  == x_end[:rdim]]
        constraints = constraints + [x[-2]  == zeros(rdim)]
        constraints = constraints + [x[-3]  == x[-2]]
    #~ f = -cp.sum(cp.log(cp.quad_form(tis[0]**2, identity(1)) - Tis[0]))
    #~ f = -cp.log(cp.quad_form(tis[0], identity(1)))
    #~ f = -cp.log(1-cp.quad_form(tis[0], identity(1)))
    #~ f = cp.max(0, cp.quad_form(tis[0], identity(1)) - 0.5)
    #~ f = cp.max(0, tis[0]- 0.5)
    #~ f = [Tis[k] - cp.quad_form(tis[k], identity(1)) for k in range(len(tis))]
    #~ f = -cp.log(tis[0])
    #~ f = cp.quad_form(tis[0], identity(1))
    #~ log_cost = [-cp.log(cp.quad_form(tis[k],identity(1))) for k in range(len(tis)) ]
    #~ log_cost = [-cp.quad_form(tis[k],identity(1)) - Tis[k] for k in range(len(tis)) ]
    
    #~ obj = cp.Minimize(sum(tis) + sum(log_cost))
    #~ obj = cp.Minimize(sum(tis) - sum(f))
    obj = cp.Minimize(sum(tis) + sum(Tis))
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS",verbose=False)
    #~ prob.solve(verbose=True)
    
    tis, Tis, xs = tovals(tis), tovals(Tis), tovals(x)
    
    tis = [max(tis[i], np.sqrt(Tis[i])) for i in range(len(tis))]
    print "tis", tis
    
    for j in range(1,len(tis)):
        idx = j*nvars
        if tis[j-1] < tis[j]:
            xs[idx-2] =  xs[idx+1] / tis[j] * tis[j-1]
            xs[idx-3] =  xs[idx+2] / tis[j] * tis[j-1]
        else:
            xs[idx+1] =  xs[idx-2] / tis[j-1] * tis[j]
            xs[idx+2] =  xs[idx-3] / tis[j-1] * tis[j]
    
    return prob, xs, tis
    
    

def bezierFromVal(xis, ti):
    wps = zeros((len(xis),dim))
    wps[0,:rdim] = array(xis[0])
    wps[1,:rdim] = array(xis[0]) +   array(xis[1])
    wps[2,:rdim] = array(xis[0]) + 2*array(xis[1]) + array(xis[2])
    wps[5,:rdim] = array(xis[5])
    wps[4,:rdim] = array(xis[5]) -   array(xis[4])
    wps[3,:rdim] = array(xis[5]) - 2*array(xis[4]) + array(xis[3])
    
    return bezier(wps.transpose(), ti)
        
def min_time(xs,num_phases,V):
    v = V[1]; V   = V[0][:,:rdim]; 
    A   = V        ; a = v
    constraints = []
    tis = [cp.Variable(1) for i in range(num_phases)]
    for j in range(num_phases):
        idx = j*nvars
        b = bezierFromVal(xs[idx:idx+nvars], 1)
        wps  = b.compute_derivate(1).waypoints()
        wpsa = b.compute_derivate(2).waypoints()
        t = tis[j]
        for i in range(wps.shape[1]):
            constraints = constraints + [ V.dot(wps[:2,i]) <= v *t]
        #~ for i in range(wpsa.shape[1]):
            #~ constraints = constraints + [ A.dot(wpsa[:2,i]) <= a *t*t ]
    obj = cp.Minimize(sum(tis))
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS",verbose=False)
    return tovals(tis)
    
    

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
            plotBezier(b, colors[-1], label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[-1],linewidth=2)
        
        from no_time import no_time, tailored_cost
        from no_time import bezierFromVal as beznotime
        
        ntis = min_time(xis,nphase,boundIneq())
        
        prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = None, tis = ntis)
        
        #~ print "xis", xis
        
        for i in range(nphase):
            b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)
            #~ plotBezier(b, colors[i], label = None, linewidth = 3.0)            
            plotBezier(b, colors[0], label = "classic", linewidth = 3.0)            
            
        prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = 3, tis = ntis)
        
        for i in range(nphase):
            b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)      
            plotBezier(b, colors[1], label = "jerk", linewidth = 3.0)   
            
        prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = 2, tis = ntis)
        
        
        for i in range(nphase):
            b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)       
            plotBezier(b, colors[2], label = "acc", linewidth = 3.0)         
            
        prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = 1, tis = ntis)
        
        
        for i in range(nphase):
            b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)
            plotBezier(b, colors[3], label = "vel", linewidth = 3.0)      
             
        #~ prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = 0, tis = ntis)
        
        
        #~ for i in range(nphase):
            #~ b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)        
            #~ plotBezier(b, colors[4], label = "distance", linewidth = 3.0)                
        plt.legend()
        plt.show()
        
        #~ plt.show()
        return b, tis
        
    b, tis = one(2)
    
