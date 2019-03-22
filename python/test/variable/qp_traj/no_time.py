
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array, identity
eps =0.000001


degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2


from scipy.special import binom as bino

binos = {}
for deg in range(2,degree+1):
    n_pts = deg+1
    binos[deg] = [bino(deg,i) for i in range(n_pts)]

def derivate(xvars_phase, t, num):
    if num <= 0:
        return xvars_phase
    deg = len(xvars_phase)-1
    print "deg", deg
    deriv = [ (deg / t) * (xvars_phase[j+1]-xvars_phase[j]) for j in range(deg) ]
    return derivate(deriv,t,num-1)
       
def eval_at_point_square(pis, t, square = True):
    if (t <0):
        print "t", t
    if (t >1.):
        print "t", t
    assert(t>=0.)
    assert(t<=1.)
    deg = len(pis)-1
    n_pts = deg+1
    res = sum([binos[deg][i] * t**i * (1-t)**(deg-i)*pis[i] for i in range(n_pts)])
    if square:
        return cp.quad_form(res,identity(2))
    else:
        return res
    
def eval_at_points(pis, T, timestep):
    res = 0
    for i in range(int(T/timestep)):
        res = res + eval_at_point_square(pis, (float)(i)*timestep / T)
    res = res + eval_at_point_square(pis, 1.)
    return res
    
def eval_at_pointszero(pis, T, timestep):
    res = 0
    lastval = 0.
    for i in range(1,int(T/timestep)):        
        lastval = (float)(i)*timestep / T
        rest = eval_at_point_square(pis, lastval, square = False) - eval_at_point_square(pis, (float)(i-1)*timestep / T, square = False)
        res += cp.quad_form(rest,identity(2))
    res = res+ cp.quad_form(eval_at_point_square(pis, 1., square = False) - eval_at_point_square(pis, lastval, square = False),identity(2))
    return res
    
def eval_cost_integral(xvars, tis, num_points, num_deriv):
    T = sum(tis)
    timestep = T / num_points
    res = 0    
    if num_deriv ==0:
        for i in range(len(xvars) / nvars - 1):
            if tis[i] > 0.:
                pis = derivate(xvars[nvars*i:nvars*(i+1)],tis[i], num_deriv)
                res += eval_at_pointszero(pis, tis[i], timestep)
    else:
        for i in range(len(xvars) / nvars):
            if tis[i] > 0.:
                pis = derivate(xvars[nvars*i:nvars*(i+1)],tis[i], num_deriv)
                res += eval_at_points(pis, tis[i], timestep)
        
    return res
#~ def cost_integral(xvars, t, num, num_eval):
    
    #~ derivars = derivate(xvars, t, num)

#~ def shortbezier(xvars, t):

#find the minimum time sequence of straight lines connecting from start to end position along union of convex sets

def no_time(Pis, V, x_start = None, x_end = None, ccost = None, tis = None): #TODO: add V and A    
    num_phases = len(Pis)
    num_vars = nvars * num_phases    
    if tis is None:
        tis = ones(len(Pis))
    tis = [el * 1.0 for el in tis]
    x   = [cp.Variable(rdim) for i in range(num_vars)] 
    Vn   = V[0][:,:rdim]*n        ; v = V[1]
    Ann1 = Vn*(n-1)        ; a = v * 10.
    
    constraints = []
    
    for j in range(len(Pis)):
        idx = j*nvars
        Pi = Pis[j]
        P    = Pi[0][:,:rdim]       ; p  = Pi  [1]
        #positions
        for k in range(nvars):            
            constraints = constraints + [
            P*(x[idx+k]) <= p 
        ]
        #velocities
        vel = [x[idx+k]-x[idx+k-1] for k in range(1,nvars)]
        constraints = constraints + [Vn*(el) <= v*tis[j] for el in vel]
        #accelerations
        acc = [vel[k+1]-vel[k] for k in range(0,nvars-2)]
        constraints = constraints + [Ann1*(el) <= a*tis[j]*tis[j] for el in acc]
        
    #continuity constraints
    for j in range(1,len(Pis)):
        idx = j*nvars
        constraints = constraints + [
        x[idx]   ==  x[idx-1], #just positions because velocities at 0
        (x[idx+1] -   x[idx]) * tis[j-1]    ==   (x[idx] -   x[idx-2])  * tis[j] ,#  v/n = (p1 - p0) / to = (p0-p_n_1) / t1
        (x[idx+3] - 2*x[idx+1] + x[idx]) * tis[j-1] * tis[j-1]   == (x[idx-3] - 2*x[idx-2] + x[idx]) * tis[j] * tis[j] #just positions because velocities at 0
        ]
            
    if x_start is not None:
        constraints = constraints + [x[0]  == x_start[:rdim]]
        constraints = constraints + [x[0]  == x[1]]
        constraints = constraints + [x[1]  == x[2]]
    if x_end is not None:
        constraints = constraints + [x[-1]  == x_end[:rdim]]
        constraints = constraints + [x[-1]  == x[-2]]
        constraints = constraints + [x[-2]  == x[-3]]
    
    #~ cost = [cp.quad_form(x[idx],identity(x[idx].shape[0])) for idx in range(1,2)]
    #~ cost = sum([cp.quad_form(x[idx],identity(x[idx].shape[0])) - cp.quad_form(x[idx+1],identity(x[idx].shape[0])) for idx in range(len(x)-1)])
    
    #~ print "tis", tis
    #~ tis = [idx%nvars for idx in range(len(x)-1)]
    #~ print "tis", tis
    tis = [tis[idx/nvars] for idx in range(len(x)-1)]
    #~ print "tis", tis
    vel = []; acc = []; jerk = []
    for i in range(num_phases):
        off  = i * nvars
        vb = [(x[off+idx+1] - x[off+idx]) / degree    * tis[i] for idx in range(nvars-1)]
        ab = [(vb[idx+1] - vb[idx])    /(degree-1) * tis[i] for idx in range(len(vb)-1)]
        jerk = jerk + [(ab[idx+1] - ab[idx])  /(degree-2) * tis[i] for idx in range(len(ab)-1)]
        vel = vel + vb
        acc = acc + ab
    #~ v = [(x[idx+1] - x[idx]) / degree * tis[idx/nvars] for idx in range(len(x)-1)]
    #~ a = [(v[idx+1] - v[idx]) / (degree-1) * tis[idx/(nvars-1)] for idx in range(len(v)-1)]
    #~ obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in vel]))    
    #~ obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in acc]))    
    obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in jerk]))    
    if ccost is not None:
        print "ccost", ccost
        if ccost == 1:
            obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in vel]))    
        elif ccost == 2:
            obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in acc]))    
        elif ccost == 3:
            obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in jerk]))    
        else:
            obj = cp.Minimize(eval_cost_integral(x, tis, num_points = 40, num_deriv = ccost))
        #~ obj = cp.Minimize(sum ([cp.quad_form(xi, ccost[k%nvars][0]) + ccost[k%nvars][1].T * xi for k, xi in enumerate(x)]))
        
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS",verbose=False)
    return prob, tovals(x)
    
#~ def makeVelContinuous(xis):
    #~ for k in range(nvars):
        #compute end vel
    

def bezierFromVal(xis, ti, numvars = nvars):
    wps = zeros((len(xis),dim))
    for i in range(numvars):
        wps[i,:rdim] = array(xis[i])
    #~ wps[1,:rdim] = wps[0,:rdim] + array(xis[1])
    #~ wps[3,:rdim] = array(xis[3])
    #~ wps[2,:rdim] = wps[3,:rdim] - array(xis[2])
    
    return bezier(wps.transpose(), ti)
    
def bezierFromValAcc(xis, ti):
    wps = zeros((len(xis),dim))
    for i in range(3):
        wps[i,:rdim] = array(xis[0])
    for i in range(3,6):
        wps[i,:rdim] = array(xis[-1])    
    return bezier(wps.transpose(), ti)
    
def min_time(xs,V):
    b = bezierFromValAcc(xs, 1)
    wps = b.compute_derivate(1).waypoints()
    t = cp.Variable(1) 
    constraints = [ t >= 0.00001]
    for i in range(wps.shape[1]):
        #~ if (np.linalg.norm(wps[:2,i]) > 0.01):
        constraints = constraints + [ V[0].dot(wps[:2,i]) <= V[1] *t ]
    obj = cp.Minimize(t)
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS",verbose=True)
    return t.value[0]
    
def tailored_cost(c, nvars):
    A = c.A; b = c.b
    return [(A[i:i+2,i:i+2],b[i:i+2]) for i in range(0,nvars*3,3)]
    
        #~ pData = generate_problem(pDef)
        #~ costPerPhase = pData.cost
        

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
        pDef, Pis = genProblemDef(nvars,nphase)
        V = boundIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))
        pDef.splits = problemDefinition().splits
        pDef.flag = NONE
        pDef.costFlag = JERK
        pData = generate_problem(pDef)
        costPerPhase = pData.cost

        #compute initial times
        from straight_lines import solve_straight_lines
        from straight_lines import bezierFromVal as bezierStraight
        prob, xis, tis = solve_straight_lines(Pis, V, x_start = x_start, x_end = x_end)
        
        #~ print "xis", len(xis)
        #~ print "tis", len(tis), tis
        
        #~ for i in range(len(tis)):
            #~ ti = abs(tis[i][0])
            #~ b = bezierStraight(xis[i*4:i*4+4], abs(ti))
            #~ plotBezier(b, "r", label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[i],linewidth=2)

        #~ prob, xis = no_time(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end, ccost = tailored_cost(costPerPhase, nvars), tis = tis)
        
        #~ for i in range(nphase):
            #~ b = bezierFromVal(xis[i*nvars:i*nvars+nvars], 1.)
            #~ plotBezier(b, colors[i], label = None, linewidth = 3.0)
        

        prob, xis = no_time(Pis[:], V, x_start=x_start, x_end=x_end, ccost = None)
        
        for i in range(nphase):
            b = bezierFromVal(xis[i*nvars:i*nvars+nvars], 1.)
            plotBezier(b, "g", label = None, linewidth = 3.0)
            
        prob, xis = no_time(Pis[:], V, x_start=x_start, x_end=x_end, ccost = 2, tis = tis)
        
        for i in range(nphase):
            b = bezierFromVal(xis[i*nvars:i*nvars+nvars], 1.)
            plotBezier(b, "b", label = None, linewidth = 3.0)
        
        print "tis, ", tis
        
        plt.show()
        return b
        
    b = one(2)
    
