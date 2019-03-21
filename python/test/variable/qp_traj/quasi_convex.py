import cvxpy as cp
from numpy.random import randn
from numpy import eye

from straight_lines import *

def wpsFromStraightLine(xis, nphase):
    res = []
    for i in range(nphase):
        xi = xis[i*4:i*4+4]
        res = res + [xis[0], xis[0] + xis[1], xis[3]-xis[2], xis[3]]
    return res
    
def bezierFromRefined(xis,ti):
    wps = zeros((len(xis),dim))
    for i in range(4):
        wps[i,:rdim] = array(xis[i])
    return bezier(wps.transpose(), ti)
    

#given 4 control points of 2 connecting bezier curves, find minimum time to achieve
#transition while achieving position and velocity continuity
# only variables are times and end velocities, endposition is fixed
def one_step(Pi, Pip1, ci, ci1, t,t1, Vv): #TODO: add V and A  
    # we have only 2 variables, position and before last control point that will set end
    # here the control points are the variables
    x  = [cp.Variable(rdim) for i in range(3)] 
    #~ t  = cp.Variable(1)
    #~ t1 = cp.Variable(1)
    
    P    = Pi  [0][:,:rdim]       ; p  = Pi  [1]
    P1   = Pip1[0][:,:rdim]       ; p1 = Pip1[1]
    V    = Vv[0]         ; v = Vv[1]
    
    #offset with respect to control point index
    
    c2_0 = cp.Variable(rdim)
    c3_0 = cp.Variable(rdim);
    c1_1 = cp.Variable(rdim)
    
    #dc0 constant
    dc0_0 = t /n * (ci[1]  - ci[0] )
    dc1_0 = t /n * (c2_0   - ci[1] )
    dc2_0 = t /n * (c3_0   - c2_0  )
    dc0_1 = t1/n * (c1_1   - c3_0  )
    dc1_1 = t1/n * (ci1[2] - c1_1  )
    dc2_1 = t1/n * (ci1[3] - ci1[2])
    
    ddc0_0 = t /(n-1) *(dc1_0  - dc0_0)
    ddc2_0 = t /(n-1) *(dc2_0  - dc1_0)
    ddc0_1 = t1/(n-1) *(dc1_1  - dc0_1 )
    ddc1_1 = t1/(n-1) *(dc2_1  - dc1_1)
    
    off = 2
    
    
    constraints = [
    P * (c2_0) <= p,
    P * (c3_0) <= p,
    P1* (c3_0) <= p1,
    P1* (c1_1) <= p1,
    #velocities
    V * (dc1_0)  <= v,
    V * (dc2_0)  <= v,
    V * (dc0_1)  <= v,
    V * (dc1_1)  <= v,
    ]
    
    #continuity constraints
    constraints = constraints +[
    dc2_0 == dc0_1
    ]
    
    #~ constraints = constraints +[
    #~ #position constraints for next phase
    #~ P1*(x[3-off]) <= p1, #start pos in constraint
    #~ P1*(x[3-off] + x[2-off]) <= p1, #vel pos in constraint
    #~ Vn*x[2-off] <= v * t1,
    #~ Vn*(xi1[3]-xi1[2]-x[2-off]-x[3-off]) <= v * t1,
    #~ ]    
    
    obj = cp.Minimize(dc1_0[0])
    prob = cp.Problem(obj, constraints)
    res = prob.solve(solver="ECOS", verbose=True, )
    print res
    #~ prob.solve(verbose=False, )
    
    
    c2_0 = cp.Variable(rdim)
    c3_0 = cp.Variable(rdim);
    c1_1 = cp.Variable(rdim)
       
    #reassigning values
    ci [2] = c2_0.value
    ci [3] = c3_0.value
    ci1[0] = c3_0.value
    ci1[1] = c1_1.value
    
    return prob, ci, ci1
    
    


from hpp_spline import bezier

if __name__ == '__main__':
    from  cord_methods import *
    from  plot_cord import *
    
    
        
    def one(nphase=6):
        plt.close()
        pDef, inequalities_per_phase = genProblemDef(nvars,nphase)
        V = boundIneq()
        print "v1", V[1]
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))

        prob, xis, tis = solve_straight_lines(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end)
        for i in range(len(tis)):
            ti = abs(tis[i])[0]
            b = bezierFromVal(xis[i*4:i*4+4], abs(ti))
            print "wtf0",  b.waypoints()
            print "wtf",  b.compute_derivate(1).waypoints()
        
        #velocity / acc variables replaced with control points variables
        xis = wpsFromStraightLine(xis,nphase)        
        tis = [el[0] for el in tis]
        
        
        #now find time for one phase
        
        
        def refinePass(start = 0, draw = False):
            #now try to improve
            for i in range(start,nphase-start,2):
                print "1", i
                print "1", xis[i*4:i*4+4]
                print "1", xis[(i+1)*4:(i+1)*4+4]
                prob, xi, xi1 = one_step(inequalities_per_phase[i], inequalities_per_phase[i+1], xis[i*4:i*4+4], xis[(i+1)*4:(i+1)*4+4], tis[i],tis[i+1], V)
                xis[i*4:i*4+4] = xi[:]
                xis[(i+1)*4:(i+1)*4+4] = xi1[:]
                
                if draw:
                    b = bezierFromRefined(xi, abs(t0))
                    plotBezier(b, colors[start], label = None, linewidth = 3.0)
                    #~ plotControlPoints(b, colors[start],linewidth=2)
                    
                    b = bezierFromRefined(xi1, abs(t1))
                    plotBezier(b, colors[start], label = None, linewidth = 3.0)
                    #~ plotControlPoints(b, colors[start],linewidth=2)
                
                #~ plt.plot(x_start[0],x_start[1],"g",linewidth=12)
                #~ plt.plot(x_end[0],x_end[1],"g",linewidth=12)
                
                
        for _ in range(10):
            refinePass(1);
            refinePass(0);
        refinePass(0,True);
        refinePass(1,True);
        
        #~ pDef.costFlag = derivative_flag.VELOCITY
        #~ print "tis ", tis
        #~ tis = addSplit(tis)
        #~ lastval = tis[-1]
        #~ tis = [el / lastval for el in tis[:-1]]  
        #~ print "tis ", tis
        #~ pDef.splits = array([tis]).T 
        #~ pDef.splits = array([genSplit(nphase)]).T 
        #~ res, cos = computeTrajectory(pDef, saveToFile=False)
        #~ plot(pDef,res, "", saveToFile=False)
        
        plt.show()
        #now retrieve the tis and tries to solve original problem
        
        
    one(6)
    
