import cvxpy as cp
from numpy.random import randn
from numpy import eye

from straight_lines import *


#given 4 control points of 2 connecting bezier curves, find minimum time to achieve
#transition while achieving position and velocity continuity
# only variables are times and end velocities, endposition is fixed
def one_step(Pi, Pip1, xi, xi1, ti,ti1, V): #TODO: add V and A  
    # we have only 2 variables, position and before last control point that will set end
    #velocity  
    x  = [cp.Variable(rdim) for i in range(2)] 
    t  = cp.Variable(1)
    t1 = cp.Variable(1)
    
    P    = Pi  [0][:,:rdim]       ; p  = Pi  [1]
    P1   = Pip1[0][:,:rdim]       ; p1 = Pip1[1]
    Vn   = V[0][:,:rdim]*n        ; v = V[1]
    
    #offset with respect to control point index
    off = 2
    
    
    dc0_0_n_t = xi[1]
    dc1_0_n_t = x[3-off]-x[2-off]-xi[1]-xi[0]
    dc2_0_n_t = x[2-off]
    
    ddc0_0_t_t = 1./(n-1) *(dc1_0_n_t  - dc0_0_n_t)
    ddc2_0_t_t = 1./(n-1) *(dc2_0_n_t  - dc1_0_n_t)
    
    constraints = [
    P*(x[3-off]) <= p,
    P*(x[3-off] - x[2-off]) <= p,
    #velocities
    #~ Vn*x[2-off] <= v * t,
    #~ Vn*(x[3-off]-x[2-off]-xi[1]-xi[0]) <= v * t,
    Vn*dc1_0_n_t <= v * t,
    Vn*dc2_0_n_t <= v * t,
    #acceleration
    Vn*ddc0_0_t_t <= v * t*ti,#this is an approximation because we premultiply by the previous ti and not the real value
    Vn*ddc2_0_t_t <= v * t*ti,#this is an approximation because we premultiply by the previous ti and not the real value
    ]
    
    constraints = constraints +[
    #position constraints for next phase
    P1*(x[3-off]) <= p1, #start pos in constraint
    P1*(x[3-off] + x[2-off]) <= p1, #vel pos in constraint
    Vn*x[2-off] <= v * t1,
    Vn*(xi1[3]-xi1[2]-x[2-off]-x[3-off]) <= v * t1,
    ]    
    
    obj = cp.Minimize(t + t1)
    prob = cp.Problem(obj, constraints)
    prob.solve(solver="ECOS", verbose=False, )
    #~ prob.solve(verbose=False, )
    
    #speed is given by minimum velocity between two solutions
    cn = x[2-off].value / t.value[0]
    if t1.value[0] > t.value[0]:
        cn = x[2-off].value / t1.value[0]
    
    
    #reassigning values
    xi [3] = x[3-off].value
    xi1[0] = xi[3]
    xi [2] = cn * t.value[0]
    # as long as ti1 / t.value[0] > 1, the problem will be feasible
    #~ xi1[1] = x[2-off].value * (ti1 / t.value[0])
    xi1[1] = cn * t1.value[0]
    
    return prob, xi, xi1, t.value[0],t1.value[0]
    
    


from hpp_spline import bezier

if __name__ == '__main__':
    from  cord_methods import *
    from  plot_cord import *
    
    
        
    def one(nphase=6):
        plt.close()
        pDef, inequalities_per_phase = genProblemDef(nvars,nphase)
        #~ pDef, inequalities_per_phase = genProblemDef(nphase,nphase)
        V = boundIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))

        prob, xis, tis = solve_straight_lines(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end)
        for i in range(len(tis)):
            ti = abs(tis[i])[0]
            b = bezierFromVal(xis[i*4:i*4+4], abs(ti))
            plotBezier(b, colors[i], label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[i],linewidth=2)
        
        tis = [el[0] for el in tis]
        
        print("status:", prob.status)
        
        def refinePass(start = 0, draw = False, add = 0):
            #now try to improve
            for i in range(start,nphase-start,2):
                prob, xi, xi1, t0,t1 = one_step(inequalities_per_phase[i], inequalities_per_phase[i+1], xis[i*4:i*4+4], xis[(i+1)*4:(i+1)*4+4], tis[i],tis[i+1], V)
                xis[i*4:i*4+4] = xi[:]
                xis[(i+1)*4:(i+1)*4+4] = xi1[:]
                tis[i] = t0
                tis[i+1] = t1
                
                if draw:
                    b = bezierFromVal(xi, abs(t0))
                    plotBezier(b, colors[1+add], label = None, linewidth = 3.0)
                    #~ plotControlPoints(b, colors[start],linewidth=2)
                    
                    b = bezierFromVal(xi1, abs(t1))
                    plotBezier(b, colors[1+add], label = None, linewidth = 3.0)
                    #~ plotControlPoints(b, colors[start],linewidth=2)
                
                #~ plt.plot(x_start[0],x_start[1],"g",linewidth=12)
                #~ plt.plot(x_end[0],x_end[1],"g",linewidth=12)
                
                
        for _ in range(2):
            refinePass(0,);
            refinePass(1,);
        refinePass(0,True,0);
        refinePass(1,True,0);
        #~ for _ in range(10):
            #~ refinePass(1);
            #~ refinePass(0);
        #~ refinePass(0,True,3);
        #~ refinePass(1,True,3);
        
        
        pDef.costFlag = derivative_flag.VELOCITY
        #~ print "tis ", tis
        tis = addSplit(tis)
        lastval = tis[-1]
        tis = [el / lastval for el in tis[:-1]]  
        #~ print "tis ", tis
        #~ pDef.splits = array([tis]).T 
        #~ pDef.splits = array([genSplit(nphase)]).T 
        
        from test_cord import gen as gen_py
        gen_py(save = False, Pis = inequalities_per_phase, start=x_start, end = x_end, tis = tis)
        #~ res, cos = computeTrajectory(pDef, saveToFile=False)
        #~ pDef.degree = 8
        #~ plot(pDef,res, "", saveToFile=False)
        
        plt.show()
        #now retrieve the tis and tries to solve original problem
        
        
    one(6)
    
