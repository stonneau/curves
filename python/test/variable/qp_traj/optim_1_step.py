import cvxpy as cp
from numpy.random import randn
from numpy import eye

from straight_lines import *


#given 4 control points of 2 connecting bezier curves, find minimum time to achieve
#transition while achieving position and velocity continuity
# only variables are times and end velocities, endposition is fixed
def one_step(Pi, Pip1, xi, xi1, ti,ti1, V): #TODO: add V and A  
    print "t0", ti
    print "t1", ti1
    if(ti[0] < ti1[0]):
        raise ValueError
    # we have only 2 variables, position and before last control point that will set end
    #velocity  
    x  = [cp.Variable(rdim) for i in range(2)] 
    t = cp.Variable(1)
    
    P    = Pi  [0][:,:rdim]       ; p  = Pi  [1]
    P1   = Pip1[0][:,:rdim]       ; p1 = Pip1[1]
    Vn   = V[0][:,:rdim]*n        ; v = V[1]
    
    #offset with respect to control point index
    off = 2
    
    constraints = [
    t >= ti1,
    P*(x[3-off]) <= p,
    P*(x[3-off] - x[2-off]) <= p,
    #end point constraint    
    P1*x[3-off] <= p1,
    #velocities
    Vn*x[2-off] <= v * t,
    Vn*(x[3-off]-x[2-off]-xi[1]-xi[0]) <= v * t,
    ]
    
    constraints = constraints +[
    #position constraints for next phase
    P1*(x[3-off]) <= p1, #start pos in constraint
    P1*(x[3-off] + x[2-off]) <= p1, #vel pos in constraint
    Vn*(xi1[3]-xi1[2]-x[1]-x[0]) <= v * ti,
    ]    
    
    obj = cp.Minimize(t)
    prob = cp.Problem(obj, constraints)
    prob.solve()
    
    #reassigning values
    xi [3] = x[3-off].value
    xi1[0] = xi[3]
    xi [2] = x[2-off].value
    xi1[1] = x[2-off].value / t.value[0] * ti1
    #~ xi1[1] = x[2-off].value 
    print "new t0", t.value[0]
    
    return prob, xi, xi1, t.value[0],ti1[0]
    
    


from hpp_spline import bezier

if __name__ == '__main__':
    from  cord_methods import *
    from  plot_cord import *
    
    def oneVelIneq():
        v = ones(4)*0.001; #v[:-3] = -v[3:]
        V = zeros((4,2)); 
        V[:2,:] = identity(2); 
        V[:2,:] =-identity(2); 
        return (V, v)
        
    def one(nphase=3):
        plt.close()
        pDef, inequalities_per_phase = genProblemDef(nvars,nphase)
        V = oneVelIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.start.reshape((-1,))

        prob, xis, tis = solve_straight_lines(inequalities_per_phase[:], V, x_start=x_start, x_end=x_end)
        for i in range(len(tis)):
            ti = abs(tis[i])[0]
            #~ print("optimal value", ti)
            b = bezierFromVal(xis[i*4:i*4+4], abs(ti))
            plotBezier(b, colors[i], label = None, linewidth = 3.0)
            plotControlPoints(b, colors[i],linewidth=2)
        
        print("status:", prob.status)
        
        #now try to improve
        prob, xi, xi1, t0,t1 = one_step(inequalities_per_phase[0], inequalities_per_phase[1], xis[0:4], xis[4:8], tis[0],tis[1], V)
        
        b = bezierFromVal(xi, abs(t0))
        plotBezier(b, colors[2], label = None, linewidth = 3.0)
        plotControlPoints(b, colors[2],linewidth=2)
        
        b = bezierFromVal(xi1, abs(t1))
        plotBezier(b, colors[3], label = None, linewidth = 3.0)
        plotControlPoints(b, colors[3],linewidth=2)
        
        plt.plot(x_start[0],x_start[1],"g",linewidth=12)
        plt.plot(x_end[0],x_end[1],"g",linewidth=12)
        plt.show()
        return b
        
    b = one()
    
