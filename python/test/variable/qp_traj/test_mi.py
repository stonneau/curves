
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array, identity
import gurobipy
eps =0.000001


degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2
totalnumvars = rdim * nvars

def integrate(init_val, all_vals):
    ints = [init_val]
    for _, j in enumerate(all_vals):
        ints = ints + [ints[-1] + j ]
    return ints[1:]

def boundIneq():
    v = ones(4)*1.; #v[:-3] = -v[3:]
    V = zeros((4,2)); 
    V[:2,:] = identity(2); 
    V[2:,:] =-identity(2); 
    return (V, v)
    
def mip1(Pis,num_steps,x_start, x_end, total_time):
    (V,v) = boundIneq()
    
    time_per_step = total_time / (float(num_steps)) * 1.2
    print "total time", total_time
    print "total time", time_per_step
    
    #~ P = P[:,:rdim]
    jerks = [cp.Variable(rdim) for i in range(num_steps)]
    accs = integrate(zeros(rdim), jerks)
    vels = integrate(zeros(rdim), accs)
    pos  = integrate(x_start[:rdim], vels)
    select = [[cp.Variable(1, boolean=True) for _ in range(len(Pis))] for _ in range(num_steps)]
    #~ tswitch = [cp.Variable(1, integer=True) for _ in range(len(Pis))] 
    #~ select = [cp.Variable(1, integer=True) for _ in range(num_steps)] 
    #~ constraints = [P * x <= p for x in pos]
    M = 5000
    eps = ones(rdim)*0.001
    #~ constraints = [pos[-1] == x_end[:rdim],vels[0] <= eps, accs[0] <= eps] 
    constraints = [pos[-1] == x_end[:rdim],  vels[-1]<= eps,  accs[-1] <= eps] 
    #~ constraints = [pos[-1] == x_end[:rdim],vels[-1] <= ones(rdim)*0.001,accs[-1] <= ones(rdim)*0.001] 
    #~ constraints = [pos[-1] == x_end[:rdim]] 
    #~ constraints = constraints + [ts >=0 for ts in tswitch] + [ts<= num_steps for ts in tswitch]
    #~ constraints = constraints + [ts >=1 for ts in select] + [ts<= num_steps for ts in select]
    #~ constraints = constraints + [ts[i+1] >= , ts<= num_steps for ts in tswitch]
    constraints = constraints + [sum(select[j]) == 1 for j in range(num_steps)]
    for j in range(num_steps-1):
        constraints = constraints + [Pis[i][0][:,:rdim] * x -Pis[i][1] <= M*(1-select[j][i]) for x in pos for i in range(len(Pis))]
        #~ constraints = constraints + [Pis[i][0][:,:rdim] * x -Pis[i][1] <= M*((tswitch[i] <= j) and (tswitch[i+1] >= j)) for x in pos for i in range(len(Pis))]
        #~ constraints = constraints + [Pis[i][0][:,:rdim] * x -Pis[i][1] <= M*(cp.quad_form(num_steps - select[j],identity(1))) for x in pos for i in range(len(Pis))]
    #~ constraints = constraints + [Pis[i][0][:,:rdim] * x -Pis[i][1] <1.2 M*((tswitch[i] <= j) * (tswitch[i+1] >= j)) for x in pos for i in range(len(Pis))]
    #~ constraints = constraints + [V * x <= v*20 for x in jerks]
    constraints = constraints + [V * x <= v * time_per_step for x in accs]
    constraints = constraints + [V * x <= v * time_per_step for x in vels]
    #~ obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for i,j in enumerate(jerks) if (i%2 == 0 or len(jerks)<20)]))
    obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in jerks]))
    #~ obj = cp.Minimize(sum([cp.quad_form(j,identity(rdim)) for j in vels]))
    prob = cp.Problem(obj, constraints)
    res = prob.solve(solver=cp.GUROBI, verbose=True )
    print "res ", res
    return pos

if __name__ == '__main__':
    from  cord_methods import *
    from  plot_cord import *
    
    def one(nphase=2, num_steps = 10):
        plt.close()
        pDef, Pis = genProblemDef(nvars,nphase)
        #~ pDef, inequalities_per_phase = genProblemDef(nphase,nphase)
        V = boundIneq()
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))
        
        
        ##### COMPARISON TRAJECTORY #################"
        #find time by solving sqp
        from sqp import solve_straight_lines, min_time
        prob, xis, tis = solve_straight_lines(Pis[:], V, x_start=x_start, x_end=x_end)
        from no_time import no_time, tailored_cost
        from no_time import bezierFromVal as beznotime
        
        ntis = min_time(xis,nphase,boundIneq())
        prob, xis = no_time(Pis[:], V, x_start=x_start, x_end=x_end, ccost = 3, tis = ntis)
        
        for i in range(nphase):
            b = beznotime(xis[i*nvars:i*nvars+nvars], 1.)
            plotBezier(b, colors[3], label = "vel", linewidth = 3.0)      
        
        ##### MIXED INTEGER PROGRAMM #######"
        
        res = mip1(Pis, num_steps, x_start, x_end, total_time = sum(ntis))
        xs = [e.value for e in res]
        xp = [x_start[0]] + [e.value[0] for e in res]
        yp = [x_start[1]] + [e.value[1] for e in res]
        plt.scatter(xp,yp,color = "k")
                
        #~ print "times", tis
        #~ for i in range(nphase):
            #~ ti = abs(tis[i][0])
            #~ b = bezierFromVal(xis[i*4:i*4+4], 1., 4)
            #~ print "init vel", b.compute_derivate(1)(0.)
            #~ print "end vel", b.compute_derivate(1)(b.max())
            #~ plotBezier(b, colors[i], label = None, linewidth = 3.0)
            #~ plotControlPoints(b, colors[i],linewidth=2)
        
        plt.show()
        #~ coeffs = np.polyfit(xp,yp,7)
        #~ ffit = np.poly1d(coeffs)
        #~ x_new = np.linspace(10, 40, num=len(xp)*10)
        #~ plt.plot(x_new, ffit(x_new)) 
        return xp, yp
        
    xp, yp = one(3,20)
    #~ x = cp.Variable(2, boolean=True)
    #~ y = cp.Variable(1)
    #~ eq2 = y >= -2
    #~ eq1 = y >= -1
    
    #~ M = 5000
    #~ m = -5000
    
    #~ constraints = [
    #~ -y - 2 <= M*(1-x[0]),
    #~ -y - 1 <= M*(1-x[1]),    
    #~ x[0] + x[1] == 1
    #~ ];
    #~ obj = cp.Minimize(y)
    #~ prob = cp.Problem(obj, constraints)
    #~ res = prob.solve(verbose=True)
    #~ print res
    
