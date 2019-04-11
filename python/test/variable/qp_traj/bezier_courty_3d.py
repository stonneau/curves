
from __future__ import division
import cvxpy as cp
from numpy.random import randn
from numpy import matrix, array, zeros, ones, diag, cross, vstack, identity
eps =0.000001

from numpy.linalg import pinv, norm

degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 3

#try to solve distance based ik for a simple manipulator with 2 DOF, each with length 1


##############" CONSTRAINTS #############
def genf2vars(dist,i,j):
    def f(x):
        xi = x[rdim*i:rdim*(i+1)]
        xj = x[rdim*j:rdim*(j+1)]
        return ( norm(xi-xj)**2- dist*dist)**2
    return f
    
def genf(dist,i):
    def f(x):
        xi = x[rdim*i:rdim*(i+1)]
        return (xi.T.dot(xi) - dist*dist)**2
    return f

def gen_target(x_end):
    def f(x):
        x1 = x[rdim:] - x_end
        return x1.dot(x1)
    return f

# A x <= b

def ineqf(A,b,x):
    blen = b.shape[0]
    lenx = (int)(x.shape[0] / rdim)
    res = zeros(blen * lenx)
    for i in range(lenx):
        val = A.dot(x[rdim*i:rdim*(i+1)]) - b
        res[i*blen:(i+1)] = np.where(val>-0.01, val, 0.)
    return res

def gen_ineq(A,b):
    def f(x):
        return ineqf(A,b,x)
        #~ res = ineqf(A,b,x)
        #~ return res.T.dot(res)
    return f

#distance from base to joint 1 
f10 = genf(1,0)

#distance from joint 2 to joint 1 
f12 = genf2vars(1,1,0)

def positions():
    def f(x):
        res = zeros(2)
        res[0] = f10(x)
        res[1] = f12(x)
        return res
    return f


##############" END CONSTRAINTS #############


##############" JACOBIANS #############


###### HARD CONSTRAINTS JACOBIANS ######

#position

def Jpositions():
    def f(x):
        J = zeros((2,2*rdim))
        x0 = x[:rdim]
        x1 = x[rdim:]
        x1mx0 = (x1-x0)
        x0mx1 = (x0-x1)    
        df1dx0 = 4*x0*(x0.T.dot(x0) - 1)
        df1dx1 = zeros(rdim)

        df2dx0 = 4*(x0mx1)*(x1mx0.T.dot(x1mx0) -1)
        df2dx1 = 4*(x1mx0)*(x1mx0.T.dot(x1mx0) -1)
        
        J[0,:rdim] = df1dx0
        J[0,rdim:] = df1dx1
        J[1,:rdim] = df2dx0
        J[1,rdim:] = df2dx1
        return J
    return f
    
    
def Jineq(A,b):
    
    blen = b.shape[0]
    Arows = A.shape[0]
    blen = b.shape[0]
    #~ fineq  = gen_ineq(A, b)
    def f(x):
        #~ ineqf(A,b,x)
        lenx = x.shape[0]
        lenxrdim = (int)(lenx / rdim)
        #~ fin = fineq(x)
        #~ fin = ineqf(A,b,x)
        J = zeros((lenxrdim*Arows,lenx))
        for i in range(lenxrdim):
            xi = x[rdim*i:rdim*(i+1)]
            fval = ineqf(A,b,x)
            if fval.T.dot(fval) > 0:
            #~ if residual > 0.:
                 J[i*Arows,i*rdim:(i+1)*rdim] = A[:]
            #~ val = A.dot(x[rdim*i:rdim*(i+1)]) - b
            #~ res[i*blen:(i+1)] = np.where(val>0., val, 0.)
                #~ if(fin[0]>0.): 
                    #~ J[0:Arows,:2] = A[:]
                #~ if(fin[1]>0.):
                    #~ J[Arows: ,2:] = A[:]
        return J
    return f
  
###### SOFT CONSTRAINTS JACOBIANS ######  

def Jtarget(x_end):  
    def fu(x):
        J = zeros((1,2*rdim)) 
        J[0,rdim:] = 2*(x[rdim:]-x_end) 
        return J
    return fu
    
##############" CONSTRAINT DIC #############
    
all_constraints = { "pos" : [positions, Jpositions], "target" : [gen_target, Jtarget], "ineq" : [gen_ineq, Jineq]}

def constraint(name,*args):
    cons = all_constraints[name]
    return [cons[0](*args), cons[1](*args)]
    
##############" CONSTRAINT DIC #############


    
##############" BEZIER STUFF  #############

def bezier2(x0,x2,rank):
    def ft(t):
        def f(x): 
            xi = x[rank:rank+rdim]
            return 6. * t * x2 * (1-t)**2 + 4.*t* xi * (1-t)**3  + x0 * (1-t)**4
    return ft
            
def gradBezier2(rank):
        def ft(t):
            def f(x):
                grad = 4.*t* (1-t)**3
                J = zeros((1,2*rdim)) 
                J[0,rank:rank+rdim] = array([grad for _ in rdim])
                return J
        return ft
        
# discretize on 4 points
def genBezier(x0,x2, rank = 0, nsteps=20):
    res = []
    fns = float(nsteps)
    
    bts = bezier2(x0,x2,rank)
    gbs = gradBezier2(rank)
    
    for i in range(nsteps+1):
        t = float(i) / fns
        res += [[bts(t),gbs(t)]]
    return res
    
    
from sympy import *
from sympy.tensor.array import derive_by_array
init_printing(use_unicode=True)
from scipy.special import binom as bino

t = symbols('t')

#for the moment start / goal are given
def symbolicbezier(numvars, pointId):
    symbs = [symbols('x'+str(pointId)+'_'+str(i))  for i in range(numvars)] 
    constants =  [symbols('x'+str(pointId)+'s')] + [symbols('x'+str(pointId)+'g')]    
    npts = len(symbs) + len(constants)
    deg = npts -1
    pis = [constants[0]] + symbs + [constants[1]]
    factors0 = [bino(deg,i) * t**i * (1-t)**(deg-i)*pis[i] for i in range(len(pis))]
    eq0 = sum(factors0);
    return eq0, pis

def makeTemporalBezierVal(symBezier, pis, xsId, xgId):    
    def ft(ti):
        fsymbt = lambdify([t],symBezier)     
        evalt = fsymbt(ti) #eval delta t and exponents 
        fsymb =  lambdify(pis,evalt)  
        def fu(x, fsymb = fsymb):
            #~ print "x", x
            xvars = tuple([xsId]+[x[i*rdim:(i+1)*(rdim)] for i in range(len(pis[1:-1]))] +[xgId])       
            #~ print "xvaras ",    xvars  
            return fsymb(*xvars)
        return fu, evalt
    return ft
        
def generateBezierSymbolic(numvars, pointId, xsId, xgId):
    symb, pis = symbolicbezier(numvars, pointId)
    return makeTemporalBezierVal(symb, pis, xsId, xgId), pis

def bezier2Dist(xstart, xend, numControlPoints = 1, nsteps=10):
    xs0 = xstart[:rdim]
    xs1 = xstart[rdim:]
    xg0 = xend[:rdim]
    xg1 = xend[rdim:]
    fns = float(nsteps)
    res = []
    
    ft0, pis0 = generateBezierSymbolic(numControlPoints, 0, xs0, xg0)
    ft1, pis1 = generateBezierSymbolic(numControlPoints, 1, xs1, xg1)
    
    for i in range(nsteps+1):
        t = float(i) / fns
        pt0t, eq0 = ft0(t)
        pt1t, eq1 = ft1(t)
        
        args = pis0[1:-1] + pis1[1:-1]
        
        #distance to base must be 1
        distBase = (eq0*eq0 - 1)**2
        #~ gradDistbase = derive_by_array(distBase, (pis0[0], pis1[1]))
        gradDistbase = derive_by_array(distBase, args)
        lambdaGradDistbase = lambdify(pis0,gradDistbase)  
        
        #~ print "gradDistbase", gradDistbase
        
        #distance between two variables must be 1
        dist10 = ((eq1-eq0)*(eq1-eq0) - 1)**2
        gradDist10 = derive_by_array(dist10, args)
        lambdaGradDist10 = lambdify(pis0+pis1,gradDist10)  
        
        def fu1(x, pt0t=pt0t, pt1t=pt1t):
            x0 = x[:rdim*numControlPoints]
            x1 = x[rdim*numControlPoints:]
            pt0 = pt0t(x0)
            pt1 = pt1t(x1)
            return ( norm(pt1-pt0)**2- 1)**2
            
        #~ def gs1(x, lambdaGradDist10x0=lambdaGradDist10x0, lambdaGradDist10x1=lambdaGradDist10x1):
        def gs1(x, lambdaGradDist10=lambdaGradDist10):
            J = zeros((1,2*rdim*numControlPoints))
            x0 = x[:rdim*numControlPoints]
            x1 = x[rdim*numControlPoints:]
            xvars = tuple([xs0]+[x0[i*rdim:(i+1)*(rdim)] for i in range(numControlPoints)] + [xg0] + 
            [xs1] +[x1[i*rdim:(i+1)*(rdim)] for i in range(numControlPoints)] + [xg1]) 
            
            grad = lambdaGradDist10(*xvars)
            if grad != 0.:
                for i, el in enumerate(grad):                
                    J[0,i*rdim:(i+1)*rdim] = el               
            return J
            
        def fu0(x,  pt0t=pt0t):
            x0 = x[:rdim*numControlPoints]
            pt0 = pt0 = pt0t(x0)
            return (pt0.T.dot(pt0) -1)**2
            
        def gs0(x, lambdaGradDistbase=lambdaGradDistbase):
            x0 = x[:rdim*numControlPoints]
            J = zeros((1,2*rdim*numControlPoints)) 
            xvars = tuple([xs0]+[x0[i*rdim:(i+1)*(rdim)] for i in range(numControlPoints)] + [xg0] )  
            grad = lambdaGradDistbase(*xvars)
            for i, el in enumerate(grad):         
                J[0,i*rdim:(i+1)*rdim] = el    
            return J
        res += [[fu0,gs0],[fu1,gs1]]
    return res
    

# first, try to generate bezier that enforces length 

#~ def discrete_bezier(b, num_steps):
    #~ fns = float(num_steps)
    #~ return [b (float(i) / fns) for i in range(num_steps + 1)]

##############" BEZIER STUFF #############


def genFJ(cons, sizejac):
    def F(x):
        F = zeros(0);
        for (fi, Ji) in cons:
            F = hstack([F,fi(x)])
        return F    
    def J(x):
        Jx = zeros((0,sizejac))
        for (fi, Ji) in cons:
            Jx = vstack([Jx,Ji(x)])
        return Jx
    return F,J
        
def genf_df(F,J):
    def f(x, Fx=None):
        if Fx is not None:
            pass
        else:
            Fx = F(x)
        return 0.5* Fx.T.dot(Fx)
        
    def df(x, Fx=None, Jx = None):
        if Fx is not None:
            pass
        else:
            Fx = F(x)
        if Jx is not None:
            pass
        else:
            Jx = J(x)
        return Fx.T.dot(Jx)
    return f, df
        
        
def dx(f, x):
    return abs(0-f(x))
 
def backtrack(f, fx, dfx, x, p, tau = 0.5, c = 0.5, alpha = 10.): #p is search direction
    m = p.T.dot(dfx)
    t = -c * m
    while not ((fx) - f(x + alpha * p) >= alpha * t) :
        alpha = tau * alpha
    #~ print "alpha", alpha
    return alpha
        
def stepC(x, eps = 1.5, hard = [constraint("pos")], soft = []):
    
    #calling appropriate constraints
    sizejac = x.shape[0]
    F,J  = genFJ(hard, sizejac)
    G,JG = genFJ(soft, sizejac)
    
    Fx = F(x)
    Jx = J(x)
    Gx = G(x)
    JGx = JG(x)
    
    f, df = genf_df(F,J)
    
    fx  =  f(x, Fx)
    dfx    = df(x, Fx, Jx)
    dfxinv = pinv(dfx.reshape((-1,1))).reshape((-1,))
    
    eps = backtrack(f, fx, dfx, x, -dfxinv * fx)
    
    #~ print "grad_FT_D", norm(F)
    #~ nfx = norm(Fx) + norm(Gx)
    #~ if nfx >= 0.0001:   
    print "fx",  fx
    if eps >= 0.0001:    
        Ji = pinv(Jx)
        JiJ = Ji.dot(Jx)
        JGi = pinv(JGx)
        idnull = identity(JiJ.shape[0])
        #~ return x -  eps * pinv(J).dot(F) - eps * (idnull - JiJ).dot(JGi).dot(G)
        #~ print "nx", x -  eps * dfxinv * fx - eps * (idnull - JiJ).dot(JGi).dot(Gx)
        return x -  eps * dfxinv * fx - eps * (idnull - JiJ).dot(JGi).dot(Gx)
    return x
    
    
def getLineFromSegment(line):
        a = line[0]; b = line[1]; c = a.copy() ; c[2] = 1.
        normal = cross((b-a),(c-a))
        normal /= norm(normal)
        # get inequality
        dire = b - a
        coeff = normal
        rhs = a.dot(normal)
        return (coeff, array([rhs]))

from hpp_spline import bezier
from  cord_methods import *
from  plot_cord import *



def bezierFromVal(xis, numvars = 3):
    wps = zeros((len(xis),3))
    for i in range(numvars):
        wps[i,:rdim] = array(xis[i])
    return bezier(wps.transpose(), 1.)

from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    A = zeros(2); A[1] = 1
    A = A.reshape((1,2))
    b = 1.26
    
    #~ A = array([ 0.707, -0.707]); A = A.reshape((1,2))
    #~ b = array([-0.354])
    
    A = array([ 0.97 ,  0.243,0.]); A = A.reshape((1,rdim))
    b = array([1.213])
    #~ b = array([1.019])
    
    #~ (array([ 0.97 ,  0.243, -0.   ]), array([1.019]))
    
    toggle3D()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    setFigure(ax)
    
    def ik(x_end):        
        hard = [constraint("pos")]
        #~ soft = [constraint("target",x_end),constraint("ineq",A,b)]
        soft = [constraint("target",x_end)]
        
        x = zeros(rdim*2); x[:rdim]= [0.,1.2,0.]
        x[rdim:]= [0.,2.,0.]
        for i in range(100):
            xn = stepC(x, 0.1,hard, soft)
            if norm(xn-x) <= 0.0001:
                print "optimium local"
                break
            x = xn
            
        xis = array([[0,0,0.],x[:rdim],x[rdim:]])
            
        plotPoints(xis, color = "b")
        plotPoints(array([x_end]), color = "r")
        plotSegment(xis[:2], color = "b")
        plotSegment(xis[1:], color = "b")
        return x
    
    xs = ik([1.2,0.2,0.])
    xg = ik([1.,1.,0.])
                     
    plt.show()
    
    
    nvars = 4
    hard = bezier2Dist(xs,xg,nvars) + [constraint("ineq",A,b)]
    #~ hard = bezier2Dist(xs,xg,nvars)
    #~ x = xs[:]
    x = zeros(nvars*rdim*2)
    for i in range(nvars):
        if i < nvars /2:
            x[i*rdim:(i+1)*rdim] = xs[:rdim]
        else:
            x[i*rdim:(i+1)*rdim] = xg[rdim:]
    for i in range(nvars,nvars*2):
        if i < nvars /2:
            x[i*rdim:(i+1)*rdim] = xs[:rdim]
        else:
            x[i*rdim:(i+1)*rdim] = xg[rdim:]
            
    print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ", x
    print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ", x.shape
    
    #~ for i in range(1000):
    print "x", x
    for i in range(500):
        xn = stepC(x, 0.1,hard, soft=[constraint("ineq",A,b)])
        #~ xn = stepC(x, 0.1,hard, soft=[])
        #~ x = stepC(x, 0.01,hard, soft=[])
        if norm(xn-x) <= 0.0001:
            print "optimium local"
            break
        x = xn
        
        
    xis = [x[i*rdim:(i+1)*rdim] for i in range(nvars)] + [x[i*rdim:(i+1)*rdim] for i in range(nvars, 2*nvars)]
    #~ x1 = [xs[rdim:]] + [x[i*rdim:(i+1)*rdim] for i in range(nvars, 2*nvars)] + [xg[rdim:]]
    #~ xis = array([[0,0],x[:rdim],x[rdim:]])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')    
    setFigure(ax1)
    plotPoints(array(xis), color = "g")
    
    
    #retrieve bezier curves
    x0 = [xs[:rdim]] + [x[i*rdim:(i+1)*rdim] for i in range(nvars)] + [xg[:rdim]]
    x1 = [xs[rdim:]] + [x[i*rdim:(i+1)*rdim] for i in range(nvars, 2*nvars)] + [xg[rdim:]]
    
    print "x0", x0
    print "x1", x1
    
    #plot bzier
    #~ b = bezierFromVal([xs[:rdim],x[:rdim],xg[:rdim]])
    b = bezierFromVal(x0,nvars+2)
    plotBezier(b, "g", label = "x0", linewidth = 2.0)
    
    #~ plt.show()
    
    b2 = bezierFromVal(x1,nvars+2)
    plotBezier(b2, "y", label = "x1", linewidth = 2.0)
    
    #~ plotPoints(array([x_end]), color = "r")
    xis = array([[1.2,0.2,0.],[1.,1.,0.]])
    plotSegment(xis, color = "y")
    #~ plotSegment(xis[1:], color = "b")
    plt.legend()
    plt.show()
    #try to animate
    
    import matplotlib.animation as animation
    
    #~ xlim(b(0.)[0], b(1.)[0])
    #~ ylim(b(0.)[1], b(1.)[1])
    
    dt = 1./30. #30 fps
    t = 0.

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-2, 2), ylim=(-2, 2))
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    line2, = ax.plot([], [], 'o-', lw=2)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    energy_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)

    def init():
        """initialize animation"""
        line.set_data([], [])
        line2.set_data([], [])
        time_text.set_text('')
        energy_text.set_text('')
        return line, line2, time_text, energy_text

    def animate(i):
        """perform animation step"""
        global t
        t = t+dt
        if(t/10.>1.11):
            t = 0.
        ti = min(t/10., 0.99)
        ti = max(ti,0.)
        #~ print "ti", t/10
        p1 = b(ti)[:rdim]
        p2 = b2(ti)[:rdim]
        x = [0,p1[0],p2[0]]
        y = [0,p1[1],p2[1]]
        #~ position = b(ti)[:rdim]
        
        #~ x2 = [0.,0.5]
        #~ y2 = [0.5,1.]
        x2 = [1.2,1.]
        y2 = [0.2,1.]
        
        line2.set_data(*(x2,y2))
        
        line.set_data(*(x,y))
        #~ time_text.set_text('time = %.1f' % pendulum.time_elapsed)
        #~ energy_text.set_text('energy = %.3f J' % pendulum.energy())
        return line, line2, time_text, energy_text

    # choose the interval based on dt and the time to animate one step
    from time import time
    t0 = time()
    animate(0)
    t1 = time()
    interval = 300 * dt - (t1 - t0)

    ani = animation.FuncAnimation(fig, animate, frames=300,
                                  #~ interval=interval, blit=True, init_func=init, repeat = True)
                                  interval=interval, blit=True, init_func=init, repeat = True, repeat_delay = 500)
                                  
    #~ ani.save('test_bezier_cons.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
                                  
    plt.show()
