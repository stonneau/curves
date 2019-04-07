
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array, identity     
eps =0.000001

from numpy.linalg import pinv, norm

degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2

#try to solve distance based ik for a simple manipulator with 2 DOF, each with length 1


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

f10 = genf(1,0)
f11 = genf2vars(1,1,0)
target = None


def gen_target(x_end):
    global target
    def f(x):
        x1 = x[rdim:] - x_end
        #~ print "x ", x
        #~ print "x_end ", x_end
        #~ print "x1 ", x1
        #~ print "x[:rdim] ", x[rdim:]
        return x1.dot(x1)
        #~ return (x1.T.dot(x1))
    target = f


def Jacob2(x,x_end):
    #jacobian is 2 * 2 matrix
    #~ J = zeros((3*rdim,2*rdim))
    J = zeros((3,2*rdim))
    # df10 / dx0
    x0 = x[:2]
    x1 = x[2:]
    x1mx0 = (x1-x0)
    x0mx1 = (x0-x1)
    
    df1dx0 = 4*x0*(x0.T.dot(x0) - 1)
    df1dx1 = zeros(rdim)

    df2dx0 = 4*(x0mx1)*(x1mx0.T.dot(x1mx0) -1)
    df2dx1 = 4*(x1mx0)*(x1mx0.T.dot(x1mx0) -1)
    
    J[0,:2] = df1dx0
    J[0,2:] = df1dx1
    J[1,:2] = df2dx0
    J[1,2:] = df2dx1
    J[2,2:] = 2*(x1-x_end)
    return J
    
#~ def Hessian(x,x_end):
    
    
def F(x):
    res = zeros(3)
    res[0] = f10(x)
    res[1] = f11(x)
    res[2] = target(x)
    return res

def wolfe1(x,alpha,f,fx,G):
    c1 = 0.1
    pk = -G
    return (f(x+alpha*pk) <= fx + c1 * alpha * pk.T.dot(G) + 0.01).all()
        
def wolfe2(x,alpha,f,G, x_end):
    c2 = 0.9
    pk = -G
    Jalpha = Jacob2(x+alpha*pk,x_end)
    g2 = pinv(Jalpha).dot(f(x+alpha*pk))
    return pk.T.dot(g2) >= c2 * pk.T.dot(G)

def linesearch(x,f,fx,G,x_end):
    c1 = 0.001
    c2 = 0.9
    alpha = 1.
    for i in range(100):
        #~ if(wolfe1(x,alpha,f,fx,G)) and (wolfe2(x,alpha,f,G,x_end)):
        if(wolfe1(x,alpha,f,fx,G)):
            return alpha
        alpha = alpha / 2.
    return 0
    
J0iJ0  = None  
J1i  = None  
idnull = None
J = None

def step(x, x_end, eps = 1.):
    global J0iJ0
    global J1i
    global idnull
    global J
    fx = F(x)
    nfx = norm(fx)
    if nfx >= 0.0009:
        J = Jacob2(x,x_end)
        J0 = J[:2,:]
        J1 = J[2:,:]
        J1i = pinv(J1)
        F0 = fx[:2]
        G0 = fx[2:]
        J0i = pinv(J0)
        J0iJ0 = J0i.dot(J0)
        idnull = identity(J0iJ0.shape[0])
        G = pinv(J).dot(fx)
        return x -  eps * pinv(J0).dot(F0) - eps * (idnull - J0iJ0).dot(J1i).dot(G0)
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

if __name__ == '__main__':
    
    #~ x_end = [-1,-1.]
    x_end = [1.2,1.2]
    gen_target(x_end)
    
    x = zeros(4); x[:rdim]= [0.,1.2]
    x[rdim:]= [0.,2.]
    #~ for i in range(5):
        #~ x = step(x,x_end, 1)
    for i in range(100):
        x = step(x,x_end, 0.1)
    #~ for i in range(10):
        #~ x = step(x,0.01)
        
    xis = array([[0,0],x[:rdim],x[rdim:]])
        
    plotPoints(xis, color = "b")
    plotPoints(array([x_end]), color = "r")
    plotSegment(xis[:2], color = "b")
    plotSegment(xis[1:], color = "b")
    
               
    #~ plt.legend()
    #~ plt.show()
    
