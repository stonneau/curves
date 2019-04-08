
import cvxpy as cp
from numpy.random import randn
from numpy import matrix, array, zeros, ones, diag, cross, vstack, identity
eps =0.000001

from numpy.linalg import pinv, norm

degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2

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
    global softC
    def f(x):
        x1 = x[rdim:] - x_end
        return x1.dot(x1)
    return f

#distance from base to joint 1 
f10 = genf(1,0)

#distance from joint 2 to joint 1 
f12 = genf2vars(1,1,0)


def positions(x):
    res = zeros(2)
    res[0] = f10(x)
    res[1] = f12(x)
    return res

#target position constraint
#~ target = None

#~ def gen

##############" END CONSTRAINTS #############


##############" JACOBIANS #############


###### HARD CONSTRAINTS JACOBIANS ######

#position

def Jpositions(x):
    J = zeros((2,2*rdim))
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
    return J
  
###### SOFT CONSTRAINTS JACOBIANS ######  

def Jtarget(x_end):  
    def fu(x):
        J = zeros((1,2*rdim)) 
        J[0,2:] = 2*(x[2:]-x_end) 
        return J
    return fu
    
##############" CONSTRAINT DIC #############
    
hardC = { "pos" : [positions, Jpositions] }

softC = {"target" : [gen_target, Jtarget] } 

def targetConstraint(x_end):
    target = softC["target"]
    return [target[0](x_end), target[1](x_end)]
    
    
def F(x):
    res = zeros(3)
    res[0] = f10(x)
    res[1] = f12(x)
    res[2] = target(x)
    return res
    
def stepC(x, x_end, eps = 1., hard = [hardC["pos"]], soft = []):
    
    #calling appropriate constraints
    F = zeros(0);  G = zeros(0); J = zeros((0,4)); JG = zeros((0,4))
    for (fi, Ji) in hard:
        F = hstack([F,fi(x)])
        J = vstack([J ,Ji(x)])
    for (gi, Ji) in soft:
        G  = hstack([G,gi(x)])
        JG = vstack([JG ,Ji(x)])    
    #evluation
    nfx = norm(F) + norm(G)
    if nfx >= 0.0009:    
        Ji = pinv(J)
        JiJ = Ji.dot(J)
        JGi = pinv(JG)
        idnull = identity(JiJ.shape[0])
        return x -  eps * pinv(J).dot(F) - eps * (idnull - JiJ).dot(JGi).dot(G)
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
    
    def ik(x_end):        
        hard = [hardC["pos"]]
        soft = [targetConstraint(x_end)]
        
        x = zeros(4); x[:rdim]= [0.,1.2]
        x[rdim:]= [0.,2.]
        for i in range(100):
            x = stepC(x,x_end, 0.1,hard, soft)
            
        xis = array([[0,0],x[:rdim],x[rdim:]])
            
        plotPoints(xis, color = "b")
        plotPoints(array([x_end]), color = "r")
        plotSegment(xis[:2], color = "b")
        plotSegment(xis[1:], color = "b")
    
    ik([1.2,1.2])
    ik([2,2])
    
               
    #~ plt.legend()
    plt.show()
    
