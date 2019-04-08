
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
    def f(x):
        x1 = x[rdim:] - x_end
        return x1.dot(x1)
    return f

# A x <= b
def gen_ineq(A,b):
    def f(x):
        x0 = x[:rdim]
        x1 = x[rdim:]
        res = zeros(rdim)
        res[0] = max(A.dot(x0) - b,0.)
        res[1] = max(A.dot(x1) - b,0.)
        return res
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
    
    
def Jineq(A,b):
    Arows = A.shape[0]
    fineq  = gen_ineq(A, b)
    def f(x):
        fin = fineq(x)
        J = zeros((2*Arows,2*rdim))
        if(fin[0]>0.): 
            J[0:Arows,:2] = A[:]
        if(fin[1]>0.):
            J[Arows: ,2:] = A[:]
        return J
    return f
  
###### SOFT CONSTRAINTS JACOBIANS ######  

def Jtarget(x_end):  
    def fu(x):
        J = zeros((1,2*rdim)) 
        J[0,2:] = 2*(x[2:]-x_end) 
        return J
    return fu
    
##############" CONSTRAINT DIC #############
    
all_constraints = { "pos" : [positions, Jpositions], "target" : [gen_target, Jtarget], "ineq" : [gen_ineq, Jineq]}

def constraint(name,*args):
    cons = all_constraints[name]
    return [cons[0](*args), cons[1](*args)]
    
##############" CONSTRAINT DIC #############
        
def stepC(x, x_end, eps = 1., hard = [all_constraints["pos"]], soft = []):
    
    #calling appropriate constraints
    F = zeros(0);  G = zeros(0); J = zeros((0,4)); JG = zeros((0,4))
    for (fi, Ji) in hard:
        F = hstack([F,fi(x)])
        J = vstack([J ,Ji(x)])
    for (gi, Ji) in soft:
        G  = hstack([G,gi(x)])
        JG = vstack([JG ,Ji(x)])    
    #evluation
    #~ print "F", F
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
    A = zeros(2); A[1] = 1
    A = A.reshape((1,2))
    b = 1.26
    ineq = ineqConstraint(A,b)
    
    def ik(x_end):        
        hard = [all_constraints["pos"]]
        soft = [constraint("target",x_end),constraint("ineq",A,b)]
        
        x = zeros(4); x[:rdim]= [0.,1.2]
        x[rdim:]= [0.,2.]
        for i in range(100):
            x = stepC(x,x_end, 0.1,hard, soft)
            
        xis = array([[0,0],x[:rdim],x[rdim:]])
            
        plotPoints(xis, color = "b")
        plotPoints(array([x_end]), color = "r")
        plotSegment(xis[:2], color = "b")
        plotSegment(xis[1:], color = "b")
    
    ik([1.2,0.2])
    ik([2,2])
    
               
    #~ plt.legend()
    plt.show()
    
