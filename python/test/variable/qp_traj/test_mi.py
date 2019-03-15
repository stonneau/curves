
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array, identity
eps =0.000001


degree = 5
n = degree
nvars = degree+1
dim = 3
rdim = 2

if __name__ == '__main__':
    
    x = cp.Variable(2, boolean=True)
    y = cp.Variable(1)
    #~ eq1 = y <= 2
    eq2 = y >= -2
    eq1 = y >= -1
    
    M = 5000
    m = -5000
    
    constraints = [
    #~ m*(1-x[0]) <= y +2 ,
    #~ m*(1-x[0]) <= y +1 ,    
    -y - 2 <= M*(1-x[0]),
    -y - 1 <= M*(1-x[1]),    
    x[0] + x[1] == 1
    ];
    
    #~ if x[0] then y >= -2,  x[0] * y >= -1, x[0] + x[1] == 1]

    obj = cp.Minimize(y)
    prob = cp.Problem(obj, constraints)
    res = prob.solve(verbose=True)
    print res
    
