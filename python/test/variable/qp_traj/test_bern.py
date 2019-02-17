from numpy import polynomial
from math import sqrt

from scipy.special import binom as bino

def bern(i, n):
    def eva(t):
        return bino(n,i) * t**i * (1-t)**(n-i)
    return eva

def aeq(a,b):
    return sqrt(a**2) - sqrt(b**2) < 1e-9
    

def productBern(i0,n0, i1, n1):
    pos =  [(i,j) for i in range(n0+1) for j in range(n1+1)]
    pos = pos.index((i0,i1))
    def eval0(t):
        return bern(i0, n0)(t) * bern(i1, n1)(t)
    bn =  bern(i0+i1, n0+n1)
    def eval1(t):
        return  ((i0+i1)/ (n0+n1)) *  bn(t)
    def eval2(t):
        return  1 - (i0+i1)/ (n0+n1) *  bn(t)
    eva = eval1
    if(pos)%2 == 0:
        eva = eval2    
    for i in range(10):            
        t = (float)(i) / 10.
        if not aeq(eval1(t), eval0(t)):
            print 'failed 0 at', t, eval1(t), eval0(t)
            return False
    return True
            
    
