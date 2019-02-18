from numpy import polynomial
from math import sqrt

from scipy.special import binom as bino

def bern(i, n):
    def eva(t):
        return bino(n,i) * t**i * (1-t)**(n-i)
    return eva

def aeq(a,b):
    return (a**2) -(b**2) < 1e-6
    

def productBern(i0,n0, i1, n1):
    pos =  [(i,j) for i in range(n0+1) for j in range(n1+1)]
    pos = pos.index((i0,i1))
    weight = 1; ratio = 0
    notfound = True
    pos1 = i0 * (n1+1) + i1
    #~ print "pos1 ", pos1
    #~ print "n0 + n1", n0 + n1
    #~ print "i0+i1", i0+i1
    #~ val =   float(i0+i1)/ float(n0+n1) 
    val = float(bino(n0,i0)*bino(n1,i0-i1)) / float(bino(n0+n1,i0))
    #~ print "top", float(bino(n0,i0)*bino(n1,i0-i1))
    #~ print "dow", float(bino(n0+n1,i0))
    #~ print "val", val
    #~ if(val < 0 ):
    if(i0 > i1):
        val = float(bino(n1,i1)*bino(n0,i1-i0)) / float(bino(n0+n1,i1))
        #~ print "val 2!!!!!!!!!!", val
    print "val", val
    assert(pos1 == pos)
    res = None
    def eval0(t):
        return bern(i0, n0)(t) * bern(i1, n1)(t)
    bn =  bern(i0+i1, n0+n1)
    def eval1(t):
        return  (val *  bn(t))
    def eval2(t):
        #~ print "val", 1 - val
        return  val2 *  bn(t)
    eva = eval1    
    res = [(i0,i1),val]
    for i in range(10):            
        t = (float)(i) / 10.
        if not aeq(eval1(t), eval0(t)):
            print "i0, i1", i0, i1
            if not aeq(eval2(t), eval0(t)):
                print 'failed 0 at', t, eva(t), eval0(t)
                return False
    return True, res
            

