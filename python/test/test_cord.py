from spline import *

from numpy import matrix, array, zeros
from numpy.linalg import norm

__EPS = 1e-6

import eigenpy
eigenpy.switchToNumpyArray()

waypointsA = array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.],[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]).transpose()
waypointsb = array([[1.,2.,3.],[1.,2.,3.]]).transpose()

#testing bezier curve
a = bezierVar(waypointsA, waypointsb, 3.)

subBeziers = a.split(array([[0.2,0.4]]).transpose())
assert(subBeziers.size == 3)
assert(subBeziers.at(0).max() - 0.2 <= __EPS)
assert(subBeziers.at(1).max() - 0.2 <= __EPS)
assert(subBeziers.at(2).max() - 2.6 <= __EPS)


_zeroMat = array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]).transpose()
_I3 = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).transpose()
_zeroVec = array([[0.,0.,0.]]).transpose()

def createControlPoint(val):
        print type(val)
        print val
        if type(val) == str:
                return (_I3, _zeroVec)
        else:
                return (_zeroMat, val.reshape([3,1]))
        

def createWaypointList(waypoints):
        mat = zeros([3, len(waypoints)*3])
        vec = zeros([3, len(waypoints)])
        for i, val in enumerate(waypoints):
                print "wt", val
                print "wt2", val
                mvar, vvar = createControlPoint(val)
                mat [:,i*3:i*3+3] = mvar
                vec [:,i:i+1] = vvar
        return mat, vec

class varBezier:
        #waypoints is a list that contains either 3d arrays (constants), or a string "variable"
        def __init__(self, waypoints=None, time=1.):
                if(waypoints != None):
                        mat, vec = createWaypointList(waypoints)
                        self.bezier = bezierVar(waypointsA, waypointsb, time)
                
        def fromBezier(self, bez):
                var = varBezier ()
                var.bezier = bez
                return var
        
        #splits into n+1 continuous curves, with n the number of time values
        def split(self, times):
                timearray = array(times).reshape([-1,1])
                subBeziers = self.bezier.split(timearray)
                dim = subBeziers.size                
                return [self.fromBezier(subBeziers.at(i)) for i in range(dim)]
                
        # for each control point of the curve, gives the linear depency 
        # of a given variable
        def waypoints(self, varId):
                if varId < 0:
                        return self.bezier.waypoints().A, self.bezier.waypoints().b
                assert self.bezier.nbWaypoints > varId
                mat = self.bezier.waypoints().A[:, varId*3:varId*3+3]
                vec = self.bezier.waypoints().b[:, varId]
                return mat, vec
                
                
a=varBezier([array([1,2,3]),""], 1.)
subs = a.split([0.4])
