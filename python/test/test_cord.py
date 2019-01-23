from spline import *

from numpy import matrix, array, zeros, ones
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
        if type(val) == str:
                return (_I3, _zeroVec)
        else:
                return (_zeroMat, val.reshape([3,1]))
        

def createWaypointList(waypoints):
        mat = zeros([3, len(waypoints)*3])
        vec = zeros([3, len(waypoints)])
        for i, val in enumerate(waypoints):
                mvar, vvar = createControlPoint(val)
                mat [:,i*3:i*3+3] = mvar
                vec [:,i:i+1] = vvar
        return mat, vec
        
def toPoint(val, mat, vec, i):
        return mat[i*3:(i+1)*3,:].dot(val)+vec[i*3:(i+1)*3]       

class varBezier:
        #waypoints is a list that contains either 3d arrays (constants), or a string "variable"
        def __init__(self, waypoints=None, time=1.):
                if(waypoints != None):
                        mat, vec = createWaypointList(waypoints)
                        self.bezier = bezierVar(mat, vec, time)
                
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
        def waypoints(self, varId=-1):
                if varId < 0:
                        return self.bezier.waypoints().A, self.bezier.waypoints().b
                assert self.bezier.nbWaypoints > varId
                mat = self.bezier.waypoints().A[:, varId*3:varId*3+3]
                vec = self.bezier.waypoints().b[:, varId]
                return mat, vec
                
        def matrixFormWaypoints(self, varId):
                assert(varId >= 0)
                mat, vec = self.waypoints(varId)
                #~ resmat = zeros([mat.shape[0],mat.shape[0]])
                resvec = zeros(3)
                for i in range(0, mat.shape[0]/3, 1):
                        #~ resmat[i*3:i*3+3, i*3:i*3+3] = mat[i*3:i*3+3, :]
                        resvec += vec[i*3:i*3+3]
                return mat.transpose(), resvec
             
        
                
        def toBezier3(self, x):
                wps = []
                for i in range(self.bezier.nbWaypoints):
                        mat, vec = self.matrixFormWaypoints(i) 
                        wps += [mat.dot(x) + vec]
                        #~ mat, vec = self.waypoints(varId)
                #~ assert(len(x_list)*3 == self.bezier.waypoints().b.shape[0])
                #~ wps = []
                #~ for i in range(self.bezier.nbWaypoints):
                        #~ mat, vec = self.waypoints(i)
                        #~ pt = array([0.,0.,0.])
                        #~ for j, el in enumerate(x_list):
                                #~ pt+=toPoint(el,mat,vec,j)
                        #~ wps += [pt]
                return bezier(array(wps).transpose(),self.bezier.max())

#~ def square(mat, vec):
        #~ resmat = zeros([mat.shape[0],mat.shape[0]])
        #~ resvec = zeros(3)
        #~ for i in range(0, mat.shape[0], 3):
                #~ resmat[i*3:i*3+3, i*3:i*3+3] = mat[i*3:i*3+3, :]
                #~ resvec += [i*3:i*3+3]
        

def segmentConstraint(varBez, a, b, wpIndex, constraintVarIndex, totalAddVarConstraints):
        mat, vec = varBez.matrixFormWaypoints(wpIndex)
        vec = b - vec   
        resmat = zeros([mat.shape[0],mat.shape[1]+totalAddVarConstraints])
        resmat[:,:mat.shape[1]]=mat
        resmat[:, mat.shape[1]+constraintVarIndex-1]= b-a
        return resmat, vec
        #mat X' + vec = (alpha) a + (1 - alpha) b = alpha(a-b) + b
        #[mat, (a-b)] [X, alpha]' + b 
                      

#try to solve the qp with one segment constraint

testConstant = varBezier([array([1,2,3]),"","",array([4,5,5])], 1.)
#~ testConstant = varBezier([""], 1.)
subs = testConstant.split([0.4])

dim = testConstant.bezier.nbWaypoints

a = array([2.,3.,2.])
b = array([0.,0.,1.])

c0 = segmentConstraint(subs[0],a,b,dim-1,1,1)

#~ a=varBezier([array([1,2,3]),""], 1.)
#~ subs = a.split([0.4])
from qp import solve_lp
q = zeros(c0[0].shape[1])
q[-1] = -1
G = zeros([2,q.shape[0]])
h = zeros(2)
G[0,-1] =  1 ; h[0]=1.
G[1,-1] = -1 ; h[1]=0.
C = c0[0]
d = c0[1]
res = solve_lp(q, G=G, h=h, C=c0[0], d=c0[1])
x_list = [res[i:i+3] for i in range(dim)]
test =testConstant.toBezier3(res[:-1])
testsub =subs[0].toBezier3(res[:-1])

#~ testConstant=varBezier([array([1,2,3]),"","",array([4,5,5])], 1.)
#~ subs = testConstant.split([0.4])

#~ zero3 = [zeros(3) for _ in range(3)]
#~ ones3 = [ones(3) for _ in range(3)]
#~ ones4 = ones(12)
#~ #check that t works

#~ b = testConstant.toBezier3(ones4)

#~ b_sub0 =  subs[0].toBezier3(ones4)
#~ b_sub1 =  subs[1].toBezier3(ones4)

#~ for i in range(10):
        #~ t = i / 10. * 0.4
        #~ assert(abs(b_sub0(t) - b(t)) <= __EPS).all()
#~ for i in range(10):
        #~ t = i / 10. * 0.6
        #~ assert(abs(b_sub1(t) - b(t + 0.4)) <= __EPS).all()




#get points
#~ beziersub1wps = 

import numpy as np
import matplotlib.pyplot as plt

step = 1000.
points =  np.array([(test(i/step)[0][0],test(i/step)[1][0]) for i in range(int(step))])

#~ points = np.array([(0, 1), (2, 4), (3, 1), (9, 3)])
# get x and y vectors
x = points[:,0]
y = points[:,1]

# calculate polynomial
#~ z = np.polyfit(x, y, 2)
#~ f = np.poly1d(z)

# calculate new x's and y's
#~ x_new = np.linspace(x[0], x[-1], 50)
#~ y_new = f(x_new)

#~ plt.plot(x,y,'o', x, y)
plt.plot(x,y)
#~ plt.xlim([x[0]-10, x[-1] + 1 ])
plt.show()


