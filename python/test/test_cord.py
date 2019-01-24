from spline import *

from numpy import matrix, array, zeros, ones, diag, cross
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
        return (resmat, vec)
        #mat X' + vec = (alpha) a + (1 - alpha) b = alpha(a-b) + b
        #[mat, (a-b)] [X, alpha]' + b 

from numpy import vstack, identity
def concat(m1,m2):
        if (type(m1) == type(None)):
                return m2
                
        return vstack([m1,m2]).reshape([-1,m2.shape[-1]])
        
def concatvec(m1,m2):
        if (type(m1) == type(None)):
                return m2
        return array(m1.tolist()+m2.tolist())

def lineConstraint(varBez, C, d, totalAddVarConstraints):
        resmat = None
        resvec = None
        for wpIndex in range(varBez.bezier.nbWaypoints):
                mat, vec = varBez.matrixFormWaypoints(wpIndex)
                mat = C.dot(mat)
                vec = d - C.dot(vec)
                resmat = concat(resmat, mat)
                resvec = concatvec(resvec, vec)
        augmented = zeros([resmat.shape[0],resmat.shape[1]+totalAddVarConstraints])
        augmented[:,:resmat.shape[1]]=resmat
        return (augmented, resvec)

#try to solve the qp with one segment constraint

testConstant = varBezier([array([1,2,0]),"","",array([4,5,0])], 1.)
#~ testConstant = varBezier([""], 1.)
subs = testConstant.split([0.4])

dim = testConstant.bezier.nbWaypoints

a = array([2.,3.,0.])
b = array([0.,0.,0.])

def inequality(v, n): 
	#the plan has for equation ax + by + cz = d, with a b c coordinates of the normal
	#inequality is then ax + by +cz -d <= 0 
	# last var is v because we need it
	return [n[0], n[1], n[2], np.array(v).dot(np.array(n)) + __EPS]

# constraint is left of line
lines0 = [[array([1.1,0.,0.]), array([1.1,4.,0.])], [array([1.,1.5,0.]), array([2.,1.5,0.])]]
def getLineFromSegment(line):
        a = line[0]; b = line[1]; c = a.copy() ; c[2] = 1.
        normal = cross((b-a),(c-a))
        normal /= norm(normal)
        # get inequality
        dire = b - a
        coeff = normal
        rhs = a.dot(normal)
        return (coeff, array([rhs]))

eq0 = segmentConstraint(subs[0],a,b,dim-1,1,1)
matineq0 = None; vecineq0 = None

from convex_hull import *

lines0 = genFromLine(lines0[0], 5, [[0,6],[0,7]])

for line in lines0:
        (mat,vec) = getLineFromSegment(line)
        (mat,vec) = lineConstraint(subs[0], mat,vec,1)
        (matineq0, vecineq0) = (concat(matineq0,mat), concatvec(vecineq0,vec))
ineq0 = (matineq0, vecineq0)

#~ lineConstraint(varBez, C, d, totalAddVarConstraints)
#~ a=varBezier([array([1,2,3]),""], 1.)
#~ subs = a.split([0.4])
from qp import solve_lp
q = zeros(eq0[0].shape[1])
q[-1] = -1
G = zeros([2,q.shape[0]])
h = zeros(2)
G[0,-1] =  1 ; h[0]=1.
G[1,-1] = -1 ; h[1]=0.
G = vstack([G,ineq0[0]])
h = concatvec(h,ineq0[1])
C = eq0[0]
d = eq0[1]
#~ C = None; d = None
res = solve_lp(q, G=G, h=h, C=C, d=d)
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

step = 100.
points1 =  np.array([(test(i/step*0.4)[0][0],test(i/step*0.4)[1][0]) for i in range(int(step))])
points2 =  np.array([(test(i/step*0.6+0.4)[0][0],test(i/step*0.6+0.4)[1][0]) for i in range(int(step))])

#~ lines0 = genFromLine(lines0[0], 50, [[0,6],[0,7]])
#~ step = 10.
for line in lines0:
        a_0 = line[0]
        b_0 = line[1]
        pointsline =  np.array([ a_0 * i / step + b_0 * (1. - i / step) for i in range(int(step))])
        xl = pointsline[:,0]
        yl = pointsline[:,1]
        plt.plot(xl,yl,'b')
#~ points = np.array([(0, 1), (2, 4), (3, 1), (9, 3)])
# get x and y vectors
x = points1[:,0]
y = points1[:,1]
x2 = points2[:,0]
y2 = points2[:,1]


#~ for line in lines:
        

# calculate polynomial
#~ z = np.polyfit(x, y, 2)
#~ f = np.poly1d(z)

# calculate new x's and y's
#~ x_new = np.linspace(x[0], x[-1], 50)
#~ y_new = f(x_new)

plt.plot(x,y,'b')
plt.plot(x2,y2,'r')

#~ plt.plot(x,y)
#~ plt.xlim([x[0]-10, x[-1] + 1 ])
plt.show()


