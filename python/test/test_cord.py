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


from convex_hull import *

lines0 = genFromLine(lines0[0], 5, [[0,6],[0,7]])

from qp import solve_lp, quadprog_solve_qp

def getRightMostLine(ptList):
        pt1 = array([0.,0.,0.])
        pt2 = array([0.,0.,0.])
        for pt in ptList:
                if pt[0] > pt1[0]:
                       pt1 =  pt
                elif pt[0] > pt2[0]:
                       pt2 =  pt
        if pt1[1] < pt2[1]:
                return [pt2,pt1]
        else:
                return [pt1,pt2]

def plotBezier(bez, color):
        step = 100.
        points1 =  np.array([(bez(i/step*bez.max())[0][0],bez(i/step*bez.max())[1][0]) for i in range(int(step))])
        x = points1[:,0]
        y = points1[:,1]
        plt.plot(x,y,color,linewidth=2.0)
        
def plotControlPoints(bez, color):
        wps = bez.waypoints()
        wps = np.array([wps[:2,i] for i in range(wps.shape[1]) ])
        x = wps[:,0]
        y = wps[:,1]
        plt.scatter(x,y,color=color)
        #~ [plt.plot(el[0],el[1]) for el in points1]
        
        
def plotPoly(lines, color):
        step = 1000.
        for line in lines:
                a_0 = line[0]
                b_0 = line[1]
                pointsline =  np.array([ a_0 * i / step + b_0 * (1. - i / step) for i in range(int(step))])
                xl = pointsline[:,0]
                yl = pointsline[:,1]
                plt.plot(xl,yl,color,linewidth=0.5)

idxFile = 0
import string
import uuid; uuid.uuid4().hex.upper()[0:6]
#~ def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        #~ return ''.join(np.random.choice(string.ascii_uppercase + string.digits) for _ in range(size))



def accelerationcost(bezVar):
        # derivate twice for acceleration
        acc = bezVar.bezier.compute_derivate(2)
        #~ prim = acc.compute_primitive(1)
        
        
        
        hacc = varBezier(); hacc.bezier = acc
        wps_i=[[],[],[]]
        for i in range(3): #integrate for each axis
                for j in range(acc.nbWaypoints):
                        A_i = hacc.matrixFormWaypoints(j)[0][i,:].reshape([1,-1])
                        b_i = hacc.matrixFormWaypoints(j)[1][i]
                        wps_i[i] += [(A_i.transpose().dot(A_i), b_i*b_i) ]
        # now integrate each bezier curve
        for i, wps in enumerate(wps_i):
              wps_i[i] = compute_primitive(wps)
        
        resmat = wps_i[0][-1][0]; resvec = wps_i[0][-1][1]
        for i in range(1,3):
                resmat = resmat + wps_i[i][-1][0]
                resvec = resvec + wps_i[i][-1][1]
        return (resmat, resvec)

        
        #~ endwp = hprim.matrixFormWaypoints(prim.nbWaypoints-1)
        
        #~ return hprim
        # integral is given by difference of end control points of primitive

def computeTrajectory(bezVar, splits, save, filename = uuid.uuid4().hex.upper()[0:6]):
        global idxFile
        colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        subs = bezVar.split(splits)
        #generate random constraints for each line
        line_current = [array([1.,1.,0.]), array([0.,1.,0.])]
        
        #qp vars
        P = accelerationcost(bezVar)[0]; P = P + identity(P.shape[0]) * 0.0001
        #~ print "P", P
        q = zeros(3*bezVar.bezier.nbWaypoints)
        q[-1] = -1
        G = zeros([2,q.shape[0]])
        h = zeros(2)
        G[0,-1] =  1 ; h[0]=1.
        G[1,-1] = -1 ; h[1]=0.
        
        dimExtra = 0
        
        for i, bez in enumerate(subs):
                color = colors[i]
                init_points = []
                if i == 0:
                        init_points = [bezVar.waypoints()[1][:3,0][:2]]
                if i == len(subs)-1:
                        init_points = [bezVar.waypoints()[1][-3:,-1][:2]]
                lines, ptList = genFromLine(line_current, 5, [[0,5],[0,5]],init_points)
                matineq0 = None; vecineq0 = None
                for line in lines:
                        (mat,vec) = getLineFromSegment(line)
                        (mat,vec) = lineConstraint(bez, mat,vec,dimExtra)
                        (matineq0, vecineq0) = (concat(matineq0,mat), concatvec(vecineq0,vec))
                ineq  = (matineq0, vecineq0)
                G = vstack([G,ineq[0]])
                h = concatvec(h,ineq[1])
                line_current = getRightMostLine(ptList)
                plotPoly  (lines, color)
        C = None; d = None
        try:
                #~ res = solve_lp(q, G=G, h=h, C=C, d=d)
                res = quadprog_solve_qp(P, q, G=G, h=h, C=C, d=d)
                #plot bezier
                for i, bez in enumerate(subs):
                        color = colors[i]
                        test = bez.toBezier3(res[:])
                        plotBezier(test, color)
                if save:
                        plt.savefig(filename+str(idxFile))
                #plot subs control points
                for i, bez in enumerate(subs):
                        color = colors[i]
                        test = bez.toBezier3(res[:])
                        plotBezier(test, color)
                        plotControlPoints(test, color)
                if save:
                        plt.savefig("subcp"+filename+str(idxFile))
                final = bezVar.toBezier3(res[:])
                plotControlPoints(final, "black")
                if save:
                        plt.savefig("cp"+filename+str(idxFile))
                idxFile += 1
                if save:
                        plt.close()
                else:
                        plt.show()
                
                return final
        except ValueError:
                #~ pass
                plt.close()
     
     
#~ bezier_curve_t compute_primitive(const std::size_t order) const
    #~ {
        #~ if(order == 0) return *this;
        #~ num_t new_degree = (num_t)(degree_+1);
        #~ t_point_t n_wp;
        #~ point_t current_sum =  point_t::Zero(Dim);
        #~ // recomputing waypoints q_i from derivative waypoints p_i. q_0 is the given constant.
        #~ // then q_i = (sum( j = 0 -> j = i-1) p_j) /n+1
        #~ n_wp.push_back(current_sum);
        #~ for(typename t_point_t::const_iterator pit =  pts_.begin(); pit != pts_.end(); ++pit)
        #~ {
            #~ current_sum += *pit;
            #~ n_wp.push_back(current_sum / new_degree);
        #~ }
        #~ bezier_curve_t integ(n_wp.begin(), n_wp.end(),T_, mult_T_*T_);
        #~ return integ.compute_primitive(order-1);
    #~ }
    
def compute_primitive(wps):
        new_degree = (len(wps)-1+1)
        current_sum = [zeros(wps[0][0].shape),0]
        n_wp = [(current_sum[0],0.)]
        for wp in wps:
                current_sum[0] = current_sum[0] + wp[0]
                current_sum[1] = current_sum[1] + wp[1]
                #~ print "current_sum[0]", current_sum[0].shape
                #~ print "current_sum[1]", current_sum[1]
                n_wp +=[(current_sum[0]/new_degree, current_sum[1]/new_degree)]
        return n_wp
        
        
        
                
def genBezierInput(numvars = 3):
        valDep = array([np.random.uniform(0., 1.), np.random.uniform(0.,5.), 0.])
        valEnd = array([np.random.uniform(5., 10.), np.random.uniform(0.,5.), 0.])
        return varBezier([valDep]+["" for _ in range(numvars)]+[valEnd], 1.)
        
def genSplit(numCurves):
        splits = []
        lastval = np.random.uniform(0., 1.)
        for i in range(numCurves):
                splits += [lastval]
                lastval += np.random.uniform(0., 1.)
        return [el / lastval for el in splits[:-1]]
                


def gen(save = False):
        testConstant = genBezierInput(4)
        splits = genSplit(4)
        return computeTrajectory(testConstant,splits, save), testConstant

res = None
for i in range(1000):
        res = gen()
        if res[0] != None:
                break


