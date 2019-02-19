from hpp_spline import *
from numpy import matrix, array, zeros, ones, diag, cross
from numpy.linalg import norm

__EPS = 1e-6


from varBezier import varBezier
from qp_cord import *
from plot_cord import *


idxFile = 0
import string
import uuid; uuid.uuid4().hex.upper()[0:6]


######################## generate problems ########################

def genSplit(numCurves):
        if (numCurves == 0):
                return [-1]
        splits = []
        lastval = np.random.uniform(0., 1.)
        for i in range(numCurves):
                splits += [lastval]
                lastval += np.random.uniform(0., 1.)
        return [el / lastval for el in splits[:-1]]        
                

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
                

#inequality constraint from line
def getLineFromSegment(line):
        a = line[0]; b = line[1]; c = a.copy() ; c[2] = 1.
        normal = cross((b-a),(c-a))
        normal /= norm(normal)
        # get inequality
        dire = b - a
        coeff = normal
        rhs = a.dot(normal)
        return (coeff, array([rhs]))
   
def genConstraintsPerPhase(pDef, numphases):        
        pData = setupControlPoints(pDef)
        vB = varBezier()
        bezVar = vB.fromBezier(pData.bezier())
        line_current = [array([1.,1.,0.]), array([0.,1.,0.]),array([0.,1.,0.])]
        dimExtra = 0
        for i in range(numphases):
                init_points = []
                if i == 0:
                        init_points = [bezVar.waypoints()[1][:3,0][:2]]
                if i == numphases-1:
                        init_points = [bezVar.waypoints()[1][-3:,-1][:2]]
                lines, ptList = genFromLine(line_current, 5, [[0,5],[0,5]],init_points)
                matineq0 = None; vecineq0 = None
                for line in lines:
                        (mat,vec) = getLineFromSegment(line)
                        (matineq0, vecineq0) = (concat(matineq0,mat), concatvec(vecineq0,vec))
                ineq  = (matineq0, vecineq0)
                pDef.addInequality(matineq0,vecineq0.reshape((-1,1)))
                line_current = getRightMostLine(ptList)
                plotPoly  (lines, colors[i])
######################## generate problems ########################


######################## save problems ########################

def writearray(f, a):
        for i in range(a.shape[0]):
                line = ""
                for j in range(a.shape[1]-1):
                        line += str(a[i][j]) + " "
                line+= str(a[i][-1])+"\n"
                f.write(line)
        f.write("\n")

def saveProblem(pDef):
        f = open("test","w")
        # degree totaltime flag
        f.write(str(pDef.degree)+"\n")
        f.write(str(pDef.totalTime)+"\n")
        f.write(str(int(pDef.flag))+"\n")
        f.write("\n")
        writearray(f, pDef.start)
        writearray(f, pDef.end)
        writearray(f, pDef.splits)
        i = 0
        while(True):
                try:
                        ineq = pDef.inequality(i)
                        writearray(f, ineq.A)
                        writearray(f, ineq.b)
                        i+=1
                except:
                        f.close()
                        return
        f.write()
        f.close()

######################## solve a given problem ########################


colors=['b', 'g', 'r', 'c', 'm', 'y', 'k']
labels=['Distance', 'Velocity', 'Acceleration', 'Jerk']
colors2=[colors[len(colors)-1-i] for i in range(len(colors))]



def genProblemDef(numvars = 3, numcurves= 4):
        valDep = array([[np.random.uniform(0., 1.), np.random.uniform(0.,5.), 0. ]]).T
        valEnd = array([[np.random.uniform(5., 10.), np.random.uniform(0.,5.), 0.]]).T
        pDef = problemDefinition()
        #~ pDef.flag =  int(constraint_flag.END_POS) | int(constraint_flag.INIT_POS) | int(constraint_flag.END_VEL)
        pDef.flag =  int(constraint_flag.END_POS) |  int(constraint_flag.END_VEL)
        #~ pDef.flag =  constraint_flag.INIT_POS
        #~ pDef.flag =  constraint_flag.END_POS
        pDef.start = valDep
        pDef.end = valEnd
        pDef.degree = numvars - 1
        pDef.splits = array([genSplit(numcurves)]).T        
        genConstraintsPerPhase(pDef, numcurves) #load random inequality constraints
        saveProblem(pDef)
        return pDef

def __getbezVar(pDef): #for plotting only
        vB = varBezier()
        pData = setupControlPoints(pDef)
        return vB.fromBezier(pData.bezier())

def __addZeroConstants(res, pDef):
        r = res.tolist()
        if (int)(constraint_flag.INIT_POS) & (int)(pDef.flag):
                r = [0.,0.,0.] + r                
                if (int)(constraint_flag.INIT_VEL) & (int)(pDef.flag):
                        r = [0.,0.,0.] + r
        if (int)(constraint_flag.END_POS) & (int)(pDef.flag):
                r = r + [0.,0.,0.]
                if (int)(constraint_flag.END_VEL) & (int)(pDef.flag):
                        r = r + [0.,0.,0.]
        return array(r)
        
import time
import math
        
def computeTrajectory(pDef, save, filename = uuid.uuid4().hex.upper()[0:6]):
        global idxFile
        global colors
        bezVar = __getbezVar(pDef)
        subs = bezVar.split(pDef.splits.reshape((-1)).tolist())
        a = time.clock()
        ineq = generate_problem(pDef);
        b = time.clock()
        print "clock ", b-a
        
        #qp vars
        dimVar = ineq.cost.A.shape[0]
        P = ineq.cost.A * 2.
        q = ineq.cost.b.flatten()
        G = zeros([2,dimVar])
        h = zeros(2)
        G[0,-1] =  1 ; h[0]=1.
        G[1,-1] = -1 ; h[1]=0.
        #~ G = vstack([G,ineq.A])
        G = ineq.A
        #~ h = concatvec(h,ineq.b.reshape((-1)))
        h = ineq.b.reshape((-1))
        
        dimExtra = 0
        
        C = None; d = None
        #~ try:
        if True:                        
                res = quadprog_solve_qp(P, q, G=G, h=h, C=C, d=d)
                res = __addZeroConstants(res, pDef)
                #plot bezier
                derivative = (int)(pDef.costFlag)
                color = colors[derivative]
                for i, bez in enumerate(subs):
                        #~ color = colors[i]
                        test = bez.toBezier3(res[:])
                        if( i == 0):
                                plotBezier(test, color, label = labels[derivative], linewidth =7 - derivative*2)
                        else:
                                plotBezier(test, color, linewidth =7 - derivative*2)
                if save:
                        plt.savefig(filename+str(idxFile))
                #plot subs control points
                for i, bez in enumerate(subs):
                        #~ color = colors[i]
                        test = bez.toBezier3(res[:])
                        #~ plotBezier(test, color)
                        #~ plotControlPoints(test, color)
                if save:
                        plt.savefig("subcp"+filename+str(idxFile))
                final = bezVar.toBezier3(res[:])
                #~ plotControlPoints(final, "black")
                if save:
                        plt.savefig("cp"+filename+str(idxFile))
                idxFile += 1
                if save:
                        plt.close()
                #~ else:
                        #~ plt.show()
                
                return final
        #~ except ValueError:
                #~ plt.close()
                return P, q
######################## solve a given problem ########################


#solve and gen problem
def gen(save = False):
        pDef = genProblemDef(5,3)
        pDef.costFlag = derivative_flag.DISTANCE
        computeTrajectory(pDef, save)
        pDef.costFlag = derivative_flag.VELOCITY
        computeTrajectory(pDef, save)
        pDef.costFlag = derivative_flag.ACCELERATION
        computeTrajectory(pDef, save)
        pDef.costFlag = derivative_flag.JERK
        res = computeTrajectory(pDef, save)
        plt.legend(loc='upper left')
        plt.show()
        return res

res = None
for i in range(1):
        res = gen(False)
        #~ if res[0] != None:
                #~ break

np.set_printoptions(precision=3, suppress=True, threshold=np.nan)
