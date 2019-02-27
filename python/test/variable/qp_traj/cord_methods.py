from hpp_spline import *
from numpy import matrix, array, zeros, ones, diag, cross, dot
from numpy.linalg import norm

__EPS = 1e-6

import numpy as np
np.set_printoptions(precision=3, suppress=True, threshold=np.nan)
from varBezier import varBezier
from qp_cord import *
from plot_cord import *


idxFile = 0
import string
import uuid; uuid.uuid4().hex.upper()[0:6]


######################## generate problems ########################

def genSplit(numCurves, minTimes = None, maxTimes = None):
        if (numCurves == 0):
                return [-1]
        splits = []
        lastval = np.random.uniform(0., 1.)
        if minTimes is None:
                minTimes = [0. for i in range(numCurves)]
        if maxTimes is None:
                maxTimes = [1. for i in range(numCurves)]
        for i in range(numCurves):
                splits += [lastval]
                #~ lastval += np.random.uniform(minTimes[i], maxTimes[i])
                lastval += np.random.uniform(0.,1.)
        return [el / lastval for el in splits[:-1]]   
             
def addSplit(times):
        res =[times[0]]
        for el in times[1:]:
                res += [res[-1] + el]    
        return res
                

def getRightMostLine(ptList):
        pt1 = array([0.,0.,0.])
        pt2 = array([0.,0.,0.])
        for pt in ptList:
                if pt[0] > pt1[0]:
                       pt1 =  pt[:]
                elif pt[0] > pt2[0]:
                       pt2 =  pt[:]
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
        ineqs = []
        for i in range(numphases):
                init_points = []
                if i == 0:
                        #~ init_points = [bezVar.waypoints()[1][:3,0][:2]]
                        init_points = [pDef.start.flatten()[:2]]; 
                        init_points = init_points + [init_points[-1] + np.array([0.1,0.1])]
                        init_points = init_points + [init_points[-1] + np.array([-0.1,0.1])]
                        init_points = init_points + [init_points[-1] + np.array([-0.1,-0.1])]
                        init_points = init_points + [init_points[-1] + np.array([0.1,0.1])]
                if i == numphases-1:
                        init_points = [bezVar.waypoints()[1][-3:,-1][:2]]
                lines, ptList = genFromLine(line_current, 5, [[0,5],[0,5]],init_points)
                matineq0 = None; vecineq0 = None
                for line in lines:
                        (mat,vec) = getLineFromSegment(line)
                        (matineq0, vecineq0) = (concat(matineq0,mat), concatvec(vecineq0,vec))
                ineq  = (matineq0, vecineq0)
                ineqs += [(matineq0[:], vecineq0[:])]
                pDef.addInequality(matineq0,vecineq0.reshape((-1,1)))
                line_current = getRightMostLine(ptList)
                plotPoly  (lines, colors[i])
        return ineqs
   
def genConstraintsPerPhaseNew(pDef, numphases):        
        pData = setupControlPoints(pDef)
        vB = varBezier()
        bezVar = vB.fromBezier(pData.bezier())
        line_current = [array([1.,1.,0.]), array([0.,1.,0.]),array([0.,1.,0.])]
        dimExtra = 0
        ineqs = []
        for i in range(numphases):
                #~ res_tmp = []
                init_points = []
                if i == 0 and int (pDef.flag) & int(constraint_flag.INIT_POS):
                        #~ init_points = [bezVar.waypoints()[1][:3,0][:2]]
                        init_points = [pDef.start.flatten()[:2]]
                if i == numphases-1 :
                        init_points = [pDef.end.flatten()[:2]]
                        #~ init_points = [bezVar.waypoints()[1][-3:,-1][:2]]
                #~ lines, ptList = genFromLine(line_current, 5, [[4*i-1,4*i+3],[0,5]],init_points)
                lines, ptList = genFromLine(line_current, 5, [[0,5],[0,5]],init_points)
                matineq0 = None; vecineq0 = None
                for line in lines:
                        (mat,vec) = getLineFromSegment(line)
                        #~ res_tmp += [(mat,vec)]
                        (matineq0, vecineq0) = (concat(matineq0,mat), concatvec(vecineq0,vec))
                ineq  = (matineq0, vecineq0)
                
                ineqs += [(matineq0[:], vecineq0[:])]
                pDef.addInequality(matineq0,vecineq0.reshape((-1,1)))
                line_current = getRightMostLine(ptList)
                plotPoly  (lines, colors[i])
        return ineqs
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
        #~ valEnd = array([[1.,1.,0.]]).T
        #~ valEnd = array([[2.,2.,2.]]).T
        pDef = problemDefinition()
        #~ pDef.flag =  int(constraint_flag.INIT_POS) | int(constraint_flag.INIT_VEL)
        pDef.flag =  int(constraint_flag.INIT_POS) | int(constraint_flag.INIT_VEL) | int(constraint_flag.END_POS)
        #~ pDef.flag =  int(constraint_flag.INIT_POS) | int(constraint_flag.INIT_VEL) | int(constraint_flag.INIT_ACC)
        #~ pDef.flag =  int(constraint_flag.INIT_POS) | int(constraint_flag.INIT_VEL) | int(constraint_flag.INIT_ACC)
        pDef.start = valDep
        pDef.end = valEnd
        pDef.degree = numvars - 1
        pDef.splits = array([genSplit(numcurves)]).T   
        inequalities_per_phase = genConstraintsPerPhase(pDef, numcurves) #load random inequality constraints
        #~ saveProblem(pDef)
        return pDef,inequalities_per_phase

def __getbezVar(pDef): #for plotting only
        vB = varBezier()
        pData = setupControlPoints(pDef)
        return vB.fromBezier(pData.bezier())

import time
import math


def qpineq(cost, ineq, eq = None, verbose = False):        
        #~ print "ineq ", ineq.A
        dimVar = cost[0].shape[0]
        P = cost[0]* 2.
        q = cost[1].flatten()
        G = ineq[0]
        h = ineq[1].reshape((-1))
        eqDim = 0
        if eq is not None:
                eqDim = eq[1].shape[0]
        C = zeros((dimVar+eqDim,dimVar)); 
        for i in range(2,dimVar,3):
                C[i,i]=1.                
        d = zeros(dimVar+eqDim)  
        if eq is not None:
                C[eqDim,:] = eq[0]
                d [-eqDim] = eq[1]  
                C=eq[0]
                d=eq[1]
                #~ G = None
                #~ h = None
        
        #~ res = quadprog_solve_qp(P, q, G=G, h=h, C=C, d=d, verbose = verbose)
        return [quadprog_solve_qp(P, q, G=G, h=h, C=C, d=d, verbose = verbose), P, q, G, h]
        #~ return [quadprog_solve_qp(P, q, G=None, h=None, C=eq[0], d=eq[1]), P, q, G, h]


def evalBez(res, pDef):
        bezVar = __getbezVar(pDef)
        return bezVar.toBezier3(res[:]);


gb = -1

def plot(pDef, res, filename, saveToFile, resIsBezier = False):        
        global idxFile
        global colors
        global gb
        bezVar = __getbezVar(pDef)
        subs = bezVar.split(pDef.splits.reshape((-1)).tolist())
        final = None
        if resIsBezier:
                final =res
        else:
                final = evalBez(res,pDef)
        derivative = (int)(pDef.costFlag)
        color = colors[(derivative + gb) % len(colors)]
        #~ gb = gb + 1
        if resIsBezier:
                plotBezier(final, "r", linewidth =7 - derivative*2)
                plotControlPoints(final, "b")
        else:
                for i, bez in enumerate(subs):
                        #~ color = colors[i]
                        test = bez.toBezier3(res[:])
                        if( i == 0):
                                plotBezier(test, color, label = labels[derivative], linewidth =7 - derivative*2)
                        else:
                                plotBezier(test, color, linewidth =7 - derivative*2)
                if saveToFile:
                        plt.savefig(filename+str(idxFile))
                #plot subs control points
                for i, bez in enumerate(subs):
                        #~ color = colors[i]
                        test = bez.toBezier3(res[:])
                        #~ plotBezier(test, color)
                        #~ plotControlPoints(test, color)
                if saveToFile:
                        plt.savefig("subcp"+filename+str(idxFile))
                #~ plotControlPoints(final, "black")
                if saveToFile:
                        plt.savefig("cp"+filename+str(idxFile))
                idxFile += 1
                if saveToFile:
                        plt.close()
        #~ else:
                #~ plt.show()
        
        return final
        
def oneVelIneq(pDef, ineq):
        bvel = 10*ones(3)
        
        bezVar, b = __getbezVar(pDef)
        #derivates
        db = b.compute_derivate(1)
        (matineq0, vecineq0) = ineq
        for i in range(db.nbWaypoints):
                A,b = db.waypoints().A[:, i*3:i*3+3].transpose() * db.mult_T , db.waypoints().b[:,i] * db.mult_T
                (matineq0, vecineq0) = (concat(matineq0, A), concatvec(vecineq0,bvel - b))    
                (matineq0, vecineq0) = (concat(matineq0,-A), concatvec(vecineq0,bvel + b))    
        return (matineq0, vecineq0)
        

        
def computeTrajectory(pDef, saveToFile, filename = uuid.uuid4().hex.upper()[0:6]):        
        ineq = generate_problem(pDef);
        #~ ineq = oneVelIneq(ineq)
        
        #~ print "INEQ " , ineq.A
        (res,P, q, G, H) = (None, None, None, None,None)
        #~ try:                  
        if True:                  
                (res,P, q, G, H) = qpineq((ineq.cost.A, ineq.cost.b),(ineq.A, ineq.b), verbose = True)
                #~ print "cost ", res[1]
                #~ plot(pDef, res[0], filename, saveToFile)
                return res[0], res[1]
        #~ except ValueError:
                #~ print "FAIl traj"
                #~ raise ValueError
######################## solve a given problem ########################
import math

def ineq2Phases(phase0, phase1):
        return (concat(phase0[0],phase1[0]), concatvec(phase0[1],phase1[1]))
        
def def_find_min_t_out(degree, phase0, phase1, bezierPrev, i  = None, relax_vel = False, last = False):
        pDef = problemDefinition()
        pDef.flag =  int(constraint_flag.INIT_POS) | int(constraint_flag.INIT_VEL)
        if relax_vel:
                print "relax ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                pDef.flag =  int(constraint_flag.INIT_POS)
        if last:
                pDef.flag =  pDef.flag | int(constraint_flag.END_POS)
        pDef.costFlag = derivative_flag.VELOCITY
        pDef.splits = array([genSplit(1)]).T   
        c = curve_constraints()
        T_prev_phase = bezierPrev.max()
        init_vel_calc = bezierPrev.degree * T_prev_phase * (bezierPrev.waypoints()[:,-1].flatten() - bezierPrev.waypoints()[:,-2].flatten())
        #~ print "init_vel_calc", init_vel_calc
        dprev = bezierPrev.compute_derivate(1)
        dprev2 = dprev.compute_derivate(1)
        #~ try:
        c.init_vel =  init_vel_calc.reshape((3,1))
        #~ c.init_acc =  dprev2 (bezierPrev.max())
        pDef.curveConstraints = c
        
        #~ print "bezierPrev max", bezierPrev.max()
        #~ print "init_vel constraint", c.init_vel
        pDef.start = bezierPrev.waypoints()[:,-1].reshape((3,1))
        pDef.degree = degree
        
        ph0 = (phase0[0][:],phase0[1][:])
        
        pDef.addInequality(ph0[0],ph0[1].reshape((-1,1)))
        ineq = generate_problem(pDef)
        
        if math.isnan(ineq.cost.b[0][0]):
                print "pb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ", ineq.cost.b
                saveProblem(pDef)
                ineq = generate_problem(pDef)
                if math.isnan(ineq.cost.b[0][0]):
                        print "pbkqjlkjqs", ineq.cost.b
                        raise Exception('I know Python!') 
                        
        (matineq0, vecineq0) = (ineq.A, ineq.b.reshape((-1)))        
        
        #add constraint on last waypoint
        lastinphase1 = zeros((phase1[1].shape[0],ineq.A.shape[1]))
        lastinphase1[:,-3:] = phase1[0]     
        (matineq0, vecineq0) = (concat(ineq.A,lastinphase1), concatvec(ineq.b.reshape((-1)),phase1[1]))           
        #add constraint on velocity waypoint
        # x+1 = 2 x_end - x-1
        # with P1 and p1 phase 1 matrix / vector
        # P1 (2 x_end - x-1) <= p1
        
        
        #add constraint on end velocity waypoint
        # x+1 = x_end + velConstraint(i) / n / T_next with T_next = 1
        # velConstraint(i) = (x_end - x-1) * n * T_prev
        # x+1 = x_end + (x_end - x-1) * T_prev 
        # x+1 = x_end(1 + T_prev) + (- x-1) * T_prev 
        # P1 (x_end(1 + T_prev) - x-1) <= p1
        #  T_prev * (P1 (x_end - x-1) ) + P1 * x_end <= p1 
        #  T_prev * (P1 (x_end - x-1) ) <= p1 - P1 * x_end-1  TODO
        #~ velConstraint = 
        #~ m1T_prev = bezierPrev.max()
        Tm1 =1.; Tp1 = 1.
        veclcon = zeros((phase1[1].shape[0],ineq.A.shape[1]))  
        veclcon[:,-6:-3] = -phase1[0] * Tm1 / Tp1
        veclcon[:,-3:] = (1+Tp1)*phase1[0]        
        (matineq0, vecineq0) = (concat(matineq0,veclcon), concatvec(vecineq0,phase1[1])) 
        #add constraint on acceleration waypoint
        # with P1 and p1 phase 1 matrix / vector
        # P1 (x-2) <= p1            
        #~ acccon = zeros((phase1[1].shape[0],ineq.A.shape[1]))    
        #~ acccon[:,-9:-6] = phase1[0]              
        #~ (matineq0, vecineq0) = (concat(matineq0,veclcon), concatvec(vecineq0,phase1[1]))
        #~ (matineq0, vecineq0) = (concat(matineq0,acccon), concatvec(vecineq0,phase1[1]))
        
        (costA,costB) = (ineq.cost.A,ineq.cost.b)    
                
        #~ if relax_vel: #TODO
                #~ Avel = zeros((costA.shape[0],ineq.A.shape[1])) 
                #compute cost least square
                #~ bve = 2* pDef.start - c.init_vel
                #~ veclcon = zeros((ineq.cost.A.shape))    
                #~ print "shapes ", costB.shape
                #~ costA[3:6,3:6] = costA[3:6,3:6] + 1000 * identity(3)   
                #~ costB[3:6] = (1000. * 2)*(bve.reshape((-1,1))) 
                #~ print "relax"
                     
                #~ pDef.flag =  int(constraint_flag.INIT_POS)
        return pDef, (costA,costB),(matineq0, vecineq0) 
        
def zeroBez():
        waypoints = array([[0.,0.,0.], [0.,0.,0.]]).transpose()
        return bezier(waypoints, 1.)
        
def Bez(pos):
        waypoints = array([pos, pos]).transpose()
        return bezier(waypoints, 1.)
        
        
def find_minimum_time_constrained_vel(bezierPrev, bez, phaseA, phaseB, Min = True):
        
        #~ bretime = bezier(bez.waypoints())
        #~ return 1., bretime
        
        (matineq0, vecineq0) = (ones(1)*(-1),ones(1)*(-0.01)) # t > 0
        #~ (matineq0, vecineq0) = (zeros(1),zeros(1)) 
        (matineq0, vecineq0) = (concatvec(matineq0, ones(1)), concatvec(vecineq0,ones(1)*2)) 
        bvel = 10*ones(3)
        #~ velConstraint = bezierPrev.compute_derivate(1)(bezierPrev.max()).flatten()
        #~ velConstraint  = bez.derivate(bez.min(),1).flatten()
        T_prev_phase = bezierPrev.max()
        Tm1 = bez.max()
        #~ print "T_prev_phase ", T_prev_phase
        #~ print "tm1 ", Tm1
        velConstraint  = bezierPrev.degree * T_prev_phase * (bezierPrev.waypoints()[:,-1].flatten() - bezierPrev.waypoints()[:,-2].flatten())
        velConstraint2 = bez.degree * Tm1 * (bez.waypoints()[:,1].flatten() - bez.waypoints()[:,0].flatten())
        
        #~ print "bezierPrev start ", bez.waypoints()[:,0].flatten()
        #~ print "bez start ", bez.waypoints()[:,0].flatten()
        
        #~ print "velConstraint ", velConstraint
        #~ print "velConstraint2 ", velConstraint2
        b = bez.compute_derivate(1)
        cost = ones(1)
        if not Min:
                cost = ones(1)*(-1)        
        #t scaling constraints
        # x(1)  = x(0) + v_final(i-1) /n / t_out(i)
        # A(x(0) + v_final(i-1) /n / t_out(i)) <= a
        # A(x(0) + v_final(i-1) /n / t_out(i)) <= a
        # A x(0) * t_out(i) + A v_final(i-1) /n <= a * t_out(i)
        # (A x(0) - a) * t_out(i) <= - A v_final(i-1) /n 
        wps = bez.waypoints()
        #~ print "init wp ", bez.waypoints()
        # make sure end vel and acceleration(t) is in the next polytope
        # TODO
        
        deg = bez.degree
        #add constraint on velocity waypoint
        # x+1 = x_0 + velConstraint(i-1) / n / T
        # P1 * x_0 * T + P1 / n * velConstraint <= p 1 * T
        # (P1 * x_0 -  p 1) * T  <= - P1 / n * velConstraint
        veclcon_t =  phaseA[0].dot(wps[:,0].flatten())  - phaseA[1] 
        veclcon_c = -phaseA[0].dot(velConstraint) / deg
        #~ print "constrants t ", veclcon_t
        #~ print "veclcon_c t ", veclcon_c
        #~ print "deg ", deg
        #~ print "velConstraint ", velConstraint
        #~ print "velConstraint now ",  bez.derivate(bez.min(),1).flatten()
        #~ print "wp0 ", wps[:,0].flatten()
        #~ print "bezprev wps ",bezierPrev.compute_derivate(1).waypoints()
        (matineq0, vecineq0) = (concatvec(matineq0, veclcon_t), concatvec(vecineq0,veclcon_c)) 
        #add constraint on acceleration waypoint
        #TODO
        
        #add constraint on end velocity waypoint
        # x+1 = x_end + velConstraint(i) / n / T_next with T_next = 1
        # velConstraint(i) = (x_end - x-1) * n * T_prev
        # x+1 = x_end + (x_end - x-1) * T_prev 
        # P2 (x_end + (x_end - x-1) * T_prev ) <= p2
        #  T * (P2 (x_end - x-1) ) + P2 * x_end <= p2 
        #  T * (P2 (x_end - x-1) ) <= p2 - P2 * x_end-1
        #~ print "AA ", deg
        veclcon_t = phaseB[0].dot((wps[:,-1] - wps[:,-2]).flatten()) 
        #~ # comment
        #~ print "A A", deg
        veclcon_c = phaseB[1] - phaseB[0].dot(wps[:,-2].flatten()) 
        #~ (matineq0, vecineq0) = (concatvec(matineq0, veclcon_t), concatvec(vecineq0,veclcon_c)) 
        #~ print "A ", deg
        
        #polytope velocity constraints
        wps = b.waypoints()
        
        #~ print "B ", deg
        for i in range(b.nbWaypoints):
                if i != 1:
                        A = wps[:,i].flatten()
                        (matineq0, vecineq0) = (concatvec(matineq0, -bvel), concatvec(vecineq0,A))   
                # TODO ADD CONSTRAINT FOR 1
        try:                
                #~ print "C ", deg
                vecineq0 = vecineq0 + ones(vecineq0.shape[0])*0.0001
                t = quadprog_solve_qp(ones((1,1))*0.001,cost,matineq0.reshape((-1,1)),vecineq0).tolist()[0]
                #~ print "D", deg
        except:
                print "dafuq !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", (matineq0.dot(1.)- vecineq0 < 0.001 ).all()
                #~ print "prev max ", bezierPrev.max()
                #~ print "Unexpected error:", sys.exc_info()[0]
        #rewrite velocity wp
        wps = bez.waypoints()
        #x(1)  = x(0) + v_final(i-1) * t_out(i)
        # x(1)  = x(0) + v_final(i-1) /n / t_out(i)
        #~ wps[:,1] = wps[:,-1] + (wps[:,-1] - wps[:,-2])
        
        #~ print "E", deg
        wps[:,1] = wps[:,0] + (velConstraint) / deg * t
        #~ print "F", deg
        #~ print "t ", t
        #~ bretime = bezier(wps, t)
        bretime = bezier(wps,t)
        
        #~ print "G", deg
        #~ print "end wp ", bretime.waypoints()
        #~ if abs(t - 1.)>0.1:
        print "t ", t
                #~ t = 1.
        return t, bretime
        
        
def solveForPhase(degree, bezierPrev, phaseA, phaseB, filename="", saveToFile=False, Min = False, relax_vel = False, last = False):  
        npDef, cost, ineq = def_find_min_t_out(degree, phaseA,phaseB, bezierPrev, relax_vel=relax_vel, last=last)
        (res,P, q, G, H) = (None, None, None, None,None)
        # now try for each constraint
        best_res = None
        best_cost = -10000;
        best_time = None;
        bretime = None
        foundOne = False
        if Min:
                best_cost = -best_cost
        
        for i in range(phaseA[0].shape[0]):
                try:                  
                #~ if True:              
                        (lineIneq,lineinec)= (phaseA[0][i,:], phaseA[1][i])
                        lastinphase1 = zeros((1,ineq[0].shape[1]))
                        lastinphase1[:,-3:] = lineIneq
                        eq = None
                        if (not last):
                                eq = (lastinphase1, array([lineinec]))
                        (res,P, q, G, H) = qpineq(cost, ineq, eq,verbose = True)
                        foundOne = True
                        b = evalBez(res[0],npDef)
                        c = approxLengthBez(b)
                        #~ c = res[1]
                        c, br = find_minimum_time_constrained_vel(bezierPrev, b, phaseA, phaseB, Min)
                        #~ if (Min and c < best_cost) or ((not Min) and c > best_cost):
                        if (c < abs(best_cost)):
                                best_res = res[0]
                                best_cost = c
                                best_time = c
                                #~ bretime = b
                                bretime = br
                                #~ print "t", best_time
                except ValueError, e:
                        #~ print e
                        #~ if relax_vel:
                                #~ print "not saved"
                        #~ else:
                                #~ "print wtd"
                        #~ vel =  npDef.start.flatten() + (npDef.curveConstraints.init_vel.flatten()) / float(npDef.degree)
                        # plt.scatter([npDef.start.flatten()[0]],[npDef.start.flatten()[1]],color="r")      
                        pass
        if best_res is not None: 
                if best_time < 0.3:
                        plot(npDef, best_res, filename, saveToFile)
                        plot(npDef, bretime, filename, saveToFile, True)
                        plt.show()
                tg = 0
        else:
                print "no min time found"
                print foundOne
                raise ValueError  
        #~ return evalBez(best_res,npDef), 1., evalBez(best_res,npDef)
        return evalBez(best_res,npDef), best_time, bretime
        
def approxLengthBez(b):
        res = 0.
        mx = b.max()
        for i in range(100):
                t1 = float(i)/100.
                t2 = float(i+1)/100.
                res = res + norm(b(t2).flatten() - b(t1).flatten())
        return res
        
#~ def solveForPhases(pDef, inequalities_per_phase, filename="", saveToFile=False, Min = True): 
        #~ bezierPrev = Bez(pDef.start.reshape((-1)))
        #~ timesMin = [];
        #~ for i in range(len(inequalities_per_phase)-1):
                #~ A = inequalities_per_phase[i]
                #~ bezierPrev = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[i], inequalities_per_phase[i+1], filename=filename, saveToFile=saveToFile, Min = Min)
                #~ timesMin += [norm(bezierPrev(bezierPrev.max()).flatten() - bezierPrev(bezierPrev.min()).flatten())]   
        #~ if Min:
                #~ timesMin += [0.]
        #~ else:
                #~ bezierPrev = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[-1], inequalities_per_phase[-1], filename=filename, saveToFile=saveToFile, Min = False)
                # timesMin += [norm(bezierPrev(bezierPrev.max()).flatten() - bezierPrev(bezierPrev.min()).flatten())]
                #~ timesMin += [approxLengthBez(bezierPrev)]
        #~ return timesMin
        
def solveForPhases(pDef, inequalities_per_phase, filename="", saveToFile=False, Min = True): 
        bezierPrev = Bez(pDef.start.reshape((-1)))
        timesMin = [];
        for i in range(len(inequalities_per_phase)-1):
                print "phase ", i
                A = inequalities_per_phase[i]
                bezierPrev, time, bezretime = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[i], inequalities_per_phase[i+1], filename=filename, saveToFile=saveToFile, Min = Min)
                timesMin += [time]   
                #~ bezierPrev = bezretime
        #~ if Min:
                #~ timesMin += array([0.01])
        #~ else:
        #~ bezierPrev, time = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[-1], inequalities_per_phase[-1], filename=filename, saveToFile=saveToFile, Min = Min, last = True)
        #~ timesMin += [norm(bezierPrev(bezierPrev.max()).flatten() - bezierPrev(bezierPrev.min()).flatten())]
        #~ timesMin += [approxLengthBez(bezierPrev)]
        timesMin += [time]
        return array(timesMin).reshape((-1)).tolist()
        
def findTimesToSplit(pDef, inequalities_per_phase, filename="", saveToFile=False):
        times2 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=False)
        times1 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=True)
        avg = times1[:]
        for i in range(len(times1)):
                if(times1[i] > times2[i]):
                        raise ValueError
        avg = [(times1[i] + times2[i]) / 2. for i in range(len(times1))]
        total = sum(avg)
        avg = addSplit([el / total for el in avg])
        #~ print "avg", avg
        return avg[:-1]

if __name__ == '__main__':
                
                
        totalScenarios = 0
        totalMinDist = 0
        totalHeuristicTimes = 0
        totalHeuristicTimesHighD = 0
        totalRandomTimes = 0
        totalRandomTimesHighD = 0
        avgcost = 0.
        avgcostRand = 0.
        avgcostHighD = 0.
        avgcostRandHighD = 0.
        #solve and gen problem
        def gen(saveToFile = False):
                global totalScenarios
                global totalMinDist
                global totalScenarios
                global totalHeuristicTimes
                global totalRandomTimes
                global totalHeuristicTimesHighD
                global totalRandomTimesHighD
                global avgcost
                global avgcostRand
                global avgcostHighD
                global avgcostRandHighD
                plt.close()
                while(True):
                        numcurves = 4
                        pDef, inequalities_per_phase = genProblemDef(5,numcurves)
                        pDef.costFlag = derivative_flag.VELOCITY
                        #~ times = findTimesToSplit(pDef, inequalities_per_phase)
                        try:
                                times = findTimesToSplit(pDef, inequalities_per_phase)
                                totalMinDist = totalMinDist + 1
                                totalScenarios = totalScenarios + 1
                                plt.show()
                                return
                                #~ plt.show()
                        
                        except ValueError:
                                pass      
                        except:
                                return
                        for i in range(10):
                                try:
                                        pDef.splits = array([genSplit(numcurves)]).T 
                                        res, cos = computeTrajectory(pDef, saveToFile)
                                        b = evalBez(res,pDef)
                                        plot(pDef, res, "", False)
                                        plt.show()
                                        totalScenarios = totalScenarios + 1
                                        return 
                                except ValueError:
                                        pass
                        return
                                #~ plt.close()

        (P, q, G,h, res) = (None,None,None,None, None)
        for i in range(100):
                #~ (P, q, G,h, res) = gen(False)
                gen(False)
                gb = -1
                #~ if res[0] != None:
                        #~ break
        print "total scenarios", totalScenarios 
        print "total min time success", totalMinDist
        #~ print "total HeuristicTimes", totalHeuristicTimes, ' cost ', avgcost / totalHeuristicTimes
        #~ print "total totalRandomTimes", totalRandomTimes, ' cost ', avgcostRand / totalRandomTimes
        #~ print "total totalHeuristicTimesHighD", totalHeuristicTimesHighD, ' cost ', avgcostHighD / totalHeuristicTimesHighD
        #~ print "total totalRandomtotalRandomTimesHighDTimes", totalRandomTimesHighD, ' cost ', avgcostRandHighD / totalRandomTimesHighD

        def cost(P, q, G,h, x):
                print (x.T.dot(P).dot(x) / 2. + q.dot(x))
                print "ineq ?",  (G.dot(x) -h <=0.0001).all()

        #~ zero = array([2.,2.,0.])
        #~ print "res, " ; cost(P,q, G,h, res)
        #~ print "zero, " ; cost(P,q, G,h, zero)
        #~ print "-zero, " ; cost(P,q, G,h, -zero)

        #~ plt.show()

