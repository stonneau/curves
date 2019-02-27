from  cord_methods import *

        
problem_gen_times = 0.
qp_times = 0.
num_qp_times = 0
######################## solve a given problem ########################

        
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
        
def findTimesToSplit(pDef, inequalities_per_phase, filename="", saveToFile=False):
        times2 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=False)
        times1 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=True)
        #~ print "avg", avg
        return [0. for el in times1], times1

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
def gen(saveToFile = False, degree = 4, numcurves= 10):
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
        #~ plt.close()
        while(True):
                #~ numcurves = 3
                pDef, inequalities_per_phase = genProblemDef(degree,numcurves)
                pDef.costFlag = derivative_flag.VELOCITY
                totalScenarios = totalScenarios + 1  
                #~ timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                try:
                #~ if True:
                        print 'a'
                        timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                        print 'b'
                        #~ if not relaxed:
                                #~ break
                        #~ if relaxed:
                                #~ break
                                #~ print "relaxed"
                                #~ totalScenarios = totalScenarios + 1  
                                #~ print "times", timesMax
                        totalMinDist = totalMinDist + 1
                        for i in range(10):
                                #~ if True:
                                try:
                                        print "heho"
                                        pDef.degree = 10
                                        pDef.splits = array([genSplit(numcurves,timesMin,timesMax)]).T 
                                        #~ timesMid = [(timesMin[i]+ timesMax[i])/2. for i in range(len(timesMin))]
                                        #~ pDef.splits = array([genSplit(numcurves,timesMid,timesMid)]).T 
                                        res1, cos1 = computeTrajectory(pDef, saveToFile)
                                        #~ pDef.splits = array([genSplit(numcurves,timesMin,timesMax)]).T 
                                        #~ cos2, res2 = computeTrajectory(pDef, saveToFile)
                                        b1 = evalBez(res1, pDef)
                                        #~ b2 = evalBez(res2, pDef)
                                        #~ plt.close()
                                        #~ plotBezier(b1, "r", label = None, linewidth = 2.0)
                                        #~ plotBezier(b2, "g", label = None, linewidth = 2.0)
                                        #~ plt.show()  
                                        totalRandomTimes = totalRandomTimes + 1 
                                        return
                                except ValueError, e:
                                        print e
                                        pass
                                        #~ plt.close()     
                        return 0
                except ValueError:
                        #~ print "total fail"
                        pass

(P, q, G,h, res) = (None,None,None,None, None)
totaltrials = 10


benchs = []
degree = []
numphases = []

#~ def stat(degree, numphases):
        #~ [gen(False, degree, numphases)]

for i in range(totaltrials):
        #~ (P, q, G,h, res) = gen(False)
        gen(False)
        gb = -1
        #~ if res[0] != None:
                #~ break
print "totalScenarios", totalScenarios 
print "total success", totalMinDist 
print "total average  random time to success", float (totalRandomTimes) / float(totalMinDist)
#~ print "avg time to generate problem", float (qp_times) / float(num_qp_times)
#~ print "avg time for successful qp", float (problem_gen_times) / float(num_qp_times)
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

