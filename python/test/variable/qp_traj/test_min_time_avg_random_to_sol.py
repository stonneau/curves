from  cord_methods import *

        
problem_gen_times = 0.
qp_times = 0.
num_qp_times = 0
        
def computeTrajectory(pDef, saveToFile, filename = uuid.uuid4().hex.upper()[0:6]):
        global problem_gen_times
        global qp_times
        global num_qp_times    
        a = time.clock()
        ineq = generate_problem(pDef);
        b = time.clock()
        
        #~ print "INEQ " , ineq.A
        (res,P, q, G, H) = (None, None, None, None,None)
        #~ try:                  
        if True:                               
                c = time.clock()   
                for i in range(1):
                        (res,P, q, G, H) = qpineq((ineq.cost.A, ineq.cost.b),(ineq.A, ineq.b), verbose = True)
                d = time.clock()   
                num_qp_times = num_qp_times + 1
                problem_gen_times = problem_gen_times + (b-a)
                qp_times = qp_times + (d-c)
                #~ print "cost ", res[1]
                plot(pDef, res[0], filename, saveToFile)
                return res[1], res[0]
        #~ except ValueError:
                #~ print "FAIl traj"
                #~ raise ValueError
######################## solve a given problem ########################

def findTimesToSplit(pDef, inequalities_per_phase, filename="", saveToFile=False):
        #~ times2 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=False)
        times1 = solveForPhases(pDef, inequalities_per_phase, filename, saveToFile, Min=True)
        #~ print "avg", avg
        return [0 for _ in times1], times1

def solveForPhases(pDef, inequalities_per_phase, filename="", saveToFile=False, Min = True): 
        bezierPrev = Bez(pDef.start.reshape((-1)))
        timesMin = [];
        for i in range(len(inequalities_per_phase)-1):
                A = inequalities_per_phase[i]
                #~ try:
                bezierPrev, time = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[i], inequalities_per_phase[i+1], filename=filename, saveToFile=saveToFile, Min = Min, relax_vel = False)
                timesMin += [time]   
                        #~ print "ok"
                #~ except:
                        #~ print "relaxvel"
                        #~ bezierPrev, time = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[i], inequalities_per_phase[i+1], filename=filename, saveToFile=
                        #~ saveToFile, Min = Min, relax_vel = True)
                        #~ print "saved"
                        #~ timesMin += [time]
        plt.show()
        if Min:
                timesMin += array([0.01])
        else:
                bezierPrev, time = solveForPhase(pDef.degree, bezierPrev, inequalities_per_phase[-1], inequalities_per_phase[-1], filename=filename, saveToFile=saveToFile, Min = False)
                #~ timesMin += [norm(bezierPrev(bezierPrev.max()).flatten() - bezierPrev(bezierPrev.min()).flatten())]
                #~ timesMin += [approxLengthBez(bezierPrev)]
                timesMin += [time]
        return array(timesMin).reshape((-1)).tolist()

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
def gen(saveToFile = False, degree = 5, numcurves= 3):
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
                #~ numcurves = 3
                pDef, inequalities_per_phase = genProblemDef(degree,numcurves)
                pDef.costFlag = derivative_flag.VELOCITY
                totalScenarios = totalScenarios + 1    
                timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                try:
                #~ if True:
                        timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                        totalMinDist = totalMinDist + 1
                        for i in range(1):
                                #~ if True:
                                try:
                                        pDef.splits = array([genSplit(numcurves,timesMin,timesMax)]).T 
                                        cos1, res1 = computeTrajectory(pDef, saveToFile)
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
                                except ValueError:
                                        pass
                                        plt.close()                        
                except ValueError:
                        #~ print "total fail"
                        pass

(P, q, G,h, res) = (None,None,None,None, None)
totaltrials = 100


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
print "avg time to generate problem", float (qp_times) / float(num_qp_times)
print "avg time for successful qp", float (problem_gen_times) / float(num_qp_times)
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

