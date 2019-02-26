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
                try:
                #~ if True:
                        timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                        print "times max 1", timesMax
                        totalMinDist = totalMinDist + 1
                        pDef.degree +=3
                        timesMin,timesMax = findTimesToSplit(pDef, inequalities_per_phase)
                        print "times max 2", timesMax
                        totalMinDist = totalMinDist + 1
                        return
                except ValueError:
                        #~ print "total fail"
                        pass
                except:
                        totalScenarios = totalScenarios -1

(P, q, G,h, res) = (None,None,None,None, None)

from cPickle import dump, load

        
totaltrials = 100    
benchs = []
#~ degrees = [5,7,9,12]
degrees = [5, 9]
numphases = [2, 4, 6]


def savebench():
        res = [benchs, degrees, numphases]
        fname = "success_rate_min_dist"
        f = open(fname, "w")
        dump(res,f)
        f.close()
        
def loadbench():
        global benchs
        global degrees
        global numphases
        fname = "success_rate_min_dist"
        f = open(fname, "r")
        res = load(f)
        f.close()
        benchs = res[0]
        degrees = res[1]
        numphases = res[2]

def gen_bench():
        global benchs
        def stat(degree, numphases):
                print "numphases ", numphases
                [gen(False, degree, numphases) for _ in range(totaltrials)]
                res = (float)( totalMinDist) / (float) (totalScenarios)       
                global totalScenarios
                global totalMinDist
                totalScenarios = 0
                totalMinDist = 0
                return res
                
        def stats(degree):
                return [[numphase, stat(degree, numphase)] for numphase in numphases]

        for i, deg in enumerate(degrees):
                res = array(stats(deg))
                benchs += [res]
                #~ savebench()
       
def plotbench():
        plt.close()
        global benchs
        for i, ben in enumerate(benchs):
                ben[:,0];
                ben[:,1];
                colors[i];
                degrees[i]
                plt.plot(ben[:,0],ben[:,1],color=colors[i],label="degree " + str(degrees[i]),linewidth=1 )  
                plt.scatter(ben[:,0],ben[:,1],color=colors[i],label="degree " + str(degrees[i]),linewidth=1 )  
                plt.legend(loc='upper left')
        
                
#~ gen_bench(); savebench();
gen_bench()
#~ loadbench()

plotbench()    
plt.show()
