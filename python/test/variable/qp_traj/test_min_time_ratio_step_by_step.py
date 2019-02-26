from  cord_methods import *

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
                        pDef, inequalities_per_phase = genProblemDef(7,numcurves)
                        pDef.costFlag = derivative_flag.VELOCITY
                        #~ times = findTimesToSplit(pDef, inequalities_per_phase)
                        try:
                                times = findTimesToSplit(pDef, inequalities_per_phase)
                                totalMinDist = totalMinDist + 1
                                totalScenarios = totalScenarios + 1
                                #~ plt.show()
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
                                        totalScenarios = totalScenarios + 1
                                        return 
                                except ValueError:
                                        pass
                        return
                                #~ plt.close()

        (P, q, G,h, res) = (None,None,None,None, None)
        for i in range(10):
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

