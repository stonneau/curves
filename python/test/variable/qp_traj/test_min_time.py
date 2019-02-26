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
                        totalScenarios = totalScenarios +1
                        pDef, inequalities_per_phase = genProblemDef(5,4)
                        #~ pDef.costFlag = derivative_flag.DISTANCE
                        #~ pDef.costFlag = derivative_flag.ACCELERATION
                        pDef.costFlag = derivative_flag.VELOCITY
                        #~ res = computeTrajectory(pDef, saveToFile)
                        #~ res = computeTrajectory(pDef, saveToFile)
                        #~ times = findTimesToSplit(pDef, inequalities_per_phase)
                        plt.close()
                        try:
                                times = findTimesToSplit(pDef, inequalities_per_phase)
                                totalMinDist = totalMinDist + 1
                                print "times", times
                                # now solve with and without times
                                #~ pDef.degree = 4
                                #~ pDef.costFlag = derivative_flag.VELOCITY
                                oldsplits = pDef.splits
                                #~ plt.show()
                                plt.close()
                                
                                
                                pDef.splits = array([times]).T 
                                global gb
                                gb = -1
                                try:
                                        res, cos = computeTrajectory(pDef, saveToFile)
                                        totalHeuristicTimes = totalHeuristicTimes +1
                                        # compute best cost 
                                        plt.close()
                                        b = evalBez(res, pDef)
                                        plotBezier(b,'r',None,3.)
                                        plt.show()
                                        avgcost = avgcost + cos
                                except ValueError:
                                        pass
                                        #~ print "failed"
                                        #~ plt.close()
                                #~ print "SUCESS RETIME"
                                try:
                                        pDef.splits = oldsplits
                                        res, cos = computeTrajectory(pDef, saveToFile)
                                        totalRandomTimes = totalRandomTimes +1
                                        avgcostRand = avgcostRand + cos
                                except ValueError:
                                        pass
                                        
                                pDef.degree = 9
                                pDef.costFlag = derivative_flag.VELOCITY
                                oldsplits = pDef.splits
                                
                                pDef.splits = array([times]).T 
                                global gb
                                gb = -1
                                try:
                                        res, cos = computeTrajectory(pDef, saveToFile)
                                        totalHeuristicTimesHighD = totalHeuristicTimesHighD +1
                                        avgcostHighD = avgcostHighD + cos
                                except ValueError:
                                        pass
                                        #~ print "failed"
                                        #~ plt.close()
                                #~ print "SUCESS RETIME"
                                try:
                                        pDef.splits = oldsplits
                                        res, cos = computeTrajectory(pDef, saveToFile)
                                        totalRandomTimesHighD = totalRandomTimesHighD +1
                                        avgcostRandHighD = avgcostRandHighD + cos
                                except ValueError:
                                        pass
                                #~ print "SUCESS RANDOM"
                        
                        
                        #~ plt.show()
                        #~ pDef.costFlag = derivative_flag.ACCELERATION
                        #~ res = computeTrajectory(pDef, saveToFile)
                        #~ pDef.costFlag = derivative_flag.JERK
                        #~ res = computeTrajectory(pDef, saveToFile)
                        #~ plt.legend(loc='upper left')
                        #~ plt.show()
                        #~ plt.close()
                                return res
                        except ValueError:
                                print "failed"
                        #~ except Exception, e:
                        except:
                                totalScenarios = totalScenarios -1
                                print "total fail "
                                #~ plt.close()

        (P, q, G,h, res) = (None,None,None,None, None)
        for i in range(10):
                #~ (P, q, G,h, res) = gen(False)
                gen(False)
                #~ plt.show()
                gb = -1
                #~ if res[0] != None:
                        #~ break
        print "total scenarios", totalScenarios 
        print "total min time success", totalMinDist
        print "total HeuristicTimes", totalHeuristicTimes, ' cost ', avgcost / totalHeuristicTimes
        print "total totalRandomTimes", totalRandomTimes, ' cost ', avgcostRand / totalRandomTimes
        print "total totalHeuristicTimesHighD", totalHeuristicTimesHighD, ' cost ', avgcostHighD / totalHeuristicTimesHighD
        print "total totalRandomtotalRandomTimesHighDTimes", totalRandomTimesHighD, ' cost ', avgcostRandHighD / totalRandomTimesHighD

        def cost(P, q, G,h, x):
                print (x.T.dot(P).dot(x) / 2. + q.dot(x))
                print "ineq ?",  (G.dot(x) -h <=0.0001).all()

        #~ zero = array([2.,2.,0.])
        #~ print "res, " ; cost(P,q, G,h, res)
        #~ print "zero, " ; cost(P,q, G,h, zero)
        #~ print "-zero, " ; cost(P,q, G,h, -zero)

        #~ plt.show()

