from cord_methods import *
from cord_methods import __getbezVar
######################## solve a given problem ########################


colors=['b', 'g', 'r', 'c', 'm', 'y', 'k']
labels=['Distance', 'Velocity', 'Acceleration', 'Jerk']
colors2=[colors[len(colors)-1-i] for i in range(len(colors))]


gb = -1

def plot(pDef, res, filename, saveToFile):        
        global idxFile
        global colors
        global gb
        bezVar = __getbezVar(pDef)
        subs = bezVar.split(pDef.splits.reshape((-1)).tolist())
        final = evalBez(res,pDef)
        derivative = (int)(pDef.costFlag)
        color = colors[(derivative + gb) % len(colors)]
        gb = gb + 1
        for i, bez in enumerate(subs):
                color = colors[i]
                test = bez.toBezier3(res[:])
                if( i == 0):
                        plotBezier(test, color, label = labels[derivative], linewidth = 2)
                else:
                        plotBezier(test, color, linewidth =2)
        if saveToFile:
                plt.savefig(filename+str(idxFile))
        #plot subs control points
        for i, bez in enumerate(subs):
                color = colors[i]
                test = bez.toBezier3(res[:])
                #~ plotBezier(test, color, linewidth = 1)
                #~ plotControlPoints(test, color)
        if saveToFile:
                plt.savefig("subcp"+filename+str(idxFile))
        #~ plotControlPoints(final, "black")
        if saveToFile:
                plt.savefig("cp"+filename+str(idxFile))
        idxFile += 1
        if saveToFile:
                plt.close()
        else:
                plt.show()
        
        return final

#solve and gen problem
def gen(saveToFile = False):
        plt.close()
        try:
        #~ if True:
                pDef, inequalities_per_phase = genProblemDef(10,5)
                pDef.costFlag = derivative_flag.VELOCITY
                res, cost = computeTrajectory(pDef, saveToFile)
                plot(pDef,res, saveToFile=False, filename = uuid.uuid4().hex.upper()[0:6])
                plt.legend(loc='upper left')
                plt.show()
                #~ plt.close()
                return res
        except ValueError:
                print "failed"
                #~ plt.close()

for i in range(1):
        gen(False)
        #~ if res[0] != None:
                #~ break

def cost(P, q, G,h, x):
        print (x.T.dot(P).dot(x) / 2. + q.dot(x))
        print "ineq ?",  (G.dot(x) -h <=0.0001).all()

#~ zero = array([2.,2.,0.])
#~ print "res, " ; cost(P,q, G,h, res)
#~ print "zero, " ; cost(P,q, G,h, zero)
#~ print "-zero, " ; cost(P,q, G,h, -zero)

#~ plt.show()

np.set_printoptions(precision=3, suppress=True, threshold=np.nan)
