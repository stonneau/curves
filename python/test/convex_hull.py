import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity, cross
from numpy.linalg import norm

from scipy.spatial import ConvexHull
import numpy as np

def genConvexHullLines(points):
         hull = ConvexHull(points)
         lineList = [points[el] for el in hull.vertices] + [points[hull.vertices[0]]]
         print "lineList", len(lineList)
         print "len(hull.vertices)", len(hull.vertices)
         #~ lineList = points[hull.vertices].tolist()
         #~ lineList = lineList + points[hull.vertices][0].tolist()
         lineList = [array(el[:2].tolist() + [0.]) for el in lineList]
         print "lineList", len(lineList)
         return [[lineList[i],lineList[i+1]] for i in range(len(hull.vertices))]
         #now compute lines


def getLineFromSegment(line):
        a = line[0]; b = line[1]; c = a.copy() ; c[2] = 1.
        normal = cross((b-a),(c-a))
        normal /= norm(normal)
        # get inequality
        dire = b - a
        coeff = normal
        rhs = a.dot(normal)
        #~ signs = [1., -1., -1.]
        #~ rhs = 0.
        #~ for i in range(3):
                #~ if dire[i] != 0.:
                        #~ coeff[i] = signs[i]/dire[i]
                        #~ rhs += coeff[i] * a[i] *(-signs[i])
        
        #~ coeff = array([1. / dire[0], -1. / dire[1], -1. / dire[2]])
        print "coeff", coeff
        #~ print "rhs", rhs
        return (coeff, array([rhs]))

import numpy as np
import matplotlib.pyplot as plt


#generate right of the line
def genFromLine(line, num_points, ranges):
        coeff, rhs = getLineFromSegment(line)
        num_gen = 0
        gen = [line[0][:2], line[1][:2]]
        while(len(gen) < num_points):
                #~ pt = array(np.random.rand(2).tolist()) # + [0.])
                pt = array([np.random.uniform(ranges[0][0], ranges[0][1]), np.random.uniform(ranges[1][0], ranges[1][1])])
                if coeff[:2].dot(pt) <= rhs :
                        gen += [pt]
        print "gen", gen
        return genConvexHullLines(gen)
        #~ step = 1000.
        #~ for line in lines:
                #~ print "one line"
                #~ a_0 = line[0]
                #~ b_0 = line[1]
                #~ pointsline =  np.array([ a_0 * i / step + b_0 * (1. - i / step) for i in range(int(step))])
                #~ xl = pointsline[:,0]
                #~ yl = pointsline[:,1]
                #~ print "heo"
                #~ print "xl", xl
                #~ print "yl", yl
                #~ plt.plot(xl,yl,'b')
                
#~ genFromLine([array([0.5,0.,0.]),array([0.5,0.5,0.])],5)
#~ genFromLine([array([0.5,0.,0.]),array([0.5,-0.5,0.])],5)
        
