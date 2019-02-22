import matplotlib.pyplot as plt
import numpy as np


def plotBezier(bez, color, label = None, linewidth = 2.0):
        step = 1000.
        points1 =  np.array([(bez(i/step*bez.max())[0][0],bez(i/step*bez.max())[1][0]) for i in range(int(step))])
        x = points1[:,0]
        y = points1[:,1]
        plt.plot(x,y,color,linewidth=linewidth, label=label)
        
def plotControlPoints(bez, color):
        wps = bez.waypoints()
        wps = np.array([wps[:2,i] for i in range(wps.shape[1]) ])
        x = wps[:,0]
        y = wps[:,1]
        plt.scatter(x,y,color=color)        
        
def plotPoly(lines, color):
        step = 1000.
        for line in lines:
                a_0 = line[0]
                b_0 = line[1]
                pointsline =  np.array([ a_0 * i / step + b_0 * (1. - i / step) for i in range(int(step))])
                xl = pointsline[:,0]
                yl = pointsline[:,1]
                plt.plot(xl,yl,color,linewidth=0.5)

def plotFun(bez, start, end, color, label = None, linewidth = 2.0):
        step = 1000.
        points1 =  np.array([[start + i/step*end,bez(start + i/step*end)] for i in range(int(step))])
        x = points1[:,0]
        y = points1[:,1]
        plt.plot(x,y,color,linewidth=linewidth, label=label)


if __name__ == '__main__':
        def cost(x):
                return x**2 + x * 2
        def cost_f(x):
                return 2./3. * x**2 - x * 2./3. * 2
             
             
        plotFun(cost,-3.,6.,'r')
        plotFun(cost_f,-3.,6.,'b')
        
