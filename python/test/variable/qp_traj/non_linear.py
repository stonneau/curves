
import cvxpy as cp
from numpy.random import randn
from numpy import eye, ones, zeros, array
import numpy as np
from scipy.optimize import minimize

from  cord_methods import *
from  plot_cord import *

from straight_lines import solve_straight_lines, add_acc, bezierFromVal

eps =0.000001

degree = 5
n = degree
nvars = degree+1
dim = 2
rdim = 2

jmp = nvars * dim


def boundIneq():
    v = ones(4)*1.; #v[:-3] = -v[3:]
    V = zeros((4,2)); 
    V[:2,:] = identity(2); 
    V[2:,:] =-identity(2); 
    return (V, v)


def get_x(x,i,xid):
    return x[i*jmp + xid*dim:i*jmp+(xid+1)*dim]

def get_ti(x,i):
    #~ return max(0.00000001, x[i*jmp+nvars*dim + i])
    return x[i*jmp+nvars*dim + i]

def c0(x, i):
    return get_x(x,i,0)
    
def c5(x, i):
    return get_x(x,i,5)
    
def c1(x, i):    
    return c0(x, i) + get_x(x,i,1) * get_ti(x,i)    
    
def c4(x, i):
    return c5(x, i) - get_x(x,i,4) * get_ti(x,i)
    
def c2(x, i):
    return get_x(x,i,2) * get_ti(x,i) * get_ti(x,i) + 2*c1(x, i) - c0(x, i)
    
def c3(x, i):
    return get_x(x,i,3) * get_ti(x,i) * get_ti(x,i) + 2*c4(x, i) - c5(x, i)
    
cjs = [c0,c1,c2,c3,c4,c5]
    
def cj(x,j, i):
    return cjs[j](x,i)
    
def dcj(x,j,i):
    #~ return n / get_ti(x,i) * (cj(x,j+1, i) - cj(x,j, i))
    return n * (cj(x,j+1, i) - cj(x,j, i))
    
def ddcj(x,j,i):
    #~ return (n-1) / get_ti(x,i) * (dcj(x,j+1, i) - dcj(x, j, i))
    return (n-1) * (dcj(x,j+1, i) - dcj(x, j, i))
#~ def position_constraints(Pis):


def posjac0(P,x,i):
    j = 0
    res = zeros((P.shape[0],x.shape[0]))
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim]
    return res
    
def posjac1(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 0
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x0
    j = 1
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] * get_ti(x,i) #x1
    res [:,i*jmp+nvars*dim + i] = -P[:,:rdim].dot(get_x(x,i,1))  # t
    return res
    
def posjac2(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 0
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x0
    j = 1
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = - 2* P[:,:rdim] * get_ti(x,i) #2*x1
    j = 2
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = - P[:,:rdim] * get_ti(x,i) * get_ti(x,i) #x2
    res [:,i*jmp+nvars*dim + i] = - 2 * P[:,:rdim].dot(get_x(x,i,1)) - P[:,:rdim].dot(get_x(x,i,2)) * get_ti(x,i)  # t1 + t1**2
    return res
    
def posjac3(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 5
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x5
    j = 4
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] =  2* P[:,:rdim] * get_ti(x,i) #-2*x4
    j = 3
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] * get_ti(x,i) * get_ti(x,i) #x2
    res [:,i*jmp+nvars*dim + i] =  2 * P[:,:rdim].dot(get_x(x,i,4)) - P[:,:rdim].dot(get_x(x,i,3)) * get_ti(x,i)  # t1 + t1**2
    return res
    
def posjac4(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 5
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x5
    j = 4
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] =  P[:,:rdim] * get_ti(x,i) #- x4
    res [:,i*jmp+nvars*dim + i] = -P[:,:rdim].dot(get_x(x,i,4))  # t
    return res
    
def posjac5(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    res [:,i*jmp + 5*dim:i*jmp+(5+1)*dim] = -P[:,:rdim]
    return res

posjacj = [posjac0, posjac1, posjac2, posjac3, posjac4, posjac5]

def position_cons(P,p,i):
    
    def evaluate(j, i):
        def ev(x):
            return p - P[:,:rdim].dot(cj(x,j,i))
        return ev
        
    def evaljac(j, i):
        def ev(x):            
            return posjacj[j](P,x,i)
        return ev
    
    return [ 
              {'type': 'ineq',
              'fun' : evaluate(j, i) , #contraint is >0 
              #~ 'jac' : evaljac(j, i) , #contraint is >0 
               } #contraint is >0 
              for j in range(nvars)
        ]
       


def veljac0(P,x,i):
    j = 0
    res = zeros((P.shape[0],x.shape[0]))
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim]
    return res
    
def veljac1(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 0
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x0
    j = 1
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] * get_ti(x,i) #x1
    res [:,i*jmp+nvars*dim + i] = -P[:,:rdim].dot(get_x(x,i,1))  # t
    return res
    
def veljac2(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 0
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x0
    j = 1
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = - 2* P[:,:rdim] * get_ti(x,i) #2*x1
    j = 2
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = - P[:,:rdim] * get_ti(x,i) * get_ti(x,i) #x2
    res [:,i*jmp+nvars*dim + i] = - 2 * P[:,:rdim].dot(get_x(x,i,1)) - P[:,:rdim].dot(get_x(x,i,2)) * get_ti(x,i)  # t1 + t1**2
    return res
    
def veljac3(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 5
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x5
    j = 4
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] =  2* P[:,:rdim] * get_ti(x,i) #-2*x4
    j = 3
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] * get_ti(x,i) * get_ti(x,i) #x2
    res [:,i*jmp+nvars*dim + i] =  2 * P[:,:rdim].dot(get_x(x,i,4)) - P[:,:rdim].dot(get_x(x,i,3)) * get_ti(x,i)  # t1 + t1**2
    return res
    
def veljac4(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    j = 5
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] = -P[:,:rdim] #x5
    j = 4
    res [:,i*jmp + j*dim:i*jmp+(j+1)*dim] =  P[:,:rdim] * get_ti(x,i) #- x4
    res [:,i*jmp+nvars*dim + i] = -P[:,:rdim].dot(get_x(x,i,4))  # t
    return res
    
def veljac5(P,x,i):
    res = zeros((P.shape[0],x.shape[0]))
    res [:,i*jmp + 5*dim:i*jmp+(5+1)*dim] = -P[:,:rdim]
    return res

veljacj = [veljac0, veljac1, veljac2, veljac3, veljac4]

def vel_cons(V,v,i):
    
    def evaluate(j, i):
        return lambda x: v * get_ti(x,i) - V.dot(dcj(x,j,i)) * get_ti(x,i)
        
    return [ 
              {'type': 'ineq',
              'fun' : evaluate(j, i) ,} #contraint is >0 
              for j in range(nvars-1)
        ]
        
def acc_cons(A,a,i):
    
    def evaluate(j, i):
        #~ def ev(x):
            #~ print "consV", (v - V.dot(dcj(x,j,i))>=0).all()
            #~ return v - V.dot(dcj(x,j,i))
        return lambda x: a * get_ti(x,i) * get_ti(x,i)- A.dot(ddcj(x,j,i)) * get_ti(x,i) * get_ti(x,i)
        #~ return ev
        
    return [ 
              {'type': 'ineq',
              'fun' : evaluate(j, i) ,} #contraint is >0 
              for j in range(nvars-2)
        ]
    
def t_pos_cons(nphases):
    return [{'type': 'ineq',
              'fun' : lambda x:  x[-nphases:] - ones(nphases)*0.01}]
    
def flatlist(l):
    return [item for sublist in l for item in sublist]
        
def all_position_constraints(Pis, nphases):
    return flatlist ([position_cons(Pis[i][0],Pis[i][1],i) for i in range(nphases)])
    
def all_velocity_constraints(Vv, nphases):
    l = [vel_cons(Vv[0],Vv[1],i) for i in range(nphases)]
    return flatlist (l)
    
def all_acc_constraints(Aa, nphases):
    l = [acc_cons(Aa[0],Aa[1],i) for i in range(nphases)]
    return flatlist (l)
    
def cn_continuity_constraints(i,m):    
    def eq(x):
        return get_x(x,i-1,nvars-1-m) - get_x(x,i,m)
        
    return   {'type': 'eq',
              'fun' : eq ,} #contraint is >0     
    
def all_cn_continuity_constraint(nphases):
    return [cn_continuity_constraints(i,0) for i in range(1,nphases)] + [cn_continuity_constraints(i,1) for i in range(1,nphases)] + [cn_continuity_constraints(i,2) for i in range(1,nphases)]
    #~ return [cn_continuity_constraints(i,0) for i in range(1,nphases)] + [cn_continuity_constraints(i,1) for i in range(1,nphases)]
    #~ return [cn_continuity_constraints(i,0) for i in range(1,nphases)]
    #~ return []
    
def start_state_constraint(start):
    def eqs(x):
        #~ print "start", start 
        #~ print "start", get_x(x,0,0) - start[:rdim] 
        #~ res = zeros(x.shape[0])
        #~ res[:rdim] = get_x(x,0,0) - start[:rdim]
        return get_x(x,0,0) - start[:rdim] 
        
    def jac(x):
        res = zeros(x.shape[0])
        res[:rdim] = ones(rdim)
        return res.T
    return   [{'type': 'eq',
              'fun' : eqs ,
              #~ 'jac' : jac
              } ]
              
def end_state_constraint(nphase, end):
    def eqe(x):
        #~ res = zeros((1,x.shape[0]))
        #~ res[:,-rdim-nphase:-nphase] = get_x(x,nphase-1,nvars-1) - end[:rdim]
        #~ return res
        return get_x(x,nphase-1,nvars-1) - end[:rdim]
        
    def jac(x):
        res = zeros((1,x.shape[0]))
        res[:,-rdim-nphase:-nphase] = ones(rdim)
        return res
        
    return   [{'type': 'eq',
              'fun' : eqe  ,
              #~ 'jac' : jac
              } ]
    
    
def all_constraints(Pis,nphases, start, end):
    Vv = boundIneq()
    return all_position_constraints(Pis,nphases) + all_velocity_constraints(Vv, nphases) + t_pos_cons(nphases)  + all_cn_continuity_constraint(nphases) + start_state_constraint(start) + end_state_constraint(nphases, end)
    #~ return all_position_constraints(Pis,nphases) + all_velocity_constraints(Vv, nphases) + t_pos_cons(nphases)  + all_cn_continuity_constraint(nphases) + start_state_constraint(start) + end_state_constraint(nphases, end)
    #~ return all_position_constraints(Pis,nphases) + all_velocity_constraints(Vv, nphases) + t_pos_cons(nphases)  + all_cn_continuity_constraint(nphases) + start_state_constraint(start) + end_state_constraint(nphases, end)
    #~ return all_position_constraints(Pis,nphases)  + all_cn_continuity_constraint(nphases) + all_velocity_constraints(Vv, nphases) +  t_pos_cons(nphases)  
    #~ return all_position_constraints(Pis,nphases)  + all_cn_continuity_constraint(nphases) + all_velocity_constraints(Vv, nphases)  + end_state_constraint(nphases, end) +  t_pos_cons(nphases) 
    
        
def objective(nphases):
    #~ obj = ones(nvars * dim * nphases + nphases) #times at the end
    obj = zeros(nvars * dim * nphases + nphases) #times at the end
    obj [-nphases:] = ones(nphases)
    def fun(x):
        return obj.dot(x) * 1000.
    return fun
    
def bounds(nphases):
    #~ return[(-1,100) for _ in range(nvars * dim * nphases)] + [(0.0001,None) for _ in range(nphases)]    
    res = []
    for i in range(nphases):
        res = res + [(0.,100.),(0.,100.),(-1.,1),(-1.,1),(-1.,1),(-1.,1),(-1.,1),(-1.,1),(-1.,1),(-1.,1),(0.,100.),(0.,100.)] 
    return res + [(0.0001,None) for _ in range(nphases)]    
 
def acc(xs,ts,nphase, V):
    xis_rec = []
    for i in range(nphase):
        xis = xs[i*4:i*4+4]
        xis = [xis[0], zeros(dim), zeros(dim)] + [zeros(dim), zeros(dim), xis[-1]]
        xis_rec = xis_rec + xis
    return xis_rec

def init_guess(pDef,Pis,nphases, plot=True):
    V = boundIneq()
    x_start = pDef.start.reshape((-1,))
    x_end = pDef.end.reshape((-1,))
    prob, xis, tis = solve_straight_lines(Pis[:], V, x_start=x_start, x_end=x_end)    
    
    #~ print "stragiht ", xis
    if plot:
        for i in range(len(tis)):
            ti = abs(tis[i][0])
            b = bezierFromVal(xis[i*4:i*4+4], abs(ti))
            plotBezier(b, colors[i], label = None, linewidth = 3.0)
            plotControlPoints(b, colors[i],linewidth=2)
    
    xis = acc(xis,tis,nphases,V)
    #~ print "acc ", xis
    l = flatlist([el.tolist() for el in xis]) + tis
    return array(l).reshape((-1,))
    
    

    
    
if __name__ == '__main__':
    
        
    def one(nphases=3):
        plt.close()
        pDef, Pis = genProblemDef(6,nphases)
        x_start = pDef.start.reshape((-1,))
        x_end = pDef.end.reshape((-1,))

        x0 = init_guess(pDef, Pis,nphases, False)   
        #~ x0 = init_guess(pDef, Pis,nphases, False)   
        
        #~ print "bounds", bounds(nphases)
        #~ print "bounds", (x0 >= bounds(nphases)[0]).all()
        #~ print "bounds", (x0 <= bounds(nphases)[1]).all()
        
        #~ print "x0", x0
        
        #~ print "c0", c0(x0, 0)
        #~ print "c1", c1(x0, 0)
        #~ print "c2", c2(x0, 0)
        #~ print "c3", c3(x0, 0)
        #~ print "c4", c4(x0, 0)
        #~ print "c5", c5(x0, 0)
        #~ print "c0", c0(x0, 1)
        #~ print "c1", c1(x0, 1)
        #~ print "c2", c2(x0, 1)
        #~ print "c3", c3(x0, 1)
        #~ print "c4", c4(x0, 1)
        #~ print "c5", c5(x0, 1)
        
        #~ print "bounds", bounds(nphases)[0].shape
        
        #~ plt.show()
        res = minimize(objective(nphases),x0, method='SLSQP',
               #~ constraints=all_constraints(Pis,nphases,x_start,x_end), options={'ftol': 1e-9, 'eps': 1.4901161193847656e-08},'disp': True, 'maxiter' : 500},
               constraints=all_constraints(Pis,nphases,x_start,x_end), options={'ftol': 1e-5,'disp': True, 'maxiter' : 500},
               bounds=bounds(nphases)
               )
               
       
               
        for i in range(nphases):
            wp1 = array([array(cj(res.x,j, i).tolist() + [0.]) for j in range(6)]).T
            b1 = bezier(wp1)
            plotBezier(b1, colors[i], label = None, linewidth = 2.0)
            #~ plotControlPoints(b1, colors[i],linewidth=2)
        #~ wp1 = array([array(cj(res.x,j, 0).tolist() + [0.]) for j in range(6)]).T
        #~ wp2 = array([array(cj(res.x,j, 1).tolist() + [0.]) for j in range(6)]).T
        #~ b1 = bezier(wp1)
        #~ b2 = bezier(wp2)
        #~ plotBezier(b1, "r", label = None, linewidth = 2.0)
        #~ plotControlPoints(b1, "r",linewidth=2)
        #~ plotBezier(b2, "b", label = None, linewidth = 2.0)
        #~ plotControlPoints(b2, "b",linewidth=2)
        plt.show()
        return res
        
    res = one()
    #now plot
