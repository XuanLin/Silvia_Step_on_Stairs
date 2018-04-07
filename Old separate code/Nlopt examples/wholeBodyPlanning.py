import nlopt
from numpy import *

def myfunc(x, grad):
    if grad.size > 0:
        grad[0] = 0.0
        grad[1] = 0.5 / sqrt(x[1])
    return sqrt(x[1])
def myconstraint(x, grad, a, b):
    if grad.size > 0:
        grad[0] = 3 * a * (a*x[0] + b)**2
        grad[1] = -1.0
    return (a*x[0] + b)**3 - x[1]


#First, need a function to calculate the body COM r by body configuration (thetas)
#Can also take derivative to get body COM velocity by theta_v, and body COM acceleration by theta_accel

#Optimization variables are body COM r, r_dot, r_dot_dot
#leg configurations theta, theta_dot, theta_dot_dot
#body posture is assumed to be fixed (get from plane fitting to stairs)

#Here, walking on flat ground, R = [1,0,0; 0,1,0; 0,0,1]
#hexagon1: [750, 1180, 0], [1250, 1180, 0], [1500, 780, 0], [1250, 380, 0], [750, 380, 0], [500, 780, 0]
#hexagon2: [750, 1326.5, 0], [1250, 1326.5, 0], [1500, 924, 0], [1250, 526.5, 0], [750, 526.5, 0], [500, 924, 0]
#Let's go from [1000, 780, 0] to the center of the second hexagon [1000, 926.5, 0]
#The switching point is precalculated and given as [900, 852, 0], by the middle point of [1500, 780, 0] and [500, 924, 0] offested a little bit

#So r0=[1000, 780, 0], rN=[1000, 926.5, 0], and rN/2=[900, 852, 0]
opt = nlopt.opt(nlopt.LD_MMA, 2)
opt.set_lower_bounds([-float('inf'), 0])
opt.set_min_objective(myfunc)
opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,2,0), 1e-8)
opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,-1,1), 1e-8)
opt.set_xtol_rel(1e-4)
x = opt.optimize([1.234, 5.678])
minf = opt.last_optimum_value()
print("optimum at ", x[0], x[1])
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())

