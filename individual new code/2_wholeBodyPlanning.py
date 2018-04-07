#First, need a function to calculate the body COM r by body configuration (thetas)
#Can also take derivative to get body COM velocity by theta_v, and body COM acceleration by theta_accel

#Optimization variables are body GC(geometrc center) px, py, pz, px_dot, py_dot, pz_dot, COM cx, cy, cz, cx_dot, cy_dot, cz_dot
#leg configurations thetax18, theta_dotx18,
#body posture is assumed to be fixed (get from plane fitting to stairs)
#constraint 1: c = function(theta, p), c_dot = function(theta, theta_dot, p, p_dot)
#How do deal with accelerations?

#Here, walking on flat ground, R = [1,0,0; 0,1,0; 0,0,1]
#hexagon1: [1250, 1180, 0], [1500, 780, 0], [1250, 380, 0], [750, 380, 0], [500, 780, 0], [750, 1180, 0]
#hexagon2: [1250, 1326.5, 0], [1500, 924, 0], [1250, 526.5, 0], [750, 526.5, 0], [500, 924, 0], [750, 1326.5, 0]
#Let's go from [1000, 780, 0] to the center of the second hexagon [1000, 926.5, 0]
#The switching point is precalculated and given as [900, 852, 0], by the middle point of [1500, 780, 0] and [500, 924, 0] offested a little bit

#So r0=[1000, 780, 0], rN=[1000, 926.5, 0], and rN/2=[900, 852, 0]

#Let's start from only positions first, np.since taking derivatives will make equation very long
#Thus this problem is 24*N dimensions: 18 motor angles and 3 CoM, 3 GC, N is the number of variables along time dimension
#constraint 1: c = function(theta, p)

#constraint 2: r0=[1000, 780, 0], rN=[1000, 926.5, 0], and rN/2=[900, 852, 0]
#constraint 2 is verified to be working
#constraint 3: toepoints = specified toepoints on the ground (with limits on each motor angle), note before N/2 only first 3 toes, after N/2 only second 3 toes, at N/2 all 6 toes

#constraint 4: let's specify a triangular trajectory for the limbs on the air

#constraint 5: CoM within support polygon
#constraint 5 is verified to be working
#constraint 6: each joint angle and geometri center should have limited amount of moving
#constraint 6 is verified to be working, don't set the limit at too low!

import nlopt
import numpy as np
import time

N = 20  #Dimension along time axis
m1 = 0.3
m2 = 0.3
m3 = 0.3
mb = 5
L_coxa = 57
L_femur = 199
L_tibia = 381
Xm1 = 0.0
Xm2 = 0.0
Ym3 = 0.0
x_body_mid_half = 224
x_body_front_half = 111
y_body_half = 194

#r0=[1000, 780, 100]
r0x = 1000.0
r0y = 780.0
r0z = 200.0

#rS=[900, 852, 100]
rSx = 900.0
rSy = 852.0
rSz = 200.0

#rN=[1000, 926.5, 100]
rNx = 1000.0
rNy = 926.5
rNz = 200.0

#------------------------------------ Constraint 1 -----------------------------------------------------------------
def CoM(result, x, grad, step):
    #Step should start from 0
    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])

        grad[0, 0+24*step] = 1.0
        grad[0, 1+24*step] = 0.0
        grad[0, 2+24*step] = 0.0
        grad[0, 3+24*step] = -1.0
        grad[0, 4+24*step] = 0.0
        grad[0, 5+24*step] = 0.0
        grad[0, 6+24*step] = -(L_coxa*m2*np.sin(x[6+24*step]) + L_coxa*m3*np.sin(x[6+24*step]) + Xm1*m1*np.sin(x[6+24*step]) + L_femur*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step]) + Xm2*m2*np.cos(x[7+24*step])*np.sin(x[6+24*step]) - Ym3*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.sin(x[8+24*step]) - Ym3*m3*np.sin(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 7+24*step] = -(L_femur*m3*np.cos(x[6+24*step])*np.sin(x[7+24*step]) + Xm2*m2*np.cos(x[6+24*step])*np.sin(x[7+24*step]) + Ym3*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.cos(x[8+24*step]) - Ym3*m3*np.cos(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 8+24*step] = -(Ym3*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.cos(x[8+24*step]) - Ym3*m3*np.cos(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 9+24*step] =  -(L_coxa*m2*np.sin(x[9+24*step]) + L_coxa*m3*np.sin(x[9+24*step]) + Xm1*m1*np.sin(x[9+24*step]) + L_femur*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step]) + Xm2*m2*np.cos(x[10+24*step])*np.sin(x[9+24*step]) - Ym3*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.sin(x[11+24*step]) - Ym3*m3*np.sin(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 10+24*step] = -(L_femur*m3*np.cos(x[9+24*step])*np.sin(x[10+24*step]) + Xm2*m2*np.cos(x[9+24*step])*np.sin(x[10+24*step]) + Ym3*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.cos(x[11+24*step]) - Ym3*m3*np.cos(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 11+24*step] = -(Ym3*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.cos(x[11+24*step]) - Ym3*m3*np.cos(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 12+24*step] = -(L_coxa*m2*np.sin(x[12+24*step]) + L_coxa*m3*np.sin(x[12+24*step]) + Xm1*m1*np.sin(x[12+24*step]) + L_femur*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step]) + Xm2*m2*np.cos(x[13+24*step])*np.sin(x[12+24*step]) - Ym3*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.sin(x[14+24*step]) - Ym3*m3*np.sin(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 13+24*step] = -(L_femur*m3*np.cos(x[12+24*step])*np.sin(x[13+24*step]) + Xm2*m2*np.cos(x[12+24*step])*np.sin(x[13+24*step]) + Ym3*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.cos(x[14+24*step]) - Ym3*m3*np.cos(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 14+24*step] = -(Ym3*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.cos(x[14+24*step]) - Ym3*m3*np.cos(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 15+24*step] = -(L_coxa*m2*np.sin(x[15+24*step]) + L_coxa*m3*np.sin(x[15+24*step]) + Xm1*m1*np.sin(x[15+24*step]) + L_femur*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step]) + Xm2*m2*np.cos(x[16+24*step])*np.sin(x[15+24*step]) - Ym3*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.sin(x[17+24*step]) - Ym3*m3*np.sin(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 16+24*step] = -(L_femur*m3*np.cos(x[15+24*step])*np.sin(x[16+24*step]) + Xm2*m2*np.cos(x[15+24*step])*np.sin(x[16+24*step]) + Ym3*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.cos(x[17+24*step]) - Ym3*m3*np.cos(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 17+24*step] = -(Ym3*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.cos(x[17+24*step]) - Ym3*m3*np.cos(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 18+24*step] = -(L_coxa*m2*np.sin(x[18+24*step]) + L_coxa*m3*np.sin(x[18+24*step]) + Xm1*m1*np.sin(x[18+24*step]) + L_femur*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step]) + Xm2*m2*np.cos(x[19+24*step])*np.sin(x[18+24*step]) - Ym3*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.sin(x[20+24*step]) - Ym3*m3*np.sin(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 19+24*step] = -(L_femur*m3*np.cos(x[18+24*step])*np.sin(x[19+24*step]) + Xm2*m2*np.cos(x[18+24*step])*np.sin(x[19+24*step]) + Ym3*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.cos(x[20+24*step]) - Ym3*m3*np.cos(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 20+24*step] = -(Ym3*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.cos(x[20+24*step]) - Ym3*m3*np.cos(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 21+24*step] = -(L_coxa*m2*np.sin(x[21+24*step]) + L_coxa*m3*np.sin(x[21+24*step]) + Xm1*m1*np.sin(x[21+24*step]) + L_femur*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step]) + Xm2*m2*np.cos(x[22+24*step])*np.sin(x[21+24*step]) - Ym3*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.sin(x[23+24*step]) - Ym3*m3*np.sin(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 22+24*step] = -(L_femur*m3*np.cos(x[21+24*step])*np.sin(x[22+24*step]) + Xm2*m2*np.cos(x[21+24*step])*np.sin(x[22+24*step]) + Ym3*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.cos(x[23+24*step]) - Ym3*m3*np.cos(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 23+24*step] = -(Ym3*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.cos(x[23+24*step]) - Ym3*m3*np.cos(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

        grad[1, 0+24*step] = 0.0
        grad[1, 1+24*step] = 1.0
        grad[1, 2+24*step] = 0.0
        grad[1, 3+24*step] = 0.0
        grad[1, 4+24*step] = -1.0
        grad[1, 5+24*step] = 0.0
        grad[1, 6+24*step] = (L_coxa*m2*np.cos(x[6+24*step]) + L_coxa*m3*np.cos(x[6+24*step]) + Xm1*m1*np.cos(x[6+24*step]) + L_femur*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step]) + Xm2*m2*np.cos(x[6+24*step])*np.cos(x[7+24*step]) - Ym3*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.sin(x[8+24*step]) - Ym3*m3*np.cos(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 7+24*step] = -(L_femur*m3*np.sin(x[6+24*step])*np.sin(x[7+24*step]) + Xm2*m2*np.sin(x[6+24*step])*np.sin(x[7+24*step]) + Ym3*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.cos(x[8+24*step]) - Ym3*m3*np.sin(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 8+24*step] = -(Ym3*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.cos(x[8+24*step]) - Ym3*m3*np.sin(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 9+24*step] = (L_coxa*m2*np.cos(x[9+24*step]) + L_coxa*m3*np.cos(x[9+24*step]) + Xm1*m1*np.cos(x[9+24*step]) + L_femur*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step]) + Xm2*m2*np.cos(x[9+24*step])*np.cos(x[10+24*step]) - Ym3*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.sin(x[11+24*step]) - Ym3*m3*np.cos(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 10+24*step] = -(L_femur*m3*np.sin(x[9+24*step])*np.sin(x[10+24*step]) + Xm2*m2*np.sin(x[9+24*step])*np.sin(x[10+24*step]) + Ym3*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.cos(x[11+24*step]) - Ym3*m3*np.sin(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 11+24*step] = -(Ym3*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.cos(x[11+24*step]) - Ym3*m3*np.sin(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 12+24*step] =  (L_coxa*m2*np.cos(x[12+24*step]) + L_coxa*m3*np.cos(x[12+24*step]) + Xm1*m1*np.cos(x[12+24*step]) + L_femur*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step]) + Xm2*m2*np.cos(x[12+24*step])*np.cos(x[13+24*step]) - Ym3*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.sin(x[14+24*step]) - Ym3*m3*np.cos(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 13+24*step] =  -(L_femur*m3*np.sin(x[12+24*step])*np.sin(x[13+24*step]) + Xm2*m2*np.sin(x[12+24*step])*np.sin(x[13+24*step]) + Ym3*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.cos(x[14+24*step]) - Ym3*m3*np.sin(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 14+24*step] =  -(Ym3*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.cos(x[14+24*step]) - Ym3*m3*np.sin(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 15+24*step] = (L_coxa*m2*np.cos(x[15+24*step]) + L_coxa*m3*np.cos(x[15+24*step]) + Xm1*m1*np.cos(x[15+24*step]) + L_femur*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step]) + Xm2*m2*np.cos(x[15+24*step])*np.cos(x[16+24*step]) - Ym3*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.sin(x[17+24*step]) - Ym3*m3*np.cos(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 16+24*step] = -(L_femur*m3*np.sin(x[15+24*step])*np.sin(x[16+24*step]) + Xm2*m2*np.sin(x[15+24*step])*np.sin(x[16+24*step]) + Ym3*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.cos(x[17+24*step]) - Ym3*m3*np.sin(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 17+24*step] = -(Ym3*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.cos(x[17+24*step]) - Ym3*m3*np.sin(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 18+24*step] = (L_coxa*m2*np.cos(x[18+24*step]) + L_coxa*m3*np.cos(x[18+24*step]) + Xm1*m1*np.cos(x[18+24*step]) + L_femur*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step]) + Xm2*m2*np.cos(x[18+24*step])*np.cos(x[19+24*step]) - Ym3*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.sin(x[20+24*step]) - Ym3*m3*np.cos(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 19+24*step] = -(L_femur*m3*np.sin(x[18+24*step])*np.sin(x[19+24*step]) + Xm2*m2*np.sin(x[18+24*step])*np.sin(x[19+24*step]) + Ym3*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.cos(x[20+24*step]) - Ym3*m3*np.sin(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 20+24*step] = -(Ym3*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.cos(x[20+24*step]) - Ym3*m3*np.sin(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 21+24*step] =  (L_coxa*m2*np.cos(x[21+24*step]) + L_coxa*m3*np.cos(x[21+24*step]) + Xm1*m1*np.cos(x[21+24*step]) + L_femur*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step]) + Xm2*m2*np.cos(x[21+24*step])*np.cos(x[22+24*step]) - Ym3*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.sin(x[23+24*step]) - Ym3*m3*np.cos(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 22+24*step] = -(L_femur*m3*np.sin(x[21+24*step])*np.sin(x[22+24*step]) + Xm2*m2*np.sin(x[21+24*step])*np.sin(x[22+24*step]) + Ym3*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.cos(x[23+24*step]) - Ym3*m3*np.sin(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 23+24*step] = -(Ym3*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.cos(x[23+24*step]) - Ym3*m3*np.sin(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

        grad[2, 0+24*step] = 0.0
        grad[2, 1+24*step] = 0.0
        grad[2, 2+24*step] = 1.0
        grad[2, 3+24*step] = 0.0
        grad[2, 4+24*step] = 0.0
        grad[2, 5+24*step] = -1.0
        grad[2, 6+24*step] = 0.0
        grad[2, 7+24*step] =  -(m3*(Ym3*np.sin(x[7+24*step] + x[8+24*step]) - L_femur*np.cos(x[7+24*step])) - Xm2*m2*np.cos(x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 8+24*step] = -(Ym3*m3*np.sin(x[7+24*step] + x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 9+24*step] = 0.0
        grad[2, 10+24*step] = -(m3*(Ym3*np.sin(x[10+24*step] + x[11+24*step]) - L_femur*np.cos(x[10+24*step])) - Xm2*m2*np.cos(x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 11+24*step] =  -(Ym3*m3*np.sin(x[10+24*step] + x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 12+24*step] = 0.0
        grad[2, 13+24*step] = -(m3*(Ym3*np.sin(x[13+24*step] + x[14+24*step]) - L_femur*np.cos(x[13+24*step])) - Xm2*m2*np.cos(x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 14+24*step] =  -(Ym3*m3*np.sin(x[13+24*step] + x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 15+24*step] = 0.0
        grad[2, 16+24*step] = -(m3*(Ym3*np.sin(x[16+24*step] + x[17+24*step]) - L_femur*np.cos(x[16+24*step])) - Xm2*m2*np.cos(x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 17+24*step] = -(Ym3*m3*np.sin(x[16+24*step] + x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 18+24*step] = 0.0
        grad[2, 19+24*step] = -(m3*(Ym3*np.sin(x[19+24*step] + x[20+24*step]) - L_femur*np.cos(x[19+24*step])) - Xm2*m2*np.cos(x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 20+24*step] = -(Ym3*m3*np.sin(x[19+24*step] + x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 21+24*step] = 0.0
        grad[2, 22+24*step] = -(m3*(Ym3*np.sin(x[22+24*step] + x[23+24*step]) - L_femur*np.cos(x[22+24*step])) - Xm2*m2*np.cos(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 23+24*step] = -(Ym3*m3*np.sin(x[22+24*step] + x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

    #This CM has nothing to do with body dimension, think about vector summation to convince yourself
    result[:] = np.array(
        [x[0+24*step] - x[3+24*step] + \
        (L_coxa*m2*np.cos(x[6+24*step]) + L_coxa*m2*np.cos(x[9+24*step]) + L_coxa*m3*np.cos(x[6+24*step]) + L_coxa*m2*np.cos(x[12+24*step]) +
        L_coxa*m3*np.cos(x[9+24*step]) + L_coxa*m2*np.cos(x[15+24*step]) + L_coxa*m3*np.cos(x[12+24*step]) + L_coxa*m2*np.cos(x[18+24*step]) +
        L_coxa*m3*np.cos(x[15+24*step]) + L_coxa*m2*np.cos(x[21+24*step]) + L_coxa*m3*np.cos(x[18+24*step]) + L_coxa*m3*np.cos(x[21+24*step]) +
        L_femur*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step]) + L_femur*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step]) + L_femur*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step]) +
        L_femur*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step]) + L_femur*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step]) + L_femur*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step]) +
        Xm1*m1*np.cos(x[6+24*step]) + Xm1*m1*np.cos(x[9+24*step]) + Xm1*m1*np.cos(x[12+24*step]) +
        Xm1*m1*np.cos(x[15+24*step]) + Xm1*m1*np.cos(x[18+24*step]) + Xm1*m1*np.cos(x[21+24*step]) +
        Xm2*m2*np.cos(x[6+24*step])*np.cos(x[7+24*step]) + Xm2*m2*np.cos(x[9+24*step])*np.cos(x[10+24*step]) + Xm2*m2*np.cos(x[12+24*step])*np.cos(x[13+24*step]) +
        Xm2*m2*np.cos(x[15+24*step])*np.cos(x[16+24*step]) + Xm2*m2*np.cos(x[18+24*step])*np.cos(x[19+24*step]) + Xm2*m2*np.cos(x[21+24*step])*np.cos(x[22+24*step]) -
        Ym3*m3*np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.sin(x[8+24*step]) - Ym3*m3*np.cos(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step]) -
        Ym3*m3*np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.sin(x[11+24*step]) - Ym3*m3*np.cos(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step]) -
        Ym3*m3*np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.sin(x[14+24*step]) - Ym3*m3*np.cos(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step]) -
        Ym3*m3*np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.sin(x[17+24*step]) - Ym3*m3*np.cos(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step]) -
        Ym3*m3*np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.sin(x[20+24*step]) - Ym3*m3*np.cos(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step]) -
        Ym3*m3*np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.sin(x[23+24*step]) - Ym3*m3*np.cos(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb),

        x[1+24*step] - x[4+24*step] + \
        (L_coxa*m2*np.sin(x[6+24*step]) + L_coxa*m2*np.sin(x[9+24*step]) + L_coxa*m3*np.sin(x[6+24*step]) + L_coxa*m2*np.sin(x[12+24*step]) +
        L_coxa*m3*np.sin(x[9+24*step]) + L_coxa*m2*np.sin(x[15+24*step]) + L_coxa*m3*np.sin(x[12+24*step]) + L_coxa*m2*np.sin(x[18+24*step]) +
        L_coxa*m3*np.sin(x[15+24*step]) + L_coxa*m2*np.sin(x[21+24*step]) + L_coxa*m3*np.sin(x[18+24*step]) + L_coxa*m3*np.sin(x[21+24*step]) +
        L_femur*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step]) + L_femur*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step]) + L_femur*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step]) +
        L_femur*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step]) + L_femur*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step]) + L_femur*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step]) +
        Xm1*m1*np.sin(x[6+24*step]) + Xm1*m1*np.sin(x[9+24*step]) + Xm1*m1*np.sin(x[12+24*step]) +
        Xm1*m1*np.sin(x[15+24*step]) + Xm1*m1*np.sin(x[18+24*step]) + Xm1*m1*np.sin(x[21+24*step]) +
        Xm2*m2*np.cos(x[7+24*step])*np.sin(x[6+24*step]) + Xm2*m2*np.cos(x[10+24*step])*np.sin(x[9+24*step]) + Xm2*m2*np.cos(x[13+24*step])*np.sin(x[12+24*step]) +
        Xm2*m2*np.cos(x[16+24*step])*np.sin(x[15+24*step]) + Xm2*m2*np.cos(x[19+24*step])*np.sin(x[18+24*step]) + Xm2*m2*np.cos(x[22+24*step])*np.sin(x[21+24*step]) -
        Ym3*m3*np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.sin(x[8+24*step]) - Ym3*m3*np.sin(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step]) -
        Ym3*m3*np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.sin(x[11+24*step]) - Ym3*m3*np.sin(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step]) -
        Ym3*m3*np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.sin(x[14+24*step]) - Ym3*m3*np.sin(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step]) -
        Ym3*m3*np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.sin(x[17+24*step]) - Ym3*m3*np.sin(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step]) -
        Ym3*m3*np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.sin(x[20+24*step]) - Ym3*m3*np.sin(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step]) -
        Ym3*m3*np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.sin(x[23+24*step]) - Ym3*m3*np.sin(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb),

        x[2+24*step] - x[5+24*step] + \
        (m3*(Ym3*np.cos(x[7+24*step] + x[8+24*step]) + L_femur*np.sin(x[7+24*step])) + m3*(Ym3*np.cos(x[10+24*step] + x[11+24*step]) + L_femur*np.sin(x[10+24*step])) +
        m3*(Ym3*np.cos(x[13+24*step] + x[14+24*step]) + L_femur*np.sin(x[13+24*step])) + m3*(Ym3*np.cos(x[16+24*step] + x[17+24*step]) + L_femur*np.sin(x[16+24*step])) +
        m3*(Ym3*np.cos(x[19+24*step] + x[20+24*step]) + L_femur*np.sin(x[19+24*step])) + m3*(Ym3*np.cos(x[22+24*step] + x[23+24*step]) + L_femur*np.sin(x[22+24*step])) +
        Xm2*m2*np.sin(x[7+24*step]) + Xm2*m2*np.sin(x[10+24*step]) + Xm2*m2*np.sin(x[13+24*step]) + Xm2*m2*np.sin(x[16+24*step]) + Xm2*m2*np.sin(x[19+24*step]) + Xm2*m2*np.sin(x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)])

    return result

#------------------------------------ Constraint 2 -----------------------------------------------------------------
def r0(result, x, grad):

    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 0] = 1.0
        grad[1, 1] = 1.0
        grad[2, 2] = 1.0

    result[:] = np.array([x[0]-r0x, x[1]-r0y, x[2]-r0z])
    return result

def rS(result, x, grad):


    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 24*(N/2-1)+0] = 1.0
        grad[1, 24*(N/2-1)+1] = 1.0
        grad[2, 24*(N/2-1)+2] = 1.0

    result[:] = np.array([x[24*(N/2-1)]-rSx, x[24*(N/2-1)+1]-rSy, x[24*(N/2-1)+2]-rSz])
    return result

def rN(result, x, grad):

    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 24*(N-1)+0] = 1.0
        grad[1, 24*(N-1)+1] = 1.0
        grad[2, 24*(N-1)+2] = 1.0

    result[:] = np.array([x[24*(N-1)]-rNx, x[24*(N-1)+1]-rNy, x[24*(N-1)+2]-rNz])
    return result

#------------------------------------ Constraint 3 -----------------------------------------------------------------
def leg1_end_point_XYZ(result, x, grad, step, px, py, pz):

    #This is RF leg, the algebra is compared with matlab to verify it's correctness. Should carryout actual number to double check it!
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[6+24*step]) + L_tibia*(np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.sin(x[8+24*step]) + np.cos(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) +
        L_femur*np.cos(x[6+24*step])*np.cos(x[7+24*step]) + x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[6+24*step]) + L_tibia*(np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.sin(x[8+24*step]) + np.sin(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) +
        L_femur*np.cos(x[7+24*step])*np.sin(x[6+24*step]) + y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(x[7+24*step]) - L_tibia*(np.cos(x[7+24*step])*np.cos(x[8+24*step]) - np.sin(x[7+24*step])*np.sin(x[8+24*step]))])

    print result    #Print this to verify convergence, should be 0 when converge!

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 6+24*step] = - L_coxa*np.sin(x[6+24*step]) - L_tibia*(np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.sin(x[8+24*step]) + np.sin(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) - \
                             L_femur*np.cos(x[7+24*step])*np.sin(x[6+24*step])
        grad[0, 7+24*step] = L_tibia*(np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.cos(x[8+24*step]) - np.cos(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step])) - \
                             L_femur*np.cos(x[6+24*step])*np.sin(x[7+24*step])
        grad[0, 8+24*step] = L_tibia*(np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.cos(x[8+24*step]) - np.cos(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 6+24*step] = L_coxa*np.cos(x[6+24*step]) + L_tibia*(np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.sin(x[8+24*step]) + np.cos(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) + \
                             L_femur*np.cos(x[6+24*step])*np.cos(x[7+24*step])
        grad[1, 7+24*step] = L_tibia*(np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.cos(x[8+24*step]) - np.sin(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step])) - \
                             L_femur*np.sin(x[6+24*step])*np.sin(x[7+24*step])
        grad[1, 8+24*step] = L_tibia*(np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.cos(x[8+24*step]) - np.sin(x[6+24*step])*np.sin(x[7+24*step])*np.sin(x[8+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 6+24*step] = 0.0
        grad[2, 7+24*step] = L_tibia * (np.cos(x[7+24*step]) * np.sin(x[8+24*step]) + np.cos(x[8+24*step]) * np.sin(x[7+24*step])) + L_femur * np.cos(x[7+24*step])
        grad[2, 8+24*step] = L_tibia * (np.cos(x[7+24*step]) * np.sin(x[8+24*step]) + np.cos(x[8+24*step]) * np.sin(x[7+24*step]))

    return result

def leg2_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is RM leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[9+24*step]) + L_tibia*(np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) +
        L_femur*np.cos(x[9+24*step])*np.cos(x[10+24*step]) + x_body_mid_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[9+24*step]) + L_tibia*(np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.sin(x[11+24*step]) + np.sin(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) +
        L_femur*np.cos(x[10+24*step])*np.sin(x[9+24*step]),

        x[2+24*step] - pz +
        L_femur*np.sin(x[10+24*step]) - L_tibia*(np.cos(x[10+24*step])*np.cos(x[11+24*step]) - np.sin(x[10+24*step])*np.sin(x[11+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 9+24*step] = - L_coxa*np.sin(x[9+24*step]) - L_tibia*(np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.sin(x[11+24*step]) + np.sin(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) - \
                             L_femur*np.cos(x[10+24*step])*np.sin(x[9+24*step])
        grad[0, 10+24*step] = L_tibia*(np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.cos(x[11+24*step]) - np.cos(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step])) - \
                              L_femur*np.cos(x[9+24*step])*np.sin(x[10+24*step])
        grad[0, 11+24*step] =  L_tibia*(np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.cos(x[11+24*step]) - np.cos(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 9+24*step] = L_coxa*np.cos(x[9+24*step]) + L_tibia*(np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) + \
                             L_femur*np.cos(x[9+24*step])*np.cos(x[10+24*step])
        grad[1, 10+24*step] = L_tibia*(np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.cos(x[11+24*step]) - np.sin(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step])) - \
                              L_femur*np.sin(x[9+24*step])*np.sin(x[10+24*step])
        grad[1, 11+24*step] =  L_tibia*(np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.cos(x[11+24*step]) - np.sin(x[9+24*step])*np.sin(x[10+24*step])*np.sin(x[11+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 9+24*step] = 0.0
        grad[2, 10+24*step] = L_tibia*(np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[11+24*step])*np.sin(x[10+24*step])) + L_femur*np.cos(x[10+24*step])
        grad[2, 11+24*step] = L_tibia*(np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[11+24*step])*np.sin(x[10+24*step]))

    return result

def leg3_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is RR leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[12+24*step]) + L_tibia*(np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) +
        L_femur*np.cos(x[12+24*step])*np.cos(x[13+24*step]) + x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[12+24*step]) + L_tibia*(np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.sin(x[14+24*step]) + np.sin(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) +
        L_femur*np.cos(x[13+24*step])*np.sin(x[12+24*step]) - y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(x[13+24*step]) - L_tibia*(np.cos(x[13+24*step])*np.cos(x[14+24*step]) - np.sin(x[13+24*step])*np.sin(x[14+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 12+24*step] =  - L_coxa*np.sin(x[12+24*step]) - L_tibia*(np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.sin(x[14+24*step]) + np.sin(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) - \
                               L_femur*np.cos(x[13+24*step])*np.sin(x[12+24*step])
        grad[0, 13+24*step] = L_tibia*(np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.cos(x[14+24*step]) - np.cos(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step])) - \
                              L_femur*np.cos(x[12+24*step])*np.sin(x[13+24*step])
        grad[0, 14+24*step] = L_tibia*(np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.cos(x[14+24*step]) - np.cos(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 12+24*step] = L_coxa*np.cos(x[12+24*step]) + L_tibia*(np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) + \
                              L_femur*np.cos(x[12+24*step])*np.cos(x[13+24*step])
        grad[1, 13+24*step] =  L_tibia*(np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.cos(x[14+24*step]) - np.sin(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step])) - \
                               L_femur*np.sin(x[12+24*step])*np.sin(x[13+24*step])
        grad[1, 14+24*step] = L_tibia*(np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.cos(x[14+24*step]) - np.sin(x[12+24*step])*np.sin(x[13+24*step])*np.sin(x[14+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 12+24*step] = 0.0
        grad[2, 13+24*step] = L_tibia*(np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[14+24*step])*np.sin(x[13+24*step])) + L_femur*np.cos(x[13+24*step])
        grad[2, 14+24*step] = L_tibia*(np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[14+24*step])*np.sin(x[13+24*step]))

    return result


def leg4_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LR leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[15+24*step]) + L_tibia*(np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) +
        L_femur*np.cos(x[15+24*step])*np.cos(x[16+24*step]) - x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[15+24*step]) + L_tibia*(np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.sin(x[17+24*step]) + np.sin(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) +
        L_femur*np.cos(x[16+24*step])*np.sin(x[15+24*step]) - y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(x[16+24*step]) - L_tibia*(np.cos(x[16+24*step])*np.cos(x[17+24*step]) - np.sin(x[16+24*step])*np.sin(x[17+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 15+24*step] = - L_coxa*np.sin(x[15+24*step]) - L_tibia*(np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.sin(x[17+24*step]) + np.sin(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) - \
                              L_femur*np.cos(x[16+24*step])*np.sin(x[15+24*step])
        grad[0, 16+24*step] = L_tibia*(np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.cos(x[17+24*step]) - np.cos(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step])) - \
                              L_femur*np.cos(x[15+24*step])*np.sin(x[16+24*step])
        grad[0, 17+24*step] = L_tibia*(np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.cos(x[17+24*step]) - np.cos(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 15+24*step] = L_coxa*np.cos(x[15+24*step]) + L_tibia*(np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) + \
                              L_femur*np.cos(x[15+24*step])*np.cos(x[16+24*step])
        grad[1, 16+24*step] = L_tibia*(np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.cos(x[17+24*step]) - np.sin(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step])) - \
                              L_femur*np.sin(x[15+24*step])*np.sin(x[16+24*step])
        grad[1, 17+24*step] = L_tibia*(np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.cos(x[17+24*step]) - np.sin(x[15+24*step])*np.sin(x[16+24*step])*np.sin(x[17+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 15+24*step] = 0.0
        grad[2, 16+24*step] = L_tibia*(np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[17+24*step])*np.sin(x[16+24*step])) + L_femur*np.cos(x[16+24*step])
        grad[2, 17+24*step] = L_tibia*(np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[17+24*step])*np.sin(x[16+24*step]))

    return result

def leg5_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LM leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[18+24*step]) +
        L_tibia*(np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) +
        L_femur*np.cos(x[18+24*step])*np.cos(x[19+24*step]) - x_body_mid_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[18+24*step]) + L_tibia*(np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.sin(x[20+24*step]) + np.sin(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) +
        L_femur*np.cos(x[19+24*step])*np.sin(x[18+24*step]),

        x[2+24*step] - pz +
        L_femur*np.sin(x[19+24*step]) - L_tibia*(np.cos(x[19+24*step])*np.cos(x[20+24*step]) - np.sin(x[19+24*step])*np.sin(x[20+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 18+24*step] = - L_coxa*np.sin(x[18+24*step]) - L_tibia*(np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.sin(x[20+24*step]) + np.sin(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) - \
                              L_femur*np.cos(x[19+24*step])*np.sin(x[18+24*step])
        grad[0, 19+24*step] = L_tibia*(np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.cos(x[20+24*step]) - np.cos(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step])) - \
                              L_femur*np.cos(x[18+24*step])*np.sin(x[19+24*step])
        grad[0, 20+24*step] = L_tibia*(np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.cos(x[20+24*step]) - np.cos(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 18+24*step] = L_coxa*np.cos(x[18+24*step]) + L_tibia*(np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) + \
                              L_femur*np.cos(x[18+24*step])*np.cos(x[19+24*step])
        grad[1, 19+24*step] = L_tibia*(np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.cos(x[20+24*step]) - np.sin(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step])) - \
                              L_femur*np.sin(x[18+24*step])*np.sin(x[19+24*step])
        grad[1, 20+24*step] = L_tibia*(np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.cos(x[20+24*step]) - np.sin(x[18+24*step])*np.sin(x[19+24*step])*np.sin(x[20+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 18+24*step] = 0.0
        grad[2, 19+24*step] = L_tibia*(np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[20+24*step])*np.sin(x[19+24*step])) + L_femur*np.cos(x[19+24*step])
        grad[2, 20+24*step] = L_tibia*(np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[20+24*step])*np.sin(x[19+24*step]))

    return result

def leg6_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LF leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(x[21+24*step]) + L_tibia*(np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) +
        L_femur*np.cos(x[21+24*step])*np.cos(x[22+24*step]) - x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(x[21+24*step]) + L_tibia*(np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.sin(x[23+24*step]) + np.sin(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) +
        L_femur*np.cos(x[22+24*step])*np.sin(x[21+24*step]) + y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(x[22+24*step]) - L_tibia*(np.cos(x[22+24*step])*np.cos(x[23+24*step]) - np.sin(x[22+24*step])*np.sin(x[23+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 21+24*step] = - L_coxa*np.sin(x[21+24*step]) - L_tibia*(np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.sin(x[23+24*step]) + np.sin(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) - \
                              L_femur*np.cos(x[22+24*step])*np.sin(x[21+24*step])
        grad[0, 22+24*step] = L_tibia*(np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.cos(x[23+24*step]) - np.cos(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step])) - \
                              L_femur*np.cos(x[21+24*step])*np.sin(x[22+24*step])
        grad[0, 23+24*step] = L_tibia*(np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.cos(x[23+24*step]) - np.cos(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))

        grad[1, 1+24*step] = 1.0
        grad[1, 21+24*step] =  L_coxa*np.cos(x[21+24*step]) + L_tibia*(np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) + \
                               L_femur*np.cos(x[21+24*step])*np.cos(x[22+24*step])
        grad[1, 22+24*step] = L_tibia*(np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.cos(x[23+24*step]) - np.sin(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step])) - \
                              L_femur*np.sin(x[21+24*step])*np.sin(x[22+24*step])
        grad[1, 23+24*step] = L_tibia*(np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.cos(x[23+24*step]) - np.sin(x[21+24*step])*np.sin(x[22+24*step])*np.sin(x[23+24*step]))

        grad[2, 2+24*step] = 1.0
        grad[2, 21+24*step] = 0.0
        grad[2, 22+24*step] = L_tibia*(np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[23+24*step])*np.sin(x[22+24*step])) + L_femur*np.cos(x[22+24*step])
        grad[2, 23+24*step] = L_tibia*(np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[23+24*step])*np.sin(x[22+24*step]))

    return result

#------------------------------------ Constraint 4 -----------------------------------------------------------------
def leg1_end_point_Z(x, grad, step, pz):
    #This is RF leg, the algebra is compared with matlab to verify it's correctness. Should carryout actual number to double check it!
    result = x[2+24*step] - pz + L_femur*np.sin(x[7+24*step]) - L_tibia*(np.cos(x[7+24*step])*np.cos(x[8+24*step]) - np.sin(x[7+24*step])*np.sin(x[8+24*step]))

    print result    #Print this to verify convergence, should be 0 when converge!

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[6+24*step] = 0.0
        grad[7+24*step] = -(L_tibia * (np.cos(x[7+24*step]) * np.sin(x[8+24*step]) + np.cos(x[8+24*step]) * np.sin(x[7+24*step])) + L_femur * np.cos(x[7+24*step]))
        grad[8+24*step] = -(L_tibia * (np.cos(x[7+24*step]) * np.sin(x[8+24*step]) + np.cos(x[8+24*step]) * np.sin(x[7+24*step])))

    return -result

def leg2_end_point_Z(x, grad, step, pz):
    # This is RM leg
    result = x[2+24*step] - pz + L_femur*np.sin(x[10+24*step]) - L_tibia*(np.cos(x[10+24*step])*np.cos(x[11+24*step]) - np.sin(x[10+24*step])*np.sin(x[11+24*step]))

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[9+24*step] = 0.0
        grad[10+24*step] = -(L_tibia*(np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[11+24*step])*np.sin(x[10+24*step])) + L_femur*np.cos(x[10+24*step]))
        grad[11+24*step] = -(L_tibia*(np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[11+24*step])*np.sin(x[10+24*step])))

    return -result

def leg3_end_point_Z(x, grad, step, pz):
    # This is RR leg
    result = x[2+24*step] - pz + L_femur*np.sin(x[13+24*step]) - L_tibia*(np.cos(x[13+24*step])*np.cos(x[14+24*step]) - np.sin(x[13+24*step])*np.sin(x[14+24*step]))

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[12+24*step] = 0.0
        grad[13+24*step] = -(L_tibia*(np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[14+24*step])*np.sin(x[13+24*step])) + L_femur*np.cos(x[13+24*step]))
        grad[14+24*step] = -(L_tibia*(np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[14+24*step])*np.sin(x[13+24*step])))

    return -result


def leg4_end_point_Z(x, grad, step, pz):
    # This is LR leg
    result = x[2+24*step] - pz + L_femur*np.sin(x[16+24*step]) - L_tibia*(np.cos(x[16+24*step])*np.cos(x[17+24*step]) - np.sin(x[16+24*step])*np.sin(x[17+24*step]))

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[15+24*step] = 0.0
        grad[16+24*step] = -(L_tibia*(np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[17+24*step])*np.sin(x[16+24*step])) + L_femur*np.cos(x[16+24*step]))
        grad[17+24*step] = -(L_tibia*(np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[17+24*step])*np.sin(x[16+24*step])))

    return -result

def leg5_end_point_Z(x, grad, step, pz):
    # This is LM leg
    result = x[2+24*step] - pz + L_femur*np.sin(x[19+24*step]) - L_tibia*(np.cos(x[19+24*step])*np.cos(x[20+24*step]) - np.sin(x[19+24*step])*np.sin(x[20+24*step]))

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[18+24*step] = 0.0
        grad[19+24*step] = -(L_tibia*(np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[20+24*step])*np.sin(x[19+24*step])) + L_femur*np.cos(x[19+24*step]))
        grad[20+24*step] = -(L_tibia*(np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[20+24*step])*np.sin(x[19+24*step])))

    return -result

def leg6_end_point_Z(x, grad, step, pz):
    # This is LF leg
    result = x[2+24*step] - pz + L_femur*np.sin(x[22+24*step]) - L_tibia*(np.cos(x[22+24*step])*np.cos(x[23+24*step]) - np.sin(x[22+24*step])*np.sin(x[23+24*step]))

    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[21+24*step] = 0.0
        grad[22+24*step] = -(L_tibia*(np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[23+24*step])*np.sin(x[22+24*step])) + L_femur*np.cos(x[22+24*step]))
        grad[23+24*step] = -(L_tibia*(np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[23+24*step])*np.sin(x[22+24*step])))

    return -result

#------------------------------------ Constraint 5 -----------------------------------------------------------------
#if step < N/2, then leg: right, right, left => leg1, leg3, leg5
#if step > N/2, then leg: left, left, right => leg2, leg4, leg6
#hexagon1: [1250, 1180, 0], [1500, 780, 0], [1250, 380, 0], [750, 380, 0], [500, 780, 0], [750, 1180, 0]
#hexagon2: [1250, 1326.5, 0], [1500, 924, 0], [1250, 526.5, 0], [750, 526.5, 0], [500, 924, 0], [750, 1326.5, 0]

# support polygon 1: [x1=1250, y1=1180], [x2=1250, y2=380], [x3=500, y3=780]

# ((y3-y2)/(x3-x2))*CoMx-CoMy+(x3*y2-y3*x2)/(x3-x2)<0
# (-(y3-y1)/(x3-x1))*CoMx+CoMy-(x3*y1-y3*x1)/(x3-x1)<0
# CoMx-x1<0

# 0.53333*CoMx+CoMy-1580<0
# 0.53333*CoMx-CoMy-20<0
# 750-CoMx<0

# support polygon 2: [x1=1500, y1=924], [x2=750, y2=526.5], [x3=750, y3=1326.5]

# ((y2-y1)/(x2-x1))*CoMx-CoMy+(x2*y1-y2*x1)/(x2-x1)<0
# (-(y3-y1)/(x3-x1))*CoMx+CoMy-(x3*y1-y3*x1)/(x3-x1)<0
# x3-CoMx<0

# -0.53667*CoMx+CoMy-655.66667<0
# -0.53*CoMx-CoMy+1189<0
# CoMx-1250<0
def supportHexagonConstraint(result, x, grad, step):
    CoMx = x[step*24+3]
    CoMy = x[step*24+4]

    if (step <= N/2):
       #Polygon 1
       x1 = 1250
       y1 = 1180

       x2 = 1250
       y2 = 380

       x3 = 500
       y3 = 780
       result[:] = np.array([((y3-y2)/(x3-x2))*CoMx-CoMy+(x3*y2-y3*x2)/(x3-x2), (-(y3-y1)/(x3-x1))*CoMx+CoMy-(x3*y1-y3*x1)/(x3-x1), CoMx-x1])
       if grad.size > 0:
            grad[:] = np.zeros([3, 24*N])

            grad[0, step*24+3] = ((y3-y2)/(x3-x2))
            grad[0, step*24+4] = -1.0

            grad[1, step*24+3] = (-(y3-y1)/(x3-x1))
            grad[1, step*24+4] = 1.0

            grad[2, step*24+3] = 1.0

    else:
       #Polygon 2
       x1 = 1500
       y1 = 924

       x2 = 750
       y2 = 526.5

       x3 = 750
       y3 = 1326.5
       result[:] = np.array([((y2-y1)/(x2-x1))*CoMx-CoMy+(x2*y1-y2*x1)/(x2-x1), (-(y3-y1)/(x3-x1))*CoMx+CoMy-(x3*y1-y3*x1)/(x3-x1), x3-CoMx])
       if grad.size > 0:
            grad[:] = np.zeros([3, 24*N])

            grad[0, step*24+3] = ((y2-y1)/(x2-x1))
            grad[0, step*24+4] = -1.0

            grad[1, step*24+3] = (-(y3-y1)/(x3-x1))
            grad[1, step*24+4] = 1.0

            grad[2, step*24+3] = -1.0

    return result

#---------------------------------- Constraint 6 -------------------------------------------------------------------
def rotSpdLimUp(result, x, grad, step):
    #Note here step should range from 0 to N-2
    rotationalSpeedLimit = np.pi/8
    if grad.size > 0:
        grad[:] = np.zeros([18, 24 * N])
        for i in range(18):
            grad[i, 6+i+step*24] = 1.0
            grad[i, 6+i+(step+1)*24] = -1.0

    for i in range(18):
        result[i] = x[6+i+step*24] - x[6+i+(step+1)*24] - rotationalSpeedLimit
    return result

def rotSpdLimLow(result, x, grad, step):
    #Note here step should range from 0 to N-2
    rotationalSpeedLimit = np.pi/8
    if grad.size > 0:
        grad[:] = np.zeros([18, 24 * N])
        for i in range(18):
            grad[i, 6+i+step*24] = -1.0
            grad[i, 6+i+(step+1)*24] = 1.0

    for i in range(18):
        result[i] = x[6+i+(step+1)*24] - x[6+i+step*24] - rotationalSpeedLimit
    return result

def GMSpdLimUp(result, x, grad, step):
    #Note here step should range from 0 to N-2
    gmSpeedLimit = 14  # Limit the speed to 12mm between steps
    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+step*24] = 1.0
        grad[0, 0+(step+1)*24] = -1.0

        grad[1, 1+step*24] = 1.0
        grad[1, 1+(step+1)*24] = -1.0

        grad[2, 2+step*24] = 1.0
        grad[2, 2+(step+1)*24] = -1.0

    result[:] = np.array([
        x[0+step*24] - x[0+(step+1)*24] - gmSpeedLimit,
        x[1+step*24] - x[1+(step+1)*24] - gmSpeedLimit,
        x[2+step*24] - x[2+(step+1)*24] - gmSpeedLimit])
    return result

def GMSpdLimLow(result, x, grad, step):
    #Note here step should range from 0 to N-2
    gmSpeedLimit = 14  # Limit the speed to 12mm between steps
    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+(step+1)*24] = 1.0
        grad[0, 0+step*24] = -1.0

        grad[1, 1+(step+1)*24] = 1.0
        grad[1, 1+step*24] = -1.0

        grad[2, 2+(step+1)*24] = 1.0
        grad[2, 2+step*24] = -1.0

    result[:] = np.array([
        x[0+(step+1)*24] - x[0+step*24] - gmSpeedLimit,
        x[1+(step+1)*24] - x[1+step*24] - gmSpeedLimit,
        x[2+(step+1)*24] - x[2+step*24] - gmSpeedLimit])
    return result

#----------------------------------print end points-----------------------------------------------------------------
def printEndPoints(step):
    #This is RF leg
    resultRF = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[6+24*step]) + L_tibia*(np.cos(x[6+24*step])*np.cos(x[7+24*step])*np.sin(x[8+24*step]) + np.cos(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) +
        L_femur*np.cos(x[6+24*step])*np.cos(x[7+24*step]) + x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(x[6+24*step]) + L_tibia*(np.cos(x[7+24*step])*np.sin(x[6+24*step])*np.sin(x[8+24*step]) + np.sin(x[6+24*step])*np.cos(x[8+24*step])*np.sin(x[7+24*step])) +
        L_femur*np.cos(x[7+24*step])*np.sin(x[6+24*step]) + y_body_half,

        x[2+24*step] + L_femur*np.sin(x[7+24*step]) - L_tibia*(np.cos(x[7+24*step])*np.cos(x[8+24*step]) - np.sin(x[7+24*step])*np.sin(x[8+24*step]))])

    #This is RM leg
    resultRM = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[9+24*step]) + L_tibia*(np.cos(x[9+24*step])*np.cos(x[10+24*step])*np.sin(x[11+24*step]) + np.cos(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) +
        L_femur*np.cos(x[9+24*step])*np.cos(x[10+24*step]) + x_body_mid_half,

        x[1+24*step] +
        L_coxa*np.sin(x[9+24*step]) + L_tibia*(np.cos(x[10+24*step])*np.sin(x[9+24*step])*np.sin(x[11+24*step]) + np.sin(x[9+24*step])*np.cos(x[11+24*step])*np.sin(x[10+24*step])) +
        L_femur*np.cos(x[10+24*step])*np.sin(x[9+24*step]),

        x[2+24*step] + L_femur*np.sin(x[10+24*step]) - L_tibia*(np.cos(x[10+24*step])*np.cos(x[11+24*step]) - np.sin(x[10+24*step])*np.sin(x[11+24*step]))])

    #This is RR leg
    resultRR = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[12+24*step]) + L_tibia*(np.cos(x[12+24*step])*np.cos(x[13+24*step])*np.sin(x[14+24*step]) + np.cos(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) +
        L_femur*np.cos(x[12+24*step])*np.cos(x[13+24*step]) + x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(x[12+24*step]) + L_tibia*(np.cos(x[13+24*step])*np.sin(x[12+24*step])*np.sin(x[14+24*step]) + np.sin(x[12+24*step])*np.cos(x[14+24*step])*np.sin(x[13+24*step])) +
        L_femur*np.cos(x[13+24*step])*np.sin(x[12]) - y_body_half,

        x[2+24*step] + L_femur*np.sin(x[13+24*step]) - L_tibia*(np.cos(x[13+24*step])*np.cos(x[14+24*step]) - np.sin(x[13+24*step])*np.sin(x[14+24*step]))])

    #This is LR leg
    resultLR = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[15+24*step]) + L_tibia*(np.cos(x[15+24*step])*np.cos(x[16+24*step])*np.sin(x[17+24*step]) + np.cos(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) +
        L_femur*np.cos(x[15+24*step])*np.cos(x[16+24*step]) - x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(x[15+24*step]) + L_tibia*(np.cos(x[16+24*step])*np.sin(x[15+24*step])*np.sin(x[17+24*step]) + np.sin(x[15+24*step])*np.cos(x[17+24*step])*np.sin(x[16+24*step])) +
        L_femur*np.cos(x[16+24*step])*np.sin(x[15+24*step]) - y_body_half,

        x[2+24*step] + L_femur*np.sin(x[16+24*step]) - L_tibia*(np.cos(x[16+24*step])*np.cos(x[17+24*step]) - np.sin(x[16+24*step])*np.sin(x[17+24*step]))])

    #This is LM leg
    resultLM = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[18+24*step]) + L_tibia*(np.cos(x[18+24*step])*np.cos(x[19+24*step])*np.sin(x[20+24*step]) + np.cos(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) +
        L_femur*np.cos(x[18+24*step])*np.cos(x[19+24*step]) - x_body_mid_half,

        x[1+24*step] +
        L_coxa*np.sin(x[18+24*step]) + L_tibia*(np.cos(x[19+24*step])*np.sin(x[18+24*step])*np.sin(x[20+24*step]) + np.sin(x[18+24*step])*np.cos(x[20+24*step])*np.sin(x[19+24*step])) +
        L_femur*np.cos(x[19+24*step])*np.sin(x[18+24*step]),

        x[2+24*step] + L_femur*np.sin(x[19+24*step]) - L_tibia*(np.cos(x[19+24*step])*np.cos(x[20+24*step]) - np.sin(x[19+24*step])*np.sin(x[20+24*step]))])

    #This is LF leg
    resultLF = np.array([
        x[0+24*step] +
        L_coxa*np.cos(x[21+24*step]) + L_tibia*(np.cos(x[21+24*step])*np.cos(x[22+24*step])*np.sin(x[23+24*step]) + np.cos(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) +
        L_femur*np.cos(x[21+24*step])*np.cos(x[22+24*step]) - x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(x[21+24*step]) + L_tibia*(np.cos(x[22+24*step])*np.sin(x[21+24*step])*np.sin(x[23+24*step]) + np.sin(x[21+24*step])*np.cos(x[23+24*step])*np.sin(x[22+24*step])) +
        L_femur*np.cos(x[22+24*step])*np.sin(x[21+24*step]) + y_body_half,

        x[2+24*step] + L_femur*np.sin(x[22+24*step]) - L_tibia*(np.cos(x[22+24*step])*np.cos(x[23+24*step]) - np.sin(x[22+24*step])*np.sin(x[23+24*step]))])

    return [resultRF, resultRM, resultRR, resultLR, resultLM, resultLF]

#----------------------------------objective------------------------------------------------------------------------
W_coxa = 5
W_femur = 5
W_tibia = 5
W_diff = 2

def objective(x, grad):
    if grad.size > 0:
        grad[:] = [0] * (24 * N)

        #Derivative of nominal posture part
        for i in range(N):
            grad[6+24*i]  = 2*W_coxa*(x[6 + 24 * i] - RF_coxa_nominal)
            grad[7+24*i]  = 2*W_femur*(x[7 + 24 * i] - femur_nominal)
            grad[8+24*i]  = 2*W_tibia*(x[8 + 24 * i] - tibia_nominal)
            grad[9+24*i]  = 2*W_coxa*(x[9 + 24 * i] - RM_coxa_nominal)
            grad[10+24*i] = 2*W_femur*(x[10 + 24 * i] - femur_nominal)
            grad[11+24*i] = 2*W_tibia*(x[11 + 24 * i] - tibia_nominal)
            grad[12+24*i] = 2*W_coxa*(x[12 + 24 * i] - RR_coxa_nominal)
            grad[13+24*i] = 2*W_femur*(x[13 + 24 * i] - femur_nominal)
            grad[14+24*i] = 2*W_tibia*(x[14 + 24 * i] - tibia_nominal)
            grad[15+24*i] = 2*W_coxa*(x[15 + 24 * i] - LR_coxa_nominal)
            grad[16+24*i] = 2*W_femur*(x[16 + 24 * i] - femur_nominal)
            grad[17+24*i] = 2*W_tibia*(x[17 + 24 * i] - tibia_nominal)
            grad[18+24*i] = 2*W_coxa*(x[18 + 24 * i] - LM_coxa_nominal)
            grad[19+24*i] = 2*W_femur*(x[19 + 24 * i] - femur_nominal)
            grad[20+24*i] = 2*W_tibia*(x[20 + 24 * i] - tibia_nominal)
            grad[21+24*i] = 2*W_coxa*(x[21 + 24 * i] - LF_coxa_nominal)
            grad[22+24*i] = 2*W_femur*(x[22 + 24 * i] - femur_nominal)
            grad[23+24*i] = 2*W_tibia*(x[23 + 24 * i] - tibia_nominal)

        # Derivative of limited difference part
        grad[0] = grad[0] + 2*W_diff*(x[0]-x[0+24])   #same as - 2*W_diff*(x[0+24]-x[0])
        grad[1] = grad[1] + 2*W_diff*(x[1]-x[1+24])
        grad[2] = grad[2] + 2*W_diff*(x[2]-x[2+24])
        grad[3] = grad[3] + 2*W_diff*(x[3]-x[3+24])
        grad[4] = grad[4] + 2*W_diff*(x[4]-x[4+24])
        grad[5] = grad[5] + 2*W_diff*(x[5]-x[5+24])
        grad[6] = grad[6] + 2*W_diff*(x[6]-x[6+24])
        grad[7] = grad[7] + 2*W_diff*(x[7]-x[7+24])
        grad[8] = grad[8] + 2*W_diff*(x[8]-x[8+24])
        grad[9] = grad[9] + 2*W_diff*(x[9]-x[9+24])
        grad[10] = grad[10] + 2*W_diff*(x[10]-x[10+24])
        grad[11] = grad[11] + 2*W_diff*(x[11]-x[11+24])
        grad[12] = grad[12] + 2*W_diff*(x[12]-x[12+24])
        grad[13] = grad[13] + 2*W_diff*(x[13]-x[13+24])
        grad[14] = grad[14] + 2*W_diff*(x[14]-x[14+24])
        grad[15] = grad[15] + 2*W_diff*(x[15]-x[15+24])
        grad[16] = grad[16] + 2*W_diff*(x[16]-x[16+24])
        grad[17] = grad[17] + 2*W_diff*(x[17]-x[17+24])
        grad[18] = grad[18] + 2*W_diff*(x[18]-x[18+24])
        grad[19] = grad[19] + 2*W_diff*(x[19]-x[19+24])
        grad[20] = grad[20] + 2*W_diff*(x[20]-x[20+24])
        grad[21] = grad[21] + 2*W_diff*(x[21]-x[21+24])
        grad[22] = grad[22] + 2*W_diff*(x[22]-x[22+24])
        grad[23] = grad[23] + 2*W_diff*(x[23]-x[23+24])

        grad[0+(N-1)*24] = grad[0+(N-1)*24] - 2*W_diff*(x[0+(N-2)*24]-x[0+(N-1)*24])   # Same as + 2*W_diff*(x[0+(N-1)*24]-x[0+(N-2)*24])
        grad[1+(N-1)*24] = grad[1+(N-1)*24] - 2*W_diff*(x[1+(N-2)*24]-x[1+(N-1)*24])
        grad[2+(N-1)*24] = grad[2+(N-1)*24] - 2*W_diff*(x[2+(N-2)*24]-x[2+(N-1)*24])
        grad[3+(N-1)*24] = grad[3+(N-1)*24] - 2*W_diff*(x[3+(N-2)*24]-x[3+(N-1)*24])
        grad[4+(N-1)*24] = grad[4+(N-1)*24] - 2*W_diff*(x[4+(N-2)*24]-x[4+(N-1)*24])
        grad[5+(N-1)*24] = grad[5+(N-1)*24] - 2*W_diff*(x[5+(N-2)*24]-x[5+(N-1)*24])
        grad[6+(N-1)*24] = grad[6+(N-1)*24] - 2*W_diff*(x[6+(N-2)*24]-x[6+(N-1)*24])
        grad[7+(N-1)*24] = grad[7+(N-1)*24] - 2*W_diff*(x[7+(N-2)*24]-x[7+(N-1)*24])
        grad[8+(N-1)*24] = grad[8+(N-1)*24] - 2*W_diff*(x[8+(N-2)*24]-x[8+(N-1)*24])
        grad[9+(N-1)*24] = grad[9+(N-1)*24] - 2*W_diff*(x[9+(N-2)*24]-x[9+(N-1)*24])
        grad[10+(N-1)*24] = grad[10+(N-1)*24] - 2*W_diff*(x[10+(N-2)*24]-x[10+(N-1)*24])
        grad[11+(N-1)*24] = grad[11+(N-1)*24] - 2*W_diff*(x[11+(N-2)*24]-x[11+(N-1)*24])
        grad[12+(N-1)*24] = grad[12+(N-1)*24] - 2*W_diff*(x[12+(N-2)*24]-x[12+(N-1)*24])
        grad[13+(N-1)*24] = grad[13+(N-1)*24] - 2*W_diff*(x[13+(N-2)*24]-x[13+(N-1)*24])
        grad[14+(N-1)*24] = grad[14+(N-1)*24] - 2*W_diff*(x[14+(N-2)*24]-x[14+(N-1)*24])
        grad[15+(N-1)*24] = grad[15+(N-1)*24] - 2*W_diff*(x[15+(N-2)*24]-x[15+(N-1)*24])
        grad[16+(N-1)*24] = grad[16+(N-1)*24] - 2*W_diff*(x[16+(N-2)*24]-x[16+(N-1)*24])
        grad[17+(N-1)*24] = grad[17+(N-1)*24] - 2*W_diff*(x[17+(N-2)*24]-x[17+(N-1)*24])
        grad[18+(N-1)*24] = grad[18+(N-1)*24] - 2*W_diff*(x[18+(N-2)*24]-x[18+(N-1)*24])
        grad[19+(N-1)*24] = grad[19+(N-1)*24] - 2*W_diff*(x[19+(N-2)*24]-x[19+(N-1)*24])
        grad[20+(N-1)*24] = grad[20+(N-1)*24] - 2*W_diff*(x[20+(N-2)*24]-x[20+(N-1)*24])
        grad[21+(N-1)*24] = grad[21+(N-1)*24] - 2*W_diff*(x[21+(N-2)*24]-x[21+(N-1)*24])
        grad[22+(N-1)*24] = grad[22+(N-1)*24] - 2*W_diff*(x[22+(N-2)*24]-x[22+(N-1)*24])
        grad[23+(N-1)*24] = grad[23+(N-1)*24] - 2*W_diff*(x[23+(N-2)*24]-x[23+(N-1)*24])

        for i in range(1,N-1):
            grad[0+i*24] = grad[0+i*24] - 2*W_diff*(x[0+(i-1)*24]-x[0+i*24]) + 2*W_diff*(x[0+i*24]-x[0+(i+1)*24])
            grad[1+i*24] = grad[1+i*24] - 2*W_diff*(x[1+(i-1)*24]-x[1+i*24]) + 2*W_diff*(x[1+i*24]-x[1+(i+1)*24])
            grad[2+i*24] = grad[2+i*24] - 2*W_diff*(x[2+(i-1)*24]-x[2+i*24]) + 2*W_diff*(x[2+i*24]-x[2+(i+1)*24])
            grad[3+i*24] = grad[3+i*24] - 2*W_diff*(x[3+(i-1)*24]-x[3+i*24]) + 2*W_diff*(x[3+i*24]-x[3+(i+1)*24])
            grad[4+i*24] = grad[4+i*24] - 2*W_diff*(x[4+(i-1)*24]-x[4+i*24]) + 2*W_diff*(x[4+i*24]-x[4+(i+1)*24])
            grad[5+i*24] = grad[5+i*24] - 2*W_diff*(x[5+(i-1)*24]-x[5+i*24]) + 2*W_diff*(x[5+i*24]-x[5+(i+1)*24])
            grad[6+i*24] = grad[6+i*24] - 2*W_diff*(x[6+(i-1)*24]-x[6+i*24]) + 2*W_diff*(x[6+i*24]-x[6+(i+1)*24])
            grad[7+i*24] = grad[7+i*24] - 2*W_diff*(x[7+(i-1)*24]-x[7+i*24]) + 2*W_diff*(x[7+i*24]-x[7+(i+1)*24])
            grad[8+i*24] = grad[8+i*24] - 2*W_diff*(x[8+(i-1)*24]-x[8+i*24]) + 2*W_diff*(x[8+i*24]-x[8+(i+1)*24])
            grad[9+i*24] = grad[9+i*24] - 2*W_diff*(x[9+(i-1)*24]-x[9+i*24]) + 2*W_diff*(x[9+i*24]-x[9+(i+1)*24])
            grad[10+i*24] = grad[10+i*24] - 2*W_diff*(x[10+(i-1)*24]-x[10+i*24]) + 2*W_diff*(x[10+i*24]-x[10+(i+1)*24])
            grad[11+i*24] = grad[11+i*24] - 2*W_diff*(x[11+(i-1)*24]-x[11+i*24]) + 2*W_diff*(x[11+i*24]-x[11+(i+1)*24])
            grad[12+i*24] = grad[12+i*24] - 2*W_diff*(x[12+(i-1)*24]-x[12+i*24]) + 2*W_diff*(x[12+i*24]-x[12+(i+1)*24])
            grad[13+i*24] = grad[13+i*24] - 2*W_diff*(x[13+(i-1)*24]-x[13+i*24]) + 2*W_diff*(x[13+i*24]-x[13+(i+1)*24])
            grad[14+i*24] = grad[14+i*24] - 2*W_diff*(x[14+(i-1)*24]-x[14+i*24]) + 2*W_diff*(x[14+i*24]-x[14+(i+1)*24])
            grad[15+i*24] = grad[15+i*24] - 2*W_diff*(x[15+(i-1)*24]-x[15+i*24]) + 2*W_diff*(x[15+i*24]-x[15+(i+1)*24])
            grad[16+i*24] = grad[16+i*24] - 2*W_diff*(x[16+(i-1)*24]-x[16+i*24]) + 2*W_diff*(x[16+i*24]-x[16+(i+1)*24])
            grad[17+i*24] = grad[17+i*24] - 2*W_diff*(x[17+(i-1)*24]-x[17+i*24]) + 2*W_diff*(x[17+i*24]-x[17+(i+1)*24])
            grad[18+i*24] = grad[18+i*24] - 2*W_diff*(x[18+(i-1)*24]-x[18+i*24]) + 2*W_diff*(x[18+i*24]-x[18+(i+1)*24])
            grad[19+i*24] = grad[19+i*24] - 2*W_diff*(x[19+(i-1)*24]-x[19+i*24]) + 2*W_diff*(x[19+i*24]-x[19+(i+1)*24])
            grad[20+i*24] = grad[20+i*24] - 2*W_diff*(x[20+(i-1)*24]-x[20+i*24]) + 2*W_diff*(x[20+i*24]-x[20+(i+1)*24])
            grad[21+i*24] = grad[21+i*24] - 2*W_diff*(x[21+(i-1)*24]-x[21+i*24]) + 2*W_diff*(x[21+i*24]-x[21+(i+1)*24])
            grad[22+i*24] = grad[22+i*24] - 2*W_diff*(x[22+(i-1)*24]-x[22+i*24]) + 2*W_diff*(x[22+i*24]-x[22+(i+1)*24])
            grad[23+i*24] = grad[23+i*24] - 2*W_diff*(x[23+(i-1)*24]-x[23+i*24]) + 2*W_diff*(x[23+i*24]-x[23+(i+1)*24])

    result = 0
    for i in range(N):
        # Nominal posture part
        result = result + (x[6+24*i]-RF_coxa_nominal)* W_coxa* (x[6+24*i]-RF_coxa_nominal) + \
                          (x[7+24*i]-femur_nominal)*   W_femur*(x[7+24*i]-femur_nominal) + \
                          (x[8+24*i]-tibia_nominal)*   W_tibia*(x[8+24*i]-tibia_nominal) + \
                          (x[9+24*i]-RM_coxa_nominal)* W_coxa* (x[9+24*i]-RM_coxa_nominal) + \
                          (x[10+24*i]-femur_nominal)*  W_femur*(x[10+24*i]-femur_nominal) + \
                          (x[11+24*i]-tibia_nominal)*  W_tibia*(x[11+24*i]-tibia_nominal) + \
                          (x[12+24*i]-RR_coxa_nominal)*W_coxa* (x[12+24*i]-RR_coxa_nominal) + \
                          (x[13+24*i]-femur_nominal)*  W_femur*(x[13+24*i]-femur_nominal) + \
                          (x[14+24*i]-tibia_nominal)*  W_tibia*(x[14+24*i]-tibia_nominal) + \
                          (x[15+24*i]-LR_coxa_nominal)*W_coxa* (x[15+24*i]-LR_coxa_nominal) + \
                          (x[16+24*i]-femur_nominal)*  W_femur*(x[16+24*i]-femur_nominal) + \
                          (x[17+24*i]-tibia_nominal)*  W_tibia*(x[17+24*i]-tibia_nominal) + \
                          (x[18+24*i]-LM_coxa_nominal)*W_coxa* (x[18+24*i]-LM_coxa_nominal) + \
                          (x[19+24*i]-femur_nominal)*  W_femur*(x[19+24*i]-femur_nominal) + \
                          (x[20+24*i]-tibia_nominal)*  W_tibia*(x[20+24*i]-tibia_nominal) + \
                          (x[21+24*i]-LF_coxa_nominal)*W_coxa* (x[21+24*i]-LF_coxa_nominal) + \
                          (x[22+24*i]-femur_nominal)*  W_femur*(x[22+24*i]-femur_nominal) + \
                          (x[23+24*i]-tibia_nominal)*  W_tibia*(x[23+24*i]-tibia_nominal)

    for i in range(N-1):
        # Limited difference part
        result = result + (x[0+i*24]-x[0+(i+1)*24])*W_diff*(x[0+i*24]-x[0+(i+1)*24]) + \
                          (x[1+i*24]-x[1+(i+1)*24])*W_diff*(x[1+i*24]-x[1+(i+1)*24]) + \
                          (x[2+i*24]-x[2+(i+1)*24])*W_diff*(x[2+i*24]-x[2+(i+1)*24]) + \
                          (x[3+i*24]-x[3+(i+1)*24])*W_diff*(x[3+i*24]-x[3+(i+1)*24]) + \
                          (x[4+i*24]-x[4+(i+1)*24])*W_diff*(x[4+i*24]-x[4+(i+1)*24]) + \
                          (x[5+i*24]-x[5+(i+1)*24])*W_diff*(x[5+i*24]-x[5+(i+1)*24]) + \
                          (x[6+i*24]-x[6+(i+1)*24])*W_diff*(x[6+i*24]-x[6+(i+1)*24]) + \
                          (x[7+i*24]-x[7+(i+1)*24])*W_diff*(x[7+i*24]-x[7+(i+1)*24]) + \
                          (x[8+i*24]-x[8+(i+1)*24])*W_diff*(x[8+i*24]-x[8+(i+1)*24]) + \
                          (x[9+i*24]-x[9+(i+1)*24])*W_diff*(x[9+i*24]-x[9+(i+1)*24]) + \
                          (x[10+i*24]-x[10+(i+1)*24])*W_diff*(x[10+i*24]-x[10+(i+1)*24]) + \
                          (x[11+i*24]-x[11+(i+1)*24])*W_diff*(x[11+i*24]-x[11+(i+1)*24]) + \
                          (x[12+i*24]-x[12+(i+1)*24])*W_diff*(x[12+i*24]-x[12+(i+1)*24]) + \
                          (x[13+i*24]-x[13+(i+1)*24])*W_diff*(x[13+i*24]-x[13+(i+1)*24]) + \
                          (x[14+i*24]-x[14+(i+1)*24])*W_diff*(x[14+i*24]-x[14+(i+1)*24]) + \
                          (x[15+i*24]-x[15+(i+1)*24])*W_diff*(x[15+i*24]-x[15+(i+1)*24]) + \
                          (x[16+i*24]-x[16+(i+1)*24])*W_diff*(x[16+i*24]-x[16+(i+1)*24]) + \
                          (x[17+i*24]-x[17+(i+1)*24])*W_diff*(x[17+i*24]-x[17+(i+1)*24]) + \
                          (x[18+i*24]-x[18+(i+1)*24])*W_diff*(x[18+i*24]-x[18+(i+1)*24]) + \
                          (x[19+i*24]-x[19+(i+1)*24])*W_diff*(x[19+i*24]-x[19+(i+1)*24]) + \
                          (x[20+i*24]-x[20+(i+1)*24])*W_diff*(x[20+i*24]-x[20+(i+1)*24]) + \
                          (x[21+i*24]-x[21+(i+1)*24])*W_diff*(x[21+i*24]-x[21+(i+1)*24]) + \
                          (x[22+i*24]-x[22+(i+1)*24])*W_diff*(x[22+i*24]-x[22+(i+1)*24]) + \
                          (x[23+i*24]-x[23+(i+1)*24])*W_diff*(x[23+i*24]-x[23+(i+1)*24])


    return result

#-------------------------------------------------------------------------------------------------------------------
opt = nlopt.opt(nlopt.LD_SLSQP, 24*N)   #Make first 3 GC, second 3 CoM, and the rest motor angles
xUpperBound = 2000.0
xLowerBound = 0.0
yUpperBound = 2000.0
yLowerBound = 0.0
zUpperBound = 400.0
zLowerBound = -10.0
# coxaUpperBound = float('inf')
# coxaLowerBound = -float('inf')
# femurUpperBound = float('inf')
# femurLowerBound = -float('inf')
# tibiaUpperBound = float('inf')
# tibiaLowerBound = -float('inf')

coxaUpperBound_RF = np.pi/2
coxaLowerBound_RF = 0.0
coxaUpperBound_RM = np.pi/4
coxaLowerBound_RM = -np.pi/4
coxaUpperBound_RR = 0.0
coxaLowerBound_RR = -np.pi/2
coxaUpperBound_LR = 3*np.pi/2
coxaLowerBound_LR = np.pi
coxaUpperBound_LM = np.pi+np.pi/4
coxaLowerBound_LM = np.pi/2+np.pi/4
coxaUpperBound_LF = np.pi
coxaLowerBound_LF = np.pi/2

femurUpperBound = np.pi/2
femurLowerBound = -np.pi/2
tibiaUpperBound = np.pi/2
tibiaLowerBound = -np.pi/2

upperBound = [xUpperBound, yUpperBound, zUpperBound,
              xUpperBound, yUpperBound, zUpperBound,
              coxaUpperBound_RF, femurUpperBound, tibiaUpperBound,
              coxaUpperBound_RM, femurUpperBound, tibiaUpperBound,
              coxaUpperBound_RR, femurUpperBound, tibiaUpperBound,
              coxaUpperBound_LR, femurUpperBound, tibiaUpperBound,
              coxaUpperBound_LM, femurUpperBound, tibiaUpperBound,
              coxaUpperBound_LF, femurUpperBound, tibiaUpperBound]*N

lowerBound = [xLowerBound, yLowerBound, zLowerBound,
              xLowerBound, yLowerBound, zLowerBound,
              coxaLowerBound_RF, femurLowerBound, tibiaLowerBound,
              coxaLowerBound_RM, femurLowerBound, tibiaLowerBound,
              coxaLowerBound_RR, femurLowerBound, tibiaLowerBound,
              coxaLowerBound_LR, femurLowerBound, tibiaLowerBound,
              coxaLowerBound_LM, femurLowerBound, tibiaLowerBound,
              coxaLowerBound_LF, femurLowerBound, tibiaLowerBound]*N

opt.set_upper_bounds(upperBound)
opt.set_lower_bounds(lowerBound)

tol = 1e-6
# Constraint 1: One 3-D vector equality constraints speaking c = function(theta, p)
for i in range(N):
    opt.add_equality_mconstraint(lambda result, x, grad, step=i: CoM(result, x, grad, step), np.array([tol, tol, tol]))


# Constraint 2: Three 3-D vector equality constraints speaking r0=[1000, 780, 0], rN=[1000, 926.5, 0], and rS=[900, 852, 0]
opt.add_equality_mconstraint(r0, np.array([tol, tol, tol]))
opt.add_equality_mconstraint(rS, np.array([tol, tol, tol]))
opt.add_equality_mconstraint(rN, np.array([tol, tol, tol]))

#constraint 3: toepoints = specified toepoints on the ground (with limits on each motor angle), note before N/2 only first 3 toes, after N/2 only second 3 toes, at N/2 all 6 toes
#if step < N/2, then leg: right, right, left => leg1, leg3, leg5
#if step > N/2, then leg: left, left, right => leg2, leg4, leg6
#hexagon1: [1250, 1180, 0], [1500, 780, 0], [1250, 380, 0], [750, 380, 0], [500, 780, 0], [750, 1180, 0]
#hexagon2: [1250, 1326.5, 0], [1500, 924, 0], [1250, 526.5, 0], [750, 526.5, 0], [500, 924, 0], [750, 1326.5, 0]
#Note those end points are on the ground
px_RF_1 = 1250.0
py_RF_1 = 1180.0
pz_RF_1 = 0.0

px_RM_2 = 1500.0
py_RM_2 = 924.0
pz_RM_2 = 0.0

px_RR_1 = 1250.0
py_RR_1 = 380.0
pz_RR_1 = 0.0

px_LR_2 = 750.0
py_LR_2 = 526.5
pz_LR_2 = 0.0

px_LM_1 = 500.0
py_LM_1 = 780.0
pz_LM_1 = 0.0

px_LF_2 = 750.0
py_LF_2 = 1326.5
pz_LF_2 = 0.0

for i in range(N/2):
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_1, py=py_RF_1, pz=pz_RF_1: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_1, py=py_RR_1, pz=pz_RR_1: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_1, py=py_LM_1, pz=pz_LM_1: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

    #z>0
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg2_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg4_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg6_end_point_Z(x, grad, step, pz), tol)

i = N/2
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_1, py=py_RF_1, pz=pz_RF_1: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_1, py=py_RR_1, pz=pz_RR_1: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_1, py=py_LM_1, pz=pz_LM_1: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_2, py=py_RM_2, pz=pz_RM_2: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_2, py=py_LR_2, pz=pz_LR_2: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_2, py=py_LF_2, pz=pz_LF_2: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

for i in range(N/2+1, N):
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_2, py=py_RM_2, pz=pz_RM_2: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_2, py=py_LR_2, pz=pz_LR_2: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_2, py=py_LF_2, pz=pz_LF_2: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

    #z>0
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg1_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg3_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=5.0: leg5_end_point_Z(x, grad, step, pz), tol)

#constraint 4: let's specify a triangular trajectory for the limbs on the air
liftHeight = 30.0  #Specify the liftheight to 50mm
i = N/4
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg2_end_point_Z(x, grad, step, pz), tol)
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg4_end_point_Z(x, grad, step, pz), tol)
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg6_end_point_Z(x, grad, step, pz), tol)
i = 3*N/4
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg1_end_point_Z(x, grad, step, pz), tol)
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg3_end_point_Z(x, grad, step, pz), tol)
opt.add_equality_constraint(lambda x, grad, step=i, pz=liftHeight: leg5_end_point_Z(x, grad, step, pz), tol)

#constraint 5: CoM within support polygon
for i in range(N):
    opt.add_inequality_mconstraint(lambda result, x, grad, step=i: supportHexagonConstraint(result, x, grad, step), np.array([tol, tol, tol]))

#constraint 6: each joint angle and geometri center should have limited amount of moving
for i in range(N-1):
    opt.add_inequality_mconstraint(lambda result, x, grad, step=i: rotSpdLimUp(result, x, grad, step), np.array([tol]*18))
    opt.add_inequality_mconstraint(lambda result, x, grad, step=i: rotSpdLimLow(result, x, grad, step), np.array([tol]*18))

    opt.add_inequality_mconstraint(lambda result, x, grad, step=i: GMSpdLimUp(result, x, grad, step), np.array([tol]*3))
    opt.add_inequality_mconstraint(lambda result, x, grad, step=i: GMSpdLimLow(result, x, grad, step), np.array([tol]*3))

#Objective function
opt.set_min_objective(objective)
opt.set_xtol_rel(1e-4)
#opt.set_ftol_rel(tol)  #This makes it stop too early!
maxtime = 1200   #stop at 20 min
opt.set_maxtime(maxtime)
ini = [0.0]*(24*N)

RF_coxa_nominal = np.pi/4
RM_coxa_nominal = 0.0
RR_coxa_nominal = -np.pi/4
LR_coxa_nominal = 5*np.pi/4
LM_coxa_nominal = np.pi
LF_coxa_nominal = 3*np.pi/4
femur_nominal = np.pi/6
tibia_nominal = -np.pi/8
for i in range(N):
    # RF
    ini[6 + 24 * i] = RF_coxa_nominal
    ini[7 + 24 * i] = femur_nominal
    ini[8 + 24 * i] = tibia_nominal

    # RM
    ini[9 + 24 * i] = RM_coxa_nominal
    ini[10 + 24 * i] = femur_nominal
    ini[11 + 24 * i] = tibia_nominal

    # RR
    ini[12 + 24 * i] = RR_coxa_nominal
    ini[13 + 24 * i] = femur_nominal
    ini[14 + 24 * i] = tibia_nominal

    # LR
    ini[15 + 24 * i] = LR_coxa_nominal
    ini[16 + 24 * i] = femur_nominal
    ini[17 + 24 * i] = tibia_nominal

    # LM
    ini[18 + 24 * i] = LM_coxa_nominal
    ini[19 + 24 * i] = femur_nominal
    ini[20 + 24 * i] = tibia_nominal

    # LF
    ini[21 + 24 * i] = LF_coxa_nominal
    ini[22 + 24 * i] = femur_nominal
    ini[23 + 24 * i] = tibia_nominal

for i in range(N/2):
    ini[0 + 24 * i] = r0x
    ini[1 + 24 * i] = r0y
    ini[2 + 24 * i] = r0z
    ini[3 + 24 * i] = r0x
    ini[4 + 24 * i] = r0y
    ini[5 + 24 * i] = r0z

for i in range(N/2, N):
    ini[0 + 24 * i] = rNx
    ini[1 + 24 * i] = rNy
    ini[2 + 24 * i] = rNz
    ini[3 + 24 * i] = rNx
    ini[4 + 24 * i] = rNy
    ini[5 + 24 * i] = rNz

ini[0] = r0x
ini[1] = r0y
ini[2] = r0z

ini[24*(N/2-1)+0] = rSx
ini[24*(N/2-1)+1] = rSy
ini[24*(N/2-1)+2] = rSz

ini[24*(N-1)+0] = rNx
ini[24*(N-1)+1] = rNy
ini[24*(N-1)+2] = rNz

#Record time
t1 = time.time()
x = opt.optimize(ini)
t2 = time.time()
print("Takes %.2f seconds"%(t2-t1))
minf = opt.last_optimum_value()

for i in range(N):
    print("Step", i)
    print("GM and CM %.2f   %.2f   %.2f   %.2f   %.2f   %.2f"%(x[0+i*24], x[1+i*24], x[2+i*24], x[3+i*24], x[4+i*24], x[5+i*24]))
    print("---------------------------------------------------")
for i in range(N):
    print("Step", i)
    print("RF motor angles %.2f   %.2f   %.2f" % (x[6 + i * 24]*180/3.14,  x[7 + i * 24]*180/3.14,  x[8 + i * 24]*180/3.14))
    print("RM motor angles %.2f   %.2f   %.2f" % (x[9 + i * 24]*180/3.14,  x[10 + i * 24]*180/3.14, x[11 + i * 24]*180/3.14))
    print("RR motor angles %.2f   %.2f   %.2f" % (x[12 + i * 24]*180/3.14, x[13 + i * 24]*180/3.14, x[14 + i * 24]*180/3.14))
    print("LR motor angles %.2f   %.2f   %.2f" % (x[15 + i * 24]*180/3.14, x[16 + i * 24]*180/3.14, x[17 + i * 24]*180/3.14))
    print("LM motor angles %.2f   %.2f   %.2f" % (x[18 + i * 24]*180/3.14, x[19 + i * 24]*180/3.14, x[20 + i * 24]*180/3.14))
    print("LF motor angles %.2f   %.2f   %.2f" % (x[21 + i * 24]*180/3.14, x[22 + i * 24]*180/3.14, x[23 + i * 24]*180/3.14))
    print("---------------------------------------------------")

print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())

#Print all endpoints
for i in range(N):
    results = printEndPoints(i)
    print("Step", i)
    print("RF Leg", results[0])
    print("RM Leg", results[1])
    print("RR Leg", results[2])
    print("LR Leg", results[3])
    print("LM Leg", results[4])
    print("LF Leg", results[5])
    print("---------------------------------------------------")