#!/usr/bin/python

import gurobipy as go
import numpy as np

#1 needs go.GRB.something
#2 Open IDE to check grammar
#3 always row vector
#4 needs model.update()

def addNormConstraint1(mod, x1, y1, z1, x2, y2, z2, th):
	mod.addQConstr(x1*x1+x2*x2-2*x1*x2 + y1*y1+y2*y2-2*y1*y2 + z1*z1+z2*z2-2*z1*z2 <= (th*th))

def addNormConstraint2(mod, x1, y1, z1, x2, y2, z2, cx, cy, cz, th):
	mod.addQConstr(x1*x1+x2*x2-2*x1*x2 + y1*y1+y2*y2-2*y1*y2 + z1*z1+z2*z2-2*z1*z2 + cx*cx + cy*cy + cz*cz + \
				   2*x2*cx + 2*y2*cy + 2*z2+cz - 2*x1*cx - 2*y1*cy - 2*z1*cz <= (th*th))

def addPolytopeConstraint(mod, A, x, b):
	dim1 = np.size(b)
	dim2 = np.size(x)
	for i in range(dim1):
		ss = (go.quicksum(A[i,j]*x[j] for j in range(dim2)) <= b[i])
		#print ss
		mod.addConstr(ss)

def setQuadraticObjective(mod, x, g, Q):
	dim = np.size(x)
	ss = go.quicksum((x[i]-g[i])*Q[i,j]*(x[j]-g[j]) for i in range(dim) for j in range(dim))
	#not checking if g and Q have correspoinding size. Make sure that they have!
	mod.setObjective(ss, go.GRB.MINIMIZE)



# Create a new model
m = go.Model("setpOnStairs")

# Step variables, first index is count of legs, second index is count of step
#Step 0: starting point -> constants
# g = np.array([200.0, 1200.0, 0.0, 200.0, 1500.0, 0.0, 700.0, 1350.0, 0.0])   #Goal position
#foot 1
p1x_0 = 750.0
p1y_0 = 0.0
p1z_0 = 0.0

#foot 2
p2x_0 = 1250.0
p2y_0 = 0.0
p2z_0 = 0.0

#foot 3
p3x_0 = 1500.0
p3y_0 = -400.0
p3z_0 = 0.0

#foot 4
p4x_0 = 1250.0
p4y_0 = -800.0
p4z_0 = 0.0

#foot 5
p5x_0 = 750.0
p5y_0 = -800.0
p5z_0 = 0.0

#foot 6
p6x_0 = 500.0
p6y_0 = -400.0
p6z_0 = 0.0

#####################Step 1###############################
# foot1
p1x_1 = m.addVar()
p1y_1 = m.addVar()
p1z_1 = m.addVar()

# foot2
p2x_1 = m.addVar()
p2y_1 = m.addVar()
p2z_1 = m.addVar()

# foot3
p3x_1 = m.addVar()
p3y_1 = m.addVar()
p3z_1 = m.addVar()

# foot4
p4x_1 = m.addVar()
p4y_1 = m.addVar()
p4z_1 = m.addVar()

# foot5
p5x_1 = m.addVar()
p5y_1 = m.addVar()
p5z_1 = m.addVar()

# foot6
p6x_1 = m.addVar()
p6y_1 = m.addVar()
p6z_1 = m.addVar()

#####################Step 2###############################
# foot1
p1x_2 = m.addVar()
p1y_2 = m.addVar()
p1z_2 = m.addVar()

# foot2
p2x_2 = m.addVar()
p2y_2 = m.addVar()
p2z_2 = m.addVar()

# foot3
p3x_2 = m.addVar()
p3y_2 = m.addVar()
p3z_2 = m.addVar()

# foot4
p4x_2 = m.addVar()
p4y_2 = m.addVar()
p4z_2 = m.addVar()

# foot5
p5x_2 = m.addVar()
p5y_2 = m.addVar()
p5z_2 = m.addVar()

# foot6
p6x_2 = m.addVar()
p6y_2 = m.addVar()
p6z_2 = m.addVar()

#####################Step 3###############################
# foot1
p1x_3 = m.addVar()
p1y_3 = m.addVar()
p1z_3 = m.addVar()

# foot2
p2x_3 = m.addVar()
p2y_3 = m.addVar()
p2z_3 = m.addVar()

# foot3
p3x_3 = m.addVar()
p3y_3 = m.addVar()
p3z_3 = m.addVar()

# foot4
p4x_3 = m.addVar()
p4y_3 = m.addVar()
p4z_3 = m.addVar()

# foot5
p5x_3 = m.addVar()
p5y_3 = m.addVar()
p5z_3 = m.addVar()

# foot6
p6x_3 = m.addVar()
p6y_3 = m.addVar()
p6z_3 = m.addVar()

# Array of steps on feasible area
# First index is leg count, second index is step, third index is feasible area
H1_11 = m.addVar(vtype=go.GRB.BINARY)
H1_12 = m.addVar(vtype=go.GRB.BINARY)
H1_13 = m.addVar(vtype=go.GRB.BINARY)

H1_21 = m.addVar(vtype=go.GRB.BINARY)
H1_22 = m.addVar(vtype=go.GRB.BINARY)
H1_23 = m.addVar(vtype=go.GRB.BINARY)

H1_31 = m.addVar(vtype=go.GRB.BINARY)
H1_32 = m.addVar(vtype=go.GRB.BINARY)
H1_33 = m.addVar(vtype=go.GRB.BINARY)

H2_11 = m.addVar(vtype=go.GRB.BINARY)
H2_12 = m.addVar(vtype=go.GRB.BINARY)
H2_13 = m.addVar(vtype=go.GRB.BINARY)

H2_21 = m.addVar(vtype=go.GRB.BINARY)
H2_22 = m.addVar(vtype=go.GRB.BINARY)
H2_23 = m.addVar(vtype=go.GRB.BINARY)

H2_31 = m.addVar(vtype=go.GRB.BINARY)
H2_32 = m.addVar(vtype=go.GRB.BINARY)
H2_33 = m.addVar(vtype=go.GRB.BINARY)

H3_11 = m.addVar(vtype=go.GRB.BINARY)
H3_12 = m.addVar(vtype=go.GRB.BINARY)
H3_13 = m.addVar(vtype=go.GRB.BINARY)

H3_21 = m.addVar(vtype=go.GRB.BINARY)
H3_22 = m.addVar(vtype=go.GRB.BINARY)
H3_23 = m.addVar(vtype=go.GRB.BINARY)

H3_31 = m.addVar(vtype=go.GRB.BINARY)
H3_32 = m.addVar(vtype=go.GRB.BINARY)
H3_33 = m.addVar(vtype=go.GRB.BINARY)

H4_11 = m.addVar(vtype=go.GRB.BINARY)
H4_12 = m.addVar(vtype=go.GRB.BINARY)
H4_13 = m.addVar(vtype=go.GRB.BINARY)

H4_21 = m.addVar(vtype=go.GRB.BINARY)
H4_22 = m.addVar(vtype=go.GRB.BINARY)
H4_23 = m.addVar(vtype=go.GRB.BINARY)

H4_31 = m.addVar(vtype=go.GRB.BINARY)
H4_32 = m.addVar(vtype=go.GRB.BINARY)
H4_33 = m.addVar(vtype=go.GRB.BINARY)

H5_11 = m.addVar(vtype=go.GRB.BINARY)
H5_12 = m.addVar(vtype=go.GRB.BINARY)
H5_13 = m.addVar(vtype=go.GRB.BINARY)

H5_21 = m.addVar(vtype=go.GRB.BINARY)
H5_22 = m.addVar(vtype=go.GRB.BINARY)
H5_23 = m.addVar(vtype=go.GRB.BINARY)

H5_31 = m.addVar(vtype=go.GRB.BINARY)
H5_32 = m.addVar(vtype=go.GRB.BINARY)
H5_33 = m.addVar(vtype=go.GRB.BINARY)

H6_11 = m.addVar(vtype=go.GRB.BINARY)
H6_12 = m.addVar(vtype=go.GRB.BINARY)
H6_13 = m.addVar(vtype=go.GRB.BINARY)

H6_21 = m.addVar(vtype=go.GRB.BINARY)
H6_22 = m.addVar(vtype=go.GRB.BINARY)
H6_23 = m.addVar(vtype=go.GRB.BINARY)

H6_31 = m.addVar(vtype=go.GRB.BINARY)
H6_32 = m.addVar(vtype=go.GRB.BINARY)
H6_33 = m.addVar(vtype=go.GRB.BINARY)

#Update the model to incorporate new variables
m.update()

#Bound on step size, 9 quadratic constraints in total
#th_stepSize = 230
th_stepSize = 850   #For 3 steps, need 450 stepsize
#Step 0-1
addNormConstraint1(m, p1x_0, p1y_0, p1z_0, p1x_1, p1y_1, p1z_1, th_stepSize)
addNormConstraint1(m, p2x_0, p2y_0, p2z_0, p2x_1, p2y_1, p2z_1, th_stepSize)
addNormConstraint1(m, p3x_0, p3y_0, p3z_0, p3x_1, p3y_1, p3z_1, th_stepSize)
addNormConstraint1(m, p4x_0, p4y_0, p4z_0, p4x_1, p4y_1, p4z_1, th_stepSize)
addNormConstraint1(m, p5x_0, p5y_0, p5z_0, p5x_1, p5y_1, p5z_1, th_stepSize)
addNormConstraint1(m, p6x_0, p6y_0, p6z_0, p6x_1, p6y_1, p6z_1, th_stepSize)

#Step 1-2
addNormConstraint1(m, p1x_1, p1y_1, p1z_1, p1x_2, p1y_2, p1z_2, th_stepSize)
addNormConstraint1(m, p2x_1, p2y_1, p2z_1, p2x_2, p2y_2, p2z_2, th_stepSize)
addNormConstraint1(m, p3x_1, p3y_1, p3z_1, p3x_2, p3y_2, p3z_2, th_stepSize)
addNormConstraint1(m, p4x_1, p4y_1, p4z_1, p4x_2, p4y_2, p4z_2, th_stepSize)
addNormConstraint1(m, p5x_1, p5y_1, p5z_1, p5x_2, p5y_2, p5z_2, th_stepSize)
addNormConstraint1(m, p6x_1, p6y_1, p6z_1, p6x_2, p6y_2, p6z_2, th_stepSize)

#Step 2-3
addNormConstraint1(m, p1x_2, p1y_2, p1z_2, p1x_3, p1y_3, p1z_3, th_stepSize)
addNormConstraint1(m, p2x_2, p2y_2, p2z_2, p2x_3, p2y_3, p2z_3, th_stepSize)
addNormConstraint1(m, p3x_2, p3y_2, p3z_2, p3x_3, p3y_3, p3z_3, th_stepSize)
addNormConstraint1(m, p4x_2, p4y_2, p4z_2, p4x_3, p4y_3, p4z_3, th_stepSize)
addNormConstraint1(m, p5x_2, p5y_2, p5z_2, p5x_3, p5y_3, p5z_3, th_stepSize)
addNormConstraint1(m, p6x_2, p6y_2, p6z_2, p6x_3, p6y_3, p6z_3, th_stepSize)


#Reachability
#Rotation matrices
rotMatrix1 = np.array([[np.cos(60.0/180.0*np.pi),  -np.sin(60.0/180.0*np.pi),  0.0],[np.sin(60.0/180.0*np.pi),  np.cos(60.0/180.0*np.pi),  0.0],[0.0, 0.0, 1.0]])
rotMatrix2 = np.array([[np.cos(120.0/180.0*np.pi), -np.sin(120.0/180.0*np.pi), 0.0],[np.sin(120.0/180.0*np.pi), np.cos(120.0/180.0*np.pi), 0.0],[0.0, 0.0, 1.0]])
rotMatrix3 = np.array([[np.cos(180.0/180.0*np.pi), -np.sin(180.0/180.0*np.pi), 0.0],[np.sin(180.0/180.0*np.pi), np.cos(180.0/180.0*np.pi), 0.0],[0.0, 0.0, 1.0]])
rotMatrix4 = np.array([[np.cos(240.0/180.0*np.pi), -np.sin(240.0/180.0*np.pi), 0.0],[np.sin(240.0/180.0*np.pi), np.cos(240.0/180.0*np.pi), 0.0],[0.0, 0.0, 1.0]])
rotMatrix5 = np.array([[np.cos(300.0/180.0*np.pi), -np.sin(300.0/180.0*np.pi), 0.0],[np.sin(300.0/180.0*np.pi), np.cos(300.0/180.0*np.pi), 0.0],[0.0, 0.0, 1.0]])

#First index: foot count, second index: direction, third index: step
################Offset vector#######################################
r11 = np.array([0,0,0])
r12 = np.array([600.0,0,0])   #nominal 600
r21 = np.array([0,0,0])
r22 = rotMatrix1.dot(r12)
r31 = np.array([0,0,0])
r32 = rotMatrix2.dot(r12)
r41 = np.array([0,0,0])
r42 = rotMatrix3.dot(r12)
r51 = np.array([0,0,0])
r52 = rotMatrix4.dot(r12)
r61 = np.array([0,0,0])
r62 = rotMatrix5.dot(r12)

radi_1 = 900   #Nominal 600
radi_2 = 580   #Nominal 275

################################### Try to do: get constraint to see what constraint does it set#############################################
#step 0 - 1
addNormConstraint1(m, p2x_1, p2y_1, p2z_1,  p1x_1+r11[0], p1y_1+r11[1], p1z_1+r11[2], radi_1)
#addNormConstraint1(m, p2x_1, p2y_1, p2z_1,  p1x_1+r12[0], p1y_1+r12[1], p1z_1+r12[2], radi_1)   #This is used to test quadratic expression
addNormConstraint2(m, p2x_1, p2y_1, p2z_1,  p1x_1,        p1y_1,        p1z_1,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p3x_1, p3y_1, p3z_1,  p2x_1+r21[0], p2y_1+r21[1], p2z_1+r21[2], radi_1)
addNormConstraint2(m, p3x_1, p3y_1, p3z_1,  p2x_1,        p2y_1,        p2z_1,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p4x_1, p4y_1, p4z_1,  p3x_1+r31[0], p3y_1+r31[1], p3z_1+r31[2], radi_1)
addNormConstraint2(m, p4x_1, p4y_1, p4z_1,  p3x_1,        p3y_1,        p3z_1,        r32[0], r32[1], r32[2], radi_2)

addNormConstraint1(m, p5x_1, p5y_1, p5z_1,  p4x_1+r41[0], p4y_1+r41[1], p4z_1+r41[2], radi_1)
addNormConstraint2(m, p5x_1, p5y_1, p5z_1,  p4x_1,        p4y_1,        p4z_1,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p6x_1, p6y_1, p6z_1,  p5x_1+r51[0], p5y_1+r51[1], p5z_1+r51[2], radi_1)
addNormConstraint2(m, p6x_1, p6y_1, p6z_1,  p5x_1,        p5y_1,        p5z_1,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p1x_1, p1y_1, p1z_1,  p6x_1+r61[0], p6y_1+r61[1], p6z_1+r61[2], radi_1)
addNormConstraint2(m, p1x_1, p1y_1, p1z_1,  p6x_1,        p6y_1,        p6z_1,        r32[0], r32[1], r32[2], radi_2)

#Step 1 - 2
addNormConstraint1(m, p2x_2, p2y_2, p2z_2,  p1x_2+r11[0], p1y_2+r11[1], p1z_2+r11[2], radi_1)
addNormConstraint2(m, p2x_2, p2y_2, p2z_2,  p1x_2,        p1y_2,        p1z_2,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p3x_2, p3y_2, p3z_2,  p2x_2+r21[0], p2y_2+r21[1], p2z_2+r21[2], radi_1)
addNormConstraint2(m, p3x_2, p3y_2, p3z_2,  p2x_2,        p2y_2,        p2z_2,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p4x_2, p4y_2, p4z_2,  p3x_2+r31[0], p3y_2+r31[1], p3z_2+r31[2], radi_1)
addNormConstraint2(m, p4x_2, p4y_2, p4z_2,  p3x_2,        p3y_2,        p3z_2,        r32[0], r32[1], r32[2], radi_2)

addNormConstraint1(m, p5x_2, p5y_2, p5z_2,  p4x_2+r41[0], p4y_2+r41[1], p4z_2+r41[2], radi_1)
addNormConstraint2(m, p5x_2, p5y_2, p5z_2,  p4x_2,        p4y_2,        p4z_2,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p6x_2, p6y_2, p6z_2,  p5x_2+r51[0], p5y_2+r51[1], p5z_2+r51[2], radi_1)
addNormConstraint2(m, p6x_2, p6y_2, p6z_2,  p5x_2,        p5y_2,        p5z_2,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p1x_2, p1y_2, p1z_2,  p6x_2+r61[0], p6y_2+r61[1], p6z_2+r61[2], radi_1)
addNormConstraint2(m, p1x_2, p1y_2, p1z_2,  p6x_2,        p6y_2,        p6z_2,        r32[0], r32[1], r32[2], radi_2)

#Step 2 - 3
addNormConstraint1(m, p2x_3, p2y_3, p2z_3,  p1x_3+r11[0], p1y_3+r11[1], p1z_3+r11[2], radi_1)
addNormConstraint2(m, p2x_3, p2y_3, p2z_3,  p1x_3,        p1y_3,        p1z_3,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p3x_3, p3y_3, p3z_3,  p2x_3+r21[0], p2y_3+r21[1], p2z_3+r21[2], radi_1)
addNormConstraint2(m, p3x_3, p3y_3, p3z_3,  p2x_3,        p2y_3,        p2z_3,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p4x_3, p4y_3, p4z_3,  p3x_3+r31[0], p3y_3+r31[1], p3z_3+r31[2], radi_1)
addNormConstraint2(m, p4x_3, p4y_3, p4z_3,  p3x_3,        p3y_3,        p3z_3,        r32[0], r32[1], r32[2], radi_2)

addNormConstraint1(m, p5x_3, p5y_3, p5z_3,  p4x_3+r41[0], p4y_3+r41[1], p4z_3+r41[2], radi_1)
addNormConstraint2(m, p5x_3, p5y_3, p5z_3,  p4x_3,        p4y_3,        p4z_3,        r12[0], r12[1], r12[2], radi_2)

addNormConstraint1(m, p6x_3, p6y_3, p6z_3,  p5x_3+r51[0], p5y_3+r51[1], p5z_3+r51[2], radi_1)
addNormConstraint2(m, p6x_3, p6y_3, p6z_3,  p5x_3,        p5y_3,        p5z_3,        r22[0], r22[1], r22[2], radi_2)

addNormConstraint1(m, p1x_3, p1y_3, p1z_3,  p6x_3+r61[0], p6y_3+r61[1], p6z_3+r61[2], radi_1)
addNormConstraint2(m, p1x_3, p1y_3, p1z_3,  p6x_3,        p6y_3,        p6z_3,        r32[0], r32[1], r32[2], radi_2)

# th_reachability = 500
# #Step 1 - 2
# addNormConstraint(m, p1x_2, p1y_2, p1z_2, p1x_1, p1y_1, p1z_1, th_reachability)
# addNormConstraint(m, p1x_2, p1y_2, p1z_2, p2x_1, p2y_1, p2z_1, th_reachability)
#
# addNormConstraint(m, p2x_2, p2y_2, p2z_2, p2x_1, p2y_1, p2z_1, th_reachability)
# addNormConstraint(m, p2x_2, p2y_2, p2z_2, p3x_1, p3y_1, p3z_1, th_reachability)
#
# addNormConstraint(m, p3x_2, p3y_2, p3z_2, p1x_1, p1y_1, p1z_1, th_reachability)
# addNormConstraint(m, p3x_2, p3y_2, p3z_2, p3x_1, p3y_1, p3z_1, th_reachability)
#
# #Step 2 - 3
# addNormConstraint(m, p1x_3, p1y_3, p1z_3, p1x_2, p1y_2, p1z_2, th_reachability)
# addNormConstraint(m, p1x_3, p1y_3, p1z_3, p2x_2, p2y_2, p2z_2, th_reachability)
#
# addNormConstraint(m, p2x_3, p2y_3, p2z_3, p2x_2, p2y_2, p2z_2, th_reachability)
# addNormConstraint(m, p2x_3, p2y_3, p2z_3, p3x_2, p3y_2, p3z_2, th_reachability)
#
# addNormConstraint(m, p3x_3, p3y_3, p3z_3, p1x_2, p1y_2, p1z_2, th_reachability)
# addNormConstraint(m, p3x_3, p3y_3, p3z_3, p3x_2, p3y_2, p3z_2, th_reachability)

#Summation of a row of H (total regions in one step) is 1
# First index is leg count, second index is step, third index is feasible area
m.addConstr(go.quicksum([H1_11, H1_12, H1_13]) == 1)
m.addConstr(go.quicksum([H1_21, H1_22, H1_23]) == 1)
m.addConstr(go.quicksum([H1_31, H1_32, H1_33]) == 1)

m.addConstr(go.quicksum([H2_11, H2_12, H2_13]) == 1)
m.addConstr(go.quicksum([H2_21, H2_22, H2_23]) == 1)
m.addConstr(go.quicksum([H2_31, H2_32, H2_33]) == 1)

m.addConstr(go.quicksum([H3_11, H3_12, H3_13]) == 1)
m.addConstr(go.quicksum([H3_21, H3_22, H3_23]) == 1)
m.addConstr(go.quicksum([H3_31, H3_32, H3_33]) == 1)

m.addConstr(go.quicksum([H4_11, H4_12, H4_13]) == 1)
m.addConstr(go.quicksum([H4_21, H4_22, H4_23]) == 1)
m.addConstr(go.quicksum([H4_31, H4_32, H4_33]) == 1)

m.addConstr(go.quicksum([H5_11, H5_12, H5_13]) == 1)
m.addConstr(go.quicksum([H5_21, H5_22, H5_23]) == 1)
m.addConstr(go.quicksum([H5_31, H5_32, H5_33]) == 1)

m.addConstr(go.quicksum([H6_11, H6_12, H6_13]) == 1)
m.addConstr(go.quicksum([H6_21, H6_22, H6_23]) == 1)
m.addConstr(go.quicksum([H6_31, H6_32, H6_33]) == 1)

#Length of a feasible area
length = 220.0
#Witdh of a feasible area
width = 2000.0
#Height of a feasible area
height = 1.0
#Step height, the height of each step
step_height = 130
#Distance between two feasible areas
distance = 0.0
#Expression of a polytope
#Remember: always row vector
A = np.array([[ 1.0,  0.0,  0.0], \
			  [-1.0,  0.0,  0.0], \
			  [ 0.0,  1.0,  0.0], \
			  [ 0.0, -1.0,  0.0], \
			  [ 0.0,  0.0,  1.0], \
			  [ 0.0,  0.0, -1.0]])
b1 = np.array([width,  0.0,    \
			   length, 5000.0, \
			   height, 0.0])
b2 = np.array([width, 0.0, \
			   length*2+distance,  -(length+distance), \
			   height+step_height, -step_height])
b3 = np.array([width, 0.0, \
			   5000.0, -(length*2+distance*2), \
			   height+2*step_height, -2*step_height])
#Big-M constants
M1 = 5000.0
M2 = 5000.0
M3 = 5000.0

#Constraints are:
# First index is leg count, second index is step, third index is feasible area
######################################## Step 1 ##########################################
#Feasible area 1

# foot1 step1 feasible_area1
xx = np.array([p1x_1, p1y_1, p1z_1])
yy = np.array([M1*(1.0-H1_11), M1*(1.0-H1_11), M1*(1.0-H1_11), M1*(1.0-H1_11), M1*(1.0-H1_11), M1*(1.0-H1_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area1
xx = np.array([p2x_1, p2y_1, p2z_1])
yy = np.array([M1*(1.0-H2_11), M1*(1.0-H2_11), M1*(1.0-H2_11), M1*(1.0-H2_11), M1*(1.0-H2_11), M1*(1.0-H2_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area1
xx = np.array([p3x_1, p3y_1, p3z_1])
yy = np.array([M1*(1.0-H3_11), M1*(1.0-H3_11), M1*(1.0-H3_11), M1*(1.0-H3_11), M1*(1.0-H3_11), M1*(1.0-H3_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area1
xx = np.array([p4x_1, p4y_1, p4z_1])
yy = np.array([M1*(1.0-H4_11), M1*(1.0-H4_11), M1*(1.0-H4_11), M1*(1.0-H4_11), M1*(1.0-H4_11), M1*(1.0-H4_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area1
xx = np.array([p5x_1, p5y_1, p5z_1])
yy = np.array([M1*(1.0-H5_11), M1*(1.0-H5_11), M1*(1.0-H5_11), M1*(1.0-H5_11), M1*(1.0-H5_11), M1*(1.0-H5_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area1
xx = np.array([p6x_1, p6y_1, p6z_1])
yy = np.array([M1*(1.0-H6_11), M1*(1.0-H6_11), M1*(1.0-H6_11), M1*(1.0-H6_11), M1*(1.0-H6_11), M1*(1.0-H6_11)])
b = np.add(b1,yy)
addPolytopeConstraint(m, A, xx, b)

#Feasible area 2

# foot1 step1 feasible_area2
xx = np.array([p1x_1, p1y_1, p1z_1])
yy = np.array([M2*(1.0-H1_12), M2*(1.0-H1_12), M2*(1.0-H1_12), M2*(1.0-H1_12), M2*(1.0-H1_12), M2*(1.0-H1_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area2
xx = np.array([p2x_1, p2y_1, p2z_1])
yy = np.array([M2*(1.0-H2_12), M2*(1.0-H2_12), M2*(1.0-H2_12), M2*(1.0-H2_12), M2*(1.0-H2_12), M2*(1.0-H2_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area2
xx = np.array([p3x_1, p3y_1, p3z_1])
yy = np.array([M2*(1.0-H3_12), M2*(1.0-H3_12), M2*(1.0-H3_12), M2*(1.0-H3_12), M2*(1.0-H3_12), M2*(1.0-H3_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area2
xx = np.array([p4x_1, p4y_1, p4z_1])
yy = np.array([M2*(1.0-H4_12), M2*(1.0-H4_12), M2*(1.0-H4_12), M2*(1.0-H4_12), M2*(1.0-H4_12), M2*(1.0-H4_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area2
xx = np.array([p5x_1, p5y_1, p5z_1])
yy = np.array([M2*(1.0-H5_12), M2*(1.0-H5_12), M2*(1.0-H5_12), M2*(1.0-H5_12), M2*(1.0-H5_12), M2*(1.0-H5_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area2
xx = np.array([p6x_1, p6y_1, p6z_1])
yy = np.array([M2*(1.0-H6_12), M2*(1.0-H6_12), M2*(1.0-H6_12), M2*(1.0-H6_12), M2*(1.0-H6_12), M2*(1.0-H6_12)])
b = np.add(b2,yy)
addPolytopeConstraint(m, A, xx, b)

#Feasible area 3

# foot1 step1 feasible_area3
xx = np.array([p1x_1, p1y_1, p1z_1])
yy = np.array([M3*(1.0-H1_13), M3*(1.0-H1_13), M3*(1.0-H1_13), M3*(1.0-H1_13), M3*(1.0-H1_13), M3*(1.0-H1_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area3
xx = np.array([p2x_1, p2y_1, p2z_1])
yy = np.array([M3*(1.0-H2_13), M3*(1.0-H2_13), M3*(1.0-H2_13), M3*(1.0-H2_13), M3*(1.0-H2_13), M3*(1.0-H2_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area3
xx = np.array([p3x_1, p3y_1, p3z_1])
yy = np.array([M3*(1.0-H3_13), M3*(1.0-H3_13), M3*(1.0-H3_13), M3*(1.0-H3_13), M3*(1.0-H3_13), M3*(1.0-H3_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area3
xx = np.array([p4x_1, p4y_1, p4z_1])
yy = np.array([M3*(1.0-H4_13), M3*(1.0-H4_13), M3*(1.0-H4_13), M3*(1.0-H4_13), M3*(1.0-H4_13), M3*(1.0-H4_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area3
xx = np.array([p5x_1, p2y_1, p5z_1])
yy = np.array([M3*(1.0-H5_13), M3*(1.0-H5_13), M3*(1.0-H5_13), M3*(1.0-H5_13), M3*(1.0-H5_13), M3*(1.0-H5_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area3
xx = np.array([p6x_1, p6y_1, p6z_1])
yy = np.array([M3*(1.0-H6_13), M3*(1.0-H6_13), M3*(1.0-H6_13), M3*(1.0-H6_13), M3*(1.0-H6_13), M3*(1.0-H6_13)])
b = np.add(b3,yy)
addPolytopeConstraint(m, A, xx, b)

######################################## Step 2 ##########################################
# Feasible area 1

# foot1 step1 feasible_area1
xx = np.array([p1x_2, p1y_2, p1z_2])
yy = np.array([M1 * (1.0 - H1_21), M1 * (1.0 - H1_21), M1 * (1.0 - H1_21), M1 * (1.0 - H1_21), M1 * (1.0 - H1_21), M1 * (1.0 - H1_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area1
xx = np.array([p2x_2, p2y_2, p2z_2])
yy = np.array([M1 * (1.0 - H2_21), M1 * (1.0 - H2_21), M1 * (1.0 - H2_21), M1 * (1.0 - H2_21), M1 * (1.0 - H2_21), M1 * (1.0 - H2_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area1
xx = np.array([p3x_2, p3y_2, p3z_2])
yy = np.array([M1 * (1.0 - H3_21), M1 * (1.0 - H3_21), M1 * (1.0 - H3_21), M1 * (1.0 - H3_21), M1 * (1.0 - H3_21), M1 * (1.0 - H3_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area1
xx = np.array([p4x_2, p4y_2, p4z_2])
yy = np.array([M1 * (1.0 - H4_21), M1 * (1.0 - H4_21), M1 * (1.0 - H4_21), M1 * (1.0 - H4_21), M1 * (1.0 - H4_21), M1 * (1.0 - H4_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area1
xx = np.array([p5x_2, p5y_2, p5z_2])
yy = np.array([M1 * (1.0 - H5_21), M1 * (1.0 - H5_21), M1 * (1.0 - H5_21), M1 * (1.0 - H5_21), M1 * (1.0 - H5_21), M1 * (1.0 - H5_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area1
xx = np.array([p6x_2, p6y_2, p6z_2])
yy = np.array([M1 * (1.0 - H6_21), M1 * (1.0 - H6_21), M1 * (1.0 - H6_21), M1 * (1.0 - H6_21), M1 * (1.0 - H6_21), M1 * (1.0 - H6_21)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)

# Feasible area 2

# foot1 step1 feasible_area2
xx = np.array([p1x_2, p1y_2, p1z_2])
yy = np.array([M2 * (1.0 - H1_22), M2 * (1.0 - H1_22), M2 * (1.0 - H1_22), M2 * (1.0 - H1_22), M2 * (1.0 - H1_22), M2 * (1.0 - H1_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area2
xx = np.array([p2x_2, p2y_2, p2z_2])
yy = np.array([M2 * (1.0 - H2_22), M2 * (1.0 - H2_22), M2 * (1.0 - H2_22), M2 * (1.0 - H2_22), M2 * (1.0 - H2_22), M2 * (1.0 - H2_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area2
xx = np.array([p3x_2, p3y_2, p3z_2])
yy = np.array([M2 * (1.0 - H3_22), M2 * (1.0 - H3_22), M2 * (1.0 - H3_22), M2 * (1.0 - H3_22), M2 * (1.0 - H3_22), M2 * (1.0 - H3_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area2
xx = np.array([p4x_2, p4y_2, p4z_2])
yy = np.array([M2 * (1.0 - H4_22), M2 * (1.0 - H4_22), M2 * (1.0 - H4_22), M2 * (1.0 - H4_22), M2 * (1.0 - H4_22), M2 * (1.0 - H4_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area2
xx = np.array([p5x_2, p5y_2, p5z_2])
yy = np.array([M2 * (1.0 - H5_22), M2 * (1.0 - H5_22), M2 * (1.0 - H5_22), M2 * (1.0 - H5_22), M2 * (1.0 - H5_22), M2 * (1.0 - H5_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area2
xx = np.array([p6x_2, p6y_2, p6z_2])
yy = np.array([M2 * (1.0 - H6_22), M2 * (1.0 - H6_22), M2 * (1.0 - H6_22), M2 * (1.0 - H6_22), M2 * (1.0 - H6_22), M2 * (1.0 - H6_22)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)

# Feasible area 3

# foot1 step1 feasible_area3
xx = np.array([p1x_2, p1y_2, p1z_2])
yy = np.array([M3 * (1.0 - H1_23), M3 * (1.0 - H1_23), M3 * (1.0 - H1_23), M3 * (1.0 - H1_23), M3 * (1.0 - H1_23), M3 * (1.0 - H1_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area3
xx = np.array([p2x_2, p2y_2, p2z_2])
yy = np.array([M3 * (1.0 - H2_23), M3 * (1.0 - H2_23), M3 * (1.0 - H2_23), M3 * (1.0 - H2_23), M3 * (1.0 - H2_23), M3 * (1.0 - H2_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area3
xx = np.array([p3x_2, p3y_2, p3z_2])
yy = np.array([M3 * (1.0 - H3_23), M3 * (1.0 - H3_23), M3 * (1.0 - H3_23), M3 * (1.0 - H3_23), M3 * (1.0 - H3_23), M3 * (1.0 - H3_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area3
xx = np.array([p4x_2, p4y_2, p4z_2])
yy = np.array([M3 * (1.0 - H4_23), M3 * (1.0 - H4_23), M3 * (1.0 - H4_23), M3 * (1.0 - H4_23), M3 * (1.0 - H4_23), M3 * (1.0 - H4_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area3
xx = np.array([p5x_2, p5y_2, p5z_2])
yy = np.array([M3 * (1.0 - H5_23), M3 * (1.0 - H5_23), M3 * (1.0 - H5_23), M3 * (1.0 - H5_23), M3 * (1.0 - H5_23), M3 * (1.0 - H5_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area3
xx = np.array([p6x_2, p6y_2, p6z_2])
yy = np.array([M3 * (1.0 - H6_23), M3 * (1.0 - H6_23), M3 * (1.0 - H6_23), M3 * (1.0 - H6_23), M3 * (1.0 - H6_23), M3 * (1.0 - H6_23)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)

######################################## Step 3 ##########################################
# Feasible area 1

# foot1 step1 feasible_area1
xx = np.array([p1x_3, p1y_3, p1z_3])
yy = np.array([M1 * (1.0 - H1_31), M1 * (1.0 - H1_31), M1 * (1.0 - H1_31), M1 * (1.0 - H1_31), M1 * (1.0 - H1_31), M1 * (1.0 - H1_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area1
xx = np.array([p2x_3, p2y_3, p2z_3])
yy = np.array([M1 * (1.0 - H2_31), M1 * (1.0 - H2_31), M1 * (1.0 - H2_31), M1 * (1.0 - H2_31), M1 * (1.0 - H2_31), M1 * (1.0 - H2_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area1
xx = np.array([p3x_3, p3y_3, p3z_3])
yy = np.array([M1 * (1.0 - H3_31), M1 * (1.0 - H3_31), M1 * (1.0 - H3_31), M1 * (1.0 - H3_31), M1 * (1.0 - H3_31), M1 * (1.0 - H3_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area1
xx = np.array([p4x_3, p4y_3, p4z_3])
yy = np.array([M1 * (1.0 - H4_31), M1 * (1.0 - H4_31), M1 * (1.0 - H4_31), M1 * (1.0 - H4_31), M1 * (1.0 - H4_31), M1 * (1.0 - H4_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area1
xx = np.array([p5x_3, p5y_3, p5z_3])
yy = np.array([M1 * (1.0 - H5_31), M1 * (1.0 - H5_31), M1 * (1.0 - H5_31), M1 * (1.0 - H5_31), M1 * (1.0 - H5_31), M1 * (1.0 - H5_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area1
xx = np.array([p6x_3, p6y_3, p6z_3])
yy = np.array([M1 * (1.0 - H6_31), M1 * (1.0 - H6_31), M1 * (1.0 - H6_31), M1 * (1.0 - H6_31), M1 * (1.0 - H6_31), M1 * (1.0 - H6_31)])
b = np.add(b1, yy)
addPolytopeConstraint(m, A, xx, b)

# Feasible area 2

# foot1 step1 feasible_area2
xx = np.array([p1x_3, p1y_3, p1z_3])
yy = np.array([M2 * (1.0 - H1_32), M2 * (1.0 - H1_32), M2 * (1.0 - H1_32), M2 * (1.0 - H1_32), M2 * (1.0 - H1_32), M2 * (1.0 - H1_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area2
xx = np.array([p2x_3, p2y_3, p2z_3])
yy = np.array([M2 * (1.0 - H2_32), M2 * (1.0 - H2_32), M2 * (1.0 - H2_32), M2 * (1.0 - H2_32), M2 * (1.0 - H2_32), M2 * (1.0 - H2_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area2
xx = np.array([p3x_3, p3y_3, p3z_3])
yy = np.array([M2 * (1.0 - H3_32), M2 * (1.0 - H3_32), M2 * (1.0 - H3_32), M2 * (1.0 - H3_32), M2 * (1.0 - H3_32), M2 * (1.0 - H3_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot1 step1 feasible_area2
xx = np.array([p4x_3, p4y_3, p4z_3])
yy = np.array([M2 * (1.0 - H4_32), M2 * (1.0 - H4_32), M2 * (1.0 - H4_32), M2 * (1.0 - H4_32), M2 * (1.0 - H4_32), M2 * (1.0 - H4_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area2
xx = np.array([p5x_3, p5y_3, p5z_3])
yy = np.array([M2 * (1.0 - H5_32), M2 * (1.0 - H5_32), M2 * (1.0 - H5_32), M2 * (1.0 - H5_32), M2 * (1.0 - H5_32), M2 * (1.0 - H5_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area2
xx = np.array([p6x_3, p6y_3, p6z_3])
yy = np.array([M2 * (1.0 - H6_32), M2 * (1.0 - H6_32), M2 * (1.0 - H6_32), M2 * (1.0 - H6_32), M2 * (1.0 - H6_32), M2 * (1.0 - H6_32)])
b = np.add(b2, yy)
addPolytopeConstraint(m, A, xx, b)

# Feasible area 3

# foot1 step1 feasible_area3
xx = np.array([p1x_3, p1y_3, p1z_3])
yy = np.array([M3 * (1.0 - H1_33), M3 * (1.0 - H1_33), M3 * (1.0 - H1_33), M3 * (1.0 - H1_33), M3 * (1.0 - H1_33), M3 * (1.0 - H1_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot2 step1 feasible_area3
xx = np.array([p2x_3, p2y_3, p2z_3])
yy = np.array([M3 * (1.0 - H2_33), M3 * (1.0 - H2_33), M3 * (1.0 - H2_33), M3 * (1.0 - H2_33), M3 * (1.0 - H2_33), M3 * (1.0 - H2_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot3 step1 feasible_area3
xx = np.array([p3x_3, p3y_3, p3z_3])
yy = np.array([M3 * (1.0 - H3_33), M3 * (1.0 - H3_33), M3 * (1.0 - H3_33), M3 * (1.0 - H3_33), M3 * (1.0 - H3_33), M3 * (1.0 - H3_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot4 step1 feasible_area3
xx = np.array([p4x_3, p4y_3, p4z_3])
yy = np.array([M3 * (1.0 - H4_33), M3 * (1.0 - H4_33), M3 * (1.0 - H4_33), M3 * (1.0 - H4_33), M3 * (1.0 - H4_33), M3 * (1.0 - H4_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot5 step1 feasible_area3
xx = np.array([p5x_3, p5y_3, p5z_3])
yy = np.array([M3 * (1.0 - H5_33), M3 * (1.0 - H5_33), M3 * (1.0 - H5_33), M3 * (1.0 - H5_33), M3 * (1.0 - H5_33), M3 * (1.0 - H5_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)
# foot6 step1 feasible_area3
xx = np.array([p6x_3, p6y_3, p6z_3])
yy = np.array([M3 * (1.0 - H6_33), M3 * (1.0 - H6_33), M3 * (1.0 - H6_33), M3 * (1.0 - H6_33), M3 * (1.0 - H6_33), M3 * (1.0 - H6_33)])
b = np.add(b3, yy)
addPolytopeConstraint(m, A, xx, b)


############################################### Set objective ##########################################################
g = np.array([750,  length*2.0+distance*2.0+800, step_height*2.0,\
			  1250, length*2.0+distance*2.0+800, step_height*2.0,\
			  1500, length*2.0+distance*2.0+400, step_height*2.0,\
			  1250, length*2.0+distance*2.0,     step_height*2.0,\
			  750,  length*2.0+distance*2.0,     step_height*2.0,\
			  500,  length*2.0+distance*2.0+400, step_height*2.0])   #Goal position
x = np.array([p1x_3, p1y_3, p1z_3, p2x_3, p2y_3, p2z_3, p3x_3, p3y_3, p3z_3, p4x_3, p4y_3, p4z_3, p5x_3, p5y_3, p5z_3, p6x_3, p6y_3, p6z_3])
Q = np.diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
setQuadraticObjective(m, x, g, Q)

m.optimize()

vars = m.getVars()
print('----------------------- Step 1 -----------------------')
print('p1x_1 =', vars[0], 'p1y_1 =', vars[1], 'p1z_1 =', vars[2])
print('p2x_1 =', vars[3], 'p2y_1 =', vars[4], 'p2z_1 =', vars[5])
print('p3x_1 =', vars[6], 'p3y_1 =', vars[7], 'p3z_1 =', vars[8])
print('p4x_1 =', vars[9], 'p4y_1 =', vars[10], 'p4z_1 =', vars[11])
print('p5x_1 =', vars[12], 'p5y_1 =', vars[13], 'p5z_1 =', vars[14])
print('p6x_1 =', vars[15], 'p6y_1 =', vars[16], 'p6z_1 =', vars[17])

print('----------------------- Step 2 -----------------------')
print('p1x_2 =', vars[18], 'p1y_2 =', vars[19], 'p1z_2 =', vars[20])
print('p2x_2 =', vars[21], 'p2y_2 =', vars[22], 'p2z_2 =', vars[23])
print('p3x_2 =', vars[24], 'p3y_2 =', vars[25], 'p3z_2 =', vars[26])
print('p4x_2 =', vars[27], 'p4y_2 =', vars[28], 'p4z_2 =', vars[29])
print('p5x_2 =', vars[30], 'p5y_2 =', vars[31], 'p5z_2 =', vars[32])
print('p6x_2 =', vars[33], 'p6y_2 =', vars[34], 'p6z_2 =', vars[35])

print('----------------------- Step 3 -----------------------')
print('p1x_3 =', vars[36], 'p1y_3 =', vars[37], 'p1z_3 =', vars[38])
print('p2x_3 =', vars[39], 'p2y_3 =', vars[40], 'p2z_3 =', vars[41])
print('p3x_3 =', vars[42], 'p3y_3 =', vars[43], 'p3z_3 =', vars[44])
print('p4x_3 =', vars[45], 'p4y_3 =', vars[46], 'p4z_3 =', vars[47])
print('p5x_3 =', vars[48], 'p5y_3 =', vars[49], 'p5z_3 =', vars[50])
print('p6x_3 =', vars[51], 'p6y_3 =', vars[52], 'p6z_3 =', vars[53])

#How to read H values?
# First index is leg count, second index is step, third index is feasible area
# So H1_11 = 1, H1_12 = 0, H1_13 = 0 means the first leg at step one is in first feasible area
print('----------------------- Binary Variables -----------------------')
print('H1_11 =', vars[54], 'H1_12 =', vars[55], 'H1_13 =', vars[56])
print('H1_21 =', vars[57], 'H1_22 =', vars[58], 'H1_23 =', vars[59])
print('H1_31 =', vars[60], 'H1_32 =', vars[61], 'H1_33 =', vars[62])
print('H2_11 =', vars[63], 'H2_12 =', vars[64], 'H2_13 =', vars[65])
print('H2_21 =', vars[66], 'H2_22 =', vars[67], 'H2_23 =', vars[68])
print('H2_31 =', vars[69], 'H2_32 =', vars[70], 'H2_33 =', vars[71])
print('H3_11 =', vars[72], 'H3_12 =', vars[73], 'H3_13 =', vars[74])
print('H3_21 =', vars[75], 'H3_22 =', vars[76], 'H3_23 =', vars[77])
print('H3_31 =', vars[78], 'H3_32 =', vars[79], 'H3_33 =', vars[80])
print('H4_11 =', vars[81], 'H4_12 =', vars[82], 'H4_13 =', vars[83])
print('H4_21 =', vars[84], 'H4_22 =', vars[85], 'H4_23 =', vars[86])
print('H4_31 =', vars[87], 'H4_32 =', vars[88], 'H4_33 =', vars[89])
print('H5_11 =', vars[90], 'H5_12 =', vars[91], 'H5_13 =', vars[92])
print('H5_21 =', vars[93], 'H5_22 =', vars[94], 'H5_23 =', vars[95])
print('H5_31 =', vars[96], 'H5_32 =', vars[97], 'H5_33 =', vars[98])
print('H6_11 =', vars[99], 'H6_12 =', vars[100], 'H6_13 =', vars[101])
print('H6_21 =', vars[102], 'H6_22 =', vars[103], 'H6_23 =', vars[104])
print('H6_31 =', vars[105], 'H6_32 =', vars[106], 'H6_33 =', vars[107])