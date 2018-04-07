#!/usr/bin/python

import gurobipy as go
import numpy as np
import nlopt
import time

class Vector3D:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def __getitem__(self):
		return self

	def __neg__(self):
		return Vector3D(-self.x, -self.y, -self.z)

	def norm(self):
		return np.sqrt(self.x^2+self.y^2+self.z^2)

def addNormConstr(model, vec1, vec2, vec_offset, th):
	#Add a constraint of the form norm(vec1+vec2+vec_offset))<=threshold
	#Input variables are Vector3D type
	x = vec1.x + vec2.x + vec_offset.x
	y = vec1.y + vec2.y + vec_offset.y
	z = vec1.z + vec2.z + vec_offset.z
	expr = go.QuadExpr(x*x+y*y+z*z)  #Still needs to test it to see if it works as expected
	model.addQConstr(expr <= (th*th))

def addPolytopeConstr(model, A, x, b):
	#Add a constraint of the form Ax<=b
	#Not checking if column(A)=dim(x), make sure they match!
	dim_row = np.size(b)
	dim_col = np.size(x)
	for i in range(dim_row):
		expr = (go.quicksum(A[i,j]*x[j] for j in range(dim_col)) <= b[i])
		model.addConstr(expr)

def setQuadraticObjective(model, x, x0, Q):
	#Add a objective function of the form minimize(x) (x-x0)'Q(x-x0)
	#not checking if g and Q have correspoinding size. Make sure that they match!
	dim = np.size(x)
	expr = go.quicksum((x[i]-x0[i])*Q[i,j]*(x[j]-x0[j]) for i in range(dim) for j in range(dim))
	model.setObjective(expr, go.GRB.MINIMIZE)

# Create a new model
m = go.Model("setpOnStairs")

########################################### Create Footstep Variables ##################################################
N = 10   #Total number of rounds
S = 6    #Number of footsteps to be planned each round

#Initial footstep condition
footstep0 = np.zeros([1,S,3])
footstep0[0,0,0] = 750.0;  footstep0[0,0,1] = 800.0;  footstep0[0,0,2] = 0.0;
footstep0[0,1,0] = 1250.0; footstep0[0,1,1] = 800.0;  footstep0[0,1,2] = 0.0;
footstep0[0,2,0] = 1500.0; footstep0[0,2,1] = 400.0;  footstep0[0,2,2] = 0.0;
footstep0[0,3,0] = 1250.0; footstep0[0,3,1] = 0.0;    footstep0[0,3,2] = 0.0;
footstep0[0,4,0] = 750.0;  footstep0[0,4,1] = 0.0;    footstep0[0,4,2] = 0.0;
footstep0[0,5,0] = 500.0;  footstep0[0,5,1] = 400.0;  footstep0[0,5,2] = 0.0;

# The footstep variables are NxS matrix of 3D vectors
# First index is round, second index is footsteps, third index is xyz
footsteps = m.addVars(N,S,3)   #Note here footsteps[0] is actually the first step

########################################### Create Footstep on Feasible Area Variables #################################
#Description of feasible areas
#Number of feasible areas
F = 3
#Length of a step area
length_step = 300.0
#Witdh of a step area
width_step = 2000.0
#Length before reaching a step
length_before = 1200.0
#Length after finishing the steps
length_after = 3000.0
#Height of a feasible box area
box_height = 1.0
#Step height, the height of each step
step_height = 130.0    #This is based on the step height of the first floor stairs outside
#Distance between two feasible areas
distance = 0.0
#Expression of a polytope
#Remember: always row vector
A = np.array([[ 1.0,  0.0,  0.0],
			  [-1.0,  0.0,  0.0],
			  [ 0.0,  1.0,  0.0],
			  [ 0.0, -1.0,  0.0],
			  [ 0.0,  0.0,  1.0],
			  [ 0.0,  0.0, -1.0]])
b = []
b.append(np.array([width_step,  0.0,
				   length_before, 0.0,
			       box_height, 0.0]))

b.append(np.array([width_step, 0.0,
				   length_before+length_step+distance,  -(length_before+distance),
			       box_height+step_height, -step_height]))

b.append(np.array([width_step, 0.0,
				   length_before+length_step+2*distance+length_after, -(length_before+length_step+2*distance),
			       box_height+2*step_height, -2*step_height]))

#Goal position
offset_after = 200  #The amound of offset that robot wants to stay after finishing the steps
x_goal = np.array([ 750,  length_before+length_step+2*distance+offset_after+800, step_height*2.0,
			   1250,  length_before+length_step+2*distance+offset_after+800, step_height*2.0,
			   1500,  length_before+length_step+2*distance+offset_after+400, step_height*2.0,
			   1250,  length_before+length_step+2*distance+offset_after,     step_height*2.0,
			    750,  length_before+length_step+2*distance+offset_after,     step_height*2.0,
			    500,  length_before+length_step+2*distance+offset_after+400, step_height*2.0])

print x_goal

#Big-M constants
M = []
for i in range(F):
	M.append(10000.0)

# Array of footsteps on feasible area
# First index is round, second index is footsteps, third index is count of feasible area
H = m.addVars(N,S,F,vtype=go.GRB.BINARY)

########################################## Create Footstep Trim Variables ##############################################

trim = m.addVars(N,vtype=go.GRB.BINARY)

#Update the model to incorporate new variables
m.update()

###################################### Constraint 1: Bounds on the step size ###########################################
#This is the constraint between this step and next step
th_stepSize = 380   #Here, not sure why the first step is so large if the initial condition is set to -800,
					#needs stepsize=800 if initial condition is set to -800
					#For some reason, if the initial condition is negative, will have issue with the first step!
					#Maybe it is because of the big-M expression

#First step: starting from initial condition
for j in range(S):
	v1 = Vector3D(footstep0[0, j, 0], footstep0[0, j, 1], footstep0[0, j, 2])
	v2 = Vector3D(footsteps[0, j, 0], footsteps[0, j, 1], footsteps[0, j, 2])
	v_offset = Vector3D(0.0, 0.0, 0.0)
	addNormConstr(m, v1, -v2, v_offset, th_stepSize)

for i in range(1,N):
	for j in range(S):
		v1 = Vector3D(footsteps[i-1, j, 0], footsteps[i-1, j, 1], footsteps[i-1, j, 2])
		v2 = Vector3D(footsteps[i, j, 0], footsteps[i, j, 1], footsteps[i, j, 2])
		v_offset = Vector3D(0.0, 0.0, 0.0)
		addNormConstr(m, v1, -v2, v_offset, th_stepSize)

########################################## Constraint 2: Leg workspace #################################################
#This is the constraint between footsteps within the same step
#For this constraint, we assume that the body orientation is unchanged, with y forward, x to the right
#i.e. not optimize over rotation matrix R
#Rotation matrices
rotMatrix1 = np.array([[np.cos(0.0 / 180.0 * np.pi),    -np.sin(0.0 / 180.0 * np.pi),    0.0],
					   [np.sin(0.0 / 180.0 * np.pi),     np.cos(0.0 / 180.0 * np.pi),    0.0], [0.0, 0.0, 1.0]])
rotMatrix2 = np.array([[np.cos(-60.0 / 180.0 * np.pi),  -np.sin(-60.0 / 180.0 * np.pi),  0.0],
					   [np.sin(-60.0 / 180.0 * np.pi),   np.cos(-60.0 / 180.0 * np.pi),  0.0], [0.0, 0.0, 1.0]])
rotMatrix3 = np.array([[np.cos(-120.0 / 180.0 * np.pi), -np.sin(-120.0 / 180.0 * np.pi), 0.0],
					   [np.sin(-120.0 / 180.0 * np.pi),  np.cos(-120.0 / 180.0 * np.pi), 0.0], [0.0, 0.0, 1.0]])
rotMatrix4 = np.array([[np.cos(-180.0 / 180.0 * np.pi), -np.sin(-180.0 / 180.0 * np.pi), 0.0],
					   [np.sin(-180.0 / 180.0 * np.pi),  np.cos(-180.0 / 180.0 * np.pi), 0.0], [0.0, 0.0, 1.0]])
rotMatrix5 = np.array([[np.cos(-240.0 / 180.0 * np.pi), -np.sin(-240.0 / 180.0 * np.pi), 0.0],
					   [np.sin(-240.0 / 180.0 * np.pi),  np.cos(-240.0 / 180.0 * np.pi), 0.0], [0.0, 0.0, 1.0]])
rotMatrix6 = np.array([[np.cos(-300.0 / 180.0 * np.pi), -np.sin(-300.0 / 180.0 * np.pi), 0.0],
					   [np.sin(-300.0 / 180.0 * np.pi),  np.cos(-300.0 / 180.0 * np.pi), 0.0], [0.0, 0.0, 1.0]])


#Center offset vector
p_left =  np.array([0.0,   0.0, 0.0])    #Left vector is from center to the left
p_right = np.array([600.0, 0.0, 0.0])   #Right vector is from center ot the right
r1_left  = rotMatrix1.dot(p_left)
r1_right = rotMatrix1.dot(p_right)
r2_left  = rotMatrix2.dot(p_left)
r2_right = rotMatrix2.dot(p_right)
r3_left  = rotMatrix3.dot(p_left)
r3_right = rotMatrix3.dot(p_right)
r4_left  = rotMatrix4.dot(p_left)
r4_right = rotMatrix4.dot(p_right)
r5_left  = rotMatrix5.dot(p_left)
r5_right = rotMatrix5.dot(p_right)
r6_left  = rotMatrix6.dot(p_left)
r6_right = rotMatrix6.dot(p_right)

radi_left = 600   #Nominal 600
radi_right = 275   #Nominal 275

for i in range(N):
	#1 to 2
	v1 = Vector3D(footsteps[i, 0, 0], footsteps[i, 0, 1], footsteps[i, 0, 2])
	v2 = Vector3D(footsteps[i, 1, 0], footsteps[i, 1, 1], footsteps[i, 1, 2])
	vec_offset_left = Vector3D(r1_left[0], r1_left[1], r1_left[2])
	vec_offset_right = Vector3D(r1_right[0], r1_right[1], r1_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left,  radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

	# 2 to 3
	v1 = Vector3D(footsteps[i, 1, 0], footsteps[i, 1, 1], footsteps[i, 1, 2])
	v2 = Vector3D(footsteps[i, 2, 0], footsteps[i, 2, 1], footsteps[i, 2, 2])
	vec_offset_left = Vector3D(r2_left[0], r2_left[1], r2_left[2])
	vec_offset_right = Vector3D(r2_right[0], r2_right[1], r2_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left, radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

	#3 to 4
	v1 = Vector3D(footsteps[i, 2, 0], footsteps[i, 2, 1], footsteps[i, 2, 2])
	v2 = Vector3D(footsteps[i, 3, 0], footsteps[i, 3, 1], footsteps[i, 3, 2])
	vec_offset_left = Vector3D(r3_left[0], r3_left[1], r3_left[2])
	vec_offset_right = Vector3D(r3_right[0], r3_right[1], r3_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left,  radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

	#4 to 5
	v1 = Vector3D(footsteps[i, 3, 0], footsteps[i, 3, 1], footsteps[i, 3, 2])
	v2 = Vector3D(footsteps[i, 4, 0], footsteps[i, 4, 1], footsteps[i, 4, 2])
	vec_offset_left = Vector3D(r4_left[0], r4_left[1], r4_left[2])
	vec_offset_right = Vector3D(r4_right[0], r4_right[1], r4_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left, radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

	#5 to 6
	v1 = Vector3D(footsteps[i, 4, 0], footsteps[i, 4, 1], footsteps[i, 4, 2])
	v2 = Vector3D(footsteps[i, 5, 0], footsteps[i, 5, 1], footsteps[i, 5, 2])
	vec_offset_left = Vector3D(r5_left[0], r5_left[1], r5_left[2])
	vec_offset_right = Vector3D(r5_right[0], r5_right[1], r5_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left,  radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

	#6 to 1
	v1 = Vector3D(footsteps[i, 5, 0], footsteps[i, 5, 1], footsteps[i, 5, 2])
	v2 = Vector3D(footsteps[i, 0, 0], footsteps[i, 0, 1], footsteps[i, 0, 2])
	vec_offset_left = Vector3D(r6_left[0], r6_left[1], r6_left[2])
	vec_offset_right = Vector3D(r6_right[0], r6_right[1], r6_right[2])
	addNormConstr(m, v1, -v2, vec_offset_left,  radi_left)
	addNormConstr(m, v1, -v2, vec_offset_right, radi_right)

############################## Constraint 3: Footsteps within feasible area ############################################
for i in range(N):
	for j in range(S):
		for k in range(F):
			xx = np.array([footsteps[i,j,0], footsteps[i,j,1], footsteps[i,j,2]])
			yy = np.array([M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k]),
				           M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k])])
			bb = np.add(b[k], yy)
			addPolytopeConstr(m, A, xx, bb)

#Within one round, each footstep can only be within one region
for i in range(S):
	for j in range(N):
		sum = go.quicksum(H[j,i,k] for k in range(F))
		m.addConstr( sum == 1)

###########################  Constraint 4: if trim(i)==1 then footstep i is set to initial footstep ####################
m_trim = -10000.0
M_trim =  10000.0

for i in range(N):
	for j in range(S):
		for k in range(F):
			m.addConstr( (footsteps[i, j, k] - footstep0[0, j, k]) >= m_trim * trim[i] )
			m.addConstr( (footsteps[i, j, k] - footstep0[0, j, k]) <= M_trim * trim[i] )


############################################### Set objective ##########################################################
x = np.array([footsteps[N - 1, 0, 0], footsteps[N - 1, 0, 1], footsteps[N - 1, 0, 2],
			  footsteps[N - 1, 1, 0], footsteps[N - 1, 1, 1], footsteps[N - 1, 1, 2],
			  footsteps[N - 1, 2, 0], footsteps[N - 1, 2, 1], footsteps[N - 1, 2, 2],
			  footsteps[N - 1, 3, 0], footsteps[N - 1, 3, 1], footsteps[N - 1, 3, 2],
			  footsteps[N - 1, 4, 0], footsteps[N - 1, 4, 1], footsteps[N - 1, 4, 2],
			  footsteps[N - 1, 5, 0], footsteps[N - 1, 5, 1], footsteps[N - 1, 5, 2]])

Q = np.diag([1.0, 1.0, 1.0,
			 1.0, 1.0, 1.0,
			 1.0, 1.0, 1.0,
			 1.0, 1.0, 1.0,
			 1.0, 1.0, 1.0,
			 1.0, 1.0, 1.0 ])

#Quadratic expression for constraints on goal footstep position
dim1 = np.size(x)
expr1 = go.quicksum((x[i]-x_goal[i])*Q[i,j]*(x[j]-x_goal[j]) for i in range(dim1) for j in range(dim1))
#Quadratic expression for constraints on footstep length
weight1 = 0.01
expr2 = go.quicksum((footsteps[i,j,k]-footsteps[i-1,j,k])*weight1*(footsteps[i,j,k]-footsteps[i-1,j,k]) for i in range(1,N) for j in range(S) for k in range(3))
#Punishment on trim varaibles
weight2 = 10
expr3 = go.quicksum(weight2*trim[i] for i in range(N))
obj = expr1 + expr2 + expr3
m.setObjective(obj, go.GRB.MINIMIZE)

m.optimize()

footstep_1_RF = [footsteps[0,1,0].X, footsteps[0,1,1].X, footsteps[0,1,2].X]
footstep_1_RM = [footsteps[0,2,0].X, footsteps[0,2,1].X, footsteps[0,2,2].X]
footstep_1_RR = [footsteps[0,3,0].X, footsteps[0,3,1].X, footsteps[0,3,2].X]
footstep_1_LR = [footsteps[0,4,0].X, footsteps[0,4,1].X, footsteps[0,4,2].X]
footstep_1_LM = [footsteps[0,5,0].X, footsteps[0,5,1].X, footsteps[0,5,2].X]
footstep_1_LF = [footsteps[0,0,0].X, footsteps[0,0,1].X, footsteps[0,0,2].X]

footstep_2_RF = [footsteps[1,1,0].X, footsteps[1,1,1].X, footsteps[1,1,2].X]
footstep_2_RM = [footsteps[1,2,0].X, footsteps[1,2,1].X, footsteps[1,2,2].X]
footstep_2_RR = [footsteps[1,3,0].X, footsteps[1,3,1].X, footsteps[1,3,2].X]
footstep_2_LR = [footsteps[1,4,0].X, footsteps[1,4,1].X, footsteps[1,4,2].X]
footstep_2_LM = [footsteps[1,5,0].X, footsteps[1,5,1].X, footsteps[1,5,2].X]
footstep_2_LF = [footsteps[1,0,0].X, footsteps[1,0,1].X, footsteps[1,0,2].X]

for i in range(N):
	print("------------Round {} ------------------".format(i+1))
	print("Foot 1 [{}, {}, {}]".format(footsteps[i, 0, 0], footsteps[i, 0, 1], footsteps[i, 0, 2]))
	print("Foot 2 [{}, {}, {}]".format(footsteps[i, 1, 0], footsteps[i, 1, 1], footsteps[i, 1, 2]))
	print("Foot 3 [{}, {}, {}]".format(footsteps[i, 2, 0], footsteps[i, 2, 1], footsteps[i, 2, 2]))
	print("Foot 4 [{}, {}, {}]".format(footsteps[i, 3, 0], footsteps[i, 3, 1], footsteps[i, 3, 2]))
	print("Foot 5 [{}, {}, {}]".format(footsteps[i, 4, 0], footsteps[i, 4, 1], footsteps[i, 4, 2]))
	print("Foot 6 [{}, {}, {}]".format(footsteps[i, 5, 0], footsteps[i, 5, 1], footsteps[i, 5, 2]))

############################################### The end of MICP part ###############################################

N = 21  #Dimension along time axis
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

#x y upper/lower bounds should be set appropriately!!!
xUpperBound = 3000.0
xLowerBound = 0.0
yUpperBound = 4000.0
yLowerBound = 0.0
zUpperBound = 500.0
zLowerBound = -10.0

#Need to calculate r0, rS, rN values
bodyHeight = 160.0 #The minimal height of body from the ground

r0x = ((footstep_1_RM[0] + footstep_1_LM[0])/2.0 + (footstep_1_RF[0] + footstep_1_LF[0])/2.0 + (footstep_1_RR[0] + footstep_1_LR[0])/2.0)/3.0
r0y = ((footstep_1_RF[1] + footstep_1_RR[1])/2.0 + (footstep_1_LF[1] + footstep_1_LR[1])/2.0)/2.0
r0z = (footstep_1_RF[2]+footstep_1_RM[2]+footstep_1_RR[2]+footstep_1_LR[2]+footstep_1_LM[2]+footstep_1_LF[2])/6.0+bodyHeight

rNx = ((footstep_2_RM[0] + footstep_2_LM[0])/2.0 + (footstep_2_RF[0] + footstep_2_LF[0])/2.0 + (footstep_2_RR[0] + footstep_2_LR[0])/2.0)/3.0
rNy = ((footstep_2_RF[1] + footstep_2_RR[1])/2.0 + (footstep_2_LF[1] + footstep_2_LR[1])/2.0)/2.0
rNz = (footstep_2_RF[2]+footstep_2_RM[2]+footstep_2_RR[2]+footstep_2_LR[2]+footstep_2_LM[2]+footstep_2_LF[2])/6.0+bodyHeight

rSx = (r0x + rNx)/2.0
rSy = (r0y + rNy)/2.0
rSz = (r0z + rNz)/2.0

deg2rad_const = np.pi/180.0

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
        grad[0, 6+24*step] = -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[6+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[6+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[6+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[8+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 7+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 8+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 9+24*step] =  -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[9+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[9+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[9+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[11+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 10+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 11+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 12+24*step] = -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[12+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[12+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[12+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[14+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 13+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 14+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 15+24*step] = -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[15+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[15+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[15+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[17+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 16+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 17+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 18+24*step] = -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[18+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[18+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[18+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[20+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 19+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 20+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 21+24*step] = -deg2rad_const*(L_coxa*m2*np.sin(deg2rad_const*x[21+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[21+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[21+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[23+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 22+24*step] = -deg2rad_const*(L_femur*m3*np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[0, 23+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

        grad[1, 0+24*step] = 0.0
        grad[1, 1+24*step] = 1.0
        grad[1, 2+24*step] = 0.0
        grad[1, 3+24*step] = 0.0
        grad[1, 4+24*step] = -1.0
        grad[1, 5+24*step] = 0.0
        grad[1, 6+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[6+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[6+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[6+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 7+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 8+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 9+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[9+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[9+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[9+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 10+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 11+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 12+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[12+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[12+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[12+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 13+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 14+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 15+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[15+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[15+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[15+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 16+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 17+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 18+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[18+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[18+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[18+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 19+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 20+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 21+24*step] = deg2rad_const*(L_coxa*m2*np.cos(deg2rad_const*x[21+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[21+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[21+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 22+24*step] = -deg2rad_const*(L_femur*m3*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]) + Ym3*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[1, 23+24*step] = -deg2rad_const*(Ym3*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

        grad[2, 0+24*step] = 0.0
        grad[2, 1+24*step] = 0.0
        grad[2, 2+24*step] = 1.0
        grad[2, 3+24*step] = 0.0
        grad[2, 4+24*step] = 0.0
        grad[2, 5+24*step] = -1.0
        grad[2, 6+24*step] = 0.0
        grad[2, 7+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[7+24*step] + deg2rad_const*x[8+24*step]) - L_femur*np.cos(deg2rad_const*x[7+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[7+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 8+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[7+24*step] + deg2rad_const*x[8+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 9+24*step] = 0.0
        grad[2, 10+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[10+24*step] + deg2rad_const*x[11+24*step]) - L_femur*np.cos(deg2rad_const*x[10+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[10+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 11+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[10+24*step] + deg2rad_const*x[11+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 12+24*step] = 0.0
        grad[2, 13+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[13+24*step] + deg2rad_const*x[14+24*step]) - L_femur*np.cos(deg2rad_const*x[13+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[13+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 14+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[13+24*step] + deg2rad_const*x[14+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 15+24*step] = 0.0
        grad[2, 16+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[16+24*step] + deg2rad_const*x[17+24*step]) - L_femur*np.cos(deg2rad_const*x[16+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[16+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 17+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[16+24*step] + deg2rad_const*x[17+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 18+24*step] = 0.0
        grad[2, 19+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[19+24*step] + deg2rad_const*x[20+24*step]) - L_femur*np.cos(deg2rad_const*x[19+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[19+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 20+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[19+24*step] + deg2rad_const*x[20+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 21+24*step] = 0.0
        grad[2, 22+24*step] = -deg2rad_const*(m3*(Ym3*np.sin(deg2rad_const*x[22+24*step] + deg2rad_const*x[23+24*step]) - L_femur*np.cos(deg2rad_const*x[22+24*step])) - Xm2*m2*np.cos(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)
        grad[2, 23+24*step] = -deg2rad_const*(Ym3*m3*np.sin(deg2rad_const*x[22+24*step] + deg2rad_const*x[23+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)

    #This CM has nothing to do with body dimension, think about vector summation to convince yourself
    result[:] = np.array(
        [x[0+24*step] - x[3+24*step] + \
        (L_coxa*m2*np.cos(deg2rad_const*x[6+24*step]) + L_coxa*m2*np.cos(deg2rad_const*x[9+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[6+24*step]) + L_coxa*m2*np.cos(deg2rad_const*x[12+24*step]) +
        L_coxa*m3*np.cos(deg2rad_const*x[9+24*step]) + L_coxa*m2*np.cos(deg2rad_const*x[15+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[12+24*step]) + L_coxa*m2*np.cos(deg2rad_const*x[18+24*step]) +
        L_coxa*m3*np.cos(deg2rad_const*x[15+24*step]) + L_coxa*m2*np.cos(deg2rad_const*x[21+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[18+24*step]) + L_coxa*m3*np.cos(deg2rad_const*x[21+24*step]) +
        L_femur*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) +
        L_femur*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) +
        Xm1*m1*np.cos(deg2rad_const*x[6+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[9+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[12+24*step]) +
        Xm1*m1*np.cos(deg2rad_const*x[15+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[18+24*step]) + Xm1*m1*np.cos(deg2rad_const*x[21+24*step]) +
        Xm2*m2*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) +
        Xm2*m2*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) - Ym3*m3*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb),

        x[1+24*step] - x[4+24*step] + \
        (L_coxa*m2*np.sin(deg2rad_const*x[6+24*step]) + L_coxa*m2*np.sin(deg2rad_const*x[9+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[6+24*step]) + L_coxa*m2*np.sin(deg2rad_const*x[12+24*step]) +
        L_coxa*m3*np.sin(deg2rad_const*x[9+24*step]) + L_coxa*m2*np.sin(deg2rad_const*x[15+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[12+24*step]) + L_coxa*m2*np.sin(deg2rad_const*x[18+24*step]) +
        L_coxa*m3*np.sin(deg2rad_const*x[15+24*step]) + L_coxa*m2*np.sin(deg2rad_const*x[21+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[18+24*step]) + L_coxa*m3*np.sin(deg2rad_const*x[21+24*step]) +
        L_femur*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]) +
        L_femur*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]) + L_femur*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) +
        Xm1*m1*np.sin(deg2rad_const*x[6+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[9+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[12+24*step]) +
        Xm1*m1*np.sin(deg2rad_const*x[15+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[18+24*step]) + Xm1*m1*np.sin(deg2rad_const*x[21+24*step]) +
        Xm2*m2*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]) +
        Xm2*m2*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]) + Xm2*m2*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[8+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[11+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[14+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[17+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[20+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step]) -
        Ym3*m3*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[23+24*step]) - Ym3*m3*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb),

        x[2+24*step] - x[5+24*step] + \
        (m3*(Ym3*np.cos(deg2rad_const*x[7+24*step] + deg2rad_const*x[8+24*step]) + L_femur*np.sin(deg2rad_const*x[7+24*step])) + m3*(Ym3*np.cos(deg2rad_const*x[10+24*step] + deg2rad_const*x[11+24*step]) + L_femur*np.sin(deg2rad_const*x[10+24*step])) +
        m3*(Ym3*np.cos(deg2rad_const*x[13+24*step] + deg2rad_const*x[14+24*step]) + L_femur*np.sin(deg2rad_const*x[13+24*step])) + m3*(Ym3*np.cos(deg2rad_const*x[16+24*step] + deg2rad_const*x[17+24*step]) + L_femur*np.sin(deg2rad_const*x[16+24*step])) +
        m3*(Ym3*np.cos(deg2rad_const*x[19+24*step] + deg2rad_const*x[20+24*step]) + L_femur*np.sin(deg2rad_const*x[19+24*step])) + m3*(Ym3*np.cos(deg2rad_const*x[22+24*step] + deg2rad_const*x[23+24*step]) + L_femur*np.sin(deg2rad_const*x[22+24*step])) +
        Xm2*m2*np.sin(deg2rad_const*x[7+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[10+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[13+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[16+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[19+24*step]) + Xm2*m2*np.sin(deg2rad_const*x[22+24*step]))/(6*m1 + 6*m2 + 6*m3 + mb)])

    return result

#------------------------------------ Constraint 2 -----------------------------------------------------------------
def r0(result, x, grad):

    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 0] = 1.0
        grad[1, 1] = 1.0
        grad[2, 2] = 1.0

    result[:] = np.array([x[0]-r0x, x[1]-r0y, x[2]-r0z])
    #print("r0 residual is %f %f %f" %(result[0], result[1], result[2]))
    return result

def rS(result, x, grad):
    # rS = r[set 11] for N=21
    # For N=21, variable set 11 is in the middle, that means (11-1)*24 -> 10*24
    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 24*((N-1)/2)+0] = 1.0
        grad[1, 24*((N-1)/2)+1] = 1.0
        grad[2, 24*((N-1)/2)+2] = 1.0

    result[:] = np.array([x[24*((N-1)/2)]-rSx, x[24*((N-1)/2)+1]-rSy, x[24*((N-1)/2)+2]-rSz])
    #print("rS residual is %f %f %f" % (result[0], result[1], result[2]))
    return result

def rN(result, x, grad):

    if grad.size > 0:
        grad[:] = np.zeros([3, 24*N])
        grad[0, 24*(N-1)+0] = 1.0
        grad[1, 24*(N-1)+1] = 1.0
        grad[2, 24*(N-1)+2] = 1.0

    result[:] = np.array([x[24*(N-1)]-rNx, x[24*(N-1)+1]-rNy, x[24*(N-1)+2]-rNz])
    #print("rN residual is %f %f %f" % (result[0], result[1], result[2]))
    return result

#------------------------------------ Constraint 3 -----------------------------------------------------------------
def leg1_end_point_XYZ(result, x, grad, step, px, py, pz):

    #This is RF leg, the algebra is compared with matlab to verify it's correctness. Should carryout actual number to double check it!
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[6+24*step]) + L_tibia*(np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) +
        L_femur*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) + x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[6+24*step]) + L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) +
        L_femur*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) + y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[7+24*step]) - L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))])

    print("leg 1 end point results px=%f py=%f pz=%f res_x=%f res_y=%f res_z=%f" % (px, py, pz, result[0], result[1], result[2]))
    #Print this to verify convergence, should be 0 when converge!

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 6+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[6+24*step]) - L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) -
                             L_femur*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]))
        grad[0, 7+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step])) -
                             L_femur*np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]))
        grad[0, 8+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.cos(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 6+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[6+24*step]) + L_tibia*(np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) +
                             L_femur*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]))
        grad[1, 7+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step])) -
                             L_femur*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step]))
        grad[1, 8+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 6+24*step] = 0.0
        grad[2, 7+24*step] = deg2rad_const*(L_tibia * (np.cos(deg2rad_const*x[7+24*step]) * np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[8+24*step]) * np.sin(deg2rad_const*x[7+24*step])) + L_femur * np.cos(deg2rad_const*x[7+24*step]))
        grad[2, 8+24*step] = deg2rad_const*(L_tibia * (np.cos(deg2rad_const*x[7+24*step]) * np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[8+24*step]) * np.sin(deg2rad_const*x[7+24*step])))

    return result

def leg2_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is RM leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[9+24*step]) + L_tibia*(np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) +
        L_femur*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) + x_body_mid_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[9+24*step]) + L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) +
        L_femur*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]),

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[10+24*step]) - L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 9+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[9+24*step]) - L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) -
                             L_femur*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]))
        grad[0, 10+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]))
        grad[0, 11+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.cos(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 9+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[9+24*step]) + L_tibia*(np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) +
                             L_femur*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]))
        grad[1, 10+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step])) -
                              L_femur*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step]))
        grad[1, 11+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 9+24*step] = 0.0
        grad[2, 10+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) + L_femur*np.cos(deg2rad_const*x[10+24*step]))
        grad[2, 11+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])))

    return result

def leg3_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is RR leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[12+24*step]) + L_tibia*(np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) +
        L_femur*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) + x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[12+24*step]) + L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) +
        L_femur*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]) - y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[13+24*step]) - L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 12+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[12+24*step]) - L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) -
                               L_femur*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step]))
        grad[0, 13+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]))
        grad[0, 14+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.cos(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 12+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[12+24*step]) + L_tibia*(np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) +
                              L_femur*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]))
        grad[1, 13+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step])) -
                               L_femur*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step]))
        grad[1, 14+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 12+24*step] = 0.0
        grad[2, 13+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) + L_femur*np.cos(deg2rad_const*x[13+24*step]))
        grad[2, 14+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])))

    return result


def leg4_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LR leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[15+24*step]) + L_tibia*(np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) +
        L_femur*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) - x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[15+24*step]) + L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) +
        L_femur*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) - y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[16+24*step]) - L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 15+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[15+24*step]) - L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]))
        grad[0, 16+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]))
        grad[0, 17+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.cos(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 15+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[15+24*step]) + L_tibia*(np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) +
                              L_femur*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]))
        grad[1, 16+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step])) -
                              L_femur*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step]))
        grad[1, 17+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 15+24*step] = 0.0
        grad[2, 16+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) + L_femur*np.cos(deg2rad_const*x[16+24*step]))
        grad[2, 17+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])))

    return result

def leg5_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LM leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[18+24*step]) +
        L_tibia*(np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) +
        L_femur*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) - x_body_mid_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[18+24*step]) + L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) +
        L_femur*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]),

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[19+24*step]) - L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 18+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[18+24*step]) - L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]))
        grad[0, 19+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]))
        grad[0, 20+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.cos(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 18+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[18+24*step]) + L_tibia*(np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) +
                              L_femur*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]))
        grad[1, 19+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step])) -
                              L_femur*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step]))
        grad[1, 20+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 18+24*step] = 0.0
        grad[2, 19+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) + L_femur*np.cos(deg2rad_const*x[19+24*step]))
        grad[2, 20+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])))

    return result

def leg6_end_point_XYZ(result, x, grad, step, px, py, pz):
    # This is LF leg
    result[:] = np.array([
        x[0+24*step] - px +
        L_coxa*np.cos(deg2rad_const*x[21+24*step]) + L_tibia*(np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) +
        L_femur*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) - x_body_front_half,

        x[1+24*step] - py +
        L_coxa*np.sin(deg2rad_const*x[21+24*step]) + L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) +
        L_femur*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) + y_body_half,

        x[2+24*step] - pz +
        L_femur*np.sin(deg2rad_const*x[22+24*step]) - L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))])

    if grad.size > 0:
        grad[:] = np.zeros([3, 24 * N])

        grad[0, 0+24*step] = 1.0
        grad[0, 21+24*step] = -deg2rad_const*(L_coxa*np.sin(deg2rad_const*x[21+24*step]) - L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]))
        grad[0, 22+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step])) -
                              L_femur*np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]))
        grad[0, 23+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.cos(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step])))

        grad[1, 1+24*step] = 1.0
        grad[1, 21+24*step] = deg2rad_const*(L_coxa*np.cos(deg2rad_const*x[21+24*step]) + L_tibia*(np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) +
                              L_femur*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]))
        grad[1, 22+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step])) -
                              L_femur*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step]))
        grad[1, 23+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step])))

        grad[2, 2+24*step] = 1.0
        grad[2, 21+24*step] = 0.0
        grad[2, 22+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) + L_femur*np.cos(deg2rad_const*x[22+24*step]))
        grad[2, 23+24*step] = deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])))

    return result

#------------------------------------ Constraint 4 -----------------------------------------------------------------
def leg1_end_point_Z(x, grad, step, pz):
    #This is RF leg, the algebra is compared with matlab to verify it's correctness. Should carryout actual number to double check it!
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[7+24*step]) - L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[7+24*step] = -deg2rad_const*(L_tibia * (np.cos(deg2rad_const*x[7+24*step]) * np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[8+24*step]) * np.sin(deg2rad_const*x[7+24*step])) + L_femur * np.cos(deg2rad_const*x[7+24*step]))
        grad[8+24*step] = -deg2rad_const*(L_tibia * (np.cos(deg2rad_const*x[7+24*step]) * np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[8+24*step]) * np.sin(deg2rad_const*x[7+24*step])))

    return result

def leg2_end_point_Z(x, grad, step, pz):
    # This is RM leg
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[10+24*step]) - L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[10+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) + L_femur*np.cos(deg2rad_const*x[10+24*step]))
        grad[11+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])))

    print("Converge diffenence along Z for leg2 pz=%f result=%f", pz, result)

    return result

def leg3_end_point_Z(x, grad, step, pz):
    # This is RR leg
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[13+24*step]) - L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[13+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) + L_femur*np.cos(deg2rad_const*x[13+24*step]))
        grad[14+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])))

    return result


def leg4_end_point_Z(x, grad, step, pz):
    # This is LR leg
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[16+24*step]) - L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[16+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) + L_femur*np.cos(deg2rad_const*x[16+24*step]))
        grad[17+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])))

    return result

def leg5_end_point_Z(x, grad, step, pz):
    # This is LM leg
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[19+24*step]) - L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[19+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) + L_femur*np.cos(deg2rad_const*x[19+24*step]))
        grad[20+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])))

    return result

def leg6_end_point_Z(x, grad, step, pz):
    # This is LF leg
    result = -(x[2+24*step] - pz + L_femur*np.sin(deg2rad_const*x[22+24*step]) - L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step])))

    if grad.size > 0:
        grad[:] = [0.0] * (24 * N)

        grad[2+24*step] = -1.0
        grad[22+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) + L_femur*np.cos(deg2rad_const*x[22+24*step]))
        grad[23+24*step] = -deg2rad_const*(L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])))

    return result

#------------------------------------ Constraint 5 -----------------------------------------------------------------
#Step 0-10 is first part
#Step 11-20 is the second part
def supportHexagonConstraint(result, x, grad, step):
    CoMx = x[step*24+3]
    CoMy = x[step*24+4]

    if (step <= (N-1)/2):
       #Polygon 1
       x1 = footstep_1_RF[0]
       y1 = footstep_1_RF[1]

       x2 = footstep_1_RR[0]
       y2 = footstep_1_RR[1]

       x3 = footstep_1_LM[0]
       y3 = footstep_1_LM[1]

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
       x1 = footstep_2_RM[0]
       y1 = footstep_2_RM[1]

       x2 = footstep_2_LR[0]
       y2 = footstep_2_LR[1]

       x3 = footstep_2_LF[0]
       y3 = footstep_2_LF[1]

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

gmSpeedLimit = abs(footstep_1_RF[1] - footstep_2_RF[1])/8.0 # Limit the speed
rotationalSpeedLimit = (0.52/np.pi*180.0) * (3.0 / ((N-1.0) / 4.0))   #Loaded speed 5rpm=0.52rad/s, let's say the whole process of lifting one leg is 1s, with 5 steps in between,
                                                  # then one step is 0.2s, give it 0.4s here. So the maximum rad rate is 0.52*0.4=0.21rad

def rotSpdLimUp(result, x, grad, step):
    #Note here step should range from 0 to N-2
    if grad.size > 0:
        grad[:] = np.zeros([18, 24 * N])
        for i in range(18):
            grad[i, 6+i+step*24] = 1.0
            grad[i, 6+i+(step+1)*24] = -1.0

    for i in range(18):
        result[i] = x[6+i+step*24] - x[6+i+(step+1)*24] - rotationalSpeedLimit
    print("rotational speed UP limit = %f %f %f" %(result[0], result[1], result[2]))
    return result

def rotSpdLimLow(result, x, grad, step):
    #This is the negative part of the rotational speed
    if grad.size > 0:
        grad[:] = np.zeros([18, 24 * N])
        for i in range(18):
            grad[i, 6+i+step*24] = -1.0
            grad[i, 6+i+(step+1)*24] = 1.0

    for i in range(18):
        result[i] = x[6+i+(step+1)*24] - x[6+i+step*24] - rotationalSpeedLimit
    print("rotational speed DOWN limit = %f %f %f" % (result[0], result[1], result[2]))
    return result

def GMSpdLimUp(result, x, grad, step):
    #Note here step should range from 0 to N-2
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
    # This is the negative part of the moving speed
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
        L_coxa*np.cos(deg2rad_const*x[6+24*step]) + L_tibia*(np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) +
        L_femur*np.cos(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[7+24*step]) + x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[6+24*step]) + L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step])*np.sin(deg2rad_const*x[8+24*step]) + np.sin(deg2rad_const*x[6+24*step])*np.cos(deg2rad_const*x[8+24*step])*np.sin(deg2rad_const*x[7+24*step])) +
        L_femur*np.cos(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[6+24*step]) + y_body_half,

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[7+24*step]) - L_tibia*(np.cos(deg2rad_const*x[7+24*step])*np.cos(deg2rad_const*x[8+24*step]) - np.sin(deg2rad_const*x[7+24*step])*np.sin(deg2rad_const*x[8+24*step]))])

    #This is RM leg
    resultRM = np.array([
        x[0+24*step] +
        L_coxa*np.cos(deg2rad_const*x[9+24*step]) + L_tibia*(np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) +
        L_femur*np.cos(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[10+24*step]) + x_body_mid_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[9+24*step]) + L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step])*np.sin(deg2rad_const*x[11+24*step]) + np.sin(deg2rad_const*x[9+24*step])*np.cos(deg2rad_const*x[11+24*step])*np.sin(deg2rad_const*x[10+24*step])) +
        L_femur*np.cos(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[9+24*step]),

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[10+24*step]) - L_tibia*(np.cos(deg2rad_const*x[10+24*step])*np.cos(deg2rad_const*x[11+24*step]) - np.sin(deg2rad_const*x[10+24*step])*np.sin(deg2rad_const*x[11+24*step]))])

    #This is RR leg
    resultRR = np.array([
        x[0+24*step] +
        L_coxa*np.cos(deg2rad_const*x[12+24*step]) + L_tibia*(np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) +
        L_femur*np.cos(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[13+24*step]) + x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[12+24*step]) + L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12+24*step])*np.sin(deg2rad_const*x[14+24*step]) + np.sin(deg2rad_const*x[12+24*step])*np.cos(deg2rad_const*x[14+24*step])*np.sin(deg2rad_const*x[13+24*step])) +
        L_femur*np.cos(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[12]) - y_body_half,

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[13+24*step]) - L_tibia*(np.cos(deg2rad_const*x[13+24*step])*np.cos(deg2rad_const*x[14+24*step]) - np.sin(deg2rad_const*x[13+24*step])*np.sin(deg2rad_const*x[14+24*step]))])

    #This is LR leg
    resultLR = np.array([
        x[0+24*step] +
        L_coxa*np.cos(deg2rad_const*x[15+24*step]) + L_tibia*(np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) +
        L_femur*np.cos(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[16+24*step]) - x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[15+24*step]) + L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step])*np.sin(deg2rad_const*x[17+24*step]) + np.sin(deg2rad_const*x[15+24*step])*np.cos(deg2rad_const*x[17+24*step])*np.sin(deg2rad_const*x[16+24*step])) +
        L_femur*np.cos(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[15+24*step]) - y_body_half,

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[16+24*step]) - L_tibia*(np.cos(deg2rad_const*x[16+24*step])*np.cos(deg2rad_const*x[17+24*step]) - np.sin(deg2rad_const*x[16+24*step])*np.sin(deg2rad_const*x[17+24*step]))])

    #This is LM leg
    resultLM = np.array([
        x[0+24*step] +
        L_coxa*np.cos(deg2rad_const*x[18+24*step]) + L_tibia*(np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) +
        L_femur*np.cos(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[19+24*step]) - x_body_mid_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[18+24*step]) + L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step])*np.sin(deg2rad_const*x[20+24*step]) + np.sin(deg2rad_const*x[18+24*step])*np.cos(deg2rad_const*x[20+24*step])*np.sin(deg2rad_const*x[19+24*step])) +
        L_femur*np.cos(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[18+24*step]),

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[19+24*step]) - L_tibia*(np.cos(deg2rad_const*x[19+24*step])*np.cos(deg2rad_const*x[20+24*step]) - np.sin(deg2rad_const*x[19+24*step])*np.sin(deg2rad_const*x[20+24*step]))])

    #This is LF leg
    resultLF = np.array([
        x[0+24*step] +
        L_coxa*np.cos(deg2rad_const*x[21+24*step]) + L_tibia*(np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) +
        L_femur*np.cos(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[22+24*step]) - x_body_front_half,

        x[1+24*step] +
        L_coxa*np.sin(deg2rad_const*x[21+24*step]) + L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step])*np.sin(deg2rad_const*x[23+24*step]) + np.sin(deg2rad_const*x[21+24*step])*np.cos(deg2rad_const*x[23+24*step])*np.sin(deg2rad_const*x[22+24*step])) +
        L_femur*np.cos(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[21+24*step]) + y_body_half,

        x[2+24*step] + L_femur*np.sin(deg2rad_const*x[22+24*step]) - L_tibia*(np.cos(deg2rad_const*x[22+24*step])*np.cos(deg2rad_const*x[23+24*step]) - np.sin(deg2rad_const*x[22+24*step])*np.sin(deg2rad_const*x[23+24*step]))])

    return [resultRF, resultRM, resultRR, resultLR, resultLM, resultLF]

#----------------------------------objective------------------------------------------------------------------------
W_coxa = 1
W_femur = 1
W_tibia = 1
W_diff = 10

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
tol = 1e-4
opt = nlopt.opt(nlopt.LD_SLSQP, 24*N)   #Make first 3 GC, second 3 CoM, and the rest motor angles

#opt = nlopt.opt(nlopt.GN_ISRES, 24*N)    #Improved Stochastic Ranking Evolution Strategy method

# opt = nlopt.opt(nlopt.AUGLAG, 24*N)           #Augmented lagrangian method
# local_opt = nlopt.opt(nlopt.LD_SLSQP, 24*N)   #Needed for augmented lagrangian method
# local_opt.set_xtol_rel(tol)                   #Needed for augmented lagrangian method
# opt.set_local_optimizer(local_opt)            #Needed for augmented lagrangian method

coxaUpperBound_RF = 90.0
coxaLowerBound_RF = 0.0
coxaUpperBound_RM = 45.0
coxaLowerBound_RM = -45.0
coxaUpperBound_RR = 0.0
coxaLowerBound_RR = -90.0
coxaUpperBound_LR = 270.0
coxaLowerBound_LR = 180.0
coxaUpperBound_LM = 225.0
coxaLowerBound_LM = 135.0
coxaUpperBound_LF = 180.0
coxaLowerBound_LF = 90.0

femurUpperBound = 90.0
femurLowerBound = -45.0
tibiaUpperBound = 45.0
tibiaLowerBound = -90.0

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

# Constraint 1: One 3-D vector equality constraints speaking c = function(theta, p)
for i in range(N):
    opt.add_equality_mconstraint(lambda result, x, grad, step=i: CoM(result, x, grad, step), np.array([tol, tol, tol]))


# Constraint 2: Three 3-D vector equality constraints speaking r0, rN, and rS
# opt.add_equality_mconstraint(r0, np.array([tol, tol, tol]))
# opt.add_equality_mconstraint(rS, np.array([tol, tol, tol]))
# opt.add_equality_mconstraint(rN, np.array([tol, tol, tol]))

#constraint 3: toepoints = specified toepoints on the ground (with limits on each motor angle), note before N/2 only first 3 toes, after N/2 only second 3 toes, at N/2 all 6 toes
#if step < (N-1)/2, then leg: right, right, left => leg1, leg3, leg5
#if step > (N-1)/2, then leg: left, left, right => leg2, leg4, leg6
#Note those end points are on the ground
px_RF_1 = footstep_1_RF[0]
py_RF_1 = footstep_1_RF[1]
pz_RF_1 = footstep_1_RF[2]

px_RF_2 = footstep_2_RF[0]
py_RF_2 = footstep_2_RF[1]
pz_RF_2 = footstep_2_RF[2]

px_RM_1 = footstep_1_RM[0]
py_RM_1 = footstep_1_RM[1]
pz_RM_1 = footstep_1_RM[2]

px_RM_2 = footstep_2_RM[0]
py_RM_2 = footstep_2_RM[1]
pz_RM_2 = footstep_2_RM[2]

px_RR_1 = footstep_1_RR[0]
py_RR_1 = footstep_1_RR[1]
pz_RR_1 = footstep_1_RR[2]

px_RR_2 = footstep_2_RR[0]
py_RR_2 = footstep_2_RR[1]
pz_RR_2 = footstep_2_RR[2]

px_LR_1 = footstep_1_LR[0]
py_LR_1 = footstep_1_LR[1]
pz_LR_1 = footstep_1_LR[2]

px_LR_2 = footstep_2_LR[0]
py_LR_2 = footstep_2_LR[1]
pz_LR_2 = footstep_2_LR[2]

px_LM_1 = footstep_1_LM[0]
py_LM_1 = footstep_1_LM[1]
pz_LM_1 = footstep_1_LM[2]

px_LM_2 = footstep_2_LM[0]
py_LM_2 = footstep_2_LM[1]
pz_LM_2 = footstep_2_LM[2]

px_LF_1 = footstep_1_LF[0]
py_LF_1 = footstep_1_LF[1]
pz_LF_1 = footstep_1_LF[2]

px_LF_2 = footstep_2_LF[0]
py_LF_2 = footstep_2_LF[1]
pz_LF_2 = footstep_2_LF[2]

#All points on the ground at step=0 and step=20
i = 0
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_1, py=py_RF_1, pz=pz_RF_1: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_1, py=py_RR_1, pz=pz_RR_1: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_1, py=py_LM_1, pz=pz_LM_1: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_1, py=py_RM_1, pz=pz_RM_1: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_1, py=py_LR_1, pz=pz_LR_1: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_1, py=py_LF_1, pz=pz_LF_1: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

i = N-1
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_2, py=py_RF_2, pz=pz_RF_2: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_2, py=py_RR_2, pz=pz_RR_2: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_2, py=py_LM_2, pz=pz_LM_2: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_2, py=py_RM_2, pz=pz_RM_2: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_2, py=py_LR_2, pz=pz_LR_2: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_2, py=py_LF_2, pz=pz_LF_2: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

groundOffset = 15.0
# 3 foot points on the ground
for i in range(1, (N-1)/2):   #This is 1-9 for N=21
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_1, py=py_RF_1, pz=pz_RF_1: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_1, py=py_RR_1, pz=pz_RR_1: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_1, py=py_LM_1, pz=pz_LM_1: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

    #z>30.0
    # Put 30.0 to be the safe offset from the ground
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=groundOffset: leg2_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=groundOffset: leg4_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=groundOffset: leg6_end_point_Z(x, grad, step, pz), tol)

i = (N-1)/2    #For N=20, we specify step 10 so that all legs are on the ground
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RF_1, py=py_RF_1, pz=pz_RF_1: leg1_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RR_1, py=py_RR_1, pz=pz_RR_1: leg3_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LM_1, py=py_LM_1, pz=pz_LM_1: leg5_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_2, py=py_RM_2, pz=pz_RM_2: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_2, py=py_LR_2, pz=pz_LR_2: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_2, py=py_LF_2, pz=pz_LF_2: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

for i in range((N-1)/2+1, N-1):   #This is 11-19
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_RM_2, py=py_RM_2, pz=pz_RM_2: leg2_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LR_2, py=py_LR_2, pz=pz_LR_2: leg4_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))
    opt.add_equality_mconstraint(lambda result, x, grad, step=i, px=px_LF_2, py=py_LF_2, pz=pz_LF_2: leg6_end_point_XYZ(result, x, grad, step, px, py, pz), np.array([tol, tol, tol]))

    #z>30.0, and for RF leg, z should be larger than step_height
    # Put 30.0 to be the safe offset from the ground
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=step_height: leg1_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=groundOffset: leg3_end_point_Z(x, grad, step, pz), tol)
    opt.add_inequality_constraint(lambda x, grad, step=i, pz=groundOffset: leg5_end_point_Z(x, grad, step, pz), tol)

#constraint 4: let's specify a triangular trajectory for the limbs on the air
liftHeight = 60.0  #Specify the liftheight
#Put as inequality constraint to be able to incorporate step height constraint
i = (N-1)/4      #This is step=5 for N=21
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg2_end_point_Z(x, grad, step, pz), tol)
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg4_end_point_Z(x, grad, step, pz), tol)
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg6_end_point_Z(x, grad, step, pz), tol)
i = 3*(N-1)/4    #This is step=15 for N=21
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg1_end_point_Z(x, grad, step, pz), tol)
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg3_end_point_Z(x, grad, step, pz), tol)
opt.add_inequality_constraint(lambda x, grad, step=i, pz=liftHeight: leg5_end_point_Z(x, grad, step, pz), tol)

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
opt.set_xtol_rel(tol)
#opt.set_xtol_abs(tol)
#opt.set_ftol_abs(tol)
#opt.set_ftol_rel(tol)  #This makes it stop too early!
maxtime = 600  #stop at 10 min
opt.set_maxtime(maxtime)
ini = [0.0]*(24*N)

RF_coxa_nominal = 60.0
RM_coxa_nominal = 0.0
RR_coxa_nominal = -60.0
LR_coxa_nominal = 240.0
LM_coxa_nominal = 180.0
LF_coxa_nominal = 120.0
femur_nominal = 30.0
tibia_nominal = -30.0
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

for i in range((N-1)/2):
    ini[0 + 24 * i] = r0x
    ini[1 + 24 * i] = r0y
    ini[2 + 24 * i] = r0z
    ini[3 + 24 * i] = r0x
    ini[4 + 24 * i] = r0y
    ini[5 + 24 * i] = r0z

for i in range((N-1)/2, N):
    ini[0 + 24 * i] = rNx
    ini[1 + 24 * i] = rNy
    ini[2 + 24 * i] = rNz
    ini[3 + 24 * i] = rNx
    ini[4 + 24 * i] = rNy
    ini[5 + 24 * i] = rNz

ini[0] = r0x
ini[1] = r0y
ini[2] = r0z

ini[24*((N-1)/2-1)+0] = rSx
ini[24*((N-1)/2-1)+1] = rSy
ini[24*((N-1)/2-1)+2] = rSz

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

motorKeyFrame = []
for i in range(N):
    print("Step", i)
    motorKeyFrame.append([[x[6 + i * 24],  x[7 + i * 24],  x[8 + i * 24]],
                          [x[9 + i * 24],  x[10 + i * 24], x[11 + i * 24]],
                          [x[12 + i * 24], x[13 + i * 24], x[14 + i * 24]],
                          [x[15 + i * 24], x[16 + i * 24], x[17 + i * 24]],
                          [x[18 + i * 24], x[19 + i * 24], x[20 + i * 24]],
                          [x[21 + i * 24], x[22 + i * 24], x[23 + i * 24]]])

    print("[[%.2f,   %.2f,   %.2f]," % (x[6 + i * 24],  x[7 + i * 24],  x[8 + i * 24]))
    print(" [%.2f,   %.2f,   %.2f]," % (x[9 + i * 24],  x[10 + i * 24], x[11 + i * 24]))
    print(" [%.2f,   %.2f,   %.2f]," % (x[12 + i * 24], x[13 + i * 24], x[14 + i * 24]))
    print(" [%.2f,   %.2f,   %.2f]," % (x[15 + i * 24], x[16 + i * 24], x[17 + i * 24]))
    print(" [%.2f,   %.2f,   %.2f]," % (x[18 + i * 24], x[19 + i * 24], x[20 + i * 24]))
    print(" [%.2f,   %.2f,   %.2f]]" % (x[21 + i * 24], x[22 + i * 24], x[23 + i * 24]))
    print("                                                                   ")

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