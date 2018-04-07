#!/usr/bin/python

import gurobipy as go
import numpy as np

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
b = []
b.append(np.array([width_step,  0.0, \
				   length_before, 0.0, \
			       box_height, 0.0]))

b.append(np.array([width_step, 0.0, \
				   length_before+length_step+distance,  -(length_before+distance), \
			       box_height+step_height, -step_height]))

b.append(np.array([width_step, 0.0, \
				   length_before+length_step+2*distance+length_after, -(length_before+length_step+2*distance), \
			       box_height+2*step_height, -2*step_height]))

#Goal position
offset_after = 200  #The amound of offset that robot wants to stay after finishing the steps
x0 = np.array([ 750,  length_before+length_step+2*distance+offset_after+800, step_height*2.0,\
			   1250,  length_before+length_step+2*distance+offset_after+800, step_height*2.0,\
			   1500,  length_before+length_step+2*distance+offset_after+400, step_height*2.0,\
			   1250,  length_before+length_step+2*distance+offset_after,     step_height*2.0,\
			    750,  length_before+length_step+2*distance+offset_after,     step_height*2.0,\
			    500,  length_before+length_step+2*distance+offset_after+400, step_height*2.0])

print x0

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
			yy = np.array([M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k]), M[k] * (1 - H[i,j,k]), \
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
x = np.array([footsteps[N - 1, 0, 0], footsteps[N - 1, 0, 1], footsteps[N - 1, 0, 2], \
			  footsteps[N - 1, 1, 0], footsteps[N - 1, 1, 1], footsteps[N - 1, 1, 2], \
			  footsteps[N - 1, 2, 0], footsteps[N - 1, 2, 1], footsteps[N - 1, 2, 2], \
			  footsteps[N - 1, 3, 0], footsteps[N - 1, 3, 1], footsteps[N - 1, 3, 2], \
			  footsteps[N - 1, 4, 0], footsteps[N - 1, 4, 1], footsteps[N - 1, 4, 2], \
			  footsteps[N - 1, 5, 0], footsteps[N - 1, 5, 1], footsteps[N - 1, 5, 2]])

Q = np.diag([1.0, 1.0, 1.0, \
			 1.0, 1.0, 1.0, \
			 1.0, 1.0, 1.0, \
			 1.0, 1.0, 1.0, \
			 1.0, 1.0, 1.0, \
			 1.0, 1.0, 1.0 ])

#Quadratic expression for constraints on goal footstep position
dim1 = np.size(x)
expr1 = go.quicksum((x[i]-x0[i])*Q[i,j]*(x[j]-x0[j]) for i in range(dim1) for j in range(dim1))
#Quadratic expression for constraints on footstep length
weight1 = 0.01
expr2 = go.quicksum((footsteps[i,j,k]-footsteps[i-1,j,k])*weight1*(footsteps[i,j,k]-footsteps[i-1,j,k]) for i in range(1,N) for j in range(S) for k in range(3))
#Punishment on trim varaibles
weight2 = 10
expr3 = go.quicksum(weight2*trim[i] for i in range(N))
obj = expr1 + expr2 + expr3
print obj
m.setObjective(obj, go.GRB.MINIMIZE)

m.optimize()

for i in range(N):
	print("------------Round {} ------------------".format(i+1))
	print("Foot 1 [{}, {}, {}]".format(footsteps[i, 0, 0], footsteps[i, 0, 1], footsteps[i, 0, 2]))
	print("Foot 2 [{}, {}, {}]".format(footsteps[i, 1, 0], footsteps[i, 1, 1], footsteps[i, 1, 2]))
	print("Foot 3 [{}, {}, {}]".format(footsteps[i, 2, 0], footsteps[i, 2, 1], footsteps[i, 2, 2]))
	print("Foot 4 [{}, {}, {}]".format(footsteps[i, 3, 0], footsteps[i, 3, 1], footsteps[i, 3, 2]))
	print("Foot 5 [{}, {}, {}]".format(footsteps[i, 4, 0], footsteps[i, 4, 1], footsteps[i, 4, 2]))
	print("Foot 6 [{}, {}, {}]".format(footsteps[i, 5, 0], footsteps[i, 5, 1], footsteps[i, 5, 2]))

for i in range(N):
	print trim[i]

# for i in range(N):
# 	for j in range(S):
# 		print(H[i,j,0], H[i,j,1], H[i,j,2])