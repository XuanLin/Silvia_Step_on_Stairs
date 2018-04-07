import gurobipy as go
import numpy as np

#1 needs go.GRB.something
#2 Open IDE to check grammar
#3 always row vector
#4 needs model.update()

def addNormConstraint(mod, x1, y1, z1, x2, y2, z2, th):
	mod.addQConstr(x1*x1+x2*x2-2*x1*x2 + y1*y1+y2*y2-2*y1*y2 + z1*z1+z2*z2-2*z1*z2 <= th^2)

def addPolytopeConstraint(mod, A, x, b):
	#Do not put a 1-D vector as column vector!
	dim1 = np.size(b)
	dim2 = np.size(x)
	for i in range(dim1):
		mod.addConstr((go.quicksum(A[i,j]*x[j] for j in range(dim2)) <= b[i]))


m = go.Model("Test")   #If import gurobipy then should be gurobipy.Model

a = m.addVar()
b = m.addVar()
c = m.addVar()
x = m.addVar()
y = m.addVar()
z = m.addVar()
H11 = m.addVar(vtype=go.GRB.BINARY, name="H11")
p1x_1 = m.addVar(name="p1x_1")
p1y_1 = m.addVar(name="p1y_1")
p1z_1 = m.addVar(name="p1z_1")

m.update()

A = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],  [0, 0, 1], [0, 0, -1]])
b1 = np.array([900, 0, 400, 1000, 10, 0])
M1 = 10000
yy = [M1*(1-H11), M1*(1-H11), M1*(1-H11), M1*(1-H11), M1*(1-H11), M1*(1-H11)]
bb = np.add(b1,yy)
xx = np.array([p1x_1, p1y_1, p1z_1])
addPolytopeConstraint(m, A, xx, bb)

#For feasibility, just don't set the objective value
#m.setObjective(a+b+c*c+x*x+y+z, GRB.MINIMIZE)
m.setObjective(a-b+c-x+y-z, go.GRB.MAXIMIZE)

addNormConstraint(m, a, b, c, x, y, z, 10)

m.optimize()

vars = m.getVars()
print('a=', vars[0])
print('b=', vars[1])
print('c=', vars[2])
print('x=', vars[3])
print('y=', vars[4])
print('z=', vars[5])



