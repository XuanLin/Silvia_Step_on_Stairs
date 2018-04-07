import gurobipy as go
import numpy as np

def setQuadraticObjective(mod, x, g, Q):
	dim = np.size(x)
	#not checking if g and Q have correspoinding size. Make sure that they have!
	s = go.quicksum((x[i]-g[i])*Q[i,j]*(x[j]-g[j]) for i in range(dim) for j in range(dim))
	mod.setObjective(s, go.GRB.MINIMIZE)

m = go.Model("Test")   

#Step 3
# foot1
p1x_3 = m.addVar(name="p1x_3")
p1y_3 = m.addVar(name="p1y_3")
p1z_3 = m.addVar(name="p1z_3")

# foot2
p2x_3 = m.addVar(name="p2x_3")
p2y_3 = m.addVar(name="p2y_3")
p2z_3 = m.addVar(name="p2z_3")

# foot3
p3x_3 = m.addVar(name="p3x_3")
p3y_3 = m.addVar(name="p3y_3")
p3z_3 = m.addVar(name="p3z_3")

m.update()

g = np.array([200, 1200, 0, 700, 1350, 0, 200, 1500, 0])   #Goal position
x = np.array([p1x_3, p1y_3, p1z_3, p2x_3, p2y_3, p2z_3, p3x_3, p3y_3, p3z_3])
Q = np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1])
setQuadraticObjective(m, x, g, Q)

m.optimize()



