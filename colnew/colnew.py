
import numpy as np
import scipy.special as special
import scikits.bvp1lg.colnew as colnew
#from scikits import colnew


nu = 3.4123
degrees = [2, 1]

def fsub(x, z):
	u, du, v = z     # it's neat to name the variables
	return np.array([-du/x + (nu**2/x**2 - 1)*u, x**(nu+1) * u])

def dfsub(x, z):
 	u, du, v = z
 	zero = np.zeros(x.shape)
 	return np.array([[(nu**2/x**2 - 1), -1/x, zero],[x**(nu+1), zero, zero]])


boundary_points = [1, 5, 10]

def gsub(z):
	u, du, v = z
	return np.array([u[0] - special.jv(nu,   1),v[1] - 5**(nu+1) * special.jv(nu+1, 5),u[2] - special.jv(nu,   10)])

def dgsub(z):
     return np.array([[1, 0, 0],
                     [0, 0, 1],
                     [1, 0, 0]])

tol = [1e-5, 0, 1e-5]

solution = colnew.solve(
     boundary_points, degrees, fsub, gsub,
     dfsub=dfsub, dgsub=dgsub,
     is_linear=True, tolerances=tol,
     vectorized=True, maximum_mesh_size=300)

solution.nmesh
x = np.linspace(1, 10, 101)
np.allclose(solution(x)[:,0], special.jv(nu, x),
            rtol=1e-4, atol=1e-8)
np.allclose(solution(x)[:,2], x**(nu+1)*special.jv(nu+1, x),
            rtol=1e-4, atol=1e-8)

import matplotlib.pyplot as plt
plt.plot(solution.mesh, solution(solution.mesh)[:,2], '.',
          x, x**(nu+1)*special.jv(nu+1, x), '-')

plt.show()
