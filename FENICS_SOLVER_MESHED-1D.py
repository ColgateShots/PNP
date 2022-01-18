# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # FENICS SOLVER

# ### Import needed packages 

# Get the libraries
import fenics as fn
import numpy as np
import sympy as sym
import scipy as sc
from scipy import constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm #Colormap
import meshio as mio
#import mshr as msr




mesh = fn.Mesh()
with fn.XDMFFile("meshing/Dmesh.xdmf") as infile:
    infile.read(mesh)
mvc = fn.MeshValueCollection("size_t", mesh, 2) 
with fn.XDMFFile("meshing/Dmf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc)


# +
M    = 2   #species

Poly    = fn.FiniteElement('Lagrange', mesh.ufl_cell(),2)
Multi   = fn.FiniteElement('Real', mesh.ufl_cell(), 0)
ElemP   = [Poly] * (M+1) 
ElemR   = [Multi] * (M)
Elem    = [ElemP + ElemR][0]
Mixed   = fn.MixedElement(Elem)
V       = fn.FunctionSpace(mesh, Mixed)

# +
# define potentials and concentrations
u_GND  = fn.Expression('0', degree=2)          #Ground
u_DD   = fn.Expression('0.5', degree=2)          #pontential
c_avg  = fn.Expression('0.0001', degree=2)    #average concentration

# set boundary conditions
bcs = []
bcs += [fn.DirichletBC(V.sub(0), u_DD, mf,3)]
bcs += [fn.DirichletBC(V.sub(0), u_GND, mf, 1)]

# define problem
UC    = fn.Function(V)
uc    = fn.split(UC)                        # trial function potential concentration lagrange multi
u, c, lam = uc[0], uc[1:M+1], uc[M+1:]

VW    = fn.TestFunctions(V)                          # test function potential concentration lagrange multi                     
v, w, mu = VW[0], VW[1:M+1], VW[M+1:]


#lets try rot
r = fn.Expression('x[0]', degree=0)

# changing concentrations charges
Rho = 0
for i in range(M):
    if i%2:
        Rho += -c[i]
    else:
        Rho += c[i]

PoissonLeft     = (fn.dot(fn.grad(u), fn.grad(v)))*fn.dx                    # weak solution Poisson left
PoissonRight    = -(Rho)*v*fn.dx                                  # weak solution Poisson right
NernstPlanck    = 0
for i in range(M):
    if i%2:
        NernstPlanck += fn.dot((-fn.grad(c[i]) + c[i]*fn.grad(u)),fn.grad(w[i]))*fn.dx     # weak solution Nernst-Planck 
    else:
        NernstPlanck += fn.dot((-fn.grad(c[i]) - c[i]*fn.grad(u)),fn.grad(w[i]))*fn.dx     # weak solution Nernst-Planck

constraint = 0
for i in range(M):
    constraint += lam[i] * w[i] * fn.dx + (c[i] - c_avg) * mu[i] * fn.dx            #constraint a la hoermann
    
        
PNP_xy     = PoissonLeft + PoissonRight + NernstPlanck + constraint        # PNP system
 
# -

# Compute solution
fn.solve(PNP_xy == 0, UC, bcs) # solve function


S1, S2, S3, S4, S5 = UC.split()

# +
fig = plt.figure()
ax = fig.add_subplot()

plt.plot(np.linspace(0,500),[S1(0,x,0) for x in np.linspace(0,500)])

ax.set_ylabel("$ \Phi\;[{U_T}]$")
ax.set_xlabel("$ z_{dl}\;[{\lambda_D}]$")
fig.tight_layout()
plt.savefig("1DPOT.png")

# -

x = 0.0327683*np.linspace(0,500)
y = [S1(0,x,0) for x in np.linspace(0,500)]
z = [S2(0,x,0) for x in np.linspace(0,500)]

np.savez("d.npz", x, y, z)


