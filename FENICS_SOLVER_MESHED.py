# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # FENICS SOLVER

# %% [markdown]
# ### Import needed packages 

# %%
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


# %%

# %%
mesh = fn.Mesh()
with fn.XDMFFile("meshing/mesh.xdmf") as infile:
    infile.read(mesh)
mvc = fn.MeshValueCollection("size_t", mesh, 2) 
with fn.XDMFFile("meshing/mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc)


# %%
M    = 2   #species

Poly    = fn.FiniteElement('Lagrange', mesh.ufl_cell(),2)
Multi   = fn.FiniteElement('Real', mesh.ufl_cell(), 0)
ElemP   = [Poly] * (M+1) 
ElemR   = [Multi] * (M)
Elem    = [ElemP + ElemR][0]
Mixed   = fn.MixedElement(Elem)
V       = fn.FunctionSpace(mesh, Mixed)

# %%
# define potentials and concentrations
u_GND  = fn.Expression('0', degree=2)          #Ground
u_DD   = fn.Expression('0.5', degree=2)          #pontential
#c_INIT = fn.Expression('0.01', degree=2)     #concentration on ground
c_avg  = fn.Expression('0.0001', degree=2)    #average concentration

# set boundary conditions
bcs = []
bcs += [fn.DirichletBC(V.sub(0), u_DD, mf, 5)]
bcs += [fn.DirichletBC(V.sub(0), u_GND, mf, 3)]
#bcs += [fn.DirichletBC(V.sub(i), c_INIT, mf, 1) for i in range(M+1)]


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
 

# %%
# Compute solution
fn.solve(PNP_xy == 0, UC, bcs) # solve function


# %%
S1, S2, S3, S4, S5 = UC.split()

# %%
fig = plt.figure()
ax = fig.add_subplot()

plt.plot(np.linspace(0,500),[S1(0,x,0) for x in np.linspace(0,500)])
plt.plot(np.linspace(0,600),[S1(100,x,0) for x in np.linspace(0,600)])
plt.plot(np.linspace(0,900),[S1(200,x,0) for x in np.linspace(0,900)])
plt.plot(np.linspace(0,1000),[S1(224,x,0) for x in np.linspace(0,1000)])

ax.set_ylabel("$ \Phi\;[\dfrac{1}{U_T}]$")
ax.set_xlabel("$ y_{dl}\;[\dfrac{l_{ref}}{\lambda_D}]$")
fig.tight_layout()
plt.savefig("POT.png")


# %%
import matplotlib.tri as tri

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:,0], xy[:,1], mesh.cells())

def plot(obj):
    plt.gca().set_aspect('equal')
    C = obj.compute_vertex_values(mesh)
    plt.tripcolor(mesh2triang(mesh), C, shading='flat')
    plt.tricontour(mesh2triang(mesh), C)
    plt.colorbar()

plot(S1)

# %%
plot(S2)

# %%
plot(S3)


# %%
def plot_tri(obj):
    plt.gca().set_aspect('equal')
    C = obj.compute_vertex_values(mesh)
    plt.tripcolor(mesh2triang(mesh), C, shading='gouraud')
    plt.tricontour(mesh2triang(mesh), C)
    plt.colorbar()


# %%
plot_tri(S1)

# %%
fig = plt.figure()
ax = fig.add_subplot()

plt.plot(0.0327683*np.linspace(223.9,500),[S1(x,1000,0) for x in np.linspace(223.9,500)])
plt.plot(0.0327683*np.linspace(100,500),[S1(x,600,0) for x in np.linspace(100,500)])
plt.plot(0.0327683*np.linspace(0,500),[S1(x,500,0) for x in np.linspace(0,500)])

ax.set_ylabel("$ \Phi\;[{U_T}]$")
ax.set_xlabel("$ x\;[{\lambda_D}]$")
fig.tight_layout()
plt.savefig("POTslahes_cart.png")

# %%
# %matplotlib notebook
fig, host = plt.subplots()
fig.subplots_adjust(right=1.1)

par1 = host.twinx()


p10, = host.plot(0.0327683*np.linspace(0,500),[S1(0,x,0) for x in np.linspace(0,500)],"-", color="#CC0000", label="Potential @ x = 0 $\lambda_D$")
p11, = host.plot(0.0327683*np.linspace(0,600),[S1(100,x,0) for x in np.linspace(0,600)],"--", color="#CC0000",label="Potential @ x = 3.25 $\lambda_D$")
p12, = host.plot(0.0327683*np.linspace(0,900),[S1(200,x,0) for x in np.linspace(0,900)], "-.",color="#CC0000", label="Potential @ x = 6.5 $\lambda_D$")
p13, = host.plot(0.0327683*np.linspace(0,1000),[S1(224,x,0) for x in np.linspace(0,1000)],".", color="#CC0000", label="Potential @ x = 7.33 $\lambda_D$")

p20, = par1.plot(0.0327683*np.linspace(0,500),[S2(0,x,0)*10000 for x in np.linspace(0,500)],"-", color="#0080FF", label="Concentration @ x = 0 $\lambda_D$")
p21, = par1.plot(0.0327683*np.linspace(0,600),[S2(100,x,0)*10000 for x in np.linspace(0,600)],"--", color="#0080FF",label="Concentration @ x = 3.25 $\lambda_D$")
p22, = par1.plot(0.0327683*np.linspace(0,900),[S2(200,x,0)*10000 for x in np.linspace(0,900)],"-.", color="#0080FF", label="Concentration @ x = 6.5 $\lambda_D$")
p23, = par1.plot(0.0327683*np.linspace(0,1000),[S2(224,x,0)*10000 for x in np.linspace(0,1000)],".", color="#0080FF", label="Concentration @ x = 7.33 $\lambda_D$")


host.set_xlim(0, 35)
host.set_ylim(0, 0.5)
par1.set_ylim(0.75, 1.35)


host.set_xlabel("Distance $ z\;[\lambda_D]$")
host.set_ylabel("Potential $ \Phi\;[U_T]$")
par1.set_ylabel("Concentration C $ [C_{ref}] $")


host.yaxis.label.set_color(p10.get_color())
par1.yaxis.label.set_color(p20.get_color())


tkw = dict(size=4, width=1.5)
host.tick_params(axis='y', colors=p10.get_color(), **tkw)
par1.tick_params(axis='y', colors=p20.get_color(), **tkw)

host.tick_params(axis='x', **tkw)

lines = [p10,p11,p12,p13,p20,p21,p22,p23]

host.legend(lines, [l.get_label() for l in lines], loc="center right", fontsize="x-small", bbox_to_anchor=(1, 0.5))
plt.savefig("CART_PLOT.png",bbox_inches = 'tight', dpi=1200)
plt.show()

# %%
x = 0.0327683*np.linspace(0,500)
y = [S1(0,x,0) for x in np.linspace(0,500)]
z = [S2(0,x,0) for x in np.linspace(0,500)]

# %%
np.savez("cart.npz", x, y, z)

# %%
