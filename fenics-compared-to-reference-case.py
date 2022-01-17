#
# Copyright 2021 Johannes Laurin Hoermann
#
# ### MIT license
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

# %% [markdown]
# # FENICS SOLVER compared with reference case
#

# %%
# if juyter notebook autocompletion won't work, this might help
# %config Completer.use_jedi = False

# %%
# Get the libraries
import fenics as fn
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# %% [markdown]
# ## reference
# (from https://github.com/libAtoms/matscipy/tree/master/examples/electrochemistry/pnp_batch/cell_1d)

# %%
import scipy.constants as const

# %%
# SI units
ref_potential_difference_SI = 0.01  # V
ref_concentration_SI = np.array([1000, 1000])  # mM or mol/m^3
ref_number_charges = np.array([1,-1]) 
ref_domain_size_SI = 2e-8  # m

average_concentration_SI = 1000
thermal_voltage_SI = const.Boltzmann * 298 / const.elementary_charge
ionic_strength_SI = 0.5*np.sum(ref_concentration_SI*np.square(ref_number_charges))
debye_length_SI = np.sqrt(82*const.epsilon_0*const.Boltzmann*298/(2*const.Avogadro*const.elementary_charge**2*ionic_strength_SI))

# %%
thermal_voltage_SI

# %%
ionic_strength_SI

# %%
debye_length_SI

# %%
# dimensionless
ref_concentration = ref_concentration_SI / average_concentration_SI
ref_potential_difference = ref_potential_difference_SI / thermal_voltage_SI
ref_domain_size = ref_domain_size_SI / debye_length_SI

# %%
ref_concentration

# %%
ref_potential_difference

# %%
ref_domain_size

# %%
# reference data (SI units)

ref_dat_SI = np.loadtxt('samples/potential_sweep/data_std/NaCl_c_1000_1000_mM_z_+1_-1_l_20e-9_m_u_0.01_V.txt')
ref_dat_x_SI = ref_dat_SI[:,0]
ref_dat_potential_SI = ref_dat_SI[:,1]
ref_dat_concentration_SI = ref_dat_SI[:,2:]

# %%
ref_dat_concentration_SI.shape

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_dat_x_SI, ref_dat_potential_SI)
ax.set_ylabel(r'potential $\varphi\, (\mathrm{V})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_dat_x_SI, ref_dat_concentration_SI[:,0], ref_dat_x_SI, ref_dat_concentration_SI[:,1])
ax.set_ylabel(r'concentration $c\, (\mathrm{mM})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %%
# reference data (dimensionless)

ref_dat_x = ref_dat_x_SI / debye_length_SI
ref_dat_potential = ref_dat_potential_SI / thermal_voltage_SI
ref_dat_concentration = ref_dat_concentration_SI / average_concentration_SI

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_dat_x, ref_dat_potential)
ax.set_ylabel(r'potential $\varphi\, (\mathrm{V})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_dat_x, ref_dat_concentration[:,0], ref_dat_x, ref_dat_concentration[:,1])
ax.set_ylabel(r'concentration $c\, (\mathrm{mM})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %% [markdown]
# ## minimalistic 

# %%
# define mesh and define function space
X    = 500  #x-limit
Y    = ref_domain_size / 2  #y-limit
NX   = 50  #x-steps
NY   = 50  #y-steps

mesh = fn.RectangleMesh(fn.Point(-X, -Y), fn.Point(X, Y), NX, NY)

Poly    = fn.FiniteElement('Lagrange', mesh.ufl_cell(),2)
Real    = fn.FiniteElement('Real', mesh.ufl_cell(), 0)
Elem    = [Poly,Poly,Poly,Real,Real] 
Mixed   = fn.MixedElement(Elem)
V       = fn.FunctionSpace(mesh, Mixed)

# %%
Elem

# %%
# define potentials and concentrations
u_GND  = fn.Constant(ref_potential_difference/2) #fn.Expression('0', degree=2)
u_DD   = fn.Constant(-ref_potential_difference/2) #fn.Expression('0.025', degree=2)

# TODO: treat reference cocentrations species-wise 
c_INIT = fn.Constant(ref_concentration[0]) #fn.Expression('0.00001', degree=2)
c_AVG  = fn.Constant(ref_concentration[0]) #fn.Expression('0.00001', degree=2)

# define boundaries
def boundaryGND(x, on_boundary):
    tol=1e-12
    # return ((x[1] < -400 + tol)) 
    return on_boundary and fn.near(x[1], -ref_domain_size/2, tol)
def boundaryHigh(x, on_boundary):
    tol=1e-12
    # return (x[1] > 400  - tol)
    return on_boundary and fn.near(x[1], ref_domain_size/2, tol)


# %%
# set boundary conditions
GammaGND  = fn.DirichletBC(V.sub(0), u_GND, boundaryGND)  # ground potential at straight electrode
GammaHigh = fn.DirichletBC(V.sub(0), u_DD, boundaryHigh)  # high potential at shaped electrode
#GammaC_GND0 = fn.DirichletBC(V.sub(0) , c_INIT, boundaryGND)
#GammaC_GND1 = fn.DirichletBC(V.sub(1) , c_INIT, boundaryGND) 
#GammaC_GND2 = fn.DirichletBC(V.sub(2) , c_INIT, boundaryGND)
bcs=[GammaGND,GammaHigh]
#bcs=[GammaGND,GammaHigh,GammaC_GND0,GammaC_GND1,GammaC_GND2]

# %%
# define problem
UC    = fn.Function(V)
uc    = fn.split(UC)
u, c1, c2, lam1, lam2  = uc[0], uc[1], uc[2], uc[3], uc[4]

VW    = fn.TestFunctions(V)                    
v, w1, w2, mu1, mu2  = VW[0], VW[1], VW[2], VW[3], VW[4]

#create rotation
#r = fn.Expression('x[0]', degree=1)

# changing concentrations charges
PoissonLeft     = (fn.dot(fn.grad(u), fn.grad(v)))*fn.dx
PoissonRight    = c1*v*fn.dx - c2*v*fn.dx
NernstPlanck1   = fn.dot((-fn.grad(c1) - (c1)*fn.grad(u)),fn.grad(w1))*fn.dx
NernstPlanck2   = fn.dot((-fn.grad(c2) + (c2)*fn.grad(u)),fn.grad(w2))*fn.dx

Constraint1     = lam1 * w1 * fn.dx + ((c1) - c_AVG) * mu1 * fn.dx
Constraint2     = lam2 * w2 * fn.dx + ((c2) - c_AVG) * mu2 * fn.dx 

PNP_xy          = PoissonLeft - PoissonRight + NernstPlanck1 + NernstPlanck2 + Constraint1 + Constraint2

# %%
# Compute solution
fn.solve(PNP_xy == 0, UC, bcs) # solve function

# %%
# # %matplotlib notebook
x,y=mesh.coordinates().T
X=x.reshape(NX+1,NY+1)
Y=y.reshape(NX+1,NY+1)

# %%
w_xy = np.array([UC(xy) for xy in mesh.coordinates()])
W_xy = w_xy.reshape((NX+1),(NY+1),5)
U_xy = W_xy[:,:,0]
C_xys = W_xy[:,:,1:]

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(Y[::,10],U_xy[::,10])
#ax.plot(Y,U_xy[::,20])
#ax.plot(Y,U_xy[::,25])

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(Y[::,10],C_xys[:,:,0][::,10])
ax.plot(Y[::,20],C_xys[:,:,0][::,20])
ax.plot(Y[::,15],C_xys[:,:,0][::,25])
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(Y[::,10],C_xys[:,:,1][::,10])
ax.plot(Y[::,20],C_xys[:,:,1][::,20])
ax.plot(Y[::,25],C_xys[:,:,1][::,25])

# %%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,U_xy)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,C_xys[:,:,0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,C_xys[:,:,1])


# %% [markdown]
# ## comparsion

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_dat_x-ref_domain_size/2, ref_dat_potential, 'b:', label='reference 1d')
ax.plot(Y[::,10], U_xy[::,10], 'r--', label='fenics 2d')
ax.set_ylabel(r'potential $\varphi$, dimensionless')
ax.set_xlabel(r'location $x$, dimensionless')
ax.legend()
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.semilogy(ref_dat_x-ref_domain_size/2, ref_dat_concentration[:,0], 
            color='moccasin', linestyle='-', linewidth=2, label='cation, reference case')
ax.semilogy(ref_dat_x-ref_domain_size/2, ref_dat_concentration[:,1], 
            color='skyblue', linestyle='-', linewidth=2, label='anion, reference case')

ax.semilogy(Y[::,10], C_xys[:,:,0][::,20], color='tab:orange', linestyle=':', label='cation, fenics 2d cross section')
ax.semilogy(Y[::,10], C_xys[:,:,1][::,20], color='tab:blue', linestyle=':', label='anion, fenics 2d cross section')

# ref_dat_x, ref_dat_concentration[:,1])
ax.set_ylabel(r'concentration $c$, dimensionless')
ax.set_xlabel(r'location $x$, dimensionless')

ax.legend() 
fig.tight_layout()

# %%
