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
#
# Copyright 2021-2022 Johannes Laurin Hoermann
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
# # Sample case 1d cell

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
ref_domain_size_SI = 2e-8  # mM

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
ref_potential_difference

# %%
ref_domain_size

# %%
# reference data (SI units)

ref_dat_SI = np.loadtxt('samples/potential_sweep/data_std/NaCl_c_1000_1000_mM_z_+1_-1_l_20e-9_m_u_0.01_V.txt')
ref_x_SI = ref_dat_SI[:,0]
ref_potential_SI = ref_dat_SI[:,1]
ref_concentration_SI = ref_dat_SI[:,2:]

# %%
ref_concentration_SI.shape

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_x_SI, ref_potential_SI)
ax.set_ylabel(r'potential $\varphi\, (\mathrm{V})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_x_SI, ref_concentration_SI[:,0], ref_x_SI, ref_concentration_SI[:,1])
ax.set_ylabel(r'concentration $c\, (\mathrm{mM})$')
ax.set_xlabel(r'location $x\, (\mathrm{m})$')
fig.tight_layout()

# %%
# reference data (dimensionless)

ref_x = ref_x_SI / debye_length_SI
ref_potential = ref_potential_SI / thermal_voltage_SI
ref_concentration = ref_concentration_SI / average_concentration_SI

# %%
ref_concentration.shape

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_x, ref_potential)
ax.set_ylabel(r'dimensionless potential $\varphi$')
ax.set_xlabel(r'dimensionless location $x$')
fig.tight_layout()

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(ref_x, ref_concentration[:,0], ref_x, ref_concentration[:,1])
ax.set_ylabel(r'dimensionless concentration $c$')
ax.set_xlabel(r'dimensionless location $x$')
fig.tight_layout()
