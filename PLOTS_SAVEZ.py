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
# Copyright 2022 Johannes Laurin Hoermann
#           2021 Maxim Kümmerle
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

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# %%
cartesian_plot = np.load("cart.npz")
rotational_plot = np.load("rot.npz")
onedimensional_plot = np.load("d.npz")

# %%
cartesian_plot["arr_1"]

# %%
# %matplotlib notebook
fig, host = plt.subplots()
fig.subplots_adjust(right=1.1)

p1, = host.plot(onedimensional_plot["arr_0"],onedimensional_plot["arr_1"],"-",label="Potential in 1D model")
p2, = host.plot(cartesian_plot["arr_0"],cartesian_plot["arr_1"],"--", label="Potential in cartesian model")
p3, = host.plot(rotational_plot["arr_0"],rotational_plot["arr_1"],"-.",label="Potential in cylindrical model")

host.set_xlim(0)
host.set_ylim(0, 0.5)
host.set_ylabel("$ \Phi\;[{U_T}]$")
host.set_xlabel("$ z_{dl}\;[{\lambda_D}]$")
plots = [p1,p2,p3]
host.legend(plots, [p.get_label() for p in plots], fontsize="small")
fig.tight_layout()
plt.savefig("ALL_POT.png", dpi=1200)


# %%
# %matplotlib notebook
fig, host = plt.subplots()
fig.subplots_adjust(right=1.1)

p1, = host.plot(onedimensional_plot["arr_0"],10000*onedimensional_plot["arr_2"],"-",label="Concentration in 1D model")
p2, = host.plot(cartesian_plot["arr_0"],10000*cartesian_plot["arr_2"],"--", label="Concentration in cartesian model")
p3, = host.plot(rotational_plot["arr_0"],10000*rotational_plot["arr_2"],"-.",label="Concentration in cylindrical model")

host.set_xlim(0)
#host.set_ylim(0.75, 1.35)
host.set_ylabel("$ c\;[c_{ref}]$")
host.set_xlabel("$ z_{dl}\;[{\lambda_D}]$")
plots = [p1,p2,p3]
host.legend(plots, [p.get_label() for p in plots], fontsize="small")
fig.tight_layout()
plt.savefig("ALL_CON.png", dpi=1200)



# %%
