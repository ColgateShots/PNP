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
