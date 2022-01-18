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

# +
from dolfin import *
#from mshr import *
import numpy as np
import ufl
from time import gmtime, strftime
import meshio

msh = meshio.read("pointy.msh")

line_cells = []
for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data
    elif  cell.type == "line":
        if len(line_cells) == 0:
            line_cells = cell.data
        else:
            line_cells = np.vstack([line_cells, cell.data])

line_data = []
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "line":
        if len(line_data) == 0:
            line_data = msh.cell_data_dict["gmsh:physical"][key]
        else:
            line_data = np.vstack([line_data, msh.cell_data_dict["gmsh:physical"][key]])
    elif key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]

triangle_mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells})
line_mesh =meshio.Mesh(points=msh.points,cells=[("line", line_cells)],cell_data={"name_to_read":[line_data]})

meshio.write("pointy_mesh.xdmf", triangle_mesh)

meshio.xdmf.write("pointy_mf.xdmf", line_mesh)
# -


