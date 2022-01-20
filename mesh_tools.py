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
#           2021 ColgateShots
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
from dolfin import *
#from mshr import *
import numpy as np
import ufl
from time import gmtime, strftime
import meshio
from gmsh import *
msh = meshio.read("pointy.msh")

"""
# 3D #############################
triangle_cells = []
for cell in msh.cells:
    if cell.type == "tetra":
        tetra_cells = cell.data
    elif  cell.type == "triangle":
        if len(triangle_cells) == 0:
            triangle_cells = cell.data
        else:
            triangle_cells = np.vstack([triangle_cells, cell.data])

triangle_data = []
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        if len(triangle_data) == 0:
            triangle_data = msh.cell_data_dict["gmsh:physical"][key]
        else:
            triangle_data = np.vstack([triangle_data, msh.cell_data_dict["gmsh:physical"][key]])
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})
triangle_mesh =meshio.Mesh(points=msh.points,
                           cells=[("triangle", triangle_cells)],
                           cell_data={"name_to_read":[triangle_data]})
meshio.write("mesh_3D.xdmf", tetra_mesh)

meshio.xdmf.write("mf_3D.xdmf", triangle_mesh)
################################
"""

# 2D ###########################
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
line_mesh = meshio.Mesh(points=msh.points,cells=[("line", line_cells)],cell_data={"name_to_read":[line_data]})
meshio.write("mesh_2D.xdmf", triangle_mesh)

meshio.xdmf.write("mf_2D.xdmf", line_mesh)
###############################
""""
mesh = Mesh()
with XDMFFile("mesh_2D.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2) 
with XDMFFile("mf_2D.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)


bc1 = (DirichletBC(W.sub(1), Constant((0.0, 0.0, 0.0)), mf, 1)) # so z.B. Randbedingung für verschiedene Ränder zuweisen über mf
bc2 = (DirichletBC(W.sub(1), Constant((0.0, 0.0, 0.0)), mf, 2))
"""
