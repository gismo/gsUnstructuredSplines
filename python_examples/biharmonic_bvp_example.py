#!/usr/bin/python

""""
    @file biharmonic_example.py

    @brief Compute biharmonic2_example using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys

path_module = "data/"
print("Module path:", path_module, "(change if needed).")
os.chdir(os.path.join(os.path.dirname(__file__), "../../../"+path_module))

gismo_path = "../build/lib"
print("G+Smo path:", os.getcwd() + "/" + gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

from enum import Enum

class Method(Enum):
    ApproxC1 =  0
    DPatch =    1
    AlmostC1 =  2
    Nitsche =   3


"""
    To run the biharmonic_bvp_example.py, we have the following options:
    
    Input: - geo_id:            [string]     A string of geometry which is stored in path_geo
           - path_geo:          [string]     The path where the geometries are stored
           - degree:            [int]        The polynomial degree for the discrete space
           - smoothness:        [int]        The smoothness for the discrete space
           - numRefine:         [int]        The number of refinements
           - method:            [Enum]       Which method should be used for solving the biharmonic problem
           
           - interpolation      [bool]       For approx C1: Compute the bf with an interpolation
           - gDdegree           [int]        For approx C1: The degree for the gluing data
           - gDsmoothness       [int]        For approx C1: The smoothness for the gluing data
           
           - penalty_init       [int]        For Nitsche's method: The initial value for Nitsche
           
           - cond               [bool]       Compute the condition number of the matrix
           - plot               [bool]       Print the result in paraview output
           - mesh               [bool]       Flag for printing the mesh lines in the output
           
    Output:- The xml file is stored in "data/biharmonic_bvp_example.xml"
    
    To run the xml file, use for this location:
    ~/gismo/extensions/gsUnstructuredSplines/python_examples $ 
    ../../../build/bin/biharmonic3_example -x "../../../data/biharmonic_bvp_example.xml"
    
           
"""
""" -------------------------------------------------------------------------------------------------- """
geo_id = "g1100"  # Without .xml extension
path_geo = "planar/geometries/"

degree = 3
smoothness = degree-1
numRefine = 2

method = Method.ApproxC1

# For Approx C1 method:
interpolation = True
gDdegree = degree-1
gDsmoothness = degree-2

# For Nitsche's method:
penalty_init = -1.0

cond = False  # Slow!

plot = False
mesh = False
""" -------------------------------------------------------------------------------------------------- """

# [!Geometry]
mp = gs.core.gsMultiPatch()
file = gs.io.gsFileData(path_geo + geo_id + ".xml")
file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch
mp.computeTopology(1e-4, False)
# [!Geometry]

# [!Right hand side]
f = gs.core.gsFunctionExpr("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)
# [!Right hand side]

# [!Exact solution]
ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
# [!Exact solution]

# [!Boundary]
dirichlet = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
neumann = gs.core.gsFunctionExpr(" -4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 " -4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2)
laplace = gs.core.gsFunctionExpr("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)

bcs = gs.pde.gsBoundaryConditions()
for bdy in mp.boundaries():
    #               patch_nr, side, boundary condition, function, unknown, parametric, component
    bcs.addCondition(bdy, gs.pde.bctype.dirichlet, dirichlet, 0, False, 0)
    bcs.addCondition(bdy, gs.pde.bctype.laplace, laplace, 0, False, 0)

# bcs.setGeoMap(mp)
# [!Boundary]

# [!Option list]
opt = gs.io.gsOptionList()
opt.addInt("numRefine", "Number of uniform h-refinement loops.", numRefine)
opt.addInt("degree", "Discrete polynomial degree.", degree)
opt.addInt("smoothness", "Discrete smoothness", smoothness)

opt.addInt("method", "Which method should be used for solving the biharmonic problem", method.value)

opt.addReal("penalty", "Fixed Penalty value for Nitsche's method", penalty_init)

opt.addSwitch("interpolation","Estimate condition number (slow!)", interpolation)
opt.addInt("gluingDataDegree", "Degree for the gluing data", gDdegree)
opt.addInt("gluingDataSmoothness", "Smoothness for the gluing data", gDsmoothness)

opt.addSwitch("cond","Estimate condition number (slow!)", cond)
opt.addSwitch("plot", "Plotting the results.", plot)
opt.addSwitch("mesh", "Plotting the mesh lines.", mesh)
# [!Option list]

# [!Save the data to the XML-file]
file = gs.io.gsFileData()
file.add(bcs)  # id=0 Boundary
file.add(f)  # id=1 Source function
file.add(opt)  # id=2 Optionlist
file.add(ms)  # id=3 Exact solution
file.add(mp)  # id=X Geometry (should be last!)
file.save("biharmonic_bvp_example", False)
print("Xml saved: biharmonic_bvp_example.xml")
# [!Save the data to the XML-file]
