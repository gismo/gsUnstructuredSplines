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

gismo_path = os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    Nitsche = 1
    DPatch = 2


def create_biharmonic_bvp(path_geo, loop, deg, reg, path_bvp, method, penalty_init=-1.0, cond=False):
    # [!Geometry]
    mp = gs.core.gsMultiPatch()
    file = gs.io.gsFileData(path_geo + ".xml")
    file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch
    mp.computeTopology(1e-4, False)

    # print(mp.nPatches())
    # boxSide = gs.core.boxSide(gs.core.side.west)
    # print(boxSide.index()) # get the side index
    # patchSide = gs.core.patchSide(1,boxSide)
    # print(patchSide.side().index()) # get the side index
    # print(patchSide.patch()) # get the patch index
    # print(mp.boundaries())
    # for bdy in mp.boundaries():
    #     print("Patch:", bdy.patch(), "Side:", bdy.side().index())
    #     print()

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
    opt.addSwitch("plot", "Plotting the results.", False)
    opt.addSwitch("info", "Plotting the results.", False)
    opt.addSwitch("cond","Estimate condition number (slow!)", cond)
    opt.addInt("refinementLoop", "Number of Uniform h-refinement loops.", loop)
    opt.addInt("discreteDegree", "Number of degree elevation steps to perform before solving (Degree 3 == 0).", deg)
    opt.addInt("discreteRegularity", "Number of degree elevation steps to perform before solving (Degree 3 == 0)", reg)
    opt.addInt("smoothing", "Which method should we use? (0 = Approx C1, 1 = Nitsche, 2 = DPatch)", method.value)
    opt.addReal("penalty", "Fixed Penalty value for Nitsche's method", penalty_init);
    # [!Option list]

    # [!Save the data to the XML-file]
    file = gs.io.gsFileData()
    file.add(bcs)  # id=0 Boundary
    file.add(f)  # id=1 Source function
    file.add(opt)  # id=2 Optionlist
    file.add(ms)  # id=3 Exact solution
    file.add(mp)  # id=X Geometry (should be last!)
    file.save(path_bvp, False)
    print("Filedata saved: " + path_bvp + ".xml")
    # [!Save the data to the XML-file]


''' ####### USER INPUT ####### '''
geo_list = ["g1501"]  # Without .xml extension
path_geo = "planar/geometries/"

deg_list = [3, 4, 5]
loop = 4

method = Method.ApproxC1
penalty_init = -1.0

# Condition number
cond = False

xml_col = "XmlCollection_input"  # Without .xml extension
path_xml = "results/XmlCollection/"
''' ##### USER INPUT END ##### '''

m_str = ""
if method == Method.ApproxC1:
    m_str = "approxC1"
elif method == Method.Nitsche:
    m_str = "nitsche"
elif method == Method.DPatch:
    m_str = "dPatch"

file = gs.io.gsXmlCollection(path_xml + xml_col + ".xml")
for geo in geo_list:
    for deg in deg_list:
        reg = deg - 1
        path_bvp = "results/" + geo + "/bvp/" + m_str + "-bvp1-p" + str(deg) + "-r" + str(reg) + "-l" + str(loop)
        create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=reg, path_bvp=path_bvp, method=method,
                              penalty_init=penalty_init, cond=cond)
        file.addFile(path_bvp + ".xml")

file.save()
print("Xml Collection saved: " + path_xml + xml_col + ".xml")
