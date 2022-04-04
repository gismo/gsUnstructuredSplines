#!/usr/bin/python

""""
    @file biharmonic_example.py

    @brief Compute biharmonic3_example using pygismo

    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os
import subprocess

path_module = "data/"
print("Module path:", path_module, "(change if needed).")
os.chdir(os.path.join(os.path.dirname(__file__), "../../../" + path_module))

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    DPatch = 1
    AlmostC1 = 2
    Nitsche = 3
    #Spline = 4  # Only for surface
    SurfASG1 = 5  # Only for surface

"""
    To run the biharmonic3_example.py, we have the following options:
    
    Input: - geo_list:          [list]       a list of geometries which we want to use
           - path_geo:          [string]     the path where the geometries are stored
           - NumRefinement:     [int]        the number of refinements
           - second:            [bool]       compute the second biharmonic problem, i.e., with Laplacian
           - deg_list           [list]       a list of integers of degrees which we want to use for each methods
           - method_list        [list]       a list of Method(Enum) which we want to use for each methods
           - compute_list       [list]       a list of booleans which methods we want to run
           
    Output:- The tikz files are stored in "tikz_files/" + domain + "geo/*"
           - A pdf and tex file is created with "geo_example.pdf" and "geo_example.tex"
           
"""
""" -------------------------------------------------------------------------------------------------- """
domain = "planar"  # "surfaces"

# Two-patch case
#geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension

# Multi-patch case
#geo_list = ["g1021", "g1121", "g1500", "g1311"]  # Without .xml extension

# Surface multi-patch
# geo_list = ["g2029", "g1021", "g2030", "g2013"] # Without .xml extension

# List for planar
geo_list = ["g1000", "g1100", "g1510", "g1400",
            "g1021", "g1121", "g1500", "g1311",
            "g1702", "g1705", "g1704", "g1023",
            "g1024", "g1122", "g1600", "g1701",
            "g1700", "g1123", "g1601", "g1801"]  # Without .xml extension

# List for surfaces
# geo_list = ["g2029", "g2024", "g2030", "g2013",
#             "g2007", "g2009", "g1031", "g1021",
#             "g2018", "g2000", "g2028", "g2008"]


path_geo = domain + "/geometries/"

NumRefinement = 4
second = False

deg_list = [
    [3, 4, 5],
    [3, 4, 5],
    [3, 4, 5],
    [3, 4],
]

# Approx C1: gluing data set to default: \tilde{p} = p-1, \tilde{r} = p-2,
# Nitsche: penalty set to default: via Eigenvalue-problem
method_list = [
    Method.ApproxC1,
    Method.DPatch,
    Method.Nitsche,
    Method.SurfASG1
]

compute_list = [
    False,
    False,
    True,
    False
]

path_example = "../build/bin/biharmonic3_example" if domain == "planar" else "../build/bin/biharmonic_surface2_example"
""" -------------------------------------------------------------------------------------------------- """
residual = True if domain == "surfaces" else False  # For surface residual computation

for idx, compute in enumerate(compute_list):
    if compute:
        for geo in geo_list:

            # Making new folder if there exist no folder.
            path_results = "results/" + domain + "/error/" + geo + "/"
            if not os.path.exists(path_results):
                os.makedirs(path_results)

            for deg in deg_list[idx]:
                m_str = ""
                reg = deg - 1
                if method_list[idx] == Method.ApproxC1:
                    m_str = "approxC1"
                elif method_list[idx] == Method.Nitsche:
                    m_str = "nitsche"
                elif method_list[idx] == Method.DPatch:
                    m_str = "dPatch"
                elif method_list[idx] == Method.AlmostC1:
                    m_str = "almostC1"
                elif method_list[idx] == Method.SurfASG1:
                    m_str = "surfASG1"
                    reg = deg - 2
                else:
                    print("METHOD NOT IMPLEMENTED!!!")

                argument_list = m_str + "-g" + geo + "-p" + str(deg) + "-s" + str(reg) + "-r" + str(NumRefinement) \
                                + "-m" + str(method_list[idx].value) + ("--second" if second else "") + ("--residual" if residual else "")

                # [!Run biharmonic2_example]
                proc = subprocess.Popen([path_example, "-g", geo, "-p", str(deg), "-s", str(reg), "-r", str(NumRefinement),
                                         "-m", str(method_list[idx].value), ("--second" if second else ""), "", ("--residual" if residual else ""), "",
                                         "-o", path_results + argument_list])
                proc.wait()
                # [!Run biharmonic2_example]
        print("Geometry: ", geo, " finished!")
print("Finished!")
