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

import numpy as np

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    DPatch = 1
    AlmostC1 = 2
    Nitsche = 3


"""
    To run the biharmonic_test.py, we have the following options:
    
    Input: - geo_list:          [list]       a list of geometries which we want to use
           - path_geo:          [string]     the path where the geometries are stored
           - caption_list:      [list]       a list of strings for the captions in the figure
           - numData:           [int]        the number of data-points
           - ms                 [gismo]      a FunctionExpr which is the exact solution
           - compute_mesh       [bool]       if we want to compute the mesh
           - compute_solution   [bool]       if we want to compute the solution
           
    Output:- The tikz files are stored in "tikz_files/geo/*"
           - A pdf and tex file is created with "geo_example.pdf" and "geo_example.tex"
           
"""
""" -------------------------------------------------------------------------------------------------- """
domain = "planar"
#domain = "surfaces"

geo_list = ["g1000", "g1100", "g1510", "g1400"]
#geo_list = ["g1021", "g1121", "g1500", "g1311"]  # Without .xml extension
#geo_list = ["g1001", "g1021", "g1030", "g1031"]  # Without .xml extension
path_geo = domain + "/geometries/"

numRefinement = 5
degree = 4

second = False

N = 100

penalty = np.linspace(-10, 20, num=N)
penalty = np.around(np.power(2, penalty), decimals=5)

path_example = "../build/bin/biharmonic3_example" if domain == "planar" else "../build/bin/biharmonic_surface2_example"
""" -------------------------------------------------------------------------------------------------- """
residual = True if domain == "surfaces" else False  # For surface residual computation

path_dir = domain + "/nitsche/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

if not os.path.exists(path_tikz):
    os.makedirs(path_tikz)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)


for geo in geo_list:

    # Making new folder if there exist no folder.
    path_results = "results/" + path_dir + geo + "/"
    if not os.path.exists(path_results):
        os.makedirs(path_results)

    for pen in penalty:
        argument_list = "nitsche-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                        + "-m" + str(Method.Nitsche.value) + ("--second" if second else "") + ("--residual" if residual else "") \
                        + "-y" + str(pen)

        # [!Run biharmonic2_example]
        proc = subprocess.Popen([path_example, "-g", geo, "-p", str(degree), "-s", str(degree - 1), "-r", str(numRefinement),
                                 "-m", str(Method.Nitsche.value), ("--second" if second else ""), "", ("--residual" if residual else ""), "",
                                 "-y", str(pen), "-o", path_results + argument_list])
        proc.wait()
        # [!Run biharmonic2_example]

    print("Geometry: ", geo, " finished!")

    # Approx C1
    argument_list = "approxC1-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                    + "-m" + str(Method.ApproxC1.value) + ("--second" if second else "") + ("--residual" if residual else "")

    # [!Run biharmonic2_example]
    proc = subprocess.Popen([path_example, "-g", geo, "-p", str(degree), "-s", str(degree - 1), "-r", str(numRefinement),
                             "-m", str(Method.ApproxC1.value), ("--second" if second else ""), "", ("--residual" if residual else ""), "",
                             "-y", str(pen), "-o", path_results + argument_list])
    proc.wait()
    # [!Run biharmonic2_example]

    # Nitsche with EW
    argument_list = "nitsche-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                    + "-m" + str(Method.Nitsche.value) + ("--second" if second else "") + ("--residual" if residual else "")\
                    + "-y" + str(-1)

    # [!Run biharmonic2_example]
    proc = subprocess.Popen([path_example, "-g", geo, "-p", str(degree), "-s", str(degree - 1), "-r", str(numRefinement),
                             "-m", str(Method.Nitsche.value), ("--second" if second else ""), "", ("--residual" if residual else ""), "",
                             "-y", str(-1), "-o", path_results + argument_list])
    proc.wait()
    # [!Run biharmonic2_example]
    print("Geometry: ", geo, " finished!")
print("Finished!")
