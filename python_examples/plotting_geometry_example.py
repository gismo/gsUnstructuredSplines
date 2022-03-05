#!/usr/bin/python

""""
    @file plotting_geometry_example.py

    @brief Print the mesh of the geometry and compute the exact solution with pylatex

    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): P. Weinmüller
"""

import os
import sys

path_module = "data-private/"
print("Module path:", path_module, "(change if needed).")
os.chdir(os.path.join(os.path.dirname(__file__), "../../../"+path_module))

gismo_path = "../build/lib"
print("G+Smo path:", os.getcwd() + "/" + gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

import matplotlib.pyplot as plt
import numpy as np

from pylatex import Document, NoEscape, TikZ, TikZOptions, Axis, Plot
import library.library as lib

"""
    To run the diss_geo_example.py, we have the following options:
    
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
#geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
#geo_list = ["g1021", "g1121", "g1500", "g1311"]  # Without .xml extension
#geo_list = ["g1000", "g1100", "g1510", "g1400","g1021", "g1121", "g1501", "g1311"]
#geo_list = [f[:f.find(".xml")] for f in os.listdir("../filedata/") ]

from glob import glob


#path_geo = "../filedata/planar/"
path_geo = "../filedata/planar/geometries/"
path_geo_global = os.path.join(os.path.dirname(__file__), path_geo)
#path_geo = "../filedata/"
#path_geo_global = path_geo
geo_list = [y for x in os.walk(path_geo_global) for y in glob(os.path.join(x[0], '*.xml'))]
geo_list_temp = []
for geo in geo_list:
    geo_list_temp.append(geo[geo.find(path_geo)+len(path_geo):geo.find(".xml")])

geo_list = geo_list_temp
caption_list = geo_list

numData = 50

compute_mesh = True
""" -------------------------------------------------------------------------------------------------- """
# Here the tikz files and the corresponding pdf are stored
# Change if necessary
path_dir = "geo/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

tikz_list = []
if compute_mesh:
    for idgeo, geo in enumerate(geo_list):
        # Making new folder if there exist no folder.
        if not os.path.exists(path_tikz):
            os.makedirs(path_tikz)
        if not os.path.exists(path_fig):
            os.makedirs(path_fig)

        mp = gs.core.gsMultiPatch()
        file = gs.io.gsFileData(path_geo_global + geo + ".xml")
        if not file.getAnyFirst(mp):  # Assume that there exist only one gsMultiPatch
            caption_list.pop(idgeo)
            continue

        if not mp.domainDim() == 2:
            print("Surfaces/Volumes are not plotted yet!")
            caption_list.pop(idgeo)
            continue

        uv_list = []

        u0 = np.zeros((2, numData))
        u0[1, :] = np.linspace(0, 1, numData)
        uv_list.append(u0)
        v1 = np.ones((2, numData))
        v1[0, :] = np.flip(np.linspace(0, 1, numData))
        uv_list.append(v1)
        u1 = np.ones((2, numData))
        u1[1, :] = np.linspace(0, 1, numData)
        uv_list.append(u1)
        v0 = np.zeros((2, numData))
        v0[0, :] = np.linspace(0, 1, numData)
        uv_list.append(v0)

        points_patchwise = []
        for idx in range(mp.nPatches()):

            points_list = []
            for uv in uv_list:
                points_list.append(mp.patch(idx).eval(uv))
            points = points_list[0]
            points_list = points_list[1:]

            while len(points_list) > 0:
                for idx, p in enumerate(points_list):
                    if points[0, -1] == p[0, 0] and points[1, -1] == p[1, 0]:
                        points = np.append(points, p, axis=1)
                        points_list = points_list[:idx] + points_list[idx + 1:]
                        break
                    elif points[0, -1] == p[0, -1] and points[1, -1] == p[1, -1]:
                        points = np.append(points, np.flip(p, axis=1), axis=1)
                        points_list = points_list[:idx] + points_list[idx + 1:]
                        break

            points = np.transpose(points)
            points_patchwise.append(points)

        # Creating the tikz figures
        # Here we can change the layout setting of the mesh
        opt_patch = [{'fill': 'gray!20', 'line width': '1pt'}, {'fill': 'gray!40', 'line width': '1pt'},
                     {'fill': 'gray!60', 'line width': '1pt'}, {'fill': 'gray!80', 'line width': '1pt'}]
        opt_patch = opt_patch * 100

        opt_mesh = ['solid']

        opt_axis = TikZOptions(width=NoEscape(r'1\textwidth'))
        doc = Document()
        with doc.create(TikZ()) as tikz:
            with doc.create(Axis(options=opt_axis)) as axis:
                for idx in range(mp.nPatches()):
                    curve = Plot(options=opt_patch[idx], coordinates=points_patchwise[idx])
                    axis.append(curve)

                    # Adding the mesh line:
                    # Works only for uniform refinement
                    for el in range(mp.basis(idx).component(0).numElements()-1):
                        u1 = np.ones((2, numData)) * (el+1) / mp.basis(idx).component(0).numElements()
                        u1[1, :] = np.linspace(0, 1, numData)
                        points = mp.patch(idx).eval(u1)
                        points = points.T
                        curve = Plot(options=opt_mesh, coordinates=points)
                        axis.append(curve)

                    # Works only for uniform refinement
                    for el in range(mp.basis(idx).component(1).numElements()-1):
                        v1 = np.ones((2, numData)) * (el+1) / mp.basis(idx).component(1).numElements()
                        v1[0, :] = np.linspace(0, 1, numData)
                        points = mp.patch(idx).eval(v1)
                        points = points.T
                        curve = Plot(options=opt_mesh, coordinates=points)
                        axis.append(curve)

        geo_name = geo.replace("/",":")
        tikz_list.append(path_dir + geo_name + "_mesh")
        tex = doc.dumps()  # The document as string in LaTeX syntax
        with open(path_tikz + geo_name + "_mesh.tikz", 'w') as f:
            begin = False
            for idx, line in enumerate(tex.splitlines()):
                if line == "\\begin{tikzpicture}%":
                    begin = True
                if begin and idx != len(tex.splitlines()) - 1:
                    f.write(line)
                    f.write("\n")
else:
    for geo in geo_list:
        mp = gs.core.gsMultiPatch()
        file = gs.io.gsFileData(path_geo_global + geo + ".xml")
        if not file.getAnyFirst(mp):  # Assume that there exist only one gsMultiPatch
           continue

        geo_name = geo.replace("/",":")
        tikz_list.append(path_dir + geo_name + "_mesh")

if not len(tikz_list) == len(caption_list):
    print("Something is wrong in tikz_list and caption_list!")
    exit()

print(tikz_list)

# Creating the tex and pdf file
doc = lib.MyDocument()
doc.addTikzFigure(tikz_list, caption_list, col=4)
doc.generate_pdf("geo_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions()
print("Finished: pdf saved to geo_example.pdf")