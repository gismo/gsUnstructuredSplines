#!/usr/bin/python

""""
    @file diss_geo_example.py

    @brief Print the mesh of the geometry and compute the exact solution with pylatex

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os
import sys

gismo_path = os.path.join(os.path.dirname(__file__), "../../../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

import matplotlib.pyplot as plt
import numpy as np

from pylatex import Document, NoEscape, TikZ, TikZOptions, Axis, Plot
import lib.diss_library as lib

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
# geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
geo_list = ["g1021", "g1121", "g1501", "g1311"]  # Without .xml extension
path_geo = "planar/geometries/"

# caption_list = ["Ex. I: geometry", "Ex. II: geometry", "Ex. III: geometry", "Ex. IV: geometry",
#                "Ex. I: solution", "Ex. II: solution", "Ex. III: solution", "Ex. IV: solution"]
caption_list = ["Ex. V: geometry", "Ex. VI: geometry", "Ex. VII: geometry", "Ex. VIII: geometry",
                "Ex. V: solution", "Ex. VI: solution", "Ex. VII: solution", "Ex. VIII: solution"]

#geo_list = ["g1000", "g1100", "g1510", "g1400","g1021", "g1121", "g1501", "g1311"]
#caption_list = geo_list + geo_list

numData = 50

ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)

compute_mesh = True
compute_solution = True
""" -------------------------------------------------------------------------------------------------- """

# Here the tikz files and the corresponding pdf are stored
# Change if necessary
path_dir = "geo/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

tikz_list = []
if compute_mesh:
    for geo in geo_list:
        # Making new folder if there exist no folder.
        if not os.path.exists("tikz_figures/" + path_dir):
            os.makedirs("tikz_figures/" + path_dir)

        mp = gs.core.gsMultiPatch()
        file = gs.io.gsFileData(path_geo + geo + ".xml")
        file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch

        v0 = np.zeros((2, numData))
        v0[0, :] = np.linspace(0, 1, numData)
        u0 = np.zeros((2, numData))
        u0[1, :] = np.linspace(0, 1, numData)

        v1 = np.ones((2, numData))
        v1[0, :] = np.flip(np.linspace(0, 1, numData))
        u1 = np.flip(np.ones((2, numData)))
        u1[1, :] = np.linspace(0, 1, numData)

        points_patchwise = []
        for idx in range(mp.nPatches()):
            points = mp.patch(idx).eval(v0)

            pointsB = mp.patch(idx).eval(u0)
            pointsC = mp.patch(idx).eval(v1)
            pointsD = mp.patch(idx).eval(u1)

            points_list = [pointsB, pointsC, pointsD]
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
        opt_patch = [{'fill': 'gray!20'}, {'fill': 'gray!40'}, {'fill': 'gray!60'}, {'fill': 'gray!80'}]
        opt_patch = opt_patch * 10

        opt_axis = TikZOptions(width=NoEscape(r'1\textwidth'))
        doc = Document()
        with doc.create(TikZ()) as tikz:
            with doc.create(Axis(options=opt_axis)) as axis:
                for idx in range(mp.nPatches()):
                    curve = Plot(options=opt_patch[idx], coordinates=points_patchwise[idx])
                    axis.append(curve)

        tikz_list.append(path_dir + geo + "_mesh")
        tex = doc.dumps()  # The document as string in LaTeX syntax
        with open(path_tikz + geo + "_mesh.tikz", 'w') as f:
            begin = False
            for idx, line in enumerate(tex.splitlines()):
                if line == "\\begin{tikzpicture}%":
                    begin = True
                if begin and idx != len(tex.splitlines()) - 1:
                    f.write(line)
                    f.write("\n")
else:
    for geo in geo_list:
        tikz_list.append(path_dir + geo + "_mesh")

if compute_solution:
    nx, ny = (500, 500)  # If you need a "smoother" picture, increase the number of points

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x, y)

    for geo in geo_list:
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        mp = gs.core.gsMultiPatch()
        file = gs.io.gsFileData(path_geo + geo + ".xml")
        file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch

        x = xv
        y = yv

        for idx in range(mp.nPatches()):
            Z = np.zeros((nx, ny))
            X = np.zeros((nx, ny))
            Y = np.zeros((nx, ny))
            for i in range(nx):
                for j in range(ny):
                    uv = [x[i, j], y[i, j]]
                    xy = mp.patch(idx).eval(uv)
                    X[i, j] = xy[0]
                    Y[i, j] = xy[1]
                    Z[i, j] = ms.eval(xy)

            surf = ax.plot_surface(X, Y, Z, edgecolor='none')
            surf._facecolors2d = surf._facecolor3d
            surf._edgecolors2d = surf._edgecolor3d

        tikz_list.append(path_dir + 'solution_' + geo + '.pdf')
        plt.savefig(path_fig + 'solution_' + geo + '.pdf')
        # plt.show()
else:
    for geo in geo_list:
        tikz_list.append(path_dir + 'solution_' + geo + '.pdf')

crop_list = []
for geo in geo_list:
    crop_list.append(path_dir + 'solution_' + geo)

# Creating the tex and pdf file
doc = lib.MyDocument()
doc.addTikzFigure(tikz_list, caption_list, row=4)
doc.generate_pdf("geo_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions(crop=crop_list)
doc.generate_pdf("geo_example", compiler="pdflatex", clean_tex=False)
lib.clean_extensions()
print("Finished: pdf saved to geo_example.pdf")
