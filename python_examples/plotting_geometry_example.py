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
os.chdir(os.path.join(os.path.dirname(__file__), "../../../" + path_module))

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
domain = "planar"
#domain = "surfaces"

# geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
# geo_list = ["g1021", "g1121", "g1500", "g1311"]  # Without .xml extension
# geo_list = ["g1000", "g1100", "g1510", "g1400","g1021", "g1121", "g1501", "g1311"]
# geo_list = ["g1121", "g1702", "g1704", "g1703"]
geo_list = ["g1001", "g1021", "g1030", "g1031"]  # Without .xml extension

from glob import glob

path_geo = "../filedata/" + domain + "/geometries/"

# path_geo_global = os.path.join(os.path.dirname(__file__), path_geo)
# geo_list = [y for x in os.walk(path_geo_global) for y in glob(os.path.join(x[0], '*.xml'))]
# geo_list_temp = []
# for geo in geo_list:
#     geo_list_temp.append(geo[geo.find(path_geo)+len(path_geo):geo.find(".xml")])
# geo_list = geo_list_temp
# print(geo_list)

caption_list = geo_list

# numData = 50
numData = 50
""" -------------------------------------------------------------------------------------------------- """
# Here the tikz files and the corresponding pdf are stored
# Change if necessary
path_dir = domain + "/geo/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

tikz_list = []
crop_list = []
for idgeo, geo in enumerate(geo_list):
    # Making new folder if there exist no folder.
    if not os.path.exists(path_tikz):
        os.makedirs(path_tikz)
    if not os.path.exists(path_fig):
        os.makedirs(path_fig)

    mp = gs.core.gsMultiPatch()
    file = gs.io.gsFileData(domain + '/geometries/' + geo + ".xml")
    if not file.getAnyFirst(mp):  # Assume that there exist only one gsMultiPatch
        caption_list.pop(idgeo)
        continue

    # if not mp.domainDim() == 2:
    #     print("Surfaces/Volumes are not plotted yet!")
    #     caption_list.pop(idgeo)
    #     continue
    if mp.targetDim() == 2:
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
                    if mp.targetDim() == 2:
                        if points[0, -1] == p[0, 0] and points[1, -1] == p[1, 0]:
                            points = np.append(points, p, axis=1)
                            points_list = points_list[:idx] + points_list[idx + 1:]
                            break
                        elif points[0, -1] == p[0, -1] and points[1, -1] == p[1, -1]:
                            points = np.append(points, np.flip(p, axis=1), axis=1)
                            points_list = points_list[:idx] + points_list[idx + 1:]
                            break
                    elif mp.targetDim() == 3:
                        if points[0, -1] == p[0, 0] and points[1, -1] == p[1, 0] and points[2, -1] == p[2, 0]:
                            points = np.append(points, p, axis=1)
                            points_list = points_list[:idx] + points_list[idx + 1:]
                            break
                        elif points[0, -1] == p[0, -1] and points[1, -1] == p[1, -1] and points[2, -1] == p[2, -1]:
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
                    for el in range(mp.basis(idx).component(0).numElements() - 1):
                        u1 = np.ones((2, numData)) * (el + 1) / mp.basis(idx).component(0).numElements()
                        u1[1, :] = np.linspace(0, 1, numData)
                        points = mp.patch(idx).eval(u1)
                        points = points.T
                        curve = Plot(options=opt_mesh, coordinates=points)
                        axis.append(curve)

                    # Works only for uniform refinement
                    for el in range(mp.basis(idx).component(1).numElements() - 1):
                        v1 = np.ones((2, numData)) * (el + 1) / mp.basis(idx).component(1).numElements()
                        v1[0, :] = np.linspace(0, 1, numData)
                        points = mp.patch(idx).eval(v1)
                        points = points.T
                        curve = Plot(options=opt_mesh, coordinates=points)
                        axis.append(curve)

        geo_name = geo.replace("/", ":")
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

    elif mp.targetDim() == 3:
        nx, ny = (numData, numData)  # If you need a "smoother" picture, increase the number of points

        x = np.linspace(0, 1, nx)
        y = np.linspace(0, 1, ny)
        xv, yv = np.meshgrid(x, y)

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        x = xv
        y = yv
        for idx in range(mp.nPatches()):
            Z = np.zeros((nx, ny))
            X = np.zeros((nx, ny))
            Y = np.zeros((nx, ny))
            for i in range(nx):
                for j in range(ny):
                    uv = [x[i, j], y[i, j]]
                    xyz = mp.patch(idx).eval(uv)
                    X[i, j] = xyz[0]
                    Y[i, j] = xyz[1]
                    Z[i, j] = xyz[2]

            surf = ax.plot_surface(X, Y, Z, edgecolor='none')
            surf._facecolors2d = surf._facecolor3d
            #surf._edgecolors2d = surf._edgecolor3d

            # add the mesh line
            # Works only for uniform refinement
            for el in range(mp.basis(idx).component(0).numElements() - 1):
                u1 = np.ones((2, numData)) * (el + 1) / mp.basis(idx).component(0).numElements()
                u1[1, :] = np.linspace(0, 1, numData)
                points = mp.patch(idx).eval(u1)
                points = points.T

                X = points[:,0]
                Y = points[:,1]
                Z = points[:,2]
                ax.plot3D(X, Y, Z, 'b')

            for el in range(mp.basis(idx).component(1).numElements() - 1):
                v1 = np.ones((2, numData)) * (el + 1) / mp.basis(idx).component(1).numElements()
                v1[0, :] = np.linspace(0, 1, numData)
                points = mp.patch(idx).eval(v1)
                points = points.T

                X = points[:,0]
                Y = points[:,1]
                Z = points[:,2]
                ax.plot3D(X, Y, Z,'b')

        geo_name = geo.replace("/", ":")
        tikz_list.append(path_dir + geo_name + "_mesh" + '.pdf')
        plt.savefig(path_fig + geo_name + "_mesh" + '.pdf')
        crop_list.append(path_dir + geo_name + "_mesh")
        #plt.show()
        plt.close(fig)

if not len(tikz_list) == len(caption_list):
    print("Something is wrong in tikz_list and caption_list!")
    exit()



# Creating the tex and pdf file
doc = lib.MyDocument()
doc.addTikzFigure(tikz_list, caption_list, col=4, row=7)
doc.generate_pdf("geo_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
if not len(crop_list) == 0:
    lib.clean_extensions(crop=crop_list)
    doc.generate_pdf("geo_example", compiler="pdflatex", clean_tex=False)
lib.clean_extensions()
print("Finished: pdf saved to geo_example.pdf")
