#!/usr/bin/python

""""
    @file diss_gluingData_example.py

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

import numpy as np

from pylatex import Document, NoEscape, TikZ, TikZOptions, Axis, Plot
import lib.diss_library as lib

"""
    To run the diss_gluingData_example.py, we have the following options:
    
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
geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
#geo_list = ["g1021", "g1121", "g1500", "g1311"]  # Without .xml extension
path_geo = "planar/geometries/"

caption_list = ["Ex. I: gluing data", "Ex. II: gluing data", "Ex. III: gluing data", "Ex. IV: gluing data"]

#geo_list = ["g1000", "g1100", "g1510", "g1400","g1021", "g1121", "g1501", "g1311"]
#caption_list = geo_list + geo_list

#geo_list = [f[:f.find(".xml")] for f in os.listdir("../filedata/" + path_geo) ]
#caption_list = geo_list + geo_list

numData = 50
""" -------------------------------------------------------------------------------------------------- """
# Here the tikz files and the corresponding pdf are stored
# Change if necessary
path_dir = "gluingData/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

list_tikz = []
list_legend = []
for geo in geo_list:
    # Making new folder if there exist no folder.
    if not os.path.exists("tikz_figures/" + path_dir):
        os.makedirs("tikz_figures/" + path_dir)
    if not os.path.exists("tikz_files/" + path_dir):
        os.makedirs("tikz_files/" + path_dir)

    mp = gs.core.gsMultiPatch()
    file = gs.io.gsFileData(path_geo + geo + ".xml")
    file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch
    mp.computeTopology(1e-4, False)

    if len(mp.interfaces()) != 1:
        print("More than one interfaces is not possible.")
        exit()

    for idx, int in enumerate(mp.interfaces()):

        dir2 = 0
        points2 = np.zeros((2, numData))
        uv = np.linspace(0, 1, numData)
        if int.second().side().index() == 1:  # West
            points2[1, :] = uv
            dir2 = 1
        elif int.second().side().index() == 2:  # East
            points2 = np.ones((2, numData))
            points2[1, :] = uv
            dir2 = 1
        elif int.second().side().index() == 3:  # South
            points2[0, :] = uv
        elif int.second().side().index() == 4:  # North
            points2 = np.ones((2, numData))
            points2[0, :] = uv

        dir1 = 0
        points1 = np.zeros((2, numData))
        if int.first().side().index() == 1:  # West
            points1[1, :] = uv
            dir1 = 1
        elif int.first().side().index() == 2:  # East
            points1 = np.ones((2, numData))
            points1[1, :] = uv
            dir1 = 1
        elif int.first().side().index() == 3:  # South
            points1[0, :] = uv
        elif int.first().side().index() == 4:  # North
            points1 = np.ones((2, numData))
            points1[0, :] = uv

        if points1[0, -1] == points2[0, 0] and points1[1, -1] == points2[1, 0]:
            points1 = np.flip(points1, axis=1)

        alpha1_vec = []
        alpha2_vec = []
        beta1_vec = []
        beta2_vec = []
        for i in range(numData):
            alpha1 = np.linalg.det(mp.patch(int.first().patch()).jacobian(points1[:,i]))

            alpha2 = np.linalg.det(mp.patch(int.second().patch()).jacobian(points2[:,i]))

            vec1 = mp.patch(int.first().patch()).jacobian(points1[:,i])[:,1-dir1]
            vec2 = mp.patch(int.second().patch()).jacobian(points2[:,i])[:,1-dir2]

            tang = mp.patch(int.second().patch()).jacobian(points2[:,i])[:,dir2]
            D1 = 1.0 / np.linalg.norm(tang)
            beta1 = - D1 * D1 * np.dot(vec1,tang)
            beta2 = - D1 * D1 * np.dot(vec2,tang)

            mat_beta = np.zeros((2,2))
            mat_beta[:,0] = vec1
            mat_beta[:,1] = vec2

            beta = np.linalg.det(mat_beta)
            if not np.linalg.norm(beta - beta1*alpha2 + beta2*alpha1) < 1e-10 \
                    and not np.linalg.norm(beta + beta1*alpha2 - beta2*alpha1) < 1e-10:
                print("Sth went wrong in the gluing data")
            else:
                alpha1_vec.append(alpha1)
                alpha2_vec.append(alpha2)
                beta1_vec.append(beta1)
                beta2_vec.append(beta2)

        M = [alpha1_vec, alpha2_vec, beta1_vec, beta2_vec]

        x = []
        for i in range(4):
            x.append(uv)

        opt_plot = [{'color': 'green', 'line width': '1pt'},
                    {'color': 'blue', 'line width': '1pt'},
                    {'color': 'red', 'line width': '1pt'},
                    {'color': 'yellow', 'line width': '1pt'},
                    {'color': 'brown', 'line width': '1pt'}]

        fig = lib.MyTikz()
        opt_axis = {'height': '4cm', 'mark options': '{solid}',
                    'xlabel': '{$u$}', #'ylabel': '{$v$}',
                    'ylabel style': '{yshift=-0.4cm}', 'xlabel style': '{yshift=0.2cm}'}
        fig.setOptions(opt_axis)
        color_list = ["red", "green", "blue", "yellow"]
        fig.setColor(color_list)
        fig.setPlotOptions(opt_plot)
        fig.create_error_plot(x, M, array=True, rate=False)  # True since M is an array
        fig.generate_tikz(path_tikz + geo + "gluing_data" + str(idx))
        list_tikz.append(path_dir + geo + "gluing_data" + str(idx))
        list_legend = fig.getLegendList()

legend_list = []
legend_image = []
legend_entry = []
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{-0.5cm}Gluing data'])
for image in list_legend:
    legend_entry.append([""])
    legend_image.append(image)
legend_entry.append([""])
legend_image.append(["empty legend"])

legend_image.append(["empty legend"])
legend_entry.append([r"\hspace{-2.3cm}$\alpha^{(k)}$"])
legend_image.append(["empty legend"])
legend_entry.append([r"\hspace{-2.3cm}$\alpha^{(l)}$"])
legend_image.append(["empty legend"])
legend_entry.append([r"\hspace{-2.3cm}$\beta^{(k)}$"])
legend_image.append(["empty legend"])
legend_entry.append([r"\hspace{-2.3cm}$\beta^{(l)}$"])
fig = lib.MyTikz()
fig.create_legend(legend_image, legend_entry)
fig.generate_tikz(path_tikz + "legend_gluing_data")
legend_list.append(path_dir + "legend_gluing_data")

doc = lib.MyDocument()
doc.addTikzFigure(list_tikz, caption_list, row=4, legend=legend_list)
doc.generate_pdf("gluingData_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions(crop=legend_list)
doc.generate_pdf("gluingData_example", compiler="pdflatex", clean_tex=False)
lib.clean_extensions()