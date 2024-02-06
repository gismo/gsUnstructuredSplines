#!/usr/bin/python

""""
    @file diss_nitsche_example.py
    @brief Compute biharmonic2_example using pygismo
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): P. Weinmüller
"""

import os
import sys
import subprocess

import numpy as np

path_module = "data/"
print("Module path:", path_module, "(change if needed).")
os.chdir(os.path.join(os.path.dirname(__file__), "../../../" + path_module))

gismo_path = "../build/lib"
print("G+Smo path:", os.getcwd() + "/" + gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

import library.library as lib

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    DPatch = 1
    AlmostC1 = 2
    Nitsche = 3


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

geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
#geo_list = ["g1001", "g1021", "g1030", "g1031"]  # Without .xml extension
path_geo = domain + "/geometries/"

numRefinement = 5
degree = 3

second = False

h = -1
N = 1

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


file_col = []
for geo in geo_list:
    path_results = "results/" + path_dir + geo + "/"
    for pen in penalty:
        argument_list = "nitsche-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                        + "-m" + str(Method.Nitsche.value) + ("--second" if second else "") + ("--residual" if residual else "") + "-y" + str(pen)
        file_col.append(path_results + argument_list)

    # Approx C1
    argument_list = "approxC1-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                    + "-m" + str(Method.ApproxC1.value) + ("--second" if second else "") + ("--residual" if residual else "")
    file_col.append(path_results + argument_list)

    # Nitsche with EW
    argument_list = "nitsche-g" + geo + "-p" + str(degree) + "-s" + str(degree - 1) + "-r" + str(numRefinement) \
                    + "-m" + str(Method.Nitsche.value) + ("--second" if second else "") + ("--residual" if residual else "") + "-y" + str(-1)
    file_col.append(path_results + argument_list)


list_dict = []
for idx in range(len(file_col)):
    my_dict = {"Matrix": gs.io.gsFileData.getMatrix(file_col[idx]+".xml"), "Deg": None, "Reg": None,
               "Geo": None, "Name": None}

    path_name = file_col[idx]
    name = path_name[path_name.find("results/" + path_dir) + len("results/" + path_dir):]
    my_dict["Name"] = name

    my_dict["Geo"] = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    my_dict["Deg"] = path_name[path_name.find('-p') + 2:path_name.find('-p') + 3]
    my_dict["Reg"] = path_name[path_name.find('-r') + 2:path_name.find('-r') + 3]
    list_dict.append(my_dict)

m_str = ""
geo_mat_list = []
geo_mat_list2 = []
name_mat_list = []
name_mat_list2 = []
m_str = "nitsche"
# Approx C1
geo_mat_list_approxC1 = []
m_str2 = "approxC1"
for geo in geo_list:
    mat_list = []
    for dict in list_dict:
        if m_str in dict["Name"]:
            if dict["Geo"] == geo and int(dict["Deg"]) == degree and "y-1" not in dict["Name"]:
                mat_list.append(dict["Matrix"])
            if dict["Geo"] == geo and int(dict["Deg"]) == degree and "y-1" in dict["Name"]:
                geo_mat_list2.append([dict["Matrix"]])

    geo_mat_list.append(mat_list)
    name_mat_list.append(geo + "-nitsche-penalty-p" + str(degree) + "-r" + str(degree-1)
                         + "-l" + str(numRefinement) + "-h" + str(h))
    name_mat_list2.append(geo + "-nitsche-penalty-jump-p" + str(degree) + "-r" + str(degree-1)
                         + "-l" + str(numRefinement) + "-h" + str(h))

    mat_list = []
    for dict in list_dict:
        if m_str2 in dict["Name"]:
            if dict["Geo"] == geo and int(dict["Deg"]) == degree:
                mat_list.append(dict["Matrix"])
    geo_mat_list_approxC1.append(mat_list)

list_tikz = []
for idx, mat_list in enumerate(geo_mat_list):  # idx = geo

    x_col = 7  # Maybe change if more then one interface

    M_list = []
    x_list = []
    M = np.zeros((len(mat_list), 3))
    x = np.zeros((1, len(mat_list)))
    for i, mat in enumerate(mat_list):
        x[0, i] = mat[h, x_col]  # Mesh size
        M[i, :] = mat[h, 2:5]  # L2 Error + H1 Error + H2 Error
    M_list.append(M)
    x_list.append(x[0])

    M = np.zeros((len(mat_list), 3))
    for i, mat in enumerate(mat_list):
        M[i, :] = geo_mat_list_approxC1[idx][0][h,2:5]
    M_list.append(M)
    x_list.append(x[0])

    a = geo_mat_list2[idx][0][h,4] * 5
    b = geo_mat_list2[idx][0][h,2] * 0.2
    yy = np.linspace(a, b, num=len(mat_list))

    M = np.zeros((len(mat_list), 1))
    x = np.zeros((1, len(mat_list)))
    for i, mat in enumerate(mat_list):
        M[i, :] = yy[i]
        x[0,i] = geo_mat_list2[idx][0][h, x_col]
    M_list.append(M)
    x_list.append(x[0])

    opt_plot = [{'color': 'green', 'line width': '1pt'},
                {'color': 'blue', 'line width': '1pt'},
                {'color': 'red', 'line width': '1pt'}]

    fig = lib.MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '4.5cm', 'mark options': '{solid}',
                'xlabel': '{Stability parameter $\eta$}', 'ylabel': '{Error}',
                'ylabel style': '{yshift=-0.4cm}', 'xlabel style': '{yshift=0.2cm}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.setPlotOptions(opt_plot)
    fig.create_error_plot(x_list, M_list, False, False)
    fig.generate_tikz(path_tikz + name_mat_list[idx])

    list_tikz.append(path_dir + name_mat_list[idx])

for idx, mat_list in enumerate(geo_mat_list):

    x_col = 7  # Maybe change if more then one interface

    M_list = []
    x_list = []
    M = np.zeros((len(mat_list), 1))
    x = np.zeros((1, len(mat_list)))
    for i, mat in enumerate(mat_list):
        x[0, i] = mat[h, x_col]  # Mesh size
        M[i, :] = mat[h, 5]  # L2 Error + H1 Error + H2 Error
    M_list.append(M)
    x_list.append(x[0])

    M = np.zeros((len(mat_list), 3))
    for i, mat in enumerate(mat_list):
        M[i, :] = geo_mat_list_approxC1[0][0][h,5]

    # For Approx C1 Jump
    #M_list.append(M)
    #x_list.append(x)

    a = geo_mat_list2[idx][0][h,5] * 5
    b = geo_mat_list2[idx][0][h,5] * 0.2
    yy = np.linspace(a, b, num=len(mat_list))

    M = np.zeros((len(mat_list), 1))
    x = np.zeros((1, len(mat_list)))
    for i, mat in enumerate(mat_list):
        M[i, :] = yy[i]
        x[0,i] = geo_mat_list2[idx][0][h, x_col]
    M_list.append(M)
    x_list.append(x[0])

    opt_plot = [{'color': 'green', 'line width': '1pt'},
                {'color': 'blue', 'line width': '1pt'},
                {'color': 'red', 'line width': '1pt'}]

    fig = lib.MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '4.5cm', 'mark options': '{solid}',
                'xlabel': '{Stability parameter $\eta$}', 'ylabel': '{Error}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.setPlotOptions(opt_plot)
    fig.create_error_plot(x_list, M_list, False, False)
    fig.generate_tikz(path_tikz + name_mat_list2[idx])

    list_tikz.append(path_dir + name_mat_list2[idx])

caption_list = []
caption_list.append('Ex. I: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. II: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. III: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. IV: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. I: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. II: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. III: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')
caption_list.append('Ex. IV: $p=' + str(degree) + '$, $r=' + str(degree-1) + '$')

doc = lib.MyDocument()
doc.addTikzFigure(list_tikz, caption_list, row=2, col=4)
doc.generate_pdf("nitsche_stabilization_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions()