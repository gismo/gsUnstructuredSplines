#!/usr/bin/python

""""
    @file plotting_error2_example.py

    @brief Plotting the results via pylatex

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

import library.library as lib

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
geo_list = ["g1000", "g1100", "g1021","g1121"]  # Without .xml extension
path_geo = "planar/geometries/"

NumRefinement = 2
second = False

deg_list = [3, 4, 5]  # For all methods the same


# Approx C1: gluing data set to default: \tilde{p} = p-1, \tilde{r} = p-2,
# Nitsche: penalty set to default: via Eigenvalue-problem
method_list = [
    Method.ApproxC1,
    Method.DPatch,
    Method.Nitsche
]

path_example = "../build/bin/biharmonic3_example"
""" -------------------------------------------------------------------------------------------------- """

path_dir = "error/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

if not os.path.exists(path_tikz):
    os.makedirs(path_tikz)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

file_coll = []
for method in method_list:
    for geo in geo_list:
        for deg in deg_list:
            m_str = ""
            if method == Method.ApproxC1:
                m_str = "approxC1"
            elif method == Method.Nitsche:
                m_str = "nitsche"
            elif method == Method.DPatch:
                m_str = "dPatch"
            elif method == Method.AlmostC1:
                m_str = "almostC1"
            else:
                print("METHOD NOT IMPLEMENTED!!!")

            path_results_geo = "results/error/" + geo + "/"
            argument_list = m_str + "-g" + geo + "-p" + str(deg) + "-s" + str(deg - 1) + "-r" + str(NumRefinement) \
                            + "-m" + str(method.value) + ("-second" if second else "")

            file_coll.append(path_results_geo + argument_list)

list_dict = []
for idx in range(len(file_coll)):

    #if np.array(gs.io.gsFileData.getMatrix(file_coll[idx] + ".xml")).shape[0] == 0:
    #    print("Empty matrix in ", file_coll[idx])
    #    continue

    my_dict = {"Matrix": gs.io.gsFileData.getMatrix(file_coll[idx] + ".xml"), "Deg": None, "Reg": None,
               "Geo": None, "Name": None}

    path_name = file_coll[idx]
    name = path_name[path_name.find('/results/') + 9:]
    my_dict["Name"] = name[:name.find('.xml')]

    my_dict["Geo"] = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    my_dict["Deg"] = path_name[path_name.find('-p') + 2:path_name.find('-p') + 3]
    my_dict["Reg"] = path_name[path_name.find('-s') + 2:path_name.find('-s') + 3]
    list_dict.append(my_dict)

# Sorting values
row_mat_coll = []
name_mat_coll = []
deg_mat_coll = []
for geo in geo_list:  # Cols
    geo_mat_list = []
    name_mat_list = []

    for deg in deg_list:  # Rows
        mat_list = []
        for dict in list_dict:
            if dict["Geo"] == geo and int(dict["Deg"]) == deg:
                mat_list.append(dict["Matrix"])
        geo_mat_list.append(mat_list)
        name_mat_list.append(geo + "-error-p" + str(deg) + "-s" + str(deg - 1) + "-r" + str(NumRefinement))

    row_mat_coll.append(geo_mat_list)
    name_mat_coll.append(name_mat_list)


def create_tikz_figure(geo_mat, name_mat, deg_list, opt_plot):
    list_tikz = []

    x_col = 0

    M_list = []
    x_list = []
    for mat in geo_mat:
        if mat.shape[0] != 0:
            x_list.append(mat[:, x_col])  # Mesh size
            M_list.append(mat[:, 2:5])  # L2 Error + H1 Error + H2 Error

    '''Computing the rates'''
    # l2error = M[:,0]
    # rate_l2 = np.log(l2error[:-1]/l2error[1:])/np.log(2.0)
    # print(np.around(rate_l2,2))
    #
    # h1error = M[:,1]
    # rate_h1 = np.log(h1error[:-1]/h1error[1:])/np.log(2.0)
    # print(np.around(rate_h1,2))
    #
    # h2error = M[:,2]
    # rate_h2 = np.log(h2error[:-1]/h2error[1:])/np.log(2.0)
    # print(np.around(rate_h2,2))

    fig = lib.MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '6cm', 'mark options': '{solid}',
                'xlabel': '{Mesh size $h$}', 'ylabel': '{Error}',
                'ylabel style': '{yshift=-0.4cm}', 'xlabel style': '{yshift=0.2cm}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.setDegree(deg_list)
    #fig.setRateShift(shift)
    fig.setPlotOptions(opt_plot)
    #fig.swapLinestyle(id)
    fig.create_error_plot(x_list, M_list, False, rate=True)  # True since M is an array
    fig.generate_tikz(path_tikz + name_mat)

    list_tikz.append(path_dir + name_mat)
    legend_list = fig.getLegendList()
    return list_tikz, legend_list


opt_plot = [{'mark': 'diamond*', 'color': 'green', 'line width': '1pt'},
            {'mark': 'square*', 'color': 'blue', 'line width': '1pt'},
            {'mark': 'triangle*', 'color': 'red', 'line width': '1pt'},
            {'mark': 'pentagon*', 'color': 'yellow', 'line width': '1pt'},
            {'mark': 'halfcircle*', 'color': 'brown', 'line width': '1pt'}]

tikz_coll = []
legend_coll = []

for idx, deg in enumerate(deg_list):
    for geo in range(len(geo_list)):
        rate_list = [deg + 1, deg, deg - 1]

        list_tikz2, legend_list2 = create_tikz_figure(row_mat_coll[geo][idx], name_mat_coll[geo][idx], rate_list, opt_plot)

        for tikz in list_tikz2:
            tikz_coll.append(tikz)
        if geo == 0:
            for legend in legend_list2:
                legend_coll.append(legend)


legend_list = []

legend_image = []
legend_entry = []
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{-0.8cm}$L^2$-error'])
for idx in range(len(method_list)):
    legend_entry.append([""])
    legend_image.append(legend_coll[idx*3])
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{-0.8cm}$H^1$-error'])
for idx in range(len(method_list)):
    legend_entry.append([""])
    legend_image.append(legend_coll[idx*3+1])
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{-0.8cm}$H^2$-error'])
for idx in range(len(method_list)):
    legend_image.append(legend_coll[idx*3+2])

legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{+0.5cm}Approx. $C^1$'])
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{+0.5cm}D-Patch'])
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{+0.5cm}Nitsche'])
# ADD MORE METHODS

fig = lib.MyTikz()
fig.create_legend(legend_image, legend_entry, col=3)
fig.generate_tikz(path_tikz + "legend_error")
legend_list.append(path_dir + "legend_error")


caption_coll = []
for deg in deg_list:
    for geo in geo_list:
        caption_coll.append("$p = " + str(deg) + "$, $r = " + str(deg - 1) + "$")


doc = lib.MyDocument()
doc.addTikzFigure(tikz_coll, caption_coll, row=3, col=4, legend=legend_list)
doc.generate_pdf("error_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions(crop=legend_list)
doc.generate_pdf("error_example", compiler="pdflatex", clean_tex=False)
lib.clean_extensions()