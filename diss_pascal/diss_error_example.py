#!/usr/bin/python

""""
    @file diss_error_example.py

    @brief Compute biharmonic2_example using pygismo

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

import lib.diss_library as lib

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    Nitsche = 1
    DPatch = 2

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
geo_list = ["g1021", "g1121", "g1501", "g1311"]  # Without .xml extension
path_geo = "planar/geometries/"

loop = 4

deg_list = [
    3, 4, 5
]
method_list = [
    Method.ApproxC1,
    Method.Nitsche
]
compute_list = [
    False,
    False
]

""" -------------------------------------------------------------------------------------------------- """

path_dir = "error/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

xml_col = "results/xml_files/xml_error_results.xml"


file_col = []
for idx, compute in enumerate(compute_list):
    if compute:
        for geo in geo_list:

            # Making new folder if there exist no folder.
            if not os.path.exists("results/" + geo + "/bvp/"):
                os.makedirs("results/" + geo + "/bvp/")
            if not os.path.exists("results/" + geo + "/results/"):
                os.makedirs("results/" + geo + "/results/")

            for deg in deg_list:
                m_str = ""
                if method_list[idx] == Method.ApproxC1:
                    m_str = "approxC1"
                elif method_list[idx] == Method.Nitsche:
                    m_str = "nitsche"
                elif method_list[idx] == Method.DPatch:
                    m_str = "dPatch"
                else:
                    print("METHOD NOT IMPLEMENTED!!!")

                path_bvp = "results/" + geo + "/bvp/" + m_str + "-bvp1-p" + str(deg) + "-r" + str(deg-1) + \
                           "-l" + str(loop)

                lib.create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=deg-1, path_bvp=path_bvp,
                                          method=method_list[idx])
                path_output = path_bvp.replace("bvp/", "results/")
                path_output = path_output + "-result.xml"
                lib.biharmonic_example(path_bvp + ".xml", path_output)
                file_col.append(path_output)
        print("Finished!")
    else:
        for geo in geo_list:
            for deg in deg_list:
                m_str = ""
                if method_list[idx] == Method.ApproxC1:
                    m_str = "approxC1"
                elif method_list[idx] == Method.Nitsche:
                    m_str = "nitsche"
                elif method_list[idx] == Method.DPatch:
                    m_str = "dPatch"
                else:
                    print("METHOD NOT IMPLEMENTED!!!")

                path_bvp = "results/" + geo + "/bvp/" + m_str + "-bvp1-p" + str(deg) + "-r" + str(deg-1) + \
                           "-l" + str(loop)
                path_output = path_bvp.replace("bvp/", "results/")
                path_output = path_output + "-result.xml"
                file_col.append(path_output)


list_dict = []
for id in range(len(file_col)):
    my_dict = {"Matrix": gs.io.gsFileData.getMatrix(file_col[id]), "Deg": None, "Reg": None,
               "Geo": None, "Name": None}

    path_name = file_col[id]
    name = path_name[path_name.find('/results/') + 9:]
    my_dict["Name"] = name[:name.find('.xml')]

    my_dict["Geo"] = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    my_dict["Deg"] = path_name[path_name.find('-p') + 2:path_name.find('-p') + 3]
    my_dict["Reg"] = path_name[path_name.find('-r') + 2:path_name.find('-r') + 3]
    list_dict.append(my_dict)

row_mat_coll = []
name_mat_coll = []
for deg in deg_list:  # Rows
    geo_mat_list = []
    name_mat_list = []

    for geo in geo_list:  # Cols
        mat_list = []
        for dict in list_dict:
            if dict["Geo"] == geo and int(dict["Deg"]) == deg:
                mat_list.append(dict["Matrix"])
        geo_mat_list.append(mat_list)
        name_mat_list.append(geo + "-error-p" + str(deg) + "-r" + str(deg - 1) + "-l" + str(loop))

    row_mat_coll.append(geo_mat_list)
    name_mat_coll.append(name_mat_list)


def create_tikz_figures(geo_mat_list, name_mat_list, deg_list, opt_plot, shift=0):
    list_tikz = []
    for idx, mat_list in enumerate(geo_mat_list):
        x_col = 0

        M_list = []
        x_list = []
        for mat in mat_list:
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
        fig.setRateShift(shift)
        fig.setPlotOptions(opt_plot)
        fig.create_error_plot(x_list, M_list, False, rate=True)  # True since M is an array
        fig.generate_tikz(path_tikz + name_mat_list[idx])

        list_tikz.append(path_dir + name_mat_list[idx])
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
    rate_list = [deg + 1, deg, deg - 1]
    list_tikz2, legend_list2 = create_tikz_figures(row_mat_coll[idx], name_mat_coll[idx], rate_list, opt_plot)
    for tikz in list_tikz2:
        tikz_coll.append(tikz)
    if idx == 0:
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

if range(len(method_list)) != 2:
    print("Need here more entries!")
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{+0.5cm}Approx. $C^1$'])
legend_image.append(["empty legend"])
legend_entry.append([r'\hspace{+0.5cm}Nitsche'])
# TODO ADD MORE METHODS

fig = lib.MyTikz()
fig.create_legend(legend_image, legend_entry, col=3)
fig.generate_tikz(path_tikz + "legend_error")
legend_list.append(path_dir + "legend_error")


caption_coll = []
for deg in deg_list:
    for geo in geo_list:
        caption_coll.append("$p=" + str(deg) + "$, $r=" + str(deg - 1) + "$")

doc = lib.MyDocument()
doc.addTikzFigure(tikz_coll, caption_coll, row=4, legend=legend_list)
doc.generate_pdf("error_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
lib.clean_extensions(crop=legend_list)
doc.generate_pdf("error_example", compiler="pdflatex", clean_tex=False)
lib.clean_extensions()
