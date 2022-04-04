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
import subprocess


path_module = "data/"
print("Module path:", path_module, "(change if needed).")
os.chdir(os.path.join(os.path.dirname(__file__), "../../../"+path_module))

gismo_path = "../build/lib"
print("G+Smo path:", os.getcwd() + "/" + gismo_path, "(change if needed).")
sys.path.append(gismo_path)

#import pygismo as gs

from enum import Enum


class Method(Enum):
    ApproxC1 =  0
    DPatch =    1
    AlmostC1 =  2
    Nitsche =   3

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
geo_list = ["g1000", "g1020", "g1702", "g1021", "g1704"]  # Without .xml extension
path_geo = "planar/geometries/"

loop = 1
second = False

deg_list = [
    [3, 4, 5],
    [3, 4, 5],
    [2]
]
method_list = [
    Method.ApproxC1,
    Method.DPatch,
    Method.AlmostC1
]
compute_list = [
    True,
    True,
    True
]

path_example = "../build/bin/biharmonic_multiPatch_example"
""" -------------------------------------------------------------------------------------------------- """

path_dir = "error/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir

file_coll = []
for idx, compute in enumerate(compute_list):
    if compute:
        for geo in geo_list:

            # Making new folder if there exist no folder.
            path_results_geo = "results/" + geo + "/results/"
            if not os.path.exists(path_results_geo):
                os.makedirs(path_results_geo)

            for deg in deg_list[idx]:
                m_str = ""
                if method_list[idx] == Method.ApproxC1:
                    m_str = "approxC1"
                elif method_list[idx] == Method.Nitsche:
                    m_str = "nitsche"
                elif method_list[idx] == Method.DPatch:
                    m_str = "dPatch"
                elif method_list[idx] == Method.AlmostC1:
                    m_str = "almostC1"
                else:
                    print("METHOD NOT IMPLEMENTED!!!")

                argument_list = m_str + "-g" + geo + "-p" + str(deg) + "-s" + str(deg-1) + "-r" + str(loop) \
                                + "-m" + str(method_list[idx].value) + ("-second" if second else "")

                # [!Run biharmonic2_example]
                proc = subprocess.Popen([path_example, "-g", geo, "-p", str(deg), "-s", str(deg-1), "-r", str(loop),
                                         "-m", str(method_list[idx].value), ("--second" if second else ""), "", "-o", path_results_geo+argument_list])
                proc.wait()
                # [!Run biharmonic2_example]

                file_coll.append(path_results_geo+argument_list)
        print("Finished!")
    else:
        for geo in geo_list:
            for deg in deg_list[idx]:
                m_str = ""
                if method_list[idx] == Method.ApproxC1:
                    m_str = "approxC1"
                elif method_list[idx] == Method.Nitsche:
                    m_str = "nitsche"
                elif method_list[idx] == Method.DPatch:
                    m_str = "dPatch"
                elif method_list[idx] == Method.AlmostC1:
                    m_str = "almostC1"
                else:
                    print("METHOD NOT IMPLEMENTED!!!")

                path_results_geo = "results/" + geo + "/results/"
                argument_list = m_str + "-g" + geo + "-p" + deg + "-s" + str(deg-1) + "-r" + loop \
                                + "-m" + str(method_list[idx].value) + ("-second" if second else "")

                file_coll.append(path_results_geo+argument_list)

print(file_coll)
#
# list_dict = []
# for id in range(len(file_col)):
#     my_dict = {"Matrix": gs.io.gsFileData.getMatrix(file_col[id]), "Deg": None, "Reg": None,
#                "Geo": None, "Name": None}
#
#     path_name = file_col[id]
#     name = path_name[path_name.find('/results/') + 9:]
#     my_dict["Name"] = name[:name.find('.xml')]
#
#     my_dict["Geo"] = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
#     my_dict["Deg"] = path_name[path_name.find('-p') + 2:path_name.find('-p') + 3]
#     my_dict["Reg"] = path_name[path_name.find('-r') + 2:path_name.find('-r') + 3]
#     list_dict.append(my_dict)
#
# row_mat_coll = []
# name_mat_coll = []
# for deg in deg_list:  # Rows
#     geo_mat_list = []
#     name_mat_list = []
#
#     for geo in geo_list:  # Cols
#         mat_list = []
#         for dict in list_dict:
#             if dict["Geo"] == geo and int(dict["Deg"]) == deg:
#                 mat_list.append(dict["Matrix"])
#         geo_mat_list.append(mat_list)
#         name_mat_list.append(geo + "-error-p" + str(deg) + "-r" + str(deg - 1) + "-l" + str(loop))
#
#     row_mat_coll.append(geo_mat_list)
#     name_mat_coll.append(name_mat_list)
#
#
# def create_tikz_figures(geo_mat_list, name_mat_list, deg_list, opt_plot, shift=0):
#     list_tikz = []
#     for idx, mat_list in enumerate(geo_mat_list):
#         x_col = 0
#
#         M_list = []
#         x_list = []
#         for mat in mat_list:
#             x_list.append(mat[:, x_col])  # Mesh size
#             M_list.append(mat[:, 2:5])  # L2 Error + H1 Error + H2 Error
#
#         '''Computing the rates'''
#         # l2error = M[:,0]
#         # rate_l2 = np.log(l2error[:-1]/l2error[1:])/np.log(2.0)
#         # print(np.around(rate_l2,2))
#         #
#         # h1error = M[:,1]
#         # rate_h1 = np.log(h1error[:-1]/h1error[1:])/np.log(2.0)
#         # print(np.around(rate_h1,2))
#         #
#         # h2error = M[:,2]
#         # rate_h2 = np.log(h2error[:-1]/h2error[1:])/np.log(2.0)
#         # print(np.around(rate_h2,2))
#
#         fig = lib.MyTikz()
#         opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '6cm', 'mark options': '{solid}',
#                     'xlabel': '{Mesh size $h$}', 'ylabel': '{Error}',
#                     'ylabel style': '{yshift=-0.4cm}', 'xlabel style': '{yshift=0.2cm}'}
#         fig.setOptions(opt_axis)
#         color_list = ["red", "green", "blue", "yellow"]
#         fig.setColor(color_list)
#         fig.setDegree(deg_list)
#         fig.setRateShift(shift)
#         fig.setPlotOptions(opt_plot)
#         fig.create_error_plot(x_list, M_list, False, rate=True)  # True since M is an array
#         fig.generate_tikz(path_tikz + name_mat_list[idx])
#
#         list_tikz.append(path_dir + name_mat_list[idx])
#         legend_list = fig.getLegendList()
#     return list_tikz, legend_list
#
#
# opt_plot = [{'mark': 'diamond*', 'color': 'green', 'line width': '1pt'},
#             {'mark': 'square*', 'color': 'blue', 'line width': '1pt'},
#             {'mark': 'triangle*', 'color': 'red', 'line width': '1pt'},
#             {'mark': 'pentagon*', 'color': 'yellow', 'line width': '1pt'},
#             {'mark': 'halfcircle*', 'color': 'brown', 'line width': '1pt'}]
#
# tikz_coll = []
# legend_coll = []
# for idx, deg in enumerate(deg_list):
#     rate_list = [deg + 1, deg, deg - 1]
#     list_tikz2, legend_list2 = create_tikz_figures(row_mat_coll[idx], name_mat_coll[idx], rate_list, opt_plot)
#     for tikz in list_tikz2:
#         tikz_coll.append(tikz)
#     if idx == 0:
#         for legend in legend_list2:
#             legend_coll.append(legend)
#
#
# legend_list = []
#
# legend_image = []
# legend_entry = []
# legend_image.append(["empty legend"])
# legend_entry.append([r'\hspace{-0.8cm}$L^2$-error'])
# for idx in range(len(method_list)):
#     legend_entry.append([""])
#     legend_image.append(legend_coll[idx*3])
# legend_image.append(["empty legend"])
# legend_entry.append([r'\hspace{-0.8cm}$H^1$-error'])
# for idx in range(len(method_list)):
#     legend_entry.append([""])
#     legend_image.append(legend_coll[idx*3+1])
# legend_image.append(["empty legend"])
# legend_entry.append([r'\hspace{-0.8cm}$H^2$-error'])
# for idx in range(len(method_list)):
#     legend_image.append(legend_coll[idx*3+2])
#
# if range(len(method_list)) != 2:
#     print("Need here more entries!")
# legend_image.append(["empty legend"])
# legend_entry.append([r'\hspace{+0.5cm}Approx. $C^1$'])
# legend_image.append(["empty legend"])
# legend_entry.append([r'\hspace{+0.5cm}Nitsche'])
# # TODO ADD MORE METHODS
#
# fig = lib.MyTikz()
# fig.create_legend(legend_image, legend_entry, col=3)
# fig.generate_tikz(path_tikz + "legend_error")
# legend_list.append(path_dir + "legend_error")
#
#
# caption_coll = []
# for deg in deg_list:
#     for geo in geo_list:
#         caption_coll.append("$p=" + str(deg) + "$, $r=" + str(deg - 1) + "$")
#
# doc = lib.MyDocument()
# doc.addTikzFigure(tikz_coll, caption_coll, row=4, legend=legend_list)
# doc.generate_pdf("error_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
# lib.clean_extensions(crop=legend_list)
# doc.generate_pdf("error_example", compiler="pdflatex", clean_tex=False)
# lib.clean_extensions()