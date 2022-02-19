#!/usr/bin/python

""""
    @file diss_jump_example.py

    @brief Compute biharmonic2_example using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os
import sys

gismo_path = os.path.join(os.path.dirname(__file__), "../../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

import diss_library as lib

# Move the following lines to diss_library.py
from python2latex.python2latex import MyDocument
from python2latex.python2tikz import MyTikz
# End move

from enum import Enum


class Method(Enum):
    ApproxC1 = 0
    Nitsche = 1
    DPatch = 2


''' ####### USER INPUT ####### '''
geo_list = ["g1000", "g1100", "g1510", "g1400"]  # Without .xml extension
path_geo = "planar/geometries/"

loop = 4

deg_approxC1 = 3  # FIX

deg_list = [
    [1, 2, 3, 4],   # p_tilde for gluing data
    [3, 4, 5]       # Degree for Nitsche
]
method_list = [
    Method.ApproxC1,
    Method.Nitsche
]
compute_list = [
    False,
    False
]
shift = [
    1,
    0
]

xml_col = "Xml_jump_results.xml"

''' ##### USER INPUT END ##### '''

path_dir = "jump/"
path_tikz = "tikz_files/" + path_dir
path_fig = "tikz_figures/" + path_dir


max_id = 0
file_col = gs.io.gsXmlCollection(xml_col)
for idx, compute in enumerate(compute_list):
    if compute:
        for geo in geo_list:
            for deg in deg_list[idx]:
                if method_list[idx] == Method.ApproxC1:
                    ptilde = deg
                    rtilde = deg-1
                    deg = deg_approxC1
                    reg = deg - 1
                    path_bvp = "results/" + geo + "/bvp/" + "approxC1" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                                                 "-P" + str(ptilde) + "-R" + str(rtilde) + "-l" + str(loop)
                else:
                    reg = deg - 1
                    path_bvp = "results/" + geo + "/bvp/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                               "-l" + str(loop)
                    ptilde = -1
                    rtilde = -1

                lib.create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=reg, path_bvp=path_bvp,
                                          method=method_list[idx], ptilde=ptilde, rtilde=rtilde)
                path_output = path_bvp.replace("bvp/", "results/")
                path_output = path_output + "-result.xml"
                lib.biharmonic_example(path_bvp + ".xml", path_output)
                file_col.addFile(path_output)
                max_id = max_id + 1
        print("Finished!")
    else:
        for geo in geo_list:
            for deg in deg_list[idx]:
                if method_list[idx] == Method.ApproxC1:
                    ptilde = int(deg)
                    rtilde = int(ptilde)-1
                    deg = deg_approxC1
                    reg = deg - 1
                    path_bvp = "results/" + geo + "/bvp/" + "approxC1" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                               "-P" + str(ptilde) + "-R" + str(rtilde) + "-l" + str(loop)
                else:
                    reg = deg - 1
                    path_bvp = "results/" + geo + "/bvp/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                               "-l" + str(loop)
                path_output = path_bvp.replace("bvp/", "results/")
                path_output = path_output + "-result.xml"
                file_col.addFile(path_output)
                max_id = max_id + 1

file_col.save()

file1 = gs.io.gsXmlCollection("")
file1.load(xml_col)
list_dict = []
for id in range(max_id):
    my_dict = {"Matrix": file1.getMatrix(id), "Deg": None, "Reg": None,
               "Geo": None, "Name": None}

    path_name = file1.getString(id)
    name = path_name[path_name.find('/results/') + 9:]
    my_dict["Name"] = name[:name.find('.xml')]

    my_dict["Geo"] = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    my_dict["Deg"] = path_name[path_name.find('-p') + 2:path_name.find('-p') + 3]
    my_dict["Reg"] = path_name[path_name.find('-r') + 2:path_name.find('-r') + 3]
    list_dict.append(my_dict)


geo_mat_coll = []
name_mat_coll = []
for method in method_list:
    m_str = ""
    geo_mat_list = []
    name_mat_list = []
    if method == Method.ApproxC1:
        m_str = "approxC1"
        for geo in geo_list:
            mat_list = []
            for dict in list_dict:
                if m_str in dict["Name"]:
                    if dict["Geo"] == geo and int(dict["Deg"]) == deg_approxC1:
                        mat_list.append(dict["Matrix"])
            geo_mat_list.append(mat_list)
            name_mat_list.append(geo + "-" + m_str + "-jump-p" + str(deg) + "-r2-l" + str(loop))
    elif method == Method.Nitsche:
        m_str = "nitsche"
        for geo in geo_list:
            mat_list = []
            for deg in deg_list[idx]:
                for dict in list_dict:
                    if m_str in dict["Name"]:
                        if dict["Geo"] == geo and int(dict["Deg"]) == deg:
                            mat_list.append(dict["Matrix"])
            geo_mat_list.append(mat_list)
            name_mat_list.append(geo + "-" + m_str + "-jump-p" + str(deg) + "-r2-l" + str(loop))
    elif method == Method.DPatch:
        m_str = "dPatch"
        print("TODO")
    geo_mat_coll.append(geo_mat_list)
    name_mat_coll.append(name_mat_list)

# def unique_values_from_list(dict_list, key):
#     all_values = set()
#     for dict in dict_list:
#         all_values.add(str(dict[key]))
#     return all_values
#
# geo_all = unique_values_from_list(list_dict, "Geo")
# deg_list = unique_values_from_list(list_dict, "Deg")

def create_tikz_figures(geo_mat_list, name_mat_list, deg_list, opt_plot, shift=0):


    list_tikz = []
    for idx, mat_list in enumerate(geo_mat_list):
        x_col = 0

        M = []
        x = []
        for mat in mat_list:
            x.append(mat[:, x_col])  # Mesh size
            M.append(mat[:, 5])  # L2 Error + H1 Error + H2 Error

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

        fig = MyTikz()
        opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '6cm', 'mark options': '{solid}',
                    'xlabel': '{Mesh size $h$}', 'ylabel' : '{Error}',
                    'ylabel style' : '{yshift=-0.4cm}', 'xlabel style':'{yshift=0.2cm}'}
        fig.setOptions(opt_axis)
        color_list = ["red", "green", "blue", "yellow"]
        fig.setColor(color_list)
        fig.setDegree(deg_list)
        fig.setRateShift(shift)
        fig.setPlotOptions(opt_plot)
        fig.create_error_plot(x, M, True, True)  # True since M is an array
        fig.generate_tikz(path_tikz + name_mat_list[idx])

        list_tikz.append(path_dir + name_mat_list[idx])
        if idx == 0:
            legend_list = fig.getLegendList()
    return list_tikz, legend_list

opt_plot = [{'mark': 'diamond*', 'color': 'green', 'line width': '1pt'},
             {'mark': 'square*', 'color': 'blue', 'line width': '1pt'},
             {'mark': 'triangle*', 'color': 'red', 'line width': '1pt'},
             {'mark': 'pentagon*', 'color': 'yellow', 'line width': '1pt'},
             {'mark': 'halfcircle*', 'color': 'brown', 'line width': '1pt'}]

tikz_coll = []
legend_coll = []
for idx in range(len(method_list)):
    list_tikz2, legend_list2 = create_tikz_figures(geo_mat_coll[idx], name_mat_coll[idx], deg_list[idx], [opt_plot[idx]], shift[idx])
    tikz_coll.append(list_tikz2)
    legend_coll.append(legend_list2)

legend_list = []
for idx, method in enumerate(method_list):
    if method == Method.ApproxC1:
        legend_image = []
        legend_entry = []
        legend_image.append(["empty legend"])
        legend_entry.append([r'\hspace{-0.5cm}Jump-error'])
        for image in legend_coll[idx]:
            legend_entry.append([""])
            legend_image.append(image)
        legend_entry.append([""])
        legend_image.append(["empty legend"])
        for deg in deg_list[idx]:
            legend_image.append(["empty legend"])
            legend_entry.append([r'\hspace{-1.7cm}$\widetilde{p}=' + str(deg) + '$, $\widetilde{r}=' + str(deg - 1) + '$'])
        fig = MyTikz()
        fig.create_legend(legend_image, legend_entry)
        fig.generate_tikz(path_tikz + "legend_jump_approxC1")
        legend_list.append([path_dir + "legend_jump_approxC1"])

    if method == Method.Nitsche:
        legend_image2 = []
        legend_entry2 = []
        legend_image2.append(["empty legend"])
        legend_entry2.append([r'\hspace{-0.5cm}Jump-error'])
        for image in legend_coll[idx]:
            legend_entry2.append([""])
            legend_image2.append(image)
        legend_entry2.append([""])
        legend_image2.append(["empty legend"])
        for deg in deg_list[idx]:
            legend_image2.append(["empty legend"])
            legend_entry2.append([r'\hspace{-1.7cm}$p=' + str(deg) + '$, $r=' + str(deg - 1) + '$'])
        fig = MyTikz()
        fig.create_legend(legend_image2, legend_entry2)
        fig.generate_tikz(path_tikz + "legend_jump_nitsche")
        legend_list.append([path_dir + "legend_jump_nitsche"])

caption_coll = []
for method in method_list:
    caption_list = []
    if method == Method.ApproxC1:
        caption_list.append('Ex. I: Approx. $C^1$')
        caption_list.append('Ex. II: Approx. $C^1$')
        caption_list.append('Ex. III: Approx. $C^1$')
        caption_list.append('Ex. IV: Approx. $C^1$')
    if method == Method.Nitsche:
        caption_list.append('Ex. I: Nitsche')
        caption_list.append('Ex. II: Nitsche')
        caption_list.append('Ex. III: Nitsche')
        caption_list.append('Ex. IV: Nitsche')
    caption_coll.append(caption_list)


doc = MyDocument()
for idx in range(len(method_list)):
    doc.addTikzFigure(tikz_coll[idx], caption_coll[idx], row=4, legend=legend_list[idx])
doc.generate_pdf("jump_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
for idx in range(len(method_list)):
    doc.clean_extensions(crop=legend_list[idx])
doc.generate_pdf("jump_example", compiler="pdflatex", clean_tex=False)
doc.clean_extensions()
