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
import glob, errno

import numpy as np

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
geo_list = ["g1000", "g1100", "g1400", "g1510"]  # Without .xml extension
path_geo = "planar/geometries/"

loop = 4

deg = 4
method = Method.Nitsche

compute = False

xml_col = "Xml_nitsche_results.xml"

h = -1
N = 500

penalty = np.linspace(-10, 20, num=N)
penalty = np.around(np.power(2, penalty), decimals=5)
''' ##### USER INPUT END ##### '''

max_id = 0
file_col = gs.io.gsXmlCollection(xml_col)
if compute:
    for geo in geo_list:
        for pen in penalty:
            reg = deg - 1
            path_bvp = "results/" + geo + "/bvp/nitsche/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                       "-l" + str(loop) + "-y" + str(pen)
            ptilde = -1
            rtilde = -1

            lib.create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=reg, path_bvp=path_bvp,
                                      method=method, ptilde=ptilde, rtilde=rtilde, penalty_init=pen)
            path_output = path_bvp.replace("bvp/", "results/")
            path_output = path_output + "-result.xml"
            lib.biharmonic_example(path_bvp + ".xml", path_output)
            file_col.addFile(path_output)
            max_id = max_id + 1
        # Approx C1
        ptilde = deg
        rtilde = deg - 1
        reg = deg - 1
        path_bvp = "results/" + geo + "/bvp/nitsche/" + "approxC1" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                   "-P" + str(ptilde) + "-R" + str(rtilde) + "-l" + str(loop)
        lib.create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=reg, path_bvp=path_bvp,
                                  method=Method.ApproxC1, ptilde=ptilde, rtilde=rtilde)
        path_output = path_bvp.replace("bvp/", "results/")
        path_output = path_output + "-result.xml"
        lib.biharmonic_example(path_bvp + ".xml", path_output)
        file_col.addFile(path_output)
        max_id = max_id + 1

        # Nitsche with EW
        pen = -1
        reg = deg - 1
        path_bvp = "results/" + geo + "/bvp/nitsche/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                   "-l" + str(loop) + "-y" + str(pen)
        ptilde = -1
        rtilde = -1

        lib.create_biharmonic_bvp(path_geo=path_geo + geo, loop=loop, deg=deg, reg=reg, path_bvp=path_bvp,
                                  method=method, ptilde=ptilde, rtilde=rtilde, penalty_init=pen)
        path_output = path_bvp.replace("bvp/", "results/")
        path_output = path_output + "-result.xml"
        lib.biharmonic_example(path_bvp + ".xml", path_output)
        file_col.addFile(path_output)
        max_id = max_id + 1
    print("Finished!")
else:
    for geo in geo_list:
        for pen in penalty:
            reg = deg - 1
            path_bvp = "results/" + geo + "/bvp/nitsche/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                       "-l" + str(loop) + "-y" + str(pen)
            path_output = path_bvp.replace("bvp/", "results/")
            path_output = path_output + "-result.xml"
            file_col.addFile(path_output)
            max_id = max_id + 1

        # Approx C1
        ptilde = deg
        rtilde = deg - 1
        reg = deg - 1
        path_bvp = "results/" + geo + "/bvp/nitsche/" + "approxC1" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                   "-P" + str(ptilde) + "-R" + str(rtilde) + "-l" + str(loop)
        path_output = path_bvp.replace("bvp/", "results/")
        path_output = path_output + "-result.xml"
        file_col.addFile(path_output)
        max_id = max_id + 1

        # Nitsche with EW
        pen = -1
        reg = deg - 1
        path_bvp = "results/" + geo + "/bvp/nitsche/" + "nitsche" + "-bvp1-p" + str(deg) + "-r" + str(reg) + \
                   "-l" + str(loop) + "-y" + str(pen)

        path_output = path_bvp.replace("bvp/", "results/")
        path_output = path_output + "-result.xml"
        file_col.addFile(path_output)
        max_id = max_id + 1

file_col.save()

path_tikz = "tikz_files/nitsche/"
path_fig = "tikz_figures/nitsche/"

file1 = gs.io.gsXmlCollection("")
file1.load(xml_col)
list_dict = []
for id in range(max_id):
    my_dict = {"Matrix": file1.getMatrix(id), "Deg": None, "Reg": None,
               "Geo": None, "Name": None}

    path_name = file1.getString(id)
    name = path_name[path_name.find('/results/nitsche/') + 17:]
    my_dict["Name"] = name[:name.find('.xml')]

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
            if dict["Geo"] == geo and int(dict["Deg"]) == deg and "y-1" not in dict["Name"]:
                mat_list.append(dict["Matrix"])
            if dict["Geo"] == geo and int(dict["Deg"]) == deg and "y-1" in dict["Name"]:
                geo_mat_list2.append([dict["Matrix"]])
    geo_mat_list.append(mat_list)
    name_mat_list.append(geo + "-nitsche-penalty-p" + str(deg) + "-r" + str(deg-1) + "-l" + str(loop) + "-h" + str(h))
    name_mat_list2.append(geo + "-nitsche-penalty-jump-p" + str(deg) + "-r" + str(deg-1) + "-l" + str(loop) + "-h" + str(h))

    mat_list = []
    for dict in list_dict:
        if m_str2 in dict["Name"]:
            if dict["Geo"] == geo and int(dict["Deg"]) == deg:
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
    x_list.append(x)

    M = np.zeros((len(mat_list), 3))
    for i, mat in enumerate(mat_list):
        M[i, :] = geo_mat_list_approxC1[idx][0][h,2:5]
    M_list.append(M)
    x_list.append(x)

    a = geo_mat_list2[idx][0][h,4] * 5
    b = geo_mat_list2[idx][0][h,2] * 0.2
    yy = np.linspace(a, b, num=len(mat_list))

    M = np.zeros((len(mat_list), 1))
    x = np.zeros((1, len(mat_list)))
    for i, mat in enumerate(mat_list):
        M[i, :] = yy[i]
        x[0,i] = geo_mat_list2[idx][0][h, x_col]
    M_list.append(M)
    x_list.append(x)

    opt_plot = [{'color': 'green', 'line width': '1pt'},
                {'color': 'blue', 'line width': '1pt'},
                {'color': 'red', 'line width': '1pt'}]

    fig = MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '4.5cm', 'mark options': '{solid}',
                'xlabel': '{Stability parameter $\eta$}', 'ylabel': '{Error}',
                'ylabel style': '{yshift=-0.4cm}', 'xlabel style': '{yshift=0.2cm}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.setPlotOptions(opt_plot)
    fig.create_error_plot(x_list, M_list, False, False)
    fig.generate_tikz(path_tikz + name_mat_list[idx])

    list_tikz.append("nitsche/" + name_mat_list[idx])

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
    x_list.append(x)

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
    x_list.append(x)

    opt_plot = [{'color': 'green', 'line width': '1pt'},
                {'color': 'blue', 'line width': '1pt'},
                {'color': 'red', 'line width': '1pt'}]

    fig = MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '4.5cm', 'mark options': '{solid}',
                'xlabel': '{Stability parameter $\eta$}', 'ylabel': '{Error}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.setPlotOptions(opt_plot)
    fig.create_error_plot(x_list, M_list, False, False)
    fig.generate_tikz(path_tikz + name_mat_list2[idx])

    list_tikz.append("nitsche/" + name_mat_list2[idx])

caption_list = []
caption_list.append('Ex. I: $p=4$, $r=3$')
caption_list.append('Ex. II: $p=4$, $r=3$')
caption_list.append('Ex. III: $p=4$, $r=3$')
caption_list.append('Ex. IV: $p=4$, $r=3$')
caption_list.append('Ex. I: $p=4$, $r=3$')
caption_list.append('Ex. II: $p=4$, $r=3$')
caption_list.append('Ex. III: $p=4$, $r=3$')
caption_list.append('Ex. IV: $p=4$, $r=3$')

doc = MyDocument()
doc.addTikzFigure(list_tikz, caption_list, row=4)
doc.generate_pdf("nitsche_stabilization_example", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
doc.clean_extensions()
