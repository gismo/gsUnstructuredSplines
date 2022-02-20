#!/usr/bin/python

""""
    @file diss_library.py

    @brief Compute biharmonic2_example using pygismo and plot the output with pylatex to latex/tikz

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import errno
import glob
import os
import subprocess
import sys

from pylatex import Document, NoEscape, Figure, SubFigure, Package, UnsafeCommand, TikZ, TikZOptions, Axis, Plot

gismo_path = os.path.join(os.path.dirname(__file__), "../../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np


def create_biharmonic_bvp(path_geo, loop, deg, reg, path_bvp, method, ptilde=-1, rtilde=-1, penalty_init=-1.0,
                          cond=False):
    # [!Geometry]
    mp = gs.core.gsMultiPatch()
    file = gs.io.gsFileData(path_geo + ".xml")
    file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch
    mp.computeTopology(1e-4, False)
    # [!Geometry]

    # [!Right hand side]
    f = gs.core.gsFunctionExpr("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)
    # [!Right hand side]

    # [!Exact solution]
    ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
    # [!Exact solution]

    # [!Boundary]
    dirichlet = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
    neumann = gs.core.gsFunctionExpr(" -4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                     " -4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2)
    laplace = gs.core.gsFunctionExpr("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)

    bcs = gs.pde.gsBoundaryConditions()
    for bdy in mp.boundaries():
        #               patch_nr, side, boundary condition, function, unknown, parametric, component
        bcs.addCondition(bdy, gs.pde.bctype.dirichlet, dirichlet, 0, False, 0)
        bcs.addCondition(bdy, gs.pde.bctype.laplace, laplace, 0, False, 0)

    # bcs.setGeoMap(mp)
    # [!Boundary]

    # [!Option list]
    opt = gs.io.gsOptionList()
    opt.addSwitch("plot", "Plotting the results.", False)
    opt.addSwitch("info", "Plotting the results.", False)
    # opt.addSwitch("cond", "Estimate condition number (slow!)", cond)

    opt.addInt("numRefine", "Number of refinement-loops.", loop)
    opt.addInt("degree", "Set the polynomial degree of the basis.", deg)
    opt.addInt("smoothness", "Set the smoothness of the basis.", reg)

    opt.addInt("gluingDataDegree", "Set the polynomial degree for the gluing data",
               ptilde)
    opt.addInt("gluingDataSmoothness", "Set the smoothness for the gluing data",
               rtilde)
    opt.addSwitch("interpolation", "Compute the basis constructions with interpolation", False)

    opt.addInt("method", "The chosen method for the biharmonic problem (0 = Approx C1, 1 = Nitsche, 2 = DPatch)",
               method.value)

    opt.addReal("penalty", "Fixed Penalty value for Nitsche's method", penalty_init)
    # [!Option list]

    # [!Save the data to the XML-file]
    file = gs.io.gsFileData()
    file.add(bcs)  # id=0 Boundary
    file.add(f)  # id=1 Source function
    file.add(opt)  # id=2 Optionlist
    file.add(ms)  # id=3 Exact solution
    file.add(mp)  # id=X Geometry (should be last!)
    file.save(path_bvp, False)
    print("Filedata saved: " + path_bvp + ".xml")
    # [!Save the data to the XML-file]


def biharmonic_example(path_bvp, path_output):
    # [!Run biharmonic2_example]
    proc = subprocess.Popen(["../../build/bin/biharmonic2_example", "-o", path_output, "-x", path_bvp])
    proc.wait()
    # [!Run biharmonic2_example]


# input: crop [list] a list of pdfs which has to be cropped
def clean_extensions(crop=None):
    extensions = ['aux', 'log', 'out', 'fls', 'auxlock',
                  'fdb_latexmk', 'md5', 'dpth', 'spl']

    for ext in extensions:
        for f in glob.glob("*" + '.' + ext):
            try:
                os.remove(f)
            except (OSError, IOError) as e:
                if e.errno != errno.ENOENT:
                    raise
    if not crop == None:
        for file in crop:
            try:
                os.system("pdfcrop tikz_figures/" + file + ".pdf tikz_figures/" + file + ".pdf")
            except:
                print("Please install pdfcrop. Cannot crop the figures.")


class MyDocument(Document):
    def __init__(self):
        super().__init__(documentclass='elsarticle', document_options=["a4paper", "3p", NoEscape(r'sort&compress')])

        self.preamble.append(Package("tikz"))
        self.preamble.append(Package("pgfplots"))
        self.preamble.append(NoEscape(r'\usepgfplotslibrary{external}'))
        self.preamble.append(NoEscape(r'\tikzexternalize[prefix=tikz_figures/]'))

        new_comm = UnsafeCommand('newcommand', '\inputtikz', options=1,
                                 extra_arguments=r'\tikzsetnextfilename{#1}'
                                                 r'\input{tikz_files/#1.tikz}')
        self.preamble.append(new_comm)

        self.preamble.append(NoEscape(r'\pgfplotsset{every axis/.append style={label style={font=\tiny},'
                                      r'tick label style={font=\tiny},'
                                      r'axis lines*=left,'
                                      r'axis line style={-},}}'))
        self.preamble.append(NoEscape(r'\pgfplotsset{every add plot/.append style={line width=1pt,}}'))

        self.preamble.append(Package("xcolor"))
        self.preamble.append(NoEscape(r'\definecolor{blue}{HTML}{95DBE5}'))
        self.preamble.append(NoEscape(r'\definecolor{red}{HTML}{078282}'))
        self.preamble.append(NoEscape(r'\definecolor{green}{HTML}{339E66}'))

    def addTikzFigure(self, tikz_list, caption_list, row=1, legend=None):
        # Input: tikz_list      [list] A list of tikz_files
        #        caption_list   [list] A list of captions for each tikz_file

        width = r'' + str(1 / row) + '\\textwidth'
        with self.create(Figure(position='h!')) as fig:
            self.append(NoEscape(r"\centering"))
            for idx, tikz in enumerate(tikz_list):
                with self.create(SubFigure(
                        position='b',
                        width=NoEscape(width))) as subfig:
                    self.append(NoEscape(r"\centering"))
                    if tikz[-3:] == "pdf":
                        self.append(NoEscape(r'\resizebox{\textwidth}{!}{'
                                             r'\includegraphics[width=\textwidth]{tikz_figures/' + tikz + '}'))
                        self.append(NoEscape(r'}'))
                    else:
                        self.append(NoEscape(r'\inputtikz{' + tikz + '}'))
                    subfig.add_caption(NoEscape(r'' + caption_list[idx]))
                if idx % row == row - 1:
                    self.append(NoEscape(r"\\"))

            if legend != None:
                self.append(NoEscape(r"\begin{center}"))
                for leg in legend:
                    with self.create(SubFigure(
                            position='b',
                            width=NoEscape(width))) as subfig:
                        self.append(NoEscape(r"\centering"))
                        self.append(NoEscape(r'\inputtikz{' + leg + '}'))
                self.append(NoEscape(r"\end{center}"))

            fig.add_caption("Two kittens")


class MyTikz(Document):
    def __init__(self):
        super().__init__(documentclass='standalone')
        self.opt = {}
        self.opt_plot = [{'mark': 'diamond*', 'color': 'green', 'line width': '1pt'},
                         {'mark': 'square*', 'color': 'blue', 'line width': '1pt'},
                         {'mark': 'triangle*', 'color': 'red', 'line width': '1pt'},
                         {'mark': 'pentagon*', 'color': 'yellow', 'line width': '1pt'},
                         {'mark': 'halfcircle*', 'color': 'brown', 'line width': '1pt'}]
        self.opt_plot = self.opt_plot * 10

        self.opt_mat = [["solid"], ["dashed"], ["dotted"], ["dashdotted"], ["densely dotted"]]

        self.legend = []

        self.deg_list = []
        self.rate_shift = 0

    def setDegree(self, deg):
        self.deg_list = deg

    def setRateShift(self, shift):
        self.rate_shift = shift

    def setOptions(self, options):
        self.opt = options

    def setPlotOptions(self, options):
        self.opt_plot = options

    def setColor(self, color):
        for col, opt in zip(color, self.opt_plot):
            opt["color"] = col

    def setMarker(self, marker):
        for mark, opt in zip(marker, self.opt_plot):
            opt["mark"] = mark

    def getLegendList(self):
        return self.legend

    def create_legend(self, legend_image, legend_entry):
        self.opt = ["hide axis", "xmin=0", "xmax=1", "ymin=0", "ymax=0.4", NoEscape("mark options={solid}"),
                    NoEscape(r"legend style={draw=white!15!black,legend cell align=left}"),
                    "transpose legend",
                    NoEscape(r"legend columns=" + str(int(len(legend_image) / 2)) +
                             ",legend style={/tikz/every even column/.append style={column sep=0.5cm}}")
                    ]
        with self.create(TikZ()) as tikz:
            with self.create(Axis(options=self.opt)) as axis:
                zip_object = zip(legend_image, legend_entry)
                for image, entry in zip_object:
                    axis.append(NoEscape(r'\addlegendimage{' + str(image[0]) + '}'))
                    axis.append(NoEscape(r'\addlegendentry{' + str(entry[0]) + '}'))

    def generate_tikz(self, filename):
        tex = self.dumps()  # The document as string in LaTeX syntax
        with open(filename + ".tikz", 'w') as f:
            begin = False
            for idx, line in enumerate(tex.splitlines()):
                if line == "\\begin{tikzpicture}%":
                    begin = True
                if begin and idx != len(tex.splitlines()) - 1:
                    f.write(line)
                    f.write("\n")

    '''
    x_list = List of x values for the x-axis
    M_list = List of matrices: rows -> number of points (== length of x)
                               cols -> number of lines 
                               
    Plotsyle: For each entry in list -> new line style [["solid"], ["dashed"], ...] see self.opt_mat 
              For each col in matrix: new color, marker, ... see self.opt_plot
    '''

    def create_error_plot(self, x_list, M_list, array=False, rate=False):
        # Number of matrices -> List of List of List:
        # [[[x_1^1,y_1^1],...[x_n^1,y_n^1]],...,[[x_1^m,y_1^m],...[x_n^m,y_n^m]]]
        points_list = []
        rates_list = []
        opt_rate = []
        for idx, mat in enumerate(M_list):  # For each entry in list
            points = []
            if array:  # if mat is an array
                points_temp = []  # Number of cols
                for i, j in zip(x_list[idx], mat[:]):
                    points_temp.append([i, j])
                points.append(points_temp)
                points_list.append(points)
                if rate:
                    rates_list.append(
                        np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[idx]) + self.rate_shift) * np.log(
                            points_temp[-1][0])))
                    opt_rate.append(
                        {NoEscape(r'domain={' + str(x_list[idx][-2]) + ':' + str(x_list[idx][-1]) + '}'), 'draw=black',
                         'samples=100',
                         'forget plot', NoEscape(r'yshift=-0.2cm'), 'line width=1pt'})

            else:  # if mat is a matrix
                for col in range(mat.shape[1]):
                    points_temp = []  # Number of cols
                    for i, j in zip(x_list[idx], mat[:, col]):
                        points_temp.append([i, j])
                    points.append(points_temp)
                    if rate and idx == 0:
                        rates_list.append(
                            np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[col]) + self.rate_shift) * np.log(
                                points_temp[-1][0])))
                        opt_rate.append(
                            {NoEscape(r'domain={' + str(x_list[idx][-2]) + ':' + str(x_list[idx][-1]) + '}'),
                             'draw=black',
                             'samples=100',
                             'forget plot', NoEscape(r'yshift=-0.2cm'), 'line width=1pt'})
                points_list.append(points)

        opt_axis = TikZOptions(**self.opt, width=NoEscape(r'1\textwidth'))
        with self.create(TikZ()) as tikz:
            with self.create(Axis(options=opt_axis)) as axis:
                for idx, mat in enumerate(points_list):  # idx == id of matrix
                    for col in range(len(mat)):  # len(mat) == number of cols

                        # Define new line style for each matrix
                        new_list = []
                        for key, val in self.opt_plot[col].items():
                            new_list.append(str(key) + "=" + str(val))
                        new_list.append(self.opt_mat[idx][0])

                        curve = Plot(options=new_list, coordinates=mat[col])
                        axis.append(curve)

                        self.legend.append([",".join(new_list)])

                        if rate and abs(mat[col][-1][1]) > 1e-12:
                            # rate_h2 = np.log(points[col][-2][1]/points[col][-1][1])/np.log(2.0)
                            # print(np.around(rate_h2,2))
                            if array:
                                curve = Plot(options=opt_rate[idx],
                                             func=str(rates_list[idx]) + '*x^' + str(
                                                 int(self.deg_list[idx]) + self.rate_shift))

                                string_curve = curve.dumps()
                                string_curve = string_curve[:string_curve.find(';%')]
                                string_curve += NoEscape(r' node[right, pos=0.75] {\tiny{$h' + str(
                                    int(self.deg_list[idx]) + self.rate_shift) + '$}};')
                                axis.append(NoEscape(r'' + string_curve))
                            elif idx == 0:  # To avoid double printing
                                curve = Plot(options=opt_rate[col],
                                             func=str(rates_list[col]) + '*x^' + str(
                                                 int(self.deg_list[col]) + self.rate_shift))

                                string_curve = curve.dumps()
                                string_curve = string_curve[:string_curve.find(';%')]
                                string_curve += NoEscape(r' node[right, pos=0.75] {\tiny{$h' + str(
                                    int(self.deg_list[col]) + self.rate_shift) + '$}};')
                                axis.append(NoEscape(r'' + string_curve))
