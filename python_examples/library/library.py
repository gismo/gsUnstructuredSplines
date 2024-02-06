#!/usr/bin/python

""""
    @file library.py
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

from pylatex import Document, NoEscape, Figure, SubFigure, Package, UnsafeCommand, TikZ, TikZOptions, Axis, Plot

import numpy as np

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
            if file[-3:] == "png":
                try:
                    os.system("convert tikz_figures/" + file + " -trim tikz_figures/" + file)
                except:
                    print("Crop png file doesn't work.")
            else:
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

    def addTikzFigure(self, tikz_list, subcaption_list, caption_list, col=1, row=7, legend=None):
        # Input: tikz_list         [list] A list of tikz_files
        #        subcaption_list   [list] A list of captions for each tikz_file

        pages = len(tikz_list)//(col*row) if len(tikz_list)%(col*row) == 0 else len(tikz_list)//(col*row) + 1

        width = r'' + str(1 / col) + '\\textwidth'
        for page in range(pages):
            with self.create(Figure(position='h!')) as fig:
                self.append(NoEscape(r"\centering"))
                if not page ==0:
                    self.append(NoEscape(r"\ContinuedFloat"))
                for idx, tikz in enumerate(tikz_list[col*row*page:col*row*(page+1)]):
                    if not tikz == "":
                        with self.create(SubFigure(
                                position='b',
                                width=NoEscape(width))) as subfig:
                            self.append(NoEscape(r"\centering"))
                            if tikz[-3:] == "pdf" or tikz[-3:] == "png":
                                self.append(NoEscape(r'\resizebox{\textwidth}{!}{'
                                                     r'\includegraphics[width=\textwidth]{tikz_figures/' + tikz + '}'))
                                self.append(NoEscape(r'}'))
                            else:
                                self.append(NoEscape(r'\inputtikz{' + tikz + '}'))
                            subfig.add_caption(NoEscape(r'' + subcaption_list[idx+col*row*page]))
                    if idx % col == col - 1:
                        self.append(NoEscape(r"\\"))

                if legend != None:
                    self.append(NoEscape(r"\begin{center}"))
                    for leg in legend:
                        with self.create(SubFigure(
                                position='b',
                                width=NoEscape(r'1\textwidth'))) as subfig:
                            self.append(NoEscape(r"\centering"))
                            self.append(NoEscape(r'\inputtikz{' + str(leg) + '}'))
                    self.append(NoEscape(r"\end{center}"))

                fig.add_caption(caption_list[page])


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

    def setDegree(self, deg):
        self.deg_list = deg

    def setOptions(self, options):
        self.opt = options

    def setPlotOptions(self, options):
        self.opt_plot = options

    def swapLinestyle(self, id):
        self.opt_mat[0] = self.opt_mat[id]

    def setColor(self, color):
        for col, opt in zip(color, self.opt_plot):
            opt["color"] = col

    def setMarker(self, marker):
        for mark, opt in zip(marker, self.opt_plot):
            opt["mark"] = mark

    def getLegendList(self):
        return self.legend

    def create_legend(self, legend_image, legend_entry, col=2):
        self.opt = ["hide axis", "xmin=0", "xmax=1", "ymin=0", "ymax=0.4", NoEscape("mark options={solid}"),
                    NoEscape(r"legend style={draw=white!15!black,legend cell align=left}"),
                    "transpose legend",
                    NoEscape(r"legend columns=" + str(int((len(legend_image)+1) / (col+1))) +
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
                        np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[idx])) * np.log(
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
                            np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[col])) * np.log(
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
                            #rate_h2 = np.log(points[col][-2][1]/points[col][-1][1])/np.log(2.0)
                            #print(np.around(rate_h2,2))

                            if array:
                                curve = Plot(options=opt_rate[idx],
                                             func=str(rates_list[idx]) + '*x^' + str(
                                                 int(self.deg_list[idx])))

                                string_curve = curve.dumps()
                                string_curve = string_curve[:string_curve.find(';%')]
                                string_curve += NoEscape(r' node[right, pos=0.75] {\tiny{$h' + str(
                                    int(self.deg_list[idx])) + '$}};')
                                axis.append(NoEscape(r'' + string_curve))
                            elif idx == 0:  # To avoid double printing
                                curve = Plot(options=opt_rate[col],
                                             func=str(rates_list[col]) + '*x^' + str(
                                                 int(self.deg_list[col])))

                                string_curve = curve.dumps()
                                string_curve = string_curve[:string_curve.find(';%')]
                                string_curve += NoEscape(r' node[right, pos=0.75] {\tiny{$h' + str(
                                    int(self.deg_list[col])) + '$}};')
                                axis.append(NoEscape(r'' + string_curve))