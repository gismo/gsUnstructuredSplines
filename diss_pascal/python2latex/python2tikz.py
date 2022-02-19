from pylatex import Document, NoEscape, TikZ, TikZOptions, Axis, Plot

import numpy as np


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

        self.opt_mat = [["solid"], ["dashed"], ["dotted"], ["dashdotted"],["densely dotted"]]

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
                if (begin and idx != len(tex.splitlines()) - 1):
                    f.write(line)
                    f.write("\n")

    # def create_error_table(self, x, result, result_firstRow):
    #
    #     l2error_col = 2
    #     h1error_col = 4
    #     h2error_col = 6
    #
    #     coord = []
    #     coord2 = []
    #     coord3 = []
    #     for i, j in zip(x, result[:, l2error_col]):
    #         coord.append([i, j])
    #     for i, j in zip(x, result[:, h1error_col]):
    #         coord2.append([i, j])
    #     for i, j in zip(x, result[:, h2error_col]):
    #         coord3.append([i, j])
    #
    #     with self.create(Tabular('c|c||c|c||c|c||c|c')) as tabular:
    #         self.append(NoEscape(r'\centering'))
    #         tabular.add_row((result_firstRow))
    #         tabular.add_hline()
    #         for xx, y0, y1, r1, y2, r2, y3, r3 in zip(x, result[:, 1], result[:, l2error_col],
    #                                                   result[:, l2error_col + 1],
    #                                                   result[:, h1error_col], result[:, h1error_col + 1],
    #                                                   result[:, h2error_col], result[:, h2error_col + 1]):
    #             tabular.add_row((xx, int(y0), y1, r1, y2, r2, y3, r3))
    #             tabular.add_hline()

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
                if rate and idx == 0:
                    rates_list.append(
                        np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[idx]) + self.rate_shift) * np.log(points_temp[-1][0])))
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
                            np.exp(np.log(points_temp[-1][1]) - (int(self.deg_list[col]) + self.rate_shift) * np.log(points_temp[-1][0])))
                        opt_rate.append(
                            {NoEscape(r'domain={' + str(x_list[idx][-2]) + ':' + str(x_list[idx][-1]) + '}'), 'draw=black',
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

                        if rate and abs(mat[col][-1][1]) > 1e-12 and idx == 0:
                            #rate_h2 = np.log(points[col][-2][1]/points[col][-1][1])/np.log(2.0)
                            #print(np.around(rate_h2,2))

                            curve = Plot(options=opt_rate[col],
                                         func=str(rates_list[col]) + '*x^' + str(int(self.deg_list[col]) + self.rate_shift),
                                         label=NoEscape(r'\tiny{$h^'+str(int(self.deg_list[col]) + self.rate_shift) + '$}'))
                            axis.append(curve)
