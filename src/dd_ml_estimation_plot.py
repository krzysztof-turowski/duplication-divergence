"""
Tool for plotting the results from dd_ml_estimation.
Run: python -B ./dd_ml_estimation_plot.py FILE [--export {pdf}]
"""

import os

import numpy
import matplotlib.pyplot as pyplot
import scipy.interpolate

import dd_plot

PLOT_STYLE = 'default'
FIGURE_SIZE_SCALE = 0.6
X_LABELS, Y_LABELS = 11, 6
RESOLUTION = 20

PASTOR_SATORRAS_PTRN = [False, True, False]
CHUNG_LU_PTRN = [False, False, True]

def read_data(data, interpolate = False):
    values = [[float(value) if value != '' else None for value in parameter.split(',')]
              for parameter in data.strip().split(' ')]
    if (numpy.equal(values[0][:len(PASTOR_SATORRAS_PTRN)], None) == PASTOR_SATORRAS_PTRN).all():
        pyplot.xlabel(r'$p$')
        pyplot.ylabel(r'$r$', rotation = 0)
        x, _, y, z = zip(*values)
    elif (numpy.equal(values[0][:len(CHUNG_LU_PTRN)], None) == CHUNG_LU_PTRN).all():
        pyplot.xlabel(r'$p$')
        pyplot.ylabel(r'$q$', rotation = 0)
        x, y, _, z = zip(*values)
    else:
        raise Exception('Unidentified type of data: {0}'.format(numpy.equal(values[0], None)))

    if interpolate:
        # TODO: fix issues with infinity - currently we replace it with 2 * minimum
        minimum = min(value for value in z if value != -float('inf'))
        z = [value if value >= minimum else 2 * minimum for value in z]
        f = scipy.interpolate.interp2d(x, y, z, kind = 'cubic')
        x_range = numpy.linspace(min(x), max(x), num = RESOLUTION, endpoint = True)
        y_range = numpy.linspace(min(y), max(y), num = RESOLUTION, endpoint = True)
        matrix = [[value if value >= minimum else -float('inf') for value in row]
                  for row in f(x_range, y_range)]
    else:
        x_range, y_range = numpy.unique(x), numpy.unique(y)
        matrix = numpy.zeros((len(y_range), len(x_range)))
        for px, py, pz in zip(x, y, z):
            matrix[numpy.searchsorted(y_range, py), numpy.searchsorted(x_range, px)] = pz
    matrix = numpy.flip(matrix, axis = 0)
    return x_range, y_range, matrix

def plot_data(data_file, filename, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    x, y, data = read_data(data_file.readlines()[0])
    image = pyplot.imshow(data, cmap = 'BuPu', interpolation = 'nearest', aspect = 'auto')
    # pcolormesh, cmap {'YlGnBu','Pastel1'}, interpolation {'nearest','bilinear'}, aspect {'auto'}

    labels_x = numpy.linspace(round(min(x)), round(max(x)), num = X_LABELS, endpoint = True)
    pyplot.gca().set_xticks(numpy.linspace(0, len(x) - 1, num = X_LABELS, endpoint = True))
    pyplot.gca().set_xticklabels([(round(float(label), 2)) for label in labels_x])
    labels_y = numpy.flip(numpy.linspace(round(min(y)), round(max(y)),
                                         num = Y_LABELS, endpoint = True))
    pyplot.gca().set_yticks(numpy.linspace(0, len(y) - 1, num = Y_LABELS, endpoint = True))
    pyplot.gca().set_yticklabels([(round(float(label), 1)) for label in labels_y])
    pyplot.gcf().colorbar(image)
    dd_plot.plot(filename + '-ML', export)

args = dd_plot.get_parser().parse_args()
with open('temp/' + args.filename + '-ML.txt') as input_file:
    plot_data(input_file, args.filename, args.export)
