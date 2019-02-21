"""
Tool for plotting the results from dd_temporal_order.
Run: python -B ./dd_temporal_order.py FILE [--export {pdf}]
Note: FILE is here *without* "-TC.txt" or "-TA.txt" suffix
"""

import os

import matplotlib.pyplot as pyplot
import numpy

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.5

COLORS = ['r', 'b', 'm', 'k', 'g', 'orange', 'crimson', 'lime', 'gray', 'lightgray', 'pink', \
          'olive', 'khaki', 'saddlebrown', 'deepskyblue']

def plot_algorithms(filename):
    if not os.path.isfile(filename):
        return
    with open(filename) as data_file:
        data = data_file.readlines()
    for line, color in zip(data, COLORS):
        values = line.strip().split(' ')
        density, precision = zip(*[[float(parameter) for parameter in value.split(',')]
                                   for value in values[1:]])
        pyplot.plot(
            [numpy.mean(density)], [numpy.mean(precision)], color = color, marker = 'o',
            linestyle = None, label = values[0], alpha = 0.7)
        pyplot.plot(
            density, precision, color = color, marker = 'o',
            linestyle = None, label = None, alpha = 0.3)

def plot_theoterical_curves(filename):
    if not os.path.isfile(filename):
        return
    with open(filename) as data_file:
        data = data_file.readlines()
    for line, color in zip(data, COLORS):
        values = line.strip().split(' ')
        density, precision = zip(*[[float(parameter) for parameter in value.split(',')]
                                   for value in values[1:]])
        pyplot.plot(density, precision, color = color, label = values[0], alpha = 0.7)

def plot_labels():
    pyplot.legend(loc = 'upper left')
    pyplot.ylabel(r'$\theta$')
    pyplot.xlabel(r'$\epsilon$')
    y_bottom, y_top = max(-0.05, sum(pyplot.ylim()) - 1.0), 1.05
    pyplot.ylim(y_bottom, y_top)
    x_scale, y_scale = 0.2, round((y_top - y_bottom) / 5, 1)
    pyplot.gca().get_xaxis().set_major_locator(pyplot.MultipleLocator(x_scale))
    pyplot.gca().get_xaxis().set_minor_locator(pyplot.MultipleLocator(x_scale / 2))
    pyplot.gca().get_yaxis().set_major_locator(pyplot.MultipleLocator(y_scale))
    pyplot.gca().get_yaxis().set_minor_locator(pyplot.MultipleLocator(y_scale / 2))

def plot_data(filename, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    plot_theoterical_curves('temp/' + filename + '-TC.txt')
    plot_algorithms('temp/' + filename + '-TA.txt')
    plot_labels()
    dd_plot.plot(filename, export)

args = dd_plot.get_parser().parse_args()
plot_data(args.filename, args.export)
