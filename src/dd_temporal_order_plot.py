"""
Tool for plotting the results from dd_temporal_order.
Run: python -B ./dd_temporal_order_plot.py FILE [--export {pdf|png}] [--detailed]
Note: FILE is here *without* "-TC.txt" or "-TA.txt" suffix
"""

import os

import matplotlib.pyplot as pyplot
import numpy

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.5

COLORS_TA = ['g', 'm', 'gray', 'lime', \
             'r', 'coral', 'orange', 'yellow', 'b', 'dodgerblue', 'skyblue', 'cyan']
COLORS_TC = ['k', 'b', 'b', 'b', 'g', 'g', 'g']
MARKERS = ['o', 'o', 'o', 'o', 's', 's', 's', 's', '^', '^', '^']
STYLES = ['-', '-', '--', ':', '-', '--', ':']
POINTS = 15

def plot_algorithms(filename, detailed):
    if not os.path.isfile(filename):
        print('Reconstruction algorithms file missing')
        return
    with open(filename) as data_file:
        data = data_file.readlines()
    for line, color, marker in zip(data, COLORS_TA, MARKERS):
        values = line.strip().split(' ')
        density, precision = zip(*[[float(parameter) for parameter in value.split(',')]
                                   for value in values[1:]])
        pyplot.plot(
            [numpy.mean(density)], [numpy.mean(precision)], color = color, marker = marker,
            linestyle = 'None', label = values[0], alpha = 0.7)
        if detailed:
            points = min(POINTS, len(density))
            pyplot.plot(
                [numpy.mean(density[offset::points]) for offset in range(points)],
                [numpy.mean(precision[offset::points]) for offset in range(points)],
                color = color, marker = marker, linestyle = 'None', label = None, alpha = 0.3)

def plot_theoterical_curves(filename):
    if not os.path.isfile(filename):
        print('Theoretical bound file missing')
        return
    with open(filename) as data_file:
        data = data_file.readlines()
    for line, color, style in zip(data, COLORS_TC, STYLES):
        values = line.strip().split(' ')
        density, precision = zip(*[[float(parameter) for parameter in value.split(',')]
                                   for value in values[1:]])
        pyplot.plot(
            density, precision, color = color, linestyle = style, label = values[0], alpha = 0.7)

def plot_labels():
    pyplot.legend(bbox_to_anchor = (0, 1.02, 1, 0.102), loc = 3, ncol = 2, mode = 'expand')
    pyplot.ylabel(r'$\theta$')
    pyplot.xlabel(r'$\epsilon$')
    y_bottom, y_top = min(max(-0.05, sum(pyplot.ylim()) - 1.05), 0.45), 1.05
    pyplot.ylim(y_bottom, y_top)
    x_scale, y_scale = 0.2, round((y_top - y_bottom) / 5, 1)
    pyplot.gca().get_xaxis().set_major_locator(pyplot.MultipleLocator(x_scale))
    pyplot.gca().get_xaxis().set_minor_locator(pyplot.MultipleLocator(x_scale / 2))
    pyplot.gca().get_yaxis().set_major_locator(pyplot.MultipleLocator(y_scale))
    pyplot.gca().get_yaxis().set_minor_locator(pyplot.MultipleLocator(y_scale / 2))

def plot_data(filename, export, detailed):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    plot_theoterical_curves('temp/' + filename + '-TC.txt')
    plot_algorithms('temp/' + filename + '-TA.txt', detailed)
    plot_labels()
    dd_plot.plot(filename + '-T', export)

parser = dd_plot.get_parser()
parser.add_argument('--detailed', action = 'store_true', help = 'plot details')
args = parser.parse_args()
plot_data(args.filename, args.export, args.detailed)
