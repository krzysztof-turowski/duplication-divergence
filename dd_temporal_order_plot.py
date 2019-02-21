"""
Tool for plotting the results from dd_temporal_order.
Run: python -B ./dd_temporal_order.py FILE [--export {pdf}]
Note: FILE is here *without* "-TC.txt" suffix
"""

import os

import matplotlib.pyplot as pyplot

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.5

COLORS = ['r', 'b', 'm', 'k', 'g']

def plot_theoterical_curves(filename):
    if not os.path.isfile(filename):
        return
    with open(filename) as file:
        data = file.readlines()
    for line, color in zip(data, COLORS):
        name, *values = line.strip().split(' ')
        X, Y = zip(*[[float(parameter) for parameter in value.split(',')] for value in values])
        pyplot.plot(X, Y, color = color, label = name, alpha = 0.7)

def plot_labels():
    pyplot.legend(loc = 'upper left')
    pyplot.ylabel(r'$\theta$')
    pyplot.xlabel(r'$\epsilon$')
    pyplot.gca().get_xaxis().set_major_locator(pyplot.MultipleLocator(0.2))
    pyplot.gca().get_xaxis().set_minor_locator(pyplot.MultipleLocator(0.1))
    y_bottom, y_top = max(-0.05, sum(pyplot.ylim()) - 1.0), 1.05
    y_scale = round((y_top - y_bottom) / 5, 1)
    pyplot.ylim(y_bottom, y_top)
    pyplot.gca().get_yaxis().set_major_locator(pyplot.MultipleLocator(y_scale))
    pyplot.gca().get_yaxis().set_minor_locator(pyplot.MultipleLocator(y_scale / 2))

def plot_data(filename, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    plot_theoterical_curves('temp/' + filename + '-TC.txt')
    plot_labels()
    dd_plot.plot(filename, export)

args = dd_plot.get_parser().parse_args()
plot_data(args.filename, args.export)
