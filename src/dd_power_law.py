import argparse

import matplotlib.pyplot as pyplot
import networkx
import powerlaw

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.75

def plot_powerlaw(filename, G, plot):
    degree_list = [degree for _, degree in dict(G.degree()).items()]
    fit = powerlaw.Fit(degree_list, discrete = True)
    print('{0}: gamma {1}, cutoff {2}, percentile {3:.2f}'.format(
        filename, fit.power_law.alpha, fit.power_law.xmin,
        100 * (1.0 - len(list(d for d in degree_list if d >= fit.power_law.xmin)) / G.order())))

    if plot:
        figure = pyplot.subplot(111)
        fit.power_law.plot_ccdf(color = 'r', linestyle = '--', ax = figure, label = 'with cutoff')
        all_fit = powerlaw.Fit(degree_list, xmin = 1, discrete = True)
        all_fit.power_law.plot_ccdf(color = 'g', linestyle = ':', ax = figure, label = 'original')
        all_fit.plot_ccdf(color = 'g', linewidth = 2, ax = figure)

def plot_labels():
    pyplot.xlabel('$deg$')
    pyplot.ylabel('CCDF($deg$)')
    pyplot.legend(loc = 'lower left', prop = {'size': 12})

def plot_data(filename, action, plot):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    if action == 'real_graph':
        G = networkx.read_edgelist('files/' + filename + '.txt', nodetype = int)
    else
        raise
    plot_powerlaw(filename, G, plot)
    if plot:
        plot_labels()
        dd_plot.plot(filename + '-P', 'pdf')

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', metavar = 'FILE', required = True, help = 'graph file')
    parser.add_argument('--action', choices = ['real_graph', 'real_seed', 'synthetic'])
    parser.add_argument('--gt', type = int, help = 'number of tries')
    parser.add_argument('--plot', action = 'store_true', help = 'plot results')
    return parser

args = get_parser().parse_args()
plot_data(args.filename, args.action, args.plot)
