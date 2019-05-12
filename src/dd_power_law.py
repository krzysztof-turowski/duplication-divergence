import argparse
import re

import matplotlib.pyplot as pyplot
import networkx
import powerlaw

import dd_header
import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.75

def get_fitness(filename, G):
    degree_list = [degree for _, degree in dict(G.degree()).items() if degree > 0]
    fit = powerlaw.Fit(degree_list, discrete = True)
    gamma = fit.power_law.alpha
    cutoff = fit.power_law.xmin
    percentile = 100 * (
        1.0 - float(len(list(d for d in degree_list if d >= fit.power_law.xmin))) / G.order())
    print('{0}: gamma {1}, cutoff {2}, percentile {3:.2f}'.format(
        filename, gamma, cutoff, percentile))
    return gamma, percentile

def plot_powerlaw(G):
    degree_list = [degree for _, degree in dict(G.degree()).items() if degree > 0]
    figure = pyplot.subplot(111)
    fit = powerlaw.Fit(degree_list, discrete = True)
    fit.power_law.plot_ccdf(color = 'r', linestyle = '--', ax = figure, label = 'with cutoff')
    all_fit = powerlaw.Fit(degree_list, xmin = 1, discrete = True)
    all_fit.power_law.plot_ccdf(color = 'g', linestyle = ':', ax = figure, label = 'original')
    all_fit.plot_ccdf(color = 'g', linewidth = 2, ax = figure)

def plot_labels():
    pyplot.xlabel('$deg$')
    pyplot.ylabel('CCDF($deg$)')
    pyplot.legend(loc = 'lower left', prop = {'size': 12})

def test_powerlaw(graph_name, action, argument_list):
    if action == 'real_graph':
        G = networkx.read_edgelist('files/' + graph_name + '.txt', nodetype = int)
        get_fitness(graph_name, G)
        if argument_list.plot:
            dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
            plot_powerlaw(G)
            plot_labels()
            dd_plot.plot(graph_name + '-P', 'pdf')
    elif action == 'real_seed':
        seed_name = re.sub('^G-', 'G0-', graph_name, flags = re.MULTILINE)
        TRIES, mode, p, r = argument_list.gt, argument_list.mode, argument_list.p, argument_list.r
        G0 = networkx.read_edgelist('files/' + seed_name + '.txt', nodetype = int)
        n = networkx.read_edgelist('files/' + graph_name + '.txt', nodetype = int).order()
        gamma_avg, percentile_avg = 0, 0
        for _ in range(TRIES):
            if mode == 'pastor_satorras':
                G = dd_header.generate_pastor_satorras(G0.copy(), n, G0.order(), p, r)
            else:
                raise ValueError('Unknown mode: {0}'.format(mode))
            gamma, percentile = get_fitness(graph_name, G)
            gamma_avg += gamma
            percentile_avg += percentile
        print('{0}: average gamma {1}, average percentile {2:.2f}'.format(
            graph_name, gamma_avg / TRIES, percentile_avg / TRIES))
    else:
        raise ValueError('Unknown action: {0}'.format(action))

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', metavar = 'FILE', help = 'graph filename')
    parser.add_argument('-action', choices = ['real_graph', 'real_seed'])
    parser.add_argument('-gt', type = int, help = 'number of tries')
    parser.add_argument('-mode', choices = ['pastor_satorras'])
    parser.add_argument('-p', type = float, help = 'parameter p')
    parser.add_argument('-r', type = float, help = 'parameter r')
    parser.add_argument('--plot', action = 'store_true', help = 'plot results')
    return parser

args = get_parser().parse_args()
test_powerlaw(args.filename, args.action, args)
