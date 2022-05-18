"""
Tool for plotting the results from dd_calculate_real_graph_characteristics.
Run: python ./dd_plot_characteristics.py FILE [--export {pdf}]
"""

import numpy
import math
import pprint
import matplotlib.pyplot as pyplot
import scipy.interpolate
import os

import dd_plot

PLOT_STYLE = 'default'
FIGURE_SIZE_SCALE = 1.3


def read_data(data, interpolate=False):
    result = {}
    for line in data:
        name, data = line.split(":")

        name = name.replace(" ", "_")
        data = data.strip().split(" ")
        if len(data) == 1:
            result[name] = float(data[0])
        else:
            result[name] = [float(x) for x in data]

    return result


def plot_single(sub, title, name, data, pre=None, numbers=True):
    pyplot.subplot(sub)
    if pre is not None:
        pre()
    if name not in data:
        return

    pyplot.title(title)
    pyplot.bar([x + 1 if numbers else str(x + 1)
                for x in range(len(data[name]))], data[name])


def plot_data(data_file, filename, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    print("Get data: ", filename)
    data = read_data(data_file)
    pyplot.subplots_adjust(hspace=0.3, wspace=0.2)

    print("Sorting: ", filename)
    data["closeness"] = sorted(data["closeness"])
    data["betweenness_centrality"] = sorted(data["betweenness_centrality"])

    def pre():
        pyplot.yscale('log')
    plot_single(221, "Graphlets", "Graphlets", data, pre, False)
    plot_single(222, "Relative Graphlet Frequency", "RGF", data, numbers=False)

    def pre():
        pyplot.yscale('log')
    plot_single(223, "Betweenness", "betweenness_centrality", data, pre)
    plot_single(224, "Closeness", "closeness", data, pre)

    print("Plotting rgf: ", filename)
    dd_plot.plot(filename + '-rgf', export)

    pyplot.clf()

    for i in range(1, 4):
        plot_single(
            220 + i,
            f"$k$-hop reachability {i+1}",
            f"khop_reachability_{i+1}",
            data)

    def pre():
        pyplot.yscale('log')

    plot_single(
        224,
        "Degree Distribution",
        "get_degree_distribution",
        data,
        pre)

    print("Plotting deg: ", filename)
    dd_plot.plot(filename + '-deg', export)
    print("Done: ", filename)


def plot_data_cmp1d(data_file, graph, synthetic, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    data = []
    target = read_data(data_file)

    def get_param_value(filename):
        param = filename.split('-')[-1][:-4]
        return float(param)
    for filename in synthetic:
        with open(f'results/synthetic/{filename}') as file:
            data.append((get_param_value(filename), read_data(file)))

    data = sorted(data, key=lambda x: x[0])

    METRICS = [("Automorphisms", "log_automorphisms_sparse"),
               ("Average shortest path", "get_average_shortest_path"),
               ("Diameter", "get_diameter"),
               ("Clustering coefficient", "clustering_coefficient_two")]

    for i, (name, metric) in enumerate(METRICS, 1):
        pyplot.subplot(220 + i)
        pyplot.title(name)
        pyplot.plot([v for v, d in data], [d[metric] for v, d in data])
        pyplot.axhline(target[metric], color="r")
    dd_plot.plot(graph + '-cmp1d', export)


def plot_data_cmp2d(data_file, graph, synthetic, export):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    data = []
    target = read_data(data_file)
    pyplot.subplots_adjust(hspace=0.4)

    def get_params_value(filename):
        params = filename.split('-')[-2:]
        x, y = params[0], params[1][:-4]
        return float(x), float(y)
    for filename in synthetic:
        with open(f'results/synthetic/{filename}') as file:
            data.append((get_params_value(filename), read_data(file)))

    data = sorted(data, key=lambda x: x[0])
    side_length = int(math.sqrt(len(data)))

    METRICS = [("Automorphisms", "log_automorphisms_sparse"),
               ("Average shortest path", "get_average_shortest_path"),
               ("Diameter", "get_diameter"),
               ("Clustering coefficient", "clustering_coefficient_two")]

    for i, (name, metric) in enumerate(METRICS, 1):
        pyplot.subplot(220 + i)
        pyplot.title(name)
        pyplot.xticks(range(side_length),
                      [*set([v[0] for v, d in data])], rotation=-90)
        pyplot.yticks(range(side_length),
                      [*set([v[1] for v, d in data])])
        pyplot.imshow(numpy.array([math.fabs(d[metric] - target[metric])
                      for v, d in data]).reshape((side_length, side_length)))
        pyplot.colorbar()
    dd_plot.plot(graph + '-cmp2d', export)


args = dd_plot.get_parser().parse_args()
with open('results/' + args.filename + '.txt') as input_file:
    if args.mode == 'single':
        plot_data(input_file, args.filename, args.export)
    elif args.mode == 'cmp1d':
        models1d = ['PD']
        for model in models1d:
            synthetic = [file for file in os.listdir(
                "results/synthetic") if file.startswith(args.filename) and model in file]

            plot_data_cmp1d(input_file, args.filename, synthetic, args.export)
    elif args.mode == 'cmp2d':
        models2d = ['PS']
        for model in models2d:
            synthetic = [file for file in os.listdir(
                "results/synthetic") if file.startswith(args.filename) and model in file]

            plot_data_cmp2d(input_file, args.filename, synthetic, args.export)
    else:
        print('Unknown mode', args.mode)
        exit(1)
