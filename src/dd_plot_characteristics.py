"""
Tool for plotting the results from dd_calculate_real_graph_characteristics.
Run: python ./dd_plot_characteristics.py FILE [--export {pdf}]
"""

import numpy
import pprint
import matplotlib.pyplot as pyplot
import scipy.interpolate

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


args = dd_plot.get_parser().parse_args()
with open('results/' + args.filename + '.txt') as input_file:
    plot_data(input_file, args.filename, args.export)
