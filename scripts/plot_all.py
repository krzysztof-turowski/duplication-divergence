import os
import json
from typing import List, Dict, Any
from dd_plot_characteristics import read_data, PLOT_STYLE, FIGURE_SIZE_SCALE
from itertools import combinations
import matplotlib.pyplot as pyplot
import dd_plot
import numpy
import math

KEYS = [
    ('log_automorphisms_sparse', "Automorphisms"),
    ('get_average_shortest_path', "Average shortest path"),
    ('clustering_coefficient_two', "Clustering coefficient 2"),
    ('clustering_coefficient_three', "Clustering coefficient 3"),
    ('clustering_coefficient_four', "Clustering coefficient 4"),
    ('get_diameter', "Diameter"),
    ('get_degree_distribution', "Degree Distribution"),
    ('betweenness_centrality', "Betweenness"),
    ('closeness', "Closeness"),
    ('khop_reachability_2', "$k$-hop reachability 2"),
    ('khop_reachability_3', "$k$-hop reachability 3"),
    ('khop_reachability_4', "$k$-hop reachability 4"),
    ('Graphlets', "Graphlets"),
    ('RGF', "Relative Graphlet Frequency"),
]


def get_all_param_pairs(model, params):
    if len(params) == 1 or model == "2STEP":
        return [(model, params[0])]
    if model == "BERG":
        params = params[2:]

    return [(f"{model}_{a[0]}_{b[0]}", a[1], b[1])
            for a, b in combinations(enumerate(params), 2)]


def get_average(results: List[dict]) -> Dict[str, float]:
    averaged = {}
    for key, _ in KEYS:
        number = 0
        summed = 0
        for result in results:
            result = result["comparison"]
            if key not in result:
                continue
            number += 1
            summed += result[key]
        if number > 0:
            averaged[key] = summed / number

    return averaged


def get_averaged_dict():
    tree = os.walk('results/compared')

    db: Dict[str, List[Any]] = {}

    grouped = {}

    for (path, _, files) in tree:
        for file in files:
            with open(path + "/" + file) as f:
                data = json.load(f)

            _, _, model = file[:-4].split("%3Amode%3A")[0].split("_")
            _, _, model, *_ = model[2:].split("-")

            param_pairs = get_all_param_pairs(model, data["params"])
            for pair in param_pairs:
                if pair not in grouped:
                    grouped[pair] = []
                grouped[pair].append(data)
    return {k: get_average(v) for k, v in grouped.items()}


def plot_data_cmp1d(key, data):
    for metric, title in KEYS:
        pyplot.title(title)
        pyplot.plot([v for v, d in data], [d[metric] for v, d in data])
        dd_plot.plot(f"{key}-cmp1d-{metric}", 'pdf')
        pyplot.clf()


def plot_data_cmp2d(key, data):
    left = [*set([v[0] for v, d in data])]
    right = [*set([v[1] for v, d in data])]
    for metric, title in KEYS:
        dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
        pyplot.title(title)
        pyplot.xticks(range(len(left)), left)
        pyplot.yticks(range(len(right)), right)
        pyplot.imshow(
            numpy.array([
                d[metric] if metric in d else 0 for v, d in data
            ]).reshape((len(right), len(left))))
        pyplot.colorbar()
        dd_plot.plot(f"{key}-cmp2d-{metric}", 'pdf')
        pyplot.clf()


def plot_all(averaged_dict):
    results = sorted(averaged_dict.items())
    last_key = None
    current = []
    for ((key, *params), data) in results:
        if last_key != key:
            if len(current) > 0:
                if len(current[-1][0]) > 1:
                    plot_data_cmp2d(last_key, current)
                else:
                    plot_data_cmp1d(last_key, current)
            current = []
            last_key = key
        current.append((params, data))
    if len(current) > 0:
        if len(current[-1][0]) > 1:
            plot_data_cmp2d(last_key, current)
        else:
            plot_data_cmp1d(last_key, current)


def calculate_mrr(results: (str, List[dict])) -> Dict[str, float]:
    scores = {name: [] for name, _ in results.items()}
    for key, _ in KEYS:
        sorted_results = sorted(
            results.items(),
            key=lambda x:
                x[1][key] if key in x[1] else -float('inf')
        )
        for (score, (name, _)) in enumerate(sorted_results, 1):
            scores[name].append(score)
    return {name: (sum(1 / score for score in scores) / len(KEYS), scores)
            for (name, scores) in scores.items()}


# plot_all(get_averaged_dict())
print("Param set", "MRR", *[name for _, name in KEYS], sep=", ")
print(*[f"\"{name}\", {', '.join([str(x) for x in value])}" for name, value in
      sorted(
          calculate_mrr(
              get_averaged_dict()).items(),
          key=lambda x: -x[1][0])],
      sep="\n")
