import json
from dd_plot_characteristics import calculate_dtw, calculate_vector_distance, read_data
from typing import List, Dict
from itertools import repeat
import os
from multiprocessing import Pool
import math


def compare_graphs(target: dict, synt: dict) -> Dict[str, float]:
    keys = [
        ('log_automorphisms_sparse', "num"),
        ('get_average_shortest_path', "num"),
        ('clustering_coefficient_two', "log"),
        ('clustering_coefficient_three', "log"),
        ('clustering_coefficient_four', "log"),
        ('get_diameter', "num"),
        ('get_degree_distribution', "list"),
        ('betweenness_centrality', "list"),
        ('closeness', "list"),
        ('khop_reachability_2', "list"),
        ('khop_reachability_3', "list"),
        ('khop_reachability_4', "list"),
        ('Graphlets', "vec"),
        ('RGF', "vec")
    ]

    result = {}
    for (key, typ) in keys:
        if key not in synt:
            continue

        if typ == 'num':
            result[key] = (synt[key] - target[key])**2
        if typ == 'log':
            result[key] = (abs(math.log(synt[key]) - math.log(target[key]))
                           if synt[key] != 0 and target[key] != 0
                           else float('inf'))
        elif typ == 'list':
            result[key] = calculate_dtw(synt[key], target[key])
        elif typ == 'vec':
            result[key] = calculate_vector_distance(
                synt[key], target[key])
    return result


real_world_data = {}
for name in os.listdir('results/sieci'):
    with open('results/sieci/' + name) as input_file:
        real_world_data[name[:-4]] = read_data(input_file)


def create_compared_file(path, file):
    print("START:", path + "/" + file)
    splitted = file[:-4].split("%3Amode%3A")
    full_name = splitted[0]
    mode = int(splitted[1]) if len(splitted) > 1 else 3
    base, _, model = full_name.split("_")
    n, n0, model, *params = model[2:].split("-")
    params = [float(p) for p in params]
    n = int(n)
    n0 = int(n0)

    with open(path + "/" + file) as f:
        target = real_world_data[base]
        comparison = compare_graphs(target, read_data(f))
        result = {
            "original_path": path + "/" + file,
            "base_graph_name": base,
            "n": n,
            "n0": n0,
            "params": params,
            "comparison": comparison,
            "mode": mode
        }

    new_path = path.replace("results/generated", "results/compared")
    os.makedirs(new_path, exist_ok=True)
    with open(new_path + "/" + file, "w") as f:
        json.dump(result, f)

    os.remove(path + "/" + file)
    print("FINISH:", new_path + "/" + file)


def main():
    tree = os.walk('results/generated')

    for (path, _, files) in tree:
        with Pool(50) as p:
            p.starmap(create_compared_file, zip(repeat(path), files))


if __name__ == "__main__":
    main()
