#!/bin/python3

import subprocess
import itertools

REAL_GRAPHS = [
    ("G-a-thaliana.txt", "G0-a-thaliana.txt", 9444, 44387),
    ("G-c-elegans.txt", "G0-c-elegans.txt", 3868, 11684),
    ("G-d-melanogaster.txt", "G0-d-melanogaster.txt", 9204, 69560),
    ("G-homo-sapiens.txt", "G0-homo-sapiens.txt", 17294, 313932),
    ("G-mus-musculus.txt", "G0-mus-musculus.txt", 6869, 25269),
    ("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", 6151, 537552),
    ("G-s-pombe.txt", "G0-s-pombe.txt", 4202, 62313),
]


def float_range(start, end, values=10):
    return [float(start + i * (end - start) / 10) for i in range(1, values)]


GENERATORS = [
    ("pure_duplication", [("p", float_range(0, 1))]),
    ("chung_lu", [("p", float_range(0, 1)), ("q", float_range(0, 1))]),
    ("pastor_satorras", [("p", float_range(0, 1)), ("r", float_range(0, 1))]),
    ("sticky", [("gamma", float_range(2, 3))]),
    ("ba", [("m", range(1, 10))]),
    (
        "copy",
        [("a", range(1, 10)), ("b", range(1, 10)), ("c", float_range(0, 1))],
    ),
    ("two_step", [("a", float_range(1, 3))]),
    (
        "berg",
        [
            ("ac", float_range(5, 7)),
            ("dr", float_range(0, 1)),
            ("lar", float_range(0, 1)),
            ("ldr", float_range(0, 1)),
        ],
    ),
    ("kumar_linear", [("d", range(1, 10)), ("alpha", range(1, 10))]),
]

SEEDS = [""]


def generate_graph(mode, stable_params, variable_params):
    generated = []
    for param_values in itertools.product(
        *(prange for _, prange in variable_params)
    ):
        args = [
            f"-mode:{mode}",
            *stable_params,
            *(
                f"-{name}:{value}"
                for (name, _), value in zip(variable_params, param_values)
            ),
        ]
        result = subprocess.run(["./dd_generate", *args], capture_output=True)
        output = str(result.stdout)
        generated.append(
            output[
                output.find("Generated file: ") + len("Generated file: ") : -3
            ]
        )
    return generated


def generate_graphs():
    generated = []
    for graph, seed, n, m in REAL_GRAPHS:
        for model, variable_params in GENERATORS:
            print(f"Generating graph using {model} with params from {graph}.")
            stable_params = [
                f"-prefix:{graph.replace('.txt', '_')}",
                f"-n:{n}",
                f"-g0:./files/{seed}",
            ]
            if model == "berg":
                stable_params.append("-tu:0.01")
                stable_params.append("-ed:2e4")
            elif model == "two_step":
                stable_params.append(f"-m:{m}")
            generated.extend(
                generate_graph(model, stable_params, variable_params)
            )
    return generated


generated = generate_graphs()
for graph in itertools.chain((g[0] for g in REAL_GRAPHS), generated):
    print(f"Calculating stats for graph {graph}.")
    subprocess.run(["./dd_calculate_real_graph_characteristics", graph])
