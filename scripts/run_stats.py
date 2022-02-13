#!/bin/python3

import subprocess
import itertools
import common
import sys


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


def get_stable_params(n, seed, model, m):
    stable_params = [
        "",
        f"-n:{n}",
        f"-g0:./files/{seed}",
    ]
    if model == "berg":
        stable_params.append("-tu:0.01")
        stable_params.append("-ed:2e4")
    elif model == "two_step":
        stable_params.append(f"-m:{m}")
    return stable_params


def generate_graphs(iters):
    generated = []
    for graph, seed, n, m in common.REAL_GRAPHS:
        for model, variable_params in GENERATORS:
            print(f"Generating graph using {model} with params from {graph}.")

            stable_params = get_stable_params(n, seed, model, m)
            common_prefix = f"-prefix:{graph.replace('.txt', '')}_"
            for i in range(iters):
                stable_params[0] = f"{common_prefix}{i}_"
                generated.extend(
                    generate_graph(
                        model,
                        stable_params,
                        variable_params,
                    )
                )
    return generated


iters = int(sys.argv[1]) if len(sys.argv) > 1 else 1

generated = generate_graphs(iters)
for graph in generated:
    common.calculate_stats(graph)
