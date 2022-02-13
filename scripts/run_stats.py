#!/bin/python3

import subprocess
import itertools
import common
import sys
import asyncio


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


async def generate_graph(mode, stable_params, variable_params, prefix):
    promises = []
    for param_values in itertools.product(
        *(prange for _, prange in variable_params)
    ):
        args = [
            f"-mode:{mode}",
            prefix,
            *stable_params,
            *(
                f"-{name}:{value}"
                for (name, _), value in zip(variable_params, param_values)
            ),
        ]
        promises.append(
            asyncio.create_subprocess_exec(
                "./dd_generate",
                *args,
                stdout=asyncio.subprocess.PIPE,
            )
        )

    generated = []
    for promise in await asyncio.gather(*promises):
        output, _ = await promise.communicate()
        output = str(output)
        generated.append(
            output[
                output.find("Generated file: ") + len("Generated file: ") : -3
            ]
        )
    return generated


def get_stable_params(n, seed, model, m):
    stable_params = [
        f"-n:{n}",
        f"-g0:./files/{seed}",
    ]
    if model == "berg":
        stable_params.append("-tu:0.01")
        stable_params.append("-ed:2e4")
    elif model == "two_step":
        stable_params.append(f"-m:{m}")
    return stable_params


async def generate_graphs(iters):
    promises = []
    for graph, seed, n, m in common.REAL_GRAPHS:
        for model, variable_params in GENERATORS:
            stable_params = get_stable_params(n, seed, model, m)
            common_prefix = f"-prefix:{graph.replace('.txt', '')}_"
            for i in range(iters):
                promises.append(
                    generate_graph(
                        model,
                        stable_params,
                        variable_params,
                        f"{common_prefix}{i}_",
                    )
                )

    generated = [item for sublist in promises for item in await sublist]

    return generated


async def main():
    iters = int(sys.argv[1]) if len(sys.argv) > 1 else 1

    generated = await generate_graphs(iters)
    for graph in generated:
        await common.calculate_stats(graph)


asyncio.run(main())
