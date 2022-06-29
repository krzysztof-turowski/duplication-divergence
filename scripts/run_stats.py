#!/bin/python3

import subprocess
import itertools
import common
import sys
import asyncio
import argparse


def int_range(start, end, values=3):
    return range(start, end, (end - start + 1) // values)


def float_range(start, end, values=4):
    return [float(start + i * (end - start) / (values + 1))
            for i in range(1, values + 1)]


GENERATORS = [
    ("pure_duplication", [("p", float_range(0, 1, 100))]),
    ("chung_lu", [("p", float_range(0, 1, 10)), ("q", float_range(0, 1, 10))]),
    ("pastor_satorras", [("p", float_range(0, 5, 10)),
     ("r", float_range(0, 1, 10))]),
    ("sticky", [("gamma", float_range(2, 3, 100))]),
    ("ba", [("m", int_range(1, 10, 10))]),
    (
        "copy",
        [("a", int_range(1, 10, 4)), ("b", int_range(
            1, 10, 4)), ("c", float_range(0, 1, 4))],
    ),
    # ("two_step", [("a", float_range(1, 3, 100))]),
    (
        "berg",
        [
            ("ac", float_range(5, 7, 3)),
            ("dr", float_range(0, 2e-4, 3)),
            ("lar", float_range(0, 1e-1, 3)),
            ("ldr", float_range(0, 1e-1, 3)),
        ],
    ),
    ("kumar_linear", [("d", int_range(1, 10, 10)),
     ("alpha", int_range(1, 10, 10))]),
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
                output.find("Generated file: ") + len("Generated file: "): -3
            ]
        )
    return generated


def get_stable_params(n, seed, model, m):
    stable_params = [
        f"-n:{n}",
        f"-g0:./files/{seed}",
    ]
    if model == "berg":
        stable_params.append("-tu:0.1")
        stable_params.append("-ed:2e4")
    elif model == "two_step":
        stable_params.append(f"-m:{m}")
    return stable_params


async def generate_graphs(start, end=None):
    if end is None:
        start, end = 0, start
    promises = []
    for graph, seed, n, m in common.REAL_GRAPHS:
        for model, variable_params in GENERATORS:
            stable_params = get_stable_params(n, seed, model, m)
            common_prefix = f"-prefix:{graph.replace('.txt', '')}_"
            for i in range(start, end):
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
    parser = argparse.ArgumentParser(description='Run statistics for graphs.')
    parser.add_argument(
        'rangeone',
        type=int,
        nargs="?",
        default=None,
        help='first arg of range')
    parser.add_argument(
        'rangetwo',
        type=int,
        nargs="?",
        default=None,
        help='second arg of range')
    parser.add_argument(
        '--mode',
        type=int,
        nargs=1,
        default=3,
        help='mode 1 - fast, 2 - slow, 3 - fast and slow')

    args = parser.parse_args()

    if args.rangeone is None:
        start, end = 0, 1
    elif args.rangetwo is None:
        start, end = 0, args.rangeone
    else:
        start, end = args.rangeone, args.rangetwo

    for i in range(start, end):
        generated = await generate_graphs(i, i + 1)
        processess = await asyncio.gather(*(common.calculate_stats(graph) for graph in generated))
        for process in processess:
            _ = await process.communicate()


asyncio.run(main())
