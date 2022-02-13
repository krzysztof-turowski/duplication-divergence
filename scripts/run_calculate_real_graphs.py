#!/bin/python3

import subprocess
import common


def float_range(start, end, values=10):
    return [float(start + i * (end - start) / 10) for i in range(1, values)]


for graph in common.REAL_GRAPHS:
    common.calculate_stats(graph)
