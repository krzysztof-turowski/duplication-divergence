import asyncio
import os

CUTOFF_SIZE = 10 * 2**20  # 10 MiB

REAL_GRAPHS = [
    # ("G-a-thaliana.txt", "G0-a-thaliana.txt", 9444, 44387),
    ("G-c-elegans.txt", "G0-c-elegans.txt", 3868, 11684),
    # ("G-d-melanogaster.txt", "G0-d-melanogaster.txt", 9204, 69560),
    # ("G-homo-sapiens.txt", "G0-homo-sapiens.txt", 17294, 313932),
    # ("G-mus-musculus.txt", "G0-mus-musculus.txt", 6869, 25269),
    # ("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", 6151, 537552),
    # ("G-s-pombe.txt", "G0-s-pombe.txt", 4202, 62313),
]


def prelude():
    os.system("mkdir -p results")


def calculate_stats(graph):
    size = os.path.getsize("files/" + graph)
    if size <= CUTOFF_SIZE:
        print(f"Starting job for graph {graph}.")
        return asyncio.create_subprocess_exec(
            "./dd_calculate_real_graph_characteristics",
            graph,
            stdout=asyncio.subprocess.PIPE,
        )
    else:
        with open(f"results/{graph}", "w") as f:
            f.write("Ignored: too big.\n")
        return asyncio.create_subprocess_exec(
            "echo", f"Graph {graph} too large. Ignoring")


prelude()
