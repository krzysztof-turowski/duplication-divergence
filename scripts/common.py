import os


REAL_GRAPHS = [
    ("G-a-thaliana.txt", "G0-a-thaliana.txt", 9444, 44387),
    ("G-c-elegans.txt", "G0-c-elegans.txt", 3868, 11684),
    ("G-d-melanogaster.txt", "G0-d-melanogaster.txt", 9204, 69560),
    ("G-homo-sapiens.txt", "G0-homo-sapiens.txt", 17294, 313932),
    ("G-mus-musculus.txt", "G0-mus-musculus.txt", 6869, 25269),
    ("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", 6151, 537552),
    ("G-s-pombe.txt", "G0-s-pombe.txt", 4202, 62313),
]


def prelude():
    os.system("mkdir -p results")


def calculate_stats(graph):
    print(f"Calculating stats for graph {graph}.")
    subprocess.run(["./dd_calculate_real_graph_characteristics", graph])


prelude()
