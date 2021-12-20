import os

graph_names = [
    "G-a-thaliana.txt",
    "G-c-elegans.txt",
    "G-college-msg.txt",
    "G-d-melanogaster.txt",
    "G-dynamic-simplewiki-10k.txt",
    "G-hep-th-citations-scc.txt",
    "G-homo-sapiens.txt",
    "G-mus-musculus.txt",
    "G-s-cerevisiae.txt",
    "G-s-pombe.txt",
]

for name in graph_names:
    os.system("nohup ./dd_calculate_real_graph_characteristics " + name + " &")
