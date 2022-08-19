import itertools
import random

import networkx

def generate_seed(n, p):
    G = networkx.Graph()
    for i in range(n):
        G.add_node(i, parent = i, ancestor = i)
    for i, j in itertools.combinations(range(n), 2):
        if random.random() <= p:
            G.add_edge(i, j)
    return G

def generate_pure_duplication(G, n, n0, p):
    for i in range(n0, n):
        parent = random.randint(0, i - 1)
        G.add_node(i, parent = parent, ancestor = G.nodes[parent].get("ancestor", parent))
        for _, j in G.edges(parent):
            if random.random() <= p:
                G.add_edge(j, i)
    return G

def generate_pastor_satorras(G, n, n0, p, r):
    for i in range(n0, n):
        parent = random.randint(0, i - 1)
        G.add_node(i, parent = parent, ancestor = G.nodes[parent].get("ancestor", parent))
        for j in range(i):
            if G.has_edge(parent, j):
                if random.random() <= p:
                    G.add_edge(j, i)
            else:
                if random.random() <= r / i:
                    G.add_edge(j, i)
    return G
