"""
Attempt to code network propagation
"""
import sys

import networkx as nx
import numpy as np
import os

seedlist = r"txt/string_seeds.txt"

def weights_from_seeds(graph,seedlist, weight=100):
    p0 = []
    # Get seed proteins and assign to them an input weight of 100 and 0 to the rest
    with open(os.path.join(seedlist), 'r') as file:
        seeds = [line.rstrip() for line in file.readlines()]

    for name in graph.nodes:
        if name in seeds:
            p0.append(100)
        else:
            p0.append(0)

    p0 = np.array(p0)

    return p0

# Load graph
graph = nx.read_gexf('graph.gexf')

# restart parameter
alpha = 0.5
# Input weights for network
p0 = weights_from_seeds(graph,seedlist,100)


# Get adjacency matrix of graph
A = nx.graphmatrix.adjacency_matrix(graph)

# Create degree matrix of graph
D = np.eye(graph.number_of_nodes())
i = 0
for name,degree in graph.degree:
    D[i,i] = degree
    i += 1

W = np.dot(A,np.linalg.inv(D))

# RWR
p = alpha*np.linalg.inv((1 - (1 - alpha)*W))*p0

