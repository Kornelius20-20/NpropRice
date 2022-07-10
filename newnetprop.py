import networkx as nx
from scipy.sparse import csr_matrix,csr_array
import os

seedlist = r"txt/string_seeds.txt"

graph = nx.read_gexf("graph_test.gexf")
seeds = []


def weights_from_seeds(graph, seedlist):
    """
    Takes an input graph and creates a list of weights for the graph nodes, where each seed node from a list
    of seeds will get assigned the starting weight value (defaults to 100) and all other nodes will start at
    0

    :param graph: input graph
    :param seedlist: list of seed nodes to assign initial weight to
    :param weight: starting weight value
    :return: List of weights corresponding to each node in the graph
    """

    p0 = []
    # Get seed proteins and assign to them an input weight of 100 and 0 to the rest
    with open(os.path.join(seedlist), 'r') as file:
        seeds = [line.rstrip() for line in file.readlines()]

    propability = float(1 / len(seeds))

    for name in graph.nodes:
        if name in seeds:
            p0.append(propability)
        else:
            p0.append(0)

    p0 = csr_array(p0)

    return p0


def netprop(graph,seeds,iter=0,alpha = 0.01,threshold = 1.0e-6):

    A = nx.adjacency_matrix(graph)

    # Normalize adjacency matrix
    for i in range(A.shape[1]-1):
        colsum = A.data[A.indices[i]:A.indices[i+1]].sum()

        if colsum > 0:
            A.data[A.indices[i]:A.indices[i + 1]] /= colsum

     p0 = weights_from_seeds(graph,seeds)

    p = p0
    if iter > 0:
        for i in range(iter):
            p = alpha*p0 + (1-alpha)*A*p

    else:
        p = alpha * p0 + (1 - alpha) * A * p
        _p = csr_array(p.shape)
        while p - _p > threshold:
            _p = p
            p = alpha * p0 + (1 - alpha) * A * p

    return p