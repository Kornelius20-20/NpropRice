import networkx as nx
import numpy as np

# load weights
weights = np.load("p.npy")

# Load graph
graph = nx.read_gexf('graph.gexf')

weights *= 63.0/weights.max()
i = 0
for node in graph:
    graph.nodes[node]['weight'] = weights[i]
    i += 1

nx.write_gexf(graph,'test.gexf')