import networkx as nx
import pandas as pd
from collections import Counter

graphfile = "testgraph.gexf"
dframe = "txt/processed_uniprot.csv"

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

# Propagate GO term from node to neighbors

for node,attrs in graph.nodes(True):
    if 'GO' in attrs.keys():
        label = attrs['GO']

        for neighbor in graph.neighbors(node):
            if 'GO' in graph.nodes[neighbor].keys():
                continue
            else:
                graph.nodes[neighbor]['GO'] = label

nx.write_gexf(graph,"testgraph2.gexf")