import networkx as nx
import pandas as pd
import re

graphfile = "testgraph.gexf"
dframe = "txt/processed_uniprot.csv"
pd.set_option('display.max_colwidth', 1000)

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

def go_label_propagate_dumb(graph):
    # Propagate GO term from node to neighbors
    visited = []
    for node,attrs in graph.nodes(True):
        if 'GO' in attrs.keys():
            label = attrs['GO']

            for neighbor in graph.neighbors(node):
                if neighbor in visited or 'GO' in graph.nodes[neighbor].keys():
                    continue
                else:
                    graph.nodes[neighbor]['GO'] = label
                    visited.append(neighbor)

    return graph

nx.write_gexf(graph,"testgraph2.gexf")