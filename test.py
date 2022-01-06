import networkx as nx
import pandas as pd
import re

graphfile = "testgraph.gexf"
dframe = "txt/processed_uniprot.csv"
pd.set_option('display.max_colwidth', 1000)

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

# Propagate GO term from node to neighbors
visited = []
for node,attrs in graph.nodes(True):
    if 'GO' in attrs.keys():
        label = attrs['GO']
        # Change node label from uniprot name to protein name
        # If node name in dataframe
        if not frame[frame['Entry'] == node].empty:
            # Get protein name
            name = frame[frame['Entry'] == node]['Protein names'].to_string()
            result = re.findall("\(.*?\)",name)

            # If there's a result, replace label with last result
            if len(result) > 0: graph.nodes[node]['label'] = result[-1][1:-1]

        for neighbor in graph.neighbors(node):
            if neighbor in visited or 'GO' in graph.nodes[neighbor].keys():
                continue
            else:
                graph.nodes[neighbor]['GO'] = label
                visited.append(neighbor)

nx.write_gexf(graph,"testgraph2.gexf")