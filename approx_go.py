import networkx as nx
import pandas as pd
from collections import Counter


graphfile = "drought_module2.gexf"
dframe = "txt/processed_uniprot.csv"

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

import re

# Add GO terms from the GO columns into a 'GO' attribute as a set
for node in graph.nodes:
    row = frame[frame['Entry'] == node]

    if not row.empty:
        # For every node who is present in the uniprot data, get GO IDs
        biogo = row["Gene ontology (biological process)"].to_string() + row["Gene ontology (GO)"].to_string()
        biogo += row["Gene ontology (molecular function)"].to_string() + row["Gene ontology IDs"].to_string()
        # Extract just the GO terms and form them into a non-redundant set
        res = re.findall("GO:\d\d\d\d\d\d\d", biogo)
        goset = set(res)

        # Add to the 'GO' attribute
        graph.nodes[node]['GO'] = goset


# Iterate over all nodes, getting nodes and their attributes
for node,attrs in graph.nodes(True):
    # For nodes with GO terms in attributes
    if 'GO' in attrs:
        # Get the GO terms
        nodego = attrs['GO']
        combinedgos = []

        # For each neighbor of node
        for neighbor in graph.neighbors(node):
            # Try to get the neighbors GO terms if it exists
            try:
                neighborgo = graph.nodes[neighbor]['GO']
            except KeyError:
                continue

            # Find GO terms that are in common with neighbor and add to combinedgos list
            intersected = nodego.intersection(neighborgo)
            combinedgos += list(intersected)

        # Count the occurences of each GO term in the list and assign 'GO' value as highest value in list
        golabels = Counter(combinedgos)
        mostcomlabel = golabels.most_common(1)
        if len(mostcomlabel) > 0:
            graph.nodes[node]['GO'] = mostcomlabel[0][0]

        else:
            graph.nodes[node]['GO'] = list(graph.nodes[node]['GO'])[0]


